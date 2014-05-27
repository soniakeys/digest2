// Package d2solver implements the digest2 algorithm.
package d2solver

import (
	"math"

	"digest2/astro"
	"digest2/bin"
	"digest2/coord"
	"digest2/obs"
)

// D2Solver contains data and parameters needed for the digest2 algorithm.
// This includes the Solar System population model, parameters indicating
// which orbit classes to compute scores for, and standard observational
// errors to apply to observations.
type D2Solver struct {
	all, unk      bin.Model
	classCompute  []int // from config file
	obsErrMap     map[string]float64
	obsErrDefault float64
}

// New creates a D2Solver object from passed parameters.
func New(all, unk bin.Model, classCompute []int,
	obsErrMap map[string]float64, obsErrDefault float64) *D2Solver {
	return &D2Solver{all, unk, classCompute, obsErrMap, obsErrDefault}
}

// Solve runs the digest2 algorithm on a single tracklet.
//
// Rms returned is based on residuals of all observations in the tracklet
// against fitted linear great circle motion.
// Digest2 scores are returned in the slice classScores.
func (s *D2Solver) Solve(otk *obs.Tracklet, vmag float64,
	rnd Rand) (rms float64, classScores []Scores) {

	tk := s.newTracklet(otk, vmag, rnd) // create workspace
	tk.score()                          // run the algorithm
	return tk.rms, tk.classScores
}

// Scores is the return type from D2Solver.Solve
type Scores struct {
	Raw, NoId float64
}

// big messy struct is the workspace for the digest2 algorithm.
type tracklet struct {
	// inputs.  tracket constructed with these
	solver *D2Solver
	otk    *obs.Tracklet
	vmag   float64
	rnd    Rand

	// result values read by digest2.solve
	rms         float64 // rms for tracklet as a whole
	classScores []Scores

	cs []*classStats

	first, last obs.VObs // obs used for motion vector

	// observational error associated with first, last of motion vector
	firstObsErr, lastObsErr float64
	noObsErr                bool

	// distance independent working variables.  computed once per tracklet.
	dt, invdt, invdtsq  float64
	soe, coe            float64
	sunObserver0        coord.Cart // vector at time t0
	sunObserver1        coord.Cart // vector at time t1
	observerObjectUnit0 coord.Cart // vector at time t0
	observerObjectUnit1 coord.Cart // vector at time t1

	// distance dependent working variables
	observerObject0       coord.Cart
	observerObject0Mag    float64 // vector magnitude (not brightness)
	sunObject0            coord.Cart
	sunObject0Mag         float64 // vector magnitude (not brightness)
	sunObject0MagSq       float64
	observer1Object0      coord.Cart
	observer1Object0Mag   float64 // vector magnitude (not brightness)
	observer1Object0MagSq float64

	tz, hmag float64
	hmagBin  int

	dAnyTag bool
	dTag    map[int]bool

	// angle dependent working variables.  recomputed many times.
	// local variables would read more easily, but structs are here
	// to reduce garbage
	hv, v coord.Cart
}

func (s *D2Solver) newTracklet(otk *obs.Tracklet, vmag float64,
	rnd Rand) *tracklet {

	t := &tracklet{
		solver:      s,
		otk:         otk,
		vmag:        vmag,
		rnd:         rnd,
		dTag:        make(map[int]bool),
		classScores: make([]Scores, len(s.classCompute)),
		cs:          make([]*classStats, len(s.classCompute)),
	}
	for c, _ := range t.cs {
		t.cs[c] = &classStats{
			dInClass:    make(map[int]bool),
			dNonClass:   make(map[int]bool),
			tagInClass:  make(map[int]bool),
			tagNonClass: make(map[int]bool)}
	}
	return t
}

// per-class workspace, allocated in newTracklet
type classStats struct {
	tagInClass, tagNonClass       map[int]bool
	sumAllInClass, sumAllNonClass float64
	sumUnkInClass, sumUnkNonClass float64
	dInClass, dNonClass           map[int]bool
}

// Rand is an interface allowing the random number generator used by the
// solver to be be swapped between the standard library function and
// a simple linear congruential generator that can be easily ported to
// other languages.  The lcg can then be used to validate that ports of
// digest2 produce the same answers.
type Rand interface {
	Float64() float64
}

// some parameters for the algorithm
const (
	min_distance    = .05 // AU
	max_distance    = 100 // AU
	minDistanceStep = .2  // AU
	minAngleStep    = .1  // radians
	// an experiment.  values > 1 proved expensive for little benefit.
	ageLimit = 1
)

func (tk *tracklet) score() {
	// synthesize or select two observations to determine motion vector
	// this also sets rms values for the two obs and the tracklet as a whole
	firstRms, lastRms := tk.twoObs()
	m1 := tk.first.Meas()
	m2 := tk.last.Meas()

	// set observational errors to use
	solver := tk.solver
	tk.firstObsErr = solver.clipErr(firstRms, m1.Qual)
	tk.lastObsErr = solver.clipErr(lastRms, m2.Qual)

	// dt derived factors handy in computations
	tk.dt = m2.Mjd - m1.Mjd
	tk.invdt = 1 / tk.dt
	tk.invdtsq = tk.invdt * tk.invdt

	// solve sun-observer vectors in ecliptic coordinates at observation times
	tk.sunObserver0 = tk.sov(tk.first)
	tk.sunObserver1 = tk.sov(tk.last)

	if tk.firstObsErr == 0 && tk.lastObsErr == 0 {
		tk.noObsErr = true
	}

	tk.searchDistance(min_distance)
	tk.searchDistance(max_distance)
	tk.dRange(min_distance, max_distance, 0)

	var score float64
	for i, s := range tk.cs {
		switch d := s.sumAllInClass + s.sumAllNonClass; {
		case d > 0:
			score = 100 * s.sumAllInClass / d
		case solver.classCompute[i] < 2:
			score = 100
		default:
			score = 0
		}
		tk.classScores[i].Raw = score

		switch d := s.sumUnkInClass + s.sumUnkNonClass; {
		case d > 0:
			score = 100 * s.sumUnkInClass / d
		case solver.classCompute[i] < 2:
			score = 100
		default:
			score = 0
		}
		tk.classScores[i].NoId = score
	}
}

// setSOV sets sun-observer vectors in ecliptic coordinates in the tracklet
// struct.  also sets soe, coe.
func (tk *tracklet) sov(o obs.VObs) (sunObserver coord.Cart) {
	var sunEarth, earthSite coord.Cart
	sunEarth, tk.soe, tk.coe = astro.Se2000(o.Meas().Mjd)
	earthSite = o.EarthObserverVect()
	sunObserver.Sub(&earthSite, &sunEarth)
	sunObserver.RotateX(&sunObserver, tk.soe, tk.coe)
	return
}

// searchDistance sets up computations and searches angle space
// for a specified distance
//
// Args:
//   tk:  tracklet
//   d:  distance for computations
//
//Algorithm:
//   - setup stuff common to all possible orbits at this distance.  this
//     includes some vectors, H magnitude, and limits on possible angles.
//     Also tags are reset for the distance.
//
//   - search angle space, which results in tags being set.
//
//   - note newly set tags, and accumulate population totals needed for
//     final computation of scores.
//
//Returns:
//   true if any new bins were tagged.
//   false if this distance saw no bins that hadn't been seen
//    at other distances.
func (tk *tracklet) searchDistance(d float64) (result bool) {
	tk.clearDTags()
	tk.dAnyTag = false
	var newTag bool

	for ri := -1.; ri <= 1; ri++ {
		for di := -1.; di <= 1; di++ {
			tk.offsetMotionVector(ri, di)
			tk.solveDistanceDependentVectors(d)
			if tk.searchAngles() {
				newTag = true
			}
			if tk.noObsErr {
				return newTag
			}
		}
	}
	return newTag
}

func (tk *tracklet) clearDTags() {
	tk.dAnyTag = false
	tk.dTag = make(map[int]bool)
	for _, s := range tk.cs {
		if len(s.dInClass) > 0 {
			s.dInClass = make(map[int]bool)
		}
		if len(s.dNonClass) > 0 {
			s.dNonClass = make(map[int]bool)
		}
	}
}

func (tk *tracklet) offsetMotionVector(rx, dx float64) {
	// solve unit vectors
	tk.observerObjectUnit0 = tk.oouv(tk.first.Meas(), tk.firstObsErr, rx, dx)
	tk.observerObjectUnit1 = tk.oouv(tk.last.Meas(), tk.lastObsErr, -rx, -dx)
}

// setOOUV solves observerObject unit vector for sky coordinates.
func (tk *tracklet) oouv(sky *obs.VMeas, obsErr float64, rx, dx float64) (observerObjectUnit coord.Cart) {
	sdec, cdec := math.Sincos(sky.Dec + dx*obsErr*.5)
	sra, cra := math.Sincos(sky.Sphr.Ra + rx*obsErr*.5*cdec)
	observerObjectUnit = coord.Cart{
		X: cra * cdec,
		Y: sra * cdec,
		Z: sdec,
	}
	observerObjectUnit.RotateX(&observerObjectUnit, tk.soe, tk.coe)
	return
}

func (tk *tracklet) solveDistanceDependentVectors(d float64) {
	// observerObject0
	tk.observerObject0Mag = d
	tk.observerObject0 = tk.observerObjectUnit0
	tk.observerObject0.MulScalar(&tk.observerObject0, d)

	// sunObject0
	tk.sunObject0.Add(&tk.sunObserver0, &tk.observerObject0)
	tk.sunObject0MagSq = tk.sunObject0.Square()
	tk.sunObject0Mag = math.Sqrt(tk.sunObject0MagSq)

	// observer1Object0
	tk.observer1Object0.Sub(&tk.sunObject0, &tk.sunObserver1)
	tk.observer1Object0MagSq = tk.observer1Object0.Square()
	tk.observer1Object0Mag = math.Sqrt(tk.observer1Object0MagSq)

	// solve H mag
	tk.hmag = astro.HMag(&tk.observerObject0, &tk.sunObject0,
		tk.vmag, tk.observerObject0Mag, tk.sunObject0Mag)
	tk.hmagBin = bin.H(tk.hmag)
}

// dRange explores possible orbit space
//
// recursive tracklet method 
//
// Args:
//   tk:          tracklet
//   d1, d2:      distance limits
//   age:         a little history.  0 means a bin was just tagged.  1 is added
//                at each level of recursion.
//
// Algorithm:
//   much like aRange
//  
//   - if new bins were found at mid point, recurse both halves
//
//   - if range is big, recurse anyway
//
//   - if young, recurse
//  
func (tk *tracklet) dRange(d1, d2 float64, age int) {
	dmid := (d1 + d2) * .5

	if tk.searchDistance(dmid) || d2-d1 > minDistanceStep {
		tk.dRange(d1, dmid, 0)
		tk.dRange(dmid, d2, 0)
		return
	}

	if age < ageLimit {
		tk.dRange(d1, dmid, age+1)
		tk.dRange(dmid, d2, age+1)
	}
}

func (tk *tracklet) searchAngles() bool {
	ang1, ang2, ok := tk.solveAngleRange()
	if !ok {
		return false
	}

	tk.aRange(ang1, ang2, 0)

	if !tk.dAnyTag {
		return false
	}

	var newTag bool

	for i, dt := range tk.dTag {
		if dt {
			for cx, c := range tk.solver.classCompute {
				s := tk.cs[cx]
				if s.dInClass[i] && !s.tagInClass[i] {
					newTag = true
					s.tagInClass[i] = true
					s.sumAllInClass += tk.solver.all.Class[c][i]
					s.sumUnkInClass += tk.solver.unk.Class[c][i]
				}
				if s.dNonClass[i] && !s.tagNonClass[i] {
					newTag = true
					s.tagNonClass[i] = true
					s.sumAllNonClass +=
						tk.solver.all.SS[i] - tk.solver.all.Class[c][i]
					s.sumUnkNonClass +=
						tk.solver.unk.SS[i] - tk.solver.unk.Class[c][i]
				}
			}
		}
	}
	return newTag
}

// parabolic limits
func (tk *tracklet) solveAngleRange() (ang1, ang2 float64, ok bool) {
	// solve angle range at this distance
	th := tk.observer1Object0.Dot(&tk.observerObjectUnit1) /
		tk.observer1Object0Mag

	tk.tz = math.Acos(th)

	aa := tk.invdtsq
	bb := -2 * tk.observer1Object0Mag * th * aa
	cc := tk.observer1Object0MagSq*aa - 2*astro.U/tk.sunObject0Mag
	dsc := bb*bb - 4*aa*cc

	// use ! > to catch cases where dsc is Inf or NaN at this point.
	if !(dsc > 0) {
		return
	}

	sd := math.Sqrt(dsc)
	sd1 := -sd
	inv2aa := .5 / aa

	for {
		d2 := (-bb + sd1) * inv2aa
		d2s := d2 * d2
		nns := d2s + tk.observer1Object0MagSq -
			2*d2*tk.observer1Object0Mag*th
		nn := math.Sqrt(nns)
		ca := (nns + tk.observer1Object0MagSq - d2s) /
			(2 * nn * tk.observer1Object0Mag)
		sa := d2 * math.Sin(tk.tz) / nn
		ang2 = 2 * math.Atan2(sa, 1+ca)

		if sd1 == sd {
			break
		}

		ang1 = ang2
		sd1 = sd
	}
	return ang1, ang2, true
}

// aRange explores the space between two angles (at a set distance)
//
// Args:
//   tk:  tracklet with distance setup already done
//   ang1, ang2:  search boundaries.
//
// Algorithm:
//   - pick mid point.  a little jiggle is thrown in to help find new bins at
//     closely adjacent distances.
//
//   - solve angle at mid point.  if it resulted in tagging a new bin, recurse
//     both halves.
//
//   - if passed angle range is sufficiently large, recurse.
//
//   - age criterion: if a passed angle recently yielded a new bin, recurse.
func (tk *tracklet) aRange(ang1, ang2 float64, age int) {
	d3 := (ang2 - ang1) / 3
	mid := ang1 + d3 + d3*tk.rnd.Float64()

	if tk.tagAngle(mid) || d3 > minAngleStep {
		tk.aRange(ang1, mid, 0)
		tk.aRange(mid, ang2, 0)
		return
	}

	if age < ageLimit {
		tk.aRange(ang1, mid, age+1)
		tk.aRange(mid, ang2, age+1)
	}
}

// tagAngle processes a single distance-angle combination.
//
// Args:
//   tk:  tracklet with distance setup already done.
//   an:  angle for orbit solution
//
// Notes:
//   solves orbit for passed angle, converts to bin indicies, sets bin tag
//   and updates tag count.
func (tk *tracklet) tagAngle(an float64) bool {
	// compute object velocity scaled by gravitational constant
	tk.v = tk.observerObjectUnit1
	s := tk.observer1Object0Mag * math.Sin(an) / math.Sin(math.Pi-an-tk.tz)
	tk.v.MulScalar(&tk.v, s)
	tk.v.Sub(&tk.v, &tk.observer1Object0)
	tk.v.MulScalar(&tk.v, tk.invdt*astro.InvK)

	// compute (some) Keplarian elements
	a, e, i, hv := astro.AeiHv(&tk.sunObject0, &tk.v, tk.sunObject0Mag)
	if hv == nil {
		return false
	}
	tk.hv = *hv

	q := a * (1 - e)
	iq, ie, ii, inModel := bin.Qei(q, e, i)
	if !inModel {
		return false
	}
	ih := tk.hmagBin
	bx := bin.Mx(iq, ie, ii, ih)

	// meaning: some class was newly tagged for this bin at this distance.
	// used as function return value, see below
	var newTag bool

	for cx, c := range tk.solver.classCompute {
		s := tk.cs[cx]
		if bin.CList[c].IsClass(q, e, i, tk.hmag) {
			if !s.dInClass[bx] {
				s.dInClass[bx] = true
				newTag = true
			}
		} else {
			if !s.dNonClass[bx] {
				s.dNonClass[bx] = true
				newTag = true
			}
		}
	}
	if newTag {
		// meaning: some orbit was found at this distance.
		// will generally only be false if beyond parabolic limit.
		tk.dAnyTag = true

		// meaning: this bin intersects 2d surface at this distance
		tk.dTag[bx] = true
	}
	// true return means "we're finding stuff, keep searching more
	// angles at this distance"
	return newTag
}
