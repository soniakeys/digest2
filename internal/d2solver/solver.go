// Public domain.

// Package d2solver implements the digest2 algorithm.
package d2solver

import (
	"math"

	"github.com/soniakeys/astro"
	"github.com/soniakeys/coord"
	"github.com/soniakeys/digest2/internal/d2bin"
	"github.com/soniakeys/observation"
)

// D2Solver contains data and parameters needed for the digest2 algorithm.
// This includes the Solar System population model, parameters indicating
// which orbit classes to compute scores for, and standard observational
// errors to apply to observations.
type D2Solver struct {
	all, unk      d2bin.Model
	classCompute  []int // from config file
	obsErrMap     map[string]float64
	obsErrDefault float64
}

// New creates a D2Solver object from passed parameters.
func New(all, unk d2bin.Model, classCompute []int,
	obsErrMap map[string]float64, obsErrDefault float64) *D2Solver {
	return &D2Solver{all, unk, classCompute, obsErrMap, obsErrDefault}
}

// Solve runs the digest2 algorithm on a single observational arc.
//
// Rms returned is based on residuals of all observations in the arc
// against fitted linear great circle motion.
// Digest2 scores are returned in the slice classScores.
func (s *D2Solver) Solve(obs *observation.Arc, vMag float64,
	rnd Rand) (rms float64, classScores []Scores) {

	a := s.newArc(obs, vMag, rnd) // create workspace
	a.score()                     // run the algorithm
	return a.rms, a.classScores
}

// Scores is the return type from D2Solver.Solve
type Scores struct {
	Raw, NoId float64
}

// Big messy struct is the workspace for the digest2 algorithm.
// The algorithm operates on a set of observations on a single object.
// --typically a arc, but not required to be all from the same observer.
type arc struct {
	// inputs.  arc constructed with these
	solver *D2Solver
	obs    *observation.Arc
	vMag   float64
	rnd    Rand

	// result values read by digest2.solve
	rms         float64 // rms for arc as a whole
	classScores []Scores

	cs []*classStats

	first, last observation.VObs // obs used for motion vector

	// observational error associated with first, last of motion vector
	firstObsErr, lastObsErr float64
	noObsErr                bool

	// distance independent working variables.  computed once per arc.
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

func (s *D2Solver) newArc(obs *observation.Arc, vMag float64,
	rnd Rand) *arc {

	a := &arc{
		solver:      s,
		obs:         obs,
		vMag:        vMag,
		rnd:         rnd,
		dTag:        make(map[int]bool),
		classScores: make([]Scores, len(s.classCompute)),
		cs:          make([]*classStats, len(s.classCompute)),
	}
	for c, _ := range a.cs {
		a.cs[c] = &classStats{
			dInClass:    make(map[int]bool),
			dNonClass:   make(map[int]bool),
			tagInClass:  make(map[int]bool),
			tagNonClass: make(map[int]bool)}
	}
	return a
}

// per-class workspace, allocated in newArc
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
	Seed(int64)
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

func (a *arc) score() {
	// synthesize or select two observations to determine motion vector
	// this also sets rms values for the two obs and the arc as a whole
	firstRms, lastRms := a.twoObs()
	m1 := a.first.Meas()
	m2 := a.last.Meas()

	// set observational errors to use
	solver := a.solver
	a.firstObsErr = solver.clipErr(firstRms, m1.Qual)
	a.lastObsErr = solver.clipErr(lastRms, m2.Qual)

	// dt derived factors handy in computations
	a.dt = m2.MJD - m1.MJD
	a.invdt = 1 / a.dt
	a.invdtsq = a.invdt * a.invdt

	// solve sun-observer vectors in ecliptic coordinates at observation times
	a.sunObserver0 = a.sov(a.first)
	a.sunObserver1 = a.sov(a.last)

	if a.firstObsErr == 0 && a.lastObsErr == 0 {
		a.noObsErr = true
	}

	a.searchDistance(min_distance)
	a.searchDistance(max_distance)
	a.dRange(min_distance, max_distance, 0)

	var score float64
	for i, s := range a.cs {
		switch d := s.sumAllInClass + s.sumAllNonClass; {
		case d > 0:
			score = 100 * s.sumAllInClass / d
		case solver.classCompute[i] < 2:
			score = 100
		default:
			score = 0
		}
		a.classScores[i].Raw = score

		switch d := s.sumUnkInClass + s.sumUnkNonClass; {
		case d > 0:
			score = 100 * s.sumUnkInClass / d
		case solver.classCompute[i] < 2:
			score = 100
		default:
			score = 0
		}
		a.classScores[i].NoId = score
	}
}

// setSOV sets sun-observer vectors in ecliptic coordinates in the arc struct.
// also sets soe, coe.
func (a *arc) sov(o observation.VObs) (sunObserver coord.Cart) {
	var sunEarth, earthSite coord.Cart
	sunEarth, a.soe, a.coe = astro.Se2000(o.Meas().MJD)
	earthSite = o.EarthObserverVect()
	sunObserver.Sub(&earthSite, &sunEarth)
	sunObserver.RotateX(&sunObserver, a.soe, a.coe)
	return
}

// searchDistance sets up computations and searches angle space
// for a specified distance
//
// Args:
//   a:  arc
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
func (a *arc) searchDistance(d float64) (result bool) {
	a.clearDTags()
	a.dAnyTag = false
	var newTag bool

	for ri := -1.; ri <= 1; ri++ {
		for di := -1.; di <= 1; di++ {
			a.offsetMotionVector(ri, di)
			a.solveDistanceDependentVectors(d)
			if a.searchAngles() {
				newTag = true
			}
			if a.noObsErr {
				return newTag
			}
		}
	}
	return newTag
}

func (a *arc) clearDTags() {
	a.dAnyTag = false
	a.dTag = make(map[int]bool)
	for _, s := range a.cs {
		if len(s.dInClass) > 0 {
			s.dInClass = make(map[int]bool)
		}
		if len(s.dNonClass) > 0 {
			s.dNonClass = make(map[int]bool)
		}
	}
}

func (a *arc) offsetMotionVector(rx, dx float64) {
	// solve unit vectors
	a.observerObjectUnit0 = a.oouv(a.first.Meas(), a.firstObsErr, rx, dx)
	a.observerObjectUnit1 = a.oouv(a.last.Meas(), a.lastObsErr, -rx, -dx)
}

// setOOUV solves observerObject unit vector for sky coordinates.
func (a *arc) oouv(
	sky *observation.VMeas,
	obsErr float64,
	rx, dx float64,
) (observerObjectUnit coord.Cart) {
	sdec, cdec := math.Sincos(sky.Dec + dx*obsErr*.5)
	sra, cra := math.Sincos(sky.Sphr.RA + rx*obsErr*.5*cdec)
	observerObjectUnit = coord.Cart{
		X: cra * cdec,
		Y: sra * cdec,
		Z: sdec,
	}
	observerObjectUnit.RotateX(&observerObjectUnit, a.soe, a.coe)
	return
}

func (a *arc) solveDistanceDependentVectors(d float64) {
	// observerObject0
	a.observerObject0Mag = d
	a.observerObject0 = a.observerObjectUnit0
	a.observerObject0.MulScalar(&a.observerObject0, d)

	// sunObject0
	a.sunObject0.Add(&a.sunObserver0, &a.observerObject0)
	a.sunObject0MagSq = a.sunObject0.Square()
	a.sunObject0Mag = math.Sqrt(a.sunObject0MagSq)

	// observer1Object0
	a.observer1Object0.Sub(&a.sunObject0, &a.sunObserver1)
	a.observer1Object0MagSq = a.observer1Object0.Square()
	a.observer1Object0Mag = math.Sqrt(a.observer1Object0MagSq)

	// solve H mag
	a.hmag = astro.HMag(&a.observerObject0, &a.sunObject0,
		a.vMag, a.observerObject0Mag, a.sunObject0Mag)
	a.hmagBin = d2bin.H(a.hmag)
}

// dRange explores possible orbit space
//
// recursive arc method
//
// Args:
//   a:          arc
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
func (a *arc) dRange(d1, d2 float64, age int) {
	dmid := (d1 + d2) * .5

	if a.searchDistance(dmid) || d2-d1 > minDistanceStep {
		a.dRange(d1, dmid, 0)
		a.dRange(dmid, d2, 0)
		return
	}

	if age < ageLimit {
		a.dRange(d1, dmid, age+1)
		a.dRange(dmid, d2, age+1)
	}
}

func (a *arc) searchAngles() bool {
	ang1, ang2, ok := a.solveAngleRange()
	if !ok {
		return false
	}

	a.aRange(ang1, ang2, 0)

	if !a.dAnyTag {
		return false
	}

	var newTag bool

	for i, dt := range a.dTag {
		if dt {
			for cx, c := range a.solver.classCompute {
				s := a.cs[cx]
				if s.dInClass[i] && !s.tagInClass[i] {
					newTag = true
					s.tagInClass[i] = true
					s.sumAllInClass += a.solver.all.Class[c][i]
					s.sumUnkInClass += a.solver.unk.Class[c][i]
				}
				if s.dNonClass[i] && !s.tagNonClass[i] {
					newTag = true
					s.tagNonClass[i] = true
					s.sumAllNonClass +=
						a.solver.all.SS[i] - a.solver.all.Class[c][i]
					s.sumUnkNonClass +=
						a.solver.unk.SS[i] - a.solver.unk.Class[c][i]
				}
			}
		}
	}
	return newTag
}

// parabolic limits
func (a *arc) solveAngleRange() (ang1, ang2 float64, ok bool) {
	// solve angle range at this distance
	th := a.observer1Object0.Dot(&a.observerObjectUnit1) /
		a.observer1Object0Mag

	a.tz = math.Acos(th)

	aa := a.invdtsq
	bb := -2 * a.observer1Object0Mag * th * aa
	cc := a.observer1Object0MagSq*aa - 2*astro.U/a.sunObject0Mag
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
		nns := d2s + a.observer1Object0MagSq -
			2*d2*a.observer1Object0Mag*th
		nn := math.Sqrt(nns)
		ca := (nns + a.observer1Object0MagSq - d2s) /
			(2 * nn * a.observer1Object0Mag)
		sa := d2 * math.Sin(a.tz) / nn
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
//   a:  arc with distance setup already done
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
func (a *arc) aRange(ang1, ang2 float64, age int) {
	d3 := (ang2 - ang1) / 3
	mid := ang1 + d3 + d3*a.rnd.Float64()

	if a.tagAngle(mid) || d3 > minAngleStep {
		a.aRange(ang1, mid, 0)
		a.aRange(mid, ang2, 0)
		return
	}

	if age < ageLimit {
		a.aRange(ang1, mid, age+1)
		a.aRange(mid, ang2, age+1)
	}
}

// tagAngle processes a single distance-angle combination.
//
// Args:
//   a:  arc with distance setup already done.
//   an:  angle for orbit solution
//
// Notes:
//   solves orbit for passed angle, converts to bin indicies, sets bin tag
//   and updates tag count.
func (a *arc) tagAngle(an float64) bool {
	// compute object velocity scaled by gravitational constant
	a.v = a.observerObjectUnit1
	s := a.observer1Object0Mag * math.Sin(an) / math.Sin(math.Pi-an-a.tz)
	a.v.MulScalar(&a.v, s)
	a.v.Sub(&a.v, &a.observer1Object0)
	a.v.MulScalar(&a.v, a.invdt*astro.InvK)

	// compute (some) Keplarian elements
	sa, e, i, hv := astro.AeiHv(&a.sunObject0, &a.v, a.sunObject0Mag)
	if hv == nil {
		return false
	}
	a.hv = *hv

	q := sa * (1 - e)
	iq, ie, ii, inModel := d2bin.Qei(q, e, i)
	if !inModel {
		return false
	}
	ih := a.hmagBin
	bx := d2bin.Mx(iq, ie, ii, ih)

	// meaning: some class was newly tagged for this bin at this distance.
	// used as function return value, see below
	var newTag bool

	for cx, c := range a.solver.classCompute {
		s := a.cs[cx]
		if d2bin.CList[c].IsClass(q, e, i, a.hmag) {
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
		a.dAnyTag = true

		// meaning: this bin intersects 2d surface at this distance
		a.dTag[bx] = true
	}
	// true return means "we're finding stuff, keep searching more
	// angles at this distance"
	return newTag
}
