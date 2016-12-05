// Public domain.

package d2solver

import (
	"math"

	"github.com/soniakeys/coord"
	"github.com/soniakeys/lmfit"
	"github.com/soniakeys/observation"
	"github.com/soniakeys/unit"
)

// twoObs computes two observations suitable for computing motion vector.
// in the process it also computes an rms of great circle residuals.
//
// at least two obs must be in  a.obs.  a.rms is set to the rms over
// the whole tracklet.  Return values firstRms and lastRms will be
// non-zero if a GC fit applies to the corresponding motion vector endpoint.
func (a *arc) twoObs() (firstRms, lastRms unit.Angle) {
	obs := a.obs.Obs
	// default obs
	a.first = obs[0]
	a.last = obs[len(obs)-1]
	if len(obs) == 2 {
		// simplest case, return the only two points given, rms = 0
		return
	}
	// > 2 obs, do a great circle fit over all obs to get rms return value.
	// Fit may also be used in some cases for synthesizing observations.
	t := make([]float64, len(obs))
	s := make(coord.EquaS, len(obs))
	for i, o := range obs {
		m := o.Meas()
		t[i] = m.MJD
		s[i] = m.Equa
	}
	lmf := lmfit.New(t, s)
	a.rms = lmf.Rms() // set tracklet rms

	// scan site infomation. determine if all obs are from same site and
	// any space based observations are present
	allSameSite := true
	var ok, spaceBased bool
	var site0 *observation.SiteObs
	var par0 *observation.ParallaxConst
	if site0, ok = a.first.(*observation.SiteObs); ok {
		par0 = site0.Par
	} else {
		spaceBased = true
	}
	for _, o := range obs[1:] {
		if s, ok := o.(*observation.SiteObs); ok {
			if s.Par != par0 {
				allSameSite = false
			}
		} else {
			spaceBased = true
		}
	}
	// find observations near 17th and 83rd percentile.
	// percentile is used rather than percent to allow for non-uniformly
	// distributed observation times.
	whole, fs := math.Modf(float64(len(obs)-1) / 6)
	is := int(whole)

	// not always sensible to interpolate when space based obs are present,
	// especially if observer position is ever taken into account.
	// in this case just return obs at 17th and 83rd percentile.
	// leave first, last rms = 0 to take the default rms for the observatory.
	if spaceBased {
		a.first = obs[is]
		a.last = obs[len(obs)-1-is]
		return
	}
	// from here on, obs are known to be all ground based.
	// type assertions to *observation.SiteObs are guaranteed to work.
	// site0 is already computed.  siteLast is handy now.
	siteLast := obs[len(obs)-1].(*observation.SiteObs)

	// compute times t17 and t83 at these points of interest.
	// the times will be used in a few different ways.
	t17 := obs[is].Meas().MJD
	t17 += (obs[is+1].Meas().MJD - t17) * fs
	is = len(obs) - 1 - is
	t83 := obs[is].Meas().MJD
	t83 -= (t83 - obs[is-1].Meas().MJD) * fs

	// next case still fairly simple:  single site, obs < 3 hrs.
	//    => use gc fit of whole obs and synthesize obs at the 17th
	//    and 83rd percentile times.  first, last rms same as arc.
	if allSameSite && siteLast.MJD-site0.MJD < .125 {
		so := &observation.SiteObs{
			VMeas: site0.VMeas,
			Par:   par0,
		}
		so.VMeas.MJD = t17
		so.VMeas.Equa = *lmf.Pos(t17)
		a.first = so

		so = &observation.SiteObs{
			VMeas: siteLast.VMeas,
			Par:   par0,
		}
		so.VMeas.MJD = t83
		so.VMeas.Equa = *lmf.Pos(t83)
		a.last = so

		return a.rms, a.rms
	}

	// remaining case is involved.  not appropriate to gc fit the entire
	// obs, but probably possible to derive better endpoints than just
	// first and last observation.

	// first step is to split off initial and final assumed tracklets, each
	// tracklet being all the same code and < 3hrs, but otherwise being as
	// long as possible.  if tracklets use all obs and are both of same code,
	// they should be as equal in dt as possible.
	o1 := 0
	o2 := len(obs) - 1
	site1 := site0 // new variables just for naming consistency
	site2 := siteLast
	par1 := par0
	par2 := site2.Par
	t1 := site1.MJD
	t2 := site2.MJD
	for {
		s1next := obs[o1+1].(*observation.SiteObs)
		dt1 := s1next.MJD - t1
		if s1next.Par != par1 || dt1 > .125 {
			// initial obs is done, just try to extend final obs
			for {
				o := o2 - 1
				if o == o1 {
					break
				}
				s2prev := obs[o].(*observation.SiteObs)
				if s2prev.Par != par2 || t2-s2prev.MJD > .125 {
					break
				}
				o2 = o
			}
			break
		}
		s2prev := obs[o2-1].(*observation.SiteObs)
		dt2 := t2 - s2prev.MJD
		if s2prev.Par != par2 || dt2 > .125 {
			// final obs is done, just try to extend initial obs
			for {
				o := o1 + 1
				if o == o2 {
					break
				}
				s1next := obs[o].(*observation.SiteObs)
				if s1next.Par != par1 || s1next.MJD-t1 > .125 {
					break
				}
				o1 = o
			}
			break
		}

		if dt1 < dt2 {
			o1++
		} else {
			o2--
		}

		if o2 == o1+1 {
			break
		}
	}
	// handle each tracklet
	a.first, firstRms = oneObs(0, o1, o2 == o1+1, t17, obs)
	a.last, lastRms = oneObs(o2, len(obs)-1, o2 == o1+1, t83, obs)
	return
}

func oneObs(
	o1, o2 int,
	obssUseAllObs bool,
	pt float64,
	obs []observation.VObs,
) (result *observation.SiteObs, rms unit.Angle) {
	// a default result. (again, obs is guaranteed to be all ground based)
	result = obs[o1].(*observation.SiteObs)

	// case 1.  simple, only a single obs for the "obs".
	if o1 == o2 {
		return
	}

	// case 2.  two points for the obs.
	if o1 == o2-1 {
		// case 2.1.  fairly simple:  if there's stuff between the obss
		// and the target percentile (17 or 83) is off of this obs,
		// just return the end of the obs.
		if !obssUseAllObs {
			var dt float64
			end := result
			if o1 == 0 {
				end = obs[1].(*observation.SiteObs)
				dt = end.MJD - pt
			} else {
				dt = pt - end.MJD
			}
			if dt < 0 {
				return end, 0
			}
		}

		// case 2.2.  linearly interpolate along the great circle
		// connecting the points.
		r2 := *(obs[o2].(*observation.SiteObs))
		t := make([]float64, 2)
		s := make(coord.EquaS, 2)
		t[0] = result.MJD
		s[0] = result.Equa
		t[1] = r2.MJD
		s[1] = r2.Equa
		lmf := lmfit.New(t, s)
		if obssUseAllObs {
			// gc midpoint
			r2.MJD = (result.MJD + r2.MJD) * .5
		} else {
			// gc at pt
			r2.MJD = pt
		}
		r2.Equa = *lmf.Pos(r2.MJD)
		return &r2, 0
	}

	// case 3.  3 or more points in obs.
	// synthesize obs at point along gc
	var tr float64
	if obssUseAllObs {
		// median time of obs
		is := (o1 + o2) / 2
		tr = obs[is].(*observation.SiteObs).MJD
		if is+is < o1+o2 {
			tr = (tr + obs[is+1].(*observation.SiteObs).MJD) * .5
		}
	} else {
		var dt float64
		if o1 == 0 {
			result = obs[o2].(*observation.SiteObs)
			dt = result.MJD - pt
		} else {
			dt = pt - result.MJD
		}
		if dt < 0 {
			tr = result.MJD
		} else {
			tr = pt
		}
	}

	// gc fit, result is computed obs at time tr
	np := o2 - o1 + 1
	t := make([]float64, np)
	s := make(coord.EquaS, np)
	for i, v := range obs[o1 : o2+1] {
		m := v.(*observation.SiteObs)
		t[i] = m.MJD
		s[i] = m.Equa
	}
	lmf := lmfit.New(t, s)
	r2 := *result
	r2.MJD = tr
	r2.Equa = *lmf.Pos(tr)
	return &r2, lmf.Rms() // return an rms for this obs
}
