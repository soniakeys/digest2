package d2solver

import (
	"math"

	"coord"
	"digest2/obs"
	"lmfit"
)

// twoObs computes two observations suitable for computing motion vector.
// in the process it also computes an rms of great circle residuals.
//
// at least two obs must be in  tk.obs.  tk.rms is set to the rms over
// the whole tracklet.  Return values firstRms and lastRms will be
// non-zero if a GC fit applies to the corresponding motion vector endpoint.
func (tk *tracklet) twoObs() (firstRms, lastRms float64) {
	arc := tk.otk.Obs
	// default obs
	tk.first = arc[0]
	tk.last = arc[len(arc)-1]
	if len(arc) == 2 {
		// simplest case, return the only two points given, rms = 0
		return
	}
	// > 2 obs, do a great circle fit over all obs to get rms return value.
	// Fit may also be used in some cases for synthesizing observations.
	t := make([]float64, len(arc))
	s := make(coord.SphrS, len(arc))
	for i, o := range arc {
		m := o.Meas()
		t[i] = m.Mjd
		s[i] = m.Sphr
	}
	lmf := lmfit.New(t, s)
	tk.rms = lmf.Rms() // set tracklet rms

	// scan site infomation. determine if all obs are from same site and
	// any space based observations are present
	allSameSite := true
	var ok, spaceBased bool
	var site0 *obs.SiteObs
	var par0 *obs.ParallaxConst
	if site0, ok = tk.first.(*obs.SiteObs); ok {
		par0 = site0.Par
	} else {
		spaceBased = true
	}
	for _, o := range arc[1:] {
		if s, ok := o.(*obs.SiteObs); ok {
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
	whole, fs := math.Modf(float64(len(arc)-1) / 6)
	is := int(whole)

	// not always sensible to interpolate when space based obs are present,
	// especially if observer position is ever taken into account.
	// in this case just return obs at 17th and 83rd percentile.
	// leave first, last rms = 0 to take the default rms for the observatory.
	if spaceBased {
		tk.first = arc[is]
		tk.last = arc[len(arc)-1-is]
		return
	}
	// from here on, obs are known to be all ground based.
	// type assertions to *obs.SiteObs are guaranteed to work.
	// site0 is already computed.  siteLast is handy now.
	siteLast := arc[len(arc)-1].(*obs.SiteObs)

	// compute times t17 and t83 at these points of interest.
	// the times will be used in a few different ways.
	t17 := arc[is].Meas().Mjd
	t17 += (arc[is+1].Meas().Mjd - t17) * fs
	is = len(arc) - 1 - is
	t83 := arc[is].Meas().Mjd
	t83 -= (t83 - arc[is-1].Meas().Mjd) * fs

	// next case still fairly simple:  single site, arc < 3 hrs.
	//    => use gc fit of whole arc and synthesize obs at the 17th
	//    and 83rd percentile times.  first, last rms same as tracklet.
	if allSameSite && siteLast.Mjd-site0.Mjd < .125 {
		so := &obs.SiteObs{
			VMeas: site0.VMeas,
			Par:   par0,
		}
		so.VMeas.Mjd = t17
		so.VMeas.Sphr = *lmf.Pos(t17)
		tk.first = so

		so = &obs.SiteObs{
			VMeas: siteLast.VMeas,
			Par:   par0,
		}
		so.VMeas.Mjd = t83
		so.VMeas.Sphr = *lmf.Pos(t83)
		tk.last = so

		return tk.rms, tk.rms
	}

	// remaining case is involved.  not appropriate to gc fit the entire
	// arc, but probably possible to derive better endpoints than just
	// first and last obs.

	// first step is to split off initial and final arcs, each arc being
	// all the same code and < 3hrs, but otherwise being as long as possible.
	// if arcs use all obs and are both of same code, they should be as equal
	// in dt as possible.
	o1 := 0
	o2 := len(arc) - 1
	site1 := site0 // new variables just for naming consistency
	site2 := siteLast
	par1 := par0
	par2 := site2.Par
	t1 := site1.Mjd
	t2 := site2.Mjd
	for {
		s1next := arc[o1+1].(*obs.SiteObs)
		dt1 := s1next.Mjd - t1
		if s1next.Par != par1 || dt1 > .125 {
			// initial arc is done, just try to extend final arc
			for {
				o := o2 - 1
				if o == o1 {
					break
				}
				s2prev := arc[o].(*obs.SiteObs)
				if s2prev.Par != par2 || t2-s2prev.Mjd > .125 {
					break
				}
				o2 = o
			}
			break
		}
		s2prev := arc[o2-1].(*obs.SiteObs)
		dt2 := t2 - s2prev.Mjd
		if s2prev.Par != par2 || dt2 > .125 {
			// final arc is done, just try to extend initial arc
			for {
				o := o1 + 1
				if o == o2 {
					break
				}
				s1next := arc[o].(*obs.SiteObs)
				if s1next.Par != par1 || s1next.Mjd-t1 > .125 {
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
	// handle each arc
	tk.first, firstRms = oneObs(0, o1, o2 == o1+1, t17, arc)
	tk.last, lastRms = oneObs(o2, len(arc)-1, o2 == o1+1, t83, arc)
	return
}

func oneObs(o1, o2 int, arcsUseAllObs bool, pt float64, arc []obs.VObs) (result *obs.SiteObs, rms float64) {
	// a default result. (again, arc is guaranteed to be all ground based)
	result = arc[o1].(*obs.SiteObs)

	// case 1.  simple, only a single obs for the "arc".
	if o1 == o2 {
		return
	}

	// case 2.  two points for the arc.
	if o1 == o2-1 {
		// case 2.1.  fairly simple:  if there's stuff between the arcs
		// and the target percentile (17 or 83) is off of this arc,
		// just return the end of the arc.
		if !arcsUseAllObs {
			var dt float64
			end := result
			if o1 == 0 {
				end = arc[1].(*obs.SiteObs)
				dt = end.Mjd - pt
			} else {
				dt = pt - end.Mjd
			}
			if dt < 0 {
				return end, 0
			}
		}

		// case 2.2.  linearly interpolate along the great circle
		// connecting the points.
		r2 := *(arc[o2].(*obs.SiteObs))
		t := make([]float64, 2)
		s := make(coord.SphrS, 2)
		t[0] = result.Mjd
		s[0] = result.Sphr
		t[1] = r2.Mjd
		s[1] = r2.Sphr
		lmf := lmfit.New(t, s)
		if arcsUseAllObs {
			// gc midpoint
			r2.Mjd = (result.Mjd + r2.Mjd) * .5
		} else {
			// gc at pt
			r2.Mjd = pt
		}
		r2.Sphr = *lmf.Pos(r2.Mjd)
		return &r2, 0
	}

	// case 3.  3 or more points in arc.
	// synthesize obs at point along gc
	var tr float64
	if arcsUseAllObs {
		// median time of arc
		is := (o1 + o2) / 2
		tr = arc[is].(*obs.SiteObs).Mjd
		if is+is < o1+o2 {
			tr = (tr + arc[is+1].(*obs.SiteObs).Mjd) * .5
		}
	} else {
		var dt float64
		if o1 == 0 {
			result = arc[o2].(*obs.SiteObs)
			dt = result.Mjd - pt
		} else {
			dt = pt - result.Mjd
		}
		if dt < 0 {
			tr = result.Mjd
		} else {
			tr = pt
		}
	}

	// gc fit, result is computed obs at time tr
	np := o2 - o1 + 1
	t := make([]float64, np)
	s := make(coord.SphrS, np)
	for i, v := range arc[o1 : o2+1] {
		m := v.(*obs.SiteObs)
		t[i] = m.Mjd
		s[i] = m.Sphr
	}
	lmf := lmfit.New(t, s)
	r2 := *result
	r2.Mjd = tr
	r2.Sphr = *lmf.Pos(tr)
	return &r2, lmf.Rms() // return an rms for this arc
}
