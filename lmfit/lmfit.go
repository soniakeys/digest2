// Package lmfit fits observations to linear motion along a great circle.
package lmfit

import (
	"math"

	"coord"
)

const (
	twoPi     = 2 * math.Pi
	threePi   = 3 * math.Pi
	arcSecRad = math.Pi / (180 * 3600)
)

// LmFit represents the linear motion resulting from a great circle fit.
// It can be queried for various statistics on the fit and used to generate
// synthetic coordinates corresponding to the linear motion.
type LmFit struct {
	// various intermediate values, useful for generating statistics after
	// the fit is done
	mRot  coord.M3    // rotation matrix
	ra0   float64     // ra offset
	rs    coord.SphrS // rotated ra and dec
	t0    float64     // time offset
	nTime []float64   // normalized times
	// fit solution parameters
	r0, rr, d0, dr float64
}

// New does the great circle fit.
//
// Args:
//   mjd  -- mjd times of observations
//   s    -- ra, dec coordinates in radians
func New(mjd []float64, s coord.SphrS) *LmFit {
	nObs := len(mjd)
	if len(s) != nObs {
		return nil
	}
	var lmf LmFit

	// convert obs to cartesian
	var c coord.CartS
	c.FromSphrS(s)

	// vector normal to motion
	var n coord.Cart
	n.Cross(&c[0], &c[nObs-1])
	nmag2 := n.Square()
	nmag := math.Sqrt(nmag2)

	// rotation axis is right angle to n projected onto z=0 plane
	// and normalized to unit vector
	nxy := 1 / math.Sqrt(n.X*n.X+n.Y*n.Y)
	gcix := n.Y * nxy
	gciy := -n.X * nxy

	// rotation angle is angle from n to +z
	sina := math.Sqrt(nmag2-n.Z*n.Z) / nmag
	cosa := n.Z / nmag

	// build mRot, rotation matrix that will rotate the coordinate system
	// to the obs, so that a cylindrical projection will have negligible
	// distortion.)
	sinagx := sina * gcix
	sinagy := sina * gciy
	onemcosa := 1 - cosa
	onemcosagx := onemcosa * gcix
	onemcosagxgy := onemcosagx * gciy
	lmf.mRot[0] = cosa + onemcosagx*gcix
	lmf.mRot[1] = onemcosagxgy
	lmf.mRot[2] = sinagy
	lmf.mRot[3] = onemcosagxgy
	lmf.mRot[4] = cosa + onemcosa*gciy*gciy
	lmf.mRot[5] = -sinagx
	lmf.mRot[6] = -sinagy
	lmf.mRot[7] = sinagx
	lmf.mRot[8] = cosa

	// rotate all of cart
	var rotated coord.CartS
	rotated.Mult3S(&lmf.mRot, c)

	// transpose rotation array  so it will derotate, after least-squares fit.
	lmf.mRot.Transpose(&lmf.mRot)
	/*
		// invariant:  first and last points should be on the z=0 plane
		// of the rotated coordinate system now.
		const rTol = arcSecRad
		if math.Abs(rotated[0].Z) > rTol || math.Abs(rotated[nObs-1].Z) > rTol {
			panic(fmt.Sprintf("%f %f %f\n*** rotation failed ***",
				rotated[0].Z, rotated[nObs-1].Z, rTol))
		}
	*/
	// convert back to spherical coordinates for adjustment.
	// (this does the cylindrical projection.)
	lmf.rs.FromCartS(rotated)

	// normalize ra to near 0 to avoid wraparound problems
	ra0 := lmf.rs[0].Ra
	lmf.ra0 = ra0
	for i, r1 := range lmf.rs {
		lmf.rs[i].Ra = math.Mod(r1.Ra+threePi-ra0, twoPi) - math.Pi
	}

	// normalize time to near 0 to maintain precision
	t0 := mjd[0]
	lmf.t0 = t0
	lmf.nTime = make([]float64, nObs)
	for i, m1 := range mjd {
		lmf.nTime[i] = m1 - t0
	}

	if nObs == 2 {
		lmf.r0 = 0
		lmf.rr = lmf.rs[1].Ra / lmf.nTime[1]
		lmf.d0 = 0
		lmf.dr = 0
	} else {
		// here's the least squares stuff
		var sumT, sumRA, sumDec, sumT2, sumTRA, sumTDec float64
		for i, t1 := range lmf.nTime {
			sumT += t1
			sumRA += lmf.rs[i].Ra
			sumDec += lmf.rs[i].Dec
			sumT2 += t1 * t1
			sumTRA += t1 * lmf.rs[i].Ra
			sumTDec += t1 * lmf.rs[i].Dec
		}
		fn := float64(nObs)
		invD := 1 / (fn*sumT2 - sumT*sumT)
		lmf.r0 = invD * (sumRA*sumT2 - sumTRA*sumT)
		lmf.rr = invD * (fn*sumTRA - sumRA*sumT)
		lmf.d0 = invD * (sumDec*sumT2 - sumTDec*sumT)
		lmf.dr = invD * (fn*sumTDec - sumDec*sumT)
	}
	return &lmf
}

// Pos returns position at time t with fitted linear motion.
func (lmf *LmFit) Pos(t float64) (p *coord.Sphr) {
	nt := t - lmf.t0
	s := &coord.Sphr{
		Ra:  lmf.r0 + lmf.rr*nt + lmf.ra0,
		Dec: lmf.d0 + lmf.dr*nt,
	}
	var c coord.Cart
	return s.FromCart(c.Mult3(&lmf.mRot, c.FromSphr(s)))
}

// Res computes and returns residuals.
func (lmf *LmFit) Res() coord.SphrS {
	// residuals are computed on rotated-and-derotated observed to
	// minimize any systematic errors from rotating and derotating the
	// computed.

	// observed positions:  copy from rotated observed and fix up ra
	so := append(coord.SphrS{}, lmf.rs...)
	for i := range so {
		so[i].Ra += lmf.ra0
	}
	// then derotate back up to original place in the sky
	var c coord.CartS
	so.FromCartS(c.Mult3S(&lmf.mRot, c.FromSphrS(so)))

	// computed positions, derotated similarlly
	sc := make(coord.SphrS, len(so))
	for i, nt := range lmf.nTime {
		sc[i] = coord.Sphr{
			Ra:  lmf.r0 + lmf.rr*nt + lmf.ra0,
			Dec: lmf.d0 + lmf.dr*nt,
		}
	}
	sc.FromCartS(c.Mult3S(&lmf.mRot, c.FromSphrS(sc)))

	// repurpose so to hold residuals
	for i, so1 := range so {
		so[i] = coord.Sphr{
			Ra:  (so1.Ra - sc[i].Ra) * math.Cos(sc[i].Dec) / arcSecRad,
			Dec: (so1.Dec - sc[i].Dec) / arcSecRad,
		}
	}
	return so
}

// RmsRes returns rms of 2d residual and the component ra and dec residuals.
//
// Note:  The 2d residual is the 2d distance between computed and observed.
// By contrast, it is not the rms of the ra and dec residuals considered as
// separate values.  (The 2d distance that this function returns is smaller
// than the rms of separate values by a factor of sqrt(2).)
func (lmf *LmFit) RmsRes() (float64, coord.SphrS) {
	res := lmf.Res()
	var s float64
	for _, r1 := range res {
		s += r1.Ra*r1.Ra + r1.Dec*r1.Dec
	}
	return math.Sqrt(s / float64(len(res))), res
}

// Rms returns just the rms, never mind the residuals.
func (lmf *LmFit) Rms() float64 {
	rms, _ := lmf.RmsRes()
	return rms
}
