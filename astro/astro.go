// Package astro, stuff generally useful in astronomy.
package astro

import (
	"math"

	"digest2/coord"
)

const (
	K    = .01720209895
	InvK = 1 / K
	U    = K * K
)

var twoPi = 2 * math.Pi

// AeiHv, solves Keplarian elements from state vectors.
//
// Actually stretching the package claim of "generally useful," this
// function has parameters and return values most efficient for Digest2.
//
// Args:
//   p = position: sun object vector, in AU
//   v = object velocity vector, scaled by gravitational constant
//   d = sun-object distance pre computed from p
//
// Returns keplarian a, e, i in AU and radians, plus the momentum vector.
func AeiHv(p, v *coord.Cart, d float64) (a, e, i float64, hvp *coord.Cart) {

	// momentum vector
	var hv coord.Cart
	hv.Cross(p, v)
	hsq := hv.Square()
	hm := math.Sqrt(hsq)

	// solve for semi-major axis
	// (and the inverse--it comes in handy)
	vsq := v.Square()
	temp := 2 - d*vsq

	// for stability, require a < 100
	if d > temp*100 {
		return
	}
	a = d / temp
	inva := temp / d

	// solve for eccentricity
	// (stability test on a (above) should keep result real)
	e = math.Sqrt(1 - hsq*inva)

	// stability test:  require e < .99
	if e > .99 {
		return
	}

	// solve for inclination.

	// reliable check for i=0.  handles loss of precision in h computation.
	iZero := hv.Z >= hm
	// combination of stability tests on a and e (above) should
	// ensure that hm is well above zero.
	if !iZero {
		i = math.Acos(hv.Z/hm) * 180 / math.Pi
	}
	return a, e, i, &hv
}

// HMag computes H from V magnitude.
func HMag(oov, sov *coord.Cart, vmag, ood, sod float64) float64 {
	rdelta := ood * sod
	cospsi := oov.Dot(sov) / rdelta

	if cospsi < -.9999 {
		// object is straight into the sun.  doesn't seem too likely,
		// but anyway, this returns a valid value.
		return 30
	}

	tanhalf := math.Sqrt(1-cospsi*cospsi) / (1 + cospsi)
	phi1 := math.Exp(-3.33 * math.Pow(tanhalf, 0.63))
	phi2 := math.Exp(-1.87 * math.Pow(tanhalf, 1.22))
	return vmag -
		5.*math.Log10(rdelta) +
		2.5*math.Log10(.85*phi1+.15*phi2)
}

// Lst, computes local sidereal time.
func Lst(j0, longitude float64) float64 {
	t := (j0 - 15019.5) / 36525
	th := (6.6460656 + (2400.051262+0.00002581*t)*t) / 24
	ut := math.Mod(1, j0-.5)
	return math.Mod(th+ut+longitude, twoPi)
}

// Se2000 computes solar ephemeris, J2000.
//
// Returns:
//   sunEarth:  sun-earth vector in equatorial coordinates.
//   soe, coe:  sine and cosine of ecciptic.
//
// Notes:
//   Approximate solar coordinates, per USNO.  Originally from
//   http://aa.usno.navy.mil/faq/docs/SunApprox.html, page now at
//   http://www.usno.navy.mil/USNO/astronomical-applications/
//   astronomical-information-center/approx-solar.
func Se2000(mjd float64) (sunEarth coord.Cart, soe, coe float64) {
	// USNO algorithm is in degrees.  To mimimize confusion, work in
	// degrees here too, only converting to radians as needed for trig
	// functions.
	d := mjd - 51544.5
	g := 357.529 + .98560028*d // mean anomaly of sun, in degrees
	q := 280.459 + .98564736*d // mean longitude of sun, in degrees
	g2 := g + g
	sg, cg := math.Sincos(g*math.Pi/180) // send radians to trig function
	sg2, cg2 := math.Sincos(g2*math.Pi/180)

	// ecliptic longitude, in degrees still
	l := q + 1.915*sg + .020*sg2

	// distance in AU
	r := 1.00014 - .01671*cg - .00014*cg2

	// obliquity of ecliptic in degrees
	e := 23.439 - .00000036*d
	soe, coe = math.Sincos(e*math.Pi/180)

	// equatorial coordinates
	sl, cl := math.Sincos(l*math.Pi/180)
	sunEarth.X = r * cl
	rsl := r * sl
	sunEarth.Y = rsl * coe
	sunEarth.Z = rsl * soe
	return
}
