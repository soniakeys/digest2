// Copyright 2012 Sonia Keys
// Use of this source code is governed by an MIT-style
// license that can be found in the LICENSE file.

package mpc

import (
	"bufio"
	"errors"
	"fmt"
	"io"
	"math"
	"strconv"
	"strings"

	"code.google.com/p/digest2/go/coord"
	"code.google.com/p/digest2/go/obs"
)

// ParseObs80 parses a single line observation in the MPC 80 column format.
// This function handles only ground observatory based observations.
//
// Input line80 must be a string of 80 characters.  Other lengths are an error.
// The observatory code in columns 78-80 must exist in the input map.
func ParseObs80(line80 string, ocm obs.ParallaxMap) (desig string,
	o obs.VObs, err error) {
	if len(line80) != 80 {
		err = errors.New("ParseObs80 requires 80 characters")
		return
	}

	// the intent of reallocating desig (and later, obscode) is 
	// to allow line80 to be garbage collected sooner. no idea it really helps.
	desig = string([]byte(strings.TrimSpace(line80[:12])))

	d := line80[15:32]
	mjd := parseDate(d)
	if mjd == 0 {
		err = fmt.Errorf("ParseObs80: Invalid date (%s)", d)
		return
	}

	var rah, ram int
	var ras float64
	rah, err = strconv.Atoi(strings.TrimSpace(line80[32:34]))
	if err == nil {
		ram, err = strconv.Atoi(strings.TrimSpace(line80[35:37]))
		if err == nil {
			ras, err =
				strconv.ParseFloat(strings.TrimSpace(line80[38:44]), 64)
		}
	}
	if err != nil {
		err = fmt.Errorf("ParseObs80: Invalid RA (%s), %v", line80[32:44], err)
		return
	}

	decg := line80[44] // minus sign
	var decd, decm int
	var decs float64
	decd, err = strconv.Atoi(strings.TrimSpace(line80[45:47]))
	if err == nil {
		decm, err = strconv.Atoi(strings.TrimSpace(line80[48:50]))
		if err == nil {
			decs, err =
				strconv.ParseFloat(strings.TrimSpace(line80[51:56]), 64)
		}
	}
	if err != nil {
		err = fmt.Errorf("ParseObs80: Invalid Dec (%s), %v", line80[44:56], err)
		return
	}

	var mag float64
	if ts := strings.TrimSpace(line80[65:70]); len(ts) != 0 {
		mag, err = strconv.ParseFloat(ts, 64)
		if err != nil {
			err = fmt.Errorf("ParseObs80: Invalid mag (%s), %v", ts, err)
			return
		}
		band := line80[70]
		switch band {
		case 'V':
			break
		case 'B':
			mag -= .8
		default:
			mag += .4
		}
	}

	c := line80[77:80]
	par, ok := ocm[c]
	if !ok {
		return "", nil,
			fmt.Errorf("ParseObs80: Unknown observatory code (%s)", c)
	}

	obscode := string([]byte(line80[77:80]))

	if par == nil || line80[14] == 'S' {
		o = &obs.SatObs{Sat: obscode}
	} else {
		o = &obs.SiteObs{Par: par}
	}
	m := o.Meas()
	m.Mjd = mjd
	m.Ra = (float64(rah*60+ram)*60 + ras) * math.Pi / (12 * 3600)
	m.Dec = (float64(decd*60+decm)*60 + decs) * math.Pi / (180 * 3600)
	if decg == '-' {
		m.Dec = -m.Dec
	}
	m.Vmag = mag
	// could be enhanced to store program code, eg.  if so, see obsErr
	// code in digest2.readConfig and make appropriate changes.
	m.Qual = obscode
	return
}

var flookup = [13]int{0, 306, 337, 0, 31, 61, 92, 122, 153, 184, 214, 245, 275}

func parseDate(line80 string) float64 {
	year, err := strconv.Atoi(line80[:4])
	if err != nil {
		return 0
	}
	month, err := strconv.Atoi(line80[5:7])
	if err != nil {
		return 0
	}
	day, err := strconv.ParseFloat(strings.TrimSpace(line80[8:]), 64)
	if err != nil {
		return 0
	}
	z := year + (month-14)/12
	m := flookup[month] + 365*z + z/4 - z/100 + z/400 - 678882
	return float64(m) + day
}

func parseMpcSat2(line80, desig string, s1 *obs.SatObs) {
	if desig != strings.TrimSpace(line80[:12]) {
		return
	}
	if s1.Mjd != parseDate(line80) {
		return
	}
	if line80[77:80] != s1.Sat {
		return
	}

	x, ok := parseMpcOffset(line80[34:46])
	if !ok {
		return
	}
	y, ok := parseMpcOffset(line80[46:58])
	if !ok {
		return
	}
	z, ok := parseMpcOffset(line80[58:70])
	if !ok {
		return
	}
	if line80[32] == '1' {
		// Scale factor = 1 / 1 AU in km.
		const sf = 1 / 149.59787e6
		x *= sf
		y *= sf
		z *= sf
	}
	s1.Offset = coord.Cart{X: x, Y: y, Z: z}
}

func parseMpcOffset(off string) (float64, bool) {
	v, err := strconv.ParseFloat(strings.TrimSpace(off[1:]), 64)
	switch {
	case err != nil:
		break
	case off[0] == '-':
		return -v, true
	case off[0] == '+' || off[0] == ' ':
		return v, true
	}
	return 0, false
}

// SplitTracklets splits an observation stream up into tracklets.
//
// The function does no sorting.  The stream, iObs, must have observations
// grouped by object and sorted chronologically within each object.
// (The logic is a little bit mousetrap.)
//
// Valid tracklets are parsed against ocdMap and retuned on channel tkCh.
// Read errors are relayed on errCh should be considered fatal.
// Parse errors are not fatal.  They are quietly ignored.
// Lines causing parse errors and lines not forming valid tracklets are
// unceremoniously dropped.
func SplitTracklets(iObs io.Reader, ocdMap obs.ParallaxMap,
	tkCh chan *obs.Tracklet, errCh chan error) {
	bf := bufio.NewReader(iObs)
	var des0 string
	var o obs.VObs
	var desig string
	obuf := make([]obs.VObs, 0, 4)
	for {
		bLine, pre, err := bf.ReadLine()
		if err == io.EOF {
			break
		}
		if err != nil {
			errCh <- err
			break
		}
		if pre {
			errCh <- errors.New("splitTracklets: unexpected long line")
			break
		}
		if len(bLine) != 80 {
			continue
		}
		line := string(bLine)
		if line[14] == 's' {
			if s, ok := o.(*obs.SatObs); ok {
				parseMpcSat2(line, desig, s)
			}
			continue
		}
		desig, o, err = ParseObs80(line, ocdMap)
		switch {
		case err != nil:
			sendValid(des0, obuf, tkCh)
			obuf = obuf[:0]
		default:
			sendValid(des0, obuf, tkCh)
			fallthrough
		case len(obuf) == 0:
			des0 = desig
			obuf = obuf[:1]
			obuf[0] = o
		case desig == des0:
			obuf = append(obuf, o)
		}
	}
	sendValid(des0, obuf, tkCh)
	close(tkCh)
}

// checks that observations make a valid tracklet,
// allocates and sends the tracklet.
func sendValid(desig string, obuf []obs.VObs, tkCh chan *obs.Tracklet) {
	if len(obuf) < 2 {
		return
	}
	// the first observation time must be positive and
	// observation times must increase after that
	var t0 float64
	for i := range obuf {
		t := obuf[i].Meas().Mjd
		if t <= t0 {
			return
		}
		t0 = t
	}
	// object must show motion over the tracklet
	first := obuf[0].Meas()
	last := obuf[len(obuf)-1].Meas()
	if first.Ra == last.Ra && first.Dec == last.Dec {
		return
	}
	tkCh <- &obs.Tracklet{
		Desig: desig,
		Obs:   append([]obs.VObs{}, obuf...),
	}
}
