// Copyright 2010-2012 Sonia Keys
// Use of this source code is governed by an MIT-style
// license that can be found in the LICENSE file.

package d2solver

import "math"

// clipErr computes the obs err to use based on defaults and on rms computed
// from observations in the tracklet.
func (s *D2Solver) clipErr(computedRms float64, qual string) (clipped float64) {
	// look for config file specified obs err for this site
	defaultErr, ok := s.obsErrMap[qual]
	if !ok {
		// not there, fall back on default (which also may been specified
		// in the config file, or may be hard coded default.)
		defaultErr = s.obsErrDefault
	}
	if defaultErr == 0 {
		// if obs err is configured to be zero, that
		// takes precedence over any computed rms
		return 0
	} else if computedRms == 0 {
		// if no rms could be computed, use the default obs err.
		return defaultErr
	}
	// finally, consider rms, except it is in arc seconds. we need err
	// in radians
	computedErr := computedRms * math.Pi / (180 * 3600)
	// then just return the greater of the two
	if defaultErr > computedErr {
		return defaultErr
	}
	return computedErr
}
