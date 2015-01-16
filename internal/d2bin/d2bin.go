// Public domain.

// Package d2bin defines the bin model used by digest2 and two utility programs
// muk and s3mbin.
package d2bin

import (
	"encoding/gob"
	"math"
	"os"
	"time"
)

// Sfn, the binned model
const Sfn = "s3m.dat"

// Mfn, the binned model combined with known population.
const Mfn = "digest2.gmodel"

// Model holds model population counts for all modeled objects in the
// solar system, and counts broken down by orbit class.  The slices of
// float64 are a flat representation of the 4-element model space.
type Model struct {
	SS    []float64
	Class [][]float64
}

// New allocates and initializes a model object.
func New() *Model {
	var m Model
	m.SS = make([]float64, MSize)
	m.Class = make([][]float64, len(CList))
	for c := range CList {
		m.Class[c] = make([]float64, MSize)
	}
	return &m
}

// ReadFile reads a population model.
//
// Argument fn is the filename of the model file created by muk.
//
// The model is returned in all and unk, also package variables QPart, EPart,
// IPart, HPart, MSize, and LastH are set.
func ReadFile(fn string) (all, unk Model, aoDate time.Time, aoLines int, err error) {
	var f *os.File
	f, err = os.Open(fn)
	if err != nil {
		return
	}
	defer f.Close()
	dec := gob.NewDecoder(f)
	if err = dec.Decode(&aoDate); err != nil {
		return
	}
	if err = dec.Decode(&aoLines); err != nil {
		return
	}
	if err = dec.Decode(&QPart); err != nil {
		return
	}
	dec.Decode(&EPart)
	dec.Decode(&IPart)
	dec.Decode(&HPart)
	dec.Decode(&MSize)
	dec.Decode(&LastH)
	dec.Decode(&all)
	err = dec.Decode(&unk)
	return
}

// Package variables that define the shape and size of the model.  They are
// constant after being set (in s3mbin) or loaded (in muk and digest2.)
var (
	QPart, EPart, IPart, HPart []float64
	MSize                      int
	LastH                      int
)

// Mx computes an index into the flat representation of a model.
func Mx(iq, ie, ii, ih int) int {
	return ((iq*len(EPart)+ie)*len(IPart)+ii)*len(HPart) + ih
}

// Qeih takes four real-valued elements and returns their bin indexes.
func Qeih(q, e, i, h float64) (qx, ex, ix, hx int, inModel bool) {
	if qx, ex, ix, inModel = Qei(q, e, i); inModel {
		hx = H(h)
	}
	return
}

// Qei takes three real-valued elements and returns their bin indexes.
func Qei(q, e, i float64) (qx, ex, ix int, inModel bool) {
	for q >= QPart[qx] {
		qx++
		if qx == len(QPart) {
			return
		}
	}
	for e >= EPart[ex] {
		ex++
		if ex == len(EPart) {
			return
		}
	}
	for i >= IPart[ix] {
		ix++
		if ix == len(IPart) {
			return
		}
	}
	return qx, ex, ix, true
}

// H takes a real-valued H magnitude and returns the corresponding bin index.
func H(h float64) (ih int) {
	for ; h >= HPart[ih] && ih < LastH; ih++ {
	}
	return
}

// Clist represents the modeled orbit classes
var CList = []struct {
	Abbr, Heading string
	IsClass       func(q, e, i, h float64) bool
}{
	{"Int", "MPC interest.", isMpcint},
	{"NEO", "NEO(q < 1.3)", isNeo},
	{"N22", "NEO(H <= 22)", isCMO},
	{"N18", "NEO(H <= 18)", isH18Neo},
	{"MC", "Mars Crosser", isMarsCrosser},
	{"Hun", "Hungaria gr.", isHungaria},
	{"Pho", "Phocaea group", isPhocaea},
	{"MB1", "Inner MB", isInnerMB},
	{"Pal", "Pallas group", isPallas},
	{"Han", "Hansa group", isHansa},
	{"MB2", "Middle MB", isMidMB},
	{"MB3", "Outer MB", isOuterMB},
	{"Hil", "Hilda group", isHilda},
	{"JTr", "Jupiter tr.", isTrojan},
	{"JFC", "Jupiter Comet", isJFC},
}

// 'MPC interesting' objects
// The definition of MPC interesting implemented here is:
// any of: q < 1.3, e > .5, i > 40, or Q > 10
func isMpcint(q, e, i, h float64) bool {
	return q < 1.3 || e >= .5 || i >= 40 || q*(1+e)/(1-e) > 10
}

// 'NEO' objects
// The definition of NEO implemented here is q < 1.3
func isNeo(q, e, i, h float64) bool {
	return q < 1.3
}

// H18 NEOs
// H rounded to nearest integer <= 18
func isH18Neo(q, e, i, h float64) bool {
	return q < 1.3 && h < 18.5
}

// H22 NEOs
// H rounded to nearest integer <= 22
func isCMO(q, e, i, h float64) bool {
	return q < 1.3 && h < 22.5
}

// Mars Crosser
// 1.3 <= q < 1.67, Q > 1.58
func isMarsCrosser(q, e, i, h float64) bool {
	return q < 1.67 && q >= 1.3 && q*(1+e)/(1-e) > 1.58
}

// Hungarias
// 1.78>a>2.0, e<.18, 16 < i < 28
// (a node test would be nice...)
func isHungaria(q, e, i, h float64) bool {
	if e > .18 || i < 16 || i > 34 {
		return false
	}
	a := q / (1 - e)
	return a < 2 && a > 1.78
}

// Phocaeas
// q>1.5, 2.2<a<2.45, 20<i<27
func isPhocaea(q, e, i, h float64) bool {
	if q < 1.5 || i < 20 || i > 27 {
		return false
	}
	a := q / (1 - e)
	return a < 2.45 && a > 2.2
}

// Inner Main Belt
// q>1.67, 2.1<a<2.5, i<7 at inner edge, <17 at outer
func isInnerMB(q, e, i, h float64) bool {
	if q < 1.67 {
		return false
	}
	a := q / (1 - e)
	return a < 2.5 && a > 2.1 && i < ((a-2.1)/.4)*10+7
}

// Hansas
// 2.55<a<2.72 e<.25, 20<i<23.5
func isHansa(q, e, i, h float64) bool {
	if e > .25 || i < 20 || i > 23.5 {
		return false
	}
	a := q / (1 - e)
	return a < 2.72 && a > 2.55
}

// Pallas group
// 2.5<a<2.8, e<.35, 24<i<37
func isPallas(q, e, i, h float64) bool {
	if e > .35 || i < 24 || i > 37 {
		return false
	}
	a := q / (1 - e)
	return a < 2.8 && a > 2.5
}

// Mid Main Belt
// 2.5<a<2.8, e<.45, i<20
func isMidMB(q, e, i, h float64) bool {
	if e > .45 || i > 20 {
		return false
	}
	a := q / (1 - e)
	return a < 2.8 && a > 2.5
}

// Outer Main Belt
//  2.8<a<3.25 e<.4, i < 20 inner edge, i < 36 outer edge
func isOuterMB(q, e, i, h float64) bool {
	if e > .4 {
		return false
	}
	a := q / (1 - e)
	return a > 2.8 && a < 3.25 && i < ((a-2.8)/.45)*16+20
}

// Hildas
// 3.9<a<4.02, e<.4, i<18
func isHilda(q, e, i, h float64) bool {
	if i > 18 || e > .4 {
		return false
	}
	a := q / (1 - e)
	return a > 3.9 && a < 4.02
}

// Trojans
// 5.05<a<5.35, e<.22, i<38
func isTrojan(q, e, i, h float64) bool {
	if e > .22 || i > 38 {
		return false
	}
	a := q / (1 - e)
	return a > 5.05 && a < 5.35
}

// Jupiter Family Comets
// 2 < Tj < 3, q >= 1.67
func isJFC(q, e, i, h float64) bool {
	if q < 1.3 {
		return false
	}
	tj := 5.2*(1-e)/q + 2*math.Sqrt(q*(1+e)/5.2)*math.Cos(i*math.Pi/180)
	return tj < 3 && tj > 2
}
