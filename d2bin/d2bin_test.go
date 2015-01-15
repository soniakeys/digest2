// Public domain.

package d2bin_test

import (
	"fmt"
	"testing"

	"digest2/d2bin"
)

func ExampleNew() {
	d2bin.MSize = 3
	fmt.Printf("%+v\n", d2bin.New())
	// Output:
	// &{SS:[0 0 0] Class:[[0 0 0] [0 0 0] [0 0 0] [0 0 0] [0 0 0] [0 0 0] [0 0 0] [0 0 0] [0 0 0] [0 0 0] [0 0 0] [0 0 0] [0 0 0] [0 0 0] [0 0 0]]}
}

func TestModel(t *testing.T) {
	all, unk, _, _, err := d2bin.ReadFile("../digest2.gmodel")
	if err != nil {
		t.Skip(err)
	}
	// echo partitions
	t.Log("QPart:", d2bin.QPart)
	t.Log("EPart:", d2bin.EPart)
	t.Log("IPart:", d2bin.IPart)
	t.Log("HPart:", d2bin.HPart)
	t.Log("MSize:", d2bin.MSize,
		len(d2bin.QPart)*len(d2bin.EPart)*
			len(d2bin.IPart)*len(d2bin.HPart))
	t.Log("LastH:", d2bin.LastH, len(d2bin.HPart)-1)

	// class list example
	e, i, h := .1, 20., 18.
	q := 1.8 * (1 - e)
	t.Logf("%-14s: %t\n",
		d2bin.CList[5].Heading, d2bin.CList[5].IsClass(q, e, i, h))
	t.Logf("%-14s: %t\n",
		d2bin.CList[6].Heading, d2bin.CList[6].IsClass(q, e, i, h))

	// bin indexes
	for _, h := range []float64{5.9, 6, 18, 25, 26} {
		t.Logf("H %4.1f = bin %d\n", h, d2bin.H(h))
	}
	qx, ex, ix, inModel := d2bin.Qei(1.6, .1, 18)
	t.Log("Q, E, I = 1.6, .1, 18 : bins", qx, ex, ix, inModel)

	qx, ex, ix, inModel = d2bin.Qei(100, .1, 18)
	t.Log("Q, E, I = 100, .1, 18 : bins", qx, ex, ix, inModel)

	qx, ex, ix, inModel = d2bin.Qei(1.6, 1.1, 18)
	t.Log("Q, E, I = 1.6, 1.1, 18 : bins", qx, ex, ix, inModel)

	qx, ex, ix, inModel = d2bin.Qei(1.6, .1, 180)
	t.Log("Q, E, I = 1.6, .1, 180 : bins", qx, ex, ix, inModel)

	qx, ex, ix, inModel = d2bin.Qei(0, 0, 0)
	t.Log("Q, E, I = 0, 0, 0 : bins", qx, ex, ix, inModel)

	t.Log("Mx(10, 1, 4, 11) =", d2bin.Mx(10, 1, 4, 11))

	// bin populations
	t.Log("Modeled population in bin: ", all.SS[16121])
	t.Log("Unknown population in bin: ", unk.SS[16121])
	t.Log("Modeled Hungarias in bin: ", all.Class[5][16121])
	t.Log("Unknown Hungarias in bin: ", unk.Class[5][16121])
	// in 2014,
	// Modeled population in bin:  402.49223594996226
	// Unknown population in bin:  147.58048651498615
	// Modeled Hungarias in bin:  380.1315561749643
	// Unknown Hungarias in bin:  143.10835055998658
}
