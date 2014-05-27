package mpc_test

import (
	"io/ioutil"
	"math"
	"os"
	"testing"

	"code.google.com/p/digest2/go/mpc"
)

var siteTestCases = []struct {
	code          string
	lon, cos, sin float64
}{
	{"000", 0, .62411, .77873},
	{"248", 0, 0, 0},
	{"250", 0, 0, 0},
	{"644", 243.14022, .836325, .546877},
	{"E12", 149.0642, .85563, -.51621},
}

func TestOcd(t *testing.T) {
	fn := testFetch(t)
	testRead(t, fn)
}

func testFetch(t *testing.T) string {
	f, err := ioutil.TempFile("", "digest2ocd")
	if err != nil {
		t.Fatal(err)
	}
	fn := f.Name()
	f.Close()
	if err = mpc.FetchOcd(fn); err != nil {
		t.Fatal(err)
	}
	return fn
}

func testRead(t *testing.T, fn string) {
	defer os.Remove(fn)
	ocd, err := mpc.ReadOcd(fn)
	if err != nil {
		t.Fatal(err)
	}
	for _, c := range siteTestCases {
		switch s, ok := ocd[c.code]; {
		case !ok:
			t.Fatal("missing", c.code)
		case s == nil:
			if c.cos != 0 || c.sin != 0 {
				t.Fatal("nil stored for code", c.code)
			}
		case c.cos == 0 && c.sin == 0:
			t.Fatal("expected nil for code", c.code)
		case math.Abs(s.Longitude*360-c.lon) > 1e-10:
			t.Fatal("bad longitude, code", c.code)
		case math.Abs(s.RhoCosPhi*149.59787e9/6.37814e6-c.cos) > 1e-10:
			t.Fatal("bad rho cos, code", c.code)
		case math.Abs(s.RhoSinPhi*149.59787e9/6.37814e6-c.sin) > 1e-10:
			t.Fatal("bad rho sin, code", c.code)
		}
	}
}
