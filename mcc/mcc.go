package main

import (
	"flag"
	"fmt"
	"io/ioutil"
	"log"
	"math"
	"os"
	"strconv"
	"strings"
)

const parentImport = "digest2"
const versionString = "mcc version 0.1"
const copyrightString = "Public domain, Smithsonian Astrophysical Observatory."

var col int
var ignored int

func main() {
	// parse command line
	flag.Usage = func() {
		os.Stderr.WriteString(
			"Usage: mcc [options] <in-class> <out-of-class> [threshold]\n")
		flag.PrintDefaults()
		os.Stderr.WriteString(`
For full documentation:
   go doc ` + parentImport + `/mcc
`)
	}
	flag.IntVar(&col, "c", 1, "column containing class score")
	vers := flag.Bool("v", false, "display version and copyright")
	flag.Parse()
	if *vers {
		fmt.Println(versionString)
		fmt.Println(copyrightString)
		os.Exit(0)
	}
	if n := flag.NArg(); n < 2 || n > 3 {
		flag.Usage()
		os.Exit(1)
	}
	// parse threshold
	threshold := 50.
	thresholdPrec := 0
	if flag.NArg() == 3 {
		tStr := flag.Arg(2)
		var err error
		threshold, err = strconv.ParseFloat(tStr, 64)
		if err != nil {
			log.Fatalln("Bad threshold:", err)
		}
		if p := strings.Index(tStr, "."); p >= 0 {
			thresholdPrec = len(tStr) - p - 1
		}
	}
	// read in-class file (arg 1)
	tp, fn, err := aboveThreshold(flag.Arg(0), threshold)
	if err != nil {
		log.Fatalln("in-class file:", err)
	}
	// read out-of-class file (arg 2)
	fp, tn, err := aboveThreshold(flag.Arg(1), threshold)
	if err != nil {
		log.Fatalln("out-of-class file:", err)
	}
	// compute mcc
	tpf := float64(tp)
	fnf := float64(fn)
	fpf := float64(fp)
	tnf := float64(tn)
	mcc := 0.
	if d := (tpf + fpf) * (tpf + fnf) * (tnf + fpf) * (tnf + fnf); d > 0 {
		mcc = (tpf*tnf - fpf*fnf) / math.Sqrt(d)
	}
	// report statistics
	fmt.Println("\nIn-class file:     ", flag.Arg(0))
	fmt.Println("Out-of-class file: ", flag.Arg(1))
	fmt.Println("Total objects:     ", tp+fn+fp+tn)
	if ignored != 0 {
		fmt.Println("Lines ignored:     ", ignored)
	}
	fmt.Printf("Threshold:          %.*f\n", thresholdPrec, threshold)
	fmt.Println()
	fmt.Println("                       digest2 prediction")
	fmt.Println("                    -----------------------")
	fmt.Println("                     in-class  out-of-class")
	fmt.Printf("Actual in-class       %7d       %7d\n", tp, fn)
	fmt.Printf("Actual out-of-class   %7d       %7d\n", fp, tn)
	fmt.Println()
	fmt.Printf("Matthews correlation coefficient: %.2f\n", mcc)
}

func aboveThreshold(fn string, threshold float64) (ge, lt int, err error) {
	var b []byte
	b, err = ioutil.ReadFile(fn)
	if err != nil {
		return
	}
	for _, line := range strings.Split(string(b), "\n") {
		f := strings.Fields(line)
		if len(f) <= col {
			ignored++
			continue
		}
		score, err := strconv.ParseFloat(f[col], 64)
		if err != nil {
			ignored++
			continue
		}
		if score >= threshold {
			ge++
		} else {
			lt++
		}
	}
	return
}
