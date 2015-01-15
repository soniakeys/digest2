/*
Command s3mbin generates a file, s3m.dat, for use by the program muk.

s3m.dat is distributed with the program muk, so you do not need to run
s3mbin at all.  The program is provided for those interested in generation
of s3m.dat.

s3mbin reads the PanSTARRS Synthetic Solar System Model (S3M) and reduces it
to a binned model.  While the PanSTARRS S3M authors have made the complete S3M
freely available on the internet, they have not licensed it for redistribution
in its raw form.  They have however, granted the digest2 authors permission to
distribute histograms, or binned models of the S3M.  This program generates
this binned model.

Usage

Usage:

   s3mbin [output file]
   s3mbin -v

The program looks in one of two places for the S3M files.  First, it checks
for an environment variable, S3M, which is set to a directory containing the
unzipped s3m files.  If the environment variable is not set, it looks for
a directory "s3m" in the current directory.

The output file, s3m.dat, by default is generated in the path "../muk" relative
to the s3mbin source directory.  Alternatively the output path or file name
can be specified as a command line argument.

-------------
Public domain.
*/
package main

import (
	"bufio"
	"flag"
	"fmt"
	"go/build"
	"io"
	"log"
	"os"
	"path/filepath"
	"runtime"
	"strconv"
	"strings"

	"digest2/d2bin"
	"github.com/soniakeys/exit"
)

const parentImport = "digest2"
const versionString = "s3mbin version 0.2"
const copyrightString = "Public domain."

// Orbits are binned in four dimensions of q, e, i, and H.
// The partitions in each dimension vary in size.

func init() {
	d2bin.QPart = []float64{.4, .7, .8, .9, 1, 1.1, 1.2, 1.3,
		1.4, 1.5, 1.67, 1.8, 2, 2.2, 2.4, 2.6, 2.8, 3,
		3.2, 3.5, 4, 4.5, 5, 5.5, 10, 20, 30, 40, 100}
	d2bin.EPart = []float64{.1, .2, .3, .4, .5, .7, .9, 1.1}
	d2bin.IPart = []float64{2, 5, 10, 15, 20, 25, 30, 40, 60, 90, 180}
	d2bin.HPart = []float64{6, 8, 10, 11, 12, 13, 14, 15,
		16, 17, 18, 19, 20, 21, 22, 23, 24, 25.5}
	d2bin.MSize = len(d2bin.QPart) * len(d2bin.EPart) *
		len(d2bin.IPart) * len(d2bin.HPart)
	d2bin.LastH = len(d2bin.HPart) - 1
}

func readPart(fCh chan string, mCh chan *d2bin.Model) {
	m := d2bin.New()
	for f := range fCh {
		binS3m(m, f, f != "S0")
	}
	mCh <- m
}

var nl bool
var nOrbits, nModel int
var nClass = make([]int, len(d2bin.CList))
var s3mPath string

func binS3m(m *d2bin.Model, fn string, clipNeo bool) {
	fn = filepath.Join(s3mPath, fn+".s3m")
	if nl {
		fmt.Println()
		nl = false
	}
	fmt.Println(fn)
	f, err := os.Open(fn)
	if err != nil {
		log.Println(err)
		return
	}
	defer f.Close()

	bf := bufio.NewReader(f)
	var line string
	for {
		if line, err = bf.ReadString('\n'); err != nil {
			log.Println(err)
			return
		}
		if !strings.HasPrefix(line, "!!") {
			break
		}
	}
	var q, e, i, h float64
loop:
	for {
		f := strings.Fields(line)
		if len(f) < 14 {
			fmt.Println("unexpected format:", len(f), "fields")
			return
		}
		if q, err = strconv.ParseFloat(f[2], 64); err != nil {
			break
		}
		if e, err = strconv.ParseFloat(f[3], 64); err != nil {
			break
		}
		if i, err = strconv.ParseFloat(f[4], 64); err != nil {
			break
		}
		if h, err = strconv.ParseFloat(f[8], 64); err != nil {
			break
		}
		if q <= 0 || e < 0 || e > 1.1 || i < 0 || i >= 180 {
			goto read // crazy data
		}
		nOrbits++
		if nOrbits%100000 == 0 {
			fmt.Print(".")
			nl = true
		}
		if clipNeo && q < 1.3 {
			goto read // bad data
		}
		if iq, ie, ii, ih, inModel := d2bin.Qeih(q, e, i, h); inModel {
			nModel++
			x := d2bin.Mx(iq, ie, ii, ih)
			m.SS[x]++
			for c, cs := range d2bin.CList {
				if cs.IsClass(q, e, i, h) {
					nClass[c]++
					m.Class[c][x]++
				}
			}
		}
	read:
		line, err = bf.ReadString('\n')
		switch err {
		case nil:
		case io.EOF:
			if len(line) == 0 {
				return // normal return
			}
		default:
			break loop
		}
	}
	log.Println(err)
}

var s3mFiles = []string{
	// MB
	"S1_00", "S1_01", "S1_02", "S1_03", "S1_04", "S1_05", "S1_06",
	"S1_07", "S1_08", "S1_09", "S1_10", "S1_11", "S1_12", "S1_13",
	"S0",  // NEO
	"St5", // Jupiter Trojan
	"SR",  // SPC
	"SJ",  // JFC
	"ST",  // TNO
	"SS",  // SDO
	//	"SL",   // LPC.  Don't include.
}

func main() {
	defer exit.Handler()
	// parent dir of s3mbin, muk, etc.
	parentDir := ""
	if pkg, err := build.Import(parentImport, "", build.FindOnly); err == nil {
		parentDir = pkg.Dir
	}
	flag.Usage = func() {
		os.Stderr.WriteString(`Usage:
   s3mbin [output file]
   s3mbin -v

For full documentation:
   godoc s3mbin
`)
	}
	vers := flag.Bool("v", false, "display version and copyright")
	flag.Parse()
	if *vers {
		fmt.Println(versionString)
		fmt.Println(copyrightString)
		os.Exit(0)
	}

	// determine output dir and file name
	var outDir, outFile string
	switch flag.NArg() {
	case 0: // default output path and file name
	case 1: // specified output path or file name
		outDir, outFile = filepath.Split(flag.Arg(0))
	default:
		flag.Usage()
		os.Exit(1)
	}
	if outDir == "" {
		outDir = filepath.Join(parentDir, "muk")
	}
	if outFile == "" {
		outFile = d2bin.Sfn
	}

	// determine s3m directory
	s3mPath = os.Getenv("S3M")
	if s3mPath == "" {
		s3mPath = filepath.Join(".", "s3m")
	}

	// a quick check that the s3mPath is there
	if _, err := os.Stat(s3mPath); err != nil {
		exit.Log(err)
	}

	// a source of file names
	fCh := make(chan string)
	go func() {
		for _, f := range s3mFiles {
			fCh <- f
		}
		close(fCh)
	}()

	// start a number of file readers in parallel.
	// each returns a data set on mCh
	mCh := make(chan *d2bin.Model)
	nProc := runtime.GOMAXPROCS(0)
	if nProc > len(s3mFiles) {
		nProc = len(s3mFiles)
	}
	for i := 0; i < nProc; i++ {
		go readPart(fCh, mCh)
	}

	// combine data sets from readers
	s3m := d2bin.New()
	for i := 0; i < nProc; i++ {
		mp := <-mCh
		for x, c := range mp.SS {
			s3m.SS[x] += c
		}
		for c, a := range mp.Class {
			ac := s3m.Class[c]
			for x, cc := range a {
				ac[x] += cc
			}
		}
	}

	// show status
	if nl {
		fmt.Println()
	}
	fmt.Println(nOrbits, "orbits")
	fmt.Println(nModel, "in model")
	for c, nc := range nClass {
		fmt.Printf("%8d %s\n", nc, d2bin.CList[c].Heading)
	}

	// write results
	f, err := os.Create(filepath.Join(outDir, outFile))
	if err != nil {
		exit.Log(err)
	}
	if _, err = f.WriteString("S3M binned"); err != nil {
		f.Close()
		exit.Log(err) // catch error on first write
	}

	f.WriteString("\nq")
	for _, part := range d2bin.QPart {
		fmt.Fprintf(f, " %g", part) // ignore errors in the middle
	}
	f.WriteString("\ne")
	for _, part := range d2bin.EPart {
		fmt.Fprintf(f, " %g", part)
	}
	f.WriteString("\ni")
	for _, part := range d2bin.IPart {
		fmt.Fprintf(f, " %g", part)
	}
	f.WriteString("\nh")
	for _, part := range d2bin.HPart {
		fmt.Fprintf(f, " %g", part)
	}
	f.WriteString("\n")

	for i := 0; i < len(s3m.SS); {
		for _ = range d2bin.HPart {
			fmt.Fprintf(f, "%g ", s3m.SS[i])
			i++
		}
		f.WriteString("\n")
	}

	for cx, class := range s3m.Class {
		fmt.Fprintf(f, "%s\n", d2bin.CList[cx].Heading)
		for i := 0; i < len(class); {
			for _ = range d2bin.HPart {
				fmt.Fprintf(f, "%g ", class[i])
				i++
			}
			_, err = f.WriteString("\n")
		}
	}
	if err != nil { // print error on last write
		f.Close()
		exit.Log(err)
	}
	if err := f.Close(); err != nil {
		exit.Log(err)
	}
}
