package main

import (
	"bufio"
	"compress/gzip"
	"encoding/gob"
	"flag"
	"fmt"
	"go/build"
	"io"
	"log"
	"math"
	"os"
	"os/exec"
	"path/filepath"
	"strconv"
	"strings"

	"d2bin"
)

const parentImport = "digest2"
const versionString = "muk version 0.1 Go source."
const copyrightString = "Public domain, Smithsonian Astrophysical Observatory."
const aofn = "astorb.dat"

type fatal struct {
	err error
}

func exit(err error) {
	panic(fatal{err})
}

func handleFatal() {
	if err := recover(); err != nil {
		if f, ok := err.(fatal); ok {
			log.Fatal(f.err)
		}
		panic(err)
	}
}

func main() {
	defer handleFatal()

	// muk package dir.  required location of s3m.dat, default location
	// of astorb.dat
	mukDir := ""
	if pkg, err := build.Import(parentImport+"/muk", "",
		build.FindOnly); err == nil {
		mukDir = pkg.Dir
	}
	// parent dir of muk and digest2.  location for LICENSE, also location
	// for output file digest2.gmodel (d2bin.Mfn)
	parentDir := ""
	if pkg, err := build.Import(parentImport, "", build.FindOnly); err == nil {
		parentDir = pkg.Dir
	}
	// default location for astorb.dat
	defPath := filepath.Join(mukDir, aofn)

	flag.Usage = func() {
		os.Stderr.WriteString(`Usage:
  muk                                  Use default location for astorb.dat.
  muk -v                               Display version and copyright.
  muk -a=<muk source path>/astorb.dat  Specify astorb.dat path or file name.

Default:
  -a=` + defPath + `

For full documentation:
   go doc ` + parentImport + `/muk
`)
	}
	clPath := flag.String("a", defPath, "astorb.dat path or file name")
	vers := flag.Bool("v", false, "display version and copyright")
	flag.Parse()
	if *vers {
		fmt.Println(versionString)
		fmt.Println(copyrightString)
		os.Exit(0)
	}
	if flag.NArg() > 0 {
		flag.Usage()
		os.Exit(1)
	}
	astorbPath := *clPath
	switch {
	case astorbPath != defPath:
		// user specified a path or file name.  see if it needs to be fixed up.
		clDir, clFile := filepath.Split(*clPath)
		if clDir == "" {
			// add default directory
			astorbPath = filepath.Join(mukDir, clFile)
			break
		}
		fi, statErr := os.Stat(*clPath)
		if statErr != nil {
			exit(statErr)
		}
		if fi.IsDir() {
			// add default file name
			astorbPath = filepath.Join(*clPath, aofn)
		}
	default:
		// user took default.  we're happy if it stats, ...
		if _, err := os.Stat(defPath); err == nil {
			break
		}
		// otherwise try wget.  (Go code would be nice here but existing
		// packages I found didn't have progress features, which are
		// important because this file is big and the site is slow.)
		fmt.Printf(`
%s not found.
Accessing ftp://ftp.lowell.edu/pub/elgb/astorb.dat.gz...
This process is often time consuming.

`, defPath)
		c := exec.Command("wget",
			"ftp://ftp.lowell.edu/pub/elgb/astorb.dat.gz",
			"-O", "-")
		c.Stderr = os.Stderr
		wOut, err := c.StdoutPipe()
		if err != nil {
			exit(err)
		}
		if err = c.Start(); err != nil {
			exit(err)
		}
		// gunzip.  This we can handle.
		uOut, err := gzip.NewReader(wOut)
		if err != nil {
			exit(err)
		}
		f, err := os.Create(defPath)
		if err != nil {
			exit(err)
		}
		if _, err = io.Copy(f, uOut); err != nil {
			exit(err)
		}
		f.Close()
	}

	// S3M file required to be in muk directory.
	sPath := filepath.Join(mukDir, d2bin.Sfn)

	// in a small conflation of features, status messages show only file names
	// if the default is used for astorb.dat, and full paths otherwise
	if astorbPath == defPath {
		fmt.Println("Reading", d2bin.Sfn)
	} else {
		fmt.Println("Reading", sPath)
	}

	f, err := os.Open(sPath)
	if err != nil {
		exit(err)
	}
	bf := bufio.NewReader(f)
	var ln int
	// close on f, ln
	corrupt := func(i interface{}) {
		if i != nil {
			log.Println(i)
		}
		f.Close()
		exit(fmt.Errorf("%s corrupt. line %d", d2bin.Sfn, ln))
	}
	mustRead := func() string {
		ln++
		line, isPre, err := bf.ReadLine()
		if err != nil {
			corrupt(err)
		}
		if isPre {
			corrupt("unexpected long line")
		}
		return string(line)
	}
	if mustRead() != "S3M binned" {
		corrupt(`"S3M binned" expected`)
	}
	readPart := func(ele byte) []float64 {
		line := mustRead()
		if len(line) < 1 || line[0] != ele {
			corrupt(fmt.Sprintf("%c line expected", ele))
		}
		flds := strings.Fields(line[1:])
		part := make([]float64, len(flds))
		for px, p := range flds {
			f, err := strconv.ParseFloat(p, 64)
			if err != nil {
				corrupt(err)
			}
			part[px] = f
		}
		return part
	}
	d2bin.QPart = readPart('q')
	d2bin.EPart = readPart('e')
	d2bin.IPart = readPart('i')
	d2bin.HPart = readPart('h')
	d2bin.LastH = len(d2bin.HPart) - 1
	d2bin.MSize = len(d2bin.QPart) * len(d2bin.EPart) * len(d2bin.IPart) * len(d2bin.HPart)
	readBins := func() []float64 {
		bins := make([]float64, d2bin.MSize)
		for bx := 0; bx < d2bin.MSize; {
			flds := strings.Fields(mustRead())
			if len(flds) != len(d2bin.HPart) {
				corrupt(len(flds))
			}
			for _, s := range flds {
				p, err := strconv.ParseFloat(s, 64)
				if err != nil {
					corrupt(err)
				}
				bins[bx] = p
				bx++
			}
		}
		return bins
	}
	var s3m d2bin.Model
	s3m.SS = readBins()
	s3m.Class = make([][]float64, len(d2bin.CList))
	for cx, class := range d2bin.CList {
		if mustRead() != class.Heading {
			corrupt(fmt.Sprintf(`class "%s" expected`, class.Heading))
		}
		s3m.Class[cx] = readBins()
	}
	f.Close()

	known := d2bin.New()
	_, aoFile := filepath.Split(astorbPath)
	if astorbPath == defPath {
		fmt.Printf("Reading %s...\n", aoFile)
	} else {
		fmt.Printf("Reading %s...\n", astorbPath)
	}

	// Note on file size:  astorb.dat is over 100M.  I found that the
	// following bufio code ran about twice as fast as equivalent code
	// using ioutil.Readfile.  I usually like bufio.ReadLine, but that
	// seems to offer a big advantage only when you can work with bytes.
	// Here we need strconv functions, so bufio.ReadString seems best.
	//
	// Note also that astorb.data is ASCII encoded.
	forb, err := os.Open(astorbPath)
	if err != nil {
		exit(err)
	}
	defer forb.Close()
	bfile := bufio.NewReaderSize(forb, 1<<10)
	line, err := bfile.ReadString('\n')
	if err != nil {
		exit(err)
	}
	var decpeuy_fails, decpeu_rejects, parsefails, outofmodel, lines, good int
	for ; err == nil; line, err = bfile.ReadString('\n') {
		lines += 1
		decpeuy, err := strconv.Atoi(line[242:246])
		if err != nil || decpeuy < 2000 {
			decpeuy_fails++
			continue
		}
		word := line[234:237] + "e" + line[238:241]
		decpeu, err := strconv.ParseFloat(word, 64)
		if err != nil {
			parsefails++
			continue
		}
		if decpeu > 60. {
			decpeu_rejects++
			continue
		}
		a, err := strconv.ParseFloat(strings.TrimSpace(line[169:181]), 64)
		if err != nil {
			parsefails++
			continue
		}
		e, err := strconv.ParseFloat(line[158:168], 64)
		if err != nil {
			parsefails++
			continue
		}
		i, err := strconv.ParseFloat(strings.TrimSpace(line[147:157]), 64)
		if err != nil {
			parsefails++
			continue
		}
		h, err := strconv.ParseFloat(strings.TrimSpace(line[42:47]), 64)
		if err != nil {
			parsefails++
			continue
		}
		q := a * (1 - e)
		iq, ie, ii, ih, inModel := d2bin.Qeih(q, e, i, h)
		if !inModel {
			outofmodel++
			continue
		}

		good++
		bx := d2bin.Mx(iq, ie, ii, ih)
		known.SS[bx]++
		for c, cs := range d2bin.CList {
			if cs.IsClass(q, e, i, h) {
				known.Class[c][bx]++
			}
		}
	}

	fmt.Println(lines, "lines in", aoFile)
	if parsefails > 0 {
		fmt.Println(parsefails, "lines failed to parse")
	}
	fmt.Println(decpeuy_fails, "lines had invalid date of peak ephemeris uncertainty")
	fmt.Println(decpeu_rejects, "oribits had excessive peak ephemeris uncertainty")
	if outofmodel > 0 {
		fmt.Println(outofmodel, "orbits out of model")
	}
	fmt.Println(good, "orbits usable")

	// from s3m and known, produce all=max(s3m, known)/sqrt(v)
	// and unk=(all-known)/sqrt(v), where v is the "volume" of the d2bin.
	all := d2bin.New()
	unk := d2bin.New()
	q0 := 0.
	x := 0
	for _, q1 := range d2bin.QPart {
		dq := q1 - q0
		q0 = q1
		e0 := 0.
		d1 := 1.
		for _, e1 := range d2bin.EPart {
			d0 := d1
			d1 = 1 - e1
			if d1 < 0 {
				d1 = 0
			}
			dae := dq * (e1 - e0) / (d0 + d1)
			e0 = e1
			i0 := 0.
			for _, i1 := range d2bin.IPart {
				daei := dae * (i1 - i0)
				i0 = i1
				h0 := 0.
				for _, h1 := range d2bin.HPart {
					isqv := 1 / math.Sqrt(daei*(h1-h0))
					h0 = h1
					if known.SS[x] > s3m.SS[x] {
						all.SS[x] = float64(known.SS[x]) * isqv
						unk.SS[x] = 0
					} else {
						all.SS[x] = float64(s3m.SS[x]) * isqv
						unk.SS[x] = float64(s3m.SS[x]-known.SS[x]) * isqv
					}
					for c, kc := range known.Class {
						sc := s3m.Class[c]
						if kc[x] > sc[x] {
							all.Class[c][x] = float64(kc[x]) * isqv
							unk.Class[c][x] = 0
						} else {
							all.Class[c][x] = float64(sc[x]) * isqv
							unk.Class[c][x] = float64(sc[x]-kc[x]) * isqv
						}
					}
					x++
				}
			}
		}
	}
	mPath := filepath.Join(parentDir, d2bin.Mfn)
	fbin, err := os.Create(mPath)
	if err != nil {
		exit(err)
	}
	if astorbPath == defPath {
		fmt.Println("Writing", d2bin.Mfn)
	} else {
		fmt.Println("Writing", mPath)
	}
	defer fbin.Close()
	enc := gob.NewEncoder(fbin)
	mustEncode := func(i interface{}) {
		if err := enc.Encode(i); err != nil {
			exit(err)
		}
	}
	mustEncode(d2bin.QPart)
	mustEncode(d2bin.EPart)
	mustEncode(d2bin.IPart)
	mustEncode(d2bin.HPart)
	mustEncode(d2bin.MSize)
	mustEncode(d2bin.LastH)
	mustEncode(all)
	mustEncode(unk)
}
