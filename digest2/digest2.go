package main

import (
	"bufio"
	"encoding/gob"
	"flag"
	"fmt"
	"go/build"
	"io"
	"log"
	"math"
	"math/rand"
	"os"
	"path/filepath"
	"regexp"
	"runtime"
	"strconv"
	"strings"
	"time"

	"digest2/bin"
	"digest2/mpc"
	"digest2/obs"
	"digest2/solver"
)

const parentImport = "digest2"
const versionString = "digest2 version 0.16 Go source."
const copyrightString = "Public domain, Smithsonian Astrophysical Observatory."

type fatal struct {
	err interface{}
}

// exit ultimately calls log.Fatal but first allows
// defered functions to be called.
func exit(err interface{}) {
	panic(fatal{err})
}

func main() {
	// fatal error handler
	defer func() {
		if err := recover(); err != nil {
			if f, ok := err.(fatal); ok {
				log.Fatal(f.err)
			}
			panic(err)
		}
	}()

	// these functions all set up package vars and terminate on error
	cl := parseCommandLine()
	ocdMap := readOcd(cl)
	all, unk := readModel(cl)
	classCompute, repeatable, obsErrMap, obsErrDefault, opt :=
		readConfig(cl, ocdMap)

	solver := d2solver.New(
		all, unk, classCompute, obsErrMap, obsErrDefault)

	// open obs file
	var f *os.File
	if cl.fnObs == "-" {
		f = os.Stdin
		cl.fnObs = "input stream"
	} else {
		var err error
		f, err = os.Open(cl.fnObs)
		if err != nil {
			exit(err)
		}
		defer f.Close()
	}

	// remainder of main constructs and starts all the concurrent parts
	// of the program.

	// tkChIn supplies trackets by reading obsfile.  It is fed by
	// splitTracklets, running as a separate goroutine.  If splitTracklets
	// encounters an error reading the file, it reports the error on errCh
	// and terminates immediately.
	tkChIn := make(chan *obs.Tracklet)
	errCh := make(chan error)
	go mpc.SplitTracklets(f, ocdMap, tkChIn, errCh)

	// prCh is used to keep processed results in submission order.
	// it is a buffered channel so that a fast worker can drop off the
	// result without waiting for workers ahead of it.  the size of
	// the buffer must be at least maxWorkers, but otherwise isn't critical.
	// Having it somewhat larger allows more results to back up behind
	// a slow worker.  We expect processing time to not vary too much
	// anyway.
	maxWorkers := runtime.GOMAXPROCS(0)
	prCh := make(chan chan string, maxWorkers*2)
	tkChSeq := make(chan *tkSeq)

	// "dispatcher," dispatches tracklets to workers.
	// for each tracklet, attach a return channel that works like a ticket
	// for picking up the result of processing the tracklet.  wait for an
	// available worker, send the tracklet to the worker and drop the
	// ticket in the queue for printing.
	go func() {
		for tk := range tkChIn { // for each tracklet to be solved
			rch := make(chan string, 1) // create return channel for tracklet
			tkChSeq <- &tkSeq{tk, rch}  // queue tracklet for solving
			prCh <- rch                 // queue return channel for printing
		}
		close(prCh)
	}()

	// this function literal, run as a separate goroutine, starts
	// the worker goroutines (solve.)  they are not all started up front,
	// but only as dispatcher (described below) calls for them through avCh.
	// after all, we may have more cores than tracklets.  once it has
	// started the maximum number of workers, it's work is done.
	go func() {
		for n := 0; n < maxWorkers; n++ {
			tk, ok := <-tkChSeq
			if !ok {
				return
			}
			go solve(solver, tk, tkChSeq, repeatable, opt)
		}
	}()

	// column headings, delayed until now to avoid printing column headings
	// only to terminate with an error message if some initialization fails.
	printHeadings(opt)

	// everything is on it's way.  just wait for results and print them
	// as they are available.  prch is our channel of result channels in
	// the correct order.
	for {
		select {
		case err := <-errCh:
			exit(err)
		// wait here for next result channel in processing order
		case rch, ok := <-prCh:
			if !ok {
				return // normal return
			}
			select {
			case err := <-errCh:
				exit(err)
			case r := <-rch:
				fmt.Println(r) // wait here for processing result
			}
		}
	}
}

type tkSeq struct {
	tk  *obs.Tracklet
	rch chan string
}

type lcgRand uint64

var lcga, lcgm lcgRand
var invLcgm float64

func init() {
	lcga = lcgRand(math.Pow(13, 13))
	lcgm = 1
	lcgm <<= 59
	invLcgm = 1 / float64(lcgm)
	lcgm--
}

func (r *lcgRand) Float64() float64 {
	*r = *r * lcga & lcgm
	return float64(*r) * invLcgm
}

// worker process, solves tracklets.
// the first tracklet to solve will be waiting in tkCh.
// additional tracklets are requested by sending tkCh back over avCh.
func solve(solver *d2solver.D2Solver,
	tk *tkSeq, // first tracklet to solve
	tkCh chan *tkSeq, // channel for getting more tracklets
	repeatable bool,
	opt *outputOptions) {
	var rnd d2solver.Rand
	if repeatable {
		rnd = new(lcgRand)
	} else {
		rnd = rand.New(rand.NewSource(time.Now().UnixNano()))
	}
	// this is an infinite loop.  it just runs until the program shuts down.
	for ; ; tk = <-tkCh {
		if repeatable {
			*rnd.(*lcgRand) = 3
		}

		// average whatever magnitudes are there.  default to V=21 if none.
		// this is here rather than d2math.go just to keep d2math.go more
		// general and, as much as practical, not specific to the kind of
		// observations being processed.  code here is specific to the case
		// of typical MPC observations varying in number, often missing
		// magnitudes, and having limiting magnitude around 21.
		var vmag, mSum, mCount float64
		for _, obs := range tk.tk.Obs {
			m := obs.Meas()
			if m.Vmag > 0 {
				mSum += m.Vmag
				mCount++
			}
		}

		if mCount > 0 {
			vmag = mSum / mCount
		} else {
			vmag = 21
		}

		rms, classScores := solver.Solve(tk.tk, vmag, rnd)

		// build output line
		ol := fmt.Sprintf("%7s", tk.tk.Desig)
		if opt.rms {
			if rs := fmt.Sprintf(" %5.2f", rms); len(rs) == 6 {
				ol += rs
			} else {
				ol += " **.**"
			}
		}
		if opt.classPossible {
			// specified columns first
			for _, c := range opt.classColumn {
				cs := classScores[c]
				if opt.raw {
					ol = fmt.Sprintf("%s %3.0f", ol, cs.Raw)
				}
				if opt.noid {
					ol = fmt.Sprintf("%s %3.0f", ol, cs.NoId)
				}
			}
			// then other possibilities
		clist:
			for c := range bin.CList {
				for _, cc := range opt.classColumn {
					if cc == c {
						continue clist // already in a column
					}
				}
				// else output if possible
				cs := classScores[c]
				var pScore float64
				if opt.noid {
					pScore = cs.NoId
				} else {
					pScore = cs.Raw
				}
				if pScore > .5 {
					ol = fmt.Sprintf("%s (%s %.0f)", ol, bin.CList[c].Abbr, pScore)
				} else if pScore > 0 {
					ol = fmt.Sprintf("%s (%s <1)", ol, bin.CList[c].Abbr)
				}
			}
		} else {
			// other possibilities not computed.
			for _, cs := range classScores {
				if opt.raw {
					ol = fmt.Sprintf("%s %3.0f", ol, cs.Raw)
				}
				if opt.noid {
					ol = fmt.Sprintf("%s %3.0f", ol, cs.NoId)
				}
			}
		}

		// processing results sent on private result channel.
		tk.rch <- ol // buffered.  just drop off results and continue
	}
}

type commandLine struct {
	dc    string // config file
	dm    string // model file
	do    string // obscode file
	dp    string // default path
	fnObs string // observations
}

func parseCommandLine() *commandLine {
	// Package path of digest2/go is used for a few things.
	// There's not actually a package here, but it's the parent directory
	// for all commands and packages associated with digest2.
	pp, ppErr := build.Import(parentImport, "", build.FindOnly)
	var cl commandLine
	if ppErr == nil {
		cl.dp = pp.Dir
	}
	dh := flag.Bool("h", false, "")
	dv := flag.Bool("v", false, "")
	flag.StringVar(&cl.dc, "c", "", "")
	flag.StringVar(&cl.dm, "m", "", "")
	flag.StringVar(&cl.do, "o", "", "")
	flag.StringVar(&cl.dp, "p", cl.dp, "")
	flag.Usage = func() {
		os.Stderr.WriteString(`
Usage: digest2 [options] <obsfile>    score observations in file
       digest2 [options] -            score observations from stdin
       digest2 -h                     display help and quick reference
       digest2 -v                     display version and copyright

Options:
       -c <config-file>
       -m <model-file>
       -o <obscode-file>
       -p <path>
`)
		if ppErr == nil {
			os.Stderr.WriteString(`
Default:
       -p=` + pp.Dir + "\n")
		}
	}
	flag.Parse()
	switch {
	case *dh:
		printHelp()
		os.Exit(0)
	case *dv:
		fmt.Println(versionString)
		fmt.Println(copyrightString)
		os.Exit(0)
	case flag.NArg() != 1:
		flag.Usage()
		os.Exit(1)
	}
	cl.fnObs = flag.Arg(0)
	return &cl
}

func readOcd(cl *commandLine) obs.ParallaxMap {
	ocdFile := cl.fixupCP(cl.do, "digest2.obscodes")
	ocdMap, readErr := mpc.ReadOcd(ocdFile)
	if readErr == nil {
		return ocdMap
	}
	// that didn't work.  try getting a fresh copy.
	if err := mpc.FetchOcd(ocdFile); err != nil {
		log.Println(readErr) // show error from read attempt,
		exit(err)            // and error from download attempt
	}
	// retry with downloaded file.  see if this copy works better
	if ocdMap, readErr = mpc.ReadOcd(ocdFile); readErr != nil {
		exit(readErr)
	}
	return ocdMap
}

type outputOptions struct {
	headings, rms, raw, noid, classPossible bool
	classColumn                             []int
}

func readConfig(cl *commandLine, ocdMap obs.ParallaxMap) (classCompute []int, repeatable bool,
	obsErrMap map[string]float64, obsErrDefault float64,
	opt *outputOptions) {
	// default observational error = 1 arc sec
	obsErrDefault = 1 * math.Pi / (180 * 3600)
	obsErrMap = make(map[string]float64)
	opt = new(outputOptions)
	// default configuration
	opt.classPossible = true
	classCompute = make([]int, len(bin.CList))
	for i := range classCompute {
		classCompute[i] = i
	}
	opt.classColumn = classCompute[:4] // MPC Int .. N18
	opt.headings = true
	opt.rms = true
	opt.noid = true
	f, err := os.Open(cl.fixupCP(cl.dc, "digest2.config"))
	if err != nil {
		if cl.dc == "" {
			return
		}
		exit(err)
	}
	defer f.Close()

	rxObserr := regexp.MustCompile(`^[ \t]*(.*?)[ \t]*=[ \t]*(.+)$`)
	parseObsErr := func(s string) (parseErr string) {
		ss := rxObserr.FindStringSubmatch(s)
		if len(ss) != 3 {
			return "Invalid format for obserr."
		}
		oe, err := strconv.ParseFloat(ss[2], 64)
		if err != nil {
			return err.Error()
		}
		if oe > 10 {
			return "Observational error > 10 arc seconds not allowed."
		}
		if ss[1] == "" {
			obsErrDefault = oe * math.Pi / (180 * 3600)
			return ""
		}
		// replace or remove this check if code is changed in
		// mpc.obs80.ParseObs80 to do something more advanced
		// with obs.VMeas.Qual
		_, ok := ocdMap[ss[1]]
		if !ok {
			return "Obscode not recognized."
		}
		obsErrMap[ss[1]] = oe * math.Pi / (180 * 3600)
		return ""
	}

	var rawSpec, classSpec bool
read:
	for lr := bufio.NewReader(f); ; {
		l, isPre, err := lr.ReadLine()
		switch {
		case err == io.EOF:
			if classSpec && !opt.classPossible {
				classCompute = opt.classColumn
			}
			return
		case err != nil:
			exit(err)
		case isPre:
			exit("Unexpected long line in config file.")
		case len(l) == 0:
			continue
		case l[0] == '#':
			continue
		}
		ls := string(l)
		switch ls {
		case "headings":
			opt.headings = true
			continue
		case "noheadings":
			opt.headings = false
			continue
		case "rms":
			opt.rms = true
			continue
		case "norms":
			opt.rms = false
			continue
		case "raw":
			if !rawSpec {
				rawSpec = true
				opt.noid = false
			}
			opt.raw = true
			continue
		case "noid":
			if !rawSpec {
				rawSpec = true
				opt.raw = false
			}
			opt.noid = true
			continue
		case "poss":
			if !classSpec {
				classSpec = true
				opt.classColumn = []int{}
			}
			opt.classPossible = true
			continue
		case "repeatable":
			repeatable = true
			continue
		case "random":
			repeatable = false
			continue
		}
		if strings.HasPrefix(ls, "obserr") {
			errStr := parseObsErr(ls[6:])
			if errStr > "" {
				exit(fmt.Sprintf("%s\nConfig file line: %s", errStr, ls))
			}
			continue
		}
		// only valid possibility left is a class name
		for cx, c := range bin.CList {
			if ls == c.Abbr || ls == c.Heading {
				if !classSpec {
					classSpec = true
					opt.classColumn = []int{}
					opt.classPossible = false
				}
				opt.classColumn = append(opt.classColumn, cx)
				continue read
			}
		}
		exit("Unrecognized line in config file: " + ls)
	}
	return
}

func printHeadings(opt *outputOptions) {
	if opt.headings {
		fmt.Println(versionString)
		// heading line 1
		if opt.raw && opt.noid && len(opt.classColumn) > 0 {
			fmt.Print("-------")
			if opt.rms {
				fmt.Print("  ----")
			}
			for _, c := range opt.classColumn {
				fmt.Printf("   %3s  ", bin.CList[c].Abbr)
			}
			if opt.classPossible {
				fmt.Println(" ---------------")
			} else {
				fmt.Println()
			}
		}
		// heading line 2
		fmt.Printf("Desig. ")
		if opt.rms {
			fmt.Printf("   RMS")
		}
		for _, c := range opt.classColumn {
			if opt.raw && opt.noid {
				fmt.Print(" Raw NID")
			} else {
				fmt.Printf(" %3s", bin.CList[c].Abbr)
			}
		}
		switch {
		case !opt.classPossible:
			fmt.Println()
		case len(opt.classColumn) == 0:
			fmt.Println(" Possibilities")
		default:
			fmt.Println(" Other Possibilities")
		}
	}
}

func (cl *commandLine) fixupCP(fnSpec, fnDefault string) string {
	if fnSpec > "" {
		return fnSpec
	}
	return filepath.Join(cl.dp, fnDefault)
}

func printHelp() {
	fmt.Println(`
Digest2 uses statistical ranging techniques on short arc astrometry to
compute probabilities that observed objects are of various orbit classes.
Input is a file of 80 column MPC-format observations, with at least two
observations per object.  Output is orbit class scores for each object.

Config file keywords:
   headings
   noheadings
   rms
   norms
   raw
   noid
   repeatable
   random
   poss
   obserr

Orbit classes:`)
	for _, c := range bin.CList {
		fmt.Printf("   %3s   %s\n", c.Abbr, c.Heading)
	}
	fmt.Println(`
For full documentation:
   go doc ` + parentImport + `/digest2`)
}

//  reads population model (created by muk)
func readModel(cl *commandLine) (all, unk bin.Model) {
	f, err := os.Open(cl.fixupCP(cl.dm, bin.Mfn))
	if err != nil {
		log.Println(err)
		exit(`Use command "muk" to generate this file.`)
	}
	defer f.Close()
	dec := gob.NewDecoder(f)
	if err = dec.Decode(&bin.QPart); err != nil {
		exit(err)
	}
	dec.Decode(&bin.EPart)
	dec.Decode(&bin.IPart)
	dec.Decode(&bin.HPart)
	dec.Decode(&bin.MSize)
	dec.Decode(&bin.LastH)
	dec.Decode(&all)
	if err = dec.Decode(&unk); err != nil {
		exit(err)
	}
	return
}