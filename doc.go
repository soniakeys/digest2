/*
Command digest2 uses statistical ranging techniques to compute chances that an
object is of various orbit classes.

Contents

Version 0.15

  Program overview
  Installing from the Internet
  Command line usage
  Configuring file locations
  File formats
  Algorithm outline


Program overview

Input is a file of 80 column MPC-format observations, with at least
two observations per object.  Output is orbit class scores for each object.

The MPC observation format is documented at
http://www.minorplanetcenter.net/iau/info/OpticalObs.html.  This is an ASCII
encoded format.  There is no allowance for non-ASCII characters.

The program is provided as Go source code, and also as a C program with
nearly identical functionality.  Comparing the two, there are some small
differences in operation but the scores produced should be practically
identical.  The Go and C sources are independent.  You can use either one
without the other.

Sample run:

Here are a few observations of known NEOs (with made up designations.)

     NE00030  C2004 09 16.15206 16 13 11.57 +20 52 23.7          21.1 Vd     291
     NE00030  C2004 09 16.15621 16 13 11.34 +20 52 16.8          20.8 Vd     291
     NE00030  C2004 09 16.16017 16 13 11.13 +20 52 09.6          20.7 Vd     291
     NE00199  C2007 02 09.24234 06 08 06.06 +43 13 26.2          20.1  c     704
     NE00199  C2007 02 09.25415 06 08 05.51 +43 13 01.7          20.1  c     704
     NE00199  C2007 02 09.26683 06 08 04.80 +43 12 37.5          19.9  c     704
     NE00269  C2003 01 06.51893 12 40 50.09 +18 27 46.9          21.4 Vd     291
     NE00269  C2003 01 06.52850 12 40 50.71 +18 27 46.1          21.8 Vd     291
     NE00269  C2003 01 06.54359 12 40 51.68 +18 27 42.5          21.9 Vd     291

You put them in a file, say fmo.obs, then type "digest2 fmo.obs" and
get the following output:

  Digest2 version 0.15 Go source
  Desig.    RMS Int NEO N22 N18 Other Possibilities
  NE00030  0.15 100 100  38   0
  NE00199  0.56  98  97  18   0 (MC 2) (JFC 1)
  NE00269  0.42  18  17   3   0 (MC 9) (Hun 6) (Pho 26) (MB1 <1) (Han <1) (MB2 30) (MB3 12) (JFC 1)

Only considering these short arcs, digest2 predicts that the first two
objects are almost certain to be NEOs.  The last one, with a NEO score of
only 17 shows little chance of being a NEO.  Pho 26, MB2 30, and MB3 12
show that it is much more likely to have a Phocaea-like or main belt orbit.
The scores may not add to 100 because they are computed independently and
because all possible orbits are not categorized.  Orbit classes are based
on orbital elements and are either known dynamic populations like "Hungarias"
or popular classifications like "NEO."

The RMS figure is a root-mean-square of residuals in arc seconds of the
observations against a great circle fit.  A high RMS could indicate either bad
observations or significant great circle departure over the observation arc.
It can thus be used as a quick check that scores are meaningful.  If the RMS
is low, scores are meaningful.  Scores are still meaningful if the RMS is high
due to great circle departure, but digest2 offers no way of distinguishing this
case from one where observations are bad.


Installing from the Internet

You need Go 1 or later installed and configured.  If you are new to Go, see
http://golang.org/ and http://golang.org/doc/install.  Highly recommended
is setting GOPATH as described at http://golang.org/doc/code.html#tmp_13,
including adding <gopath>/bin to your PATH environment variable.
The section http://golang.org/doc/code.html#remote explains this in more
detail.  The Go programs digest2 and muk are remote packages.

Then type the following command:

    go get code.google.com/p/digest2/go/{digest2,muk}

This downloads the source code for the two commands, and also all of their
subordinate packages.  It also compiles and installs the commands and
packages.

Muk is a program that initializes files for digest2.  It requires astorb.dat
for input and will attempt to download a copy if it is not found.  Because
the download is time consuming, you may wish to use an existing copy of
astorb.dat.

   muk -h

shows the default value for the -a option, where muk looks for astorb.dat
by default.  To avoid the download, you can copy or link astorb.dat to this
location, or you can specify a different path or filename with -a.

Run muk by typing,

	muk

or

	muk -a <path>

Command line usage

The main executable is digest2.  Invoking the program without command line
arguments (or with invalid arguments) shows this usage prompt.

  Usage: digest2 [options] <obsfile>    Score observations in file.
         digest2 [options] -            Score observations from stdin.
         digest2 -h                     Display help and quick reference.
         digest2 -v                     Display version and copyright.

  Options:
       -c <config-file>
       -m <model-file>
       -o <obscode-file>
       -p <path>

The help information lists a quick reference to keywords and orbit classes
allowed in the configuration file.  The configuration file is explained
below under File Formats.


Configuring file locations

When digest2 runs, it reads observations either from a file specified on the
command line or from stdin.

It also reads from two other required data files and an optional configuration
file.  The default location for these three files is different for the C and Go
programs.  For the Go program, their initial location is determined by GOPATH.
The full path to these files is shown at the end of the usage message.
Type digest2 without command line arguments and read the default location as
the -p default.

You can maintain these three files in their default location or you can
relocate them and specify their locations with command line options.

	File              Command line option
	digest2.obscodes  -o
	digest2.gmodel    -m
	digest2.config    -c

A configuration file is required to be present if -c is used.

You can use the -p option to specify a common path to the three files,
overriding the default location.  They will be accessed with their default
names but in the specified location.

If you specify -p in combination with -c, -o, or -m, the path specified
with the -c, -o, or -m option takes precedence.  That is, the path specified
with -p is not joined with with a file name specified with -c, -o, or -m.


File formats

Observations, whether supplied in a file or through stdin, should contain
observations in the MPC 80 column observation format and nothing else.
The observations should be sorted first by designation and then by time
of observation, and there should be at least two observations of each object.

digest2.obscodes is a text file containing observatory codes in the standard
MPC format.  If the file is missing, digest2 will access the Minor Planet
Center web site and download a copy.  This normally happens the first time
you run digest2 and is normally quick and is not noticable.

digest2.gmodel is a binary file generated by the program muk, as described
above.  See the full documentation on muk with,

	go doc code.google.com/p/digest2/go/muk

digest2.config, the optional configuration file, is a text file with a simple
format.  Empty lines and lines beginning with # are ignored.  Other lines must
contain either a keyword or an orbit class.

Allowable keywords:

   headings
   noheadings
   rms
   norms
   raw
   noid
   repeatable
   random
   obserr
   poss

Headings and the rms column can be turned off if desired.

Keywords raw and noid determine the score produced as described below under
Algorithm Outline. The default is noid.  If both keywords are present, both
scores are output.

The keywords repeatable and random determine if program output is
strictly repeatable or can vary slightly from one run to the next.
The program uses a Monte Carlo method.  By default, the pseudo random
number generator is seeded randomly.  When the keyword repeatable is
used, it is reseeded with a constant value for each tracklet, yielding
repeatable scores.

Keyword obserr specifies the amount of observational error that the algorithm
should allow for.  It is specified in arc seconds as in,

  obserr=0.7

The default, if no obserr is specified, is 1.0 arc seconds.
Obserr may be specified for individual observatory codes as in,

  obserrF51=.3
  obserr 704 = 1

As shown, white space is optional.

The keyword poss specifies to output the "Other Possibilities" column.
By default, other possibilities are suppressed if orbit classes are
explicitly specified.

Orbit classes:

   Abbr.  Long Form
   ---    -------------
   Int    MPC interest.
   NEO    NEO(q < 1.3)
   N22    NEO(H <= 22)
   N18    NEO(H <= 18)
   MC     Mars Crosser
   Hun    Hungaria gr.
   Pho    Phocaea group
   MB1    Inner MB
   Pal    Pallas group
   Han    Hansa group
   MB2    Middle MB
   MB3    Outer MB
   Hil    Hilda group
   JTr    Jupiter tr.
   JFC    Jupiter Comet

Listing an orbit class limits scoring to only the listed classes.
Other possibilities are not computed or listed.  Either abbreviations or
long forms may be used.  In any case they must be spelled exactly as
shown.

Example 1:

  Int
  Neo
  N22
  N18
  poss

This is equivalent to default program behavior without a config file.

Example 2:

  # just three
  NEO
  Hun
  JTr

Program output is,

  Digest2 version 0.15 Go source
  Desig.    RMS NEO Hun JTr
  NE00030  0.15 100   0   0
  NE00199  0.56  97   0   0
  NE00269  0.42  17   6   0

The program runs considerably faster in this case, as it computes scores for
only these three classes and not all possible classes.

Example 3:

  noheadings
  norms
  N22

Output:

  NE00030  38
  NE00199  18
  NE00269   3

This might be useful for generating results to be analyzed by another program,
See the program mcc, for example, included with digest2.


Algorithm outline

1.  For each object, the program computes a motion vector from the
first and last observation, and computes a V magnitude from whatever
photometry is provided.

2.  It then generates a number of orbits that are consistent with
the computed motion, complete with absolute magnitudes consistent with
the computed V magnitude.

3.  Each orbit is located within a binned, or histogram, model of the
solar system.  The model is binned in the four dimensions of q, e, i, and H.
As the bin is determined for each orbit, a tag is set for that bin.
Additionally, each orbit is tested for each configured class and a separate
tag is set, by bin, for each class.

4.  A search algorithm is used to generate orbits covering the entire range
of possible orbits, and tag corresponding possible bins.  As the algorithm
generates variant orbits, it checks if the orbits are yielding new bin
taggings, either in general or of specific orbit classes.  The algorithm
terminates after reaching a point of diminishing returns in finding bins.

5.  The histogram contains object population counts by bin.  For each bin
there is a population of each orbit class, and also a complete population
count.  After orbit search, the sum of complete population of tagged bins
represents the number of possible candidate objects in the solar system.
The population sum of tagged bins of a specified class represents the number
of possible candidates of that class.  The class sum divided by the complete
sum represents the fraction of candidate objects that are of that class.
This fraction is multiplied by 100 and output as the "raw" percentage.

6.  No-ID scores are computed similarly with a parallel histogram.
In it, population counts are not of the expected complete population of the
solar system, but of the expected yet-undiscovered population.  This
population is computed by reducing the modeled complete population by known
discoveries.  As the intended context of the no-ID score is after attempted
object identification, the selected known population is that which is readily
identifiable.  The current criteria used is sky uncertainty < 1' arc.
The uncertainty parameter selected for this comparison is field 24 of
astorb.dat, which is a peak ephemeris uncertainty over a 10 year period.

-------------
Public domain 2014, Smithsonian Astrophysical Observatory.
*/
package main
