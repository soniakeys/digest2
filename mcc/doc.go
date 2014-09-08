/*
Command mcc computes Matthews correlation coefficient on digest2 results.

Matthews correlation coefficient is a statistic indicating how well
a classifier works.  Here, we are testing how successfully digest2
classifies objects by orbit class.  Compared to similar statistics, MCC
produces a meaningful measure even when the relative number of objects
in the classes (NEOs and non-NEOS, for example) is greatly different.

  Usage: mcc [options] <in-class> <out-of-class> [threshold]
    -c=1: column containing class score
    -v=false: display version and copyright

The command line arguments <in-class> and <out-of-class> are files containing
captured output of digest2.  Prepare these two files as follows:

1.  Pick an orbit class of interest, NEO, for example, and collect a set of
observations for which you know the orbit class truth for each object.
That is, you must have a set of observations which have all been identified
or otherwise used in orbit solutions of sufficient quality to tell if the
objects are in the class of interest or not.

2.  Select just one "tracklet" for each object in this set of observations.
You can select the discovery observations or some other tracklet.  A tracklet
is a few observations from a single observer.  It is an arc long enough that
the instantaneous sky motion is well determined, but hardly ever long enough
for orbit determination.  It spans more than a few minutes but no more than
a few hours.

3.  Partition the selected tracklets into two files, one file for objects
known to be in the orbit class and the other file for objects known to be
out of the orbit class.

4.  Create or edit your digest2.config file so that the class of interest
will be output in a column.  You probably want to make these changes:

- Turn off headings.  mcc can usually ignore the headings, but there is
the strange case where it would try to interpret the digest2 version number
as a class score.

- Turn off the rms column.  mcc has no need for it.  It's best to keep
the output as clean as possible.

- Specify the orbit class of interest.  This minimizes output as well,
but more importantly it allows digest2 to run much faster.  By default
digest2 will compute scores for all orbit classes so that it can display
"Other Possibilities."  If you specify only a single class however, it
computes scores for only that class and consequently runs much faster.

An example config file would look like this,

  noheadings
  norms
  NEO

5.  Run digest2 on the two tracklet files, capturing an output file for each.
These output files become the ones you specify on the mcc command line.

The mcc -c option identifies the output column containing the class of
interest.  The default is 1, corresponding to the minimized output just
described.  (To process the default output of digest2 for NEOs, you would
specify -c=3.)

The optional threshold argument specifies the threshold you use for predicting
if an object is in the orbit class or not.  The default is 50, meaning that
that you interpret a score of 50 or higher as a prediction that the object
is in the class.

mcc allows scores containing a decimal point but it otherwise ignores lines
where it does not find a numeric score in the specified column.  This should
allow it to accept not only output from the current version of digest2, but
also output of other versions or even other programs.
*/
package main
