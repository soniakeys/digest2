/*
Command muk prepares a data file for digest2.

Usage

Command line options:

  muk                                  Use default location for astorb.dat.
  muk -v                               Display version and copyright.
  muk -a=<muk source path>/astorb.dat  Specify astorb.dat path or file name.

Input

The program reads two files:

    s3m.dat, the S3M binned model.
    astorb.dat, the Lowell orbit catalog.

s3m.dat is required to be in the muk package source directory.  A copy of
this file is included with the source code, but it can also be regenerated
from the original S3M data files by the the program s3mbin.

astorb.dat will also be taken from the muk package source directory by default,
but if it is not found an attempt will be made to fetch it with wget.
Fetching a copy of astorb.dat is time consuming but only has to be done once.

If you happen to have a copy of astorb.dat in another location or with another
file name, you can specify this with the -a option.  If the file is not found
in this case, the program fails with an error message and does not attempt
to download astorb.dat.

Output

The output is a single file, digest2.gmodel, containing a merging of the inputs
in a format readily useful to digest2.  This format is the Go "gob" format, a
binary format that is not human readable.

-------------
Public domain.
*/
package main
