# BusyFit

## Introduction

BusyFit is a stand-alone software for fitting the so-called Busy Function
(Westmeier et al. 2014, MNRAS, 438, 1176) to the integrated spectrum of a
galaxy for the purpose of parametrisation.  It is entirely written in C++
and licensed under the GNU General Public License (version 3 or later).

Another Busy Function fitting software was developed by Russell Jurek and
is available from https://github.com/RussellJurek/busy-function-fitting/.


## Requirements

The following software packages  are required to successfully compile and
run BusyFit:

* Linux/Unix operating system (e.g. Ubuntu, Mac OS X) with graphical user
  interface (optional; required for plotting the fitted spectrum) and
  shell (e.g. bash, tcsh)
* GNU Compiler Collection (GCC), including the C++ compiler (g++)
* GNU Scientific Library (GSL), including development (dev) packages
* GNU less (optional; required for help system)
* gnuplot (optional; required for plotting the fitted spectrum)

If you are running a GNU/Linux system  on your computer,  you will likely
have all of these software packages either already installed or available
through your system’s package manager.  In this case,  there should be no
need to manually download and compile any of them.

The GNU less and gnuplot requirements are optional. The former is used to
display a detailed  user manual  when typing  `busyfit -help`,  while the
latter is needed to automatically plot the fitted spectrum on the screen.
Plotting can be disabled by invoking BusyFit with the `-noplot` option.


## Installation

To install BusyFit on your computer, open a terminal window and change to
the folder where the downloaded archive was saved. To unpack the archive,
type

    tar -zxvf busyfit-[ver].tar.gz

where `[ver]` is the version number of the package you downloaded,  e.g.
0.3. This will unpack all files into a new directory called `BusyFit`.

Assuming that the GNU Scientific Library is installed at a location where
it can be found  by the compiler  and linker,  BusyFit can then simply be
compiled  by  running the  `compile.sh`  shell  script  coming  with this
package. Change into the new `BusyFit` sub-directory and type

    chmod 764 compile.sh
    ./compile.sh

This should compile the source code and generate an executable file named
`busyfit`.  To make the software accessible across the entire system, you
may want to create a symbolic link in `/usr/bin`.  If using Linux, simply
type

    sudo ln -s [installation directory]/busyfit /usr/bin/busyfit

where `[installation directory]` is the full path  of the directory where
the `busyfit` executable file is located. Alternatively, you could create
an alias in the  `.alias`,  `.cshrc`,  or  `.bashrc`  file  (wherever you
usually store your aliases) in your home directory.

The default installation directory for the GSL on Linux systems should be
`/usr/lib` for `libgsl.a` and `libgslcblas.a`, and `/usr/include/gsl` for
all header files.  If that  is not  the case  on your computer,  then the
correct paths must be passed on to the GCC by modifying the  `compile.sh`
script.  This can be achieved by adding the  `-L /path1`  and `-I /path2`
options to the three `g++` calls,  where `/path1` is the full path to the
directory where  `libgsl.a`  is located,  and  `/path2`  is the full path
(without the final  `gsl`  folder)  to the location  where the GSL header
files are stored (e.g., `gsl_rng.h`).  This should enable GCC to find all
required GSL files.


## Documentation

The BusyFit software package  comes with its own internal help system and
user manual. A simple overview of all available options will be displayed
by calling  `busyfit`  without  any options  or input  file name.  A more
detailed user manual can be displayed by typing `busyfit -help`.

The following  example illustrates  how to use BusyFit  in its most basic
form:

    busyfit -c 3 4 spectrum.txt

The spectrum to be fitted is read from a simple text file (`spectrum.txt`
in the example above)  where the velocity and flux values of the spectrum
must be located in two columns the positions of which  are specified with
the `-c` option  (the third column for velocities,  the fourth column for
fluxes in the example above).  Without the `-c` option, BusyFit will read
the first two columns  of the input file by default.  Please see the user
manual for a  detailed description  of all  available options,  including
specifying the noise level, providing initial estimates of the free para-
meters, fixing individual parameters, etc.


## Version history


* BusyFit 0.3.1
  * Released 29/11/2017
  * Fixed a bug that resulted in SoFiA crashing under GSL 2.0 or higher
    due to the Jacobian matrix being assigned the wrong dimensions.
* BusyFit 0.3
  * Released 12/09/2017
  * It is now possible to specify a spectral channel range to be fitted.
  * BusyFit is now compatible with GSL 2.0 or higher.
* BusyFit 0.2
  * Released 10/09/2015
  * BusyFit now fits the exponent of the polynomial component as well.
  * BusyFit can now be made to automatically measure the rms noise of the
    spectrum under certain conditions.
  * It is now possible to specify the column containing the spectral axis
    of the spectrum. BusyFit will use this information to provide all
    physical parameters in their native spectral units (e.g. velocity in
    km/s) rather than just channels.
* BusyFit 0.1
  * Released 25/09/2013
  * Initial, stable release of the package with basic functionality and
    built-in help system


## Copyright and licence

If you publish  results  that have been  derived  with  the help  of this
software, I would appreciate a reference to the original paper describing
the Busy Function:

> The busy function: a new analytic function for describing the
> integrated 21-cm spectral profile of galaxies.
> Westmeier, T., Jurek, R., Obreschkow, D., Koribalski, B. S.,
> Staveley-Smith L., 2014, MNRAS, 438, 1176

> http://adsabs.harvard.edu/abs/2014MNRAS.438.1176W

BusyFit is free software:  you can redistribute it and/or modify it under
the terms  of the  GNU General Public License  as published  by the  Free
Software Foundation, either version 3 of the License, or (at your option)
any later version.

BusyFit is distributed  in the hope  that it will be useful,  but without
any warranty;  without even the  implied warranty  of merchantability  or
fitness for a particular purpose.  See the GNU General Public License for
more details.

You should have received a copy  of the GNU General Public License  along
with BusyFit. If not, see http://www.gnu.org/licenses/.

© 2017 Tobias Westmeier
