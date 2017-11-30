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

