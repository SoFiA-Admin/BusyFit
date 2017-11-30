/// ____________________________________________________________________ ///
///                                                                      ///
/// BusyFit 0.3.1 (BusyFit.h) - Busy Function fitting programme          ///
/// Copyright (C) 2017 Tobias Westmeier                                  ///
/// ____________________________________________________________________ ///
///                                                                      ///
/// Address:  Tobias Westmeier                                           ///
///           ICRAR M468                                                 ///
///           The University of Western Australia                        ///
///           35 Stirling Highway                                        ///
///           Crawley WA 6009                                            ///
///           Australia                                                  ///
///                                                                      ///
/// E-mail:   tobias.westmeier (at) uwa.edu.au                           ///
/// ____________________________________________________________________ ///
///                                                                      ///
/// This program is free software: you can redistribute it and/or modify ///
/// it under the terms of the GNU General Public License as published by ///
/// the Free Software Foundation, either version 3 of the License, or    ///
/// (at your option) any later version.                                  ///
///                                                                      ///
/// This program is distributed in the hope that it will be useful,      ///
/// but WITHOUT ANY WARRANTY; without even the implied warranty of       ///
/// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the         ///
/// GNU General Public License for more details.                         ///
///                                                                      ///
/// You should have received a copy of the GNU General Public License    ///
/// along with this program. If not, see http://www.gnu.org/licenses/.   ///
/// ____________________________________________________________________ ///
///                                                                      ///

#ifndef BUSYFIT_H
#define BUSYFIT_H

#include <gsl/gsl_version.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_multifit_nlin.h>
#include <gsl/gsl_linalg.h>

#define BUSYFIT_FREE_PARAM 8

class BusyFit
{
public:
	BusyFit();
	~BusyFit();
	
	int     setup(size_t n, double *newData, double *newSigma, bool noPlot = false, bool relax = false, bool verbose = false, size_t offset = 0);
	int     setFreeParameters(bool *newMask);
	int     fit();
	int     fit(double new_a, double new_b1, double new_b2, double new_c, double new_xe0, double new_xp0, double new_w, double new_n);
	int     getResult(double *parValues, double *parUncert, double &chi2, double &chi2dof, size_t &degFree);
	int     getParameters(double &posX, double &w50, double &w20, double &Fpeak, double &Fint);
	int     getParametersUncertainties(double &posX, double &w50, double &w20, double &Fpeak, double &Fint, double &stddev_posX, double &stddev_w50, double &stddev_w20, double &stddev_Fpeak, double &stddev_Fint, int method = 0);
	
private:
	struct data
	{
		size_t  nChannels;
		double *values;
		double *sigma;
		bool    mask[BUSYFIT_FREE_PARAM];   // Free parameter mask; 'true' if parameter free, 'false' if fixed.
	} spectrum;
	// WARNING: A copy of this structure is defined in some of the static functions declared below!
	//          Any change must be applied to all copies, otherwise a segmentation fault will occur!
	
	size_t  offset;          // Note that this is only relevant for plotting and doesn't affect the fit!
	
	double  freeParameters[BUSYFIT_FREE_PARAM];
	double  initialEstimates[BUSYFIT_FREE_PARAM];
	double  parameterUncertainties[BUSYFIT_FREE_PARAM];
	
	gsl_matrix *covar;
	gsl_matrix *cholesky;
	
	double  chi;
	double  chiSquare;       // Reduced chi-squared
	size_t  dof;             // Number of degrees of freedom
	//double  aicc;          // Akaike information criterion with correction for finite sample size
	
	bool    flagNoPlot;
	bool    flagRelax;
	bool    flagVerbose;
	
	// LM fitting functions:
	int     fitWithEstimates();
	int     LMSolver();
	void    print_state(size_t p, size_t iter, gsl_multifit_fdfsolver *s);
	int     printResult();
	int     plotResult();
	
	// Busy Function:
	double  B(double x);
	
	// The following functions have to be static so they can be pointed to:
	static int expb_f   (const gsl_vector *x, void *d, gsl_vector *f);
	static int expb_df  (const gsl_vector *x, void *d, gsl_matrix *J);
	static int expb_fdf (const gsl_vector *x, void *d, gsl_vector *f, gsl_matrix *J);
};

#endif
