/// ____________________________________________________________________ ///
///                                                                      ///
/// BusyFit 0.3.1 (BusyFit.cpp) - Busy Function fitting programme        ///
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

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <cstdlib>      // Needed to invoke gnuplot via system().

#include "helperFunctions.hpp"
#include "BusyFit.hpp"



// ----------- //
// Constructor //
// ----------- //

BusyFit::BusyFit()
{
	std::cout << "\nInitialising BusyFit using GSL " << GSL_MAJOR_VERSION << "." << GSL_MINOR_VERSION << "\n";
	
	offset             = 0;
	spectrum.nChannels = 0;
	spectrum.values    = 0;
	spectrum.sigma     = 0;
	chi                = 0.0;
	chiSquare          = 0.0;
	dof                = 0;
	//aicc               = 0.0;
	flagNoPlot         = false;
	flagRelax          = false;
	flagVerbose        = false;
	
	for(int i = 0; i < BUSYFIT_FREE_PARAM; i++)
	{
		freeParameters[i]         = 0.0;
		initialEstimates[i]       = 0.0;
		parameterUncertainties[i] = 0.0;
		spectrum.mask[i]          = true;
	}
	
	// Allocate memory for covariance matrix:
	covar = gsl_matrix_alloc(BUSYFIT_FREE_PARAM, BUSYFIT_FREE_PARAM);
	
	return;
}



// ---------- //
// Destructor //
// ---------- //

BusyFit::~BusyFit()
{
	// Free memory used for covariance matrix:
	if(covar != 0) gsl_matrix_free(covar);
	
	return;
}



// ------------------------------ //
// Function to initialise BusyFit //
// ------------------------------ //

int BusyFit::setup(size_t n, double *newData, double *newSigma, bool noPlot, bool relax, bool verbose, size_t specOffset)
{
	if(n < BUSYFIT_FREE_PARAM)
	{
		printErrorMessage("Not enough spectral channels to fit.", MSG_ERROR);
		return 1;
	}
	
	spectrum.nChannels = n;
	spectrum.values    = newData;
	spectrum.sigma     = newSigma;
	offset             = specOffset;
	
	flagNoPlot  = noPlot;
	flagRelax   = relax;
	flagVerbose = verbose;
	
	return 0;
}



// ----------------------------------- //
// Function to set up fixed parameters //
// ----------------------------------- //

int BusyFit::setFreeParameters(bool *newMask)
{
	for(int i = 0; i < BUSYFIT_FREE_PARAM; i++) spectrum.mask[i] = newMask[i];
	
	return 0;
}



// ------------------------------------------------- //
// Function to find initial estimates before fitting //
// ------------------------------------------------- //

int BusyFit::fit()
{
	if(spectrum.nChannels == 0 or spectrum.values == 0 or spectrum.sigma == 0)
	{
		printErrorMessage("No valid data found.", MSG_ERROR);
		return 1;
	}
	
	// Define initial estimates:
	initialEstimates[0] = 0.0;
	initialEstimates[1] = 0.5;
	initialEstimates[2] = 0.5;
	initialEstimates[3] = 0.01;
	initialEstimates[4] = 0.0;
	initialEstimates[5] = 0.0;
	initialEstimates[6] = 0.0;
	initialEstimates[7] = 2.0;
	
	size_t posMax   = 0.0;
	
	// Calculate maximum to use as A:
	for(size_t i = 0; i < spectrum.nChannels; i++)
	{
		if(spectrum.values[i] > initialEstimates[0])
		{
			initialEstimates[0] = spectrum.values[i];
			posMax = i;
		}
	}
	
	initialEstimates[0] /= 1.5;      // Make amplitude a little bit smaller, as it actually is the value at the profile centre, not the peak.
	
	// Find flanks to determine position and width:
	size_t leftFlank  = 0;
	size_t rightFlank = 0;
	size_t currentPos = posMax;
	
	while(currentPos > 0 and spectrum.values[currentPos] > 0.2 * initialEstimates[0]) currentPos--;
	
	leftFlank  = currentPos;
	currentPos = posMax;
	
	while(currentPos < spectrum.nChannels - 1 and spectrum.values[currentPos] > 0.2 * initialEstimates[0]) currentPos++;
	
	rightFlank = currentPos;
	
	if(leftFlank >= 0 and rightFlank - leftFlank > 2)
	{
		initialEstimates[4] = static_cast<double>(rightFlank + leftFlank) / 2.0;
		initialEstimates[5] = initialEstimates[4];
		initialEstimates[6] = static_cast<double>(rightFlank - leftFlank - 2) / 2.0;
		
		// xp0 must not be integer to ensure that Jacobian is well-behaved:
		if(initialEstimates[5] == 0.0) initialEstimates[5] = 1.0e-10;
		else if(initialEstimates[5] == floor(initialEstimates[5])) initialEstimates[5] *= 1.000001;
	}
	else
	{
		printErrorMessage("Failed to find line flanks; calculating moments instead.", MSG_WARNING);
		
		// Calculate first moment to use as xe0 and xp0:
		double sum = 0.0;
		
		for(size_t i = 0; i < spectrum.nChannels; i++)
		{
			if(spectrum.values[i] > 0.1 * initialEstimates[0])
			{
				initialEstimates[4] += spectrum.values[i] * static_cast<double>(i);
				sum += spectrum.values[i];
			}
		}
		
		initialEstimates[4] /= sum;
		initialEstimates[5] = initialEstimates[4];    // This is unlikely to be integer, so shouldn't mess up Jacobian.
		
		// Calculate second moment to use as w:
		sum = 0.0;
		
		for(size_t i = 0; i < spectrum.nChannels; i++)
		{
			if(spectrum.values[i] > 0.1 * initialEstimates[0])
			{
				initialEstimates[6] += spectrum.values[i] * (static_cast<double>(i) - initialEstimates[4]) * (static_cast<double>(i) - initialEstimates[4]);
				sum += spectrum.values[i];
			}
		}
		
		initialEstimates[6] = sqrt(initialEstimates[6] / sum);
	}
	
	// Print initial estimates:
	std::cout << "\nUsing the following initial estimates:\n\n";
	
	std::cout << "A   =\t" << initialEstimates[0] << "\n";
	std::cout << "B₁  =\t" << initialEstimates[1] << "\n";
	std::cout << "B₂  =\t" << initialEstimates[2] << "\n";
	std::cout << "C   =\t" << initialEstimates[3] << "\n";
	std::cout << "XE₀ =\t" << initialEstimates[4] << "\n";
	std::cout << "XP₀ =\t" << initialEstimates[5] << "\n";
	std::cout << "W   =\t" << initialEstimates[6] << "\n";
	std::cout << "N   =\t" << initialEstimates[7] << std::endl;
	
	// Copy initial estimates into parameter variables:
	for(int i = 0; i < BUSYFIT_FREE_PARAM; i++) freeParameters[i] = initialEstimates[i];
	
	return fitWithEstimates();
}



// ---------------------------------------------------- //
// Function to ask for initial estimates before fitting //
// ---------------------------------------------------- //

int BusyFit::fit(double new_a, double new_b1, double new_b2, double new_c, double new_xe0, double new_xp0, double new_w, double new_n)
{
	if(spectrum.nChannels == 0 or spectrum.values == 0 or spectrum.sigma == 0)
	{
		printErrorMessage("No valid data found.", MSG_ERROR);
		return 1;
	}
	
	// Ensure that xp0 is not integer as otherwise Jacobian will not be well-behaved:
	if(new_xp0 == 0.0) new_xp0 = 1.0e-10;
	else if(new_xp0 == floor(new_xp0)) new_xp0 *= 1.000001;
	
	// The following is actually possible and allowed in C++:
	freeParameters[0] = initialEstimates[0] = new_a;
	freeParameters[1] = initialEstimates[1] = new_b1;
	freeParameters[2] = initialEstimates[2] = new_b2;
	freeParameters[3] = initialEstimates[3] = new_c;
	freeParameters[4] = initialEstimates[4] = new_xe0;
	freeParameters[5] = initialEstimates[5] = new_xp0;
	freeParameters[6] = initialEstimates[6] = new_w;
	freeParameters[7] = initialEstimates[7] = new_n;
	
	return fitWithEstimates();
}



// -------------------------------------------- //
// Function to fit data using initial estimates //
// -------------------------------------------- //

int BusyFit::fitWithEstimates()
{
	int status = 0;
	
	// Initial fit:
	status = LMSolver();
	
	// Ensure that polynomial order is positive:
	if(freeParameters[7] < 0.0) freeParameters[7] *= -1.0;
	
	// Check: polynomial negative or offset?
	// ### WARNING: This is useless if C and XP0 are fixed!!! (Will need to be checked as well below.)
	if(flagRelax == false and (freeParameters[3] < 0.0 or fabs(freeParameters[4] - freeParameters[5]) > 2.0 * initialEstimates[6]))
	{
		freeParameters[0] = initialEstimates[0];
		freeParameters[1] = initialEstimates[1];
		freeParameters[2] = initialEstimates[2];
		freeParameters[3] = 0.0;
		freeParameters[4] = initialEstimates[4];
		freeParameters[5] = 0.1;   // This must not be integer to ensure that Jacobian is well-behaved!
		freeParameters[6] = initialEstimates[6];
		freeParameters[7] = 0.0;
		
		spectrum.mask[3]  = false;
		spectrum.mask[5]  = false;
		spectrum.mask[7]  = false;
		
		printErrorMessage("Repeating fit without polynomial component (C = 0).", MSG_WARNING);
		
		// Fit again with c and xp0 fixed to 0:
		status = LMSolver();
	}
	
	// Plot and print results:
	plotResult();
	printResult();
	
	if(status == 0)
	{
		// Success:
		if(freeParameters[0] <= 0.0 or freeParameters[1] <= 0.0 or freeParameters[2] <= 0.0 or freeParameters[6] <= 0.0 or freeParameters[3] < 0)
		{
			// Warn that some parameters are negative:
			printErrorMessage("Fit successful, but some parameters zero or negative.", MSG_WARNING);
			return 1;
		}
	}
	else
	{
		// Failure:
		printErrorMessage("Failed to fit spectrum.", MSG_ERROR);
		return 2;
	}
	
	return 0;
}



// --------------------------------------------------- //
// Function to write results to file and generate plot //
// --------------------------------------------------- //

int BusyFit::plotResult()
{
	if(flagNoPlot == true) return 0;   // Check whether plotting is desired.
	
	std::ofstream outputFileSpectrum("busyfit_output_spectrum.txt");
	std::ofstream outputFileFit("busyfit_output_fit.txt");
	
	if(outputFileSpectrum.is_open() and outputFileFit.is_open())
	{
		for(size_t i = 0; i < spectrum.nChannels; i++) outputFileSpectrum << (i + offset) << '\t' << spectrum.values[i] << '\t' << spectrum.values[i] - B(static_cast<double>(i)) << '\n';
		
		for(double x = 0.0; x < static_cast<double>(spectrum.nChannels - 1); x += static_cast<double>(spectrum.nChannels) / 1000.0) outputFileFit << (x + static_cast<double>(offset)) << '\t' << B(x)  << '\n';
		
		outputFileSpectrum.close();
		outputFileFit.close();
	}
	else
	{
		printErrorMessage("Unable to write output files.", MSG_ERROR);
		return 1;
	}
	
	// Invoke gnuplot for real-time plotting of fitted spectrum:
	std::string gnuplotCommand = "echo \'unset key; set grid; set xlabel \"Spectral Channel\"; set ylabel \"Flux\"; plot \"busyfit_output_spectrum.txt\" using 1:3 with points linecolor rgb \"#4060C0\", \"busyfit_output_spectrum.txt\" using 1:2 with histeps linecolor rgb \"#808080\", \"busyfit_output_fit.txt\" using 1:2 with lines linecolor rgb \"#FF0000\"\' | gnuplot -persist";
	
	return system(gnuplotCommand.c_str());
}



// --------------------------------------- //
// Function to print fit results on screen //
// --------------------------------------- //

int BusyFit::printResult()
{
	std::cout << std::setprecision(4) << "No. chan. = " << spectrum.nChannels << "\n";
	std::cout << std::setprecision(4) << "dof       = " << dof       << "\n";
	std::cout << std::setprecision(4) << "chi²      = " << chi * chi << "\n";
	std::cout << std::setprecision(4) << "chi²/dof  = " << chiSquare << "\n";
	//std::cout << std::setprecision(4) << "AICc      = " << aicc      << "\n";
	
	std::cout << std::setprecision(4) << "\nA    =\t" << freeParameters[0]   << "\t±  " << parameterUncertainties[0];
	if(spectrum.mask[0] == false) std::cout << " \033[1;33m(fixed)\033[0m";
	
	std::cout << std::setprecision(4) << "\nB1   =\t" << freeParameters[1]  << "\t±  " << parameterUncertainties[1];
	if(spectrum.mask[1] == false) std::cout << " \033[1;33m(fixed)\033[0m";
	
	std::cout << std::setprecision(4) << "\nB2   =\t" << freeParameters[2]  << "\t±  " << parameterUncertainties[2];
	if(spectrum.mask[2] == false) std::cout << " \033[1;33m(fixed)\033[0m";
	
	std::cout << std::setprecision(4) << "\nC    =\t" << freeParameters[3]   << "\t±  " << parameterUncertainties[3];
	if(spectrum.mask[3] == false) std::cout << " \033[1;33m(fixed)\033[0m";
	
	std::cout << std::setprecision(4) << "\nXE0  =\t" << freeParameters[4] << "\t±  " << parameterUncertainties[4];
	if(spectrum.mask[4] == false) std::cout << " \033[1;33m(fixed)\033[0m";
	
	std::cout << std::setprecision(4) << "\nXP0  =\t" << freeParameters[5] << "\t±  " << parameterUncertainties[5];
	if(spectrum.mask[5] == false) std::cout << " \033[1;33m(fixed)\033[0m";
	
	std::cout << std::setprecision(4) << "\nW    =\t" << freeParameters[6]   << "\t±  " << parameterUncertainties[6];
	if(spectrum.mask[6] == false) std::cout << " \033[1;33m(fixed)\033[0m";
	
	std::cout << std::setprecision(4) << "\nN    =\t" << freeParameters[7]   << "\t±  " << parameterUncertainties[7];
	if(spectrum.mask[7] == false) std::cout << " \033[1;33m(fixed)\033[0m";
	
	std::cout << '\n' << std::endl;
	
	return 0;
}



// -------------------------------------- //
// Function to return fit results to user //
// -------------------------------------- //

int BusyFit::getResult(double *parValues, double *parUncert, double &chi2, double &chi2dof, size_t &degFree)
{
	for(int i = 0; i < BUSYFIT_FREE_PARAM; i++)
	{
		parValues[i] = freeParameters[i];
		parUncert[i] = parameterUncertainties[i];
	}
	
	chi2    = chi * chi;
	chi2dof = chiSquare;
	degFree = dof;
	
	return 0;
}



// ----------------------------------------------- //
// Function to determine uncertainties of physical //
// parameters and return results to user           //
// ----------------------------------------------- //

int BusyFit::getParametersUncertainties(double &posX, double &w50, double &w20, double &Fpeak, double &Fint, double &stddev_posX, double &stddev_w50, double &stddev_w20, double &stddev_Fpeak, double &stddev_Fint, int method)
{
	if(freeParameters[0] == 0.0) return 1;
	
	// Save original parameters:
	double origPar[BUSYFIT_FREE_PARAM];
	for(int i = 0; i < BUSYFIT_FREE_PARAM; i++) origPar[i] = freeParameters[i];
	
	double tmp_posX  = 0.0;
	double tmp_w50   = 0.0;
	double tmp_w20   = 0.0;
	double tmp_Fpeak = 0.0;
	double tmp_Fint  = 0.0;
	
	stddev_posX  = 0.0;
	stddev_w50   = 0.0;
	stddev_w20   = 0.0;
	stddev_Fpeak = 0.0;
	stddev_Fint  = 0.0;
	
	// Determine parameters:
	getParameters(posX, w50, w20, Fpeak, Fint);
	
	
	if(method == 0)
	{
		// Method 1: Linear propagation of covariance matrix
		
		// Allocate memory:
		gsl_vector *P_orig = gsl_vector_alloc(5);
		gsl_vector *P_var  = gsl_vector_alloc(5);
		gsl_matrix *F      = gsl_matrix_alloc(BUSYFIT_FREE_PARAM, 5);
		gsl_matrix *covarD = gsl_matrix_alloc(5, 5);
		gsl_matrix *mtxtmp = gsl_matrix_alloc(BUSYFIT_FREE_PARAM, 5);
		
		gsl_vector_set(P_orig, 0, posX);
		gsl_vector_set(P_orig, 1, w50);
		gsl_vector_set(P_orig, 2, w20);
		gsl_vector_set(P_orig, 3, Fpeak);
		gsl_vector_set(P_orig, 4, Fint);
		
		for(int i = 0; i < BUSYFIT_FREE_PARAM; i++)
		{
			// Add relative offset epsilon to current parameter:
			double epsilon = 1.0e-5 * freeParameters[i];
			if(freeParameters[i] == 0.0) epsilon = 1.0e-5;
			freeParameters[i] += epsilon;
			
			getParameters(tmp_posX, tmp_w50, tmp_w20, tmp_Fpeak, tmp_Fint);
			
			gsl_vector_set(P_var, 0, tmp_posX);
			gsl_vector_set(P_var, 1, tmp_w50);
			gsl_vector_set(P_var, 2, tmp_w20);
			gsl_vector_set(P_var, 3, tmp_Fpeak);
			gsl_vector_set(P_var, 4, tmp_Fint);
			
			// Calculate Jacobian matrix elements:
			for(int j = 0; j < 5; j++) gsl_matrix_set(F, i, j, (gsl_vector_get(P_var, j) - gsl_vector_get(P_orig, j)) / epsilon);
			
			// Restore original parameter prior to offset:
			freeParameters[i] = origPar[i];
		}
		
		// Calculate covariance matrix propagation, C' = F^T * C * F:
		gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, covar, F,      0.0, mtxtmp);
		gsl_blas_dgemm(CblasTrans,   CblasNoTrans, 1.0, F,     mtxtmp, 0.0, covarD);
		
		// ### NOTE: This should automatically take care of fixed parameters, because the 
		// ###       covariance matrix rows and columns of fixed parameters are automatically 
		// ###       set to zero during the fit, including the diagonal elements, so no 
		// ###       further action will be required in this case.
		
		// Extract standard deviations from diagonal elements:
		stddev_posX  = sqrt(gsl_matrix_get(covarD, 0, 0));
		stddev_w50   = sqrt(gsl_matrix_get(covarD, 1, 1));
		stddev_w20   = sqrt(gsl_matrix_get(covarD, 2, 2));
		stddev_Fpeak = sqrt(gsl_matrix_get(covarD, 3, 3));
		stddev_Fint  = sqrt(gsl_matrix_get(covarD, 4, 4));
		
		// Release memory:
		gsl_vector_free(P_orig);
		gsl_vector_free(P_var);
		gsl_matrix_free(F);
		gsl_matrix_free(covarD);
		gsl_matrix_free(mtxtmp);
	} // End method 1
	else
	{
		// Method 2: Variation of parameters and Cholesky decomposition
		
		const int numberIterations = (method >= 1000) ? method : 1000;    // Run at least 1000 iterations.
		int successfulIterations   = 0;
		
		initRand();    // Initialise random number generator.
		
		for(int i = 0; i < numberIterations; i++)
		{
			// Create GSL parameter vector and parameter uncertainty vector (with original, unscaled uncertainties):
			gsl_vector *parameters    = gsl_vector_alloc(BUSYFIT_FREE_PARAM);
			gsl_vector *uncertainties = gsl_vector_alloc(BUSYFIT_FREE_PARAM);
			
			for(int i1 = 0; i1 < BUSYFIT_FREE_PARAM; i1++)
			{
				// Create random parameters, but with mean of 0 and variance of 1 for 
				// the Cholesky decomposition to work. This will be corrected again further down.
				freeParameters[i1] = randGauss(1.0);
				
				// Set GSL parameter and uncertainty vector:
				gsl_vector_set(parameters, i1, freeParameters[i1]);
				gsl_vector_set(uncertainties, i1, parameterUncertainties[i1]);
			}
			
			// Create Cholesky matrix, initially holding a copy of the covariance matrix:
			cholesky = gsl_matrix_alloc(BUSYFIT_FREE_PARAM, BUSYFIT_FREE_PARAM);
			gsl_matrix_memcpy(cholesky, covar);
			
			// Turn covariance matrix (the copy stored in cholesky) into correlation matrix:
			for(int i1 = 0; i1 < BUSYFIT_FREE_PARAM; i1++)
			{
				for(int i2 = 0; i2 < BUSYFIT_FREE_PARAM; i2++)
				{
					// Is any of the parameters fixed or their uncertainty zero?
					if(spectrum.mask[i1] == true and spectrum.mask[i2] == true and gsl_vector_get(uncertainties, i1) > 0.0 and gsl_vector_get(uncertainties, i2) > 0.0)
					{
						// No -> calculate correlation:
						gsl_matrix_set(cholesky, i1, i2, gsl_matrix_get(cholesky, i1, i2) / (gsl_vector_get(uncertainties, i1) * gsl_vector_get(uncertainties, i2)));
					}
					else
					{
						// Yes -> avoid division by zero and set auto-correlation to 1 and cross-correlations to 0:
						if(i1 == i2) gsl_matrix_set(cholesky, i1, i2, 1.0);
						else         gsl_matrix_set(cholesky, i1, i2, 0.0);
					}
				}
			}
			
			// Calculate Cholesky decomposition:
			if(gsl_linalg_cholesky_decomp(cholesky))
			{
				printErrorMessage("Cholesky decomposition failed; correlation matrix not positive-definite.", MSG_ERROR);
			}
			else
			{
				// Multiply Cholesky matrix with parameter vector to obtain correlated parameters:
				gsl_blas_dtrmv(CblasLower, CblasNoTrans, CblasNonUnit, cholesky, parameters);
				
				// ### NOTE: Only the lower triangular section of the matrix must be multiplied (hence _dtrmv), because 
				// ###       for some strange reason the GSL fills the upper triangular section with the transposed 
				// ###       Cholesky matrix values instead of zero!
				
				// Copy vector back into parameter variables and translate/scale to original mean/variance again:
				for(int i1 = 0; i1 < BUSYFIT_FREE_PARAM; i1++) freeParameters[i1] = gsl_vector_get(uncertainties, i1) * gsl_vector_get(parameters, i1) + origPar[i1];
			}
			
			getParameters(tmp_posX, tmp_w50, tmp_w20, tmp_Fpeak, tmp_Fint);
			
			// Ensure that results are meaningful:
			if(!std::isnan(tmp_posX) and !std::isinf(tmp_posX) and !std::isnan(tmp_w50) and !std::isinf(tmp_w50) and !std::isnan(tmp_w20) and !std::isinf(tmp_w20) and !std::isnan(tmp_Fpeak) and !std::isinf(tmp_Fpeak) and !std::isnan(tmp_Fint) and !std::isinf(tmp_Fint) and freeParameters[0] > 0.0 and freeParameters[1] > 0.0 and freeParameters[2] > 0.0 and freeParameters[3] >= 0.0 and freeParameters[6] > 0.0)
			{
				stddev_posX  += (tmp_posX  - posX)  * (tmp_posX  - posX);
				stddev_w50   += (tmp_w50   - w50)   * (tmp_w50   - w50);
				stddev_w20   += (tmp_w20   - w20)   * (tmp_w20   - w20);
				stddev_Fpeak += (tmp_Fpeak - Fpeak) * (tmp_Fpeak - Fpeak);
				stddev_Fint  += (tmp_Fint  - Fint)  * (tmp_Fint  - Fint);
				
				successfulIterations++;
			}
			
			gsl_vector_free(parameters);
			gsl_vector_free(uncertainties);
			gsl_matrix_free(cholesky);
		}
		
		if(successfulIterations > 1)
		{
			stddev_posX  = sqrt(stddev_posX  / static_cast<double>(successfulIterations - 1));
			stddev_w50   = sqrt(stddev_w50   / static_cast<double>(successfulIterations - 1));
			stddev_w20   = sqrt(stddev_w20   / static_cast<double>(successfulIterations - 1));
			stddev_Fpeak = sqrt(stddev_Fpeak / static_cast<double>(successfulIterations - 1));
			stddev_Fint  = sqrt(stddev_Fint  / static_cast<double>(successfulIterations - 1));
		}
		else
		{
			stddev_posX  = 0.0;
			stddev_w50   = 0.0;
			stddev_w20   = 0.0;
			stddev_Fpeak = 0.0;
			stddev_Fint  = 0.0;
			
			printErrorMessage("Calculation of uncertainties failed.", MSG_WARNING);
		}
		
		// Restore original parameters:
		for(int i = 0; i < BUSYFIT_FREE_PARAM; i++) freeParameters[i] = origPar[i];
	} // End method 2
	
	return 0;
}



// -------------------------------------- //
// Measure and return physical parameters //
// -------------------------------------- //

int BusyFit::getParameters(double &posX, double &w50, double &w20, double &Fpeak, double &Fint)
{
	if(freeParameters[0] == 0.0) return 1;
	
	double limitLo = 0.0;
	double limitHi = static_cast<double>(spectrum.nChannels);
	double step    = 0.2;
	
	Fpeak = 0.0;
	
	// Find maximum:
	for(double x = limitLo; x < limitHi; x += step)
	{
		if(B(x) > Fpeak) Fpeak = B(x);
	}
	
	// Determine w50 and w20:
	double current_x  = limitLo;
	double current_x0 = limitLo;
	double w20Lo;
	double w50Lo;
	
	while(current_x < limitHi and B(current_x) < 0.2 * Fpeak) current_x += step;
	current_x0 = current_x - step;
	w20Lo = (0.2 * Fpeak - B(current_x0)) * step / (B(current_x) - B(current_x0)) + current_x0;
	
	while(current_x < limitHi and B(current_x) < 0.5 * Fpeak) current_x += step;
	current_x0 = current_x - step;
	w50Lo = (0.5 * Fpeak - B(current_x0)) * step / (B(current_x) - B(current_x0)) + current_x0;
	
	current_x = limitHi;
	
	while(current_x > w20Lo and B(current_x) < 0.2 * Fpeak) current_x -= step;
	current_x0 = current_x + step;
	w20 = -1.0 * (0.2 * Fpeak - B(current_x0)) * step / (B(current_x) - B(current_x0)) + current_x0 - w20Lo;
	
	while(current_x > w50Lo and B(current_x) < 0.5 * Fpeak) current_x -= step;
	current_x0 = current_x + step;
	w50 = -1.0 * (0.5 * Fpeak - B(current_x0)) * step / (B(current_x) - B(current_x0)) + current_x0 - w50Lo;
	
	// Determine integral and centroid:
	Fint = 0.0;
	posX = 0.0;
	
	for(double x = limitLo; x < limitHi; x += step)
	{
		Fint += B(x);
		posX += B(x) * x;
	}
	
	posX /= Fint;
	Fint *= step;
	
	return 0;
}



// ----------------------------------------------- //
// Actual Levenberg-Marquardt (LM) fitting routine //
// ----------------------------------------------- //

int BusyFit::LMSolver()
{
	int                status;
	unsigned int       iter    = 0;                   // Iteration
	const unsigned int iterMax = 10000;               // Maximum number of iterations
	double             x_init[BUSYFIT_FREE_PARAM];    // Initial estimates of free parameters
	for(size_t i = 0; i < BUSYFIT_FREE_PARAM; i++) x_init[i] = freeParameters[i];
	
	// Setup of GSL solver:
	const gsl_multifit_fdfsolver_type *T;
	gsl_multifit_fdfsolver    *s;
	gsl_multifit_function_fdf  f;
	gsl_vector_view            x = gsl_vector_view_array(x_init, BUSYFIT_FREE_PARAM);
	const gsl_rng_type        *type;
	gsl_rng                   *r;
	
	#if GSL_MAJOR_VERSION >= 2
	gsl_matrix                *J = gsl_matrix_alloc(spectrum.nChannels, BUSYFIT_FREE_PARAM);
	#endif
	
	gsl_rng_env_setup();
	
	type     = gsl_rng_default;
	r        = gsl_rng_alloc(type);
	
	f.f      = &expb_f;
	f.df     = &expb_df;
	f.fdf    = &expb_fdf;
	
	f.n      = spectrum.nChannels;
	f.p      = BUSYFIT_FREE_PARAM;
	f.params = &spectrum;
	
	T        = gsl_multifit_fdfsolver_lmsder;
	s        = gsl_multifit_fdfsolver_alloc(T, spectrum.nChannels, BUSYFIT_FREE_PARAM);
	gsl_multifit_fdfsolver_set(s, &f, &x.vector);
	
	std::cout << "\nFitting \"Busy Function\" B(x) to spectrum:\nB(x) = (A / 4) × [erf(B1 (W + X − XE₀)) + 1]\n       × [erf(B2 (W − X + XE₀)) + 1] × (C |X − XP₀|ᴺ + 1)\n" << std::endl;
	
	if(flagVerbose == true) print_state(BUSYFIT_FREE_PARAM, iter, s);
	
	// Solve the χ² minimisation problem:
	do
	{
		iter++;
		status = gsl_multifit_fdfsolver_iterate(s);
		
		if(flagVerbose == true)
		{
			//std::cout << "Status: " << gsl_strerror(status) << std::endl;
			print_state(BUSYFIT_FREE_PARAM, iter, s);
		}
		
		if(status) break;
		
		status = gsl_multifit_test_delta(s->dx, s->x, 1.0e-4, 1.0e-3);
		// Here, the first epsilon is absolute, the second is relative, 
		// such that |dx_i| < epsabs + epsrel × |x_i|. GSL suggests 1.0e-4 
		// as the default for both, but that leads to difficulties in 
		// conversion for some spectra, so let's be less strict about 
		// the relative error.
	}
	while(status == GSL_CONTINUE and iter < iterMax);
	
	if(flagVerbose == false) std::cout << "Iterations completed: " << iter << std::endl;
	
	#if GSL_MAJOR_VERSION < 2
	gsl_multifit_covar(s->J, 0.0, covar);
	#else 
	gsl_multifit_fdfsolver_jac(s, J);
	gsl_multifit_covar(J, 0.0, covar);
	#endif
	
	int nFreePar = 0;
	
	for(int i = 0; i < BUSYFIT_FREE_PARAM; i++)
	{
		if(spectrum.mask[i] == true) nFreePar++;
	}
	
	chi = gsl_blas_dnrm2(s->f);
	dof = spectrum.nChannels - nFreePar;
	double cc  = GSL_MAX_DBL(1.0, chi / sqrt(static_cast<double>(dof)));
	// ### Why is this? This may cause problems in other parts of the software, 
	// ### because σ(i) != sqrt(Cov(i,i)) any longer! Possible motivation: 
	// ### χ²/dof > 1 in the ideal situation means that the model does not 
	// ### exactly describe the data due to intrinsic variations in the data 
	// ### that are not modelled. Multiplication with cc will take this additional 
	// ### uncertainty into account. Alternatively, one could set cc = 1.0.
	
	// Extract results:
	chiSquare = chi * chi / static_cast<double>(dof);
	//aicc      = chi * chi + static_cast<double>(2 * nFreePar * spectrum.nChannels) / static_cast<double>(dof - 1);
	
	for(int i = 0; i < BUSYFIT_FREE_PARAM; i++)
	{
		freeParameters[i]         = gsl_vector_get(s->x, i);
		parameterUncertainties[i] = cc * sqrt(gsl_matrix_get(covar, i, i));
	}
	
	std::cout << "Status: " << gsl_strerror(status) << "\n" << std::endl;
	
	gsl_multifit_fdfsolver_free(s);
	gsl_rng_free(r);
	#if GSL_MAJOR_VERSION >= 2
	gsl_matrix_free(J);
	#endif
	
	return status;
}



// ----------------------------------- //
// Definition of function to be fitted //
// ----------------------------------- //

// (This function is static!)
int BusyFit::expb_f(const gsl_vector *x, void *d, gsl_vector *f)
{
	// For B(x):
	double a   = gsl_vector_get(x, 0);
	double b1  = gsl_vector_get(x, 1);
	double b2  = gsl_vector_get(x, 2);
	double c   = gsl_vector_get(x, 3);
	double xe0 = gsl_vector_get(x, 4);
	double xp0 = gsl_vector_get(x, 5);
	double w   = gsl_vector_get(x, 6);
	double n   = gsl_vector_get(x, 7);
	
	struct data
	{
		size_t  nChannels;
		double *values;
		double *sigma;
		bool    mask[BUSYFIT_FREE_PARAM];
	};
	
	size_t  nChan = ((struct data *)d)->nChannels;
	double *y     = ((struct data *)d)->values;
	double *sigma = ((struct data *)d)->sigma;
	
	for(size_t i = 0; i < nChan; i++)
	{
		double x = static_cast<double>(i);
		double Yi = (a / 4.0) * (erf(b1 * (w + x - xe0)) + 1.0) * (erf(b2 * (w + xe0 - x)) + 1.0) * (c * pow(fabs(x - xp0), n) + 1.0);
		
		gsl_vector_set(f, i, (Yi - y[i]) / sigma[i]);
	}
	
	return GSL_SUCCESS;
}



// ------------------------------------------------------ //
// Definition of Jacobian matrix of function to be fitted //
// ------------------------------------------------------ //

// (This function is static!)
int BusyFit::expb_df(const gsl_vector *x, void *d, gsl_matrix *J)
{
	// For B(x):
	double a   = gsl_vector_get(x, 0);
	double b1  = gsl_vector_get(x, 1);
	double b2  = gsl_vector_get(x, 2);
	double c   = gsl_vector_get(x, 3);
	double xe0 = gsl_vector_get(x, 4);
	double xp0 = gsl_vector_get(x, 5);
	double w   = gsl_vector_get(x, 6);
	double n   = gsl_vector_get(x, 7);
	
	struct data
	{
		size_t  nChannels;
		double *values;
		double *sigma;
		bool    mask[BUSYFIT_FREE_PARAM];
	};
	
	size_t  nChan = ((struct data *)d)->nChannels;
	double *sigma = ((struct data *)d)->sigma;
	int     mask[BUSYFIT_FREE_PARAM];
	
	for(int i = 0; i < BUSYFIT_FREE_PARAM; i++)
	{
		mask[i] = ((struct data *)d)->mask[i];
	}
	
	for(size_t i = 0; i < nChan; i++)
	{
		// Jacobian matrix J(i,j) = dfi / dxj,
		// where fi = (Yi - yi)/sigma[i],
		//       Yi = value of function to be fitted
		// and the xj are the free parameters
		
		double x = static_cast<double>(i);
		
		double erf1   = erf(b1 * (w + x - xe0)) + 1.0;
		double erf2   = erf(b2 * (w + xe0 - x)) + 1.0;
		double para   = c * pow(fabs(x - xp0), n) + 1.0;
		double exp1   = exp(-1.0 * (w + x - xe0) * (w + x - xe0) * b1 * b1);
		double exp2   = exp(-1.0 * (w - x + xe0) * (w - x + xe0) * b2 * b2);
		double sqrtPi = sqrt(M_PI);
		
		gsl_matrix_set(J, i, 0, static_cast<double>(mask[0]) * (1.0 / 4.0) *           erf1 * erf2 * para / sigma[i]);
		gsl_matrix_set(J, i, 1, static_cast<double>(mask[1]) * (a / (2.0 * sqrtPi)) *  erf2 * para * (w + x - xe0) * exp1 / sigma[i]);
		gsl_matrix_set(J, i, 2, static_cast<double>(mask[2]) * (a / (2.0 * sqrtPi)) *  erf1 * para * (w - x + xe0) * exp2 / sigma[i]);
		gsl_matrix_set(J, i, 3, static_cast<double>(mask[3]) * (a / 4.0) *             erf1 * erf2 * pow(fabs(x - xp0), n) / sigma[i]);
		gsl_matrix_set(J, i, 4, static_cast<double>(mask[4]) * (a / (2.0 * sqrtPi)) * (erf1 * para * b2 * exp2 - erf2 * para * b1 * exp1) / sigma[i]);
		gsl_matrix_set(J, i, 5, static_cast<double>(mask[5]) * (-1.0 * a / 4.0) *      erf1 * erf2 * n * c * pow(fabs(x - xp0), n - 1.0) * ((x - xp0) / fabs(x - xp0)) / sigma[i]);
		gsl_matrix_set(J, i, 6, static_cast<double>(mask[6]) * (a / (2.0 * sqrtPi)) * (erf1 * para * b2 * exp2 + erf2 * para * b1 * exp1) / sigma[i]);
		gsl_matrix_set(J, i, 7, static_cast<double>(mask[7]) * (a / 4.0) *             erf1 * erf2 * (para - 1.0) * log(fabs(x - xp0)) / sigma[i]);
	}
	
	return GSL_SUCCESS;
}



// --------------------------------------- //
// Function to call expb_f() and expb_df() //
// --------------------------------------- //

// (This function is static!)
int BusyFit::expb_fdf(const gsl_vector *x, void *d, gsl_vector *f, gsl_matrix *J)
{
	// For B(x):
	expb_f(x, d, f);
	expb_df(x, d, J);
	
	return GSL_SUCCESS;
}



// ------------------------------------------- //
// Function to print results of each iteration //
// ------------------------------------------- //

void BusyFit::print_state(size_t p, size_t iter, gsl_multifit_fdfsolver *s)
{
	std::cout << "It " << iter << ": " << std::setprecision(4);
	for(size_t i = 0; i < p; i++) std::cout << gsl_vector_get(s->x, i) << " ";
	std::cout << "|F(x)| = " << gsl_blas_dnrm2(s->f) << std::endl;
	
	return;
}



// ------------------- //
// Busy Function, B(x) //
// ------------------- //

double BusyFit::B(double x)
{
	return (freeParameters[0] / 4.0) * (erf(freeParameters[1] * (freeParameters[6] + x - freeParameters[4])) + 1.0) * (erf(freeParameters[2] * (freeParameters[6] + freeParameters[4] - x)) + 1.0) * (freeParameters[3] * pow(fabs(x - freeParameters[5]), freeParameters[7]) + 1.0);
}
