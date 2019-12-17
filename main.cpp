/// ____________________________________________________________________ ///
///                                                                      ///
/// BusyFit 0.3.2 (main.cpp) - Busy Function fitting programme           ///
/// Copyright (C) 2019 Tobias Westmeier                                  ///
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
#include <string>
#include <fstream>
#include <vector>
#include <algorithm>

#include "helperFunctions.hpp"
#include "BusyFit.hpp"

int main(int argc, char *argv[]);
int printHelpMessage();



// ------------- //
// Main function //
// ------------- //

int main(int argc, char *argv[])
{
	std::string filename;
	std::vector <double> frequency;
	std::vector <double> spectrum;
	double init_par[BUSYFIT_FREE_PARAM] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};           // Initial estimates of parameters; will be initialised later.
	bool   mask[BUSYFIT_FREE_PARAM]     = {true, true, true, true, true, true, true, true};   // Which parameters are free (true) or fixed (false)?
	bool   initialGuess  = false;
	bool   noPlot        = false;
	bool   relax         = false;
	bool   verbose       = false;
	bool   flagNoise     = false;
	int    iterations    = 0;          // Iterations for parameter variation method; 0 = use error propagation method instead (default).
	int    columnSpec    = 0;          // Spectral axis. Default: first column.
	int    column        = 1;          // Flux values. Default: second column.
	// NOTE: columnSpec and column must be signed to be able to handle the situation where a user specifies a value of ≤ 0!
	size_t windowLower   = 0;          // Lower and upper boundary of spectral window to be fitted.
	size_t windowUpper   = 0;          // 0 = use full spectrum
	double noise         = 0.01;       // Default rms = 0.01; must later be specified by user.
	double specInc       = 1.0;        // Channel width in spectral units.
	std::string fluxUnit = "flux";
	std::string specUnit = "spec";
	
	initRand();         // Initialise random number generator.
	
	if(argc < 2)
	{
		std::string message = "\n\033[1mUsage:\033[0m busyfit [OPTIONS] FILE\n\n";
		
		message.append("  Note that there is no specific order for the options or file name. The\n");
		message.append("  spectrum must be stored on a \033[1;4mREGULAR GRID\033[0m and in the \033[1;4mCORRECT ORDER\033[0m, without\n");
		message.append("  any gaps or jumps, as otherwise all results will be wrong and meaningless!\n");
		message.append("  Redshift-dependent effects will \033[1;4mNOT\033[0m be taken into account.\n\n");
		
		message.append("\033[1mAvailable options\033[0m\n\n");
		message.append("  -c C1 [C2]   Spectral (C1) and flux (C2) column numbers in input file.\n");
		message.append("               Default: C1 = 1. The 2nd value is optional and defaults\n");
		message.append("               to C2 = C1 + 1.\n\n");
		message.append("  -h, -help    Display a more detailed user manual and exit.\n\n");
		message.append("  -n rms       The rms noise level of the spectrum. Default is 0.01. Use a\n");
		message.append("               value of ≤ 0 to automatically determine the rms noise.\n\n");
		message.append("  -noplot      Do not plot the fit results on the screen. Useful when \n");
		message.append("               used in batch mode or when gnuplot is unavailable.\n\n");
		message.append("  -p P1 … P8   Initial estimates of the free parameters, A, B₁, B₂, C, \n");
		message.append("               XE₀, XP₀, W and N, separated by spaces. If not specified, \n");
		message.append("               initial parameters will automatically be determined. \n");
		message.append("               Adding ‘f’ at the end will fix a parameter, e.g. ‘2.5f’.\n\n");
		message.append("  -relax       Do not check for negative or shifted polynomial component.\n\n");
		message.append("  -u iter      Apply Monte-Carlo method to calculate uncertainties for \n");
		message.append("               observational parameters, using ‘iter’ iterations. Default \n");
		message.append("               is to use a much faster error propagation method.\n\n");
		message.append("  -v           Use verbose mode; show progress information during fitting.\n\n");
		message.append("  -w C1 C2     Define channel range to be fitted (zero-based). Default is\n");
		message.append("               all channels. Set to 0 for first/last channel, respectively.\n");
		
		std::cerr << message << std::endl;
		return 0;
	}
	
	// Check command line arguments:
	for(int i = 1; i < argc; i++)
	{
		// ----------------
		// HELP (-h, -help)
		// ----------------
		if(std::string(argv[i]) == "-h" or std::string(argv[i]) == "-help")
		{
			return printHelpMessage();
		}
		// ------------
		// COLUMNS (-c)
		// ------------
		else if(std::string(argv[i]) == "-c")
		{
			if(i < argc - 1)
			{
				// First read spectral column:
				i++;
				columnSpec = stringToNumber<int>(argv[i]) - 1;
				
				if(columnSpec < 0)
				{
					columnSpec = 0;
					printErrorMessage("Invalid spectral column number; using default value.", MSG_WARNING);
				}
				
				// Then read/determine flux column:
				if(i < argc - 1)
				{
					i++;
					column = stringToNumber<int>(argv[i]) - 1;
					
					if(column < 0)
					{
						i--;
						column = columnSpec + 1;
					}
				}
				else column = columnSpec + 1;
			}
			else printErrorMessage("Column numbers missing; using default values.", MSG_WARNING);
		}
		// ----------
		// NOISE (-n)
		// ----------
		else if(std::string(argv[i]) == "-n")
		{
			if(i < argc - 1)
			{
				i++;
				noise = stringToNumber<double>(argv[i]);
				flagNoise = true;
			}
			else printErrorMessage("Noise value missing; using default value of 0.01.", MSG_WARNING);
		}
		// -----------------------
		// INITIAL PARAMETERS (-p)
		// -----------------------
		else if(std::string(argv[i]) == "-p")
		{
			if(i < argc - BUSYFIT_FREE_PARAM)
			{
				for(int j = 0; j < BUSYFIT_FREE_PARAM; j++)
				{
					std::string parStr = argv[i + j + 1];
					init_par[j] = stringToNumber<double>(parStr);
					if(parStr.substr(parStr.size() - 1, 1) == "f") mask[j] = 0;
				}
				
				initialGuess = true;
				
				i += BUSYFIT_FREE_PARAM;
			}
			else printErrorMessage("Incomplete list of initial parameters specified.", MSG_WARNING);
		}
		// ---------------------
		// NO PLOTTING (-noplot)
		// ---------------------
		else if(std::string(argv[i]) == "-noplot")
		{
			noPlot = true;
		}
		// ---------------------
		// RELAX CHECKS (-relax)
		// ---------------------
		else if(std::string(argv[i]) == "-relax")
		{
			relax = true;
		}
		// ------------------
		// UNCERTAINTIES (-u)
		// ------------------
		else if(std::string(argv[i]) == "-u")
		{
			if(i < argc - 1)
			{
				i++;
				iterations = stringToNumber<int>(argv[i]);
				
				if(iterations < 1000)
				{
					iterations = 1000;
					printErrorMessage("Too few iterations; adjusting value to 1000.", MSG_WARNING);
				}
			}
			else
			{
				iterations = 10000;
				printErrorMessage("Iterations missing; using default value of 10000.", MSG_WARNING);
			}
		}
		// ------------
		// VERBOSE (-v)
		// ------------
		else if(std::string(argv[i]) == "-v")
		{
			verbose = true;
		}
		// --------------------
		// SPECTRAL WINDOW (-w)
		// --------------------
		else if(std::string(argv[i]) == "-w")
		{
			if(i < argc - 1)
			{
				// First read lower window boundary:
				i++;
				windowLower = stringToNumber<size_t>(argv[i]);
				
				// Then read upper window boundary:
				if(i < argc - 1)
				{
					i++;
					windowUpper = stringToNumber<size_t>(argv[i]);
				}
				else
				{
					printErrorMessage("Upper window boundary missing; using default value.", MSG_WARNING);
				}
			}
			else printErrorMessage("Window boundaries missing; using default values.", MSG_WARNING);
		}
		// ---------
		// FILE NAME
		// ---------
		else
		{
			filename = argv[i];
		}
	}// End command line arguments
	
	std::cout << "\nRunning BusyFit:\n----------------\n";
	
	// Read in data from input file:
	if(filename.empty())
	{
		printErrorMessage("No input file name specified.", MSG_ERROR);
		return 1;
	}
	
	std::ifstream fp(filename.c_str(), std::ifstream::in);
	
	if(fp)
	{
		std::string line;
		std::vector <std::string> tokens;
		std::string firstChar;
		
		double addNoise = 0.0;
		
		while(fp.good())
		{
			getline(fp, line);               // Read a line
			stringTrim(line);                // Trim spaces
			stringTok(line, tokens, " \t");  // Tokenise line
			
			if(static_cast<int>(tokens.size()) > column and static_cast<int>(tokens.size()) > columnSpec and column >= 0 and columnSpec >= 0)
			{
				firstChar = (tokens[0]).substr(0, 1);
				
				if(firstChar != "#" and firstChar != "/" and firstChar != "$" and firstChar != "%")
				{
					// Add some extra noise, if desired (normally disabled):
					if(false)
					{
						//double r1 = static_cast<double>(rand()) / static_cast<double>(RAND_MAX);
						//double r2 = static_cast<double>(rand()) / static_cast<double>(RAND_MAX);
						noise    = 0.05;
						addNoise = randGauss(noise);
					}
					
					frequency.push_back(stringToNumber<double>(tokens[columnSpec]));
					spectrum.push_back(stringToNumber<double>(tokens[column]) + addNoise);         // addNoise = 0 unless for testing purposes.
				}
			}
			
			tokens.clear();
		}
		
		fp.close();
	}
	else
	{
		printErrorMessage("Failed to open input file.", MSG_ERROR);
		return 1;
	}
	
	if(spectrum.size() == 0)
	{
		printErrorMessage("No valid data found in input file.", MSG_ERROR);
		return 1;
	}
	
	// Shrink spectrum to spectral window if requested:
	if(not windowUpper) windowUpper = frequency.size() - 1;
	if(windowLower >= frequency.size() or windowUpper >= frequency.size() or windowUpper - windowLower < BUSYFIT_FREE_PARAM - 1)
	{
		printErrorMessage("Spectral window too small or outside range. Using full spectrum.", MSG_WARNING);
		windowLower = 0;
		windowUpper = frequency.size() - 1;
	}
	
	std::cout << "Channel range to be fitted: " << windowLower << "-" << windowUpper << "\n";
	
	if(windowLower)
	{
		frequency.erase(frequency.begin(), frequency.begin() + windowLower);
		spectrum.erase(spectrum.begin(), spectrum.begin() + windowLower);
	}
	if(windowUpper)
	{
		frequency.erase(frequency.begin() + windowUpper - windowLower + 1, frequency.end());
		spectrum.erase(spectrum.begin() + windowUpper - windowLower + 1, spectrum.end());
	}
	
	// Automatically determine rms noise level if requested:
	if(noise <= 0.0)
	{
		std::vector <double> spectrumCopy = spectrum;                // Create copy of spectrum.
		std::sort(spectrumCopy.begin(), spectrumCopy.end());         // Sort the copy into ascending order.
		spectrumCopy.resize(spectrumCopy.size() / 3);                // Extract the lowest 1/3 of the values.
		
		if(spectrumCopy.size() < 10)
		{
			printErrorMessage("Failed to measure rms. Baseline regions too small.", MSG_ERROR);
			return 1;
		}
		
		double mad = statMad(spectrumCopy);                          // Determine MAD of lowest 1/3 values.
		spectrumCopy.clear();                                        // Clear the copy of the spectrum again.
		
		for(size_t i = 0; i < spectrum.size(); i++)
		{
			// Copy all values smaller than 5 × MAD (≈ 3.4 × rms) into copy of spectrum:
			if(fabs(spectrum[i]) < 5.0 * mad) spectrumCopy.push_back(spectrum[i]);
		}
		
		if(spectrumCopy.size() < 10)
		{
			printErrorMessage("Failed to measure rms. Baseline regions too small.", MSG_ERROR);
			return 1;
		}
		
		noise = 1.4826 * statMad(spectrumCopy);  // NOTE: The conversion factor of 1.4826 between MAD and rms is for Gaussian noise only.
		//noise = statStdDev(spectrumCopy);
		
		std::cout << "\nResult of automatic noise measurement:\n";
		std::cout << "  Channels used: " << spectrumCopy.size() << " of " << spectrum.size() << " (" << std::setprecision(3) << 100.0 * static_cast<double>(spectrumCopy.size()) / static_cast<double>(spectrum.size()) << "%)\n";
		std::cout << "  Measured rms:  " << std::setprecision(3) << noise << std::endl;
	}
	
	// Create vector of uncertainties:
	std::vector <double> uncertainties(spectrum.size(), noise);
	
	// Declare some more variables:
	double parValues[BUSYFIT_FREE_PARAM];
	double parUncert[BUSYFIT_FREE_PARAM];
	size_t dof     = 0;
	double chi2    = 0.0;
	double chi2dof = 0.0;
	double posX    = 0.0;
	double w50     = 0.0;
	double w20     = 0.0;
	double Fint    = 0.0;
	double Fpeak   = 0.0;
	double stddev_posX  = 0.0;
	double stddev_w50   = 0.0;
	double stddev_w20   = 0.0;
	double stddev_Fint  = 0.0;
	double stddev_Fpeak = 0.0;
	int    success = 0;
	
	// Run the fitting procedure:
	BusyFit busyFit;
	busyFit.setup(spectrum.size(), &spectrum[0], &uncertainties[0], noPlot, relax, verbose, windowLower);
	busyFit.setFreeParameters(mask);
	
	if(initialGuess == true) success = busyFit.fit(init_par[0], init_par[1], init_par[2], init_par[3], init_par[4], init_par[5], init_par[6], init_par[7]);
	else success = busyFit.fit();
	
	busyFit.getResult(&parValues[0], &parUncert[0], chi2, chi2dof, dof);
	busyFit.getParametersUncertainties(posX, w50, w20, Fpeak, Fint, stddev_posX, stddev_w50, stddev_w20, stddev_Fpeak, stddev_Fint, iterations);
	
	// Print parameters in original channel units:
	std::cout << "Physical parameters in channels (zero-based):\n---------------------------------------------\n";
	
	std::cout << std::setprecision(4) << "Pos  =\t" << (posX + static_cast<double>(windowLower))  << "\t±   " << stddev_posX  << "\t [chan]\n";
	std::cout << std::setprecision(4) << "w50  =\t" << w50   << "\t±   " << stddev_w50   << "\t [chan]\n";
	std::cout << std::setprecision(4) << "w20  =\t" << w20   << "\t±   " << stddev_w20   << "\t [chan]\n";
	std::cout << std::setprecision(4) << "Fpk  =\t" << Fpeak << "\t±   " << stddev_Fpeak << "\t [" << fluxUnit << "]\n";
	std::cout << std::setprecision(4) << "Fint =\t" << Fint  << "\t±   " << stddev_Fint  << "\t [" << fluxUnit << " × chan]\n" << std::endl;
	
	// Convert spectral axis from channels to spectral units (if specified):
	std::cout << "Physical parameters in spectral units (as read from file):\n----------------------------------------------------------\n";
	
	if(static_cast<size_t>(posX + 1.0) < frequency.size())
	{
		specInc = frequency[static_cast<size_t>(posX + 1.0)] - frequency[static_cast<size_t>(posX)];
	}
	else
	{
		specInc = frequency[static_cast<size_t>(posX)] - frequency[static_cast<size_t>(posX - 1.0)];
	}
	
	posX         = frequency[static_cast<size_t>(posX)] + specInc * (posX - floor(posX));
	w50         *= fabs(specInc);
	w20         *= fabs(specInc);
	Fint        *= fabs(specInc);
	stddev_posX *= fabs(specInc);
	stddev_w50  *= fabs(specInc);
	stddev_w20  *= fabs(specInc);
	stddev_Fint *= fabs(specInc);
	
	std::cout << std::setprecision(4) << "Pos  =\t" << posX  << "\t±   " << stddev_posX  << "\t [" << specUnit << "]\n";
	std::cout << std::setprecision(4) << "w50  =\t" << w50   << "\t±   " << stddev_w50   << "\t [" << specUnit << "]\n";
	std::cout << std::setprecision(4) << "w20  =\t" << w20   << "\t±   " << stddev_w20   << "\t [" << specUnit << "]\n";
	std::cout << std::setprecision(4) << "Fpk  =\t" << Fpeak << "\t±   " << stddev_Fpeak << "\t [" << fluxUnit << "]\n";
	std::cout << std::setprecision(4) << "Fint =\t" << Fint  << "\t±   " << stddev_Fint  << "\t [" << fluxUnit << " × " << specUnit << "]\n" << std::endl;
	
	if(flagNoise == false) printErrorMessage("No noise level specified; all uncertainties are arbitrary!\n         Use option '-n rms' to set baseline noise level.", MSG_WARNING);
	else printErrorMessage("Parameter uncertainties are purely statistical and may be inaccurate.", MSG_INFO);
	
	// Append results to output file:
	bool historyFileExists = false;
	std::ifstream tmpFile("busyfit_history.txt");
	
	if(tmpFile.is_open())
	{
		historyFileExists = true;
		tmpFile.close();
	}
	
	std::ofstream outputFile;
	outputFile.open("busyfit_history.txt", std::ios::app);
	
	if(outputFile.is_open())
	{
		if(historyFileExists == false)
		{
			outputFile << "# This is the local BusyFit history file in which all fit results are stored.\n";
			outputFile << "# In order to delete the history, e.g. prior to a batch run of BusyFit on a\n";
			outputFile << "# large number of spectra, simply delete the file; it will be recreated \n";
			outputFile << "# during the next launch of the software. Alternatively, you may choose to \n";
			outputFile << "# remove individual entries from the file. The results of new fits will be \n";
			outputFile << "# appended at the end of the file.\n#\n";
			outputFile << "# The parameters in the respective columns (separated by tabs) are:\n#\n";
			outputFile << "# Filename\tSuccess\tNchan\tdof\tchi^2\tchi^2/dof\trms\tA\tdA\tB_1\tdB_1\tB_2\tdB_2\tC\tdC\tXE_0\tdXE_0\tXP_0\tdXP_0\tW\tdW\tN\tdN\tX\tdX\tW_50\tdW_50\tW_20\tdW_20\tF_peak\tdF_peak\tF_int\tdF_int\n#\n";
		}
		
		size_t subStrPos = filename.find_last_of("/");
		
		if(subStrPos < filename.size() - 1) filename = filename.substr(subStrPos + 1);
		
		outputFile << filename << '\t' << success << '\t' << spectrum.size() << '\t' << dof << '\t' << chi2 << '\t' << chi2dof << '\t' << noise;
		
		for(size_t i = 0; i < BUSYFIT_FREE_PARAM; i++) outputFile << '\t' << parValues[i] << '\t' << parUncert[i];
		
		outputFile << '\t' << posX << '\t' << stddev_posX << '\t' << w50 << '\t' << stddev_w50 << '\t' << w20 << '\t' << stddev_w20 << '\t' << Fpeak << '\t' << stddev_Fpeak << '\t' << Fint << '\t' << stddev_Fint << '\n';
		
		outputFile.close();
	}
	else
	{
		printErrorMessage("Failed to write output file.", MSG_ERROR);
		return 1;
	}
	
	return 0;
}



// ------------------------------ //
// Function to print help message //
// ------------------------------ //

int printHelpMessage()
{
	std::string message;
	
	message.append("\n");
	message.append("                            \033[35m┓           \033[36m ┏━┓ ┓ \033[0m\n");
	message.append("                            \033[35m┣━┓┓ ┓┏━┓┓ ┓\033[36m╺╋━┓╺╋━\033[0m\n");
	message.append("                            \033[35m┃ ┃┃ ┃┗━┓┗━┫\033[36m ┃ ┃ ┃ \033[0m\n");
	message.append("                            \033[35m┻━┛┗━┻┗━┛┗━┛\033[36m ┻ ┻ ┗━\033[0m\n");
	message.append("                            \033[1;4mBusyFit User Manual\033[0m\n");
	message.append("\n");
	message.append("                               \033[0mVersion 0.3.2\033[0m\n");
	message.append("\n");
	message.append("  \033[3mBusyFit is a stand-alone software for fitting the so-called Busy Function\033[0m\n");
	message.append("  \033[3m(Westmeier et al. 2014, MNRAS, 438, 1176) to the integrated spectrum of a\033[0m\n");
	message.append("  \033[3mgalaxy for the purpose of parameterisation. It is entirely written in C++\033[0m\n");
	message.append("  \033[3mand licensed under the GNU General Public License.\033[0m\n");
	message.append("\n");
	message.append("\n");
	message.append("\033[1mUsage:\033[0m busyfit [OPTIONS] FILE\n");
	message.append("\n");
	message.append("  There is  no specific order  for the options  or file name.  Note that  the\n");
	message.append("  spectrum must be stored on a \033[1;4mREGULAR GRID\033[0m and in the \033[1;4mCORRECT ORDER\033[0m, without\n");
	message.append("  any gaps or jumps,  as otherwise all results will be wrong and meaningless!\n");
	message.append("  Redshift-dependent effects will \033[1;4mNOT\033[0m be taken into account.\n");
	message.append("\n");
	message.append("\n");
	message.append("\033[1mAvailable options\033[0m\n");
	message.append("\n");
	message.append("  -c C1 [C2]   Column numbers of spectral and flux values  in input file.\n");
	message.append("               Default columns are  C1 = 1  and  C2 = 2.  Note that C2 is\n");
	message.append("               optional and will default to C2 = C1 + 1 if unspecified.\n");
	message.append("\n");
	message.append("  -h, -help    Print this detailed user manual and exit.\n");
	message.append("\n");
	message.append("  -n rms       The rms noise level of the spectrum.  The default is 0.01.\n");
	message.append("               It is  important  to provide  an accurate  estimate of the\n");
	message.append("               noise,  as otherwise  the calculated  uncertainties and χ²\n");
	message.append("               values of the fit will be arbitrary and meaningless.  If a\n");
	message.append("               value of ≤ 0  is specified,  BusyFit will attempt to auto-\n");
	message.append("               matically  determine  the rms noise level  using emission-\n");
	message.append("               free regions of the spectrum.\n");
	message.append("\n");
	message.append("  -noplot      Do not plot  the fit results  on the screen.  Useful  when\n");
	message.append("               used in batch mode or when gnuplot is unavailable.\n");
	message.append("\n");
	message.append("  -p P1 … P8   Initial estimates  of the free  parameters,  A, B₁, B₂, C,\n");
	message.append("               XE₀, XP₀, W and N, separated by spaces.  If not specified,\n");
	message.append("               initial parameters will be  automatically determined.  Ad-\n");
	message.append("               ding ‘f’ at the end will fix a parameter, e.g. ‘2.5f’.\n");
	message.append("\n");
	message.append("  -relax       Do not check for negative or shifted polynomial component.\n");
	message.append("               If this option is set, the fit will not be repeated if the\n");
	message.append("               polynomial is found to be either inverted or significantly\n");
	message.append("               shifted relative to the error functions.\n");
	message.append("\n");
	message.append("  -u iter      Determine uncertainties for observational parameters (such\n");
	message.append("               as flux, line width, etc.),  using  a Monte-Carlo  method.\n");
	message.append("               This will generate ‘iter’ realisations  of the fitted Busy\n");
	message.append("               Function and can therefore be very slow. Default is to use\n");
	message.append("               a significantly faster error propagation method to compute\n");
	message.append("               uncertainties.\n");
	message.append("\n");
	message.append("  -v           Verbose mode, showing progress information during fitting.\n");
	message.append("\n");
	message.append("  -w C1 C2     Define channel range to be fitted (zero-based). Default is\n");
	message.append("               to fit all channels.  Set to 0 for first/last channel, re-\n");
	message.append("               spectively.\n");
	message.append("\n");
	message.append("\n");
	message.append("\033[1mIntroduction\033[0m\n");
	message.append("\n");
	message.append("  The BusyFit programme  implements a Levenberg–Marquardt χ² minimisation al-\n");
	message.append("  gorithm to fit an asymmetric version of the \033[4mBusy Function\033[0m to the integrated\n");
	message.append("  HI spectrum of a galaxy. The fitted Busy Function is of the form\n");
	message.append("\n");
	message.append("      B(x) = (A/4) × [erf(B₁ (W + x − XE₀)) + 1]\n");
	message.append("           × [erf(B₂ (W − x + XE₀)) + 1] × [C |x − XP₀|ᴺ + 1]\n");
	message.append("\n");
	message.append("  with the free parameters A, B₁, B₂, C, XE₀, XP₀, W and N.  This form of the\n");
	message.append("  Busy Function is well suited to fitting a wide range of symmetric and asym-\n");
	message.append("  metric double-horn profiles commonly found in galaxies.  In addition to the\n");
	message.append("  fit parameters, the programme will also measure commonly used observational\n");
	message.append("  parameters, including w₅₀ and w₂₀ line widths  as well as peak flux and in-\n");
	message.append("  tegrated flux of the galaxy.\n");
	message.append("\n");
	message.append("  If not directed otherwise, the programme will first attempt to fit the full\n");
	message.append("  version of the Busy Function. If the initial fit results in a negative amp-\n");
	message.append("  litude of the polynomial component (C < 0)  or a significant shift with re-\n");
	message.append("  spect to the position of the error functions (|XP₀ − XE₀| > 2 W₀,  where W₀\n");
	message.append("  is the initial  estimate of the  profile width),  the fit will be discarded\n");
	message.append("  and automatically repeated  without  the polynomial component,  i.e. C, XP₀\n");
	message.append("  and N  will be fixed to zero and the fit repeated with only five free para-\n");
	message.append("  meters.  This behaviour can be changed  with the ‘-relax’ option,  in which\n");
	message.append("  case negative or shifted polynomial solutions  will be accepted without re-\n");
	message.append("  peating the fit.\n");
	message.append("\n");
	message.append("\n");
	message.append("\033[1mInput data\033[0m\n");
	message.append("\n");
	message.append("  The HI spectrum  will be read  from a text file,  and the user will have to\n");
	message.append("  specify both the input file name as well as the numbers of the columns that\n");
	message.append("  contain the spectral and flux values of the spectrum. For example,\n"); 
	message.append("\n");
	message.append("      busyfit -c 3 4 spectrum.dat\n");
	message.append("\n");
	message.append("  would read the third and fourth columns from the file spectrum.dat  and as-\n");
	message.append("  sume that the third column contains the spectral axis and the fourth column\n");
	message.append("  contains the flux values.  The parser assumes that the tabulated values are\n");
	message.append("  separated by either a space or a tab  and will ignore lines starting with a\n");
	message.append("  comment character (#, /, $, and %)  as well as empty lines  containing only\n");
	message.append("  spaces.\n");
	message.append("\n");
	message.append("  ╔═════════════════════════════════════════════════════════════════════════╗\n");
	message.append("  ║ Note that  the programme  makes the  implicit assumption  that the flux ║\n");
	message.append("  ║ values in the input file are listed in the \033[1;4mCORRECT ORDER\033[0m and on a \033[1;4mREGU-\033[0m ║\n");
	message.append("  ║ \033[1;4mLAR GRID\033[0m in  either frequency  or velocity,  without any gaps or jumps. ║\n");
	message.append("  ║ All results will otherwise be wrong and entirely meaningless!           ║\n");
	message.append("  ╚═════════════════════════════════════════════════════════════════════════╝\n");
	message.append("\n");
	message.append("  Redshift-dependent correction factors will  \033[4mNOT\033[0m  be taken into account when\n");
	message.append("  calculating physical parameters of the fitted lines,  e.g. velocity widths.\n");
	message.append("  If in doubt,  provide your spectra  in channels  rather than physical units\n");
	message.append("  and then manually apply the correct velocity conversion factors afterwards.\n");
	message.append("\n");
	message.append("  BusyFit will fit the entire spectrum by default. This could create problems\n");
	message.append("  in cases where multiple signals are present,  only one of which is meant to\n");
	message.append("  be fitted. In this case, a spectral channel range can be specified with the\n");
	message.append("  ‘-w C1 C2’ option,  where C1 and C2 are the lower  and upper boundary,  re-\n");
	message.append("  spectively. Either one can be set to 0 to signify the first/last channel of\n");
	message.append("  the spectrum.\n");
	message.append("\n");
	message.append("  In order for the uncertainties estimated by the algorithm to be meaningful,\n");
	message.append("  it is necessary to specify  the \033[4mrms noise level\033[0m  of the spectrum  using the\n");
	message.append("  ‘-n’ option  followed by  a numerical value  that specifies  the rms in the\n");
	message.append("  native flux units  of the data.  If a value of ≤ 0  is given,  then BusyFit\n");
	message.append("  will attempt  to automatically measure  the rms  in regions of the spectrum\n");
	message.append("  that are deemed  free of emission.  Note that this will  only work if there\n");
	message.append("  are sufficiently large regions of baseline  around the actual spectral line\n");
	message.append("  (at least 1/3 of the spectrum must be pure noise).\n");
	message.append("\n");
	message.append("\n");
	message.append("\033[1mInitial parameter estimates\033[0m\n");
	message.append("\n");
	message.append("  Without any further specifications,  BusyFit will automatically determine a\n");
	message.append("  set of initial estimates for the free parameters.  In many cases, this will\n");
	message.append("  be sufficient for the fit to converge. Alternatively, initial estimates can\n");
	message.append("  be provided with the ‘-p’ option in the order A, B₁, B₂, C, XE₀, XP₀, W and\n");
	message.append("  N. For example,\n");
	message.append("\n");
	message.append("      busyfit -c 3 4 spectrum.dat -p 1 0.5 0.5 0.1 20 20 5 2\n");
	message.append("\n");
	message.append("  will set initial parameter estimates of A = 1, B₁ = B₂ = 0.5, C = 0.1, etc.\n");
	message.append("  In addition, parameters can be fixed by adding the letter ‘f’ to the value.\n");
	message.append("  For example,  ‘2.5f’  will set the respective parameter to a fixed value of\n");
	message.append("  2.5 and exclude it from the fitting procedure.\n");
	message.append("\n");
	message.append("  Note that the \033[4munits\033[0m used for the free parameters  are the native flux units\n");
	message.append("  of the spectrum and spectral channels, not native spectral units! The first\n");
	message.append("  channel of the spectrum is channel number 0.\n");
	message.append("\n");
	message.append("\n");
	message.append("\033[1mPlotting the results\033[0m\n");
	message.append("\n");
	message.append("  The programme will automatically  plot the result of the fit on the screen,\n");
	message.append("  using \033[4mgnuplot\033[0m.  This will only work if gnuplot is installed and can be pre-\n");
	message.append("  vented by using the ‘-noplot’ option. The plot will show the original spec-\n");
	message.append("  trum, the fitted Busy Function, and the resulting residuals. Note that each\n");
	message.append("  call of BusyFit will open a \033[4mseparate gnuplot window\033[0m, allowing results to be\n");
	message.append("  compared to previous fits. Hence, when automatically fitting a large number\n");
	message.append("  of spectra,  it is important to use the  ‘-noplot’  option,  as otherwise a\n");
	message.append("  large number of gnuplot windows  would be crated,  potentially crashing the\n");
	message.append("  window manager.\n");
	message.append("\n");
	message.append("\n");
	message.append("\033[1mOutput files\033[0m\n");
	message.append("\n");
	message.append("  Every successful run of BusyFit  will create \033[4mthree output files\033[0m.  The first\n");
	message.append("  two, ‘busyfit_output_spectrum.txt’ and  ‘busyfit_output_fit.txt’ store  the\n");
	message.append("  original spectrum + residuals and a numerical copy of the fitted Busy Func-\n");
	message.append("  tion, respectively.  These files are used  to plot  the fit results  on the\n");
	message.append("  screen using gnuplot.  The previous contents of the files  will be replaced\n");
	message.append("  whenever the programme is run again.\n");
	message.append("\n");
	message.append("  The third output file, ‘busyfit_history.txt’,  will store a copy of the fit\n");
	message.append("  results. Previous results will not be deleted, but instead new results from\n");
	message.append("  the most recent run  of the programme  will be appended  at the end  of the\n");
	message.append("  file.  This behaviour is useful  when running BusyFit repeatedly, e.g. on a\n");
	message.append("  large number of spectra.  Before the run, the output file can be deleted or\n");
	message.append("  cleared,  and the results of all runs will then be available for processing\n");
	message.append("  after the last run has finished.\n");
	message.append("\n");
	message.append("  The file ‘busyfit_history.txt’ will store each set of results in a separate\n");
	message.append("  line, and each line contains the following entries,  separated by tabs,  in\n");
	message.append("  the given order:\n");
	message.append("\n");
	message.append("   • Input file name\n");
	message.append("   • Success flag (0, 1, or 2)\n");
	message.append("   • Number of data samples and degrees of freedom\n");
	message.append("   • Original and reduced χ²\n");
	message.append("   • Baseline rms\n");
	message.append("   • Values and uncertainties of parameters A, B₁, B₂, C, XE₀, XP₀, W and N\n");
	message.append("   • Values and uncertainties of centroid, w₅₀, w₂₀, peak flux, and integra-\n");
	message.append("     ted flux\n");
	message.append("\n");
	message.append("  The success flag will tell  whether the fit was fully successful (0),  con-\n");
	message.append("  verged successfully  but with some parameters  negative (1),  or completely\n");
	message.append("  failed to converge (2).  The Busy Function parameters are provided in units\n");
	message.append("  of spectral channels,  not the physical units of the spectrum.  Counting of\n");
	message.append("  channels starts at 0, not 1.  The values of centroid,  w₅₀ and w₂₀  will be\n");
	message.append("  given in physical spectral units,  the peak flux  will be given in the ori-\n");
	message.append("  ginal flux units of the input file, and the integrated flux  is measured in\n");
	message.append("  flux units × spectral units.\n");
	message.append("\n");
	message.append("\n");
	message.append("\033[1mUncertainties of observational parameters\033[0m\n");
	message.append("\n");
	message.append("  By default,  uncertainties of the observational parameters  will be derived\n");
	message.append("  using an \033[4merror propagation\033[0m method. This method uses the standard error pro-\n");
	message.append("  pagation law under the assumption that uncertainties  in the  Busy Function\n");
	message.append("  parameters  linearly propagate  into uncertainties  in the derived observa-\n");
	message.append("  tional parameters.  While this method is very fast, its accuracy can be li-\n");
	message.append("  mited in situations where non-linear effects become relevant.\n");
	message.append("\n");
	message.append("  If accuracy is essential,  an alternative \033[4mMonte-Carlo method\033[0m can be invoked\n");
	message.append("  with the ‘-u iter’ option. In this case, BusyFit will generate ‘iter’ real-\n");
	message.append("  isations  of the fitted  Busy Function  by randomly varying  the function’s\n");
	message.append("  free  parameters  and recalculating  the observational  parameters  in each\n");
	message.append("  iteration.  The uncertainty of each parameter is then taken as the standard\n");
	message.append("  deviation about the mean across all iterations.  A sensible number of iter-\n");
	message.append("  ations is 10,000,  which will determine  uncertainties  with an accuracy of\n");
	message.append("  typically a few percent.  The disadvantage  of the  Monte-Carlo approach is\n");
	message.append("  that it can be  very slow if a large number of spectra  needs to be fitted,\n");
	message.append("  and the default error propagation method may be more suitable in this case.\n");
	message.append("  It is also  strongly advisable  to fix  the polynomial  exponent, N, of the\n");
	message.append("  Busy Function in this case,  as this parameter is non-Gaussian, and its in-\n");
	message.append("  clusion would result in unrealistically large uncertainties  due to the al-\n");
	message.append("  gorithm operating outside its intended linear regime.\n");
	message.append("\n");
	message.append("  ╔═════════════════════════════════════════════════════════════════════════╗\n");
	message.append("  ║ \033[1mWarning 1:\033[0m For the error analysis method to work,  it is absolutely es- ║\n");
	message.append("  ║            sential that the correct  baseline \033[4mnoise level\033[0m  be specified ║\n");
	message.append("  ║            with the ‘-n rms’ option, as otherwise the calculated uncer- ║\n");
	message.append("  ║            tainties of all parameters will be entirely arbitrary and    ║\n");
	message.append("  ║            meaningless.                                                 ║\n");
	message.append("  ╚═════════════════════════════════════════════════════════════════════════╝\n");
	message.append("\n");
	message.append("  ╔═════════════════════════════════════════════════════════════════════════╗\n");
	message.append("  ║ \033[1mWarning 2:\033[0m The error analysis methods implemented in BusyFit  are based ║\n");
	message.append("  ║            on several  \033[4massumptions\033[0m,  some of which  will break down  in ║\n");
	message.append("  ║            certain situations, leading to unrealistic uncertainty esti- ║\n");
	message.append("  ║            mates.  In particular,  the calculated uncertainties  of the ║\n");
	message.append("  ║            derived observational parameters  will  almost certainly  be ║\n");
	message.append("  ║            wrong if the uncertainty, σ, of any of the non-additive free ║\n");
	message.append("  ║            parameters, P, of the fit  does not fulfil  the condition of ║\n");
	message.append("  ║            σ << P.                                                      ║\n");
	message.append("  ╚═════════════════════════════════════════════════════════════════════════╝\n");
	message.append("\n");
	message.append("  ╔═════════════════════════════════════════════════════════════════════════╗\n");
	message.append("  ║ \033[1mWarning 3:\033[0m Most observational data  are affected  by strong  \033[4msystematic\033[0m ║\n");
	message.append("  ║            \033[4merrors\033[0m that are usually much larger than the statistical er- ║\n");
	message.append("  ║            rors.  The statistical uncertainties  calculated  by BusyFit ║\n");
	message.append("  ║            will be entirely unrealistic and far too small in this case. ║\n");
	message.append("  ║            Instead, Monte Carlo or bootstrapping techniques will be re- ║\n");
	message.append("  ║            quired to accurately  characterise  the uncertainties of any ║\n");
	message.append("  ║            derived parameters.                                          ║\n");
	message.append("  ╚═════════════════════════════════════════════════════════════════════════╝\n");
	message.append("\n");
	message.append("\n");
	message.append("\033[1mCopyright and licence\033[0m\n");
	message.append("\n");
	message.append("  © 2019 Tobias Westmeier\n");
	message.append("\n");
	message.append("  This programme is free software:  you can redistribute it  and/or modify it\n");
	message.append("  under the terms of the \033[4mGNU General Public License\033[0m as published  by the Free\n");
	message.append("  Software Foundation,  either version 3 of the License,  or (at your option)\n");
	message.append("  any later version.\n");
	message.append("\n");
	message.append("  This programme is distributed in the hope that it will be useful, but \033[1mWITH-\033[0m\n");
	message.append("  \033[1mOUT ANY WARRANTY\033[0m;  without even the implied warranty of  \033[1mMERCHANTABILITY\033[0m or\n");
	message.append("  \033[1mFITNESS FOR A PARTICULAR PURPOSE\033[0m.  See the GNU General Public  License  for\n");
	message.append("  more details.\n");
	message.append("\n");
	message.append("  You should have  received  a copy of the  GNU General Public License  along\n");
	message.append("  with this programme. If not, see http://www.gnu.org/licenses/.\n");
	message.append("\n");
	message.append("  Contact: Tobias Westmeier\n");
	message.append("           ICRAR M468\n");
	message.append("           The University of Western Australia\n");
	message.append("           35 Stirling Highway\n");
	message.append("           Crawley WA 6009\n");
	message.append("           Australia\n");
	message.append("\n");
	message.append("  E-mail:  tobias.westmeier (at) uwa.edu.au\n");
	message.append("\n");
	message.append("\n");
	message.append("\033[1mLinks and references\033[0m\n");
	message.append("\n");
	message.append("  • \033[3mBusyFit website:\033[0m\n");
	message.append("    http://www.atnf.csiro.au/people/Tobias.Westmeier/tools_software_busyfit.php\n");
	message.append("\n");
	message.append("  • \033[3mBusy Function paper:\033[0m\n");
	message.append("    ADS   - http://adsabs.harvard.edu/abs/2014MNRAS.438.1176W\n");
	message.append("    arXiv - http://arxiv.org/abs/1311.5308\n");
	message.append("\n");
	message.append("  • \033[3mRussell Jurek’s Busy Function fitting code:\033[0m\n");
	message.append("    https://github.com/RussellJurek/busy-function-fitting\n");
	message.append("\n");
	message.append("\n");
	message.append("\n");
	message.append("\033[1;7m End of document. Press ‘Q’ to exit. \033[0m\n");
	
	std::string systemCall = "echo \"";
	systemCall.append(message);
	systemCall.append("\" | less -R");
	
	return system(systemCall.c_str());
}
