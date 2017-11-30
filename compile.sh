#!/bin/sh
#/// ____________________________________________________________________ ///
#///                                                                      ///
#/// BusyFit 0.3.1 (compile.sh) - Busy Function fitting programme         ///
#/// Copyright (C) 2017 Tobias Westmeier                                  ///
#/// ____________________________________________________________________ ///
#///                                                                      ///
#/// Address:  Tobias Westmeier                                           ///
#///           ICRAR M468                                                 ///
#///           The University of Western Australia                        ///
#///           35 Stirling Highway                                        ///
#///           Crawley WA 6009                                            ///
#///           Australia                                                  ///
#///                                                                      ///
#/// E-mail:   tobias.westmeier (at) uwa.edu.au                           ///
#/// ____________________________________________________________________ ///
#///                                                                      ///
#/// This program is free software: you can redistribute it and/or modify ///
#/// it under the terms of the GNU General Public License as published by ///
#/// the Free Software Foundation, either version 3 of the License, or    ///
#/// (at your option) any later version.                                  ///
#///                                                                      ///
#/// This program is distributed in the hope that it will be useful,      ///
#/// but WITHOUT ANY WARRANTY; without even the implied warranty of       ///
#/// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the         ///
#/// GNU General Public License for more details.                         ///
#///                                                                      ///
#/// You should have received a copy of the GNU General Public License    ///
#/// along with this program. If not, see http://www.gnu.org/licenses/.   ///
#/// ____________________________________________________________________ ///
#///                                                                      ///

echo "Compiling BusyFit..."
g++ -pedantic -Wall -O3 -c helperFunctions.cpp
g++ -pedantic -Wall -O3 -c BusyFit.cpp
g++ -pedantic -Wall -O3 -o busyfit helperFunctions.o BusyFit.o main.cpp -lgsl -lgslcblas
#rm -f BusyFit.o
#rm -f helperFunctions.o
echo "Done."
