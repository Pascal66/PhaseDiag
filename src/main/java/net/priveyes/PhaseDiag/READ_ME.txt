

                                ****************************************
                                *   ___           o                    *
                                *   |    |   |  | | |\                 *
                                *   |_   |   |  | | | \        o       *
                                *   |    |   |  | | | / |\  /| | |\ |  *
                                *   |    |__ |__| | |/  | \/ | | | \|  *
                                *                                      *
                                ****************************************

This program calculates the density and viscosity of water and CO2 as a function of pressure and temperature,
using published physical models.
The program also calculates the composition of co-existing phases in H2O-CO2-NaCl,
which includes CO2 solubility in water and brine. It relies on minimization of the energy of the system.
Error propagation from the uncertainty on the mixing parameters is carried out with the help of a Monte-Carlo simulation.
The program is written in Fortran 90, therefore it must be compiled
(for example, using gfortran from the GNU Compiler Collection).
It is organised in several different modules, which must also be compiled.
A bash script (compil.sh) is provided as an example. For Windows users, the easiest option is to use Cygwin
for scripting and compilation.

Obviously the aim of compiling sources and using the off-line version of the program is to modify it to suit your needs.
Therefore there is no need for a lengthy description of the program in its current form,
as it is described on-line with most if not all of the required equations referenced there. A quick summary is given below.

The program is called with 7 input arguments, such as:
>> FLUIDmin.exe min max nCalc var1 var2 typeCalc typeOutput

	typeOutput: 1 | 2 | 3 (integer)
		1: short (composition of water-rich and CO2-rich phases)
		2: extended (adds density and viscosity)
		3: all (adds thermodynamic parameters)

	typeCalc: 1 | 2 | 3 (integer)
		1: Temperature is varying, pressure and salinity are constant
		2: Pressure is varying, temperature and salinity are constant
		3: Salinity is varying, pressure and temperature are constant
		
	min and max: (positive real, with min < max) define the range of the varying variable (see typeCalc above)
	
	nCalc: (integer) number of iterations from min to max (see min and max above)
	
	var1 and var2: set values for the fixed variables.
		for typeCalc == 1: var1 is pressure, var2 is salinity
		for typeCalc == 2: var1 is temperature, var2 is salinity
		for typeCalc == 3: var1 is pressure, var2 is temperature

Units: 
	Pressure is in MPa
	Temperature in Celsius degrees
	Salinity in mol./kg(H2O)
	
Example:
>> FLUIDmin.exe 20 50 31 8 0 1 2
	will produce 31 calculations, evenly spaced from 20 to 50 Celsius degrees,
	at the fixed pressure of 8 MPa in H2O-CO2 (salinity of 0 mol./kg(H2O)) and output compositions + density and viscosity.
	
The output is a table containing results in columns, with one row par iteration,
in the same order as the csv export function on the website.

For questions, do not hesitate to contact the author.