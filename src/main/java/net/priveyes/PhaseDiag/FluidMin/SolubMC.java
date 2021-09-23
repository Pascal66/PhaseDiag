package net.priveyes.PhaseDiag.FluidMin;

import static net.priveyes.PhaseDiag.EVALW.EVALW.returnSigma.SIGMA;
import static net.priveyes.PhaseDiag.EVALW.EVALW.returnSigma.StDevAlpha;
import static net.priveyes.PhaseDiag.EVALW.EVALW.returnSigma.StDevW;
import static net.priveyes.PhaseDiag.EVALW.EVALW.returnSigma.varAlpha;
import static net.priveyes.PhaseDiag.EVALW.EVALW.returnSigma.varW;
import static net.priveyes.PhaseDiag.EVALW.EVALW.returnWCO2NACL.StDevwCO2NaCl0;
import static net.priveyes.PhaseDiag.EVALW.EVALW.returnWCO2NACL.WCO2NACL;
import static net.priveyes.PhaseDiag.EVALW.EVALW.returnWCO2NACL.sigwCO2NaCl0;
import static net.priveyes.PhaseDiag.EVALW.EVALW.returnWCO2NACL.wCO2NaCl0;
import static net.priveyes.PhaseDiag.EVALW.EVALW.returnWCO2NACL.wCO2NaClpm;
import static net.priveyes.PhaseDiag.EVALW.EVALW.returnWH2OCO2.WH2OCO2;
import static net.priveyes.PhaseDiag.EVALW.EVALW.returnWH2OCO2.aCO2;
import static net.priveyes.PhaseDiag.EVALW.EVALW.returnWH2OCO2.aH2O;
import static net.priveyes.PhaseDiag.EVALW.EVALW.returnWH2OCO2.wWCO2;
import static net.priveyes.PhaseDiag.EVALW.EVALW.returnWH2ONACL.StDevaNaCl0;
import static net.priveyes.PhaseDiag.EVALW.EVALW.returnWH2ONACL.StDevaNaClpm;
import static net.priveyes.PhaseDiag.EVALW.EVALW.returnWH2ONACL.StDevwWNaClo;
import static net.priveyes.PhaseDiag.EVALW.EVALW.returnWH2ONACL.StDevwWNaClpm;
import static net.priveyes.PhaseDiag.EVALW.EVALW.returnWH2ONACL.WH2ONACL;
import static net.priveyes.PhaseDiag.EVALW.EVALW.returnWH2ONACL.aNaCl0;
import static net.priveyes.PhaseDiag.EVALW.EVALW.returnWH2ONACL.aNaClpm;
import static net.priveyes.PhaseDiag.EVALW.EVALW.returnWH2ONACL.varaNaCl0;
import static net.priveyes.PhaseDiag.EVALW.EVALW.returnWH2ONACL.varaNaClpm;
import static net.priveyes.PhaseDiag.EVALW.EVALW.returnWH2ONACL.varwWNaClo;
import static net.priveyes.PhaseDiag.EVALW.EVALW.returnWH2ONACL.varwWNaClpm;
import static net.priveyes.PhaseDiag.EVALW.EVALW.returnWH2ONACL.wNaCloNaClpm;
import static net.priveyes.PhaseDiag.EVALW.EVALW.returnWH2ONACL.wWNaClo;
import static net.priveyes.PhaseDiag.EVALW.EVALW.returnWH2ONACL.wWNaClpm;
import static net.priveyes.PhaseDiag.FluidMin.Constants.Rm;
import static net.priveyes.PhaseDiag.FluidMin.Constants.mCO2;
import static net.priveyes.PhaseDiag.FluidMin.Constants.mH2O;
import static net.priveyes.PhaseDiag.FluidMin.solutionProperty.Wlist;
import static net.priveyes.PhaseDiag.FluidMin.solutionProperty.alphas;
import static net.priveyes.PhaseDiag.FluidMin.solutionProperty.nH2Oapp;
import static net.priveyes.PhaseDiag.H2OCO2.CompH2OCO2.PhaseDiagH2OCO2;
import static net.priveyes.PhaseDiag.THERMO.thermo.RTlnGamASF;
import static net.priveyes.PhaseDiag.THERMO.thermo.cacheThermo;
import static net.priveyes.PhaseDiag.THERMO.thermo.diecIAPWS97;
import static net.priveyes.PhaseDiag.Utils.csv_write;
import static net.priveyes.PhaseDiag.VISCO.visco.Bando2004;
import static net.priveyes.PhaseDiag.VISCO.visco.ratioVisco;
import static net.priveyes.PhaseDiag.VISCO.visco.viscoCO2;
import static net.priveyes.PhaseDiag.VISCO.visco.viscoWater;
import static net.priveyes.PhaseDiag.FluidMin.Minimizer.ITMAX;
import static java.lang.Double.NaN;

import android.util.Log;

import net.priveyes.PhaseDiag.CO2EOS.CO2EOS;
import net.priveyes.PhaseDiag.EVALW.EVALW.returnWCO2NACL;
import net.priveyes.PhaseDiag.IAPWS.WATEREOS;
import net.priveyes.PhaseDiag.THERMO.thermo;
import net.priveyes.PhaseDiag.neural.Net;
import net.priveyes.PhaseDiag.num_rep.VecMath;

import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.stat.StatUtils;
import org.apache.commons.math3.util.FastMath;

import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ForkJoinPool;
import java.util.concurrent.Future;
import java.util.concurrent.ThreadLocalRandom;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.concurrent.atomic.AtomicReference;

/**
 *Copyright 2015 Benoît Dubacq, CNRS.
 *The present program is distributed under the terms of the GNU General Public License.
 *http://baobab.istep.upmc.fr/model/index.php
* https://espace.curtin.edu.au/bitstream/handle/20.500.11937/8063/196580_196580.pdf?sequence=2
   *
 *  1: H2O, 2: CO2, 3: NaCl0, 4: NaCl+-
     *
*
 *For documentation, see the README file, the website and our original publication.
 *Note that some modules and fonctions used here come from sources freely available online with references cited in the online documentation of the program.
*/
public class SolubMC {
	static ForkJoinPool executor = new ForkJoinPool();
	//final static int nProcessors = Runtime.getRuntime().availableProcessors();
	//static ExecutorService executor = Executors.newWorkStealingPool(nProcessors);
	public static RealMatrix /*double[][]*/ densities = null;
	public static RealMatrix /*double[][]*/ viscosities = null;


	/**
	 * @param minCalc
	 * @param maxCalc     min and max: (positive real, with min < max) define the range of the varying variable (see typeCalc above)
	 * @param nCalc       (integer) number of iterations from min to max (see min and max above)
	 * @param var1
	 * @param var2        var1 and var2: set values for the fixed variables.
	 *                    for typeCalc == 1: var1 is pressure, var2 is salinity
	 *                    for typeCalc == 2: var1 is temperature, var2 is salinity
	 *                    for typeCalc == 3: var1 is pressure, var2 is temperature
	 * @param typeCalc:   1 | 2 | 3 (integer)
	 *                    1: Temperature is varying, pressure and salinity are constant
	 *                    2: Pressure is varying, temperature and salinity are constant
	 *                    3: Salinity is varying, pressure and temperature are constant
	 * @param typeOutput: 1 | 2 | 3 (integer)
	 *                    1: short (composition of water-rich and CO2-rich phases)
	 *                    2: extended (adds density and viscosity)
	 *                    3: all (adds thermodynamic parameters)
	 *                    Units:
	 *                    Pressure is in MPa
	 *                    Temperature in Celsius degrees
	 *                    Salinity in mol./kg(H2O)
	 */

	public static void FLUIDmin(int fast, double minCalc, double maxCalc, int nCalc, double var1, double var2, int typeCalc, int typeOutput) throws ExecutionException, InterruptedException {
		long start = System.nanoTime();

		Log.i("Scalp", "Initialising Neural...");
		Net net = new Net(new int[]{4, 5, 4}, true);
//USE MINIMIZE //! Anonymized Nelder-Mead

		//moved to thermo
		//double[] /*REAL(DP), DIMENSION(15) ::*/ propNaCl0, propNaClpm, propHalite;
		double   mu0NaCl0, mu0NaClpm, nan, arg = -1.0;
//		double   pkb, tK /*,wWCO2, aH2O, aCO2*/;
		//double   wWNaClpm, wWNaClo, wNaCloNaClpm, aNaClpm, aNaCl0, varW, varAlpha;
		double   /*wCO2NaCl0,*/ /*wCO2NaClpm,*/ ionic, ftol;
		double   DHNaClpm, sCO2inW = 0, sWinCO2 = 0;
		double   molefracCO2, fracNaClapp, xH2Oo, xNaClapp, coef, nNaCl0, nNaClpm;
		double   /*minCalc, maxCalc,*/ delta, /*var1, var2,*/ densityWater, densityCO2;
		//double   nNaCl0, nNaClpm, molalNaCl;
		double   ntotal, uncertM/*, sigwCO2NaCl0*/;
		//double   varwWNaClpm = 0, varwWNaClo = 0, varaNaClpm = 0, varaNaCl0 = 0;
		//double[][] /*REAL(DP), DIMENSION(4,4)::*/ Wlist = new double[4][4]; //MatrixUtils.createRealMatrix(4 /*+ 1*/, 4 /*+ 1*/);
		double[] /*REAL(DP), DIMENSION(100) ::*/  compoCO2 = new double[100 /*+ 1*/], xCO2 = new double[100 /*+ 1*/], nCO2 = new double[100 /*+ 1*/], XNapmc = new double[100 /*+ 1*/];
		//double[] /*REAL(DP), DIMENSION(4)::*/ alphas;
		double[] proportions;
		final double[] /*RealVector /*REAL(DP), DIMENSION(27) ::*/ optionsSPECNaCl = new double[29];
		final RealVector /*double[]*/ /*REAL(DP), DIMENSION(29) ::*/ optionsSolubCO2 = MatrixUtils.createRealVector(new double[29 /*+ 1*/]);
		int /*INTEGER(I4B) ::*/ i/*, niter*/ /*nCalc, typeOutput, typeCalc,*/ /*niter1*/;
		final int niter = ITMAX;
		double[]/*[]*/ /*REAL(DP), DIMENSION(2,1) ::*/ param = new double[2 /*+ 1*/]/*[1+1*/;
		double[] /*REAL(DP), DIMENSION(2) ::*/ y, uncerty = new double[2 /*+ 1*/];
		RealVector /*REAL(DP), DIMENSION(3) ::*/ brineEM, gasEM, compo;
		final double[] /*REAL(DP), DIMENSION(7) ::*/ randos = new double[7];
		//double /*REAL(DP), PARAMETER ::*/ /*MCO2 = 44.0095,*/ /*R = 8.314462e-3*/;
		RealMatrix /*double[][]*/ /*real(DP), dimension(:,:), ALLOCATABLE ::*/ condition;
		RealMatrix /*double[][]*/ solub;
//		RealMatrix /*double[][]*/ densities = null;
//		RealMatrix /*double[][]*/ viscosities = null;

		Map<Integer, String>/*double[][]*/ thermoString = new HashMap<>();
		double[][] thermoVal = new double[25][25];

		double[] /*real(DP), dimension(:), ALLOCATABLE ::*/ actCO2 = null, actH2O = null, actNaCl0 = null, actNaClpm = null;
		AtomicReference<Double> wWCO2MC = new AtomicReference<>((double) 0);
		AtomicReference<Double> aCO2MC = new AtomicReference<>(0.0);
		AtomicReference<Double> wWNaClpmMC = new AtomicReference<>(0.0);
		AtomicReference<Double> wWNaCloMC = new AtomicReference<>(0.0);
		AtomicReference<Double> aNaClpmMC = new AtomicReference<>(0.0); //, compoCO2mean, nCO2mean, xCO2mean;
		AtomicReference<Double> aNaCl0MC = new AtomicReference<>(0.0);
		AtomicReference<Double> wCO2NaCl0MC = new AtomicReference<>(0.0);
		//double  /*real(DP) ::*/ meanactH2O, stdevactH2O, meanactCO2, stdevactCO2, meanactNaClpm, stdevactNaClpm, meanactNaCl0, stdevactNaCl0;
		//double   StDevW, StDevwWNaClo, StDevwWNaClpm, StDevwCO2NaCl0, StDevAlpha, StDevaNaCl0, StDevaNaClpm;
		double   kk1, kk2, kk3, kk4, pb, tC;
		boolean /*LOGICAL ::*/ valid;
//! for command line arg
		int /*integer(I4B) ::*/ num_args;
		int ix;
		AtomicInteger mc = new AtomicInteger();
//String /*character(len=12), dimension(7) ::*/ args;

//		nan = Math.sqrt(arg);

//! reading command line arguments
// num_args = command_argument_count();
//		if (/*num_args*/args.length != /*/=-*/ 7) {
//			throw new Error("input problem//! need 7 input arguments: minCalc maxCalc nCalc var1 var2 typeCalc typeOutput");
//		}

//		for (ix = 1; ix <= args.length; ix++) { //do ix = 1, num_args
		// /*call*/ get_command_argument(ix,args(ix));
		//! write(*,*) "arg. ", ix, " is ", args(ix)
//		}//end do

//		minCalc = Double.parseDouble(args[1 - 1]); //read (args[1],*) minCalc;
//		maxCalc = Double.parseDouble(args[2 - 1]); //read (args[2],*) maxCalc;
//		nCalc = Integer.parseInt(args[3 - 1]); //read (args[3],'(I12)') nCalc;
//		var1 = Double.parseDouble(args[4 - 1]); //read (args[4],*) var1;
//		var2 = Double.parseDouble(args[5 - 1]); //read (args[5],*) var2;
//		typeCalc = Integer.parseInt(args[6 - 1]); //read (args[6],'(I12)') typeCalc;
//		typeOutput = Integer.parseInt(args[7 - 1]); //read (args[7],'(I12)') typeOutput;

		//allocate (conditions(nCalc, 3), solub(nCalc,6));
		condition = MatrixUtils.createRealMatrix(nCalc /*+ 1*/, 3 /*+ 1*/)/*[nCalc][3+1]*/;
		solub = MatrixUtils.createRealMatrix(nCalc /*+ 1*/, 6 /*+ 1*/);
		if (typeOutput != /*/=*/ 1) { 
			//allocate (densities(nCalc, 2), viscosities(nCalc, 4));
			densities = MatrixUtils.createRealMatrix(nCalc /*+ 1*/, 3 /*+ 1*/);
			viscosities = MatrixUtils.createRealMatrix(nCalc /*+ 1*/, 4 /*+ 1*/);
			if (typeOutput == 3) { 
				//allocate (thermoVal(nCalc,23), actH2O(100), actNaCl0(100), actNaClpm(100), actCO2(100));
				//thermoVal = MatrixUtils.createRealMatrix(nCalc /*+ 1*/, 23 /*+ 1*/);
				actH2O = new double[100 /*+ 1*/];
				actNaCl0 = new double[100 /*+ 1*/];
				actNaClpm = new double[100 /*+ 1*/];
				actCO2 = new double[100 /*+ 1*/];
				//!prepare MonteCarlo: allocate arrays
			}//end if
		}//end if
//! formatting conditions

		delta = (maxCalc - minCalc) / (nCalc /*- 1.0*/);

		RealVector vPressure = condition.getColumnVector(0);
		RealVector vTemperature = condition.getColumnVector(1);
		RealVector vSalinity = condition.getColumnVector(2);

		if (typeCalc == 1) { //then //! varying T
			vPressure.set(var1 * 10);
//			conditions.getColumnVector(1-1).set(/*conditions(:, 1) =*/ var1 * 10.0); //! MPa -> bar
			vSalinity.set(var2);
//			conditions.getColumnVector(3-1).set(/*conditions(:, 3) =*/ var2);
			condition.setEntry(0, 1,/*;[1][2] =*/ minCalc);
			condition.setEntry(nCalc - 1, 1, /*);[nCalc][2] =*/ maxCalc);
			for (i = 1; i < nCalc - 1; i++) { //do i = 2, nCalc-1
				condition.setEntry(i, 1, /*[i][2] =*/ minCalc + delta * (i /*- 1.0*/));
			} //end do//! i

		} else if (typeCalc == 2) { //then //! varying P
			vTemperature.set(var1);
//			conditions.getColumnVector(2-1).set(/**conditions(:, 2) =*/var1);
			vSalinity.set(var2);
//			conditions.getColumnVector(3-1).set(/*/*conditions(:, 3) =*/ var2);
			condition.setEntry(0, 0, /*);[1][1] = */minCalc * 10.0); //! MPa -> bar
			condition.setEntry(nCalc - 1, 0, /*);[nCalc][1] =*/ maxCalc * 10.0); //! MPa -> bar
			for (i = 1; i < nCalc - 1; i++) { //do i = 2, nCalc-1
				condition.setEntry(i, 0, /*);[i][1] =*/ (minCalc + delta * (i /*- 1.0*/)) * 10.0); //! MPa -> bar
			} //end do//! i

		} else {//! varying S
			vTemperature.set(var2);
//			RealVector svar1 = conditions.getColumnVector(2-1); svar1.set(var1);
//			conditions.setColumnVector(2-1, svar1);
			vPressure.set(var1 * 10);
			//conditions.getColumnVector(2-1).set(/*conditions(:, 2) =*/ var1);
//			RealVector svar2 = conditions.getColumnVector(1-1); svar2.set(var2 * 10);
//			conditions.setColumnVector(1-1, svar2);
			//.set(/*conditions(:, 1) =*/ var2 * 10.0); //! MPa -> bar
			for (i = 0; i < nCalc; i++) vSalinity.setEntry(i, minCalc + delta * (i/* - 1.0*/));
//			conditions.setEntry(1-1, 3-1, /*);[1][3] =*/ minCalc);
//			conditions.setEntry(nCalc-1, 3-1, /*);[nCalc][3] =*/ maxCalc);
//			for (i = 2-1; i < nCalc - 1; i++) { //do i = 2, nCalc-1
//				conditions.setEntry(i, 3-1, /*);[i][3] =*/ minCalc + delta * (i/* - 1.0*/));
//			}//end do//! i
		}//end if //! typeCalc
		condition.setColumnVector(0, vPressure);
		condition.setColumnVector(1, vTemperature);
		condition.setColumnVector(2, vSalinity);

		System.out.println((condition.toString()));

		//nH2Oapp = 1000.0 / mH2O/*18.0153d0*/;

		/* START CALCULATION LOOP*/
		for (i = 0; i < nCalc; i++) { //do i = 1, nCalc
			solutionProperty conditions = new solutionProperty(
			/*pkb =*/ vPressure/*conditions*/.getEntry(i/*, 1-1*/)/*[i][1]*/ / 1000.0,
			/*tK =*/ vTemperature/*conditions*/.getEntry(i/*, 2-1*/)/*[i][2]*/ + 273.15,
			/*molalNaCl =*/ vSalinity/*conditions*/.getEntry(i/*, 3-1*/)/*[i][3]*/);

			System.out.println("Conditions Round" + i + " pkb " + conditions.pkb + " tK " + conditions.tK + " NaClapp " + conditions.molalNaCl);
//			get_random_seed();         //! new seed for MC (function included below)
			
			if ((conditions.isValid()) ||/*.OR.*/ (conditions.molalNaCl > 25.0)) {

//!//! Eval WH2O-CO2
				WH2OCO2(conditions /*pkb, tK*/, wWCO2, aH2O, aCO2);
//				System.out.println("WH2OCO2 is:     13.572535088673398        1.0000000000000000       0.89558076691256017      ");
//				System.out.println(returnWH2OCO2.toDebug());
				SIGMA(conditions /*new double[]{pkb, tK}*/, varAlpha, varW);
//				System.out.println("SIGMA is:      8.2567163577351721E-006   7.1977716089577370E-006 ");
//				System.out.println(returnSigma.toDebug());
//!//! Eval WH2O-NaCl
				WH2ONACL(conditions /*pkb, tK*/, wWNaClpm, wWNaClo, wNaCloNaClpm, aNaClpm, aNaCl0, varwWNaClpm, varwWNaClo, varaNaClpm, varaNaCl0);
//				System.out.println("WH2ONACL is:   -9.8533015872077971       -3.6730869468974059        0.0000000000000000       0.87000002836309576       0.99999570890471678       0.11742366148414315        1.1151441481968496        2.9525668734210888E-004   2.3318698437243223E-004 ");
//				System.out.println(returnWH2ONACL.toDebug());
//!//! Eval WCO2-NaCl
				WCO2NACL(conditions /*pkb, tK*/, wCO2NaCl0, sigwCO2NaCl0);
//				System.out.println("WCO2NACL is:        25.166138530221716       0.25793864920292037      ");
				returnWCO2NACL.wCO2NaClpm = returnWCO2NACL.wCO2NaCl0;
//				System.out.println(returnWCO2NACL.toDebug());

				/* randos is an array of 7 random numbers:
				 for W(H2O-CO2), for W(H2O-NaCl+-), for W(H2O-NaCl0),  for W(CO2-NaCl), with wNaCloNaClpm = 0
				 for alpha(CO2), for alpha(NaCl+-), for alpha(NaCl0)
				 First calculation is 'best' value*/
				Arrays.fill(randos, 0.5);
				//randos = new double[]{0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5};

				/*Runnable no input no ouput, just run */
				Runnable updateWsAlpha = () -> {
					/*draw Ws and alphas from their one-sigma uncertainty*/
					wWCO2MC.set(wWCO2 + (1.0 - 2.0 * randos[0]) * StDevW);
					wWNaCloMC.set(wWNaClo + (1.0 - 2.0 * randos[1]) * StDevwWNaClo);
					wWNaClpmMC.set(wWNaClpm + (1.0 - 2.0 * randos[2]) * StDevwWNaClpm);
					wCO2NaCl0MC.set(wCO2NaCl0 + (1.0 - 2.0 * randos[3]) * StDevwCO2NaCl0);
					aCO2MC.set(aCO2 + (1.0 - 2.0 * randos[4]) * StDevAlpha);
					aNaCl0MC.set(aNaCl0 + (1.0 - 2.0 * randos[5]) * StDevaNaCl0);
					aNaClpmMC.set(aNaClpm + (1.0 - 2.0 * randos[6]) * StDevaNaClpm);

					alphas = new double[]{1.0, aCO2MC.get(), aNaCl0MC.get(), aNaClpmMC.get()};

					VecMath.setcolumn(Wlist, 0, /*Wlist(:,1) =*/ new double[]{NaN, wWCO2MC.get(), wWNaCloMC.get(), wWNaClpmMC.get()});
					VecMath.setcolumn(Wlist, 1, /*Wlist(:,2) =*/ new double[]{wWCO2MC.get(), NaN, wCO2NaCl0MC.get(), wCO2NaCl0MC.get()});
					VecMath.setcolumn(Wlist, 2, /*Wlist(:,3) =*/ new double[]{wWNaCloMC.get(), wCO2NaCl0MC.get(), NaN, wNaCloNaClpm/*0.0*/});
					VecMath.setcolumn(Wlist, 3, /*Wlist(:,4) =*/ new double[]{wWNaClpmMC.get(), wCO2NaCl0MC.get(), wNaCloNaClpm/*0.0*/, NaN});

					/* for next iteration*/
					//randos = RANDOM_NUMBER(/*randos*/);
						/* randos is an array of 7 random numbers:
						 for W(H2O-CO2), for W(H2O-NaCl+-), for W(H2O-NaCl0),  for W(CO2-NaCl), with wNaCloNaClpm = 0
						 for alpha(CO2), for alpha(NaCl+-), for alpha(NaCl0)
						 this is at the end so that first calculation is not affected ( => gives "best" value)*/
					Arrays.parallelSetAll(randos, (r)-> get_random_seed().nextDouble());
//					System.out.println("randos is:   0.47619299537462145       0.37099532787991407       0.49380364883928840        4.5601366600334758E-002  0.56301741366335156       0.55730516130607943       0.60735933353067983      ");
//OK					System.out.println("this "+ Arrays.toString(randos));
				};

				for (mc.set(0); mc.get() < fast; mc.getAndIncrement()) {
					 {
						updateWsAlpha.run();

						// uncertainty is estimated from MonteCarlo only
					Future<double[]> yy = executor.submit(() -> PhaseDiagH2OCO2(conditions /*pkb, tK*/, wWCO2MC.get(), 1.0, aCO2MC.get()/*, y*/));
//Executor Ni plus rapide ni plus lent on garde

					//y = PhaseDiagH2OCO2(pkb, tK, wWCO2MC, 1.0, aCO2MC/*, y*/);
//					System.out.println("PhaseDiagH2OCO2 is:       13.572535088673398       0.89558076691256017        5.9672719645658132E-003  0.99678904240924293      ");
//					System.out.println("PhaseDiag 297 "+wWCO2MC+" "+ aCO2MC +" "+VecMath.toString(y));

							sCO2inW = yy.get()[0];
							sWinCO2 = yy.get()[1];

					compoCO2[mc.get()] = sWinCO2;

					if (false/*y[1 - 1] != /= yy[1 - 1]*/) { //THEN TODO VERIFY//! Fell in one-phase domain during error propagation

						nCO2[mc.get()] = NaN;
						xCO2[mc.get()] = NaN;
						compoCO2[mc.get()] = NaN;

					} else { //ELSE

						if (conditions.molalNaCl != /*/=*/ 0.0) {

							thermo.calcAllThermo(conditions);

							optionsSolubCO2.setSubVector(0, /*[1]:7) =*/ MatrixUtils.createRealVector(new double[]{conditions.pkb, conditions.tK, conditions.molalNaCl, nCO2[mc.get()],
							                                                                                       cacheThermo.DH /*DHNaClpm*/, cacheThermo.Gibbs /*mu0NaClpm*/, cacheThermo.Mu /*mu0NaCl0*/}));

							optionsSolubCO2.setEntry(27, /*[28] =*/ sCO2inW);
							optionsSolubCO2.setEntry(28, /*[29] =*/ sWinCO2);
//OK							System.out.println("(SolubCO2Brine optionsSolubCO2:    1.0132500000000001E-002   298.14999999999998        1.0000000000000000        1.9427451699511160E-318 -0.10669606044223041       -205.30571258344008       -417.15252696176361        1.0000000000000000       0.89558076691256017       0.99999570890471678       0.87000002836309576                            NaN   13.572535088673398       -3.6730869468974059       -9.8533015872077971        13.572535088673398                            NaN   25.166138530221716        25.166138530221716       -3.6730869468974059        25.166138530221716                            NaN   0.0000000000000000       -9.8533015872077971        25.166138530221716        0.0000000000000000                            NaN   5.9672719645658132E-003  0.99678904240924293      ");
//OK							System.out.println("SolubMC 337 "+VecMath.toString(optionsSolubCO2.toArray()));
////TODO Inversion des initialGuess
							/*1.0d-2//!initialGuess*/
							param[0]/*[1]*/ = sCO2inW * ((/*18.0153d0*/1000.0 / mH2O) / (1.0 - sCO2inW));
							param[1]/*[1]*/ = param[0] * 1.1 /*0.9*/; //!1.1d-2
//							//		/*(sCO2inW * ((1000.0 / mH2O/*18.0153d0*/) / (1.0 - sCO2inW)))*/

//							y[1 - 1] = SolubCO2Brine.apply(/*MatrixUtils.createRealVector*/new double[]{param[0], 0}, optionsSolubCO2.toArray());
//							y[2 - 1] = SolubCO2Brine.apply(/*MatrixUtils.createRealVector*/new double[]{param[1], 0}, optionsSolubCO2.toArray());

//							System.out.println("Solub2Brine fortran init is:    17.285305526025152        15.871609812810256      ");
//							System.out.println("SolubMC 322  initialGuessMin " + VecMath.toString(y));

							ftol = 1.0e-8;
							//niter = ITMAX;

							SolubCO2BrineMinimize xt = new SolubCO2BrineMinimize(param, optionsSolubCO2.toArray());
							molefracCO2 = xt.invoke()[0];

							if (false) {
								//molefracCO2 = param[1-1]/*[1]*/;
								//molefracCO2 = x.get()/*minimizer.x[0]*/[0];

								System.out.println("Recalculate mole per kg of water from mole fractions");
//!  nH2Oo = 1000.0/18.0153
								fracNaClapp = conditions.molalNaCl / (conditions.molalNaCl + nH2Oapp); //! in end-member...
								brineEM = MatrixUtils.createRealVector(new double[]{1.0 - fracNaClapp, fracNaClapp, 0.0}); //! H2O, NaCl, CO2
								gasEM = MatrixUtils.createRealVector(new double[]{1.0 - sWinCO2, 0.0, sWinCO2});
								compo = gasEM.mapMultiply(molefracCO2).add(brineEM.mapMultiply(1 - molefracCO2)); //molefracCO2 * gasEM + (1.0 - molefracCO2) * brineEM;
								xH2Oo = compo.getEntry(0)/*[1]*/;
								xNaClapp = compo.getEntry(1)/*[2]*/;
								xCO2[mc.get()] = compo.getEntry(2)/*[3]*/;
								coef = conditions.molalNaCl / xNaClapp; //!nH2Oo / xH2Oo
								conditions.molalNaCl = xNaClapp * coef; //! per kg of water
								nCO2[mc.get()] = xCO2[mc.get()] * coef; //! per kg of water
								System.out.println("oldCO2 " + nCO2[mc.get()] + " " + xCO2[mc.get()]);
							} else {
								nCO2[mc.get()] = molefracCO2; //x.get()/*minimizer.x[0]*/[0];
								//nCO2[mc] = param[1-1]/*[1]*/;
								xCO2[mc.get()] = molefracCO2/*nCO2[mc.get()]*/ / (molefracCO2/*nCO2[mc.get()]*/ + nH2Oapp + conditions.molalNaCl);
								//xCO2[mc] = param[1-1]/*[1]*/ / (param[1-1]/*[1]*/ + nH2Oapp + molalNaCl);
							} //end if//!0
//							System.out.println("newCO2 "+nCO2[mc]+" "+xCO2[mc]);
//							(nCO2 xCO2 init is:   0.26241231463338555        4.6223125354680821E-003 )

//! calc XNap for proportions...
							//On se fiche de ne copier qu'une partie d'optionsSolubCO2
							VecMath.copyvec(optionsSPECNaCl, optionsSolubCO2.toArray());
							//optionsSolubCO2.getSubVector(1 - 1, 27 - 1)/*1:27)*/;
//							optionsSPECNaCl[3-1] = molalNaCl;
							optionsSPECNaCl[3] = nCO2[mc.get()]; //! per kg of water at saturation

//							param[1 - 1]/*[1]*/ = 0.985; //!initialGuess
//							param[2 - 1]/*[1]*/ = param[1 - 1]/*[1]*/ * 1.001;

//							y[1 - 1] = SPECNaCl.apply(/*MatrixUtils.createRealVector*/(new double[]{param[0], 0}/*[1]/*,:)*/), optionsSPECNaCl.toArray());
//							y[2 - 1] = SPECNaCl.apply(/*MatrixUtils.createRealVector*/(new double[]{param[1], 0}/*[2]/*,:)*/), optionsSPECNaCl.toArray());
//							initialDelta = y[1] - y[0]; //PP
//							System.out.println("SolubMC 376 SPECNaCl initialGuessMin " + y[0] + " initiaGuessMax " + y[1]);

							ftol = 1.0e-6;
							//niter = ITMAX;

							SPECNaClMinimize xy = new SPECNaClMinimize(param, optionsSPECNaCl);
							XNapmc[mc.get()] = xy.invoke()[0];

//							XNapmc[mc.get()] = xx.get()[0]; //(xx/*minimizer.x[0]*/[1]/*[1]*/ + xx/*minimizer.x[0]*/[0]/*[1]*/) / 2.0;
							//XNapmc[mc] = (param[2-1]/*[1]*/ + param[1-1]/*[1]*/) / 2.0; //TODO verify index
							nNaClpm = 2.0 * XNapmc[mc.get()] * conditions.molalNaCl;
							nNaCl0 = (1.0 - XNapmc[mc.get()]) * conditions.molalNaCl;

//							System.out.println("XNapmc nNaClpm nNaCl0 is:   0.99712296185731997        1.9942459237146399        2.8770381426800329E-003 ");
//							System.out.println("this "+ XNapmc[mc] + " nNaClpm "+ molalNaCl);

						} else { //ELSE
							nNaCl0 = 0.0;
							nNaClpm = 0.0;
							XNapmc[mc.get()] = 0.0;
							xH2Oo = 1.0 - sCO2inW;
							xCO2[mc.get()] = sCO2inW;
							coef = nH2Oapp / xH2Oo;
							nCO2[mc.get()] = xCO2[mc.get()] * coef; //! per kg of water
						}//END IF//!molalNaCl == 0

						if (typeOutput == 3) {
							/* Calc Activity in water-rich phase*/
							ntotal = nH2Oapp + nCO2[mc.get()] + nNaCl0 + nNaClpm;
							proportions = VecMath.vxs(new double[]{nH2Oapp, nCO2[mc.get()], nNaCl0, nNaClpm}, 1/ntotal); // / ntotal;
							actH2O[mc.get()] = Math.exp(RTlnGamASF(proportions, Wlist, alphas, 0) / (Rm * conditions.tK)) * proportions[0];
							actNaCl0[mc.get()] = Math.exp(RTlnGamASF(proportions, Wlist, alphas, 2) / (Rm * conditions.tK)) * proportions[2];
							actNaClpm[mc.get()] = Math.exp(RTlnGamASF(proportions, Wlist, alphas, 3) / (Rm * conditions.tK)) * proportions[3];
							actCO2[mc.get()] = Math.exp(RTlnGamASF(proportions, Wlist, alphas, 1) / (Rm * conditions.tK)) * proportions[1];

						}//end if//! (typeOutput == 3 )

					}//END IF//! One phase domain during error prop
				}
					}//END DO//! mc = 1, 100

//! estimate uncertainty from spread of MonteCarlo
// PP just stdev OK pour simplification
//				nCO2mean = VecMath.sum(nCO2/*(:)*/)/100.0; //! should fall very close to nCO2(1) - NB: median is spot on, but too computer-intensive for MC simulation
//				uncertM = new StandardDeviation(false).evaluate(nCO2, nCO2mean); //*/Math.sqrt(sum(Math.pow((nCO2/*(:)*/-nCO2mean), 2.0))/100.0);
				uncertM = FastMath.sqrt(StatUtils.variance(nCO2)); //
//				xCO2mean = VecMath.sum(xCO2/*(:)*/)/100.0; //! should fall very close to xCO2(1)
//				uncerty[1 - 1] = new StandardDeviation(false).evaluate(xCO2, xCO2mean); //*/Math.sqrt(sum(Math.pow((xCO2/*(:)*/-xCO2mean), 2.0))/100.0);
				uncerty[0] = FastMath.sqrt(StatUtils.variance(xCO2)); //
//				compoCO2mean = VecMath.sum(compoCO2/*(:)*/)/100.0; //! should fall very close to compoCO2(1)
//				uncerty[2 - 1] = new StandardDeviation(false).evaluate(compoCO2, compoCO2mean); //*/Math.sqrt(sum(Math.pow((compoCO2/*(:)*/-compoCO2mean), 2.0))/100.0);
				uncerty[1] = FastMath.sqrt(StatUtils.variance(compoCO2)); //

				// mole fractions - note that xCO2 uncertainty is not strictly valid yet for brine
				solub.setRow(i/*,:)*/, new double[]{nCO2[0], uncertM, xCO2[0], uncerty[0], compoCO2[0], uncerty[1]});

				if (typeOutput == 3) { 

//! save in array
					thermoString.put(0, "a(H2O) 0.95609 0.8712");thermoVal[i][0] = actH2O[0];
					thermoString.put(1, "σ aH2O 0.00022 0.0015");thermoVal[i][1]= FastMath.sqrt(StatUtils.variance(actH2O)); //stdevactH2O;
					thermoString.put(2, "a(CO2) 1.0024 1.0247");thermoVal[i][2]= actCO2[0];
					thermoString.put(3, "σ aCO2 0.0011 0.0051");thermoVal[i][3]= FastMath.sqrt(StatUtils.variance(actCO2)); // stdevactCO2;
					thermoString.put(4, "a(NaCl0) 1.37E-5 0.000223");thermoVal[i][4]= actNaCl0[0];
					thermoString.put(5, "σ aNaCl0 4.1E-6 5.1E-5");thermoVal[i][5]= FastMath.sqrt(StatUtils.variance(actNaCl0)); // stdevactNaCl0;
					thermoString.put(6, "a(NaCl± 0.0011 0.00444");thermoVal[i][6]= actNaClpm[0];
					thermoString.put(7, "σ aNaCl± 0.00016 0.0005");thermoVal[i][7]= FastMath.sqrt(StatUtils.variance(actNaClpm)); // stdevactNaClpm;
					thermoString.put(8, "X(Na+) 0.99712 0.98742");thermoVal[i][8]= XNapmc[0];
					thermoString.put(9, "W(H2O-CO2) kJ/mol 13.5728 ...");thermoVal[i][9]= wWCO2;
					thermoString.put(10,"σ WCO2 0.0054 ...");thermoVal[i][10]= StDevW;
					thermoString.put(11, "W(CO2-NaCl) kJ/mol 25.2 ...");thermoVal[i][11]= wCO2NaClpm;
					thermoString.put(12, "σ wCO2NaCl 1 ...");thermoVal[i][12]= StDevwCO2NaCl0;
					thermoString.put(13, "W(H2O-NaClo) kJ/mol -3.67309 ...");thermoVal[i][13]= wWNaClo;
					thermoString.put(14, "σ wNaClo 2.11201 ...");thermoVal[i][14]= StDevwWNaClo;
					thermoString.put(15, "W(H2O-NaCl±) kJ/mol -9.85 ...");thermoVal[i][15]= wWNaClpm;
					thermoString.put(16, "σ wNaCl± 0.69 ...");thermoVal[i][16]= StDevwWNaClpm;
					thermoString.put(17, "α CO2 0.8966 ...");thermoVal[i][17]= aCO2;
					thermoString.put(18, "σ αCO2 0.0058 ...");thermoVal[i][18]= StDevAlpha;
					thermoString.put(19, "α NaCl0 1 ...");thermoVal[i][19]= aNaCl0;
					thermoString.put(20, "σ αNaCl0 0.031 ...");thermoVal[i][20]= StDevaNaCl0;
					thermoString.put(21, "α NaCl± 0.87 ...");thermoVal[i][21]= aNaClpm;
					thermoString.put(22, "σ αNaCl± 0.034 ...");thermoVal[i][22]= StDevaNaClpm;

				}//end if//! (typeOutput == 3 )

			} else { //ELSE //! one-phase domain (not valid)
				if (typeOutput == 3) VecMath.fillvec(thermoVal[i],NaN); //! originally, nan
				RealVector setnan = solub.getRowVector(i);
				setnan.set(NaN);
				solub.setRowVector(i, setnan);

			}//End If //! (valid) == one phase domain

			if (typeOutput !=/*/=*/1) { 

				densityCO2 = CO2EOS.CO2DENSITY(conditions/*new double[]{pkb, tK}*//*, densityCO2*/);
				//!       where conditions = (/pkb, tK/)
				//!       returns density in mol/cm³
				densities.setEntry(i, 0, densityCO2 * mCO2);
				viscosities.setEntry(i, 0, viscoCO2(conditions /*pkb, tK*/));

				densityWater = WATEREOS.waterdensity(conditions/*new double[]{pkb, tK}*//*, densityWater*/);
				//!   where conditions is (/ pkb, tK/) density is returned in kg/m³ so have to divide by 10³
				densities.setEntry(i, 1, densityWater / 1000.0);
				densities.setEntry(i, 2, diecIAPWS97(conditions/*new double[]{pkb, tK}*//*, densityWater*/));
				viscosities.setEntry(i, 1, 1.0e6 * viscoWater(conditions /*pkb, tK*/)); //! micro Pa.s

//TODO
//				fugwater(conditions);
//viscoH2ONaClCO2(pkb, tK, xCO2[0]);
//viscoSystem(tK, 2);

				if (conditions.molalNaCl !=/*/=*/ 0.0) {   //! viscosity ratio
					viscosities.setEntry(i, 2, ratioVisco(conditions /*tK*//*, molalNaCl*/));
					viscosities.setEntry(i, 3, Bando2004(conditions /*pkb, tK*/, nCO2[0]/*, molalNaCl*/));
				} else { //ELSE
					viscosities.setEntry(i, 2, 0.0);
					viscosities.setEntry(i, 3, 0.0);
				}//END IF
			}//end if//! typeOutput not 1

		}//end do//! i = 1, ncalc
//! END CALCULATION LOOP

		//conditions(:,1) = conditions(:,1) / 10.0; //!bar -> MPa
		condition.setColumn(0, condition.getColumnVector(0).mapDivideToSelf(10.0).toArray());

		System.out.println("Your input");
		csv_write("Pressure 1.01325 MPa", condition.getColumnVector(0)/*(:,1)*/, true);
		csv_write("Temperature 25°C", condition.getColumnVector(1)/*(:,2)*/, true);
		csv_write("Salinity 1 3 mol/kgH2O", condition.getColumnVector(2)/*(:,3)*/, true);
		System.out.println("Calculated solubility and associated uncertainty.");
		System.out.println("CO2 solubility in water-rich phase");
		csv_write("molality 0.2615 0.1764", solub.getColumnVector(0)/*(:,1)*/, true);
		csv_write("σ 0.0038 0.0048", solub.getColumnVector(1)/*(:,2)*/, true);
		csv_write("mol. frac.0.004607 0.003006", solub.getColumnVector(2)/*(:,3)*/, true);
		csv_write("σ 6.7E-5 8.2E-5", solub.getColumnVector(3)/*(:,4)*/, true);
		System.out.println("CO2 content of CO2-rich phase");
		csv_write("mol. frac. 0.996779", solub.getColumnVector(4)/*(:,5)*/, true);
		csv_write("σ 3.5E-5", solub.getColumnVector(5)/*(:,6)*/, true);

		if (typeOutput != /*/=*/ 1) { 
			System.out.println("Above values apply to the \"pure\" phases, not to their mixture (see [7, 8, 9])");
			csv_write("CO2 density 0.0189017 g.cm-3", densities.getColumnVector(0)/*(:,1)*/, true);
			csv_write("H2O density 0.997457 g.cm-3", densities.getColumnVector(1)/*(:,2)*/, true);
			csv_write("CO2 viscosity 15.027 μPa.s", viscosities.getColumnVector(0)/*(:,1)*/, true);
			csv_write("H2O viscosity 889.9 μPa.s", viscosities.getColumnVector(1)/*(:,2)*/, true);
			System.out.println("Calc. with ref. [10] and [11].");
			System.out.println("Viscosity ratio ηr");
			csv_write("CO2-free 1.09582 1.35786", viscosities.getColumnVector(2)/*(:,3)*/, true);
			csv_write("CO2 & NaCl 0.96723 0.99806", viscosities.getColumnVector(3)/*(:,4)*/, true);
			if (typeOutput == 3) { 
				System.out.println("Activity of components and associated uncertainties and association of NaCl. nd: dimensionless");
				for (i = 0; i < 23; i++) { //DO i = 1, 23
					csv_write(thermoString.get(i), VecMath.getcolumn(thermoVal, i)/*(:,i)*/, true);
				} //END DO
			}//END IF
		}//END IF
		Log.i("Cache CO2 Hit", "Density "+ CO2EOS.cacheCO2.ask);
		Log.i("Cache H2O Hit", "Density "+ WATEREOS.cacheH2O.ask);
		Log.i("Cache Thermo Hit", "All "+ cacheThermo.ask);
		cacheFunction.printCalled();
		cacheThermo.printCalled();

		System.out.println("Time " + (System.nanoTime()-start) );
		//23092018 Time 219274644013
		//Time          219654760321 sans minimize
		//Time          219867555321 sans preminimize
		//Time          108497616776 100 Loop au lieu de 150
		//06102018 Time  30800061079 Minimize au lieu d'amoebaCSharp L 362
		//09102018 Time  15029345462 class solutionProperty pour tout le monde
	}

static class cacheFunction {
	static int callSPECNaCl = 0;
	static int callSoluBCO2Brine = 0;

	static void printCalled() { Log.w("Cache","SPECNaCl "+callSPECNaCl + " SolubCO2Brine " + callSoluBCO2Brine); }
}


	private static ThreadLocalRandom get_random_seed() {
		return ThreadLocalRandom.current();
//		//use iso_fortran_env, only: int64
//		//implicit none
//		integer, allocatable::seed (:)
//		integer::i, n, un, istat, dt(8), pid integer(int64)::t
//
//				/*call*/ */random_seed(size = n) allocate(seed(n))
//		//! First try if the OS provides a random number generator
//		open(newunit = un, file = "/dev/urandom", access = "stream", & form = "unformatted", action = "read", status = "old", iostat = istat)
//		if (istat == 0) then read(un) seed close(un)
//   else{
//			//! Fallback to XOR:ing the current time and pid. The PID is
//			//! useful in case one launches multiple instances of the same
//			//! program in parallel.
//			/*call*/
//			system_clock(t) if (t == 0) {//then
//				/*call*/
//				date_and_time(values = dt) t = (dt(1) - 1970) * 365_
//				int64 * 24 * 60 * 60 * 1000 & +dt(2) * 31_
//				int64 * 24 * 60 * 60 * 1000 & +dt(3) * 24_
//				int64 * 60 * 60 * 1000 & +dt(5) * 60 * 60 * 1000 & +dt(6) * 60 * 1000 + dt(7) * 1000 & +dt(8)
//				end if pid = getpid() t = ieor(t, int(pid, kind(t)))
//				do i = 1, n seed(i) = lcg(t) end do end if
//				/*call*/
//				random_seed(put = seed);
//				//!        contains
//				//! This simple PRNG might not be good enough for real work, but is
//				//! sufficient for seeding a better PRNG.
//
	}//end subroutine init_random_seed
//
//			int /*function*/ lcg ( int s){
//				//use iso_fortran_env, only: int64
//				//integer :: lcg
//				//integer(int64) :: s
//				if (s == 0) {//then
//					s = 104729;
//				} else {
//					s = mod(s, 4294967296/*_int64*/);
//				}//end if
//				s = mod(s * 279470273/*_int64*/, 4294967291/*_int64*/);
//				lcg = int(mod(s, int(huge(0), int64)),kind(0));
//				return lcg;
//			}//end function lcg


//	public static AmoebaCSharp.BiVectorFunction/*BiFunction<double[], double[], Double>*/ testGauss = (simplex, options) -> 2.0 - gaussian(simplex[0]) - gaussian(simplex[1]);
//
//	private static double gaussian(double x) {
//		double d = x - 0.5;
//		return exp(-d * d);
//	}
//	public static double[] test() {
//		return minimize(testGauss, new double[]{0, 0}, 0.1, ITMAX, null, null);
//	}

}//END PROGRAM FLUIDmin