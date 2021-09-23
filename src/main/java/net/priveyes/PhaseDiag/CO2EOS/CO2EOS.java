package net.priveyes.PhaseDiag.CO2EOS;

import static net.priveyes.PhaseDiag.FluidMin.Constants.Rm;

import net.priveyes.PhaseDiag.num_rep.VecMath;
import net.priveyes.PhaseDiag.FluidMin.solutionProperty;
import net.priveyes.PhaseDiag.FluidMin.Minimizer;
import net.priveyes.PhaseDiag.FluidMin.Minimizer.BiVectorFunction;

import java.util.concurrent.ExecutionException;

/**the CO2 equation of state of Sterner and Pitzer (1994)[7] has been chosen for the volume/density of CO2
  Benoit.Dubacq@upmc.fr  30/10/2013*/
public class  CO2EOS {
	static double[][] csCO2 = new double[10][6];
	static double[][] csAtT = new double[10][6];
	/**
	 * (ii) the SUBROUTINE volumeCO2(conditions,volumeJBar)
	 * where conditions = (/pkb, tK/)
	 * returns volume in J/bar
	 */
	public static double /*SUBROUTINE*/ volumeCO2(solutionProperty /*double[]*/ conditions/*, Double volumeJBar*/) throws ExecutionException, InterruptedException {
//USE nrtype

// Double /*REAL(DP), DIMENSION(2), INTENT(IN) ::*/ conditions;
		double /*REAL(DP), INTENT(OUT) ::*/ volumeJBar;
		double rho = CO2DENSITY(conditions/*, rho*/);

		volumeJBar = 1.0 / (10.0 * rho); /*rho is in mol/cm³*/
		return volumeJBar;
	}//END void /*SUBROUTINE*/ volumeCO2 {

	/**
	 * (iii) the SUBROUTINE fugCO2(conditions, fugacity)
	 * where conditions = (/pkb, tK/)
	 * returns log(fugacity)
	 */
	public static double /*SUBROUTINE*/ fugCO2(solutionProperty/*double[]*/ conditions/*, Double fugacity*/) throws ExecutionException, InterruptedException {

//	Double /*REAL(DP), DIMENSION(2), INTENT(IN) ::*/ conditions;
		double /*REAL(DP), INTENT(OUT) ::*/ fugacity;
		double  rho, AresdivRT;
		double  c1, c2, c3, c4, c5, c6, c7, c8, c9, c10;//!, pdivRT
		double  /*tK,*/ pbar;
		//RealMatrix /*REAL(DP), DIMENSION(10,6) ::*/ csCO2 = MatrixUtils.createRealMatrix(10 /*+ 1*/, 6 /*+ 1*/),
		//		csAtT = MatrixUtils.createRealMatrix(10 /*+ 1*/, 6 /*+ 1*/);
		//double[][] csCO2 = new double[10][6], csAtT = new double[10][6];

		rho = CO2DENSITY(conditions/*, rho*/); //! calc density

		pbar = conditions.pkb /*[1 - 1]*/ * 1000;
		//tK = conditions[2 - 1];

//		csCO2[1-1] = new double[]{0.0, 0.0, 0.18261340e7, 0.79224365e2, 0.0, 0.0};
//		csCO2[2-1] = new double[]{ 0.0, 0.0, 0.0, 0.66560660e-4, 0.57152798e-5, 0.3022236e-9};
//		csCO2[3-1] = new double[]{ 0.0, 0.0, 0.0, 0.59957845e-2, 0.71669631e-4, 0.62416103e-8};
//		csCO2[4-1] = new double[]{ 0.0, 0.0, -0.13270279e1, -0.15210731, 0.53654244e-3, -0.71115142e-7};
//		csCO2[5-1] = new double[]{ 0.0, 0.0, 0.12456776, 0.49045367e1, 0.98220560e-2, 0.55962121e-5};
//		csCO2[6-1] = new double[]{ 0.0, 0.0, 0.0, 0.75522299, 0.0, 0.0};
//		csCO2[7-1] = new double[]{ -0.39344644e12, 0.90918237e8, 0.42776716e6, -0.22347856e2, 0.0, 0.0};
//		csCO2[8-1] = new double[]{ 0.0, 0.0, 0.40282608e3, 0.11971627e3, 0.0, 0.0};
//		csCO2[9-1] = new double[]{ 0.0, 0.22995650e8, -0.78971817e5, -0.63376456e2, 0.0, 0.0};
//		csCO2[10-1] = new double[]{ 0.0, 0.0, 0.95029765e5, 0.18038071e2, 0.0, 0.0};
//
		VecMath.setcolumn(csAtT/*[][1]*/, 1 - 1, VecMath.vxs(VecMath.getcolumn(csCO2, /*[][1]*/1 - 1), (Math.pow(conditions.tK, -4))));
		VecMath.setcolumn(csAtT/*[][2]*/, 2 - 1, VecMath.vxs(VecMath.getcolumn(csCO2, /*[][2]*/2 - 1), (Math.pow(conditions.tK, -2))));
		VecMath.setcolumn(csAtT/*[][3]*/, 3 - 1, VecMath.vxs(VecMath.getcolumn(csCO2, /*[][3]*/3 - 1), (Math.pow(conditions.tK, -1))));
		VecMath.setcolumn(csAtT/*[][4]*/, 4 - 1, VecMath.getcolumn(csCO2, /*[][4]*/4 - 1));
		VecMath.setcolumn(csAtT/*[][5]*/, 5 - 1, VecMath.vxs(VecMath.getcolumn(csCO2, /*[][5]*/5 - 1), conditions.tK));
		VecMath.setcolumn(csAtT/*[][6]*/, 6 - 1, VecMath.vxs(VecMath.getcolumn(csCO2, /*[][6]*/6 - 1), (Math.pow(conditions.tK, 2))));

		c1 = VecMath.sum(csAtT[1 - 1]);
		c2 = VecMath.sum(csAtT[2 - 1]);
		c3 = VecMath.sum(csAtT[3 - 1]);
		c4 = VecMath.sum(csAtT[4 - 1]);
		c5 = VecMath.sum(csAtT[5 - 1]);
		c6 = VecMath.sum(csAtT[6 - 1]);
		c7 = VecMath.sum(csAtT[7 - 1]);
		c8 = VecMath.sum(csAtT[8 - 1]);
		c9 = VecMath.sum(csAtT[9 - 1]);
		c10 = VecMath.sum(csAtT[10 - 1]);

		AresdivRT = c1 * rho + 1.0 / (c2 + c3 * rho + c4 * Math.pow(rho, 2) + c5 * Math.pow(rho, 3) + c6 * Math.pow(rho, 4)) - (1 / c2) - (c7 / c8) * (Math.exp(-c8 * rho) - 1) - (c9 / c10) * (Math.exp(-c10 * rho) - 1);

		fugacity = Math.log(rho) + AresdivRT + pbar / (10 * 1000 * Rm /*8.314472*/ * conditions.tK * rho)
				+ Math.log(10 * 1000* Rm /*8.314472*/ * conditions.tK) - 1;
		return fugacity;
	}//END void /*SUBROUTINE*/ fugCO2 {


	public static class cacheCO2 {
		static solutionProperty conditions = null;
		static Double density = null;
		public static int ask = 0;
		static boolean inCacheDensity(solutionProperty askConditions) {ask +=1;return askConditions == conditions;}
	}
	/**
	 * (i) the SUBROUTINE CO2DENSITY(conditions, density)
	 * where conditions = (/pkb, tK/)
	 * returns density in mol/cm³
	 */
	public static double /*SUBROUTINE*/ CO2DENSITY(solutionProperty /*double[]*/ conditions/*, Double density*/) throws ExecutionException, InterruptedException {
		if (cacheCO2.inCacheDensity(conditions)) return cacheCO2.density;
		cacheCO2.conditions = conditions;
		//USE nrtype
//USE MINIMIZE //! Anonymized Nelder-Mead

//Double /*REAL(DP), DIMENSION(2), INTENT(IN) ::*/ conditions;
		double /*REAL(DP), INTENT(OUT) ::*/ density;
//INTERFACE
		double[] /*REAL(DP), DIMENSION(2) ::*/ y = new double[2 /*+ 1*/];
		double[] /*RealVector*/ conditionsDim;
		double[]/*[]*/ /*REAL(DP), DIMENSION(2,1) ::*/ param = new double[2 /*+ 1*/]/*[1 + 1]*/;
		double  ftol, /*pkb, tK,*/ pvap, pbar, initialGuess;
		int /*INTEGER(I4B) ::*/ niter = Minimizer.ITMAX;

		//pkb = conditions[1 - 1];
		//tK = conditions[2 - 1];
		pbar = 1000.0 * conditions.pkb;

		conditionsDim = (new double[]{pbar, conditions.tK}); //(/pbar, tK/);
		ftol = 1.0e-7;

		if (pbar > 73.8 ||/*.OR.*/ conditions.tK > 304.2) {//THEN
			initialGuess = 0.02;
		} else {//ELSE
			pvap = 73.825 * Math.exp(11.377371 * Math.pow((1.0 - conditions.tK / 304.21), 1.935) - 6.8849249 * (-1.0 + 304.21 / conditions.tK) - 9.5924263 * Math.pow((-1.0 + 304.21 / conditions.tK), 2) + 13.679755 * Math.pow((-1.0 + 304.21 / conditions.tK), 3.0) - 8.6056439 * Math.pow((-1.0 + 304.21 / conditions.tK), 4.0));
			if (pbar < pvap) {//THEN
				initialGuess = 0.0001; //! vapour phase
				//OK System.out.println("CO2EOS Vapour Pressure pbar " + pbar + " pvap " + pvap + " tK " + tK + " initialGuess " + initialGuess);
			} else {//ELSE
				initialGuess = 2.0e-2; //!liquid
				//OK System.out.println("CO2EOS Liquid Pressure pbar " + pbar + " pvap " + pvap + " tK " + tK + " initialGuess " + initialGuess);
			}//END IF//! if pbar<pvap
		}//END IF//! if

		param[1 - 1]/*[1]*/ = initialGuess;
		param[2 - 1]/*[1]*/ = param[1 - 1]/*[1]*/ * 1.01;
//		long start = System.nanoTime();
//		y[1 - 1] = DENSCO2OBJ.apply((new double[]{param[0], 0}/*[1]/*[]*/), conditionsDim);
//		y[2 - 1] = DENSCO2OBJ.apply((new double[]{param[1], 0}/*[2]/*[]*/), conditionsDim);
		//System.out.println("CO2EOS 148 DENSCO2OBJ initialGuessMin  " + y[0] + " initiaGuessMax " + y[1]);
//		double initialDelta = y[1] - y[0]; //PP

		double[] x = Minimizer.minimize/*NELDERMEAD*/(DENSCO2OBJ, param, 0.01/*initialDelta*/ /*,ftol*/, /*DENSCO2OBJ,*/ niter, null, conditionsDim/*.toArray()*/);
//		Log.w("Minimize", "Old " + (System.nanoTime() - start));
//		System.out.println("CO2EOS 161 DENSCO2OBJ Solution: x = " + VecMath.toString(x));

		density = x/*minimizer.x[0]*/[0];
//		density = param[1-1]/*[1]*/;
//		start = System.nanoTime();
//		Amoeba amo = new Amoeba(2, 1, param[0], param[1], 100, DENSCO2OBJ, conditionsDim /*new double[]{1013.25, 298.15}*/);
//		Solution sol = amo.Solve();
//		Log.w("Minimize", "New " + (System.nanoTime() - start));
//		Log.e("Amoeaba CO2EOS", VecMath.toString(sol.vector));

//		System.out.println(bestSolution.toString()); // print the best solution
//		density = sol.vector[0];

		cacheCO2.density = density;
		return density;
	}//END void /*SUBROUTINE*/ CO2DENSITY {

//END MODULE CO2EOS
//INTERFACE
//interface double  DENSCO2OBJ(double[] rho, double[] options) {
	//! pressure in bar, temperature in K, density in mol/cm³
	//USE nrtype

//	double /*REAL(DP), DIMENSION(:) ::*/ rho;
//	double  DENSCO2OBJ;
//	double /*REAL(DP), DIMENSION(:) ::*/ options;
//}//END FUNCTION DENSCO2OBJ
//END INTERFACE

	/**
	 * pressure in bar, temperature in K, density in mol/cm³
	 */
	private static BiVectorFunction DENSCO2OBJ
			/*(double[] rho, double[] options)*/ = (rho, options) -> {

//Double /*REAL(DP), DIMENSION(:) ::*/ rho;
		double  DENSCO2OBJ;
//Double /*REAL(DP), DIMENSION(:) ::*/ options;
		double  c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, pdivRT;
		double[] /*REAL(DP), DIMENSION(1) ::*/ tmp = new double[1 /*+ 1*/];
		double tK, pbar;
//		RealMatrix /*REAL(DP), DIMENSION(10,6) ::*/ csCO2 = MatrixUtils.createRealMatrix(10 /*+ 1*/, 6 ),
//				csAtT = MatrixUtils.createRealMatrix(10 /*+ 1*/, 6 /*+ 1*/);

		pbar = options[1 - 1];
		tK = options[2 - 1];
//
//		csCO2.setRow(1-1, /*[1] =*/ new double[]{ 0.0, 0.0, 0.18261340e7, 0.79224365e2, 0.0, 0.0});
//		csCO2.setRow(2-1, /*[2] =*/ new double[]{ 0.0, 0.0, 0.0, 0.66560660e-4, 0.57152798e-5, 0.3022236e-9});
//		csCO2.setRow(3-1, /*[3] =*/ new double[]{ 0.0, 0.0, 0.0, 0.59957845e-2, 0.71669631e-4, 0.62416103e-8});
//		csCO2.setRow(4-1, /*[4] =*/ new double[]{ 0.0, 0.0, -0.13270279e1, -0.15210731, 0.53654244e-3, -0.71115142e-7});
//		csCO2.setRow(5-1, /*[5] =*/ new double[]{ 0.0, 0.0, 0.12456776, 0.49045367e1, 0.98220560e-2, 0.55962121e-5});
//		csCO2.setRow(6-1, /*[6] =*/ new double[]{ 0.0, 0.0, 0.0, 0.75522299, 0.0, 0.0});
//		csCO2.setRow(7-1, /*[7] =*/ new double[]{ -0.39344644e12, 0.90918237e8, 0.42776716e6, -0.22347856e2, 0.0, 0.0});
//		csCO2.setRow(8-1, /*[8] =*/ new double[]{ 0.0, 0.0, 0.40282608e3, 0.11971627e3, 0.0, 0.0});
//		csCO2.setRow(9-1, /*[9] =*/ new double[]{ 0.0, 0.22995650e8, -0.78971817e5, -0.63376456e2, 0.0, 0.0});
//		csCO2.setRow(10-1, /*[10] =*/ new double[]{ 0.0, 0.0, 0.95029765e5, 0.18038071e2, 0.0, 0.0});

		VecMath.setcolumn(csAtT/*[][1]*/, 1 - 1, VecMath.vxs(VecMath.getcolumn(csCO2, /*[][1]*/1 - 1), (Math.pow(tK, -4))));
		VecMath.setcolumn(csAtT/*[][2]*/, 2 - 1, VecMath.vxs(VecMath.getcolumn(csCO2, /*[][2]*/2 - 1), (Math.pow(tK, -2))));
		VecMath.setcolumn(csAtT/*[][3]*/, 3 - 1, VecMath.vxs(VecMath.getcolumn(csCO2, /*[][3]*/3 - 1), (Math.pow(tK, -1))));
		VecMath.setcolumn(csAtT/*[][4]*/, 4 - 1, VecMath.getcolumn(csCO2, /*[][4]*/4 - 1));
		VecMath.setcolumn(csAtT/*[][5]*/, 5 - 1, VecMath.vxs(VecMath.getcolumn(csCO2, /*[][5]*/5 - 1), tK));
		VecMath.setcolumn(csAtT/*[][6]*/, 6 - 1, VecMath.vxs(VecMath.getcolumn(csCO2, /*[][6]*/6 - 1), (Math.pow(tK, 2))));

		c1 = VecMath.sum(csAtT[1 - 1]);
		c2 = VecMath.sum(csAtT[2 - 1]);
		c3 = VecMath.sum(csAtT[3 - 1]);
		c4 = VecMath.sum(csAtT[4 - 1]);
		c5 = VecMath.sum(csAtT[5 - 1]);
		c6 = VecMath.sum(csAtT[6 - 1]);
		c7 = VecMath.sum(csAtT[7 - 1]);
		c8 = VecMath.sum(csAtT[8 - 1]);
		c9 = VecMath.sum(csAtT[9 - 1]);
		c10 = VecMath.sum(csAtT[10 - 1]);

		pdivRT = pbar / (10.0 * 8.314472 * tK);

		tmp[1 - 1] = Math.pow((-pdivRT + rho[1 - 1] + c1 * (Math.pow(rho[1 - 1], 2.0)) - (Math.pow(rho[1 - 1], 2.0)) * ((c3 + 2.0 * rho[1 - 1] * c4 + 3.0 * c5 * (Math.pow(rho[1 - 1], 2.0)) + 4.0 * c6 * (Math.pow(rho[1 - 1], 3.0))) / (Math.pow((c2 + c3 * rho[1 - 1] + c4 * (Math.pow(rho[1 - 1], 2.0)) + c5 * (Math.pow(rho[1 - 1], 3.0)) + c6 * (Math.pow(rho[1 - 1], 4.0))), 2.0))) + c7 * (Math.pow(rho[1 - 1], 2.0)) * Math.exp(-c8 * rho[1 - 1]) + c9 * (Math.pow(rho[1 - 1], 2.0)) * Math.exp(-c10 * rho[1 - 1])), 2.0);

		DENSCO2OBJ = tmp[1 - 1];
		return DENSCO2OBJ;
	};//END FUNCTION DENSCO2OBJ

	static {
		csCO2[1 - 1] = new double[]{0.0, 0.0, 0.18261340e7, 0.79224365e2, 0.0, 0.0};
		csCO2[2 - 1] = new double[]{0.0, 0.0, 0.0, 0.66560660e-4, 0.57152798e-5, 0.3022236e-9};
		csCO2[3 - 1] = new double[]{0.0, 0.0, 0.0, 0.59957845e-2, 0.71669631e-4, 0.62416103e-8};
		csCO2[4 - 1] = new double[]{0.0, 0.0, -0.13270279e1, -0.15210731, 0.53654244e-3, -0.71115142e-7};
		csCO2[5 - 1] = new double[]{0.0, 0.0, 0.12456776, 0.49045367e1, 0.98220560e-2, 0.55962121e-5};
		csCO2[6 - 1] = new double[]{0.0, 0.0, 0.0, 0.75522299, 0.0, 0.0};
		csCO2[7 - 1] = new double[]{-0.39344644e12, 0.90918237e8, 0.42776716e6, -0.22347856e2, 0.0, 0.0};
		csCO2[8 - 1] = new double[]{0.0, 0.0, 0.40282608e3, 0.11971627e3, 0.0, 0.0};
		csCO2[9 - 1] = new double[]{0.0, 0.22995650e8, -0.78971817e5, -0.63376456e2, 0.0, 0.0};
		csCO2[10 - 1] = new double[]{0.0, 0.0, 0.95029765e5, 0.18038071e2, 0.0, 0.0};
	}
}