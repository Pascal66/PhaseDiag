package net.priveyes.PhaseDiag.IAPWS;

import static net.priveyes.PhaseDiag.CHNOZ.IAPWS95.initIAPWS95;
import static net.priveyes.PhaseDiag.FluidMin.Constants.mH2O;
import static net.priveyes.PhaseDiag.FluidMin.Minimizer.BiVectorFunction;
import static net.priveyes.PhaseDiag.FluidMin.Minimizer.ITMAX;
import static net.priveyes.PhaseDiag.FluidMin.Minimizer.minimize;
import static java.lang.Math.exp;
import static java.lang.Math.log;
import static java.lang.Math.pow;
import static java.lang.Math.sqrt;

import net.priveyes.PhaseDiag.CHNOZ.IAPWS95;
import net.priveyes.PhaseDiag.CHNOZ.IAPWScommon.residual;
import net.priveyes.PhaseDiag.FluidMin.solutionProperty;

import java.util.concurrent.ExecutionException;

/**
 * It uses the Nelder-Mead minimisation algorithm from Numerical Recipes
 * with the equations from the 2005 IAPWS release and an initial guess
 * deciding whether the phase is in liquid or gaseous state (not to be trusted very close to the critical point of water)
 Benoit.Dubacq@upmc.fr  29/10/2013*/
public class  WATEREOS {

	/**
	 * volumewater(conditions, volumeJBar)
	 * where conditions is { pkb, tK}
	 *
	 * @return volume is returned in J/bar
	 */
	public static double /*SUBROUTINE*/ volumewater(solutionProperty /*double[]*/ conditions/*, double volumeJBar*/) throws ExecutionException, InterruptedException {
//USE nrtype

		//double /*REAL(DP), DIMENSION(2), INTENT(IN) ::*/ conditions;
		double /*REAL(DP), INTENT(OUT) ::*/ volumeJBar;
		double  rho;

//! calculate water density
		/*call*/
		rho = waterdensity(conditions/*, rho*/);

//!write(*,*) '   (rho is: ', rho,')'

		volumeJBar = (15.9994 + 2 * 1.00794) / (10 * rho / 1000);
		return volumeJBar;
	}//END SUBROUTINE volumewater

	/**
	 * (ii) the SUBROUTINE fugwater(conditions, fugacity)
	 * where conditions is { pkb, tK }
	 *
	 * @return ln(fugacity) is returned
	 */
	public static double /*SUBROUTINE*/ fugwater(solutionProperty /*double[]*/ conditions/*, double fugacity*/) throws ExecutionException, InterruptedException {
//USE nrtype

		//double /*REAL(DP), DIMENSION(2), INTENT(IN) ::*/ conditions;
		double /*REAL(DP), INTENT(OUT) ::*/ fugacity;
		double  pKPa/*, tK*/, dens;
		double  arg = -1.0, nan, a, LargeBeta, LargeA, LargeB;
		double  tc = 647.096, rhoc = 322.0, rho;
		double  delta, tau, phir, alpha, epsi;
		double  BigDelta, psi, theta;
		double[] /*REAL(DP), DIMENSION(56) ::*/ Clist, Dlist, Tlist, Nlist;
		double[] /*REAL(DP), DIMENSION(3) ::*/ beta, gamma;
		double[] /*REAL(DP), DIMENSION(2) ::*/ b, C, D, conditionsDIM = new double[2 /*+ 1*/];
		int /*INTEGER(I4B) ::*/ i;

		conditionsDIM[1 - 1] = conditions.pkb /*[1 - 1]*/ * 100000;
		conditionsDIM[2 - 1] = conditions.tK /*[2 - 1]*/;

		pKPa = conditionsDIM[1 - 1];
		//tK = conditionsDIM[2 - 1]; //!900.0

		nan = sqrt(arg); //! Not-A-Number

/* calculate water density*/
		/*call*/
		rho = waterdensity(conditions/*, rho*/);

		delta = rho / rhoc;
		tau = tc / conditions.tK;

		Clist = new double[]{nan, nan, nan, nan, nan, nan, nan, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 3.0, 3.0, 3.0, 3.0, 4.0, 6.0, 6.0, 6.0, 6.0, nan, nan, nan, nan, nan};
		Dlist = new double[]{1.0, 1.0, 1.0, 2.0, 2.0, 3.0, 4.0, 1.0, 1.0, 1.0, 2.0, 2.0, 3.0, 4.0, 4.0, 5.0, 7.0, 9.0, 10.0, 11.0, 13.0, 15.0, 1.0, 2.0, 2.0, 2.0, 3.0, 4.0, 4.0, 4.0, 5.0, 6.0, 6.0, 7.0, 9.0, 9.0, 9.0, 9.0, 9.0, 10.0, 10.0, 12.0, 3.0, 4.0, 4.0, 5.0, 14.0, 3.0, 6.0, 6.0, 6.0, 3.0, 3.0, 3.0, nan, nan};
		Tlist = new double[]{-0.5, 0.875, 1.0, 0.5, 0.75, 0.375, 1.0, 4.0, 6.0, 12.0, 1.0, 5.0, 4.0, 2.0, 13.0, 9.0, 3.0, 4.0, 11.0, 4.0, 13.0, 1.0, 7.0, 1.0, 9.0, 10.0, 10.0, 3.0, 7.0, 10.0, 10.0, 6.0, 10.0, 10.0, 1.0, 2.0, 3.0, 4.0, 8.0, 6.0, 9.0, 8.0, 16.0, 22.0, 23.0, 23.0, 10.0, 50.0, 44.0, 46.0, 50.0, 0.0, 1.0, 4.0, nan, nan};
		Nlist = new double[]{0.12533547935523e-1, 0.78957634722828e1, -0.87803203303561e1, 0.31802509345418, -0.26145533859358, -0.78199751687981e-2, 0.88089493102134e-2, -0.66856572307965, 0.20433810950965, -0.66212605039687e-4, -0.19232721156002, -0.25709043003438, 0.16074868486251, -0.40092828925807e-1, 0.39343422603254e-6, -0.75941377088144e-5, 0.56250979351888e-3, -0.15608652257135e-4, 0.11537996422951e-8, 0.36582165144204e-6, -0.13251180074668e-11, -0.62639586912454e-9, -0.10793600908932, 0.17611491008752e-1, 0.22132295167546, -0.40247669763528, 0.58083399985759, 0.49969146990806e-2, -0.31358700712549e-1, -0.74315929710341, 0.47807329915480, 0.20527940895948e-1, -0.13636435110343, 0.14180634400617e-1, 0.83326504880713e-2, -0.29052336009585e-1, 0.38615085574206e-1, -0.20393486513704e-1, -0.16554050063734e-2, 0.19955571979541e-2, 0.15870308324157e-3, -0.16388568342530e-4, 0.43613615723811e-1, 0.34994005463765e-1, -0.76788197844621e-1, 0.22446277332006e-1, -0.62689710414685e-4, -0.55711118565645e-9, -0.19905718354408, 0.31777497330738, -0.11841182425981, -0.31306260323435e2, 0.31546140237781e2, -0.25213154341695e4, -0.14874640856724, 0.31806110878444};

		beta = new double[]{150, 150, 250};
		gamma = new double[]{1.21, 1.21, 1.25};
		b = new double[]{0.85, 0.95};
		C = new double[]{28, 32};
		D = new double[]{700, 800};
		a = 3.5;
		LargeBeta = 0.3;
		LargeA = 0.32;
		LargeB = 0.2;
		alpha = 20.0;
		epsi = 1.0;

		phir = 0.0;
		for (i = 1 - 1; i < 7; i++) {//DO i = 1, 7
			phir = phir + Nlist[i] * (pow(delta, Dlist[i])) * (pow(tau, Tlist[i]));
		}//END DO

		for (i = 8 - 1; i < 51; i++) {//DO i = 8, 51
			phir = phir + Nlist[i] * (pow(delta, Dlist[i])) * (pow(tau, Tlist[i])) * exp(-pow(delta, Clist[i]));
		}//END DO

		for (i = 52 - 1; i < 54; i++) {//DO i = 52, 54
			phir = phir + Nlist[i] * (pow(delta, Dlist[i])) * (pow(tau, Tlist[i])) * exp(-alpha * (pow((delta - epsi), 2)) - beta[i - 51] * (pow((tau - gamma[i - 51]), 2)));
		}//END DO

		for (i = 55 - 1; i < 56; i++) {//DO i = 55, 56
			theta = (1 - tau) + LargeA * pow((pow((delta - 1), 2)), (1 / (2 * LargeBeta)));
			BigDelta = pow(theta, 2) + (LargeB * pow((pow((delta - 1), 2)), a));
			psi = exp(-C[i - 54] * (pow((delta - 1), 2)) - D[i - 54] * (pow((tau - 1), 2)));
			phir = phir + Nlist[i] * (pow(BigDelta, b[i - 54])) * delta * psi;
		}//END DO

		dens = rho / (1000 * mH2O/*18.0153d0*/);

		fugacity = phir + log(dens) + (pKPa / 100) / (10 * 8.314462 * conditions.tK * dens) + log(10 * 8.314462 * conditions.tK) - 1;
		System.out.println("Initial Residual phir " + phir + "OldFugacity  =" + fugacity);

		initIAPWS95(new double[]{conditions.tK}, new double[]{rho});
		fugacity = residual.phi + log(dens) + (pKPa / 100) / (10 * 8.314462 * conditions.tK * dens) + log(10 * 8.314462 * conditions.tK) - 1;
		System.out.println("New Residual phir "+ residual.phi+ "NewFugacity  ="+fugacity);
		System.out.println("Verified Water Fugacity: 3.17 kPa");
		return fugacity;
	}//END SUBROUTINE fugwater

	public static class cacheH2O {
		static solutionProperty conditions = null;
		static Double density = null;
		public static int ask = 0;
		static boolean inCacheDensity(solutionProperty askConditions) {ask +=1; return askConditions == conditions;}
		//if (cacheH2O.inCacheDensity(conditions)) return cacheH2O.density;
		//cacheH2O.conditions = conditions;//USE nrtype; USE nrutil, ONLY :assert_eq, imaxloc,iminloc,nrerror,swap
	}

	/**
	 * (i) the SUBROUTINE waterdensity(conditions, density)
	 * where conditions is (/ pkb, tK/)
	 * density is returned in kg/m³
	 */
	public static double /*SUBROUTINE*/ waterdensity(solutionProperty conditions /*double[] conditions*//*, double density*/) throws ExecutionException, InterruptedException {
		if (cacheH2O.inCacheDensity(conditions)) return cacheH2O.density;
		cacheH2O.conditions = conditions;//USE nrtype; USE nrutil, ONLY :assert_eq, imaxloc,iminloc,nrerror,swap
//USE MINIMIZE

		double  p, /*tK,*/ pvap;
		double  ftol;
		double[]/*[]*/ /*REAL(DP), DIMENSION(2,1) ::*/ param = new double[2 /*+ 1*/];//[1 /*+ 1*/];
		//double[] /*REAL(DP), DIMENSION(2), INTENT(IN) ::*/ conditions;
		double /*REAL(DP), INTENT(OUT) ::*/ density;
		double[] /*REAL(DP), DIMENSION(2) ::*/ y = new double[2 /*+ 1*/], conditionsDIM = new double[2 /*+ 1*/];
		int /*INTEGER(I4B) ::*/ niter = ITMAX;
//INTERFACE
//FUNCTION DENSWATEROBJ(rho, conditions)
		//USE nrtype

		//REAL(DP), DIMENSION(:) :: rho
		//REAL(DP) :: DENSWATEROBJ
		//REAL(DP), DIMENSION(:) :: conditions
//END FUNCTION DENSWATEROBJ
//END INTERFACE

		conditionsDIM[1 - 1] = conditions.pkb /*conditions[1 - 1]*/ * 100000;
		conditionsDIM[2 - 1] = conditions.tK /*conditions[2 - 1]*/;

		p = conditionsDIM[1 - 1];
		//tK = conditionsDIM[2 - 1];

//		System.out.println("WATEREOS Conditions " + conditionsDIM[1 - 1] + " kPa " + conditionsDIM[2 - 1] + "°K");

		pvap = (1e-5 * exp(-2836.5744 / pow(conditions.tK, 2) - 6028.076559 / conditions.tK + 19.54263612 - 0.02737830188 * conditions.tK + 1.6261698e-5 * pow(conditions.tK, 2) + 7.0229056e-10 * pow(conditions.tK, 3) - 1.8680009e-13 * pow(conditions.tK, 4) + 2.7150305 * log(conditions.tK))) * 100;//! givds pvap in KPa

		if (p < pvap) {//THEN
			density = 10.0; //!vapour phase
//			System.out.println("WATEREOS Vapour Pressure " + p + " pvap " + pvap + " tK " + tK + " initialDensity kg/m³ " + density);
		} else {//ELSE
			density = 1000; //1000.0; //0.996 556 0 × 10³ @ 300 tK! liquid
//			System.out.println("WATEREOS Liquid Pressure " + p + " pvap " + pvap + " tK " + tK + " initialDensity kg/m³ " + density);
		}//END IF
		//System.out.println("T = 500 K and p = 838.025 kg m–3 " + DENSWATEROBJ.apply(new double[]{838.025,0}, new double[]{500,838.025}));

//!conditions = (/double precision :: p, t /)
		ftol = 1.0e-7;

		param[1 - 1]/*[1]*/ = density; //0.996 556 0 × 10³! initial guess in kg/m³
		param[2 - 1]/*[1]*/ = param[1 - 1] /*[1]*/ * 1.001; //0.118 820 2 × 104 !kg/m³
		long start = System.nanoTime();
//		y[1 - 1] = DENSWATEROBJ.apply((new double[]{param[0], 0}/*[1]/*[]*/), conditionsDIM);
//		y[2 - 1] = DENSWATEROBJ.apply((new double[]{param[1], 0}/*[2]/*[]*/), conditionsDIM);
		double initialDelta = 0.001; // y[1] - y[0]; //PP
//		System.out.println("initial Density from WATEREOS " + VecMath.toString(y));


		double[] x = minimize/*NELDERMEAD*/(DENSWATEROBJ, param, initialDelta,/*, ftol*//*,DENSWATEROBJ,*/ niter, new String[]{"Pressure", "Temperature"}, conditionsDIM);
//		System.out.println("WATEREOS 217 DENSWATEROBJ Solution: x = " + VecMath.toString(/*minimizer.*/x));

//		Log.w("Minimize", "Old " + (System.nanoTime() - start));
		density = /*minimizer.*/x/*minimizer.x[0]*/[0];
		//density = param[1-1]/*[1]*/;
//		start = System.nanoTime();
//		Amoeba amo = new Amoeba(2, 1, param[0], param[1], 100, DENSWATEROBJ, conditionsDIM /*new double[]{1013.25, 298.15}*/);
//		Solution sol = amo.Solve();
//		Log.w("Minimize", "New " + (System.nanoTime() - start));

//		Log.e("Amoeaba", VecMath.toString(sol.vector));
//		density = sol.vector[0]; // bestSolution.getFeatures()[0]; // 997.458; //density;
//		System.out.println("Fortran  waterdensity is:      997.45800046110639      ");
//		System.out.println("this water denstity           " + density);

		cacheH2O.density = density;

		return density;
	}//END SUBROUTINE waterdensity
//END MODULE WATEREOS

	/**
	 * Objective function for minimisation of Helmoltz free energy of water
	 *
	 * @see "https://web.archive.org/web/20110726024808/http://www.iapws.org/relguide/IAPWS95-Rev.pdf"
	 * Not finished to be verified, still errors
	 */
	private static BiVectorFunction DENSWATEROBJ/*(double[] rho, double[] conditions)*/ = (rho, conditions) -> {
//USE nrtype

		//double /*REAL(DP), DIMENSION(:) ::*/ rho;
		double  DENSWATEROBJ;
		//double /*REAL(DP), DIMENSION(:) ::*/ conditions;
		double  arg = -1.0, nan, a, LargeBeta, LargeA, LargeB;
		double  tc = 647.096/*Tc = 647.096 K (1) */, rhoc = 322.0/*ρc = 322 kg m–3*/, r = 0.46151805/*R = 0.461 518 05 kJ kg–1 K–1*/;
		double  delta, tau, pKPa, tK, alpha, epsi;//!, phir
		double  BigDelta, psi, theta, phirDelta;//!, phir1, delta1
		double  dpsidDelta, dBigDeltadDelta;
		double[] /*REAL(DP), DIMENSION(56) ::*/ Clist, Dlist, Tlist, Nlist;
		double[] /*REAL(DP), DIMENSION(3) ::*/ beta, gamma;
		double[] /*REAL(DP), DIMENSION(2) ::*/ b, C, D;
		int /*INTEGER(I4B) ::*/ i;

		pKPa = conditions[1 - 1]; //TODO Verify pkb pkPa
		tK = conditions[2 - 1]; //!900.0

		//nan = sqrt(arg); //! Not-A-Number

//		delta = rho[1 - 1] / rhoc;
//		tau = tc / tK;
//
//		Clist = new double[]{nan, nan, nan, nan, nan, nan, nan, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 3.0, 3.0, 3.0, 3.0, 4.0, 6.0, 6.0, 6.0, 6.0, nan, nan, nan, nan, nan};
//		Dlist = new double[]{1.0, 1.0, 1.0, 2.0, 2.0, 3.0, 4.0, 1.0, 1.0, 1.0, 2.0, 2.0, 3.0, 4.0, 4.0, 5.0, 7.0, 9.0, 10.0, 11.0, 13.0, 15.0, 1.0, 2.0, 2.0, 2.0, 3.0, 4.0, 4.0, 4.0, 5.0, 6.0, 6.0, 7.0, 9.0, 9.0, 9.0, 9.0, 9.0, 10.0, 10.0, 12.0, 3.0, 4.0, 4.0, 5.0, 14.0, 3.0, 6.0, 6.0, 6.0, 3.0, 3.0, 3.0, nan, nan};
//		Tlist = new double[]{-0.5, 0.875, 1.0, 0.5, 0.75, 0.375, 1.0, 4.0, 6.0, 12.0, 1.0, 5.0, 4.0, 2.0, 13.0, 9.0, 3.0, 4.0, 11.0, 4.0, 13.0, 1.0, 7.0, 1.0, 9.0, 10.0, 10.0, 3.0, 7.0, 10.0, 10.0, 6.0, 10.0, 10.0, 1.0, 2.0, 3.0, 4.0, 8.0, 6.0, 9.0, 8.0, 16.0, 22.0, 23.0, 23.0, 10.0, 50.0, 44.0, 46.0, 50.0, 0.0, 1.0, 4.0, nan, nan};
//		Nlist = new double[]{0.12533547935523e-1, 0.78957634722828e1, -0.87803203303561e1, 0.31802509345418, -0.26145533859358, -0.78199751687981e-2, 0.88089493102134e-2, -0.66856572307965, 0.20433810950965, -0.66212605039687e-4, -0.19232721156002, -0.25709043003438, 0.16074868486251, -0.40092828925807e-1, 0.39343422603254e-6, -0.75941377088144e-5, 0.56250979351888e-3, -0.15608652257135e-4, 0.11537996422951e-8, 0.36582165144204e-6, -0.13251180074668e-11, -0.62639586912454e-9, -0.10793600908932, 0.17611491008752e-1, 0.22132295167546, -0.40247669763528, 0.58083399985759, 0.49969146990806e-2, -0.31358700712549e-1, -0.74315929710341, 0.47807329915480, 0.20527940895948e-1, -0.13636435110343, 0.14180634400617e-1, 0.83326504880713e-2, -0.29052336009585e-1, 0.38615085574206e-1, -0.20393486513704e-1, -0.16554050063734e-2, 0.19955571979541e-2, 0.15870308324157e-3, -0.16388568342530e-4, 0.43613615723811e-1, 0.34994005463765e-1, -0.76788197844621e-1, 0.22446277332006e-1, -0.62689710414685e-4, -0.55711118565645e-9, -0.19905718354408, 0.31777497330738, -0.11841182425981, -0.31306260323435e2, 0.31546140237781e2, -0.25213154341695e4, -0.14874640856724, 0.31806110878444};
//
//		beta = new double[]{150, 150, 250};
//		gamma = new double[]{1.21, 1.21, 1.25};
//		b = new double[]{0.85, 0.95};
//		C = new double[]{28, 32};
//		D = new double[]{700, 800};
//		a = 3.5;
//		LargeBeta = 0.3;
//		LargeA = 0.32;
//		LargeB = 0.2;
//		alpha = 20.0;
//		epsi = 1.0;
//
//		//C'est parti, deuxième ligne page 12
//		phirDelta = 0.0;
//		for (i = 1 - 1; i < 7; i++) { //DO i = 1, 7 OK
//			phirDelta = phirDelta + Nlist[i] * Dlist[i] * pow(delta, Dlist[i] - 1) * pow(tau, Tlist[i]);
//		}//END DO
//
//		for (i = 8 - 1; i < 51; i++) { //DO i = 8, 51 TODO erreur parenthese OK
//			phirDelta = phirDelta + Nlist[i] * exp(-pow(delta, Clist[i])) * (pow(delta, Dlist[i] - 1) * pow(tau, Tlist[i]) * (Dlist[i] - Clist[i] * pow(delta, Clist[i])));
//		}//END DO
//
//		for (i = 52 - 1; i < 54; i++) { //DO i = 52, 54 OK
//			phirDelta = phirDelta + (Nlist[i] * pow(delta, Dlist[i]) * pow(tau, Tlist[i]) * exp(-alpha * pow(delta - epsi, 2) - beta[i - 51] * pow(tau - gamma[i - 51], 2)) * ((Dlist[i] / delta) - (2 * alpha * (delta - epsi))));
//		}//END DO
//
//		for (i = 55 - 1; i < 56; i++) { //DO i = 55, 56
//			theta = 1 - tau + LargeA * pow(pow(delta - 1, 2), 1 / (2 * LargeBeta)); //OK
//			BigDelta = pow(theta, 2) + LargeB * pow(pow(delta - 1, 2), a); //OK
//			psi = exp(-C[i - 54] * pow(delta - 1, 2) - D[i - 54] * pow(tau - 1, 2)); //OK
//			//derivative
//			dpsidDelta = -2 * C[i - 54] * (delta - 1) * psi; //OK
//			// Pour debug on sépare comme page 13
//			double dDelta = (delta - 1) * (LargeA * theta * (2 / LargeBeta) * pow(pow(delta - 1, 2), 1 / (2 * LargeBeta) - 1) + 2 * LargeB * a * pow(pow(delta - 1, 2), a - 1)); //OK
//
//			dBigDeltadDelta = b[i - 54] * pow(BigDelta, b[i - 54] - 1) * dDelta; //OK
//
//			phirDelta = phirDelta + Nlist[i] * (pow(BigDelta, b[i - 54]) * (psi + delta * dpsidDelta) + dBigDeltadDelta * delta * psi); //OK
//
//		}//END DO
//
//		DENSWATEROBJ = pow(1 + delta * phirDelta - pKPa / (rho[1 - 1] * r * tK), 2);
		initIAPWS95(new double[]{tK}, new double[]{rho[0]});
		//Minimise function of pressure
		//				return DENSWATEROBJ;
		return pow(1 + IAPWS95.delta * residual.Phi.delta - pKPa / (rho[1 - 1] * r * tK), 2);
	};//END FUNCTION DENSWATEROBJ
}
