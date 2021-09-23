package net.priveyes.PhaseDiag.VISCO;

import static net.priveyes.PhaseDiag.CO2EOS.CO2EOS.CO2DENSITY;
import static net.priveyes.PhaseDiag.FluidMin.Constants.mCO2;
import static net.priveyes.PhaseDiag.FluidMin.Constants.mClminus;
import static net.priveyes.PhaseDiag.FluidMin.Constants.mH2O;
import static net.priveyes.PhaseDiag.FluidMin.Constants.mNaplus;
import static net.priveyes.PhaseDiag.IAPWS.WATEREOS.waterdensity;
import static java.lang.Math.exp;
import static java.lang.Math.pow;

import net.priveyes.PhaseDiag.CHNOZ.IAPWScommon.Rho;
import net.priveyes.PhaseDiag.CHNOZ.IAPWScommon.Temperature;
import net.priveyes.PhaseDiag.num_rep.VecMath;
import net.priveyes.PhaseDiag.FluidMin.solutionProperty;

import java.util.concurrent.ExecutionException;

public class visco {
	public static double Bando2004(solutionProperty conditions /*double pkb, double tK*/, double molalCO2/*, double molalNaCl*/) {
		double Bando2004;
//double  pkb, tK, molalCO2, molalNaCl
		double pMPa, massFracNaCl, xcs, xc;
		//double /*REAL(DP), PARAMETER ::*/ /*mCO2 = 44.009,*/ mH2O = 18.015;

		pMPa = 100.0 * conditions.pkb;
//PP NaCl 58.443 g/mol 22,98976928 ± 2,0E-8 g/mol
		massFracNaCl = (mNaplus + mClminus) * conditions.molalNaCl
				/ ((mNaplus + mClminus) * conditions.molalNaCl + mCO2 * molalCO2 + 1.0e3);

		xc = molalCO2 / (conditions.molalNaCl + molalCO2 + 1000.0 / mH2O);

		xcs = pMPa / (36.1 * pMPa + 3.87 * conditions.tK - 1097.1 + massFracNaCl
				* (196.0 * pMPa + 26.9 * conditions.tK - 8810.0));
		//Bando2004 = 1.0d0 + xc * (0.2531d0 - (tK - 273.15d0) * ((4.069d0)**(-3))) / xcs
		Bando2004 = 1.0 + xc * (0.2531 - (conditions.tK - 273.15) * (pow((4.069), -3))) / xcs;
		return Bando2004;
	}//END FUNCTION Bando2004

	/**
	 * This is the equation of Mao, S.  Duan, Z.INTERNATIONAL JOURNAL OF THERMOPHYSICS, 2009, 30, 1510 - 1523
	 * The range is 273Kto 573 K, 1 bar to 1000 bar, and 0 to 6 mol\[CenterDot]kg\[Minus]1
	 TODO ionic or molalNaCl ?? */
	public static double ratioVisco(solutionProperty conditions /*double tK*//*, double ionic*/) {
		double ratioVisco;
//double  tK, ionic
		double  a0 = -0.21319213, a1 = 0.13651589e-2, a2 = -0.12191756e-5,
				b0 = 0.69161945e-1, b1 = -0.27292263e-3, b2 = 0.20852448e-6,
				c0 = -0.25988855e-2, c1 = 0.77989227e-5;
		double Av, Bv, Cv;
		Av = a0 + a1 * conditions.tK + a2 * pow(conditions.tK, 2);
		Bv = b0 + b1 * conditions.tK + b2 * pow(conditions.tK, 2);
		Cv = c0 + c1 * conditions.tK;

		//ratioVisco = Exp(Av * ionic + Bv * ionic**2 + Cv * Ionic**3)
		ratioVisco = exp(Av * conditions.molalNaCl/*ionic*/ + Bv * pow(conditions.molalNaCl/*ionic*/, 2)
				+ Cv * pow(conditions.molalNaCl/*ionic*/, 3));
		return ratioVisco;
	}//END FUNCTION ratioVisco

	/**
	 * VISCOSITY OF WATER - Gives the viscosity of water in Pa.s from the IAPWS,
	 * not implemented for the small region around the critical point for which they give an enhancement equation
	 * - this equation is off by 5 % at T + - 3 K around the critical point.
	 * Note that 1 Pa.s = 1 kg.m-1.s-1 = 10 Poises
	 * tstar: K
	 * rhostar: kg/m^3
	 * pstar = 22.064 MPa unused
	 * mustar: Pa.s
	 */

	public static double viscoWater(solutionProperty conditions /*double pkb, double tK*/) throws ExecutionException, InterruptedException {

		//double  pkb, tK;
		double viscoWater;
		double densityWater, tKdim, sum1, mu0dim;
		double si, sj, mu1dim, mu2dim, rhodim;
		double /*REAL(DP), PARAMETER ::*/ /*tstar = 647.096,*/ /*rhostar = 322.0,*/ mustar = 1.0e-6;
		double[] /*REAL(DP), DIMENSION(4) ::*/ largeH;
		double[][] /*REAL(DP), DIMENSION(7,6) ::*/ largeHbis = new double[7 /*+ 1*/][6 /*+ 1*/];
		int /*INTEGER(I4B) ::*/ i, j;

		/*   where conditions is (/ pkb, tK/)
       density is returned in kg/m³
*/
		densityWater = waterdensity(conditions /*new double[]{pkb, tK}*//*, densityWater*/);

		largeH = new double[]{1.67752, 2.20462, 0.6366564, -0.241605};

		largeHbis[1 - 1]/*,:)*/ = new double[]{5.20094e-1, 8.50895e-2, -1.08374, -2.89555e-1, 0.0, 0.0};
		largeHbis[2 - 1]/*,:)*/ = new double[]{2.22531e-1, 9.99115e-1, 1.88797, 1.26613, 0.0, 1.20573e-1};
		largeHbis[3 - 1]/*,:)*/ = new double[]{-2.81378e-1, -9.06851e-1, -7.72479e-1, -4.89837e-1, -2.5704e-1, 0.0};
		largeHbis[4 - 1]/*,:)*/ = new double[]{1.61913e-1, 2.57399e-1, 0.0, 0.0, 0.0, 0.0};
		largeHbis[5 - 1]/*,:)*/ = new double[]{-3.25372e-2, 0.0, 0.0, 6.98452e-2, 0.0, 0.0};
		largeHbis[6 - 1]/*,:)*/ = new double[]{0.0, 0.0, 0.0, 0.0, 8.72102e-3, 0.0};
		largeHbis[7 - 1]/*,:)*/ = new double[]{0.0, 0.0, 0.0, -4.35673e-3, 0.0, -5.93264e-4};

		tKdim = conditions.tK / Temperature.critical/*tstar*/;

		sum1 = 0.0;

		for (i = 0; i <= 3; i++) { //DO i = 0, 3
			sum1 = sum1 + largeH[i/* + 1*/] / (pow(tKdim, i));
		} //END DO//!i

		mu0dim = 100.0 * Math.sqrt(tKdim) / sum1;

		rhodim = densityWater / Rho.critical/*rhostar*/;

		si = 0.0;
		for (i = 0; i <= 5; i++) { //DO i = 0, 5
			sj = 0.0;
			for (j = 0; j <= 6; j++) { //DO j = 0, 6
				sj = sj + largeHbis[j/* + 1*/][i/* + 1*/] * (pow((rhodim - 1), j));
			} //END DO//!j
			si = si + (pow(((1 / tKdim) - 1), i)) * sj;
		} //END DO//! i

		mu1dim = exp(rhodim * si);

		mu2dim = 1.0;

		viscoWater = mustar * (mu0dim * mu1dim * mu2dim);
		//System.out.println("Verified Water Viscosity: 890 microPa-s\t" + viscoWater);
		return viscoWater;
	}//END FUNCTION viscoWater

	/**
	 * gives result in \[Mu]Pa.s, critical enhancement neglected
	 * so less accurate by ~2% within 5K of the critical point of CO2
	 */
	public static double viscoCO2(solutionProperty conditions /*double pkb, double tK*/) throws ExecutionException, InterruptedException {
		//double  pkb, tK;
		double viscoCO2;
		double densityCO2;
		double[] /*REAL(DP), DIMENSION(5) ::*/ a;
		double tKstar, gt, Rho, Eta0t, DeltaEtaRhot;
		//double /*REAL(DP), PARAMETER ::*/ mCO2 = 44.009;
		int /*INTEGER(I4B) ::*/ i;

		densityCO2 = CO2DENSITY(conditions /*new double[]{pkb, tK}*//*, densityCO2*/);
//!       where conditions = (/pkb, tK/) returns density in mol/cm³

		a = new double[]{0.235156, -0.491266, 0.0521116, 0.0534791, -0.015371};

		tKstar = conditions.tK / 251.196;

		Rho = 1.0e3 * mCO2 * densityCO2;

		gt = 0.0;
		for (i = 0; i <= 4; i++) { //DO i = 0, 4
			gt = gt + a[i/* + 1*/] * (pow(Math.log(tKstar), i));
		}//END DO//!i
		gt = exp(gt);

		Eta0t = 1.00697 * (Math.sqrt(conditions.tK)) / gt;

		//DeltaEtaRhot = 0.004071119d0 * Rho + 0.00007198037d0 * (Rho**2) + &
		//		((2.4117d-17 * (Rho**6)) / (tKstar**3)) + &
		//		(2.97107d-23 - (1.62788d-23/tKstar) ) * (Rho**8)
		DeltaEtaRhot = 0.004071119 * Rho + 0.00007198037 * (pow(Rho, 2))
				+ ((2.4117e-17 * (pow(Rho, 6))) / (pow(tKstar, 3)))
				+ (2.97107e-23 - (1.62788e-23 / tKstar)) * (pow(Rho, 8));

		viscoCO2 = Eta0t + DeltaEtaRhot;
		//System.out.println("Verified CO2 Viscosity: 14.9 microPa-s\t" + viscoCO2);

		return viscoCO2;
	}//END FUNCTION viscoCO2
	public static double muH2O;
	public static void viscoH2ONaClCO2(double pkb, double tK, double xCO2) {

		//http://pubs.geothermal-library.org/lib/grc/1030394.pdf
		double[][] table_1 = new double[4][4];
		//i                         a      b              c               d
		table_1[0] = new double[]{9.03591045e+01, 0, 0, -1.22757462e-01};
		table_1[1] = new double[]{0, 3.40285740e+04, 1.40090092e-02, 2.15995021e-02};
		table_1[2] = new double[]{0, 8.23556123e+08, 4.86126399e-02, -3.65253919e-04};
		table_1[3] = new double[]{0, -9.28022905e+08, 5.26696663e-02, 1.97270835e-06};

		double ai = VecMath.sum(VecMath.getcolumn(table_1, 0));
		double[] bi = VecMath.getcolumn(table_1, 1);
		double[] ciT = VecMath.vxs(VecMath.getcolumn(table_1, 2), tK);
		double[] di = VecMath.getcolumn(table_1, 3);

		//µH2O = a0 +(i=1 3)∑bi . exp(−ci . T ) + P . (i=0 3)∑di . (T − 293.15)^i
		muH2O = 0;
		for (int i = 0; i < 4; i++) {
			muH2O += bi[i] * exp(-ciT[i]) + pkb * di[i] * pow((tK - 293.15), i);
		}
		muH2O += ai;
		System.out.println("New Water Viscosity: " + muH2O);

		ai = 1.34136579e+02;
		bi = VecMath.zerovec(4);
		bi[1] = -4.07743800e+03;
		bi[2] = 1.63192756e+04;
		bi[3] = 1.37091355e+03;
		double[] ci = VecMath.zerovec(4);
		ci[1] = -5.56126409e-03;
		ci[2] = -1.07149234e-02;
		ci[3] = -5.46294495e-04;
		di = VecMath.zerovec(4);
		di[1] = 4.45861703e-01;
		di[2] = -4.51029739e-04;
//ρH2O = a0 + (i=1 3)∑ bi . 10^(ci . T) + (i=1 2)∑di . P^i  (2)
		ciT = VecMath.vxs(ci, tK);

		double rhoH2ONaCl = 0;
		for (int i = 0; i < 4; i++) {
			rhoH2ONaCl += bi[i] * pow(10, ciT[i]) + di[i] * pow(pkb, i);
		}
		rhoH2ONaCl += ai;
		System.out.println("New WaterNaCl Density: " + rhoH2ONaCl);

		//La suite utilise mu WaterNaCl mais le papier ne donne pas la formule pour mu
		// Il faut utiliser celle de :
		// Mao, S., and Z. Duan, 2009. “The Viscosity of aqueous Alkali-Chloride
		// http://www.cugb.edu.cn/upload/20600/papers_upload/news_20101129162026.pdf
		// Avec cette la nouvelle densité à la place

//Marche pas
//		di = new double[] {
//		0.28853170e7, //1-0
//		-0.11072577e5, //2-1
//		-0.90834095e1, //3-2
//		0.30925651e-1, //4-3
//		-0.27407100e-4, //5-4
//		-0.19283851e7, //6-5
//		0.56216046e4, //7-6
//		0.13827250e2, //8-7
//		-0.47609523e-1, //9-8
//		0.35545041e-4}; //10-9
//
//		// ln ηH2O =(i=1 5)∑di T^(i−3) + (i=6 10)∑di ρH2O T^(i−8)
//		double sum1 = 0, sum2 = 0;
//		for (int i = 0; i < 5; i++) sum1 += di[i] * pow(tK, (i-3));
//		for (int i = 5; i < 10; i++) sum2 += di[i] * rhoH2ONaCl * pow(tK, (i-8));
//
//		System.out.println("Other Water Viscosity: " + exp(sum1 + sum2));

//		//Table 3. Coefficients of Eq. 3.
//		bi = VecMath.zerovec(4);
//		bi[0] = 7.632609119e+02;
//		bi[1] = - 9.46077673946e+03;
//		di = VecMath.zerovec(4);
//		di[0] = -1.047187396332e+04;
//		di[1] = 3.68325597e+01; //org with error 3.6.8325597e+01;
//
//		//µr = 1 + [ ((i=1 2)∑ai . xCO2^i) / ( (i=0 1)∑bi . T^i) ]
//		//µH2O+CO2 = µr * µH2O
//		double mur = 0;
//		for (int i = 0; i < 2; i++) {
//			mur += (bi[i] * pow(xCO2, i)) / (di[i] * pow(tK, i));
//		}
//		mur += 1; mur *= muH2O;
//		System.out.println("New H2O+CO2 Viscosity: " + mur);
//
//		//µH2O+NaCl+CO2 = µH2O+NaCl . (1 + 4.65 . xCO2^1.0134)
//		double muHNC =  rhoH2ONaCl * (1 + 4.65 * pow(xCO2, 1.0134));
//		System.out.println("New H2ONaClCO2 Viscosity: " + muHNC);
	}

	/**
	 * where ηr denotes the relative viscosity; ηsol refers to the viscosity (Pa ·s) of solutions;
	 * ηH2O is the viscosity of pure water in Pa ·s; and m is the molality (mol · kg−1) of salts
	 * (LiCl, NaCl, or KCl). A, B, and C are polynomial functions of temperature T (in K):
	 */
	public static void viscoSystem(double tK, double molality) {
		//Table 3 Parameters of Eqs. 3–5
		//Parameters Systems
		//0) LiCl–H2O 1) NaCl–H2O 2) KCl–H2O
		double[][] system = new double[8][3];
		system[0]/*a0*/ = new double[]{0.62204136e-2, -0.21319213, -0.42122934};
		system[1]/*a1*/ = new double[]{0.54436974e-3, 0.13651589e-2, 0.18286059e-2};
		system[2]/*a2*/ = new double[]{-0.40443190e-6, -0.12191756e-5, -0.13603098e-5};
		system[3]/*b0*/ = new double[]{0.14987325e-1, 0.69161945e-1, 0.11380205e-1};
		system[4]/*b1*/ = new double[]{-0.66617390e-4, -0.27292263e-3, 0.47541391e-5};
		system[5]/*b2*/ = new double[]{0.52113332e-7, 0.20852448e-6, -0.99280575e-7};
		system[6]/*c0*/ = new double[]{0.12101624e-5, -0.25988855e-2, 0};
		system[7]/*c1*/ = new double[]{0.17772678e-6, 0.77989227e-5, 0};

		double[] NaClH2O = VecMath.getcolumn(system, 1);

		double A = NaClH2O[0/*a0*/] + NaClH2O[1/*a1*/] * tK + NaClH2O[2/*a2*/] * pow(tK, 2);
		double B = NaClH2O[3/*b0*/] + NaClH2O[4/*b1*/] * tK + NaClH2O[5/*b2*/] * pow(tK, 2);
		double C = NaClH2O[6/*c0*/] + NaClH2O[7/*c1*/] * tK;

		//ln ηr = Am + Bm2 + Cm3
		double viscR = exp(A * molality + B * pow(molality, 2) + C * pow(molality, 3));
		//ηr = ηsol / ηH2O
		double viscSol = viscR * muH2O;

		System.out.println("New Solution Viscosity: " + viscSol);
	}
}//END MODULE VISCO
