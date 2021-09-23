package net.priveyes.PhaseDiag.THERMO;

import static net.priveyes.PhaseDiag.CO2EOS.CO2EOS.fugCO2;
import static net.priveyes.PhaseDiag.CO2EOS.CO2EOS.volumeCO2;
import static net.priveyes.PhaseDiag.EVALW.EVALW.returnWH2OCO2.WH2OCO2;
import static net.priveyes.PhaseDiag.EVALW.EVALW.returnWH2OCO2.aCO2;
import static net.priveyes.PhaseDiag.EVALW.EVALW.returnWH2OCO2.aH2O;
import static net.priveyes.PhaseDiag.EVALW.EVALW.returnWH2OCO2.wWCO2;
import static net.priveyes.PhaseDiag.EVALW.EVALW.returnWH2ONACL.aNaClpm;
import static net.priveyes.PhaseDiag.EVALW.EVALW.returnWH2ONACL.wWNaClpm;
import static net.priveyes.PhaseDiag.FluidMin.Constants.Rm;
import static net.priveyes.PhaseDiag.FluidMin.Constants.avo;
import static net.priveyes.PhaseDiag.FluidMin.Constants.bolts;
import static net.priveyes.PhaseDiag.FluidMin.Constants.kelsq;
import static net.priveyes.PhaseDiag.FluidMin.Constants.mH2O;
import static net.priveyes.PhaseDiag.FluidMin.solutionProperty.nH2Oapp;
import static net.priveyes.PhaseDiag.IAPWS.WATEREOS.fugwater;
import static net.priveyes.PhaseDiag.IAPWS.WATEREOS.volumewater;
import static net.priveyes.PhaseDiag.IAPWS.WATEREOS.waterdensity;
import static java.lang.Math.pow;

import android.util.Log;

import net.priveyes.PhaseDiag.CHNOZ.IAPWScommon.Rho;
import net.priveyes.PhaseDiag.CHNOZ.IAPWScommon.Temperature;
import net.priveyes.PhaseDiag.num_rep.VecMath;
import net.priveyes.PhaseDiag.FluidMin.solutionProperty;

import java.util.concurrent.ExecutionException;

/**
 * //! This module defines:
 * //! (i)   FUNCTION VOLUME(pkb, tK, xCO2)
 * //! (ii)  FUNCTION vex(pkb, tK, xCO2)
 * //! (iii) FUNCTION gex(pkb, tK, xCO2)
 * //! (iv)  FUNCTION awgf(pkb, tK, xco2)
 * //! (v)   FUNCTION dielectric(pkb, tK, xco2)
 *
 * //! (vi)  Function axfunction(pkb, tK, xco2)
 * //! (vii) Function bxfunction(pkb, tK, xco2)
 * //! (iix) FUNCTION debyehuckelsolute(pkb, tK, xco2, ionic, za, zc)
 * //! (ix)  FUNCTION debyehuckelsolvent(pkb, tK, xco2, ionic)
 * //! (x)   FUNCTION RTlnGamASF(prop,Ws,alphas,l)
 * //! (xi)  FUNCTION CalcMu(pkb, tK, properties)
 * //! (xii) FUNCTION AQGIBBS(pkb, tK, properties, wCO2NaCl0, alphaAq)
 * //!
 * //! Benoit.Dubacq@upmc.fr 15-01-2014
 */
public class  thermo {

	//PUBLIC :: VOLUME, vex, gex, awgf, dielectric, axfunction, bxfunction, debyehuckelsolute, debyehuckelsolvent, RTlnGamASF, CalcMu, AQGIBBS
	public final static double[] propNaCl0 = new double[]{-393.42, 79.7, 2.965, 0.072, -0.3223, 0, 0, 50, 66, 0, 0, 0, 0, 0, 1}; //! Halite liquid from EP07
	public final static double[] propNaClpm = new double[]{-203.69, 57.565, 0.834, -0.0417156, -0.000208763, 0, 0, 0, 0, 0, 0, 0, 1, 1, 2}; //! NaClpm from Helgeson et al with refitted heat capacity
	public final static double[] propHalite = new double[]{-411.30, 72.10, 2.7020, 0.0452, 1.797, 0, 0, 26.90, 240, 0, 0, 0, 0, 0, 0}; //!halite from RH95, uncertainty from THERMOCALC

	//public static final double /*REAL(DP), PARAMETER ::*/ R = 8.314462e-3;
	//public static final double avo = 6.02214e23; //! (*Avogadros number, moles -1*)
	//public static final double kelsq = 2.307e-28; //! (*value of unit electric charge, in force units, kg m3 s-2*)
	//public static final double bolts = 1.381e-23; //! (*boltsmanns constant, kg m2 s-1 K-1*)

	static double volume(solutionProperty conditions /*double pkb, double tK*/, double xCO2) throws ExecutionException, InterruptedException {

		//double /*REAL(DP)::*/ pkb, tK, xCO2;
		double /*REAL(DP)::*/ VOLUME, vH2O, vCO2, vxs;

		vH2O = volumewater(conditions /*new double[]{pkb, tK}*//*, vH2O*/); //! vH2O in J/bar

		if (xCO2 == 0) {

			vCO2 = 0.0;
			vxs = 0.0;

		} else {

			vCO2 = volumeCO2(conditions /*new double[]{pkb, tK}*//*,vCO2*/);
			vxs = vex(conditions /*pkb, tK*/, xCO2);

		}//END IF

		VOLUME = (1.0 - xCO2) * vH2O + xCO2 * vCO2 + vxs;
		return VOLUME;
	}//END FUNCTION VOLUME

	/**numerical derivative, delta = 1 bar, vex is in kJ/kbar*/
	private static double vex(solutionProperty conditions /*double pkb, double tK*/, double xCO2) throws ExecutionException, InterruptedException {

		double /*REAL(DP)::*/ vex/*, pkb, tK, xCO2*/;

		vex = (gex(new solutionProperty(conditions.pkb + 0.001, conditions.tK, conditions.molalNaCl), xCO2) - gex(conditions /*pkb, tK*/, xCO2)) / 0.001;
		return vex;
	}//END FUNCTION vex

	private static double gex(solutionProperty conditions /*double pkb, double tK*/, double xCO2) throws ExecutionException, InterruptedException {

		double /*REAL(DP)::*/ gex/*, pkb, tK, xCO2*/;
		double  aA, aB, wAB;

		// xCO2 is nCO2 / nH2O

		/* returns wCO2NaCl0,aH2O,aCO2*/
		WH2OCO2(conditions /*pkb, tK*/,/*,wWCO2,aH2O,aCO2*/wWCO2, aH2O, aCO2);

		gex = 2. * aH2O * aCO2 * wWCO2 * xCO2 * (1.0 - xCO2)
				/ ((aH2O + aCO2) * (aH2O * (1.0 - xCO2)
				+ aCO2 * xCO2));
		return gex;
	}//END FUNCTION gex

	/**archer and wang (1990) g function
	 * https://www.nist.gov/sites/default/files/documents/srd/jpcrd383.pdf page 5*/
	private static double  awgf(solutionProperty conditions /*double pkb, double tK*/, double xco2) throws ExecutionException, InterruptedException {

		double /*REAL(DP)::*/ awgf/*, pkb, tK, xco2*/;
		double  vol;

		vol = volume(conditions /*pkb, tK*/, xco2);

		awgf = 1.0 + (((1.0 - xco2) * 1.8015 + xco2 * 4.401) / vol)
				* (-4.044525 * conditions.pkb / conditions.tK + 103.618 / Math.sqrt(conditions.tK) + 75.32165
				/ (conditions.tK - 215) - 23.23778 / Math.sqrt(conditions.tK - 215.0) - 3.548184
				/ Math.sqrt(Math.sqrt(conditions.tK - 215.0)) + Math.exp(-1246.311 / conditions.tK + 263307.7
				/ Math.pow(conditions.tK, 2) - 69.28953 * conditions.pkb / conditions.tK - 20444.73 * conditions.pkb / Math.pow(conditions.tK, 2)));
		return awgf;
	}//END FUNCTION awgf

	static double calcIonic(solutionProperty conditions) {
		if (cacheThermo.inCacheIonic(conditions)) return cacheThermo.Ionic;
		cacheThermo.conditionsIonic = conditions;

		//! max is 0.1 / nH2Oapp
		double ionic = 0.1 / nH2Oapp;
		if (conditions.molalNaCl < 0.1) {
			ionic = conditions.molalNaCl / nH2Oapp;
		}
		//else {//ELSE
		//	ionic = 0.1 / nH2Oapp;
		//} //END IF
		cacheThermo.Ionic = ionic;
		return cacheThermo.Ionic;
	}

	public static void calcAllThermo(solutionProperty conditions) throws ExecutionException, InterruptedException {
		/*ionic = */thermo.calcIonic(conditions);

		/*DHNaClpm = */debyehuckelsolute(conditions /*pkb, tK*/, 0.0, cacheThermo.Ionic/* ionic*/, 1.0, 1.0); //! USE ONLY WITH IONIC FORCE <= 0.1/55.556 //!//!
//							System.out.println("debyehuckelsolute fortran is:     -0.10669606044223041      ");
//							System.out.println("this                      "+DHNaClpm);
		/*mu0NaClpm = */AQGIBBS(conditions /*pkb, tK*/, propNaClpm, wWNaClpm, aNaClpm);
//							System.out.println("AQGIBBS fortran is:       -205.30571258344008      ");
//							System.out.println("this              "+mu0NaClpm);
		/*mu0NaCl0 = */CalcMu(conditions /*pkb, tK*/, propNaCl0);
//							System.out.println("CalcMu fortran is:       -417.15252696176361      ");
//							System.out.println("this             "+mu0NaCl0);

	}

	public static class cacheThermo {
		static solutionProperty conditions = null;
		static solutionProperty conditionsGibbs = null;
		static solutionProperty conditionsDH = null;
		static solutionProperty conditionsMu = null;
		static solutionProperty conditionsIonic = null;

		public static Double Ionic = null;
		public static Double Dielectric = null;
		public static Double Gibbs = null;
		public static Double DH = null;
		public static Double Mu = null;

		static int callRTlnGamASFin = 0;
		static int callRTlnGamASFno = 0;
		static double[] conditionsRt = null;
		public static void printCalled() { Log.w("Cache"," RTlnGamASF hit "+callRTlnGamASFin +" RTlnGamASF miss "+callRTlnGamASFno ); }

		public static int ask = 0;
		static boolean inCacheDielectric(solutionProperty askConditions) {ask +=1;return askConditions == conditions;}
		static boolean inCacheGibbs(solutionProperty askConditions) {ask +=1;return askConditions == conditionsGibbs;}
		static boolean inCacheDebyeHuckel(solutionProperty askConditions) {ask +=1;return askConditions == conditionsDH;}
		static boolean inCacheMu(solutionProperty askConditions) {ask +=1;return askConditions == conditionsMu;}
		static boolean inCacheRt(double[] askConditions) {return VecMath.equalsExactly(askConditions, conditionsRt);}
		static boolean inCacheIonic(solutionProperty askConditions) {return askConditions == conditionsIonic;}
	}

	public static double  diecIAPWS97(solutionProperty conditions /*double pkb, double tK*/) throws ExecutionException, InterruptedException {
		//double  pkb, tK;
		double  diecIAPWS97;
		double  densityWater, eps0, rhoc, rho, g, a, b;
		final double muDi = 6.138e-30, alphaMean = 1.636e-40
				, k = 1.380658e-23, /*avo = 6.0221367e23,*/ mw = 0.018015268/*, tc = 647.096*/;
		double[] /*REAL(DP), DIMENSION(12) ::*/ n;
		double[] /*REAL(DP), DIMENSION(11) ::*/ i, j;
		int /*INTEGER(I4B) ::*/ h;

		// permittivity of free space C^2.J^-1.m^-1
		eps0 = 1.0 / (4e-7 * Math.PI * (pow((299792458.0), 2)));
		rhoc = Rho.critical/*322.0*/ / mw/*0.018015268*/; //!mol/m^3
//!alphaMean = 1.636e-40 //!mean molecular polarizability C^2.J^-1.m^-2
//!muDi = 6.138e-30 C.m Molecular dipole moment
//!k = 1.380658e-23//! J/K Boltzmann's constant
//!avo = 6.0221367e23//! mol^-1 Avogadro
//!mw = 0.018015268//!kg/mol molar mass of water
//!tc = 647.096(*critical T K*),
		n = new double[]{0.978224486826, -0.957771379375, 0.237511794148, 0.714692244396, -0.298217036956, -0.108863472196, 0.949327488264e-1, -0.980469816509e-2, 0.165167634970e-4, 0.937359795772e-4, -0.123179218720e-9, 0.196096504426e-2};
		i = new double[]{1.0, 1.0, 1.0, 2.0, 3.0, 3.0, 4.0, 5.0, 6.0, 7.0, 10.0};
		j = new double[]{0.25, 1.0, 2.5, 1.5, 1.5, 2.5, 2.0, 2.0, 5.0, 0.5, 10.0};

		densityWater = waterdensity(conditions /*new double[]{pkb, tK}*//*, densityWater*/);
//!   where conditions is (/ pkb, tK/)
//!       density is returned in kg/m³

		rho = densityWater / mw;

		g = 1.0 + n[12-1] * (rho / rhoc) * pow((conditions.tK / 228.0 - 1.0), -1.2);
		for (h = 0; h < 11; h++) { //DO h = 1, 11
			g = g + n[h] * (pow((rho / rhoc), i[h])) * (pow((Temperature.critical / conditions.tK), j[h]));
		} //END DO

		a = avo * (pow(muDi, 2)) * rho * g / (eps0 * k * conditions.tK);
		b = avo * alphaMean * rho / (3 * eps0);

		diecIAPWS97 = (1.0 + a + 5.0 * b + Math.sqrt(9.0 + 2.0 * a + 18.0 * b + pow(a, 2)
				+ 10.0 * a * b + 9.0 * pow(b, 2))) / (4.0 - 4.0 * b);
		Log.w("Visco", "DieletricIAPWS97 "+ diecIAPWS97);
		return diecIAPWS97;
	}//END FUNCTION diecIAPWS97

	/**akinfiev and zotov (1999) dielectric constant*/
	private static double  dielectric(solutionProperty conditions /*double pkb, double tK*/, double xco2) throws ExecutionException, InterruptedException {
		if (cacheThermo.inCacheDielectric(conditions)) return cacheThermo.Dielectric;
		cacheThermo.conditions = conditions;

//		double  dielectric/*, pkb, tK, xco2*/;
//		double  vol, g, b;

//		vol = volume(pkb, tK, xco2);
//		g = awgf(pkb, tK, xco2);

//		b = 1.0 + (22.7024 / vol) * (xco2 * 0.265 + (1.0 - xco2) * (0.147 + 835.25 * g / tK));

//		dielectric = (b + Math.sqrt(8.0 + Math.pow(b, 2))) / 4.0;

//		Log.w("Thermo", "Dieletric "+ dielectric);

//		Log.w("Water_AW90", "Dieletric " + UtilWater.water_AW90(/*double Temp=*/tK, /*final double rho=*/997.4586327973304, /*double Pres=*/1.01325));

		cacheThermo.Dielectric = diecIAPWS97(conditions /*pkb, tK*/)/*dielectric*/;
		return cacheThermo.Dielectric;
	}//END FUNCTION dielectric

	private static double axfunction(solutionProperty conditions /*double pkb, double tK*/, double xco2) throws ExecutionException, InterruptedException {

		double  axfunction/*, pkb, tK, xco2*/;
		double  /*avo, kelsq, bolts,*/ e, v;

		//avo = 6.02214e23; //! (*Avogadros number, moles -1*)
		//kelsq = 2.307e-28; //! (*value of unit electric charge, in force units, kg m3 s-2*)
		//bolts = 1.381e-23; //! (*boltsmanns constant, kg m2 s-1 K-1*)
		e = dielectric(conditions /*pkb, tK*/, xco2); //! (*dielectric constant, dimensionless*)
		v = 1e-5 * volume(conditions /*pkb, tK*/, xco2); //! (*volume in m3 mole-1*)

// A = e^2 * B / (2.303 * 8 * pi * epsilon0 * epsilonRho * k * T')
		axfunction = (Math.sqrt(8 * Math.PI * avo) * (Math.pow(kelsq, 1.5)))
				/ (2 * Math.log(10)/*2.303*//*4.606*/ * (Math.pow(bolts, 1.5)) * Math.sqrt(v) * (Math.pow(e, 1.5))
				* (Math.pow(conditions.tK, 1.5)));
		return axfunction;
	}//END FUNCTION axfunction

	private static double bxfunction(solutionProperty conditions /*double pkb, double tK*/, double xco2) throws ExecutionException, InterruptedException {

		double  bxfunction/*, pkb, tK, xco2*/;
		double  /*avo, kelsq, bolts,*/ e, v;//!mh2o, mco2,

		//avo = 6.02214e23; //! (*Avogadros number, moles -1*)
		//kelsq = 2.307e-28; //! (*value of unit electric charge, in force units, kg m3 s-2*)
		//bolts = 1.381e-23; //! (*boltsmanns constant, kg m2 s-1 K-1*)
		e = dielectric(conditions /*pkb, tK*/, xco2); //! (*dielectric constant, dimensionless*)
		v = 1e-5 * volume(conditions /*pkb, tK*/, xco2); //! (*volume in m3 mole-1*)
// B = sqrt(2 * e^2 * N / (epsilon0 * epsilonRho * k * T)
		bxfunction = Math.sqrt(8 * Math.PI * avo * kelsq) / (Math.sqrt(bolts * v * e * conditions.tK) * 1e10);
		return bxfunction;
	}//END FUNCTION bxfunction

	public static double debyehuckelsolute(solutionProperty conditions /*double pkb, double tK*/, double xco2, double ionic, double za, double zc) throws ExecutionException, InterruptedException {
		//if (cacheThermo.inCacheDebyeHuckel(conditions) return cacheThermo.DH;
		cacheThermo.conditionsDH = conditions;

		double  debyehuckelsolute/*, pkb, tK, xco2, ionic, za, zc*/;
		double  ax, bx;
		//! USE ONLY IF IONIC FORCE IS <= 0.1/55.556 !
		ax = axfunction(conditions /*pkb, tK*/, xco2);
		bx = bxfunction(conditions /*pkb, tK*/, xco2);
		// -A * z^2 * sqrt(ionic) / (1 + B * a0 * sqrt(I) zi = charge of specie
		// A = e^2 * B / (2.303 * 8 * pi * epsilon0 * epsilonRho * k * T')
		// B = sqrt(2 * e^2 * N / (epsilon0 * epsilonRho * k * T)
		double ionsize = 4.5;
		debyehuckelsolute = -ax * za * zc * Math.sqrt(ionic) / (1.0 + bx * ionsize * Math.sqrt(ionic));

		cacheThermo.DH = debyehuckelsolute;
		return debyehuckelsolute;
	}//END FUNCTION debyehuckelsolute

	public static double  debyehuckelsolvent(solutionProperty conditions /*double pkb, double tK*/, double xco2, double ionic) throws ExecutionException, InterruptedException {

		double  debyehuckelsolvent/*, pkb, tK, xco2, ionic*/;
		double  v, /*avo, kelsq, bolts,*/ rad, e, kappa, sigma;

		//avo = 6.02214e23; //! (*Avogadros number, moles -1*)
		//kelsq = 2.307e-28; //! (*value of unit electric charge, in force units, kg m3 s-2*)
		//bolts = 1.381e-23; //! (*boltsmanns constant, kg m2 s-1 K-1*)
		rad = 4.5e-10;
		e = dielectric(conditions /*pkb, tK*/, xco2); //! (*dielectric constant, dimensionless*)
		v = 1e-5 * volume(conditions /*pkb, tK*/, xco2); //! (*volume in m3 mole-1*)

		kappa = Math.sqrt(8.0 * avo * Math.PI * kelsq * ionic / (bolts * v * e * conditions.tK));
		sigma = (3.0 / (Math.pow((kappa * rad), 3.0))) * (1 + kappa * rad - (1 / (1 + kappa * rad)) - 2 * Math.log(1 + kappa * rad));

		debyehuckelsolvent = v * Math.pow(kappa, 3) * sigma / (8 * Math.PI * avo);
		return debyehuckelsolvent;
	}//END FUNCTION debyehuckelsolvent

	public static double  RTlnGamASF(double[] prop, double[][] Ws, double[] alphas, int ix) {
		if (cacheThermo.inCacheRt(prop)) cacheThermo.callRTlnGamASFin++; else cacheThermo.callRTlnGamASFno++;
		cacheThermo.conditionsRt = prop;
		//double /*REAL(DP), DIMENSION(:) ::*/ prop, alphas;
		//double /*REAL(DP), DIMENSION(:,:) ::*/ Ws;
//int /*INTEGER(I4B) ::*/ l;
		double  RTlnGamASF;
		int /*INTEGER(I4B) ::*/ i, j;
		double[] /*REAL(DP), DIMENSION(size(prop)) ::*/ filist = new double[prop.length];
		double  tempact, Wstar, fisum;
//!for end-member l (integer)

		if (prop[ix] == 0) {
			RTlnGamASF = 0;
		} else {//ELSE
//! First calculate phis

//!   filist = Table[0, {i, Length[prop]}];
			fisum = 0.0;
			for (i = 0; i < prop.length; i++) { //DO i = 1, size(prop)
				fisum = fisum + prop[i] * alphas[i];
			}//END DO

			for (i = 0; i < prop.length; i++) { //DO i = 1, size(prop) //!i <=  Length[prop], i++,
				filist[i] = prop[i] * alphas[i] / fisum;
			}//END DO

/*Then calculate the sum. */
/*    NB: tempact == +RTlnGam because we have qj = phi and not qj = -phi */

			tempact = 0.0;
			for (i = 0; i < prop.length - 1; i++) { //DO i = 1, size(prop)-1;
				for (j = i + 1; j < prop.length; j++) { //DO j = i+1, size(prop)//!     For[j = i + 1, j <= Length[prop], j++,
					Wstar = 2.0 * Ws[i][j] * alphas[ix] / (alphas[i] + alphas[j]);

					if (i == ix) tempact = tempact + (1.0 - filist[i]) * filist[j] * Wstar;
					if (j == ix) tempact = tempact + (1.0 - filist[j]) * filist[i] * Wstar;
					if (i !=/*/=*/ ix &&/*.AND.*/ j !=/*/=*/ ix) {
						tempact = tempact - filist[i] * filist[j] * Wstar;
					}
				}//END DO //!j
			}//END DO //!i

			RTlnGamASF = tempact;

		}//END IF
		return RTlnGamASF;
	}//END FUNCTION RTlnGamASF

	public static double CalcMu(solutionProperty conditions /*double pkb, double tK*/, double[] properties) throws ExecutionException, InterruptedException {
//! HP98
//		if (cacheThermo.inCacheMu(new double[]{ pkb, tK})) return cacheThermo.Mu;
		cacheThermo.conditionsMu = conditions;

		double  CalcMu;
		//double  pkb, tK;
		//double /*REAL(DP), DIMENSION(12) ::*/ properties;
		double /*REAL(DP), PARAMETER ::*/ af = 20.0, dkdt = -0.00015, T0 = 298.15;
		double  kappa, tc0, q20, tc, q2, entropy, cpa, cpb, cpc, cpd;
		double  alpha, kappa298, landaut, smax, vmax, hp20, gland, addG;
		double  enthalpy, volume, add, lnfc, lnfh;
		//double /*REAL(DP), PARAMETER ::*/ R = 8.314462e-3;

		enthalpy = properties[1-1];
		entropy = properties[2-1] / 1000.0;
		volume = properties[3-1];
		cpa = properties[4-1];
		cpb = properties[5-1] * 1e-5;
		cpc = properties[6-1];
		cpd = properties[7-1];
		alpha = properties[8-1] * 1e-5;
		kappa298 = properties[9-1];
		landaut = properties[10-1];
		smax = properties[11-1] / 1000.0;
		vmax = properties[12-1];

		kappa = kappa298 * (1 + dkdt * (conditions.tK - T0));
		tc0 = landaut;
		if (vmax == 0) {
			add = 0;
		} else { //ELSE
			add = conditions.pkb * vmax / smax;
		}//END IF
		tc = tc0 + add;

		if (smax == 0) {
			q20 = 0;
		} else { //ELSE
			q20 = Math.sqrt(1 - T0 / tc0);
		}//END IF

		if (conditions.tK > tc) {
			q2 = 0;
		} else { //ELSE
			q2 = Math.sqrt(1 - conditions.tK / tc);
		}//END IF

		hp20 = smax * q20 * tc0 * (1 - (Math.pow(q20, 2)) / 3);
		gland = smax * q2 * ((conditions.tK - tc) + (tc * (Math.pow(q2, 2))) / 3);

		if (volume == 0) {  //! water or CO2
			if (enthalpy == -241.81) {  //!WATER
				/*Call*/
				lnfh = fugwater(conditions /*new double[]{pkb, tK}*//*, lnfh*/);
				addG = Rm * conditions.tK * lnfh;
			} else { //ELSE //!CO2
				/*Call*/
				lnfc = fugCO2(conditions /*new double[]{pkb, tK}*//*, lnfc*/);
				addG = Rm * conditions.tK * lnfc;
			} //END IF
		} else {//ELSE
			addG = (volume + vmax * q20) * (1.0 + alpha * ((conditions.tK - T0) - af * (Math.sqrt(conditions.tK) - Math.sqrt(T0)))) * (Math.pow((1 + 4 * conditions.pkb / kappa), 3. / 4.) - 1) * kappa / 3;
		} //END IF

		CalcMu = (enthalpy + hp20 + cpa * (conditions.tK - T0) + cpb * (Math.pow(conditions.tK, 2) - Math.pow(T0, 2)) / 2 - cpc * (1 / conditions.tK - 1 / T0) + 2 * cpd * (Math.pow(conditions.tK, 1. / 2.) - Math.pow(T0, 1. / 2.)) - conditions.tK * (entropy + smax * q20 + 2 * cpa * (Math.log(Math.pow(conditions.tK, 1. / 2.)) - Math.log(Math.pow(T0, 1. / 2.))) + cpb * (conditions.tK - T0) - (cpc / 2) * (Math.pow(conditions.tK, -2) - Math.pow(T0, -2)) - 2 * cpd * (Math.pow(conditions.tK, -1. / 2.) - Math.pow(T0, -1. / 2.))) + addG + gland);

		cacheThermo.Mu = CalcMu;
		return CalcMu;
	}//END FUNCTION CalcMu

	public static double AQGIBBS(solutionProperty conditions /*double pkb, double tK*/, double[] properties, double w, double alphaAq) throws ExecutionException, InterruptedException {
//! adapted from EP07
//		if (cacheThermo.inCacheGibbs(new double[]{ pkb, tK, w, alphaAq})) return cacheThermo.Gibbs;
		cacheThermo.conditionsGibbs = conditions; //new double[]{pkb, tK, w, alphaAq};

		double  AQGIBBS;
		//double  pkb, tK, wCO2NaCl0, alphaAq;
		//double /*REAL(DP), DIMENSION(15) ::*/ properties;
		double /*REAL(DP), PARAMETER ::*/ af = 20.0, dkdt = -0.00015, T0 = 298.15;
		double  entropy, cpa, cpb, cpc, cpd, nmol;
		double  alpha, kappa298, landaut, smax, vmax;
		double  dalphadT, alpha0, rhoW, rhoW0, Beta, cPstar, pA, transChem;
		double  enthalpy, volume, tKprime;
		//double /*REAL(DP), PARAMETER ::*/ R = 8.314462e-3;

		enthalpy = properties[1 - 1];
		entropy = properties[2 - 1] / 1000.0;
		volume = properties[3 - 1];
		cpa = properties[4 - 1];
		cpb = properties[5 - 1] * 1e-5;
		cpc = properties[6 - 1];
		cpd = properties[7 - 1];
		alpha = properties[8 - 1] * 1e-5;
		kappa298 = properties[9 - 1];
		landaut = properties[10 - 1];
		smax = properties[11 - 1] / 1000.0;
		vmax = properties[12 - 1];
		nmol = properties[15 - 1]; //Verify index OK

		dalphadT = 9.62656;
		alpha0 = 25.728;

		rhoW = waterdensity(conditions /*new double[]{pkb, tK}*//*, rhoW*/);
		rhoW0 = 997.0470390177107; //!d water in mol/ cm^3 calculated with IAPWS at 1 bar 298.15 K

		//![Beta]= 10^6*(Log[findDensWater[0.001+0.001, tk]/rhozero]-Log[   findDensWater[0.001, tk]/rhozero])/0.001; at tk=298.15*)

		Beta = 45240.510063909365;
		cPstar = cpa - 298.15 * cpb;

		if (conditions.tK <= 500) {
			tKprime = conditions.tK;
		} else { //ELSE
			tKprime = 500.0;
		} //END IF

/*   Add Energy required to transform chemical potential on molal scale to mole fraction scale*/
/*   This is NOT the same expression as Evans and Powell for RT ln G//! Theirs is wrong*/
		pA = Math.pow((1.0 / (1.0 + 1000.0 / mH2O/*18.0153d0*/)), 1.0 / nmol);
/*   NOTA BENE = correct expression includes nmol, the number of moles b+d produced when dissolving one mol of pure phase AbCd into phase component (AbCd)pm_ 1/(b+d)*/

		transChem = -Rm * conditions.tK * Math.log(pA) - (2.0 * w) / (1.0 + alphaAq);

		AQGIBBS = enthalpy - conditions.tK * entropy + conditions.pkb * volume + cpb * (298.15 * conditions.tK
				- (Math.pow(298.15, 2)) / 2 - (Math.pow(conditions.tK, 2)) / 2)
				+ (cPstar / (298.15 * (dalphadT * 1e-6))) * ((alpha0 * 1e-5) * (conditions.tK - 298.15)
				+ (conditions.tK / tKprime) * Math.log(rhoW / rhoW0)) - (cPstar / (298.15 * (dalphadT * 1e-6)))
				* (Beta * (conditions.pkb - 0.001) * 1e-6) + transChem;

		cacheThermo.Gibbs = AQGIBBS;
		return AQGIBBS;
	}//END FUNCTION AQGIBBS
}//END MODULE THERMO
