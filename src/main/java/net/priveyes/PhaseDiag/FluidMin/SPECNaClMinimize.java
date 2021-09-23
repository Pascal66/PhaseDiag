package net.priveyes.PhaseDiag.FluidMin;

import static net.priveyes.PhaseDiag.FluidMin.Constants.Rm;
import static net.priveyes.PhaseDiag.FluidMin.solutionProperty.Wlist;
import static net.priveyes.PhaseDiag.FluidMin.solutionProperty.alphas;
import static net.priveyes.PhaseDiag.FluidMin.solutionProperty.nH2Oapp;
import static net.priveyes.PhaseDiag.THERMO.thermo.RTlnGamASF;
import static net.priveyes.PhaseDiag.FluidMin.Minimizer.BiVectorFunction;
import static net.priveyes.PhaseDiag.FluidMin.Minimizer.ITMAX;
import static net.priveyes.PhaseDiag.FluidMin.Minimizer.minimize;

import net.priveyes.PhaseDiag.FluidMin.SolubMC.cacheFunction;

import java.util.concurrent.ExecutionException;
import java.util.concurrent.RecursiveTask;

public class SPECNaClMinimize extends RecursiveTask<double[]> {
	private final double[] optionsSPECNaCl;
	private final double[] param;
	private double[] result = new double[2];

	SPECNaClMinimize(double[] param, double[] optionsSPECNaCl) {
		this.param = param;

		param[0]/*[1]*/ = 0.985; //!initialGuess
		param[1]/*[1]*/ = param[0]/*[1]*/ * 1.001;
		this.optionsSPECNaCl = optionsSPECNaCl;
	}
	@Override
	protected double[] compute() {

		try {
			result = minimize/*NELDERMEAD*/(SPECNaCl, param, 0.001/*initialDelta*/, /*, ftol,*/ /*SPECNaCl,*/ ITMAX, null/*new String[]{"SPECNaCl", "XNapmc"}*/, optionsSPECNaCl);
		} catch (ExecutionException | InterruptedException e) {
			e.printStackTrace();
		}
		return result;
	}
	/**
	 * Use the Gibbs energy of the reaction Na+ + Cl- = NaClo as a cost function to estimate Nacl association
	 */
	static BiVectorFunction/*BiFunction<double[], double[], Double>*/  SPECNaCl
			/*(double[] XNapDIM, double[] optionsSPECNaCl)*/ = (XNapDIM, optionsSPECNaCl) -> {
		cacheFunction.callSPECNaCl++;

		double   SPECNaCl;
		double   pkb, tK, molalNaCl, nCO2, DHNaClpm;
		double   mu0NaClpm, XNap, muHltLiq, muNaClpm, mu0NaCl0;
		double[] proportions;
		double   /*nH2Oapp, nH2Oo,*/ nNaClpm, nHltLiq, total, pCO2, pH2Oo;
		double   pHltLiq, pNaClpm, aidH2Oo, aidHltLiq, aidNaClpm;
		double   rtlnGamHltLiq, rtlnGamNaClpm;
		//double /*REAL(DP), PARAMETER ::*/ R = 8.314462d - 3;
		//double[] /*REAL(DP), DIMENSION(1) ::*/ tmp = new double[1 /*+ 1*/];

		XNap = XNapDIM[0];

		if (XNap < 0.000001 ||/*.OR.*/ XNap > 0.99999999999) {
			SPECNaCl = 1.0e20;
		} else { //ELSE
			pkb = optionsSPECNaCl[0];
			tK = optionsSPECNaCl[1];
			molalNaCl = optionsSPECNaCl[2];
			nCO2 = optionsSPECNaCl[3];
			DHNaClpm = optionsSPECNaCl[4];
			mu0NaClpm = optionsSPECNaCl[5];
			mu0NaCl0 = optionsSPECNaCl[6];
			//alphas = Arrays.copyOfRange(optionsSPECNaCl, 8 - 1, 11/*-1*/)/*8:11)*/; //!alphas = (/1.0, aCO2, aNaCl0, aNaClpm  /)
//C'est Wlist ...
//			VecMath.setcolumn(Ws, 1 - 1, /*(:,1) =*/ Arrays.copyOfRange(optionsSPECNaCl, 12 - 1, 15/*-1*/))/*12:15)*/; //! (/ nan,      wWCO2,      wWNaClo,      wWNaClpm/)
//			VecMath.setcolumn(Ws, 2 - 1, /*(:,2) =*/ Arrays.copyOfRange(optionsSPECNaCl, 16 - 1, 19/*-1*/))/*16:19)*/; //! (/ wWCO2,    nan,        wCO2NaCl0,    wCO2NaClpm/)
//			VecMath.setcolumn(Ws, 3 - 1, /*(:,3) =*/ Arrays.copyOfRange(optionsSPECNaCl, 20 - 1, 23/*-1*/))/*20:23)*/; //! (/ wWNaClo,  wCO2NaCl0,  nan,          wNaCloNaClpm/)
//			VecMath.setcolumn(Ws, 4 - 1, /*(:,4) =*/ Arrays.copyOfRange(optionsSPECNaCl, 24 - 1, 27/*-1*/))/*24:27)*/; //! (/ wWNaClpm, wCO2NaClpm, wNaCloNaClpm, nan/)

			//nH2Oapp = 1000.0 / mH2O /*18.0153*/;
			//nH2Oo = nH2Oapp;

			nNaClpm = 2.0 * XNap *  molalNaCl;
			nHltLiq = (1.0 - XNap) * molalNaCl;

			total = nH2Oapp + nNaClpm + nHltLiq + nCO2;

			pCO2 = nCO2 / total;
			pH2Oo = nH2Oapp /*nH2Oo*/ / total;
			pHltLiq = nHltLiq / total;
			pNaClpm = nNaClpm / total;

			aidH2Oo = pH2Oo;
			aidHltLiq = pHltLiq;
			aidNaClpm = pNaClpm;

			proportions = new double[]{pH2Oo, pCO2, pHltLiq, pNaClpm};

			rtlnGamHltLiq = RTlnGamASF(proportions, Wlist /*Ws*/, alphas, 2);

			muHltLiq = mu0NaCl0 + Rm *  tK * Math.log(aidHltLiq) + rtlnGamHltLiq; //! R tK dhHltLiq

			rtlnGamNaClpm = RTlnGamASF(proportions, Wlist /*Ws*/, alphas, 3);

			muNaClpm = mu0NaClpm + Rm *  tK * Math.log(aidNaClpm) + Rm *  tK * DHNaClpm + rtlnGamNaClpm;

			SPECNaCl = Math.pow(2.0 * muNaClpm - muHltLiq, 2);

		}//END IF
		return SPECNaCl;
	};//END FUNCTION SPECNaCl

}
