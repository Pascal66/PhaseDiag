package net.priveyes.PhaseDiag.FluidMin;

import static net.priveyes.PhaseDiag.FluidMin.Constants.Rm;
import static net.priveyes.PhaseDiag.FluidMin.SPECNaClMinimize.SPECNaCl;
import static net.priveyes.PhaseDiag.FluidMin.solutionProperty.Wlist;
import static net.priveyes.PhaseDiag.FluidMin.solutionProperty.alphas;
import static net.priveyes.PhaseDiag.FluidMin.solutionProperty.nH2Oapp;
import static net.priveyes.PhaseDiag.THERMO.thermo.RTlnGamASF;

import net.priveyes.PhaseDiag.AmoebaCSharp.Amoeba;
import net.priveyes.PhaseDiag.AmoebaCSharp.Amoeba.Solution;
import net.priveyes.PhaseDiag.FluidMin.SolubMC.cacheFunction;
import net.priveyes.PhaseDiag.FluidMin.Minimizer.BiVectorFunction;
import net.priveyes.PhaseDiag.num_rep.VecMath;

import java.util.concurrent.ExecutionException;
import java.util.concurrent.RecursiveTask;
import java.util.concurrent.atomic.AtomicReference;

public class SolubCO2BrineMinimize extends RecursiveTask<double[]> {
	private final double[] optionsSPECNaCl;
	private final double[] param;
	private double[] result = new double[2];

	SolubCO2BrineMinimize(double[] param, double[] optionsSPECNaCl) {
		this.param = param;
		this.optionsSPECNaCl = optionsSPECNaCl;
	}
	@Override
	protected double[] compute() {
	//	param[0]/*[1]*/ = 0.985; //!initialGuess
	//	param[1]/*[1]*/ = param[0]/*[1]*/ * 1.001;

		try {
			result = Minimizer.minimize/*NELDERMEAD*/(SolubCO2Brine, param, 0.001/*initialDelta*/, /*, ftol,*/ /*SPECNaCl,*/ Minimizer.ITMAX, null/*new String[]{"SPECNaCl", "XNapmc"}*/, optionsSPECNaCl);
		} catch (ExecutionException | InterruptedException e) {
			e.printStackTrace();
		}

		return result;
	}
	/**
	 * Objective function to calculate CO2 solubility in brine:
	 * output Gibbs energy of the system, to be minimized.
	 * <p>
	 * Order of the components:
	 * 1: WATER
	 * 2: CO2
	 * 3: NaClo
	 * 4: NaClpm
	 */

	private static final BiVectorFunction/*BiFunction<double[], double[], Double>*/  SolubCO2Brine
			/*(double[] nCO2, RealVector/ optionsSolubCO2)*/ = (nCO2, optionsSolubCO2) -> {
		cacheFunction.callSoluBCO2Brine++;
// IMPLICIT NONE
//!REAL(DP), DIMENSION(:) :: xCO2rich
//double[] /*REAL(DP), DIMENSION(:) ::*/ nCO2;
//double[] /*REAL(DP), DIMENSION(:) ::*/ optionsSolubCO2;
		double SolubCO2Brine;
		double pkb, tK, molalNaCl, DHNaClpm, mu0NaClpm, mu0NaCl0, solubCO2W;
		double solubWCO2, muGasPhase, xNaClapp, xCO2, /*nH2Oo,*/ xH2Oo;
		double coef;/*nH2Oapp,*/
		double ftol;
		AtomicReference<Double> XNap = new AtomicReference<>(0.);
		double nNaClpm;
		double nNaCl0, total, pCO2, pH2Oo, pNaCl0, pNaClpm, Gid, Gxs;
		double muPoint, Gsys, fracNaClapp, xCO2gas, xCO2rich;
		int /*INTEGER(I4B) ::*/ niter;
		double[] /*REAL(DP), DIMENSION(27) ::*/ optionsSPECNaCl;
		//double /*REAL(DP), PARAMETER ::*/ R = 8.314462e-3;
		//double[] /*REAL(DP), DIMENSION(4) ::*/ alphas;
		double[] proportions = new double[4];
		double[] propNaClfree;
		//RealVector proportions = MatrixUtils.createRealVector(new double[4 /*+ 1*/]), propNaClfree;
//		double[][] Ws = new double[4][4];
//		RealMatrix /*double[][]*/ /*REAL(DP), DIMENSION(4,4) ::*/ Ws = MatrixUtils.createRealMatrix(4 /*+ 1*/, 4 /*+ 1*/);
		double[]/*[]*/ /*REAL(DP), DIMENSION(2,1) ::*/ param = new double[2 /*+ 1*/]/*[1 + 1]*/;
		double[] /*REAL(DP), DIMENSION(2) ::*/ y = new double[2 /*+ 1*/];
		double[] /*REAL(DP), DIMENSION(3) ::*/ brineEM, gasEM;
		double[] /*RealVector*/ compo;

//! Input check
		if (nCO2[0] <= 0.0 ||/*.OR.*/ nCO2[0] > 50.0) {
			SolubCO2Brine = 1.0e100;
//! ABORT
		} else { //ELSE
			pkb = optionsSolubCO2[0];
			tK = optionsSolubCO2[1];
			molalNaCl = optionsSolubCO2[2];

			DHNaClpm = optionsSolubCO2[4];
			mu0NaClpm = optionsSolubCO2[5];
			mu0NaCl0 = optionsSolubCO2[6];
			//alphas = Arrays.copyOfRange(optionsSolubCO2, 8 - 1, 11/*-1*/)/*8:11)*/; //!alphas = (/1.0, aCO2, aNaCl0, aNaClpm  /)
// C'est Wlist ...
//			VecMath.setcolumn(Ws, 1 - 1, /*);:,1) =*/ Arrays.copyOfRange(optionsSolubCO2, 12 - 1, 15/*-1*/)/*12:15*/); //! (/ nan,      wWCO2,      wWNaClo,      wWNaClpm/)
//			VecMath.setcolumn(Ws, 2 - 1, /*);:,2) =*/ Arrays.copyOfRange(optionsSolubCO2, 16 - 1, 19/*-1*/)/*16:19*/); //! (/ wWCO2,    nan,        wCO2NaCl0,    wCO2NaClpm/)
//			VecMath.setcolumn(Ws, 3 - 1, /*);:,3) =*/ Arrays.copyOfRange(optionsSolubCO2, 20 - 1, 23/*-1*/)/*20:23*/); //! (/ wWNaClo,  wCO2NaCl0,  nan,          wNaCloNaClpm/)
//			VecMath.setcolumn(Ws, 4 - 1, /*);:,4) =*/ Arrays.copyOfRange(optionsSolubCO2, 24 - 1, 27/*-1*/)/*24:27*/); //! (/ wWNaClpm, wCO2NaClpm, wNaCloNaClpm, nan/)

			/*! Need to know a) solubility of CO2 in water, b) solub of water in CO2 and c) associated chemical potentials*/
			solubCO2W = optionsSolubCO2[27]; //! CO2 in water-rich phase (mole fraction of CO2)
			solubWCO2 = optionsSolubCO2[28]; //! CO2 in CO2-rich phase (mole fraction of CO2, water = 1-CO2//!)

			propNaClfree = /*MatrixUtils.createRealVector*/new double[]{solubCO2W /*1.0 - solubWCO2*/, solubWCO2, 0.0, 0.0};
			/*!       chem. pot. of NaCl-free water-saturated CO2 phase
			  Note that absolute chem potentials are unneeded
			  Because they cancel analytically in the second derivative of the solvus calculation*/

			muGasPhase = Rm * tK * (solubWCO2 * Math.log(solubWCO2) + (solubCO2W/*1.0 - solubWCO2*/) * Math.log(solubCO2W/*1.0 - solubWCO2*/)) + propNaClfree[0] * RTlnGamASF(propNaClfree, Wlist /*Ws*/, alphas, 0) + propNaClfree[1] * RTlnGamASF(propNaClfree, Wlist /*Ws*/, alphas, 1);

			//System.out.println("SolubCO2Brine 763 muGasPhase " + muGasPhase);

			//nH2Oo = 1000.0 / mH2O/*18.0153d0*/;
			total = nH2Oapp /*nH2Oo*/ + molalNaCl + nCO2[0];
			compo = VecMath.vxs(new double[]{nH2Oapp /*nH2Oo*/, molalNaCl, nCO2[0]}, 1 / total);
			//compo = MatrixUtils.createRealVector(new double[]{nH2Oo, molalNaCl, nCO2[1 - 1]}).mapDivideToSelf(total);// / total;
			/* for Gsys calculation: fraction of gas end-member versus CO2-free brine end-member*/
			xCO2rich = compo[2] / solubWCO2;

			xH2Oo = compo[0];
			xNaClapp = compo[1]; //!
			xCO2 = compo[2];
			//nH2Oapp = nH2Oo;

			/* Calculate NaCl speciation (nCO2-dependent)*/

			optionsSPECNaCl = VecMath.copyvec(optionsSolubCO2);
			//optionsSPECNaCl = Arrays.copyOfRange(optionsSolubCO2, 1 - 1, 27/*-1*/)/*[1:27]*/;
			optionsSPECNaCl[2] = molalNaCl;
			optionsSPECNaCl[3] = nCO2[0];

			param[0]/*[1]*/ = 0.985; //!initialGuess
			param[1]/*[1]*/ = param[0]/*[1]*/ * 1.001;

//			y[1 - 1] = SPECNaCl.apply((new double[]{param[0], 0}/*[1]/*,:]*/), optionsSPECNaCl);
//			y[2 - 1] = SPECNaCl.apply((new double[]{param[1], 0}/*[2]/*,:]*/), optionsSPECNaCl);
//			double initialDelta = y[1] - y[0]; //PP
			//System.out.println("SolubCO2Brine 789 SPECNaCl initialGuessMin " + VecMath.toString(y));
			//System.out.println("Solub2Brine Fortran is:    17.076941356500200        15.671982657304325      ");
			ftol = 1.0e-6;
			niter = Minimizer.ITMAX;
//			SPECNaClMinimize xt = new SPECNaClMinimize(param, optionsSPECNaCl);
//			XNap.set(xt.invoke()[0]);

//TODO BUG Minimizer 4x plus lent meme en thread
//			System.out.println("SolubMC 728 SPECNaCl Solution: x0 = " + VecMath.toString(x));
//			System.out.println("Solub2Brine Fortran final is:   0.99709340120552148       0.99709340120554968      ");

			Amoeba amo = new Amoeba(2, 1, param[0], param[1], 100, SPECNaCl, optionsSPECNaCl /*new double[]{1013.25, 298.15}*/);
			Solution sol = amo.Solve();
//			Log.e("Amoeaba", VecMath.toString(sol.vector));
			XNap.set(sol.vector[0]);

//			try {
//          XNap = x.get()[0]; //(x/*minimizer.x[0]*/[1]/*[1]*/ + x/*minimizer.x[0]*/[0]/*[1]*/) / 2.0;
//			} catch (ExecutionException e) {
//				e.printStackTrace();
//			} catch (InterruptedException e) {
//				e.printStackTrace();
//			}

			nNaClpm = 2.0 * XNap.get() * molalNaCl;
			nNaCl0 = (1.0 - XNap.get()) * molalNaCl;

			total = nH2Oapp + nNaClpm + nNaCl0 + nCO2[0];

//			System.out.println("SolubCO2Brine 824 XNap " + XNap + " nNaClpm " + nNaClpm + " nNaCl0 " + nNaCl0 + " total " + total);

			pCO2 = nCO2[0] / total;
			pH2Oo = nH2Oapp /*nH2Oo*/ / total;
			pNaCl0 = nNaCl0 / total;
			pNaClpm = nNaClpm / total;

// calculate mu point
// Gmix = Gid + Gnonid

			proportions[0] = pH2Oo;
			proportions[1] = pCO2;
			proportions[2] = pNaCl0;
			proportions[3] = pNaClpm;
//			System.out.println("SolubCO2Brine 838 proportions pCO2 " + pCO2 + " pH2Oo " + pH2Oo + " pNaCl0 " + pNaCl0 + " pNaClpm " + pNaClpm);
			//Le rapport adimensionnel am = Pm / P0 est appelé son activité.
			//µm = µ0 + RT⋅ln am
			//G (du mélange) = sum( ni µi ) ni est le nombre de moles et µ0 le potentiel chimique standard du constituant.
			Gid = pH2Oo * Rm * tK * Math.log(pH2Oo) + pCO2 * Rm * tK * Math.log(pCO2) + pNaCl0 * Rm * tK * Math.log(pNaCl0) + pNaClpm * Rm * tK * Math.log(pNaClpm);
			/* Gnonid = use RTlnGamASF(prop,Ws,alphas,l) = Sum xi . RTlnGamASFi*/
			Gxs = proportions[0] * RTlnGamASF(proportions, Wlist /*Ws*/, alphas, 0) + proportions[1] * RTlnGamASF(proportions, Wlist /*Ws*/, alphas, 1) + proportions[2] * RTlnGamASF(proportions, Wlist /*Ws*/, alphas, 2) + proportions[3] * RTlnGamASF(proportions, Wlist /*Ws*/, alphas, 3);
			/* ensure that muPoint = muGasPhase for xCO2 = 1 (this should be okay if muGasPhase is calculated with Gid+Gnonid from RTlnGamASF)*/
			muPoint = Gid + Gxs;

			//!write(*,*) ' mu point = ', muPoint
//			System.out.println("Solub2Brine fortran muPoint is:  -0.68991862542504956      ");
//			System.out.println("Solub2Brine MuPoint " + muPoint + " Gid " + Gid + " Gxs " + Gxs);

			/* calculate mu system at fixed composition (greater than CO2 solub...) with lever rule
 We shall use the solubily of CO2 in (pure) water as a boundary so fixed composition at 1.1*solub CO2 or 50%
 Note that as long as negative proportions are allowed, the choice of system composition will have no effect on the result of the minimisation
 We will choose 50% NaCl-free water-saturated CO2 and 50% CO2-free water saturated */

			Gsys = muGasPhase * ((0.5 - xCO2rich) / (1.0 - xCO2rich)) + muPoint * ((1.0 - 0.5) / (1.0 - xCO2rich));

			SolubCO2Brine = Gsys;

		}//END IF //! input check
		return SolubCO2Brine;
	};//END FUNCTION SolubCO2Brine
}
