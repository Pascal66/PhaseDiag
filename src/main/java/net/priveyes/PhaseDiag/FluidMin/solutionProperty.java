package net.priveyes.PhaseDiag.FluidMin;

import static net.priveyes.PhaseDiag.FluidMin.Constants.mH2O;

public class solutionProperty {
public  double tK = 0.;
		double tC = 0.;
public  double pkb = 0.;
public  double molalNaCl = 0.;
public static final double nH2Oapp = 1000.0 / mH2O/*18.0153d0*/; //TODO Move to constant
public static double[] /*REAL(DP), DIMENSION(4)::*/ alphas;
public static double[][] /*REAL(DP), DIMENSION(4,4)::*/ Wlist = new double[4][4]; //MatrixUtils.createRealMatrix(4 /*+ 1*/, 4 /*+ 1*/);
	public solutionProperty(double pkb, double tK, double NaCl) {
		this.pkb = pkb; this.tK = tK; this.molalNaCl = NaCl;
		this.tC = tK - 273.15;
	}
	public boolean isValid() {
		//tC = vTemperature/*conditions*/.getEntry(i/*, 2-1*/)/*[i][2]*/;
		//pb = vPressure/*conditions*/.getEntry(i/*, 1-1*/)/*[i][1]*/;

		final double kk1 = 195.27910870306826;
		final double kk2 = 90.36072277568843;
		final double kk3 = -8.944991890967263;
		final double kk4 = -107.92860602097733;

			/*First test is critical mixing curve above critical point of water
			  second test for boiling curve*/
		boolean valid = (((this.pkb > 0.22) &&/*.AND.*/ (tC < (kk1 + kk2 * this.pkb + kk3 * (Math.pow(this.pkb, 2)) + kk4 * Math.log(this.pkb) - 5.0)))
				||/*.OR.*/ ((this.pkb < 0.22) &&/*.AND.*/ (this.pkb > 0.999
				* ((1e-8) * Math.exp(-2836.5744 / (Math.pow(this.tK, 2)) - 6028.076559 / this.tK + 19.54263612 - 0.02737830188 * this.tK
				+ 1.6261698e-5 * (Math.pow(this.tK, 2)) + 7.0229056e-10 * (Math.pow(this.tK, 3)) - 1.8680009e-13 * (Math.pow(this.tK, 4))
				+ 2.7150305 * Math.log(this.tK)))))); //! single precision - not important

		return valid;
	}
}
