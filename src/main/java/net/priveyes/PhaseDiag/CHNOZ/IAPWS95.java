package net.priveyes.PhaseDiag.CHNOZ;

import static java.lang.Double.NaN;
import static java.lang.Math.exp;
import static java.lang.Math.log;
import static java.lang.Math.pow;
import static java.lang.Math.sqrt;

import android.util.Log;

import net.priveyes.PhaseDiag.CHNOZ.IAPWScommon.Rho;
import net.priveyes.PhaseDiag.CHNOZ.IAPWScommon.Temperature;
import net.priveyes.PhaseDiag.CHNOZ.IAPWScommon.idealgas;
import net.priveyes.PhaseDiag.CHNOZ.IAPWScommon.idealgas.Phi;
import net.priveyes.PhaseDiag.CHNOZ.IAPWScommon.idealgas.d2Delta;
import net.priveyes.PhaseDiag.CHNOZ.IAPWScommon.idealgas.d2Psi;
import net.priveyes.PhaseDiag.CHNOZ.IAPWScommon.idealgas.dDelta;
import net.priveyes.PhaseDiag.CHNOZ.IAPWScommon.idealgas.dPsi;
import net.priveyes.PhaseDiag.CHNOZ.IAPWScommon.residual;
import net.priveyes.PhaseDiag.FluidMin.Constants;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.function.Supplier;
import java.util.stream.IntStream;

/**
 * @use "https://r-forge.r-project.org/scm/viewvc.php/pkg/CHNOSZ/R/IAPWS95.R?view=markup&revision=176&root=chnosz"
 * functions for properties of water using the IAPWS-95 formulation (Wagner and Pruss, 2002)
 */
public class IAPWS95 {
	public static double delta;
	static double tau;

	/**Test @see http://twt.mpei.ac.ru/mcs/worksheets/iapws/IAPWS95.xmcd */
	public static void initIAPWS95(final double[] tK/*=298.15*/, final double[] Rh/*=1000*/) {
		Rho.defaut = Rh[0]; Temperature.defaut = tK[0];
		delta = Rho.defaut / Rho.critical;
		tau = Temperature.critical / Temperature.defaut;
		IAPWS95_idealgas();
/*		Log.i("Initial", "tau "+tau+" delta "+ delta);
		Log.i("Ideal", "Phi "+idealgas.phi+
				"Phi delta " + idealgas.Phi.delta+
				"Phi delta delta" + idealgas.Phi.Delta.delta+
				"Phi tau " + idealgas.Phi.tau+
				"Phi tau tau " + idealgas.Phi.Tau.tau+
				"Phi delta tau " + idealgas.Phi.Delta.tau);
*/
		IAPWS95_residual();
		/*Log.i("Residual", "Phi "+residual.phi+
				"Phi delta " + residual.Phi.delta+
				"Phi delta delta" + residual.Phi.Delta.delta+
				"Phi tau " + residual.Phi.tau+
				"Phi tau tau " + residual.Phi.Tau.tau+
				"Phi delta tau " + residual.Phi.Delta.tau);
*/
		//IAPWS95get(new String[]{"P","Z","U","S","H","G", "CV","CP","W","MU"}, tK, Rh);
		//return (p.get() /
	}
	/**
	 * the IAPWS-95 formulation for ordinary water substance Wagner and Pruss, 2002
	 */
	public static double[] IAPWS95get /*<- function*/(String[]/*Supplier<Double>[]*/ property, final double[] tK/*=298.15*/, final double[] Rh/*=1000*/) {
		Temperature.T = tK;
		Rho.rho = Rh;
		//property = tolower(property)
		// specific and molar gas constants
		double R = 0.46151805; // kJ kg-1 K-1
		// R.M = 8.314472 // J mol-1 K-1
		// molar mass
		double M = Constants.mH2O; //18.015286; // g mol-1

		Map<String, Supplier<Double>> giveThisProperty = new HashMap<>();
		// relation of thermodynamic properties to Helmholtz free energy
		Supplier<Double> a = () -> {
			double x = idealgas.phi + residual.phi;
			return (x * R * Temperature.defaut);
		};
		giveThisProperty.putIfAbsent("A", a);
		// Table 6.3
		Supplier<Double> p = () -> {
			double x = 1 + delta * residual.Phi.delta;
			return (x * Rho.defaut * R * Temperature.defaut / 1000);  // for MPa
		};
		giveThisProperty.putIfAbsent("P", p); //Pressure
		Supplier<Double> s = () -> {
			double x = tau * (Phi.tau + residual.Phi.tau) - idealgas.phi - residual.phi;
			return (x * R);
		};
		Supplier<Double> u = () -> {
			double x = tau * (Phi.tau + residual.Phi.tau);
			return (x * R * Temperature.defaut);
		};
		giveThisProperty.putIfAbsent("U", u); //Internal Energy
		Supplier<Double> h = () -> {
			double x = 1 + tau * (Phi.tau + residual.Phi.tau) + delta * residual.Phi.delta;
			return (x * R * Temperature.defaut);
		};
		giveThisProperty.putIfAbsent("H", h); // Specific Enthalpy
		Supplier<Double> g = () -> {
			double x = 1 + idealgas.phi + residual.phi + delta * residual.Phi.delta;
			return (x * R * Temperature.defaut);
		};
		giveThisProperty.putIfAbsent("G", g);
		Supplier<Double> cv = () -> {
			double x = -pow(tau, 2) * (Phi.Tau.tau + residual.Phi.Tau.tau);
			return (x * R);
		};
		giveThisProperty.putIfAbsent("CV", cv); //Isochoric heat capacity
		Supplier<Double> cp = () -> {
			double x = -pow(tau, 2) * (Phi.Tau.tau + residual.Phi.Tau.tau) + pow((1 + delta * residual.Phi.delta - delta * tau * residual.Phi.Delta.tau), 2) / (1 + 2 * delta * residual.Phi.delta + pow(delta, 2) * residual.Phi.Delta.delta);
			return (x * R);
		};
		giveThisProperty.putIfAbsent("CP", cp); //Isobaric heat capacity
// 20090420 speed of sound calculation is incomplete
// (delta.liquid and drhos.dT not visible)
//  cs = function() {
//    x = -tau^2*(idealgas('phi.tau.tau')+residual('phi.tau.tau')) +
//         (1+delta*residual('phi.delta')-delta*tau*residual('phi.delta.tau'))^2 /
//         (1+2*delta*residual('phi.delta')+delta^2*residual('phi.delta.delta')) *
//         ((1+delta.liquid*residual('phi.delta')-delta.liquid*tau*residual('phi.tau.tau'))-rho.critical/(R*delta.liquid)*drhos.dT)
//    return(x*R)
//  }
		Supplier<Double> w = () -> {
			double x = 1 + 2 * delta * residual.Phi.delta + pow(delta, 2) * residual.Phi.Delta.delta - pow((1 + delta * residual.Phi.delta - delta * tau * residual.Phi.Delta.tau), 2) / pow(tau, 2) * (Phi.Tau.tau + residual.Phi.Tau.tau);
			return (sqrt(x * R * Temperature.defaut));
		};
		giveThisProperty.putIfAbsent("W", w);
		Supplier<Double> mu = () -> {
			double x = -(delta * residual.Phi.delta + pow(delta, 2) * residual.Phi.Delta.delta + delta * tau * residual.Phi.Delta.tau) / ((1 + delta * residual.Phi.delta - delta * tau * pow(residual.Phi.Delta.tau, 2)) - pow(tau, 2) * (Phi.Tau.tau + residual.Phi.Tau.tau) * (1 + 2 * delta * residual.Phi.delta + pow(delta, 2) * residual.Phi.Delta.delta));
			return (x / (R * Rho.defaut));
		};
		giveThisProperty.putIfAbsent("MU", mu); //Joule Thomson coefficient

		//// run the calculations
		//Double ww = null;
		//double my.T = T;
		//double my.rho = rho;
		double[] t = new double[Temperature.T.length];
		for (String aProperty : property) {
			//t = numeric()
			for (int i = 0; i < Temperature.T.length; i++) /*in 1:length(my.T))*/ {
				//Rho.defaut = /*myT =*/Rho.rho[i] ;
				//Temperature.defaut = /*myRho*/Temperature.T[i] ;
				// Equation 6.4
				//delta = Rho.defaut / Rho.critical;
				//tau = Temperature.critical / Temperature.defaut;
				initIAPWS95(new double[]{Temperature.T[i]}, new double[]{Rho.rho[i]});
				if (giveThisProperty.containsKey(aProperty)) {
					t[i] = giveThisProperty.get(aProperty).get();
				} else t[i] = NaN;
			}
			//t = data.frame(t);
			//if (j == 1) {
			//	ww = t
			//} else {
			//	ww = cbind(ww, t);
			//}
			Log.w("Property "+aProperty, "Is "+ Arrays.toString(t) +" for Tk " + Temperature.defaut+ " @Rho "+Rho.defaut);
		} //colnames(ww) = property;
		//return (ww);
		return t;
	}

////// unexported functions //////

	/**
	 * IAPWS95.idealgas and IAPWS95.residual are supporting functions to IAPWS95 for calculating
	 * the ideal-gas and residual parts in the IAPWS-95 formulation.
	 * The value of p can be one of phi, phi.delta, phi.delta.delta, phi.tau, phi.tau.tau, or phi.delta.tau,
	 * to calculate the specific dimensionless Helmholtz free energy (phi) or one of its derivatives.
	 */

	private static void /*double*/ IAPWS95_idealgas(/*Supplier<Double> p, Double delta, double tau*/) {
		//// the ideal gas part in the IAPWS-95 formulation
		// from Table 6.1 of Wagner and Pruss, 2002
		double[] n = new double[]{-8.32044648201, 6.6832105268, 3.00632, 0.012436, 0.97315, 1.27950, 0.96956, 0.24873};
		double[] gamma = new double[]{NaN, NaN, NaN, 1.28728967, 3.53734222, 7.74073708, 9.24437796, 27.5075105};
		// Equation 6.5

		idealgas.phi = log(delta) + n[1-1] + n[2-1] * tau + n[3-1] * log(tau) + IntStream.range(4 - 1, 8).mapToDouble(i -> n[i] * log(1 - exp(-gamma[i] * tau))).sum();
		// derivatives from Table 6.4
		Phi.delta = 1 / delta + 0 + 0 + 0 + 0;
		Phi.Delta.delta = -1 / pow(delta, 2. + 0. + 0. + 0. + 0.);
		Phi.tau = 0 + 0 + n[2-1] + n[3-1] / tau + IntStream.range(4 - 1, 8).mapToDouble(i -> n[i] * gamma[i] * (pow((1 - exp(-gamma[i] * tau)), -1) - 1)).sum();
		Phi.Tau.tau = 0 + 0 + 0 - n[3-1] / pow(tau, 2) - IntStream.range(4 - 1, 8).mapToDouble(i -> n[i] * pow(gamma[i], 2) * exp(-gamma[i] * tau) * pow((1 - exp(-gamma[i] * tau)), -2)).sum();
		Phi.Delta.tau = 0. + 0. + 0. + 0. + 0.;
		//return (p.get() /*get(p)()*/);
	}

	/**
	 * the residual part in the IAPWS-95 formulation
	 */
	private static void /*double*/ IAPWS95_residual(/*Supplier<Double> p, double delta, double tau*/) {
		// from Table 6.2 of Wagner and Pruss, 2002

		 List<Double> ctemp = new ArrayList<>();
		rep(ctemp, NaN, 7);
		rep(ctemp, 1, 15);
		rep(ctemp, 2, 20);
		rep(ctemp, 3, 4);
		rep(ctemp, 4, 1);
		rep(ctemp, 6, 4);
		rep(ctemp, NaN, 5);
		final Double[] c = ctemp.toArray(new Double[0]);
		//System.out.println("c len "+c.length+" "+ Arrays.toString(c));
		double[] d = new double[]/*c*/{1, 1, 1, 2, 2, 3, 4, 1, 1, 1, 2, 2, 3, 4, 4, 5, 7, 9, 10, 11, 13, 15, 1, 2, 2, 2, 3, 4, 4, 4, 5, 6, 6, 7, 9, 9, 9, 9, 9, 10, 10, 12, 3, 4, 4, 5, 14, 3, 6, 6, 6, 3, 3, 3, NaN, NaN};
		double[] t = /*c*/{-0.5, 0.875, 1, 0.5, 0.75, 0.375, 1, 4, 6, 12, 1, 5, 4, 2, 13, 9, 3, 4, 11, 4, 13, 1, 7, 1, 9, 10, 10, 3, 7, 10, 10, 6, 10, 10, 1, 2, 3, 4, 8, 6, 9, 8, 16, 22, 23, 23, 10, 50, 44, 46, 50, 0, 1, 4, NaN, NaN};
		double[] n = new double[]/*c*/{0.12533547935523E-1, 0.78957634722828E1, -0.87803203303561E1, 0.31802509345418, -0.26145533859358, -0.78199751687981E-2, 0.88089493102134E-2, -0.66856572307965, 0.20433810950965, -0.66212605039687E-4, -0.19232721156002, -0.25709043003438, 0.16074868486251, -0.40092828925807E-1, 0.39343422603254E-6, -0.75941377088144E-5, 0.56250979351888E-3, -0.15608652257135E-4, 0.11537996422951E-8, 0.36582165144204E-6, -0.13251180074668E-11, -0.62639586912454E-9, -0.10793600908932, 0.17611491008752E-1, 0.22132295167546, -0.40247669763528, 0.58083399985759, 0.49969146990806E-2, -0.31358700712549E-1, -0.74315929710341, 0.47807329915480, 0.20527940895948E-1, -0.13636435110343, 0.14180634400617E-1, 0.83326504880713E-2, -0.29052336009585E-1, 0.38615085574206E-1, -0.20393486513704E-1, -0.16554050063734E-2, 0.19955571979541E-2, 0.15870308324157E-3, -0.16388568342530E-4, 0.43613615723811E-1, 0.34994005463765E-1, -0.76788197844621E-1, 0.22446277332006E-1, -0.62689710414685E-4, -0.55711118565645E-9, -0.19905718354408, 0.31777497330738, -0.11841182425981, -0.31306260323435E2, 0.31546140237781E2, -0.25213154341695E4, -0.14874640856724, 0.31806110878444};
		ctemp.clear(); rep(ctemp, NaN, 51); ctemp.addAll(Arrays.asList(20., 20., 20., NaN, NaN));
		Double[] alpha = ctemp.toArray(new Double[0]);
		ctemp.clear(); rep(ctemp, NaN, 51); ctemp.addAll(Arrays.asList(150., 150., 250., 0.3, 0.3));
		Double[] beta = ctemp.toArray(new Double[0]);
		ctemp.clear(); rep(ctemp, NaN, 51); ctemp.addAll(Arrays.asList(1.21, 1.21, 1.25, NaN, NaN));
		Double[] gamma = ctemp.toArray(new Double[0]);
		ctemp.clear(); rep(ctemp, NaN, 51); ctemp.addAll(Arrays.asList(1., 1., 1., NaN, NaN));
		Double[] epsilon = ctemp.toArray(new Double[0]);
		ctemp.clear(); rep(ctemp, NaN, 54); ctemp.addAll(Arrays.asList(3.5, 3.5));
		Double[] a = ctemp.toArray(new Double[0]);
		ctemp.clear(); rep(ctemp, NaN, 54); ctemp.addAll(Arrays.asList(0.85, 0.95));
				Double[] b = ctemp.toArray(new Double[0]);
		ctemp.clear(); rep(ctemp, NaN, 54); ctemp.addAll(Arrays.asList(0.2, 0.2));
				Double[] B = ctemp.toArray(new Double[0]);
		ctemp.clear(); rep(ctemp, NaN, 54); ctemp.addAll(Arrays.asList(28., 32.));
				Double[] C = ctemp.toArray(new Double[0]);
		ctemp.clear(); rep(ctemp, NaN, 54); ctemp.addAll(Arrays.asList(700., 800.));
				Double[] D = ctemp.toArray(new Double[0]);
		ctemp.clear(); rep(ctemp, NaN, 54); ctemp.addAll(Arrays.asList(0.32, 0.32));
				Double[] A = ctemp.toArray(new Double[0]);
		// from Table 6.5
		// Supplier Necessaire Caused by: java.lang.IllegalStateException: stream has already been operated upon or closed
		Supplier<IntStream> i1i = () -> IntStream.range(1 - 1, 7);  // 1:7;
		Supplier<IntStream> i2i = () -> IntStream.range(8 - 1, 51);  //8:51;
		Supplier<IntStream> i3i = () -> IntStream.range(52 - 1, 54); //52:54;
		Supplier<IntStream> i4i = () -> IntStream.range(55 - 1, 56); //55:56;
		// derivatives of distance function
		double[] Theta = new double[56];
		Arrays.setAll(Theta, i -> (1 - tau) + A[i] * pow((pow((delta - 1), 2)), (1 / (2 * beta[i]))));
		double[] Delta = new double[56];
		Arrays.setAll(Delta, i -> pow(Theta[i], 2) + B[i] * pow((pow((delta - 1), 2)), a[i]));
		double[] Psi = new double[56];
		Arrays.setAll(Psi, i -> exp(-C[i] * pow((delta - 1), 2) - D[i] * pow((tau - 1), 2)));

		/*dDelta.Bi.ddelta =*/ /*function(i) */
		Arrays.setAll(dDelta.Bi.ddelta, i -> b[i] * pow(Delta[i], (b[i] - 1)) * dDelta.ddelta[i]);
		/*d2Delta.Bi.ddelta2 =*/ /*function(i) */
		Arrays.setAll(d2Delta.Bi.ddelta2, i -> b[i] * (pow(Delta[i], (b[i] - 1)) * d2Delta.ddelta2[i] + (b[i] - 1) * pow(Delta[i], (b[i] - 2)) * pow(dDelta.ddelta[i], 2)));
		/*dDelta.Bi.dtau =*/ /*function(i) */
		Arrays.setAll(dDelta.Bi.dtau, i -> -2 * Theta[i] * b[i] * pow(Delta[i], (b[i] - 1)));
		/*d2Delta.Bi.dtau2 =*/ /*function(i) */
		Arrays.setAll(d2Delta.Bi.dtau2, i -> 2 * b[i] * pow(Delta[i], (b[i] - 1)) + 4 * pow(Theta[i], 2) * b[i] * (b[i] - 1) * pow(Delta[i], (b[i] - 2)));
		/*d2Delta.Bi.Ddelta.dtau =*/ /*function(i) */
		Arrays.setAll(d2Delta.Bi.Ddelta.dtau, i -> -A[i] * b[i] * 2 / beta[i] * pow(Delta[i], (b[i] - 1)) * (delta - 1) * pow((pow((delta - 1), 2)), (1 / (2 * beta[i]) - 1)) - 2 * Theta[i] * b[i] * (b[i] - 1) * pow(Delta[i], (b[i] - 2)) * dDelta.ddelta[i]);
		/*dDelta.ddelta =*/ /*function(i) */
		Arrays.setAll(dDelta.ddelta, i -> (delta - 1) * (A[i] * Theta[i] * 2 / beta[i] * pow((pow((delta - 1), 2)), (1 / (2 * beta[i]) - 1)) + 2 * B[i] * a[i] * pow(pow((delta - 1), 2), (a[i] - 1))));
		/*d2Delta.ddelta2 =*/ /*function(i) */
		Arrays.setAll(d2Delta.ddelta2, i -> 1 / (delta - 1) * dDelta.ddelta[i] + pow((delta - 1), 2) * (4 * B[i] * a[i] * (a[i] - 1) * pow((pow((delta - 1), 2)), (a[i] - 2)) + 2 * pow(A[i], 2) * pow((1 / beta[i]), 2) * pow((pow((pow((delta - 1), 2)), (1 / (2 * B[i]) - 1))), 2) + A[i] * Theta[i] * 4 / beta[i] * (1 / (2 * B[i]) - 1) * pow((pow((delta - 1), 2)), (1 / (2 * beta[i]) - 2))));
		// derivatives of exponential function
		/*dPsi.ddelta =*/ /*function(i) */
		Arrays.setAll(dPsi.ddelta, i -> -2 * C[i] * (delta - 1) * Psi[i]);
		/*d2Psi.ddelta2 =*/ /*function(i) */
		Arrays.setAll(d2Psi.ddelta2, i -> (2 * C[i] * pow((delta - 1), 2) - 1) * 2 * C[i] * Psi[i]);
		/*dPsi.dtau =*/ /*function(i) */
		Arrays.setAll(dPsi.dtau, i -> -2 * D[i] * (tau - 1) * Psi[i]);
		/*d2Psi.dtau2 =*/ /*function(i) */
		Arrays.setAll(d2Psi.dtau2, i -> (2 * D[i] * pow((tau - 1), 2) - 1) * 2 * D[i] * Psi[i]);
		/*d2Psi.ddelta.dtau =*/ /*function(i) */
		Arrays.setAll(d2Psi.Ddelta.dtau, i -> 4 * C[i] * D[i] * (delta - 1) * (tau - 1) * Psi[i]);
		// dimensionless Helmholtz free energy and derivatives

		residual.phi = i1i.get().mapToDouble(i1 -> n[i1] * pow(delta, d[i1]) * pow(tau, t[i1])).sum()
				+ i2i.get().mapToDouble(i2 -> n[i2] * pow(delta, d[i2]) * pow(tau, t[i2]) * exp(-pow(delta, c[i2]))).sum()
				+ i3i.get().mapToDouble(i3 -> n[i3] * pow(delta, d[i3]) * pow(tau, t[i3]) * exp(-alpha[i3] * pow(delta - epsilon[i3], 2) - beta[i3] * pow(tau - gamma[i3], 2))).sum()
				+ i4i.get().mapToDouble(i4 -> n[i4] * pow(Delta[i4], b[i4]) * delta * Psi[i4]).sum();
		residual.Phi.delta = i1i.get().mapToDouble(i1 -> n[i1] * d[i1] * pow(delta, (d[i1] - 1)) * pow(tau, t[i1])).sum()
				+ i2i.get().mapToDouble(i2 -> n[i2] * exp(-pow(delta, c[i2])) * (pow(delta, (d[i2] - 1)) * pow(tau, t[i2]) * (d[i2] - c[i2] * pow(delta, c[i2])))).sum()
				+ i3i.get().mapToDouble(i3 -> n[i3] * pow(delta, d[i3]) * pow(tau, t[i3]) * exp(-alpha[i3] * pow((delta - epsilon[i3]), 2) - beta[i3] * pow((tau - gamma[i3]), 2)) * (d[i3] / delta - 2 * alpha[i3] * (delta - epsilon[i3]))).sum()
				+ i4i.get().mapToDouble(i4 -> n[i4] * (pow(Delta[i4], b[i4]) * (Psi[i4] + delta * dPsi.ddelta[i4]) + dDelta.Bi.ddelta[i4] * delta * Psi[i4])).sum();
		residual.Phi.Delta.delta = i1i.get().mapToDouble(i1 -> n[i1] * d[i1] * (d[i1] - 1) * pow(delta, (d[i1] - 2)) * pow(tau, t[i1])).sum()
				+ i2i.get().mapToDouble(i2 -> n[i2] * exp(-pow(delta, c[i2])) * (pow(delta, (d[i2] - 2)) * pow(tau, t[i2]) * ((d[i2] - c[i2] * pow(delta, c[i2])) * (d[i2] - 1 - c[i2] * pow(delta, c[i2])) - pow(c[i2], 2) * pow(delta, c[i2])))).sum()
				+ i3i.get().mapToDouble(i3 -> n[i3] * pow(tau, t[i3]) * exp(-alpha[i3] * pow((delta - epsilon[i3]), 2) - beta[i3] * pow((tau - gamma[i3]), 2)) * (-2 * alpha[i3] * pow(delta, d[i3]) + 4 * pow(alpha[i3], 2) * pow(delta, d[i3]) * pow((delta - epsilon[i3]), 2) - 4 * d[i3] * alpha[i3] * pow(delta, (d[i3] - 1)) * (delta - epsilon[i3]) + d[i3] * (d[i3] - 1) * pow(delta, (d[i3] - 2)))).sum()
				+ i4i.get().mapToDouble(i4 -> n[i4] * (pow(Delta[i4], b[i4]) * (2 * dPsi.ddelta[i4] + delta * d2Psi.ddelta2[i4]) + 2 * dDelta.Bi.ddelta[i4] * (Psi[i4] + delta * dPsi.ddelta[i4]) + d2Delta.Bi.ddelta2[i4] * delta * Psi[i4])).sum();
		residual.Phi.tau = i1i.get().mapToDouble(i1 -> n[i1] * t[i1] * pow(delta, d[i1]) * pow(tau, (t[i1] - 1))).sum()
				+ i2i.get().mapToDouble(i2 -> n[i2] * t[i2] * pow(delta, d[i2]) * pow(tau, (t[i2] - 1)) * exp(-pow(delta, c[i2]))).sum()
				+ i3i.get().mapToDouble(i3 -> n[i3] * pow(delta, d[i3]) * pow(tau, t[i3]) * exp(-alpha[i3] * pow((delta - epsilon[i3]), 2) - beta[i3] * pow((tau - gamma[i3]), 2)) * (t[i3] / tau - 2 * beta[i3] * (tau - gamma[i3]))).sum()
				+ i4i.get().mapToDouble(i4 -> n[i4] * delta * (dDelta.Bi.dtau[i4] * Psi[i4] + pow(Delta[i4], b[i4]) * dPsi.dtau[i4])).sum();
		residual.Phi.Tau.tau = i1i.get().mapToDouble(i1 -> n[i1] * t[i1] * (t[i1] - 1) * pow(delta, d[i1]) * pow(tau, (t[i1] - 2))).sum()
				+ i2i.get().mapToDouble(i2 -> n[i2] * t[i2] * (t[i2] - 1) * pow(delta, d[i2]) * pow(tau, (t[i2] - 2)) * exp(-pow(delta, c[i2]))).sum()
				+ i3i.get().mapToDouble(i3 -> n[i3] * pow(delta, d[i3]) * pow(tau, t[i3]) * exp(-alpha[i3] * pow((delta - epsilon[i3]), 2) - beta[i3] * pow((tau - gamma[i3]), 2)) * (pow(((t[i3] / tau) - 2 * beta[i3] * (tau - gamma[i3])), 2) - t[i3] / pow(tau, 2) - 2 * beta[i3])).sum()
				+ i4i.get().mapToDouble(i4 -> n[i4] * delta * (d2Delta.Bi.dtau2[i4] * Psi[i4] + 2 * dDelta.Bi.dtau[i4] * dPsi.dtau[i4] + pow(Delta[i4], b[i4]) * d2Psi.dtau2[i4])).sum();
		residual.Phi.Delta.tau = i1i.get().mapToDouble(i1 -> n[i1] * d[i1] * t[i1] * pow(delta, (d[i1] - 1)) * pow(tau, (t[i1] - 1))).sum()
				+ i2i.get().mapToDouble(i2 -> n[i2] * t[i2] * pow(delta, (d[i2] - 1)) * pow(tau, (t[i2] - 1)) * (d[i2] - c[i2] * pow(delta, c[i2])) * exp(-pow(delta, c[i2]))).sum()
				+ i3i.get().mapToDouble(i3 -> n[i3] * pow(delta, d[i3]) * pow(tau, t[i3]) * exp(-alpha[i3] * pow((delta - epsilon[i3]), 2) - beta[i3] * pow((tau - gamma[i3]), 2)) * ((d[i3] / delta) - 2 * alpha[i3] * (delta - epsilon[i3])) * (t[i3] / tau - 2 * beta[i3] * (tau - gamma[i3]))).sum()
				+ i4i.get().mapToDouble(i4 -> n[i4] * (pow(Delta[i4], b[i4]) * (dPsi.dtau[i4] + delta * d2Psi.Ddelta.dtau[i4]) + delta * dDelta.Bi.ddelta[i4] * dPsi.dtau[i4] + dDelta.Bi.dtau[i4] * (Psi[i4] + delta * dPsi.ddelta[i4]) + d2Delta.Bi.Ddelta.dtau[i4] * delta * Psi[i4])).sum();
//return (p.get()/*get(p)()*/);
	}

	private static void rep(List<Double> oldArray, double toFill, int i) {
		int len = oldArray.size();

		for (int ix = len; ix < len + i; ix++)
			oldArray.add(ix, toFill);
		//double[] newArray = VecMath.fillvec(i, toFill);
		//System.arraycopy(oldArray, i, newArray, 1, 2);
		//return;
	}

	/** define functions idealgas and residual, supplying arguments delta and tau*/
//	double idealgas(Supplier<Double> p) { return IAPWS95_idealgas(p, this.delta, this.tau); }
//	double residual(Supplier<Double> p) { return IAPWS95_residual(p, this.delta, this.tau); }
}
