package net.priveyes.PhaseDiag.H2OCO2;

import static net.priveyes.PhaseDiag.FluidMin.Constants.Rm;
import static net.priveyes.PhaseDiag.Utils.buildVector;
import static net.priveyes.PhaseDiag.Utils.resize;
import static net.priveyes.PhaseDiag.FluidMin.Minimizer.BiVectorFunction;
import static net.priveyes.PhaseDiag.FluidMin.Minimizer.ITMAX;

import net.priveyes.PhaseDiag.FluidMin.solutionProperty;
import net.priveyes.PhaseDiag.FluidMin.Minimizer;

import java.util.concurrent.ExecutionException;
import java.util.function.BiFunction;
import java.util.function.DoubleFunction;

public class  CompH2OCO2 {
/** SUBROUTINE PhaseDiagH2OCO2(pkb, tK, compo)
 returns compositions in mole fraction of H2O- and CO2-rich phases

 PhaseDiagH2OCO2SIGMA(pkb, tK, bestcompo, uncert)
 returns compositions in mole fraction of H2O- and CO2-rich phases
 and associated uncertainties on each value, also in mole fraction

 pkb, tK, w->wCO2NaCl0, aA->aH2O, aB->aCO2*/
	public static double[] PhaseDiagH2OCO2(solutionProperty conditions /*double pkb, double tK*/, double w, double aA, double aB/*, double compo*/) throws ExecutionException, InterruptedException {
//USE nrtype
//USE EVALW

//REAL(DP), INTENT(IN) :: pkb, tK, wCO2NaCl0, aH2O, aCO2
		double[] /*REAL(DP), DIMENSION(2), INTENT(OUT) ::*/ compo = new double[2 /*+ 1*/];
//INTERFACE
//}
//int /*SUBROUTINE*/ zbrak(func,x1,x2,n,xb1,xb2,nb,options) {
		//USE nrtype; USE nrutil, ONLY :arth

		//INTEGER(I4B), INTENT(IN) :: n
		//INTEGER(I4B), INTENT(OUT) :: nb
		//REAL(DP), INTENT(IN) :: x1,x2
		//REAL(DP), DIMENSION(:), POINTER :: xb1,xb2
		//REAL(DP), DIMENSION(4), INTENT(IN) :: options
		//INTERFACE
		//FUNCTION func(x,options)
		//USE nrtype

		//REAL(DP) :: x
		//REAL(DP), DIMENSION(4) :: options
		//REAL(DP) :: func
		//END FUNCTION func
		//END INTERFACE
		//}//END SUBROUTINE zbrak
//}//END INTERFACE

//INTERFACE
		//FUNCTION GBinMix(x,options)
		//USE nrtype

		//REAL(DP):: x
		//REAL(DP), DIMENSION(4) :: options
		//REAL(DP) :: GBinMix
		//END FUNCTION GBinMix
//END INTERFACE

//INTERFACE
		//FUNCTION D2GBinMix(x,options)
		//USE nrtype

		//REAL(DP) :: x
		//REAL(DP), DIMENSION(4) :: options
		//REAL(DP) :: D2GBinMix
		//END FUNCTION D2GBinMix
//END INTERFACE

//INTERFACE
		//FUNCTION rtbis(func,x1,x2,xacc,options)
		//USE nrtype; USE nrutil, ONLY :nrerror

		//REAL(DP) :: x1,x2,xacc
		//REAL(DP), DIMENSION(4) :: options
		//REAL(DP) :: rtbis
		//INTERFACE
		//FUNCTION func(x,options)
		//USE nrtype

		//REAL(DP):: x
		//REAL(DP), DIMENSION(4) :: options
		//REAL(DP) :: func
		//END FUNCTION func
		//END INTERFACE
		//END FUNCTION rtbis
//END INTERFACE

		int /*INTEGER(I4B) ::*/ niter = ITMAX;
		double  arg = -1.0; //! pb, tC,
		double[] /*REAL(DP), DIMENSION(4) ::*/ condition;
		double[] /*REAL(DP), DIMENSION(7) ::*/ optionsGBinSyst;
		int /*INTEGER(I4B) ::*/ nb, nSegments;
		double  nan, x1, x2, Gtol; //!,kk1, kk2, kk3, kk4
		double /*REAL(DP),PARAMETER ::*/ eps = 1.0e-12;
		double[] /*REAL(DP), DIMENSION(:), POINTER ::*/ xb1, xb2;
		double[][] /*REAL(DP), DIMENSION(3,2) ::*/ xs = new double[3 /*+ 1*/][2 /*+ 1*/];
		double[] /*REAL(DP), DIMENSION(3) ::*/ ys = new double[3 /*+ 1*/];

		nan = Math.sqrt(arg);

		condition = new double[]{ conditions.tK, w, aA, aB};

		x1 = eps;
		x2 = 1.0 - x1;
		nSegments = 50;

		Zbrak.zbrak(D2GBinMix, x1, x2, nSegments/*, xb1, xb2, nb*/, condition);
		//! gives nb roots to the equation, located between xb1 and xb2
//		System.out.println("CompZbrak 130 xb2"+ VecMath.toString(Zbrak.xb1));
//		System.out.println("CompZbrak 131 xb2"+ VecMath.toString(Zbrak.xb2));
//		System.out.println("CompZbrak 132 nroot "+ (Zbrak.nroot));

		if (Zbrak.nroot == 2) {//THEN
			//! proceed with Gibbs energy minimization
			//! using the Nelder-Mead simplex method
			//! with high penalties if searching out of values constrained by
			//! 0 < x1 < spinodal (=xb2(1)) and spinodal (=xb1(2)) > x2 > 1
//
			optionsGBinSyst = new double[]{ conditions.tK, w, aA, aB, Zbrak.xb2[1-1] + (Zbrak.xb1[2-1] - Zbrak.xb2[1-1]) / 2.0, Zbrak.xb2[1-1], Zbrak.xb1[2-1]};
//			long start = System.nanoTime();
			//!initial guesses ensuring that initial simplex reflections take place within the allowed boundaries
			xs[1-1]/*,:)*/ = new double[]{2 * Zbrak.xb2[1-1] / 5, Zbrak.xb1[2-1] + (2 * (1 - Zbrak.xb1[2-1]) / 5)}; //! (x1_init_min, x2_init_min)
			xs[2-1]/*,:)*/ = new double[]{3 * Zbrak.xb2[1-1] / 5, Zbrak.xb1[2-1] + (3 * (1 - Zbrak.xb1[2-1]) / 5)}; //! (x1_init_max, x2_init_max)
			xs[3-1]/*,:)*/ = new double[]{2 * Zbrak.xb2[1-1] / 5, Zbrak.xb1[2-1] + (2.5 * (1 - Zbrak.xb1[2-1]) / 5)}; //! (x1_init_min, x2_init_max)

			ys[1-1] = GBinSyst.apply(xs[1 - 1]/*,:)*/, optionsGBinSyst);
			ys[2-1] = GBinSyst.apply(xs[2-1]/*,:)*/, optionsGBinSyst);
			ys[3-1] = GBinSyst.apply(xs[3-1]/*,:)*/, optionsGBinSyst);

			Gtol = 1e-11; //! accuracy
			Minimizer.amoeba/*GBin*/(xs, ys, Gtol, GBinSyst, niter, null, optionsGBinSyst);

//Log.w("Minimize", "Old " + (System.nanoTime() - start));
//			System.out.println("Old                     " + VecMath.toString(xs));
//			start = System.nanoTime();
//			Amoeba amo = new Amoeba(3, 2, xs[2][0], xs[2][1], 150, GBinSyst, optionsGBinSyst);
//			Solution sol = amo.Solve();
//          Log.w("Minimize", "New " + (System.nanoTime() - start));
//			System.out.println("New                     " + VecMath.toString(sol.vector));
//			double[] x = sol.vector;

			compo[1-1] = (xs[1-1][1-1] + xs[2-1][1-1] + xs[3-1][1-1]) / 3.0;
			compo[2-1] = (xs[1-1][2-1] + xs[2-1][2-1] + xs[3-1][2-1]) / 3.0;

			//!The (N + 1) × N matrix p is input. Its N + 1 rows are N-dimensional
			//! vectors that are the vertices of the starting simplex. Also input is the vector y of length
			//! N + 1, whose components must be preinitialized to the values of func evaluated at the
			//! N + 1 vertices (rows) of p;
		} else {//ELSE
			if (Zbrak.nroot > 2) {//THEN
				System.out.println("More than two roots to d²Gmix/dx² - not supposed to happen //!?");
			} else {//ELSE

			}//END IF
			compo/*(:)*/ = new double[]{nan, nan};
		}//END IF//! nb == 2
		return compo;
	}//END SUBROUTINE PhaseDiagH2OCO2
	//INTERFACE
	//void /*SUBROUTINE*/ amoebaGBin(double p, double y, double ftol, double func, double iter, double[] options) {
//! from Numerical Recipes in FORTRAN 90: The Art of PARALLEL Scientific Computing (ISBN 0-521-57439-0)
//! Copyright (C) 1986-1996 by Cambridge University Press.
//! Programs Copyright (C) 1986-1996 by Numerical Recipes Software.
	//USE nrtype; USE nrutil, ONLY :assert_eq, imaxloc,iminloc,nrerror,swap

	//int /*INTEGER(I4B), INTENT(OUT) ::*/ iter;
	//double /*REAL(DP), INTENT(IN) ::*/ ftol;
	//double /*REAL(DP), DIMENSION(:), INTENT(IN) ::*/ options;
	//double /*REAL(DP), DIMENSION(:), INTENT(INOUT) ::*/ y;
	//double /*REAL(DP), DIMENSION(:,:), INTENT(INOUT) ::*/ p;
	//INTERFACE
	//FUNCTION func(xs, options)
	//USE nrtype

	//REAL(DP), DIMENSION(2) :: xs
	//REAL(DP) :: func
	//REAL(DP), DIMENSION(7) :: options
	//END FUNCTION func
//END INTERFACE
	//END SUBROUTINE amoebaGBin
//END INTERFACE

//END MODULE CompH2OCO2

	//static void /*SUBROUTINE*/ zbrak(BiFunction<Double, double[], Double> func, double x1, double x2, int n, double xb1, double xb2, int nb, double[] options) {
////USE nrtype; USE nrutil, ONLY :arth
//
////int /*INTEGER(I4B), INTENT(IN) ::*/ n;
//int /*INTEGER(I4B), INTENT(OUT) ::*/ nb;
////double /*REAL(DP), INTENT(IN) ::*/ x1,x2;
//double[] /*REAL(DP), DIMENSION(:), POINTER ::*/ xb1,xb2;
////double[] /*REAL(DP), DIMENSION(4), INTENT(IN) ::*/ options;
////INTERFACE
////FUNCTION func(x,options)
////USE nrtype
//
////REAL(DP) :: x
////REAL(DP), DIMENSION(4) :: options
////REAL(DP) :: func
////END FUNCTION func
////END INTERFACE
////!Given a function func defined on the interval from x1-x2 subdivide the interval into n
////!equally spaced segments, and search for zero crossings of the function. nb is returned as
////!the number of bracketing pairs xb1(1:nb), xb2(1:nb) that are found. xb1 and xb2 are
////!pointers to arrays of length nb that are dynamically allocated by the routine.
//int /*INTEGER(I4B) ::*/ i;
//double  dx;
//double[] /*REAL(DP), DIMENSION(0:n) ::*/ f,x;
//boolean[] /*LOGICAL(LGT), DIMENSION(1:n) ::*/ mask;
//boolean /*LOGICAL(LGT), SAVE ::*/ init=true/*.true.*/;
//if (init) {//then
//init=false/*.false.*/;
//nullify(xb1,xb2);
//}//end if
//if (associated(xb1)) deallocate(xb1);
//if (associated(xb2)) deallocate(xb2);
//dx=(x2-x1)/n; //! Determine the spacing appropriate to the mesh.
//x=x1+dx*arth(0,1,n+1);
//	for (i = 0; i <= n; i++) {//do i=0,n; //! Evaluate the function at the mesh points.
//f[i] = func.apply(x[i],options);
////!write(*,*) 'x = ',x(i),'f(x) = ', f(i)
//}//end do
//mask=f[1]/*:n)*/*f[0]/*:n-1)*/ <= 0.0; //! Record where the sign changes occur.
//nb=count(mask); //! Number of sign changes.
//allocate(xb1(nb),xb2(nb));
//xb1(1:nb)=pack(x(0:n-1),mask); //! Store the bounds of each bracket.
//xb2(1:nb)=pack(x(1:n),mask);
//}//END SUBROUTINE zbrak
	public static class Zbrak {
		static double[] xb1;
		static double[] xb2;
		static int nroot;

		/**
		 * Given a function or functor fx defined on the interval [x1,x2], subdivide
		 * the interval into n equally spaced segments, and search for zero crossings
		 * of the function. nroot will be set to the number of bracketing pairs found.
		 * If it is positive, the arrays xb1[0..nroot-1] and xb2[0..nroot-1] will be
		 * filled sequentially with any bracketing pairs that are found. On input,
		 * these vectors may have any size, including zero; they will be resized to
		 * nroot.
		 *
		 * @param fx
		 * @param x1
		 * @param x2
		 * @param n
		 * @return Zbrak.xb1
		 * @return Zbrak.xb2
		 * @return Zbrak.nroot
		 */
		static void zbrak(final BiFunction<Double, double[], Double>/*UniVarRealValueFun*/ fx, final double x1, final double x2, final int n, double[] funOptions) {
			int nb = 20;
			xb1 = new double[nb];
			xb2 = new double[nb];
			nroot = 0;
			double dx = (x2 - x1) / n;
			double x = x1;
			double fp = fx.apply(x1, funOptions);
			for (int i = 0; i < n; i++) {
				double fc = fx.apply(x += dx, funOptions);
				if (fc * fp <= 0.0) {
					xb1[nroot] = x - dx;
					xb2[nroot++] = x;
					if (nroot == nb) {
						double[] tempvec1 = buildVector(xb1);
						double[] tempvec2 = buildVector(xb2);
						xb1 = resize(xb1, 2 * nb);
						xb2 = resize(xb2, 2 * nb);
						for (int j = 0; j < nb; j++) {
							xb1[j] = tempvec1[j];
							xb2[j] = tempvec2[j];
						}
						nb *= 2;
					}
				}
				fp = fc;
			}
		}
	}

//void  rtbis(BiFunction<Double, double[], Double> func, double x1, double x2, double xacc, double[] options) {
////USE nrtype; USE nrutil, ONLY :nrerror
//
//			//double  x1,x2,xacc;
//			//double /*REAL(DP), DIMENSION(4) ::*/ options;
//			double  rtbis;
////INTERFACE
////FUNCTION func(x,options)
////USE nrtype
//
////REAL(DP):: x
////REAL(DP), DIMENSION(4) :: options
////REAL(DP) :: func
////END FUNCTION func
////END INTERFACE
//int /*INTEGER(I4B), PARAMETER ::*/ MAXIT=40;
////!Using bisection, find the root of a function func known to lie between x1 and x2. The
////!root, returned as rtbis, will be refined until its accuracy is ±xacc.
////!Parameter: MAXIT is the maximum allowed number of bisections.
//int /*INTEGER(I4B) ::*/ j;
//double  dx,f,fmid,xmid;
//fmid = func.apply(x2,options);
//f = func.apply(x1,options);
//	System.out.println("f(x2) = "+fmid);
//	System.out.println("f(x1) = "+f);
//if (f*fmid >= 0.0) /*call*/ nrerror("rtbis:root must be bracketed");
//if (f < 0.0) {//then //!Orient the search so that f>0 lies at x+dx.
//rtbis=x1;
//dx=x2-x1;
//}else{
//rtbis=x2;
//dx=x1-x2;
//}//end if
//for (j = 1; j <= MAXIT; j++) {//do j=1,MAXIT;//! Bisection loop.
//dx=dx*0.5;
//xmid=rtbis+dx;
//fmid=func.apply(xmid,options);
////!write(*,*) 'fmid = ',fmid
//if (fmid <= 0.0) rtbis=xmid;
//if (Math.abs(dx) < xacc ||/*.or.*/ fmid == 0.0) return; //RETURN
//}//end do
///*call*/ nrerror("rtbis:too many bisections");
//}//END FUNCTION rtbis

	/**
	 * Using bisection, return the root of a function or functor func known to lie
	 * between x1 and x2. The root will be refined until its accuracy is +/-xacc.
	 *
	 * @param func
	 * @param x1
	 * @param x2
	 * @param xacc
	 * @return
	 */
	public static double rtbis(final DoubleFunction/*UniVarRealValueFun*/<Double> func, final double x1, final double x2, final double xacc) {
		final int JMAX = 50;
		double dx, xmid, rtb;
		double f = func.apply(x1);
		double fmid = func.apply(x2);
		if (f * fmid >= 0.0) {
			throw new IllegalArgumentException("Root must be bracketed for bisection in rtbis");
		}
		//rtb = f < 0.0 ? (dx=x2-x1,x1) : (dx=x1-x2,x2);
		rtb = 0;
		if (f < 0.0) {
			dx = x2 - x1;
			rtb = x1;
		} else {
			dx = x1 - x2;
			rtb = x2;
		}
		for (int j = 0; j < JMAX; j++) {
			fmid = func.apply(xmid = rtb + (dx *= 0.5));
			if (fmid <= 0.0) rtb = xmid;
			if (Math.abs(dx) < xacc || fmid == 0.0) return rtb;
		}
		throw new IllegalArgumentException("Too many bisections in rtbis");
	}

	/** Gibbs mixing energy with the ASF formalism along an AB binary
 GBinMix(x,options)
 where x is the mole fraction of X
 options is {tK, wWCO2, alphaA, alphaB }
*/
	BiFunction<Double, double[], Double>  GBinMix/*(double x, double[] options) {*/ = (x, options) -> {
//USE nrtype

		//double  x;
		//double /*REAL(DP), DIMENSION(4) ::*/ options;
		double  GBinMix, tK, w, aA, aB;
		//double /*REAL(DP), PARAMETER ::*/ r = 8.314462e-3;

		tK = options[1-1];
		w = options[2-1];
		aA = options[3-1];
		aB = options[4-1];

		GBinMix = (1 - x) * Rm * tK * Math.log(1 - x) + x * Rm * tK * Math.log(x)
				+ 2 * aA * aB * w * x * (1 - x) / ((aA + aB) * (aA * (1 - x) + aB * x));
		return GBinMix;
	};//END FUNCTION GBinMix

	/** Gibbs energy of the system when unmixing*/
	static BiVectorFunction/*BiFunction<double[], double[], Double>*/  GBinSyst/*(double[] xs, double[] options) {*/
			= (xs, options) -> {

		//double /*REAL(DP), DIMENSION(2)  ::*/ xs;
		//double /*REAL(DP), DIMENSION(7) ::*/ options;
		double  GBinSyst, tK, w, aA, aB, Gx1, Gx2, f1, f2, x1, x2;
		double  xsys, xspin1, xspin2;
		//double /*REAL(DP), PARAMETER ::*/ r = 8.314462e-3;

		tK = options[1-1];
		w = options[2-1];
		aA = options[3-1];
		aB = options[4-1];
		xsys = options[5-1];
		xspin1 = options[6-1];
		xspin2 = options[7-1];

		x1 = xs[1-1];
		x2 = xs[2-1];

		if ((x1 <= 0.0) ||/*.OR.*/(x1 > xspin1) ||/*.OR.*/(x2 < xspin2) ||/*.OR.*/(x2 >= 1.0)) {//THEN
			GBinSyst = 1.0e10; //! large cost for out of bounds
		} else { //ELSE

			f1 = (x2 - xsys) / (x2 - x1);
			f2 = (xsys - x1) / (x2 - x1);

//PP Found Error t -> tK
			Gx1 = (1.0 - x1) * Rm * tK * Math.log(1.0 - x1) + x1 * Rm * tK * Math.log(x1) + 2.0 * aA * aB * w * x1 * (1.0 - x1) / ((aA + aB) * (aA * (1.0 - x1) + aB * x1));
			Gx2 = (1.0 - x2) * Rm * tK * Math.log(1.0 - x2) + x2 * Rm * tK * Math.log(x2) + 2.0 * aA * aB * w * x2 * (1.0 - x2) / ((aA + aB) * (aA * (1.0 - x2) + aB * x2));

			GBinSyst = f1 * Gx1 + f2 * Gx2;
		}//END IF
		return GBinSyst;
	};//END FUNCTION GBinSyst
/** Second derivatice of the Gibbs mixing energy with the ASF formalism along an AB binary
	 D2GBinMix(x,options)
 where x is the mole fraction of X
 options is (/tK, wWCO2, alphaA, alphaB /)*/
	static BiFunction<Double, double[], Double>  D2GBinMix/*(double x, double[] options)*/ = (x, options) -> {

		//double  x;
		double  D2GBinMix, tK, w, aA, aB;
		//double /*REAL(DP), DIMENSION(4) ::*/ options;
		//double /*REAL(DP), PARAMETER ::*/ r = 8.314462e-3;

		tK = options[1-1];
		w = options[2-1];
		aA = options[3-1];
		aB = options[4-1];

		D2GBinMix = (Rm * tK / (x - Math.pow(x, 2.0))) + (4.0 * w * (Math.pow(aA, 2.0)) * (Math.pow(aB, 2.0)) / ((aA + aB) * (Math.pow((aA * (x - 1.0) - aB * x), 3.0))));
		return D2GBinMix;
	};//END FUNCTION D2GBinMix

//	static int iter;
//	static double[] y;
//	static RealMatrix p;
//	static final int /*INTEGER(I4B), PARAMETER ::*/ ITMAX = 50000;
//	static final double /*REAL(DP), PARAMETER ::*/ TITINY = 1.0e-15;
//	static int /*INTEGER(I4B) ::*/ ihi, ndim; //! Global variables.
//	static RealVector /*REAL(DP), DIMENSION(size(p,2)) ::*/ psum;
//	static double ftol;

//static void /*SUBROUTINE*/ amoebaGBin(/*double p, double y,*/ double ftol, BiFunction<Double, double[], Double> func, double iter, double[] options) {
////! from Numerical Recipes in FORTRAN 90: The Art of PARALLEL Scientific Computing (ISBN 0-521-57439-0)
////! Copyright (C) 1986-1996 by Cambridge University Press.
////! Programs Copyright (C) 1986-1996 by Numerical Recipes Software.
////USE nrtype; USE nrutil, ONLY :assert_eq, imaxloc,iminloc,nrerror,swap
//
////int /*INTEGER(I4B), INTENT(OUT) ::*/ iter;
////double /*REAL(DP), INTENT(IN) ::*/ ftol;
////double /*REAL(DP), DIMENSION(:), INTENT(IN) ::*/ options;
////double /*REAL(DP), DIMENSION(:), INTENT(INOUT) ::*/ y;
////double /*REAL(DP), DIMENSION(:,:), INTENT(INOUT) ::*/ p;
////INTERFACE
//   //FUNCTION func(xs, options)
//   //USE nrtype
//
//   //REAL(DP), DIMENSION(2) :: xs
//   //REAL(DP) :: func
//   //REAL(DP), DIMENSION(7) :: options
//   //END FUNCTION func
////END INTERFACE
////! Minimization of the function func in N dimensions by the downhill simplex method of
////! Nelder and Mead. The (N + 1) × N matrix p is input. Its N + 1 rows are N-dimensional
////! vectors that are the vertices of the starting simplex. Also input is the vector y of length
////! N + 1, whose components must be preinitialized to the values of func evaluated at the
////! N + 1 vertices (rows) of p; and ftol the fractional convergence tolerance to be achieved
////! in the function value (n.b.//!). On output, p and y will have been reset to N +1 new points
////! all within ftol of a minimum function value, and iter gives the number of function evaluations taken.
////! Parameters: The maximum allowed number of function evaluations, and a small number.
////int /*INTEGER(I4B) ::*/ ihi,ndim; //! Global variables.
////double[] /*REAL(DP), DIMENSION(size(p,2)) ::*/ psum;
///*call*/ amoebaGBin_private(func);
////CONTAINS
//}//END SUBROUTINE amoebaGBin

//			static void /*SUBROUTINE*/ amoebaGBin_private(BiFunction<Double, double[], Double> func) {
//
//			int /*INTEGER(I4B) ::*/ i,ilo,inhi;
//			double  rtol,ysave,ytry,ytmp;
//			ndim=assert_eq(size(p,2),size(p,1)-1,size(y)-1,"amoeba");
//			iter=0;
//			//psum(:)=sum(p(:,:),dim=1);
//				p.walkInRowOrder(new DefaultRealMatrixPreservingVisitor() {
//					@Override
//					public void visit(int row, int column, double value) {
//						super.visit(row, column, value);
//						psum.addToEntry(row, value);;
//					}
//				});
//			do {//!Iteration loop.
////!write(*,*) '- iter = ', iter
//				ilo=iminloc(y/*(:)*/);   //! Determine which point is the highest (worst),
//			ihi=imaxloc(y/*(:)*/);   //! next-highest, and lowest (best).
//			ytmp=y[ihi];
//			y[ihi]=y[ilo];
//			inhi=imaxloc(y/*(:)*/);
//			y[ihi]=ytmp;
//			rtol=2.0/*_dp*/*Math.abs(y[ihi]-y[ilo])/(Math.abs(y[ihi])+Math.abs(y[ilo])+TITINY);
////!      Compute the fractional range from highest to lowest and return if satisfactory.
//			if (rtol < ftol) {//then      //! If returning, put best point and value in slot
//			/*call*/ swap(y[1],y[ilo]);   //! 1.
//			/*call*/ swap(p[1]/*,:)*/,p[ilo]/*,:)*/);
//			return; //RETURN
//			}//end if
//			if (iter >= ITMAX) /*call*/ nrerror("ITMAX exceeded in amoeba - see GBinsyst//!");
////!   Begin a new iteration. First extrapolate by a factor −1 through the face of the simplex
////!   across from the high point, i.e., reflect the simplex from the high point.
//			ytry=amotry(-1.0/*_dp*/, func);
////!   write(*,*) 'try reflecting'
//			iter=iter+1;
//			if (ytry <= y[ilo]) {//then    //!Gives a result better than the best point, so
//					ytry=amotry(2.0/*_dp*/, func);      //!try an additional extrapolation by a factor of 2.
//			iter=iter+1;
////!      write(*,*) '  reflecting successful'
//   }else if (ytry >= y[inhi]) {//then //!The reflected point is worse than the second
////!   highest, so look for an intermediate lower point, i.e., do a one-dimensional
////!   contraction.
//					ysave=y[ihi];
//			ytry=amotry(0.5/*_dp*/, func);
//			iter=iter+1;
////!      write(*,*) 'try 1D contracting'
//
//			if (ytry >= ysave) {//then
////!      Can't seem to get rid of that high point. Better contract around the lowest
////!      (best) point.
////!         write(*,*) '  failed, contract around lowest'
//			p(:,:)=0.5/*_dp*/*(p/*(:,:)*/+spread(p[ilo]/*,:)*/ ,1,size(p,1)));
//			for (i = 1; i <= ndim+1; i++) { //do i=1,ndim+1
//			if (i != /*/=*/ ilo) y[i]=func.apply(p[i]/*,:)*/, options);
//			}//end do
//				iter=iter+ndim; //!Keep track of function evaluations.
//
//			//psum(:)=sum(p(:,:),dim=1);
//				p.walkInRowOrder(new DefaultRealMatrixPreservingVisitor() {
//					@Override
//					public void visit(int row, int column, double value) {
//						super.visit(row, column, value);
//						psum.addToEntry(row, value);;
//					}
//				});
//			}//end if
//			}//end if
//			}//end do //!Go back for the test of doneness and the next iteration
//		}//END SUBROUTINE amoebaGBin_private

//		static double  amotry(BiFunction<Double, double[], Double> func, double fac) {
//
//				//double /*REAL(DP), INTENT(IN) ::*/ fac;
//				double  amotry;
////!Extrapolates by a factor fac through the face of the simplex across from the high point,
////!tries it, and replaces the high point if the new point is better.
//				double  fac1,fac2;
//				double  ytry;
//				double[] /*REAL(DP), DIMENSION(size(p,2)) ::*/ ptry;
//				fac1=(1.0-fac)/ndim;
//				fac2=fac1-fac;
//				//ptry/*(:)*/=psum/*(:)*/*fac1-p(ihi,:)*fac2;
//			ptry = psum.mapMultiply(fac1).subtract(/*:)*fac1-*/p.getRowVector(ihi).mapMultiply(fac2)/*ihi,:)*fac2*/);
////!write(*,*) 'ptry =',ptry
//				ytry=func.apply(ptry,options);       //!Evaluate the function at the trial point.
//		if (ytry < y[ihi]) {//then   //! If it's better than the highest, then replace
//		y[ihi]=ytry;          //!the highest.
//		psum/*(:)*/=psum/*(:)*/-p(ihi,:)+ptry/*(:)*/;
//		p[ihi]/*,:)*/=ptry/*(:)*/;
//		}//end if
//		amotry=ytry;
//		return amotry;
//	}//END FUNCTION amotry
}
