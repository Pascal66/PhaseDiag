package net.priveyes.PhaseDiag.EVALW;

import static net.priveyes.PhaseDiag.CO2EOS.CO2EOS.CO2DENSITY;
import static net.priveyes.PhaseDiag.FluidMin.Constants.Rm;
import static net.priveyes.PhaseDiag.FluidMin.Constants.mCO2;
import static net.priveyes.PhaseDiag.IAPWS.WATEREOS.waterdensity;

import android.telecom.TelecomManager;

import net.priveyes.PhaseDiag.CHNOZ.IAPWScommon;
import net.priveyes.PhaseDiag.num_rep.VecMath;
import net.priveyes.PhaseDiag.FluidMin.solutionProperty;

import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;

import java.util.concurrent.ExecutionException;
import java.util.stream.IntStream;

/** Module containing functions to calculate wCO2NaCl0=f(p,t,etc)*/
public class  EVALW {

//! SUBROUTINE WH2OCO2(pkb,tK,wCO2NaCl0,aH2O,aCO2) returns wCO2NaCl0,aH2O,aCO2
//! SUBROUTINE SIGMA((/pkb, tK/), varAlpha, varW) gives variance values for aCO2 and wCO2NaCl0

	static public class returnWCO2NACL {
		public static double wCO2NaCl0 = 0, sigwCO2NaCl0 = 0, wCO2NaClpm = 0;
		public static double StDevwCO2NaCl0;

		/**INTENT(OUT):: w->wCO2NaCl0, varWCO2NACL->sigwCO2NaCl0*/
			public static void /*returnWCO2NACL*/ /*SUBROUTINE*/ WCO2NACL(solutionProperty conditions /*double pkb, double tK*/,/*, double wCO2NaCl0, double sigwCO2NaCl0*/double wCO2NaCl0, double sigwCO2NaCl0) {
				//double /*REAL(DP), INTENT(IN)::*/ pkb, tK
				//double /*REAL(DP), INTENT(OUT)::*/ w, varWCO2NACL;
				double /*REAL(DP), PARAMETER ::*/ /*r = 8.314462e-3,*/ done = 50.46293666528664, dtwo = -3.532723355037698, dthree = -0.1525467514747062, dfour = 0.02306707057609801, dfive = -0.4801235974042842;
				double  pbar, x, sumVarDs, sumcovarDs;
				RealMatrix /*double[][]*/ /*REAL(DP), DIMENSION(5,5) ::*/ covDs = MatrixUtils.createRealMatrix(5 /*+ 1*/, 5 /*+ 1*/);
				double[] /*REAL(DP), DIMENSION(5) ::*/ varDS = new double[5 /*+ 1*/], derivParam;
				int /*INTEGER(I4B) ::*/ i, j;

				pbar = 1000.0 * conditions.pkb;
				x = Math.log(pbar);

				returnWCO2NACL.wCO2NaCl0 = Rm * conditions.tK * (done + dtwo * x + dthree * conditions.tK + dfour * conditions.tK * x + dfive * (x * x));

				/* derivative of wCO2NaCl0 with respect to each parameter*/
				derivParam = new double[]{1.0, x, conditions.tK, x * conditions.tK, x * x};

				//! uncertainty
				covDs.setColumn(1-1, /*covDs[][1] =*/ new double[]{8.148403839600984, -0.9419712311046312, -0.028920192801239503, 0.0057604099533277294, -0.1406780439300984});
				covDs.setColumn(2-1, /*covDs[][2] =*/ new double[]{-0.9419712311046312, 0.2686021778583372, 0.0016363211616156454, -0.00039095413120941773, -0.00927680120486655});
				covDs.setColumn(3-1, /*covDs[][3] =*/ new double[]{-0.028920192801239503, 0.0016363211616156454, 0.00012189293050311037, -0.000023752880915650645, 0.0007866971311358896});
				covDs.setColumn(4-1, /*covDs[][4] =*/ new double[]{0.0057604099533277294, -0.00039095413120941773, -0.000023752880915650645, 4.7251812408010245e-6, -0.00015228377217278442});
				covDs.setColumn(5-1, /*covDs[][5] =*/ new double[]{-0.1406780439300984, -0.00927680120486655, 0.0007866971311358896, -0.00015228377217278442, 0.007016249875572129});

				for (i = 0; i < 5; i++) {//DO i = 1, 5
					varDS[i] = covDs.getEntry(i, i)/*[i][i]*/; //! Variance
				}//END DO

		//sumVarDs = sum(varDS*(Math.pow(derivParam, 2)));
				sumVarDs = IntStream.range(0, Math.min(varDS.length, derivParam.length)).mapToDouble(id -> varDS[id] * (Math.pow(derivParam[id], 2))).sum();

				sumcovarDs = 0.0;

				for (i = 0; i < 4; i++) {//DO i = 1, 4
					for (j = i + 1; j < 5; j++) {//DO j = i+1, 5
						sumcovarDs = sumcovarDs + 2.0 * covDs.getEntry(i, j)/*[i][j]*/ * derivParam[i] * derivParam[j];
					}//END DO //!j
				}//END DO //!i

				returnWCO2NACL.sigwCO2NaCl0 = (sumcovarDs + sumVarDs) * (Rm * conditions.tK); //! timed by R*T just like the model
				StDevwCO2NaCl0 = Math.sqrt(returnWCO2NACL.sigwCO2NaCl0) * 2.0;
				//return returnWCO2NACL;
			}//END SUBROUTINE WCO2NACL
		//public static double wCO2NaCl0; //PP
		//public static double sigwCO2NaCl0; //PP
		public static String toDebug() {
			// wCO2NaCl0, sigwCO2NaCl0
			return "returnWCO2NACL wCO2NaCl0 "+ wCO2NaCl0+ " sigwCO2NaCl0 "+sigwCO2NaCl0;
		}
	}

	static public class returnWH2OCO2 {
		//public static double wWCO2 = 0;
		public static double aH2O = 0; //aA? Bizarre
		public static double aCO2 = 0; //aB
		public static double wWCO2 = 0; //wAB
		//public static double aH2O = 0;
		//public static double aCO2 = 0;
		//public static double wWCO2; // PP

		/**
		 * INTENT(OUT):: w, aH2O, aCO2;
		 */
		public static void /*returnWH2OCO2*/ /*SUBROUTINE*/ WH2OCO2(solutionProperty conditions /*double pkb, double tK*/,
				/*, double wCO2NaCl0, double aH2O, double aCO2*/double wWCO2, double aH2O, double aCO2) throws ExecutionException, InterruptedException {
			//double /*REAL(DP), INTENT(IN)::*/ pkb, tK;
			//double /*REAL(DP), INTENT(OUT)::*/ w, aH2O, aCO2;
			double /*REAL(DP)::*/ dc, rhoCO2, dh, pc, ph;
			//double /*REAL(DP), PARAMETER::*/ /*mCO2 = 44.009,*/ r = 8.314462e-3;
			double[] /*REAL(DP), DIMENSION(22) ::*/ ks;
			double /*REAL(DP)::*/ k1, k2, k3, k4, k5, k6, k7, k8, k9, k10, k11, k12, k13;
			double /*REAL(DP)::*/ k14, k15, k16, k17, k18, k19, k20, k21, k22;

			rhoCO2 = CO2DENSITY(conditions /*new double[]{pkb, tK}*//*, rhoCO2*/);
			dc = 1000.0 * mCO2 * rhoCO2;
			//System.out.println("Fortran EVALW WH2OCO2 CO2DENSITY is:        4.3095165252685422E-004 ");
			//System.out.println("This "+rhoCO2);

			dh = waterdensity(conditions /*new double[]{pkb, tK}*//*, dh*/);
			//System.out.println("Fortran EVALW WH2OCO2 waterdensity is:      997.45800046110639      ");
			//System.out.println("this xater denstity");
			pc = dc / (Rm * conditions.tK);
			ph = dh / (Rm * conditions.tK);

			ks = new double[]{ /*double precision::*/1.0, 0.0, 0.00003523673755105464, -1.3771896743353e-7, 2.2365684513434593e-10, -0.027483670839268518, 0.00028638119023642215, -1.393182728887865e-6, 1.7538136320291042e-9, 0.0, 0.0, 0.025909180044638686, 0.0007859663424532188, 9.708440316637061e-6, -4.134776731764802e-8, 5.2823424024797224e-11, 0.0, 0.0, -9.753934004224651e-8, 1.827689653226229e-10, -5.762657149786166e-6, 2.2035377878790666e-11};

			k1 = ks[1-1];
			k2 = ks[2-1];
			k3 = ks[3-1];
			k4 = ks[4-1];
			k5 = ks[5-1];
			k6 = ks[6-1];
			k7 = ks[7-1];
			k8 = ks[8-1];
			k9 = ks[9-1];
			k10 = ks[10-1];
			k11 = ks[11-1];
			k12 = ks[12-1];
			k13 = ks[13-1];
			k14 = ks[14-1];
			k15 = ks[15-1];
			k16 = ks[16-1];
			k17 = ks[17-1];
			k18 = ks[18-1];
			k19 = ks[19-1];
			k20 = ks[20-1];
			k21 = ks[21-1];
			k22 = ks[22-1];

			returnWH2OCO2.wWCO2 = ((k1 + k2 * pc + k3 * (Math.pow(pc, 2)) + k4 * (Math.pow(pc, 3)) + k5 * (Math.pow(pc, 4)) + k6 * ph + k10 * pc * ph + k7 * (Math.pow(ph, 2)) + k11 * (Math.pow(pc, 2)) * (Math.pow(ph, 2)) + k8 * (Math.pow(ph, 3)) + k9 * (Math.pow(ph, 4))) / (k12 + k13 * pc + k14 * (Math.pow(pc, 2)) + k15 * (Math.pow(pc, 3)) + k16 * (Math.pow(pc, 4)) + k17 * ph + k21 * pc * ph + k18 * (Math.pow(ph, 2)) + k22 * (Math.pow(pc, 2)) * (Math.pow(ph, 2)) + k19 * (Math.pow(ph, 3)) + k20 * (Math.pow(ph, 4)))) * Rm * conditions.tK;

			returnWH2OCO2.aH2O = 1.0;

			ks = new double[]{ /*double precision::*/ 1.0, -0.04594561654443757, -0.0007628317926453319, -6.229876106797394e-7, 6.175231929410153e-10, 0.0, -0.00004565879115063656, 1.480961757694035e-7, -1.3100616086506735e-10, 0.00004374471772733658, 2.7152488857567085e-9, -0.6398440569314919, 0.055016465822903, -0.0006789516426367684, 0.0, 0.0, .005294388998244712, -0.000011929398741765865, 0.0, 1.4401961559251803e-11, -0.0002536019541259353, 0.0};

			k1 = ks[1-1];
			k2 = ks[2-1];
			k3 = ks[3-1];
			k4 = ks[4-1];
			k5 = ks[5-1];
			k6 = ks[6-1];
			k7 = ks[7-1];
			k8 = ks[8-1];
			k9 = ks[9-1];
			k10 = ks[10-1];
			k11 = ks[11-1];
			k12 = ks[12-1];
			k13 = ks[13-1];

			k14 = ks[14-1];

			k15 = ks[15-1];
			k16 = ks[16-1];
			k17 = ks[17-1];
			k18 = ks[18-1];
			k19 = ks[19-1];
			k20 = ks[20-1];
			k21 = ks[21-1];
			k22 = ks[22-1];

			returnWH2OCO2.aCO2 = ((k1 + k2 * pc + k3 * (Math.pow(pc, 2)) + k4 * (Math.pow(pc, 3)) + k5 * (Math.pow(pc, 4)) + k6 * ph + k10 * pc * ph + k7 * (Math.pow(ph, 2)) + k11 * (Math.pow(pc, 2)) * (Math.pow(ph, 2)) + k8 * (Math.pow(ph, 3)) + k9 * (Math.pow(ph, 4))) / (k12 + k13 * pc + k14 * (Math.pow(pc, 2)) + k15 * (Math.pow(pc, 3)) + k16 * (Math.pow(pc, 4)) + k17 * ph + k21 * pc * ph + k18 * (Math.pow(ph, 2)) + k22 * (Math.pow(pc, 2)) * (Math.pow(ph, 2)) + k19 * (Math.pow(ph, 3)) + k20 * (Math.pow(ph, 4))));

			//return returnWH2OCO2;
		}//END SUBROUTINE WH2OCO2
		public static String toDebug() {
			//wWCO2, aH2O, aCO2
			return "returnWH2OCO2 wWCO2 "+wWCO2+" aH2O "+ aH2O + " aCO2 "+ aCO2 ;
		}
	}

	static public final class returnWH2ONACL {
		public static double wWNaClpm = 0;
		public static double wWNaClo = 0;
		public static double wNaCloNaClpm = 0;
		public static double aNaClpm = 0;
		public static double aNaCl0 = 0;
		public static double varwWNaClpm = 0;
		public static double varwWNaClo = 0;
		public static double varaNaClpm = 0;
		public static double varaNaCl0 = 0;
		public static double StDevwWNaClo;
		public static double StDevwWNaClpm;
		public static double StDevaNaCl0;
		public static double StDevaNaClpm;

		public static void /*returnWH2ONACL*/ /*SUBROUTINE*/ WH2ONACL(solutionProperty conditions /*double pkb, double tK*/,
				/*,  double wWNaClpm,  double wWNaClo,  double wNaCloNaClpm, double aNaClpm,  double aNaCl0*/
				/*, double varwWNaClpm,  double varwWNaClo,  double varaNaClpm,  double varaNaCl0*/
				double wWNaClpm, double wWNaClo, double wNaCloNaClpm, double aNaClpm, double aNaCl0
				, double varwWNaClpm, double varwWNaClo, double varaNaClpm, double varaNaCl0) {
			//double /*REAL(DP), INTENT(IN)::*/ pkb, tK;
			//double /*REAL(DP), INTENT(OUT)::*/ wWNaClpm, wWNaClo,
			//double wNaCloNaClpm;
			//double /*REAL(DP), INTENT(OUT)::*/ aNaClpm, aNaCl0;
			//double /*REAL(DP)::*/ varwWNaClpm, varwWNaClo, varaNaClpm, varaNaCl0;
			double /*REAL(DP)::*/ alphaNaCloSTP, TdepaNaClo, PdepaNaClo, alphaNaClpmSTP;
			double /*REAL(DP)::*/ TdepaNaClpm, PdepaNaClpm298K, tPdepaNaClpm, wWNaCloSTP;
			double /*REAL(DP)::*/ PdepwWNaClo, PdepwWNaClo2, TdepwWNaClo, wWNaClpmSTP;
			double /*REAL(DP)::*/ PdepwWNaClpm, PdepwWNaClpm2, TdepwWNaClpm, TdepwWNaClpm2;
			double /*REAL(DP)::*/ TDepTerm, sumcovar, sumVar;
			double[][] /*REAL(DP), DIMENSION(16,16)::*/ cov = new double[16][16]; //MatrixUtils.createRealMatrix(16 /*+ 1*/, 16 /*+ 1*/);
			double[] /*REAL(DP), DIMENSION(4) ::*/ derivwWNaClo, derivaNaClpm;
			double[] /*REAL(DP), DIMENSION(5) ::*/ derivwWNaClpm, var = new double[5 /*+ 1*/];
			double[] /*REAL(DP), DIMENSION(3) ::*/ derivaNaCl0;
			int /*INTEGER(I4B) ::*/ i, j;

			alphaNaCloSTP = 1.0;//!C5
			TdepaNaClo = -0.0015267175572519084;//!C7
			PdepaNaClo = -0.0004234981774667679;//!C6
			alphaNaClpmSTP = 0.87;//!C13
			TdepaNaClpm = -0.000559796437659033;//!C14
			PdepaNaClpm298K = 2.799219910543493e-6;//!C15
			tPdepaNaClpm = 6.478357893119607e-8;//!C16
			wWNaCloSTP = -3.6730730758236425; //!C1
			PdepwWNaClo = -0.001376259419628594; //!C2
			PdepwWNaClo2 = 0.0007195535642035902; //!C4
			TdepwWNaClo = 0.008923302365576268;//!C3
			wWNaClpmSTP = -9.847083476502576;//!C8
			PdepwWNaClpm = -0.614; //!C11
			PdepwWNaClpm2 = 0.0316; //!C12
			TdepwWNaClpm = -0.04180222813908087;//!C9
			TdepwWNaClpm2 = 0.00019898894805747683;//!C10
			wNaCloNaClpm = 0.0;

			//! PUBLISHED COVARIANCE MATRIX
			VecMath.setcolumn(cov,1-1, /*cov[][1] =*/ new double[]{1.1151493060848392, -0.00022577182408532612, -0.002837507852331102, -0.00008466360489179576, 0.0005163594406482213, 0.00012212068415476304, -1.3754900002474407e-6, 0.1809695782091573, 0.0017388953455562677, 0.000023686045007566874, -0.04556734335203906, -0.0006688984561209402, -0.0035738783216561305, 9.139892975144099e-6, -8.387691090992425e-7, 3.043054145303593e-8});
			VecMath.setcolumn(cov,2-1, /*cov[][2] =*/ new double[]{-0.00022577182408532612, 1.6312991518687573e-7, 5.744686731049313e-7, 7.064307800442965e-8, 9.226060237154157e-7, -1.235971714603126e-7, -2.22594034369395e-9, -0.00006061700536367417, -3.1396542968772737e-7, -5.161187521947547e-9, 8.118968123116216e-6, 1.206888454968378e-7, 2.4404580513834016e-6, -6.301072605524422e-9, -8.19342880343794e-11, 1.535541082243574e-11});
			VecMath.setcolumn(cov,3-1, /*cov[][3] =*/ new double[]{-0.002837507852331102, 5.744686731049313e-7, 7.220065464697081e-6, 2.1542244307325374e-7, -1.314077067193679e-6, -3.107271600013972e-7, 3.5005176977592553e-9, -0.0004604825936668293, -4.424486634825189e-6, -6.027201557080935e-8, 0.00011595265926285716, 1.7020059106222914e-6, 9.09386473320159e-6, -2.3256846820729423e-8, 2.1342625097899845e-9, -7.743579296932729e-11});
			VecMath.setcolumn(cov,4-1, /*cov[][4] =*/ new double[]{-0.00008466360489179576, 7.064307800442965e-8, 2.1542244307325374e-7, 3.10974473554911e-8, 5.338180474308305e-7, -5.4878472113385136e-8, -1.3027489324833588e-9, -0.000023752023388136722, -1.0908899009369157e-7, -1.978559858941007e-9, 1.2777924247915054e-6, 5.195687622776674e-8, 1.0047873300395265e-6, -2.598398944633851e-9, -4.120886573007705e-11, 7.383767332788938e-12});
			VecMath.setcolumn(cov,5-1, /*cov[][5] =*/ new double[]{0.0005163594406482213, 9.226060237154157e-7, -1.314077067193679e-6, 5.338180474308305e-7, 0.0002332015810276684, -7.208021422924908e-7, -5.840106699604749e-7, -0.0008564971127667994, 0.00004807384165217395, -6.871333873517793e-7, -0.0015987368577075113, -1.618881422924902e-6, 0.000053754940711462544, -1.4381397035573117e-7, -5.519440316205539e-9, -2.6926483882971806e-10});
			VecMath.setcolumn(cov,6-1, /*cov[][6] =*/ new double[]{0.00012212068415476304, -1.235971714603126e-7, -3.107271600013972e-7, -5.4878472113385136e-8, -7.208021422924908e-7, 1.0129955143317724e-7, 1.7283921587809876e-9, 0.00003668140348807714, 1.9986360209041606e-7, 1.9836671716462044e-9, -3.20026563994269e-6, -9.086778368920936e-8, -1.6591203992094879e-6, 4.301002501336556e-9, 8.697762538291302e-11, -2.324207678021581e-11});
			VecMath.setcolumn(cov,7-1, /*cov[][7] =*/ new double[]{-1.3754900002474407e-6, -2.22594034369395e-9, 3.5005176977592553e-9, -1.3027489324833588e-9, -5.840106699604749e-7, 1.7283921587809876e-9, 1.4755942977843876e-9, 2.5893676443583933e-6, -1.2213588832350644e-7, 1.7513266403925887e-9, 4.059240744576285e-6, 1.8024340143280729e-9, -1.6039211067193676e-7, 4.1595098472233117e-10, 1.4345759146537547e-11, 7.222456868084071e-13});
			VecMath.setcolumn(cov,8-1, /*cov[][8] =*/ new double[]{0.1809695782091573, -0.00006061700536367417, -0.0004604825936668293, -0.000023752023388136722, -0.0008564971127667994, 0.00003668140348807714, 2.5893676443583933e-6, 0.11741114105793271, 0.00034543498911485495, 5.120651859180137e-6, 0.00028892842271258345, -0.00039660344457726225, -0.005563003823646252, 0.000013847714449621944, 1.0288173147473821e-7, -5.212516906286768e-10});
			VecMath.setcolumn(cov,9-1, /*cov[][9] =*/ new double[]{0.0017388953455562677, -3.1396542968772737e-7, -4.424486634825189e-6, -1.0908899009369157e-7, 0.00004807384165217395, 1.9986360209041606e-7, -1.2213588832350644e-7, 0.00034543498911485495, 0.00004940380102658622, -6.199035865752444e-7, -0.00019410577347542688, -5.102063790862862e-6, -7.1345282233201594e-6, 1.801204647800555e-8, -2.4436831309622068e-9, -1.904124558060742e-12});
			VecMath.setcolumn(cov,10-1, /*cov[][10] =*/ new double[]{0.000023686045007566874, -5.161187521947547e-9, -6.027201557080935e-8, -1.978559858941007e-9, -6.871333873517793e-7, 1.9836671716462044e-9, 1.7513266403925887e-9, 5.120651859180137e-6, -6.199035865752444e-7, 9.840140312271186e-9, 1.1077914784247023e-6, 3.79646766817984e-8, -1.7198709683794483e-7, 4.3544455994169917e-10, 5.0221865834249e-12, 1.731236341309418e-12});
			VecMath.setcolumn(cov,11-1, /*cov[][11] =*/ new double[]{-0.04556734335203906, 8.118968123116216e-6, 0.00011595265926285716, 1.2777924247915054e-6, -0.0015987368577075113, -3.20026563994269e-6, 4.059240744576285e-6, 0.00028892842271258345, -0.00019410577347542688, 1.1077914784247023e-6, 0.06571427116580632, 2.2446783446638363e-6, -0.0002676252371541505, 6.876251640088909e-7, 5.530628803017786e-8, 1.802229564369773e-9});
			VecMath.setcolumn(cov,12-1, /*cov[][12] =*/ new double[]{-0.0006688984561209402, 1.206888454968378e-7, 1.7020059106222914e-6, 5.195687622776674e-8, -1.618881422924902e-6, -9.086778368920936e-8, 1.8024340143280729e-9, -0.00039660344457726225, -5.102063790862862e-6, 3.79646766817984e-8, 2.2446783446638363e-6, 3.3740323387351776e-6, 0.00001810848418972334, -4.4340242630316196e-8, 5.618385763300387e-10, 2.148308085744496e-11});
			VecMath.setcolumn(cov,13-1, /*cov[][13] =*/ new double[]{-0.0035738783216561305, 2.4404580513834016e-6, 9.09386473320159e-6, 1.0047873300395265e-6, 0.000053754940711462544, -1.6591203992094879e-6, -1.6039211067193676e-7, -0.005563003823646252, -7.1345282233201594e-6, -1.7198709683794483e-7, -0.0002676252371541505, 0.00001810848418972334, 0.000295256916996048, -7.335814130434789e-7, -1.1332550988142303e-8, 2.897379644771648e-10});
			VecMath.setcolumn(cov,14-1, /*cov[][14] =*/ new double[]{9.139892975144099e-6, -6.301072605524422e-9, -2.3256846820729423e-8, -2.598398944633851e-9, -1.4381397035573117e-7, 4.301002501336556e-9, 4.1595098472233117e-10, 0.000013847714449621944, 1.801204647800555e-8, 4.3544455994169917e-10, 6.876251640088909e-7, -4.4340242630316196e-8, -7.335814130434789e-7, 1.8291094099271934e-9, 2.8609925848116596e-11, -7.651015704470928e-13});
			VecMath.setcolumn(cov,15-1, /*cov[][15] =*/ new double[]{-8.387691090992425e-7, -8.19342880343794e-11, 2.1342625097899845e-9, -4.120886573007705e-11, -5.519440316205539e-9, 8.697762538291302e-11, 1.4345759146537547e-11, 1.0288173147473821e-7, -2.4436831309622068e-9, 5.0221865834249e-12, 5.530628803017786e-8, 5.618385763300387e-10, -1.1332550988142303e-8, 2.8609925848116596e-11, 2.0131785205450593e-12, -5.28517455182299e-14});
			VecMath.setcolumn(cov,16-1, /*cov[][16] =*/ new double[]{3.043054145303593e-8, 1.535541082243574e-11, -7.743579296932729e-11, 7.383767332788938e-12, -2.6926483882971806e-10, -2.324207678021581e-11, 7.222456868084071e-13, -5.212516906286768e-10, -1.904124558060742e-12, 1.731236341309418e-12, 1.802229564369773e-9, 2.148308085744496e-11, 2.897379644771648e-10, -7.651015704470928e-13, -5.28517455182299e-14, 2.6611228538790855e-14});

	//! ==================  W H2O - NaCl0  ==================
	//! wWNaCloSTP = C1
	//! PdepwWNaClo = C2
	//! PdepwWNaClo2 = C4
	//! TdepwWNaClo = C3

			returnWH2ONACL.wWNaClo = wWNaCloSTP + conditions.pkb * PdepwWNaClo + PdepwWNaClo2 * (Math.pow(conditions.pkb, 2)) + (conditions.tK - IAPWScommon.Temperature.defaut) * TdepwWNaClo;

			derivwWNaClo = new double[]{1.0, conditions.pkb, conditions.pkb * conditions.pkb, conditions.tK - IAPWScommon.Temperature.defaut}; //! d(W)/dC5 to d(W)/dC16 == 0

			for (i = 0; i < 4; i++) {//DO i = 1, 4
				var[i] = cov[i][i]; //! Variance
			}//END DO

	//sumVar = sum(var(1:4)*(Math.pow(derivwWNaClo, 2)));
			sumVar = IntStream.range(1-1, 4-1).mapToDouble(id -> var[id] * (Math.pow(derivwWNaClo[id], 2))).sum();

			sumcovar = 0.0;

			for (i = 0; i < 3; i++) {//DO i = 1, 3
				for (j = i + 1; j < 4; j++) {//DO j = i+1, 4
					sumcovar = sumcovar + 2.0 * cov[i][j] * derivwWNaClo[i] * derivwWNaClo[j];
				}//END DO //!j
			}//END DO //!i

			returnWH2ONACL.varwWNaClo = sumcovar + sumVar;
			StDevwWNaClo = Math.sqrt(returnWH2ONACL.varwWNaClo) * 2.0;

	//! ==================  W H2O - NaClpm  ==================
	//! wWNaClpmSTP = C8
	//! TdepwWNaClpm = C9
	//! TdepwWNaClpm2 = C10
	//! PdepwWNaClpm = C11
	//! PdepwWNaClpm2 = C12

			if (conditions.tK < 400) {//THEN
				TDepTerm = (conditions.tK - IAPWScommon.Temperature.defaut) * TdepwWNaClpm + TdepwWNaClpm2 * (Math.pow((conditions.tK - IAPWScommon.Temperature.defaut), 2));
			} else {//ELSE
				TDepTerm = (400 - IAPWScommon.Temperature.defaut) * TdepwWNaClpm + TdepwWNaClpm2 * (Math.pow((400 - IAPWScommon.Temperature.defaut), 2));
			}//END IF

			returnWH2ONACL.wWNaClpm = wWNaClpmSTP + TDepTerm + conditions.pkb * PdepwWNaClpm + PdepwWNaClpm2 * Math.pow(conditions.pkb, 2);

			derivwWNaClpm = new double[]{1.0, conditions.tK - IAPWScommon.Temperature.defaut, (conditions.tK - IAPWScommon.Temperature.defaut) * (conditions.tK - IAPWScommon.Temperature.defaut), conditions.pkb, conditions.pkb * conditions.pkb};

			for (i = 0; i < 5; i++) {//DO i = 1, 5
				var[i] = cov[i+7][i+7]; //! Variance from i = 8 to 12
			}//END DO

	//sumVar = sum(var(1:5)*(Math.pow(derivwWNaClpm, 2)));
			sumVar = IntStream.range(1-1, 5-1).mapToDouble(id -> var[id] * (Math.pow(derivwWNaClpm[id], 2))).sum();

			sumcovar = 0.0;

			for (i = 0; i < 4; i++) {//DO i = 1, 4
				for (j = i + 1; j < 5; j++) {//DO j = i+1, 5
					sumcovar = sumcovar + 2.0 * cov[i+7][j+7] * derivwWNaClpm[i] * derivwWNaClpm[j];
				}//END DO //!j
			}//END DO //!i

			returnWH2ONACL.varwWNaClpm = sumcovar + sumVar;
			StDevwWNaClpm = Math.sqrt(returnWH2ONACL.varwWNaClpm) * 2.0;

	//! ==================  alpha NaCl0  ==================
	//! alphaNaCloSTP = C5
	//! TdepaNaClo = C7
	//! PdepaNaClo = C6

			returnWH2ONACL.aNaCl0 = alphaNaCloSTP + PdepaNaClo * conditions.pkb + TdepaNaClo * (conditions.tK - IAPWScommon.Temperature.defaut);

			derivaNaCl0 = new double[]{1.0, conditions.pkb, conditions.tK - IAPWScommon.Temperature.defaut};

			for (i = 0; i < 3; i++) {//DO i = 1, 3
				var[i] = cov[i+4][i+4]; //! Variance from i = 5 to 7;
			}//END DO

	//sumVar = sum(var(1:3)*(derivaNaCl0**2));
			sumVar = IntStream.range(1-1, 3-1).mapToDouble(id -> var[id] * (Math.pow(derivaNaCl0[id], 2))).sum();

			sumcovar = 0.0;

			for (i = 0; i < 2; i++) {//DO i = 1, 2
				for (j = i + 1; j < 3; j++) {//DO j = i+1, 3
					sumcovar = sumcovar + 2.0 * cov[i+4][j+4] * derivaNaCl0[i] * derivaNaCl0[j];
				}//END DO //!j
			}//END DO //!i

			returnWH2ONACL.varaNaCl0 = sumcovar + sumVar;
			StDevaNaCl0 = Math.sqrt(returnWH2ONACL.varaNaCl0) * 2.0;

	//! ==================  alpha NaClpm  ==================
	//! alphaNaClpmSTP = C13
	//! TdepaNaClpm = C14
	//! PdepaNaClpm298K = C15
	//! tPdepaNaClpm = C16

			returnWH2ONACL.aNaClpm = alphaNaClpmSTP + (PdepaNaClpm298K + tPdepaNaClpm * (conditions.tK - IAPWScommon.Temperature.defaut)) * conditions.pkb + TdepaNaClpm * (conditions.tK - 298.15);

			derivaNaClpm = new double[]{1.0, conditions.tK - IAPWScommon.Temperature.defaut, conditions.pkb, conditions.pkb * (conditions.tK - IAPWScommon.Temperature.defaut)};

			for (i = 0; i < 4; i++) {//DO i = 1, 4
				var[i] = cov[i+12][i+12]; //! Variance from i = 13 to 16
			}//END DO

	//sumVar = sum(var(1:4)*(Math.pow(derivaNaClpm, 2)));
			sumVar = IntStream.range(1-1, 4-1).mapToDouble(id -> var[id] * (Math.pow(derivaNaClpm[id], 2))).sum();

			sumcovar = 0.0;

			for (i = 0; i < 3; i++) {//DO i = 1, 3
				for (j = i + 1; j < 4; j++) {//DO j = i+1, 4
					sumcovar = sumcovar + 2.0 * cov[i+12][j+12] * derivaNaClpm[i] * derivaNaClpm[j];
				}//END DO //!j
			}//END DO //!i

			returnWH2ONACL.varaNaClpm = sumcovar + sumVar;
			StDevaNaClpm = Math.sqrt(returnWH2ONACL.varaNaClpm) * 2.0;
			//return returnWH2ONACL;
		}//END SUBROUTINE WH2ONACL

		public static String toDebug() {
			//wWNaClpm, wWNaClo, wNaCloNaClpm, aNaClpm, aNaCl0, varwWNaClpm, varwWNaClo, varaNaClpm, varaNaCl0
			return "returnWH2ONACL wWNaClpm "+wWNaClpm+ " wWNaClo "+wWNaClo+" wNaCloNaClpm "+wNaCloNaClpm+" aNaClpm "+aNaClpm
					+" aNaCl0 "+aNaCl0+" varwWNaClpm "+varwWNaClpm+" varwWNaClo "+varwWNaClo+" varaNaClpm "+varaNaClpm+" varaNaCl0 "+varaNaCl0;
		}
	}

	public static final class returnSigma {
		public static double varAlpha = 0;
		public static double varW = 0;
		public static double StDevAlpha;
		public static double StDevW;

		public static void/*returnSigma*/ /*SUBROUTINE*/ SIGMA(solutionProperty conditions,/*,  double varAlpha,  double varW*/double varAlpha, double varW) throws ExecutionException, InterruptedException {
			//double /*REAL(DP), DIMENSION(2), INTENT(IN)::*/ conditions;
			//double /*REAL(DP), INTENT(OUT)::*/ varAlpha, varW;
			double /*REAL(DP)::*/ a, b, densCO2, densH2O/*, pkb, tK*/;
			//double /*REAL(DP), PARAMETER::*/ /*mCO2 = 44.009/,*/ r = 8.314462e-3;
			double[] /*REAL(DP), DIMENSION(22) ::*/ derivW = new double[22 /*+ 1*/], derivAlpha = new double[22 /*+ 1*/],
					varJS = new double[22 /*+ 1*/], varKS = new double[22 /*+ 1*/];
			double[][] /*REAL(DP), DIMENSION(22,22) ::*/ covJs = new double[22 /*+ 1*/][22 /*+ 1*/],
					covKs = new double[22 /*+ 1*/][22 /*+ 1*/];
			double /*REAL(DP)::*/ k1, k2, k3, k4, k5, k6, k7, k8, k9, k10, k11, k12, k13;
			double /*REAL(DP)::*/ k14, k15, k16, k17, k18, k19, k20, k21, k22;
			double /*REAL(DP)::*/ sumVarJs, sumcovarJs;
			double /*REAL(DP)::*/ sumVarKs, sumcovarKs;
			int /*INTEGER(I4B) ::*/ n, i, j;

			//pkb = conditions[1-1];
			//tK = conditions[2-1];

			densCO2 = CO2DENSITY(conditions /*new double[]{pkb, tK}*//*, densCO2*/);

			densH2O = waterdensity(conditions/*new double[]{pkb, tK}*//*, densH2O*/);

			a = 1000.0 * mCO2 * densCO2 / (Rm * conditions.tK);
			b = densH2O / (Rm * conditions.tK);

	//! Alpha

	//! below are k values to calculate sigma(alpha)
			k1 = 1.0;
			k2 = -0.04594561654443757;
			k3 = -0.0007628317926453319;
			k4 = -6.229876106797394e-7;
			k5 = 6.175231929410153e-10;
			k6 = 0.0;
			k7 = -0.00004565879115063656;
			k8 = 1.480961757694035e-7;
			k9 = -1.3100616086506735e-10;
			k10 = 0.00004374471772733658;
			k11 = 2.7152488857567085e-9;
			k12 = -0.6398440569314919;
			k13 = 0.055016465822903;
			k14 = -0.0006789516426367684;
			k15 = 0.0;
			k16 = 0.0;
			k17 = 0.005294388998244712;
			k18 = -0.000011929398741765865;
			k19 = 0.0;
			k20 = 1.4401961559251803e-11;
			k21 = -0.0002536019541259353;
			k22 = 0.0;
	//! above are k values to calculate sigma(alpha)

			covJs[1-1] = new double[]{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
			covJs[2-1] = new double[]{0.0, 0.00008700442722848148, 1.187917239454692e-6, 9.893292652125718e-10, -9.972739628304394e-13, 0.0, 1.8242913882158343e-8, -6.7968603778517e-11, 6.687757648570212e-14, -9.859425509158558e-8, -4.144626416916596e-12, 0.0005241033457666806, -0.00008269490874915328, 1.0780593227727057e-6, 0.0, 0.0, -2.754472869557125e-6, 3.0765264831456628e-9, 0.0, 6.684297080754853e-15, 4.0177289557592427e-7, 0.0};
			covJs[3-1] = new double[]{0.0, 1.187917239454692e-6, 2.5391004855737534e-8, 1.8953170104228833e-11, -1.7873722517316272e-14, 0.0, 1.3484778873749931e-10, -3.184569401531231e-13, 1.305061171568384e-16, -9.899203735150656e-10, -9.056496205063034e-14, 0.00002418644101421857, -1.7780131814299088e-6, 2.226219008907051e-8, 0.0, 0.0, -1.615518163425486e-7, 2.920567260201366e-10, 0.0, -1.9182723267793767e-16, 7.609117486001093e-9, 0.0};
			covJs[4-1] = new double[]{0.0, 9.893292652125718e-10, 1.8953170104228833e-11, 1.542978212431385e-14, -1.5368891802153635e-17, 0.0, 1.109468701321308e-13, -2.6809879074313077e-16, 1.171843094699598e-19, -8.979744847444579e-13, -6.776134697924772e-17, 1.7588296554722317e-8, -1.3422093962032605e-9, 1.686274323797913e-11, 0.0, 0.0, -1.1930686964136449e-10, 2.1947393740955287e-13, 0.0, -1.5127914241322698e-19, 5.86832538722592e-12, 0.0};
			covJs[5-1] = new double[]{0.0, -9.972739628304394e-13, -1.7873722517316272e-14, -1.5368891802153635e-17, 1.5922663038060778e-20, 0.0, -1.1680566952293606e-16, 3.0058985965519014e-19, -1.5937409181119185e-22, 9.451640446008518e-16, 6.389510011947817e-20, -1.624250395185382e-11, 1.2795934819339366e-12, -1.6054930997030154e-14, 0.0, 0.0, 1.1035454136471097e-13, -2.0310130743495036e-16, 0.0, 1.3685132169681811e-22, -5.6702710627776545e-15, 0.0};
			covJs[6-1] = new double[]{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
			covJs[7-1] = new double[]{0.0, 1.8242913882158343e-8, 1.3484778873749931e-10, 1.109468701321308e-13, -1.1680566952293606e-16, 0.0, 9.012413552506807e-12, -4.1002700764047836e-14, 4.803133553639355e-17, -2.5372488781639454e-11, -4.3900204071950254e-16, -2.925512909082208e-7, -7.96732812963571e-9, 1.2527667748442956e-10, 0.0, 0.0, 2.816309121094401e-9, -7.165968064530996e-12, 0.0, 1.2679169156451492e-17, 5.5043896310157805e-11, 0.0};
			covJs[8-1] = new double[]{0.0, -6.7968603778517e-11, -3.184569401531231e-13, -2.6809879074313077e-16, 3.0058985965519014e-19, 0.0, -4.1002700764047836e-14, 1.9238784946548836e-16, -2.3057342413618422e-19, 1.0359983804839943e-13, 9.463817886730917e-19, 1.7591957302790836e-9, 1.3825799551860286e-11, -3.0761489206251044e-13, 0.0, 0.0, -1.615388314266034e-11, 3.9800162836785765e-14, 0.0, -6.63187886826714e-20, -1.5274654646301703e-13, 0.0};
			covJs[9-1] = new double[]{0.0, 6.687757648570212e-14, 1.305061171568384e-16, 1.171843094699598e-19, -1.5937409181119185e-22, 0.0, 4.803133553639355e-17, -2.3057342413618422e-19, 2.8091268992448577e-22, -1.1129527233655664e-16, -2.433040986560048e-22, -2.4187328214211455e-12, 2.305277027735754e-15, 1.445912890144321e-16, 0.0, 0.0, 2.1776235086735574e-14, -5.289777272678801e-17, 0.0, 8.553761208356334e-23, 9.732423261657484e-17, 0.0};
			covJs[10-1] = new double[]{0.0, -9.859425509158558e-8, -9.899203735150656e-10, -8.979744847444579e-13, 9.451640446008518e-16, 0.0, -2.5372488781639454e-11, 1.0359983804839943e-13, -1.1129527233655664e-16, 1.287197236276622e-10, 3.3817143569708048e-15, 2.511038748932649e-7, 6.508716941981025e-8, -9.26630925231334e-10, 0.0, 0.0, -3.1791007502584974e-9, 9.580343096141975e-12, 0.0, -2.3276913259254226e-17, -3.5820274614665497e-10, 0.0};
			covJs[11-1] = new double[]{0.0, -4.144626416916596e-12, -9.056496205063034e-14, -6.776134697924772e-17, 6.389510011947817e-20, 0.0, -4.3900204071950254e-16, 9.463817886730917e-19, -2.433040986560048e-22, 3.3817143569708048e-15, 3.2401735039561317e-19, -8.684234750493175e-11, 6.313383068936992e-12, -7.927865943734428e-14, 0.0, 0.0, 5.838763570395035e-13, -1.0648121052704859e-15, 0.0, 7.346525603836833e-22, -2.6901809155307685e-14, 0.0};
			covJs[12-1] = new double[]{0.0, 0.0005241033457666806, 0.00002418644101421857, 1.7588296554722317e-8, -1.624250395185382e-11, 0.0, -2.925512909082208e-7, 1.7591957302790836e-9, -2.4187328214211455e-12, 2.511038748932649e-7, -8.684234750493175e-11, 0.05950096110975415, -0.0020292122177833385, 0.00002107312449121733, 0.0, 0.0, -0.00045109314316743045, 9.476155302991713e-7, 0.0, -1.1454442725441085e-12, 7.511538412892692e-6, 0.0};
			covJs[13-1] = new double[]{0.0, -0.00008269490874915328, -1.7780131814299088e-6, -1.3422093962032605e-9, 1.2795934819339366e-12, 0.0, -7.96732812963571e-9, 1.3825799551860286e-11, 2.305277027735754e-15, 6.508716941981025e-8, 6.313383068936992e-12, -0.0020292122177833385, 0.00013034340576153892, -1.5665486607907473e-6, 0.0, 0.0, 0.000013872419001950095, -2.5879036036862326e-8, 0.0, 2.0253343882976754e-14, -5.534567090872433e-7, 0.0};
			covJs[14-1] = new double[]{0.0, 1.0780593227727057e-6, 2.226219008907051e-8, 1.686274323797913e-11, -1.6054930997030154e-14, 0.0, 1.2527667748442956e-10, -3.0761489206251044e-13, 1.445912890144321e-16, -9.26630925231334e-10, -7.927865943734428e-14, 0.00002107312449121733, -1.5665486607907473e-6, 1.959643822661632e-8, 0.0, 0.0, -1.4066944498283406e-7, 2.5392015116971645e-10, 0.0, -1.6373967907899004e-16, 6.746302507049011e-9, 0.0};
			covJs[15-1] = new double[]{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
			covJs[16-1] = new double[]{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
			covJs[17-1] = new double[]{0.0, -2.754472869557125e-6, -1.615518163425486e-7, -1.1930686964136449e-10, 1.1035454136471097e-13, 0.0, 2.816309121094401e-9, -1.615388314266034e-11, 2.1776235086735574e-14, -3.1791007502584974e-9, 5.838763570395035e-13, -0.00045109314316743045, 0.000013872419001950095, -1.4066944498283406e-7, 0.0, 0.0, 3.5145121278898784e-6, -7.581212007985856e-9, 0.0, 9.711318644712393e-15, -5.022126002081105e-8, 0.0};
			covJs[18-1] = new double[]{0.0, 3.0765264831456628e-9, 2.920567260201366e-10, 2.1947393740955287e-13, -2.0310130743495036e-16, 0.0, -7.165968064530996e-12, 3.9800162836785765e-14, -5.289777272678801e-17, 9.580343096141975e-12, -1.0648121052704859e-15, 9.476155302991713e-7, -2.5879036036862326e-8, 2.5392015116971645e-10, 0.0, 0.0, -7.581212007985856e-9, 1.6759240207449485e-11, 0.0, -2.2566627239526905e-17, 9.092406017414963e-11, 0.0};
			covJs[19-1] = new double[]{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
			covJs[20-1] = new double[]{0.0, 6.684297080754853e-15, -1.9182723267793767e-16, -1.5127914241322698e-19, 1.3685132169681811e-22, 0.0, 1.2679169156451492e-17, -6.63187886826714e-20, 8.553761208356334e-23, -2.3276913259254226e-17, 7.346525603836833e-22, -1.1454442725441085e-12, 2.0253343882976754e-14, -1.6373967907899004e-16, 0.0, 0.0, 9.711318644712393e-15, -2.2566627239526905e-17, 0.0, 3.342163704080494e-23, -5.861281960919388e-17, 0.0};
			covJs[21-1] = new double[]{0.0, 4.0177289557592427e-7, 7.609117486001093e-9, 5.86832538722592e-12, -5.6702710627776545e-15, 0.0, 5.5043896310157805e-11, -1.5274654646301703e-13, 9.732423261657484e-17, -3.5820274614665497e-10, -2.6901809155307685e-14, 7.511538412892692e-6, -5.534567090872433e-7, 6.746302507049011e-9, 0.0, 0.0, -5.022126002081105e-8, 9.092406017414963e-11, 0.0, -5.861281960919388e-17, 2.4298907719232355e-9, 0.0};
			covJs[22-1] = new double[]{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

			derivAlpha[1-1] = 1.0 / (k12 + a * k13 + Math.pow(a, 2) * k14 + Math.pow(a, 3) * k15 + Math.pow(a, 4) * k16 + b * k17 + Math.pow(b, 2) * k18 + Math.pow(b, 3) * k19 + Math.pow(b, 4) * k20 + a * b * k21 + Math.pow(a, 2) * Math.pow(b, 2) * k22);
			derivAlpha[2-1] = a / (k12 + a * k13 + Math.pow(a, 2) * k14 + Math.pow(a, 3) * k15 + Math.pow(a, 4) * k16 + b * k17 + Math.pow(b, 2) * k18 + Math.pow(b, 3) * k19 + Math.pow(b, 4) * k20 + a * b * k21 + Math.pow(a, 2) * Math.pow(b, 2) * k22);
			derivAlpha[3-1] = Math.pow(a, 2) / (k12 + a * k13 + Math.pow(a, 2) * k14 + Math.pow(a, 3) * k15 + Math.pow(a, 4) * k16 + b * k17 + Math.pow(b, 2) * k18 + Math.pow(b, 3) * k19 + Math.pow(b, 4) * k20 + a * b * k21 + Math.pow(a, 2) * Math.pow(b, 2) * k22);
			derivAlpha[4-1] = Math.pow(a, 3) / (k12 + a * k13 + Math.pow(a, 2) * k14 + Math.pow(a, 3) * k15 + Math.pow(a, 4) * k16 + b * k17 + Math.pow(b, 2) * k18 + Math.pow(b, 3) * k19 + Math.pow(b, 4) * k20 + a * b * k21 + Math.pow(a, 2) * Math.pow(b, 2) * k22);
			derivAlpha[5-1] = Math.pow(a, 4) / (k12 + a * k13 + Math.pow(a, 2) * k14 + Math.pow(a, 3) * k15 + Math.pow(a, 4) * k16 + b * k17 + Math.pow(b, 2) * k18 + Math.pow(b, 3) * k19 + Math.pow(b, 4) * k20 + a * b * k21 + Math.pow(a, 2) * Math.pow(b, 2) * k22);
			derivAlpha[6-1] = b / (k12 + a * k13 + Math.pow(a, 2) * k14 + Math.pow(a, 3) * k15 + Math.pow(a, 4) * k16 + b * k17 + Math.pow(b, 2) * k18 + Math.pow(b, 3) * k19 + Math.pow(b, 4) * k20 + a * b * k21 + Math.pow(a, 2) * Math.pow(b, 2) * k22);
			derivAlpha[7-1] = Math.pow(b, 2) / (k12 + a * k13 + Math.pow(a, 2) * k14 + Math.pow(a, 3) * k15 + Math.pow(a, 4) * k16 + b * k17 + Math.pow(b, 2) * k18 + Math.pow(b, 3) * k19 + Math.pow(b, 4) * k20 + a * b * k21 + Math.pow(a, 2) * Math.pow(b, 2) * k22);
			derivAlpha[8-1] = Math.pow(b, 3) / (k12 + a * k13 + Math.pow(a, 2) * k14 + Math.pow(a, 3) * k15 + Math.pow(a, 4) * k16 + b * k17 + Math.pow(b, 2) * k18 + Math.pow(b, 3) * k19 + Math.pow(b, 4) * k20 + a * b * k21 + Math.pow(a, 2) * Math.pow(b, 2) * k22);
			derivAlpha[9-1] = Math.pow(b, 4) / (k12 + a * k13 + Math.pow(a, 2) * k14 + Math.pow(a, 3) * k15 + Math.pow(a, 4) * k16 + b * k17 + Math.pow(b, 2) * k18 + Math.pow(b, 3) * k19 + Math.pow(b, 4) * k20 + a * b * k21 + Math.pow(a, 2) * Math.pow(b, 2) * k22);
			derivAlpha[10-1] = (a * b) / (k12 + a * k13 + Math.pow(a, 2) * k14 + Math.pow(a, 3) * k15 + Math.pow(a, 4) * k16 + b * k17 + Math.pow(b, 2) * k18 + Math.pow(b, 3) * k19 + Math.pow(b, 4) * k20 + a * b * k21 + Math.pow(a, 2) * Math.pow(b, 2) * k22);
			derivAlpha[11-1] = (Math.pow(a, 2) * Math.pow(b, 2)) / (k12 + a * k13 + Math.pow(a, 2) * k14 + Math.pow(a, 3) * k15 + Math.pow(a, 4) * k16 + b * k17 + Math.pow(b, 2) * k18 + Math.pow(b, 3) * k19 + Math.pow(b, 4) * k20 + a * b * k21 + Math.pow(a, 2) * Math.pow(b, 2) * k22);
			derivAlpha[12-1] = (-k1 - a * b * k10 - Math.pow(a, 2) * Math.pow(b, 2) * k11 - a * k2 - Math.pow(a, 2) * k3 - Math.pow(a, 3) * k4 - Math.pow(a, 4) * k5 - b * k6 - Math.pow(b, 2) * k7 - Math.pow(b, 3) * k8 - Math.pow(b, 4) * k9) / Math.pow((k12 + a * k13 + Math.pow(a, 2) * k14 + Math.pow(a, 3) * k15 + Math.pow(a, 4) * k16 + b * k17 + Math.pow(b, 2) * k18 + Math.pow(b, 3) * k19 + Math.pow(b, 4) * k20 + a * b * k21 + Math.pow(a, 2) * Math.pow(b, 2) * k22), 2);
			derivAlpha[13-1] = -((a * (k1 + a * b * k10 + Math.pow(a, 2) * Math.pow(b, 2) * k11 + a * k2 + Math.pow(a, 2) * k3 + Math.pow(a, 3) * k4 + Math.pow(a, 4) * k5 + b * k6 + Math.pow(b, 2) * k7 + Math.pow(b, 3) * k8 + Math.pow(b, 4) * k9)) / Math.pow((k12 + a * k13 + Math.pow(a, 2) * k14 + Math.pow(a, 3) * k15 + Math.pow(a, 4) * k16 + b * k17 + Math.pow(b, 2) * k18 + Math.pow(b, 3) * k19 + Math.pow(b, 4) * k20 + a * b * k21 + Math.pow(a, 2) * Math.pow(b, 2) * k22), 2));
			derivAlpha[14-1] = -((Math.pow(a, 2) * (k1 + a * b * k10 + Math.pow(a, 2) * Math.pow(b, 2) * k11 + a * k2 + Math.pow(a, 2) * k3 + Math.pow(a, 3) * k4 + Math.pow(a, 4) * k5 + b * k6 + Math.pow(b, 2) * k7 + Math.pow(b, 3) * k8 + Math.pow(b, 4) * k9)) / Math.pow((k12 + a * k13 + Math.pow(a, 2) * k14 + Math.pow(a, 3) * k15 + Math.pow(a, 4) * k16 + b * k17 + Math.pow(b, 2) * k18 + Math.pow(b, 3) * k19 + Math.pow(b, 4) * k20 + a * b * k21 + Math.pow(a, 2) * Math.pow(b, 2) * k22), 2));
			derivAlpha[15-1] = -((Math.pow(a, 3) * (k1 + a * b * k10 + Math.pow(a, 2) * Math.pow(b, 2) * k11 + a * k2 + Math.pow(a, 2) * k3 + Math.pow(a, 3) * k4 + Math.pow(a, 4) * k5 + b * k6 + Math.pow(b, 2) * k7 + Math.pow(b, 3) * k8 + Math.pow(b, 4) * k9)) / Math.pow((k12 + a * k13 + Math.pow(a, 2) * k14 + Math.pow(a, 3) * k15 + Math.pow(a, 4) * k16 + b * k17 + Math.pow(b, 2) * k18 + Math.pow(b, 3) * k19 + Math.pow(b, 4) * k20 + a * b * k21 + Math.pow(a, 2) * Math.pow(b, 2) * k22), 2));
			derivAlpha[16-1] = -((Math.pow(a, 4) * (k1 + a * b * k10 + Math.pow(a, 2) * Math.pow(b, 2) * k11 + a * k2 + Math.pow(a, 2) * k3 + Math.pow(a, 3) * k4 + Math.pow(a, 4) * k5 + b * k6 + Math.pow(b, 2) * k7 + Math.pow(b, 3) * k8 + Math.pow(b, 4) * k9)) / Math.pow((k12 + a * k13 + Math.pow(a, 2) * k14 + Math.pow(a, 3) * k15 + Math.pow(a, 4) * k16 + b * k17 + Math.pow(b, 2) * k18 + Math.pow(b, 3) * k19 + Math.pow(b, 4) * k20 + a * b * k21 + Math.pow(a, 2) * Math.pow(b, 2) * k22), 2));
			derivAlpha[17-1] = -((b * (k1 + a * b * k10 + Math.pow(a, 2) * Math.pow(b, 2) * k11 + a * k2 + Math.pow(a, 2) * k3 + Math.pow(a, 3) * k4 + Math.pow(a, 4) * k5 + b * k6 + Math.pow(b, 2) * k7 + Math.pow(b, 3) * k8 + Math.pow(b, 4) * k9)) / Math.pow((k12 + a * k13 + Math.pow(a, 2) * k14 + Math.pow(a, 3) * k15 + Math.pow(a, 4) * k16 + b * k17 + Math.pow(b, 2) * k18 + Math.pow(b, 3) * k19 + Math.pow(b, 4) * k20 + a * b * k21 + Math.pow(a, 2) * Math.pow(b, 2) * k22), 2));
			derivAlpha[18-1] = -((Math.pow(b, 2) * (k1 + a * b * k10 + Math.pow(a, 2) * Math.pow(b, 2) * k11 + a * k2 + Math.pow(a, 2) * k3 + Math.pow(a, 3) * k4 + Math.pow(a, 4) * k5 + b * k6 + Math.pow(b, 2) * k7 + Math.pow(b, 3) * k8 + Math.pow(b, 4) * k9)) / Math.pow((k12 + a * k13 + Math.pow(a, 2) * k14 + Math.pow(a, 3) * k15 + Math.pow(a, 4) * k16 + b * k17 + Math.pow(b, 2) * k18 + Math.pow(b, 3) * k19 + Math.pow(b, 4) * k20 + a * b * k21 + Math.pow(a, 2) * Math.pow(b, 2) * k22), 2));
			derivAlpha[19-1] = -((Math.pow(b, 3) * (k1 + a * b * k10 + Math.pow(a, 2) * Math.pow(b, 2) * k11 + a * k2 + Math.pow(a, 2) * k3 + Math.pow(a, 3) * k4 + Math.pow(a, 4) * k5 + b * k6 + Math.pow(b, 2) * k7 + Math.pow(b, 3) * k8 + Math.pow(b, 4) * k9)) / Math.pow((k12 + a * k13 + Math.pow(a, 2) * k14 + Math.pow(a, 3) * k15 + Math.pow(a, 4) * k16 + b * k17 + Math.pow(b, 2) * k18 + Math.pow(b, 3) * k19 + Math.pow(b, 4) * k20 + a * b * k21 + Math.pow(a, 2) * Math.pow(b, 2) * k22), 2));
			derivAlpha[20-1] = -((Math.pow(b, 4) * (k1 + a * b * k10 + Math.pow(a, 2) * Math.pow(b, 2) * k11 + a * k2 + Math.pow(a, 2) * k3 + Math.pow(a, 3) * k4 + Math.pow(a, 4) * k5 + b * k6 + Math.pow(b, 2) * k7 + Math.pow(b, 3) * k8 + Math.pow(b, 4) * k9)) / Math.pow((k12 + a * k13 + Math.pow(a, 2) * k14 + Math.pow(a, 3) * k15 + Math.pow(a, 4) * k16 + b * k17 + Math.pow(b, 2) * k18 + Math.pow(b, 3) * k19 + Math.pow(b, 4) * k20 + a * b * k21 + Math.pow(a, 2) * Math.pow(b, 2) * k22), 2));
			derivAlpha[21-1] = -((a * b * (k1 + a * b * k10 + Math.pow(a, 2) * Math.pow(b, 2) * k11 + a * k2 + Math.pow(a, 2) * k3 + Math.pow(a, 3) * k4 + Math.pow(a, 4) * k5 + b * k6 + Math.pow(b, 2) * k7 + Math.pow(b, 3) * k8 + Math.pow(b, 4) * k9)) / Math.pow((k12 + a * k13 + Math.pow(a, 2) * k14 + Math.pow(a, 3) * k15 + Math.pow(a, 4) * k16 + b * k17 + Math.pow(b, 2) * k18 + Math.pow(b, 3) * k19 + Math.pow(b, 4) * k20 + a * b * k21 + Math.pow(a, 2) * Math.pow(b, 2) * k22), 2));
			derivAlpha[22-1] = -((Math.pow(a, 2) * Math.pow(b, 2) * (k1 + a * b * k10 + Math.pow(a, 2) * Math.pow(b, 2) * k11 + a * k2 + Math.pow(a, 2) * k3 + Math.pow(a, 3) * k4 + Math.pow(a, 4) * k5 + b * k6 + Math.pow(b, 2) * k7 + Math.pow(b, 3) * k8 + Math.pow(b, 4) * k9)) / Math.pow((k12 + a * k13 + Math.pow(a, 2) * k14 + Math.pow(a, 3) * k15 + Math.pow(a, 4) * k16 + b * k17 + Math.pow(b, 2) * k18 + Math.pow(b, 3) * k19 + Math.pow(b, 4) * k20 + a * b * k21 + Math.pow(a, 2) * Math.pow(b, 2) * k22), 2));

			n = covJs.length; //size(covJs,1);
			for (i = 0; i < n; i++) {//DO i = 1, n
				varJS[i] = covJs[i][i]; //! Variance
			}//END DO

	//sumVarJs = sum(varJS*(Math.pow(derivAlpha, 2)));
			sumVarJs = IntStream.range(0, Math.min(varJS.length, derivAlpha.length)).mapToDouble(id -> varJS[id] * (Math.pow(derivAlpha[id], 2))).sum();

			sumcovarJs = 0.0;

			for (i = 0; i < n - 1; i++) {//DO i = 1, n-1
				for (j = i+1; j < n; j++) {//DO j = i+1, n
					sumcovarJs = sumcovarJs + 2.0 * covJs[i][j] * derivAlpha[i] * derivAlpha[j];
				}//END DO //!j
			}//END DO //!i

			returnSigma.varAlpha = sumcovarJs + sumVarJs;
			StDevAlpha = Math.sqrt(returnSigma.varAlpha) * 2.0;

	//!//! below are k values to calculate sigma(wCO2NaCl0)
			k1 = 1.0;
			k2 = 0.0;
			k3 = 0.00003523673755105464;
			k4 = -1.3771896743353e-7;
			k5 = 2.2365684513434593e-10;
			k6 = -0.027483670839268518;
			k7 = 0.00028638119023642215;
			k8 = -1.393182728887865e-6;
			k9 = 1.7538136320291042e-9;
			k10 = 0.0;
			k11 = 0.0;
			k12 = 0.025909180044638686;
			k13 = 0.0007859663424532188;
			k14 = 9.708440316637061e-6;
			k15 = -4.134776731764802e-8;
			k16 = 5.2823424024797224e-11;
			k17 = 0.0;
			k18 = 0.0;
			k19 = -9.753934004224651e-8;
			k20 = 1.827689653226229e-10;
			k21 = -5.762657149786166e-6;
			k22 = 2.2035377878790666e-11;
	//! above are k values to calculate sigma(wCO2NaCl0)

			covKs[1-1] = new double[]{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
			covKs[2-1] = new double[]{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
			covKs[3-1] = new double[]{0.0, 0.0, 3.144870210050728e-11, -1.3197892581892327e-13, 2.09051465999424e-16, -4.739152053168983e-9, 7.136978509113833e-11, -4.4787903938670775e-13, 4.902468572150991e-16, 0.0, 0.0, -4.105811390748958e-8, 6.683977097930173e-10, 9.47626179264513e-12, -3.839228251424929e-14, 4.983354020034719e-17, 0.0, 0.0, -4.579749843838573e-14, 6.920648809638936e-17, -4.629021799465371e-12, 1.1484160747248122e-17};
			covKs[4-1] = new double[]{0.0, 0.0, -1.3197892581892327e-13, 5.688570312002272e-16, -9.100583862915963e-19, 1.7848305711875114e-11, -2.673562464827037e-13, 1.6798865680528183e-15, -1.7571259776073632e-18, 0.0, 0.0, 1.635381344242722e-10, -2.914776021943031e-12, -3.9284558456144976e-14, 1.6194833592958327e-16, -2.144435437655438e-19, 0.0, 0.0, 1.7669214617158198e-16, -2.5680863781614134e-19, 1.9268238632927595e-14, -4.543054503340974e-20};
			covKs[5-1] = new double[]{0.0, 0.0, 2.09051465999424e-16, -9.100583862915963e-19, 1.479581563315502e-21, -2.7910833218691427e-14, 4.188714633664097e-16, -2.6304010195292156e-18, 2.7264132424336082e-21, 0.0, 0.0, -2.510818523271467e-13, 4.580115633405086e-15, 6.216936522910385e-17, -2.570941527633963e-19, 3.470446960827511e-22, 0.0, 0.0, -2.7797084324218077e-19, 4.0099054200211814e-22, -3.031478501460643e-17, 6.899608830407701e-23};
			covKs[6-1] = new double[]{0.0, 0.0, -4.739152053168983e-9, 1.7848305711875114e-11, -2.7910833218691427e-14, 1.1275746725229698e-6, -1.6816347437494097e-8, 1.044509436019393e-10, -1.3139019294221302e-13, 0.0, 0.0, 9.144379994344241e-6, -6.917792727628069e-8, -1.4839198853067661e-9, 5.6534784835244446e-12, -7.035861440116395e-15, 0.0, 0.0, 9.590310815526327e-12, -1.6626895595114772e-14, 6.637028045115861e-10, -2.0108371146758794e-15};
			covKs[7-1] = new double[]{0.0, 0.0, 7.136978509113833e-11, -2.673562464827037e-13, 4.188714633664097e-16, -1.6816347437494097e-8, 2.5759440157696286e-10, -1.6050214089900044e-12, 2.0071064779963992e-15, 0.0, 0.0, -1.1584825319894709e-7, 1.1558160384456412e-9, 2.233763280666157e-11, -8.494403234910911e-14, 1.0493447849360293e-16, 0.0, 0.0, -1.4845432328078644e-13, 2.55609879651849e-16, -1.0507796219304502e-11, 3.1852842168666644e-17};
			covKs[8-1] = new double[]{0.0, 0.0, -4.4787903938670775e-13, 1.6798865680528183e-15, -2.6304010195292156e-18, 1.044509436019393e-10, -1.6050214089900044e-12, 1.0010990000071953e-14, -1.2494751697423876e-17, 0.0, 0.0, 7.078968211196926e-10, -7.370476847599884e-12, -1.3994986640659553e-13, 5.324983048097947e-16, -6.573142525190056e-19, 0.0, 0.0, 9.286201190578022e-16, -1.5953358224608643e-18, 6.623273070649685e-14, -2.0020855205270424e-19};
			covKs[9-1] = new double[]{0.0, 0.0, 4.902468572150991e-16, -1.7571259776073632e-18, 2.7264132424336082e-21, -1.3139019294221302e-13, 2.0071064779963992e-15, -1.2494751697423876e-17, 1.6138668130039806e-20, 0.0, 0.0, -9.037382649981693e-13, 6.790814062016557e-15, 1.5412326583240766e-16, -5.74036120315106e-19, 6.957778417738394e-22, 0.0, 0.0, -1.1258415849165494e-18, 2.009223532454372e-21, -7.04899638913527e-17, 2.309651923112212e-22};
			covKs[10-1] = new double[]{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
			covKs[11-1] = new double[]{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
			covKs[12-1] = new double[]{0.0, 0.0, -4.105811390748958e-8, 1.635381344242722e-10, -2.510818523271467e-13, 9.144379994344241e-6, -1.1584825319894709e-7, 7.078968211196926e-10, -9.037382649981693e-13, 0.0, 0.0, 0.00014032606935112047, -3.73833059911942e-7, -1.2444768135940305e-8, 4.8878997635641125e-11, -6.332703560904003e-14, 0.0, 0.0, 6.338839865583167e-11, -1.120270140805199e-13, 4.3001921621235136e-9, -1.2319212929794598e-14};
			covKs[13-1] = new double[]{0.0, 0.0, 6.683977097930173e-10, -2.914776021943031e-12, 4.580115633405086e-15, -6.917792727628069e-8, 1.1558160384456412e-9, -7.370476847599884e-12, 6.790814062016557e-15, 0.0, 0.0, -3.73833059911942e-7, 1.9900018932701613e-8, 1.897077735103052e-10, -7.991984368545677e-13, 1.0308076563381624e-15, 0.0, 0.0, -8.329585964984207e-13, 1.0972591912499802e-15, -1.1181172732813755e-10, 2.6588219168237756e-16};
			covKs[14-1] = new double[]{0.0, 0.0, 9.47626179264513e-12, -3.9284558456144976e-14, 6.216936522910385e-17, -1.4839198853067661e-9, 2.233763280666157e-11, -1.3994986640659553e-13, 1.5412326583240766e-16, 0.0, 0.0, -1.2444768135940305e-8, 1.897077735103052e-10, 2.953545921673349e-12, -1.1799915435367197e-14, 1.517624199670661e-17, 0.0, 0.0, -1.42350064317376e-14, 2.162866091330022e-17, -1.4005806107949425e-12, 3.4934651180985065e-18};
			covKs[15-1] = new double[]{0.0, 0.0, -3.839228251424929e-14, 1.6194833592958327e-16, -2.570941527633963e-19, 5.6534784835244446e-12, -8.494403234910911e-14, 5.324983048097947e-16, -5.74036120315106e-19, 0.0, 0.0, 4.8878997635641125e-11, -7.991984368545677e-13, -1.1799915435367197e-14, 4.7831714645017034e-17, -6.215535961632931e-20, 0.0, 0.0, 5.489386496846857e-17, -8.185062227630339e-20, 5.6488114062771746e-15, -1.3832692300124735e-20};
			covKs[16-1] = new double[]{0.0, 0.0, 4.983354020034719e-17, -2.144435437655438e-19, 3.470446960827511e-22, -7.035861440116395e-15, 1.0493447849360293e-16, -6.573142525190056e-19, 6.957778417738394e-22, 0.0, 0.0, -6.332703560904003e-14, 1.0308076563381624e-15, 1.517624199670661e-17, -6.215535961632931e-20, 8.305713437352831e-23, 0.0, 0.0, -6.84792304720761e-20, 1.0053754860508961e-22, -7.1799381817577e-18, 1.668336015545776e-23};
			covKs[17-1] = new double[]{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
			covKs[18-1] = new double[]{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
			covKs[19-1] = new double[]{0.0, 0.0, -4.579749843838573e-14, 1.7669214617158198e-16, -2.7797084324218077e-19, 9.590310815526327e-12, -1.4845432328078644e-13, 9.286201190578022e-16, -1.1258415849165494e-18, 0.0, 0.0, 6.338839865583167e-11, -8.329585964984207e-13, -1.42350064317376e-14, 5.489386496846857e-17, -6.84792304720761e-20, 0.0, 0.0, 8.83348583296838e-17, -1.4715200044519891e-19, 6.8947626728178935e-15, -1.97504369730598e-20};
			covKs[20-1] = new double[]{0.0, 0.0, 6.920648809638936e-17, -2.5680863781614134e-19, 4.0099054200211814e-22, -1.6626895595114772e-14, 2.55609879651849e-16, -1.5953358224608643e-18, 2.009223532454372e-21, 0.0, 0.0, -1.120270140805199e-13, 1.0972591912499802e-15, 2.162866091330022e-17, -8.185062227630339e-20, 1.0053754860508961e-22, 0.0, 0.0, -1.4715200044519891e-19, 2.552595966381773e-22, -1.01602967867454e-17, 3.1296508951839165e-23};
			covKs[21-1] = new double[]{0.0, 0.0, -4.629021799465371e-12, 1.9268238632927595e-14, -3.031478501460643e-17, 6.637028045115861e-10, -1.0507796219304502e-11, 6.623273070649685e-14, -7.04899638913527e-17, 0.0, 0.0, 4.3001921621235136e-9, -1.1181172732813755e-10, -1.4005806107949425e-12, 5.6488114062771746e-15, -7.1799381817577e-18, 0.0, 0.0, 6.8947626728178935e-15, -1.01602967867454e-17, 7.391750975307157e-13, -1.858730148698278e-18};
			covKs[22-1] = new double[]{0.0, 0.0, 1.1484160747248122e-17, -4.543054503340974e-20, 6.899608830407701e-23, -2.0108371146758794e-15, 3.1852842168666644e-17, -2.0020855205270424e-19, 2.309651923112212e-22, 0.0, 0.0, -1.2319212929794598e-14, 2.6588219168237756e-16, 3.4934651180985065e-18, -1.3832692300124735e-20, 1.668336015545776e-23, 0.0, 0.0, -1.97504369730598e-20, 3.1296508951839165e-23, -1.858730148698278e-18, 5.315630087352819e-24};

			derivW[1-1] = 1.0 / (k12 + a * k13 + Math.pow(a, 2) * k14 + Math.pow(a, 3) * k15 + Math.pow(a, 4) * k16 + b * k17 + Math.pow(b, 2) * k18 + Math.pow(b, 3) * k19 + Math.pow(b, 4) * k20 + a * b * k21 + Math.pow(a, 2) * Math.pow(b, 2) * k22);
			derivW[2-1] = a / (k12 + a * k13 + Math.pow(a, 2) * k14 + Math.pow(a, 3) * k15 + Math.pow(a, 4) * k16 + b * k17 + Math.pow(b, 2) * k18 + Math.pow(b, 3) * k19 + Math.pow(b, 4) * k20 + a * b * k21 + Math.pow(a, 2) * Math.pow(b, 2) * k22);
			derivW[3-1] = Math.pow(a, 2) / (k12 + a * k13 + Math.pow(a, 2) * k14 + Math.pow(a, 3) * k15 + Math.pow(a, 4) * k16 + b * k17 + Math.pow(b, 2) * k18 + Math.pow(b, 3) * k19 + Math.pow(b, 4) * k20 + a * b * k21 + Math.pow(a, 2) * Math.pow(b, 2) * k22);
			derivW[4-1] = Math.pow(a, 3) / (k12 + a * k13 + Math.pow(a, 2) * k14 + Math.pow(a, 3) * k15 + Math.pow(a, 4) * k16 + b * k17 + Math.pow(b, 2) * k18 + Math.pow(b, 3) * k19 + Math.pow(b, 4) * k20 + a * b * k21 + Math.pow(a, 2) * Math.pow(b, 2) * k22);
			derivW[5-1] = Math.pow(a, 4) / (k12 + a * k13 + Math.pow(a, 2) * k14 + Math.pow(a, 3) * k15 + Math.pow(a, 4) * k16 + b * k17 + Math.pow(b, 2) * k18 + Math.pow(b, 3) * k19 + Math.pow(b, 4) * k20 + a * b * k21 + Math.pow(a, 2) * Math.pow(b, 2) * k22);
			derivW[6-1] = b / (k12 + a * k13 + Math.pow(a, 2) * k14 + Math.pow(a, 3) * k15 + Math.pow(a, 4) * k16 + b * k17 + Math.pow(b, 2) * k18 + Math.pow(b, 3) * k19 + Math.pow(b, 4) * k20 + a * b * k21 + Math.pow(a, 2) * Math.pow(b, 2) * k22);
			derivW[7-1] = Math.pow(b, 2) / (k12 + a * k13 + Math.pow(a, 2) * k14 + Math.pow(a, 3) * k15 + Math.pow(a, 4) * k16 + b * k17 + Math.pow(b, 2) * k18 + Math.pow(b, 3) * k19 + Math.pow(b, 4) * k20 + a * b * k21 + Math.pow(a, 2) * Math.pow(b, 2) * k22);
			derivW[8-1] = Math.pow(b, 3) / (k12 + a * k13 + Math.pow(a, 2) * k14 + Math.pow(a, 3) * k15 + Math.pow(a, 4) * k16 + b * k17 + Math.pow(b, 2) * k18 + Math.pow(b, 3) * k19 + Math.pow(b, 4) * k20 + a * b * k21 + Math.pow(a, 2) * Math.pow(b, 2) * k22);
			derivW[9-1] = Math.pow(b, 4) / (k12 + a * k13 + Math.pow(a, 2) * k14 + Math.pow(a, 3) * k15 + Math.pow(a, 4) * k16 + b * k17 + Math.pow(b, 2) * k18 + Math.pow(b, 3) * k19 + Math.pow(b, 4) * k20 + a * b * k21 + Math.pow(a, 2) * Math.pow(b, 2) * k22);
			derivW[10-1] = (a * b) / (k12 + a * k13 + Math.pow(a, 2) * k14 + Math.pow(a, 3) * k15 + Math.pow(a, 4) * k16 + b * k17 + Math.pow(b, 2) * k18 + Math.pow(b, 3) * k19 + Math.pow(b, 4) * k20 + a * b * k21 + Math.pow(a, 2) * Math.pow(b, 2) * k22);
			derivW[11-1] = (Math.pow(a, 2) * Math.pow(b, 2)) / (k12 + a * k13 + Math.pow(a, 2) * k14 + Math.pow(a, 3) * k15 + Math.pow(a, 4) * k16 + b * k17 + Math.pow(b, 2) * k18 + Math.pow(b, 3) * k19 + Math.pow(b, 4) * k20 + a * b * k21 + Math.pow(a, 2) * Math.pow(b, 2) * k22);
			derivW[12-1] = (-k1 - a * b * k10 - Math.pow(a, 2) * Math.pow(b, 2) * k11 - a * k2 - Math.pow(a, 2) * k3 - Math.pow(a, 3) * k4 - Math.pow(a, 4) * k5 - b * k6 - Math.pow(b, 2) * k7 - Math.pow(b, 3) * k8 - Math.pow(b, 4) * k9) / Math.pow((k12 + a * k13 + Math.pow(a, 2) * k14 + Math.pow(a, 3) * k15 + Math.pow(a, 4) * k16 + b * k17 + Math.pow(b, 2) * k18 + Math.pow(b, 3) * k19 + Math.pow(b, 4) * k20 + a * b * k21 + Math.pow(a, 2) * Math.pow(b, 2) * k22), 2);
			derivW[13-1] = -((a * (k1 + a * b * k10 + Math.pow(a, 2) * Math.pow(b, 2) * k11 + a * k2 + Math.pow(a, 2) * k3 + Math.pow(a, 3) * k4 + Math.pow(a, 4) * k5 + b * k6 + Math.pow(b, 2) * k7 + Math.pow(b, 3) * k8 + Math.pow(b, 4) * k9)) / Math.pow((k12 + a * k13 + Math.pow(a, 2) * k14 + Math.pow(a, 3) * k15 + Math.pow(a, 4) * k16 + b * k17 + Math.pow(b, 2) * k18 + Math.pow(b, 3) * k19 + Math.pow(b, 4) * k20 + a * b * k21 + Math.pow(a, 2) * Math.pow(b, 2) * k22), 2));
			derivW[14-1] = -((Math.pow(a, 2) * (k1 + a * b * k10 + Math.pow(a, 2) * Math.pow(b, 2) * k11 + a * k2 + Math.pow(a, 2) * k3 + Math.pow(a, 3) * k4 + Math.pow(a, 4) * k5 + b * k6 + Math.pow(b, 2) * k7 + Math.pow(b, 3) * k8 + Math.pow(b, 4) * k9)) / Math.pow((k12 + a * k13 + Math.pow(a, 2) * k14 + Math.pow(a, 3) * k15 + Math.pow(a, 4) * k16 + b * k17 + Math.pow(b, 2) * k18 + Math.pow(b, 3) * k19 + Math.pow(b, 4) * k20 + a * b * k21 + Math.pow(a, 2) * Math.pow(b, 2) * k22), 2));
			derivW[15-1] = -((Math.pow(a, 3) * (k1 + a * b * k10 + Math.pow(a, 2) * Math.pow(b, 2) * k11 + a * k2 + Math.pow(a, 2) * k3 + Math.pow(a, 3) * k4 + Math.pow(a, 4) * k5 + b * k6 + Math.pow(b, 2) * k7 + Math.pow(b, 3) * k8 + Math.pow(b, 4) * k9)) / Math.pow((k12 + a * k13 + Math.pow(a, 2) * k14 + Math.pow(a, 3) * k15 + Math.pow(a, 4) * k16 + b * k17 + Math.pow(b, 2) * k18 + Math.pow(b, 3) * k19 + Math.pow(b, 4) * k20 + a * b * k21 + Math.pow(a, 2) * Math.pow(b, 2) * k22), 2));
			derivW[16-1] = -((Math.pow(a, 4) * (k1 + a * b * k10 + Math.pow(a, 2) * Math.pow(b, 2) * k11 + a * k2 + Math.pow(a, 2) * k3 + Math.pow(a, 3) * k4 + Math.pow(a, 4) * k5 + b * k6 + Math.pow(b, 2) * k7 + Math.pow(b, 3) * k8 + Math.pow(b, 4) * k9)) / Math.pow((k12 + a * k13 + Math.pow(a, 2) * k14 + Math.pow(a, 3) * k15 + Math.pow(a, 4) * k16 + b * k17 + Math.pow(b, 2) * k18 + Math.pow(b, 3) * k19 + Math.pow(b, 4) * k20 + a * b * k21 + Math.pow(a, 2) * Math.pow(b, 2) * k22), 2));
			derivW[17-1] = -((b * (k1 + a * b * k10 + Math.pow(a, 2) * Math.pow(b, 2) * k11 + a * k2 + Math.pow(a, 2) * k3 + Math.pow(a, 3) * k4 + Math.pow(a, 4) * k5 + b * k6 + Math.pow(b, 2) * k7 + Math.pow(b, 3) * k8 + Math.pow(b, 4) * k9)) / Math.pow((k12 + a * k13 + Math.pow(a, 2) * k14 + Math.pow(a, 3) * k15 + Math.pow(a, 4) * k16 + b * k17 + Math.pow(b, 2) * k18 + Math.pow(b, 3) * k19 + Math.pow(b, 4) * k20 + a * b * k21 + Math.pow(a, 2) * Math.pow(b, 2) * k22), 2));
			derivW[18-1] = -((Math.pow(b, 2) * (k1 + a * b * k10 + Math.pow(a, 2) * Math.pow(b, 2) * k11 + a * k2 + Math.pow(a, 2) * k3 + Math.pow(a, 3) * k4 + Math.pow(a, 4) * k5 + b * k6 + Math.pow(b, 2) * k7 + Math.pow(b, 3) * k8 + Math.pow(b, 4) * k9)) / Math.pow((k12 + a * k13 + Math.pow(a, 2) * k14 + Math.pow(a, 3) * k15 + Math.pow(a, 4) * k16 + b * k17 + Math.pow(b, 2) * k18 + Math.pow(b, 3) * k19 + Math.pow(b, 4) * k20 + a * b * k21 + Math.pow(a, 2) * Math.pow(b, 2) * k22), 2));
			derivW[19-1] = -((Math.pow(b, 3) * (k1 + a * b * k10 + Math.pow(a, 2) * Math.pow(b, 2) * k11 + a * k2 + Math.pow(a, 2) * k3 + Math.pow(a, 3) * k4 + Math.pow(a, 4) * k5 + b * k6 + Math.pow(b, 2) * k7 + Math.pow(b, 3) * k8 + Math.pow(b, 4) * k9)) / Math.pow((k12 + a * k13 + Math.pow(a, 2) * k14 + Math.pow(a, 3) * k15 + Math.pow(a, 4) * k16 + b * k17 + Math.pow(b, 2) * k18 + Math.pow(b, 3) * k19 + Math.pow(b, 4) * k20 + a * b * k21 + Math.pow(a, 2) * Math.pow(b, 2) * k22), 2));
			derivW[20-1] = -((Math.pow(b, 4) * (k1 + a * b * k10 + Math.pow(a, 2) * Math.pow(b, 2) * k11 + a * k2 + Math.pow(a, 2) * k3 + Math.pow(a, 3) * k4 + Math.pow(a, 4) * k5 + b * k6 + Math.pow(b, 2) * k7 + Math.pow(b, 3) * k8 + Math.pow(b, 4) * k9)) / Math.pow((k12 + a * k13 + Math.pow(a, 2) * k14 + Math.pow(a, 3) * k15 + Math.pow(a, 4) * k16 + b * k17 + Math.pow(b, 2) * k18 + Math.pow(b, 3) * k19 + Math.pow(b, 4) * k20 + a * b * k21 + Math.pow(a, 2) * Math.pow(b, 2) * k22), 2));
			derivW[21-1] = -((a * b * (k1 + a * b * k10 + Math.pow(a, 2) * Math.pow(b, 2) * k11 + a * k2 + Math.pow(a, 2) * k3 + Math.pow(a, 3) * k4 + Math.pow(a, 4) * k5 + b * k6 + Math.pow(b, 2) * k7 + Math.pow(b, 3) * k8 + Math.pow(b, 4) * k9)) / Math.pow((k12 + a * k13 + Math.pow(a, 2) * k14 + Math.pow(a, 3) * k15 + Math.pow(a, 4) * k16 + b * k17 + Math.pow(b, 2) * k18 + Math.pow(b, 3) * k19 + Math.pow(b, 4) * k20 + a * b * k21 + Math.pow(a, 2) * Math.pow(b, 2) * k22), 2));
			derivW[22-1] = -((Math.pow(a, 2) * Math.pow(b, 2) * (k1 + a * b * k10 + Math.pow(a, 2) * Math.pow(b, 2) * k11 + a * k2 + Math.pow(a, 2) * k3 + Math.pow(a, 3) * k4 + Math.pow(a, 4) * k5 + b * k6 + Math.pow(b, 2) * k7 + Math.pow(b, 3) * k8 + Math.pow(b, 4) * k9)) / Math.pow((k12 + a * k13 + Math.pow(a, 2) * k14 + Math.pow(a, 3) * k15 + Math.pow(a, 4) * k16 + b * k17 + Math.pow(b, 2) * k18 + Math.pow(b, 3) * k19 + Math.pow(b, 4) * k20 + a * b * k21 + Math.pow(a, 2) * Math.pow(b, 2) * k22), 2));

			n = covKs.length; //size(covKs,1);
			for (i = 0; i < n; i++) {//DO i = 1, n
				varKS[i] = covKs[i][i]; //! Variance
			}//END DO

	//sumVarKs = sum(varKS*(Math.pow(derivW, 2)));
			sumVarKs = IntStream.range(0, Math.min(varKS.length, derivW.length)).mapToDouble(id -> varKS[id] * (Math.pow(derivW[id], 2))).sum();

			sumcovarKs = 0.0;

			for (i = 0; i < n - 1; i++) {//DO i = 1, n-1
				for (j = i + 1; j < n; j++) {//DO j = i+1, n
					sumcovarKs = sumcovarKs + 2.0 * covKs[i][j] * derivW[i] * derivW[j];
				}//END DO //!j
			}//END DO //!i

			returnSigma.varW = sumcovarKs + sumVarKs;
			StDevW = Math.sqrt(returnSigma.varW) * 2.0;

			//return returnSigma;
		}// END SUBROUTINE SIGMA

		public static String toDebug() {
			//varAlpha, varW
			return "returnSIGMA varAlpha "+varAlpha+" varW "+varW;
		}
	}
}//END MODULE EVALW
