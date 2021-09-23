package net.priveyes.PhaseDiag.CHNOZ;

import net.priveyes.PhaseDiag.num_rep.VecMath;

public class IAPWScommon {
	// normal boiling point
	// triple point
	// critical point constants

	public static class Temperature {
		static double[] T;
		public static double defaut = 298.15; // K 25°C
		static double triple = 273.1600; // K
		static double boiling = 373.124;
		public static double critical = 647.096; // K
	}

	public static class Rho {
		static double[] rho;
		static double defaut = 1000.0;
		//static Phase triple;
		//static Phase boiling;
		public static double critical = 322.0; // kg m-3

		static class Phase {
			static double[] liquid = new double[100];
			static double[] vapor = new double[100];
		}
		static class triple{static double liquid; static double vapor;}

		static class boiling{static double liquid;static double vapor;}

		static {
			Rho.triple.liquid = 999.793; Rho.triple.vapor = 0.00485458;
			Rho.boiling.liquid = 958.367; Rho.boiling.vapor = 0.597657;
		}
		public static double sat;
	}

	static class Pressure {
		public static double[] init;
		public static double defaut;
		static double triple = 611.657; // Pa
		static double boiling = 0.101325;
		public static double critical = 22.064; // MPa
		public static double[] MPa;
		public static double sigma;
		static double[] convert(double[] toConvert, String mPa) {
			// 1 bar = 0.1 MPa
			return VecMath.vxs(toConvert, 0.1);
		}
	}

	static class G {
		public static double[] liquid;
		public static double[] vapor;
	}

	public static class idealgas {
		static double phi;

		static class Phi {
			static double tau;
			static double delta;

			static class Tau {
				static double tau;
			}

			static class Delta {
				static double tau;
				static double delta;
			}
		}

		static class dDelta {
			static double[] ddelta = new double[56];

			static class Bi {
				static double[] dtau = new double[56];
				static double[] ddelta = new double[56];
			}
		}

		static class d2Delta {
			static double[] ddelta2 = new double[56];

			static class Bi {
				static class Ddelta {
					static double[] dtau = new double[56];
				}

				static double[] dtau2 = new double[56];
				static double[] ddelta2 = new double[56];
			}
		}

		static class dPsi {
			static double[] Delta = new double[56];
			static double[] dtau = new double[56];
			static double[] ddelta = new double[56];

			static class Ddelta {
				static double dtau;
			}
		}

		static class d2Psi {
			static class Ddelta {
				static double[] dtau = new double[56];
			}

			static double[] dtau2 = new double[56];
			static double[] ddelta2 = new double[56];
		}
	}

	public static class residual {
		public static double phi;

		public static class Phi {
			public static double tau;
			public static double delta;

			public static class Tau {
				public static double tau;
			}

			public static class Delta {
				public static double tau;
				public static double delta;
			}
		}
	}
}
