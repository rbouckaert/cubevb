package cubevb.transform;

import java.util.List;

public class TransformUtils {

	static double cor(List<double[]> trace, double mean1, double mean2, int i1, int i2) {
		if (trace.size() == 1) {
			return 0;
		}
		double sum = 0;
		int burnin = 0;//x1.length/4;
		for (int i = burnin; i < trace.size(); i++) {
			double [] x = trace.get(i);
			sum += (x[i1] - mean1) * (x[i2] - mean2);
		}
		return sum / (trace.size() - burnin - 1.0);
	}

	static double mean(List<double[]> trace, int i) {
		double sum = 0;
		for (double [] d : trace) {
			sum += d[i];
		}
		return sum/trace.size();
	}

	static double var(List<double[]> trace, int i, double mean) {
		double sum = 0;
		for (double [] d : trace) {
			sum += (d[i] - mean)* (d[i] - mean);
		}
		return sum/(trace.size() - 1);
	}

	static public void initMeanStdevCovar(List<double[]> trace, int n, double [] mean, double [] stdev, double [][] covar0) {
		for (int i = 0; i < mean.length; i++) {
			mean[i] = mean(trace, i);
			double var = var(trace, i, mean[i]);
			
			double sigma = var != 0 ?
					Math.sqrt(var) : 
					1e-10;// rather have Double.MAX_VALUE, but that equals Double.POSITIVE_INFINITY
			stdev[i] = sigma; // Transform.softplusInverse(sigma);
		}
		for (int i = 0; i < n; i++) {
			for (int j = i + 1; j < n; j++) {
				covar0[i][j] = cor(trace, mean[i], mean[j], i, j);
				covar0[j][i] = covar0[i][j];
			}
			covar0[i][i] = cor(trace, mean[i], mean[i], i, i);//stdev[i] * stdev[i];
		}
	}

}
