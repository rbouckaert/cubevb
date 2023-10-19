package cubevb;


import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import javax.xml.parsers.ParserConfigurationException;

import org.xml.sax.SAXException;

import beast.base.core.BEASTInterface;
import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.inference.CompoundDistribution;
import beast.base.inference.Distribution;
import beast.base.inference.Logger;
import beast.base.inference.MCMC;
import beast.base.inference.Operator;
import beast.base.inference.OperatorSchedule;
import beast.base.inference.Runnable;
import beast.base.inference.State;
import beast.base.inference.StateNode;
import beast.base.inference.StateNodeInitialiser;
import beast.base.core.Log;
import beast.base.util.Randomizer;
import beast.base.inference.parameter.RealParameter;
import beast.base.math.matrixalgebra.CholeskyDecomposition;
import beast.base.math.matrixalgebra.IllegalDimension;
import beast.base.evolution.branchratemodel.BranchRateModel;
import cubevb.transform.CubeTransform;
import cubevb.transform.StickBreaking;
import cubevb.transform.Transform;

@Description("Variational inference by initialising through simulated annealing, "
		+ "finding variances through MCMC")
public class VIBySAandMCMC extends Runnable {

    final public Input<State> stateInput = new Input<>("state", "elements of the state space");

    final public Input<List<Transform>> transformsInput = new Input<>("transform", "transformations on elements of the state space", new ArrayList<>());

    final public Input<List<StateNodeInitialiser>> initialisersInput =
            new Input<>("init", "one or more state node initilisers used for determining " +
                    "the start state of the chain", new ArrayList<>());

    final public Input<Distribution> posteriorInput =
            new Input<>("distribution", "probability distribution to sample over (e.g. a posterior)", Input.Validate.REQUIRED);

    final public Input<List<Logger>> loggersInput =
            new Input<>("logger", "loggers for reporting progress of MCMC chain", new ArrayList<>(), Input.Validate.REQUIRED);

    
    final public Input<Integer> sampleCountInput = new Input<>("samples", "number of samples to generate from posterior", 1000);

    final public Input<Double> alphaInput = new Input<>("learningRate", "learning rate for SGD", 1e-8);
    
    final public Input<Integer> chainLengthInput = new Input<>("chainLength", "expected number of MCMC samples per parameter for initialisation", 64);
    
	final public Input<Boolean> useCorrelationsInput = new Input<>("useCorrelations", "if true, use correlations between " +
			"consecutive parameter entries", true);


	final public Input<Boolean> useMCMCInTransformedSpaceInput = new Input<>("useMCMCInTransformedSpace", "if true, use MCMC in transformed space "
			+ "(operators and operatorschedule will be ignored if specified) "
			+ "if false MCMC is used in state space (operators must be specified and operatorschedule can be specified)", false);

	final public Input<List<Operator>> operatorsInput =
            new Input<>("operator", "operator for generating proposals in MCMC state space",
                    new ArrayList<>());	//, Input.Validate.REQUIRED);
    final public Input<OperatorSchedule> operatorScheduleInput = new Input<>("operatorschedule", "specify operator selection and optimisation schedule");

	
    /** mean and standard deviation for multivariate distribution before transform **/
	double[] mean;
	double[] stdev;
	
	/** random draw **/
	double [] s;

	/** gradient wrt mean and stdev **/
	double[] dMean;
	double[] dStdev;
	// L = Cholesky decomposition of covariance matrix
	// used only when useCorrelations = true
	double [][] L;
	
	/** sum of gradients squared for AdaGrad weight updating **/
	double [] G;

	
	double jacobian;

	List<Transform> transforms;
    
	State state;
	Distribution posterior;
	List<Logger> loggers;
    
	@Override
	public void initAndValidate() {
		transforms = transformsInput.get();
		int n = 0;
		for (Transform t : transforms) {
			n += t.getDimension();
		}
		mean = new double[n];
		stdev = new double[n];
		Arrays.fill(stdev, -1);
		s = new double[n];
		Arrays.fill(s, 1e-3);
		dMean = new double[n];
		dStdev = new double[n];
		G = new double[2*n];
		
		posterior = posteriorInput.get();
		
		state = stateInput.get();
		state.initialise();
		state.setPosterior(posterior);
		
		loggers = loggersInput.get();

		
		// make sure all parameters in transforms are in the state
		for (Transform t : transforms) {
			StateNode transformNode = t.getStateNode();
			boolean found = false;
			for (StateNode stateNode : state.stateNodeInput.get()) {
				if (transformNode == stateNode) {
					found = true;
					break;
				}
			}
			if (!found) {
				throw new IllegalArgumentException("Found a stateNode " + transformNode.getID() + " in transform " + t.getID() + " that is not in the state");
			}
		}
		
		
		// make sure all parameters in state are associated with a transforms
		for (StateNode stateNode : state.stateNodeInput.get()) {
			boolean found = false;
			for (Transform t : transforms) {
				StateNode transformNode = t.getStateNode();
				if (transformNode == stateNode) {
					found = true;
					break;
				}
			}
			if (!found) {
				Log.warning("Found a stateNode " + stateNode.getID() + " in the state that has no transform");
			}
		}
	}

	@Override
	public void run() throws Exception {
        long startTime = System.currentTimeMillis();
		
		if (restoreFromFile) {
            try {
                state.setStateFileName(stateFileName);
				state.restoreFromFile();
			} catch (SAXException | IOException | ParserConfigurationException e) {
				e.printStackTrace();
				throw new RuntimeException(e);
			}
			
		} else {
			for (final StateNodeInitialiser initialiser : initialisersInput.get()) {
				initialiser.initStateNodes();
			}
		}

		initState();

		for (Logger logger : loggers) {
			logger.init();
		}
		
        long endTime = System.currentTimeMillis();
        Log.info.println("Total initialisation time: " + (endTime - startTime) / 1000.0 + " seconds");

		// epsilon prevents learning rate to collapse in AdaGrad
		double epsilon = 1e-3;
        // learning rate
		double alpha = alphaInput.get();
		// RMSprop weight
		double beta = 0.9;
		
		// change in mean or stdev for calculating derivatives
		double DELTA = 1e-6;
		

		


        endTime = System.currentTimeMillis();
        Log.info.println("Total optimisation time: " + (endTime - startTime) / 1000.0 + " seconds");
        
        sample();

		for (Logger logger : loggers) {
			logger.close();
		}

		endTime = System.currentTimeMillis();
        Log.info.println("Total runtime: " + (endTime - startTime) / 1000.0 + " seconds");
	}
	


	public void initState() {
		int k = 0;
		int clockIndex = -1;
		RealParameter clockParameter = null;
		Transform clockTransform = null;
		for (Transform t : transforms) {
			if (t.getDimension() == 1 && !(t instanceof StickBreaking) && !(t instanceof CubeTransform)) {
				// is this the mean clock rate parameter?
				if (t.getParameter().getID().startsWith("precision")) {
					// TODO: replace the above hack
					clockIndex = k;
					clockParameter = t.getParameter();
					clockTransform = t;
				} else 
				for (BEASTInterface o : t.getParameter().getOutputs()) {
					if (o instanceof BranchRateModel.Base) {
						if (((BranchRateModel.Base)o).meanRateInput.get() == t.parameterInput.get()) {
							clockIndex = k;
							clockParameter = t.getParameter();
							clockTransform = t;
						}
					}
				}
				
				// get transformed value
				double d = t.forward(t.getTransformNode().getDoubleValues()[0]);
				mean[k++] = d;
			} else {
				double [] d = t.forward();
				for (int i = 0; i < d.length; i++) {
					mean[k++] = d[i];
				}
			}
		}

		state.initialise();
		
		// only mean
		Arrays.fill(s, 0);

		
		
		int chainLength = chainLengthInput.get();
		if (useMCMCInTransformedSpaceInput.get()) {
			initStdevByMCMC(chainLength);
		} else {
			initStdevByMCMC1(chainLength);
		}
	}
	
	private void initStdevByMCMC1(int chainLength) {
		double [][] trace = new double[mean.length][chainLength];

		MCMC mcmc = new MCMC() {
			protected void callUserFunction(long sample) {
				if (sample > 0 && sample % mean.length == 0) {
					int x = (int)(sample/mean.length-1);
					int k = 0;
					for (Transform t : transforms) {
						if (t.getDimension() == 1 && !(t instanceof StickBreaking)&& !(t instanceof CubeTransform)) {
							// get transformed value
							double d = t.forward(t.getTransformNode().getDoubleValues()[0]);
							trace[k++][x] = d;
						} else {
							double [] d = t.forward();
							for (int i = 0; i < d.length; i++) {
								trace[k++][x] = d[i];
							}
						}
					}
				}
//				if (sample % chainLength == 0) {
//					if (sample % (10*chainLength) != 0) {
//						Log.warning.print('.');
//					} else {
//						Log.warning.print('|');
//					}
//				}
			};
		};
		
		Logger screenLogger = null;
		for (Logger logger : loggers) {
			if (logger.getFileName() == null || logger.getFileName().trim().length() == 0) {
				screenLogger = logger;
			}
		}
		screenLogger.everyInput.setValue(mean.length, screenLogger);
		OperatorSchedule operatorschedule = operatorScheduleInput.get();
		if (operatorschedule == null) {
			operatorschedule = new CubeOperatorSchedule();
			operatorschedule.initByName("autoOptimize", true, "weightSchedule", "AUTO");
		}
		mcmc.initByName("state", state,
				"distribution", posterior,
				"operator", operatorsInput.get(),
				"logger", screenLogger,
				"chainLength", (long)(mean.length * chainLength),
				"operatorschedule", operatorschedule);
		
		try {
			mcmc.run();
		} catch (IOException | SAXException | ParserConfigurationException e1) {
			e1.printStackTrace();
		}
		Log.warning.print('\n');
		

		// estimate values of mean and stdev based on trace
		for (int i = 0; i < mean.length; i++) {
			mean[i] = mean(trace[i]);
			double var = var(trace[i], mean[i]);
			
			double sigma = var != 0 ?
					Math.sqrt(var) : 
					1e-10;// rather have Double.MAX_VALUE, but that equals Double.POSITIVE_INFINITY
			stdev[i] = Transform.softplusInverse(sigma);
		}
		
		
		if (useCorrelationsInput.get()) {
			try {
				L = getL(trace);
			} catch (IllegalDimension e) {
				e.printStackTrace();
			}
		}
		
		// for debugging, make sure logP3 == logP0:
		double logP3 = calcLogP();
		System.err.println(logP3 + " " + Arrays.toString(stdev));
		//System.err.println(logP0 + " =?= " + logP3);
		//System.err.println(Arrays.toString(stdev).replaceAll(",","\n"));
	}
	

	private double var(double[] trace, double mean) {
		double sum = 0;
		for (double d : trace) {
			sum += (d - mean)* (d - mean);
		}
		return sum/(trace.length - 1);
	}

	private double mean(double[] trace) {
		double sum = 0;
		for (double d : trace) {
			sum += d;
		}
		return sum/trace.length;
	}

	private void initStdevByMCMC(int chainLength) {
		double logP0 = calcLogP();
		
		double [] values = new double[mean.length * 2];
//		int [] sampleCount = new int[mean.length];
//		double [] sum2 = new double[stdev.length];
//		double [] sum = new double[mean.length];
		double windowSize = 0.05;

		
		double [][] trace = new double[mean.length][chainLength];
		
		System.arraycopy(mean, 0, values, 0, mean.length);
		double oldLogP = logP0;
		double oldJacobian = jacobian;
		for (long j = 1; j <= chainLength * mean.length; j++) {
			int i = Randomizer.nextInt(mean.length);
			double s = Randomizer.nextGaussian();
			double oldValue = values[i];			
			double newValue = oldValue + s * windowSize;
	        state.store(j);
	        double newJacobian = setParams(i, oldValue, newValue);
			
	        if (newJacobian != Double.NEGATIVE_INFINITY) {
	            state.storeCalculationNodes();
	            state.checkCalculationNodesDirtiness();
	        	double newLogP = posterior.calculateLogP();
				double logAlpha = newLogP - oldLogP - oldJacobian + newJacobian;// + Math.log(values[i]/(values[i] - s * windowSize));
				if (logAlpha > 0 || Randomizer.nextDouble() < Math.exp(logAlpha)) {
					// accept
					oldLogP = newLogP;
					oldJacobian = newJacobian;
					values[i] = newValue; 
	                state.acceptCalculationNodes();
					jacobian = newJacobian;
				} else {
					// reject
	                state.restore();
	                state.restoreCalculationNodes();
				}
	        } else {
	        	// direct reject
	            state.restore();	        	
				jacobian = oldJacobian;
	        }
            state.setEverythingDirty(false);
//			sum2[i] += values[i] * values[i];
//			sum[i] += values[i];
//			sampleCount[i]++;

			
			if (j % mean.length == 0) {
				int x = (int) (j/mean.length-1);
				for (int k = 0; k < mean.length; k++) {
					trace[k][x] = values[k];
				}
//				for (Logger log : loggers) {
//					log.log(log.everyInput.get() * (j / mean.length - 1));
//				}
			}
			if (j % chainLength == 0) {
				if (j % (10*chainLength) != 0) {
					Log.warning.print('.');
				} else {
					Log.warning.print('|');
				}
			}
		}
		Log.warning.print('\n');

		// estimate values of mean and stdev based on trace
		for (int i = 0; i < mean.length; i++) {
			mean[i] = mean(trace[i]);
			double var = var(trace[i], mean[i]);
			
			double sigma = var != 0 ?
					Math.sqrt(var) : 
					1e-10;// rather have Double.MAX_VALUE, but that equals Double.POSITIVE_INFINITY
			stdev[i] = Transform.softplusInverse(sigma);
		}
		
		
		if (useCorrelationsInput.get()) {
			try {
				L = getL(trace);
			} catch (IllegalDimension e) {
				e.printStackTrace();
			}
		}

		
		// for debugging, make sure logP3 == logP0:
		double logP3 = calcLogP();
		System.err.println(logP3 + " " + Arrays.toString(stdev));
		//System.err.println(logP0 + " =?= " + logP3);
		//System.err.println(Arrays.toString(stdev).replaceAll(",","\n"));
	}

	/** get Cholesky decomposition cleaned up for NaN and Infinities **/
	private double[][] getL(double[][] trace) throws IllegalDimension {
		int n = mean.length;
		double[][] covar0 = new double[n][n];
		for (int i = 0; i < n; i++) {
			for (int j = i + 1; j < n; j++) {
				covar0[i][j] = cor(mean[i], mean[j], trace[i], trace[j]);
				covar0[j][i] = covar0[i][j];
			}
			covar0[i][i] = cor(mean[i], mean[i], trace[i], trace[i]);//stdev[i] * stdev[i];
		}

//		int nonzero = 0, zero = 0;
//		for (int i = 0; i < covar0.length; i++) {
//			for (int j = 0; j < i; j++) {
//				if (covar0[i][j] < 0.01) {
//					covar0[i][j] = 0;
//					zero++;
//				} else {
//					nonzero++;
//				}
//			}
//			for (int j = i+1; j < covar0[i].length; j++) {
//				if (covar0[i][j] < 0.01) {
//					covar0[i][j] = 0;
//					zero++;
//				} else {
//					nonzero++;
//				}
//			}
//		}
		System.out.println("mean = " + Arrays.toString(mean));
		
		DecimalFormat df = new DecimalFormat("#.####");
		for (int i = 0; i < covar0.length; i++) {
			System.out.print("[");
			for (int j = 0; j < covar0.length; j++) {
				String str = df.format(covar0[i][j]/(Transform.softplus(stdev[i])*Transform.softplus(stdev[j])));
				System.out.print("       ".substring(str.length()) + str);
			}
			System.out.print(" ]\n");
		}
//		System.out.println(zero + " zeros, " + nonzero +" non zeros");

		
		CholeskyDecomposition chol = new CholeskyDecomposition(covar0);
		double [][] L = chol.getL();
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				if (Double.isNaN(L[i][j])||Double.isInfinite(L[i][j])) {
					L[i][j] = 0;
				}
			}
		}
		
//		DenseCholesky chol = DenseCholesky.factorize(new DenseMatrix(covar0));
//		UpperTriangDenseMatrix L0 = chol.getU();
//		//CholeskyDecomposition chol = new CholeskyDecomposition(covar0);
//		double [][] L = new double[n][n];
//		int nanCount = 0, infCount = 0;
//		for (int i = 0; i < n; i++) {
//			for (int j = 0; j < n; j++) {
//				double d = L0.get(i,  j);
//				L[i][j] = d;
//				if (Double.isNaN(L[i][j])||Double.isInfinite(L[i][j])) {
//					L[i][j] = 0;
//					if (Double.isNaN(d)) {
//						nanCount++;
//					} else {
//						infCount++;
//					}
//				}
//			}
//		}		
//		System.err.println(nanCount + " NaNs " + infCount + " Inftys");
		
		return L;
	}
	
	private double cor(double mean1, double mean2, double[] x1, double[] x2) {
		if (x1.length == 1) {
			return 0;
		}
		double sum = 0;
		int burnin = 0;//x1.length/4;
		for (int i = burnin; i < x1.length; i++) {
			sum += (x1[i] - mean1) * (x2[i] - mean2);
		}
		return sum / (x1.length - burnin - 1.0);
	}

		
	private void initStdev(int chainLength) {
		double logP0 = calcLogP();
		
		double [] values = new double[mean.length * 2];
		double [] stdevEstimate = new double[stdev.length];
		double [] meanEstimate = new double[mean.length];
		double windowSize = 0.05;
		for (int i = 0; i < mean.length; i++) {
			// do short MCMC to determine stdev
			System.arraycopy(mean, 0, values, 0, mean.length);
			double oldLogP = logP0;
			double var = 0;
			double m0 = mean[i];
			double sum = 0;
			int acceptCount = 0;
			for (int j = 0; j < chainLength; j++) {
				double s = Randomizer.nextGaussian();
				double m = values[i];
				values[i] += s * windowSize;
				double newLogP = calcLogP(values);
				if (newLogP > oldLogP || Randomizer.nextDouble() < Math.exp(newLogP - oldLogP)) {
					// accept
					oldLogP = newLogP;
					acceptCount++;
				} else {
					// reject
					values[i] = m;
				}
				var += (values[i] - m0) * (values[i] - m0);
				sum += values[i];
			}
			if (acceptCount == 0) {
				int h = 3;
				h--;
			}
			// restore original value of mean,
			// this got changed with calls to calcLogP(values);
			mean[i] = m0;
			meanEstimate[i] = sum / chainLength;
			var /= chainLength;
			double sigma = var != 0 ?
					Math.sqrt(var) : 
					1e-10;// rather have Double.MAX_VALUE, but that equals Double.POSITIVE_INFINITY
			stdevEstimate[i] = Transform.softplusInverse(sigma);
			Log.warning.print(acceptCount + " " + var + " " + stdevEstimate[i] + "\n");

			
			if (i % 10 != 0) {
				Log.warning.print('.');
			} else {
				Log.warning.print('|');				
			}
		}
		System.arraycopy(stdevEstimate, 0, stdev, 0, stdev.length);
		System.arraycopy(meanEstimate,  0, mean,  0, mean.length);

		
		// for debugging, make sure logP3 == logP0:
		double logP3 = calcLogP();
		System.err.println(logP3 + " " + Arrays.toString(stdev));
		//System.err.println(logP0 + " =?= " + logP3);
		//System.err.println(Arrays.toString(stdev).replaceAll(",","\n"));
	}



    protected void reportLogLikelihoods(final Distribution distr, final String tabString) {
        final double full =  distr.getCurrentLogP(), last = distr.getStoredLogP();
        final String changed = full == last ? "" : "  **";
        Log.info.println(tabString + "P(" + distr.getID() + ") = " + full + " (was " + last + ")" + changed);
        if (distr instanceof CompoundDistribution) {
            for (final Distribution distr2 : ((CompoundDistribution) distr).pDistributions.get()) {
                reportLogLikelihoods(distr2, tabString + "\t");
            }
        }
    }

	public void sample() {
		long base = 0; // iterationsInput.get();
		long delta = loggers.get(0).everyInput.get();
		for (int i = 0; i < sampleCountInput.get(); i++) {
			// draw sample from variational distribution
			for (int k = 0; k < s.length; k++) {			
				s[k] = Randomizer.nextGaussian();
			}
			
			if (useCorrelationsInput.get()) {
				double [] s0 = new double[s.length];
				for (int j = 0; j < s.length; j++) {
					s0[j] = mult(L, j, s) / Transform.softplus(stdev[j]);
				}
				System.arraycopy(s0,  0,  s, 0, s.length);
			}
			setParams();

	        state.store(-1);
	        state.setEverythingDirty(true);
	        state.storeCalculationNodes();
	        state.checkCalculationNodesDirtiness();
        	double logP = posterior.calculateLogP();
			state.acceptCalculationNodes();
        	if (Double.isInfinite(logP)) {
        		i--;
        		// MAPInitialiser.reportLogLikelihoods(posterior, "");
        	} else {
        		log(base, false);
    			base += delta;
        	}
		}
		
	}
	
	
	
	private double mult(double[][] L, final int index, double[] sample) {
		double sum = 0;
		for (int i = 0; i < L.length; i++) {
			sum += L[index][i] * sample[i];
		}
		return sum;
	}

	
	protected void log(long base, boolean useScreenLogOnly) {
		for (Logger logger : loggers) {
			if (useScreenLogOnly) {
				if (logger.fileNameInput.get() == null) {
					logger.log(base);
				}
			} else if (logger.fileNameInput.get() != null) {
				logger.log(base);
			}
		}
	}

	/** return true if all transforms are finite **/
	protected boolean setParams() {
		int k = 0;
		jacobian = 0;
		boolean allFinite = true;
		for (Transform t : transforms) {
			if (t.getDimension() == 1 && !(t instanceof StickBreaking) && !(t instanceof CubeTransform)) {
				double v = Transform.softplus(stdev[k]) * s[k] + mean[k];
				t.backward(v);
				if (Double.isInfinite(v)) {
					allFinite = false;
				}
				jacobian += t.jacobian_det(v);
				k++;
			} else {
				int d = t.getDimension();
				double [] v = new double[d];
				for (int i = 0; i < d; i++) {
					v[i] = Transform.softplus(stdev[k]) * s[k] + mean[k];
					if (Double.isInfinite(v[i])) {
						allFinite = false;
					}
					k++;
				}
				t.backward(v);
				jacobian += t.jacobian_det(v);
			}
		}
		return allFinite;
	}

	/** return true if all transforms are finite **/
	protected double setParams(int dim, double oldValue, double newValue) {
		int k = 0;
		double jacobian = this.jacobian;
		for (Transform t : transforms) {
			if (t.getDimension() == 1 && !(t instanceof StickBreaking) && !(t instanceof CubeTransform)) {
				if (dim == k) {
					t.backward(0, oldValue, newValue);
					jacobian += -t.jacobian_det(oldValue) + t.jacobian_det(newValue);
					return jacobian;
				}
				k++;
			} else {
				int d = t.getDimension();
				if (k <= dim && dim < k + d) {
					double logHR = t.backward(dim - k, oldValue, newValue);
					if (Double.isInfinite(logHR)) {
						jacobian = Double.NEGATIVE_INFINITY;
						return Double.NEGATIVE_INFINITY;
					}
					jacobian += -t.jacobian_det(dim-k, oldValue) + t.jacobian_det(dim-k, newValue);
					return jacobian;
				}
				k += d;
			}
		}
		throw new IllegalArgumentException(getClass().getName() + ".setParams() programmer error: dimension too large");
	}
	
	private double calcLogQ() {
    	double logQ = 0;
    	for (int i = 0; i < s.length; i++) {
    		logQ += Math.log(Transform.softplus(stdev[i]));
    	}
    	return jacobian - logQ;
	}

	public int getDimension() {
		return s.length * 2;
	}
	
	
	private double calcLogP() {
		if (!setParams()) {
			return Double.NEGATIVE_INFINITY;
		}
		
        state.store(-1);
        state.setEverythingDirty(true);
        state.storeCalculationNodes();
        state.checkCalculationNodesDirtiness();
    	double logP = posterior.calculateLogP();
    	if (Double.isInfinite(logP) || Double.isNaN(logP)) {
    		reportLogLikelihoods(posterior, "");
        	if (Double.isNaN(logP)) {
	    		Distribution d = ((CompoundDistribution) posterior).pDistributions.get().get(1);
	    		d = ((CompoundDistribution) d).pDistributions.get().get(0);
	    		d.calculateLogP();
        	}
        	logP = -1e200;
    	}
		state.acceptCalculationNodes();
		return logP;
	}
	
	
	public double calcLogP(double [] x) {
		System.arraycopy(x, 0, mean, 0, mean.length);
		System.arraycopy(x, stdev.length, stdev, 0, stdev.length);

		Arrays.fill(s, 0);		
		setParams();
		
        //state.store(-1);
        //state.setEverythingDirty(true);
        state.storeCalculationNodes();
        state.checkCalculationNodesDirtiness();
    	double logP = posterior.calculateLogP();
		state.acceptCalculationNodes();
		
		return logP;
	}
	
	public double calcELBO(double [] x) {
		double ELBO = 0;
//		Arrays.fill(s, -1);		
//		ELBO += calcELBO0(x);
		Arrays.fill(s, 0);		
		ELBO += calcELBO0(x);
//		Arrays.fill(s, 0.5);		
//		ELBO += calcELBO0(x);
		return ELBO;
	}
		
	private double calcELBO0(double [] x) {
		System.arraycopy(x, 0, mean, 0, mean.length);
		System.arraycopy(x, stdev.length, stdev, 0, stdev.length);
		setParams();
		
        state.store(-1);
        state.setEverythingDirty(true);
        state.storeCalculationNodes();
        state.checkCalculationNodesDirtiness();
    	double logP = posterior.calculateLogP();
		state.acceptCalculationNodes();
		
		return logP + jacobian -calcLogQ(); 
	}
}
