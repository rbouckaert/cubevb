package cubevb.tool;




import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintStream;
import java.util.*;

import org.json.JSONArray;
import org.json.JSONException;
import org.json.JSONObject;

import beastfx.app.inputeditor.BeautiDoc;
import beastfx.app.tools.Application;
import beastfx.app.treeannotator.TreeAnnotator;
import beastfx.app.treeannotator.TreeAnnotator.FastTreeSet;
import beastfx.app.util.OutFile;
import beastfx.app.util.TreeFile;
import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.core.Log;
import beast.base.evolution.alignment.TaxonSet;
import beast.base.evolution.alignment.Alignment;
import beast.base.evolution.alignment.Sequence;
import beast.base.evolution.alignment.Taxon;
import beast.base.evolution.distance.Distance;
import beast.base.evolution.tree.CladeSet;
import beast.base.evolution.tree.ClusterTree;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.evolution.tree.MRCAPrior;
import beast.base.math.matrixalgebra.CholeskyDecomposition;
import beast.base.math.matrixalgebra.IllegalDimension;
// import beast.base.evolution.tree.ClusterTree;
import beast.base.util.Randomizer;

@Description("Convert posterior tree set to matrix with correlation structure and sample tree set from matrix + structure")
public class Trees2Cube extends beast.base.inference.Runnable {

	final public Input<TreeFile> treesInput = new Input<>("trees", "beast.trees on which this operation is performed",
			new TreeFile("[[none]]"));
	final public Input<OutFile> outputInput = new Input<>("out", "output file. Fails if not specified", Validate.REQUIRED);
	final public Input<Integer> sampleCountInput = new Input<>("n", "number of trees to sample", 1000);
	final public Input<Integer> burnInPercentageInput = new Input<>("burnin",
			"percentage of trees to used as burn-in (and will be ignored). Set >= 100 if only a cube is required", 10);

	final public Input<File> orderInput = new Input<>("order",
			"file with taxon ordering, can be space, tab, comma or new line separated", Validate.REQUIRED);

	final public Input<Boolean> useCorrelationsInput = new Input<>("useCorrelations", "if true, use correlations between " +
			"matrix entries", true);

	final public Input<Boolean> useLogTransformInput = new Input<>("useLogTransform", "if true, multivariate log-normal "
			+ "is used instead of normal", true);
	final public Input<Boolean> useMonoConstraintsInput = new Input<>("useMonoConstraints", "if true, monophyletic "
			+ "constraints are applied to any clade with 100% support in input tree set", true);
	final public Input<Long> seedInput = new Input<>("seed", "random number generator seed");
	
	
	private String[] taxa;
	private boolean useLogTransform = true, useMonoConstraints = true;
	
	
	private class MultivariateNormal {
		int [][] pairs;
		double [] means;
		double [][] covar0;
		double [] tipHeights;
		String [] taxaNames;
		
		List<MRCAPrior> monophylyConstraints;
		
		MultivariateNormal(int [][] pairs,
				double [] means,
				double [][] covar0,
				double [] tipHeights,
				String [] taxaNames,
				List<MRCAPrior> monophylyConstraints) {
			this.pairs = pairs;
			this.means = means;
			this.covar0 = covar0;
			this.tipHeights = tipHeights;
			this.taxaNames = taxaNames;
			this.monophylyConstraints = monophylyConstraints;
		}
		
		MultivariateNormal(File matrixFile) throws IOException, JSONException {
			fromJSON(matrixFile);
		}
		
		void toJSON(PrintStream out) {
			out.print("{");
			int n = taxaNames.length;
			int np = pairs.length;
			// first the np means
			out.println("\"n\":" + n + ",");
			
			out.print("\"taxanames\":[");
			for (int i = 0; i < n; i++) {
				out.print("\"" + taxaNames[i] + "\"" + (i < n-1 ?",":""));
			}
			out.println("],");
			
			out.print("\"tipheights\":[");
			for (int i = 0; i < n; i++) {
				out.print(tipHeights[i] + (i < n-1 ?",":""));
			}
			out.println("],");

			out.print("\"monophylyConstraints\":[");
			int k = monophylyConstraints.size();
			for (int i = 0; i < k; i++) {
				TaxonSet taxonset = monophylyConstraints.get(i).taxonsetInput.get();
				out.print("[");
				int l = taxonset.getTaxonCount();
				for (int j = 0; j < l; j++) {
					out.print(taxonset.getTaxonId(j) + (j < l-1 ?",":""));
				}
				out.println("]" + (i < k-1 ?",":""));
			}
			out.println("],");
			
			out.println("\"np\":" + np + ",");
			out.println("\"pairs\":[");
			for (int i = 0; i < np; i++) {
				out.println("[" +pairs[i][0] + "," + pairs[i][1] + "," + means[i]+"]" + (i < np-1 ?",":""));
			}
			out.println("],");
			
			out.println("\"covars\":[");
			boolean first = true;
			for (int i = 0; i < np; i++) {
				for (int j = i; j < np; j++) {
					if (covar0[i][j] != 0) {
						out.print((first ? "":",\n") +"[" + i + "," + j + "," + covar0[i][j] +"]");
						first = false;
					}
				}
			}
			out.println("]");
			out.println("}");
		}
		
		void fromJSON(File matrixFile) throws IOException, JSONException {
			String json = BeautiDoc.load(matrixFile);
			JSONObject doc = new JSONObject(json);
			
			
	        int n = doc.getInt("n");
	        int np = doc.getInt("np");
			taxaNames = new String[n];
			tipHeights = new double[n];

			pairs = new int[np][2];
			means = new double[np];
			covar0 = new double[np][np];
			
			monophylyConstraints = new ArrayList<>();

			JSONArray names = doc.getJSONArray("taxanames");
			for (int i = 0; i < n; i++) {
				taxaNames[i] = names.getString(i);
			}
	        
			JSONArray tipheights_ = doc.getJSONArray("tipheights");
			for (int i = 0; i < n; i++) {
				tipHeights[i] = tipheights_.getDouble(i);
			}

			JSONArray monophylyConstraints_ = doc.getJSONArray("monophylyConstraints");
			for (int i = 0; i < monophylyConstraints_.length(); i++) {
				Object o = monophylyConstraints_.get(i);
				if (o instanceof JSONArray) {
					JSONArray names_ = (JSONArray) o;
					
					List<Taxon> taxa2 = new ArrayList<>();
					for (int j = 0; j < names_.length(); j++) {
						taxa2.add(new Taxon(names_.getString(j)));
					}
					MRCAPrior prior = new MRCAPrior();
					TaxonSet taxonset = new TaxonSet(taxa2);
					prior.taxonsetInput.setValue(taxonset, prior);
					monophylyConstraints.add(prior);
				} else {
					throw new IllegalArgumentException("Expected JSON array in monophylyConstraints");
				}
			}

			
			JSONArray pairs_ = doc.getJSONArray("pairs");
			for (int i = 0; i < np; i++) {
				Object o = pairs_.get(i);
				if (o instanceof JSONArray) {
					JSONArray pair = (JSONArray) o;
					pairs[i][0] = pair.getInt(0);
					pairs[i][1] = pair.getInt(1);
					means[i] = pair.getDouble(2);
				} else {
					throw new IllegalArgumentException("Expected JSON array in pairs");
				}
			}
			
			JSONArray covars_ = doc.getJSONArray("covars");
			for (int i = 0; i < covars_.length(); i++) {
				Object o = covars_.get(i);
				if (o instanceof JSONArray) {
					JSONArray pair = (JSONArray) o;
					int i1 = pair.getInt(0);
					int i2 = pair.getInt(1);
					double var = pair.getDouble(2);
					covar0[i1][i2] = var;
					covar0[i2][i1] = var;					
				} else {
					throw new IllegalArgumentException("Expected JSON array in covars");
				}
				
			}
		}
		
	} // class MultivariateNormal
	
	@Override
	public void run() throws Exception {
        FastTreeSet treeSet = new TreeAnnotator().new FastTreeSet(treesInput.get().getAbsolutePath(), burnInPercentageInput.get());
        treeSet.reset();
        Tree tree = treeSet.next();
        treeSet.reset();
        taxa = tree.getTaxaNames();


		File orderFile = orderInput.get();
		String str = BeautiDoc.load(orderFile);
		String [] strs = (str.split("[\\s,]+"));
		if (strs.length != taxa.length) {
			throw new IllegalArgumentException("order file " + orderFile.getName() + " contains " + strs.length + " taxa, but tree contains " + taxa.length + " taxa");
		}
		int [] order = new int[taxa.length];
		for (int i = 0; i < order.length; i++) {
			order[i] = indexOf(taxa, strs[i]);
		}
        
        List<Tree> treeList = new ArrayList<>();
        while (treeSet.hasNext()) {
            tree = treeSet.next();
            treeList.add(tree);
            System.err.print('<');
        }
        Tree[] trees = treeList.toArray(new Tree[]{});
		
		
		
		
		
		if (seedInput.get() != null) {
			Randomizer.setSeed(seedInput.get());
		}
		useLogTransform = useLogTransformInput.get();
		useMonoConstraints = useMonoConstraintsInput.get();
		if ((outputInput.get() == null || outputInput.get().getName().equals("[[none]]"))) {
			throw new IllegalArgumentException("out flag (for output file) must be specified");
		}

		Log.warning(trees.length + " trees fetched");
		process(outputInput.get(), trees, sampleCountInput.get(), order);

		
		Log.warning("Done");
	}
	


	private int indexOf(String[] taxa, String name) {
		for (int i = 0; i < taxa.length; i++) {
			if (taxa[i].equals(name)) {
				return i;
			}
		}
		throw new IllegalArgumentException("Could not find taxon name \"" + name + "\"");
	}

	public void process(OutFile output, Tree [] trees2, int sampleCount, int [] order) throws IOException, IllegalDimension {
		int [][] pairs = getPairs(order);
		Tree tree = trees2[0];
		
		int n = tree.getLeafNodeCount();
		double[][][] fDist = new double[n][n][trees2.length];

		int k = 0;
		for (Tree tree2: trees2) {
			calcDistance(tree2.getRoot(), k++, fDist, 1.0, new ArrayList<Integer>(), new ArrayList<Double>());
		}

		int np = pairs.length;
		
		double[][] mean = new double[n][n];
		double[][] sd = new double[n][n];
		for (int i = 0; i < n; i++) {
			Arrays.fill(mean[i], 1e100);
		}
		for (int i = 0; i < np; i++) {
			int [] pair = pairs[i];
			int i1 = pair[0];
			int j1 = pair[1];
			mean[i1][j1] = mean(fDist[i1][j1]);
			mean[j1][i1] = mean[i1][j1];
			sd[i1][j1] = sd(fDist[i1][j1]);
			sd[j1][i1] = sd[i1][j1];
		}
		double[][] covar0 = new double[np][np];
		double [] means = new double[np];

		for (int i = 0; i < np; i++) {
			int i1 = pairs[i][0];
			int i2 = pairs[i][1];
			for (int j = i + 1; j < np; j++) {
				int j1 = pairs[j][0];
				int j2 = pairs[j][1];
				covar0[i][j] = cor(mean[i1][i2], mean[j1][j2], fDist[i1][i2], fDist[j1][j2]);
				covar0[j][i] = covar0[i][j];
			}
 			means[i] = mean[i1][i2];
 			System.out.println("mean[" + i +"] = " + means	[i] + " = mean["+i1+"][" + i2 +"]");
			covar0[i][i] = sd[i1][i2]*sd[i1][i2];
		}
		
		
		double [] tipheights = new double[tree.getLeafNodeCount()];
		for (int i = 0; i < tipheights.length; i++) {
			tipheights[i] = tree.getNode(i).getHeight();
		}
		
		List<MRCAPrior> monophylyconstraints = getMonophylyConstraints(trees2);

		double [][] L = getL(covar0);
		
		System.out.println(Arrays.toString(tree.getTaxaNames()));
//		for (int i = 0; i < n; i++) {
//			System.out.println(Arrays.toString(mean[i]));
//			System.out.println(Arrays.toString(sd[i]));
//		}


		sampleToFile(tree.getTaxaNames(), tipheights, output, sampleCount, means, L, pairs, monophylyconstraints);
	}
	



	private List<MRCAPrior> getMonophylyConstraints(Tree[] trees) {
		if (!useMonoConstraints) {
			// return empty list
			return new ArrayList<>();
		}
		
		String[] taxonNames = trees[0].getTaxaNames();
		Taxon[] taxa = new Taxon[taxonNames.length];
		for (int i = 0; i < taxa.length; i++) {
			taxa[i] = new Taxon(taxonNames[i]);
		}
		
		List<MRCAPrior> monophylyconstraints = new ArrayList<>();
		// pick up all clades with 100% support
		CladeSet set = new CladeSet(trees[0]);
		for (int i = 1; i < trees.length; i++) {
			set.add(trees[i]);
		}
		for (int i = 0; i < set.getCladeCount(); i++) {
			int freq = set.getFrequency(i);
			if (freq >= trees.length *0.995) {
				BitSet cladeSet = set.get(i);
				List<Taxon> taxa2 = new ArrayList<>();
				for (int j = 0; j < taxa.length; j++) {
					if (cladeSet.get(j)) {
						taxa2.add(taxa[j]);
					}
				}
				if (taxa2.size() < taxa.length) {
					MRCAPrior prior = new MRCAPrior();
					TaxonSet taxonset = new TaxonSet(taxa2);
					prior.taxonsetInput.setValue(taxonset, prior);
					monophylyconstraints.add(prior);
				}
			}
		}
		
		// initialise mono constraints
		TaxonSet taxonset = new TaxonSet();
		List<Sequence> sequences = new ArrayList<>();
		for (String s : taxonNames) {
			taxonset.taxonsetInput.get().add(new Taxon(s));
			sequences.add(new Sequence(s, "?"));
		}
		Alignment data = new Alignment(sequences, "nucleotide");
		for (MRCAPrior prior : monophylyconstraints) {
			Tree tree = new Tree();
			tree.initByName("taxonset", new TaxonSet(data));
			prior.treeInput.setValue(tree, prior);
			prior.isMonophyleticInput.setValue(true, prior);
			prior.initAndValidate();
		}

		
		Log.warning("Found " + monophylyconstraints.size() + " monophyly constraints");
		return monophylyconstraints;
	}


	void sampleToFile(String [] taxaNames, double [] tipheights, OutFile output, int sampleCount, double [] means, double [][] L, int [][] pairs) throws FileNotFoundException {
		sampleToFile(taxaNames, tipheights, output, sampleCount, means, L, pairs, new ArrayList<>());
	}
	
	void sampleToFile(String [] taxaNames, double [] tipheights, OutFile output, int sampleCount, double [] means, double [][] L, int [][] pairs, List<MRCAPrior> monophylyconstraints) throws FileNotFoundException {
		PrintStream out = new PrintStream(output);
		Log.warning("Writing to file " + output.getPath());

		ByteArrayOutputStream baos = new ByteArrayOutputStream();
		Log.info = new PrintStream(baos);
		

		int N1 = sampleCount / 10;
		int N2 = sampleCount / 50;
		
		this.monophylyconstraints = new ArrayList<>();
		monophylyconstraints.clear();
		for (MRCAPrior prior : monophylyconstraints) {
			List<Integer> taxonIndices = new ArrayList<>();
			for (String taxonName : prior.taxonsetInput.get().asStringList()) {
				int i = indexOf(taxonName, taxaNames);
				taxonIndices.add(i);
			}
			this.monophylyconstraints.add(taxonIndices);
		}
		
		
		for (int m = 0; m < sampleCount; m++) {	
			
			//String newick = getNewickSample(taxaNames, tipheights, means, L, pairs);
			//out.println(newick);
			Tree clusterTree = getSample(taxaNames, tipheights, means, L, pairs, monophylyconstraints);
			out.println(clusterTree.getRoot().toNewick());
			if (m % N1 == 0) {
				System.err.print('|');
			} else if (m % N2 == 0) {
				System.err.print('.');
			}
		}
		System.err.println("| "+ (useMonoConstraints ?  rejects + " rejects" : ""));
		out.close();
	}

	private int indexOf(String taxonName, String[] taxaNames) {
		for (int i = 0; i < taxaNames.length; i++) {
			if (taxonName.equals(taxaNames[i])) {
				return i;
			}
		}
		throw new RuntimeException("Could not find taxon " + taxonName + " in " + Arrays.toString(taxaNames));
	}

	double logWeight;

	
	Tree getSample(String [] taxaNames, double [] tipheights, double[] means, double[][] L, int[][] pairs) {
		return getSample(taxaNames, tipheights, means, L, pairs, new ArrayList<>());
	}
	
	int rejects = 0;
	String [] taxaNames;
	double [][] distance;
	List<List<Integer>> monophylyconstraints;
	Node [] nodes;
	Node root;
	int lastNode = -1;
	
//	String getNewickSample(String [] taxaNames, double [] tipheights, double[] means, double[][] L, int[][] pairs) {
//		if (nodes == null) {
//			nodes = new Node[taxaNames.length * 2 - 1];
//			for (int i = 0; i < nodes.length; i++) {
//				nodes[i] = new Node();
//				nodes[i].setNr(i);
//			}
//			for (int i = 0; i < taxaNames.length; i++) {
//				nodes[i].setID(taxaNames[i]);
//				nodes[i].setHeight(tipheights[i]);
//			}
//		} else {
//			for (int i = 0; i < nodes.length; i++) {
//				nodes[i].removeAllChildren(false);
//			}
//		}
//    	lastNode = 2*taxaNames.length - 2;
//		
//		while (true) {
//			int n = taxaNames.length;
//			double[][] matrix = new double[n][n];
//
//			logWeight = sample(matrix, means, L, pairs);
//	
//			this.distance = matrix;
//			this.taxaNames = taxaNames;
//			String newick = buildClusterer();
//			if (newick != null) {
//				// if non-null, the clustering succeeded and is conform monophylyconstraints
//				return newick;
//			}
//		}
//	}
//	
	
	Tree getSample(String [] taxaNames, double [] tipheights, double[] means, double[][] L, int[][] pairs, List<MRCAPrior> monophylyconstraints) {
		while (true) {
			int n = taxaNames.length;
			double[][] matrix = new double[n][n];
			//sample(matrix, distr);
			logWeight = sample(matrix, means, L, pairs);
			//sample(matrix, mean, sd);
	
			Distance distance = new Distance() {
				@Override
				public double pairwiseDistance(int taxon1, int taxon2) {
					return matrix[taxon1][taxon2];
				}
			};
	
			TaxonSet taxonset = new TaxonSet();
			List<Sequence> sequences = new ArrayList<>();
			for (String s : taxaNames) {
				taxonset.taxonsetInput.get().add(new Taxon(s));
				sequences.add(new Sequence(s, "?"));
			}
			Alignment data = new Alignment(sequences, "nucleotide");
	
			
//			if (monophylyconstraints.size() > 0) {
//				clusterTree = new ConstrainedClusterTree();
//				clusterTree.initByName("clusterType", "single", "taxa", data, "distance", distance, "constraint", monophylyconstraints);
//			} else {
				if (clusterTree == null) {
					clusterTree = new ClusterTree();
				}
				clusterTree.initByName("clusterType", "single", "taxa", data, "distance", distance);
				//clusterTree.initByName("clusterType", "average2", "taxa", data, "distance", distance);
//			}
			for (int i = 0; i < clusterTree.getLeafNodeCount(); i++) {
				clusterTree.getNode(i).setHeight(tipheights[i]);
			}
			boolean fitsConstraints = true;
			for (MRCAPrior prior : monophylyconstraints) {
				if (!isMonophyletic(clusterTree, prior)) {
					fitsConstraints = false;
					break;
				}
			}
			if (fitsConstraints) {
				return clusterTree;
			}
			rejects++;
		}
	}

	Tree clusterTree = null;

	
	private boolean isMonophyletic(Tree tree, MRCAPrior prior) {
		prior.treeInput.setValue(tree, prior);
		prior.isMonophyleticInput.setValue(true, prior);
		prior.initAndValidate();
		return Double.isFinite(prior.calculateLogP());
	}


	/** get Cholesky decomposition cleaned up for NaN and Infinities **/
	public double[][] getL(double[][] covar0) throws IllegalDimension {
		if (useCorrelationsInput.get()) {
			
			int nonzero = 0, zero = 0;
			double zeroLevel = 0; // was 0.01
			for (int i = 0; i < covar0.length; i++) {
				for (int j = 0; j < i; j++) {
					if (covar0[i][j] < zeroLevel) {
						covar0[i][j] = 0;
						zero++;
					} else {
						nonzero++;
					}
				}
				for (int j = i+1; j < covar0[i].length; j++) {
					if (covar0[i][j] < zeroLevel) {
						covar0[i][j] = 0;
						zero++;
					} else {
						nonzero++;
					}
				}
			}
			for (int i = 0; i < covar0.length; i++) {
				System.out.println(Arrays.toString(covar0[i]));
			}
			System.out.println(zero + " zeros, " + nonzero +" non zeros");
			
//			if (matrixOutputInput.get() != null && 
//				!matrixOutputInput.get().getName().equals("[[none]]")) {
//				return null;
//			}

			
			CholeskyDecomposition chol = new CholeskyDecomposition(covar0);
			double [][] L = chol.getL();
			for (int i = 0; i < L.length; i++) {
				for (int j = 0; j < L[0].length; j++) {
					if (Double.isNaN(L[i][j])||Double.isInfinite(L[i][j])) {
						L[i][j] = 0;
					}
				}
			}
			return L;
		} else {
			// remove all correlations
			for (int i = 0; i < covar0.length; i++) {
				for (int j = 0; j < i; j++) {
					covar0[i][j] = 0;
				}
				for (int j = i+1; j < covar0[i].length; j++) {
					covar0[i][j] = 0;
				}
			}
			
//			if (matrixOutputInput.get() != null && 
//				!matrixOutputInput.get().getName().equals("[[none]]")) {
//				return null;
//			}

			
			CholeskyDecomposition chol = new CholeskyDecomposition(covar0);
			double [][] L = chol.getL();
			for (int i = 0; i < L.length; i++) {
				for (int j = 0; j < L[0].length; j++) {
					if (Double.isNaN(L[i][j])||Double.isInfinite(L[i][j])) {
						L[i][j] = 0;
					}
				}
			}

			return L;
		}
	}
		


	/** return pairs of nodes (=matrix entries) in cube orders compatible with trees2 **/
	private int[][] getPairs(int [] order) throws IOException {
		
		int [][] pairs = new int[order.length - 1][2];
		for (int i = 0; i < order.length - 1; i++) {
			pairs[i][0] = order[i];
			pairs[i][1] = order[i+1];
		}
		Log.warning("\n" + pairs.length + " pairs (= entries in matrix), " + (order.length-1)+ " from cube");
		
		return pairs;
	}


	private double cor(double mean1, double mean2, double[] x1, double[] x2) {
		if (x1.length == 1) {
			return 0;
		}
		double sum = 0;
		for (int i = 0; i < x1.length; i++) {
			sum += (x1[i] - mean1) * (x2[i] - mean2);
		}
		return sum / (x1.length - 1.0);
	}

	private double cor(double mean1, double mean2, double[] x1, double[] x2, double [] weights) {
		if (x1.length == 1) {
			return 0;
		}
		double sum = 0;
		for (int i = 0; i < x1.length; i++) {
			sum += (x1[i] - mean1) * (x2[i] - mean2) * weights[i];
		}
		return sum / (x1.length - 1.0);
	}

	//	private double logMean(double[] ds) {
//		double[] logD = new double[ds.length];
//		for (int i = 0; i < ds.length; i++) {
//			logD[i] = Math.log(ds[i]);
//		}
//		return mean(logD);
//	}
//
//	private double logSd(double[] ds) {
//		double[] logD = new double[ds.length];
//		for (int i = 0; i < ds.length; i++) {
//			logD[i] = Math.log(ds[i]);
//		}
//		return sd(logD);
//	}
//
//	private void sample(double[][] matrix, MultivariateNormalDistribution distr) {
//		int n = matrix.length;
//		double [] sample = distr.sample();
//		int k = 0;
//		for (int i = 0; i < n; i++) {
//			for (int j = i + 1; j < n; j++) {
//				matrix[i][j] = sample[k];
//				matrix[j][i] = matrix[i][j];
//				k++;
//			}
//		}
//		
//	}
//
	double sample(double[][] matrix, double[] means, double[][] l, int [][] pairs) {
		int n = matrix.length;
		int np = pairs.length;
		
		// select n-1 pairs forming a cube
		int [] done = new int[n];
		int [] pairIndex = new int[n-1];
		// select random pair
		int x = Randomizer.nextInt(np);
		pairIndex[0] = x;
		done[pairs[x][0]] = 1;
		done[pairs[x][1]] = 1;
		
		// select n-2 adjacent pairs
		for (int i = 1; i < n-1; i++) {
			do {
				x = Randomizer.nextInt(np);
			} while ((done[pairs[x][0]]>0 && done[pairs[x][1]]>0) ||
					(done[pairs[x][0]] == 0 && done[pairs[x][1]] == 0));
			done[pairs[x][0]]++;
			done[pairs[x][1]]++;
			pairIndex[i] = x;
		}
		
		double [] sample = new double[np];
		double logP = 0;
		for (int i = 0; i < pairIndex.length; i++) {
			int k = pairIndex[i];
			sample[k] = Randomizer.nextGaussian();
			logP += -0.5 * sample[k] * sample[k];
		}
		logP += sample.length * -Math.log(2 * Math.PI);
		
		for (int i = 0; i < n; i++) {
			Arrays.fill(matrix[i], Double.POSITIVE_INFINITY);
		}
		for (int i = 0; i < pairIndex.length; i++) {
			int k = pairIndex[i];
			int i1 = pairs[k][0];
			int i2 = pairs[k][1];
			double h = mult(l, k, sample) + means[k];
			if (useLogTransform) {
				h = Math.exp(h);
			}
			matrix[i1][i2] = h;
			matrix[i2][i1] = h;
		}
		return logP;

//		double [] sample = new double[np];
//		double logP = 0;
//		for (int i = 0; i < sample.length; i++) {
//			sample[i] = Randomizer.nextGaussian();
//			logP += -0.5 * sample[i] * sample[i];
//		}
//		logP += sample.length * -Math.log(2 * Math.PI);
//		
//		for (int i = 0; i < n; i++) {
//			Arrays.fill(matrix[i], 1e100);
//		}
//		for (int i = 0; i < np; i++) {
//			int i1 = pairs[i][0];
//			int i2 = pairs[i][1];
//			matrix[i1][i2] = mult(l, i, sample) + means[i];
//			matrix[i2][i1] = matrix[i1][i2];
//		}
//		return logP;
	}

	
	private double mult(double[][] ds, final int index, double[] sample) {
		double sum = 0;
		for (int i = 0; i < ds.length; i++) {
			sum += ds[index][i] * sample[i];
		}
		return sum;
	}

//	private void sample(double[][] matrix, double[][] mean, double[][] sd) {
//		int n = matrix.length;
//		for (int i = 0; i < n; i++) {
//			for (int j = i + 1; j < n; j++) {
//				matrix[i][j] = sample(mean[i][j], sd[i][j]);
//				matrix[j][i] = matrix[i][j];
//			}
//		}
//	}

//	private double sample(double mean, double sd) {
//		return Randomizer.nextGaussian() * sd + mean;
//		// double logSigma = sd;
//		// double logMean = Math.log(mean) - (0.5 * logSigma * logSigma);
//		// double r = Math.exp(Randomizer.nextGaussian() * logSigma + logMean);
//		// return r;
//	}

	private double mean(double[] ds, double [] weights) {
		double sum = 0;
		for (int i = 0;  i < ds.length; i++) {
			sum += ds[i] * weights[i];
		}
		return sum / ds.length;
	}

	private double mean(double[] ds) {
		double sum = 0;
		for (int i = 0;  i < ds.length; i++) {
			sum += ds[i];
		}
		return sum / ds.length;
	}

	private double sd(double[] ds, double [] weights) {
		double mean = mean(ds, weights);
		double mean2 = mean * mean;

		double sum = 0;
		for (int i = 0;  i < ds.length; i++) {
			sum += (ds[i] * ds[i] - mean2) * weights[i];
		}
		double sd = Math.sqrt(sum / (ds.length - 1));
		return sd;
	}

	private double sd(double[] ds) {
		double mean = mean(ds);
		double mean2 = mean * mean;

		double sum = 0;
		for (int i = 0;  i < ds.length; i++) {
			sum += ds[i] * ds[i] - mean2;
		}
		double sd = Math.sqrt(sum / (ds.length - 1));
		return sd;
	}

	/**
	 * calculate the distance between leafs in a consensus tree and update the
	 * distance matrix weighted with the relative frequency of the tree
	 * 
	 * @param node:
	 *            current node
	 * @param nOrder:
	 *            mapping of node label to [0,...,NrOfLeafs-1]
	 * @param fDistMatrix:
	 *            distance matrix to be updated
	 * @param fWeight:
	 *            relative consensus tree frequency
	 * @param iLabel:
	 *            used to report set of leafs in sub tree below node
	 * @param fLength:
	 *            used to report set of lengths from current node to leafs in
	 *            iLabel
	 */
	private void calcDistance(Node node, int index, double[][][] fDistMatrix, double fWeight, List<Integer> iLabel,
			List<Double> fLength) {
		if (node == null) {
			return;
		}
		if (node.isLeaf()) {
			// iLabel.add(nOrder[node.m_iLabel]);
			iLabel.add(node.getNr());
			// fLength.add(node.m_fLength);
			fLength.add(node.getLength() + node.getHeight());
		} else {
			iLabel.clear();
			List<Integer> iLeft = new ArrayList<>();
			List<Integer> iRight = new ArrayList<>();
			List<Double> fLeft = new ArrayList<>();
			List<Double> fRight = new ArrayList<>();
			calcDistance(node.getChild(0), index, fDistMatrix, fWeight, iLeft, fLeft);
			calcDistance(node.getChild(1), index, fDistMatrix, fWeight, iRight, fRight);
			for (int i = 0; i < iLeft.size(); i++) {
				int i1 = iLeft.get(i);
				double f1 = fWeight * fLeft.get(i);
				for (int j = 0; j < iRight.size(); j++) {
					int i2 = iRight.get(j);
					double f2 = f1 + fWeight * fRight.get(j);
					if (useLogTransform) {
						f2 = Math.log(f2);
					}
					fDistMatrix[i1][i2][index] = f2;
					fDistMatrix[i2][i1][index] = f2;
				}
			}
			for (int i = 0; i < fLeft.size(); i++) {
				iLabel.add(iLeft.get(i));
				// fLength.add(fLeft.get(i) + node.m_fLength);
				fLength.add(fLeft.get(i) + node.getLength());
				// fLength.add(fLeft.get(i) + 1.0f);
			}
			for (int i = 0; i < fRight.size(); i++) {
				iLabel.add(iRight.get(i));
				// fLength.add(fRight.get(i) + node.m_fLength);
				fLength.add(fRight.get(i) + node.getLength());
				// fLength.add(fRight.get(i) + 1.0f);
			}
		}
	} // calcDistance


	public Trees2Cube() throws IOException {
	}

	@Override
	public void initAndValidate() {
		// nothing to do
	}

	public static void main(String[] args) throws Exception {
		new Application(new Trees2Cube(), "Trees2Cube", args);
	}

}
