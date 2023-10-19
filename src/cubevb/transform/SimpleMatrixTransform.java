package cubevb.transform;


import java.io.FileNotFoundException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import org.apache.commons.math.MathException;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.evolution.distance.Distance;
import beast.base.evolution.tree.ClusterTree;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.evolution.tree.TreeInterface;
import beast.base.evolution.tree.TreeUtils;
import beast.base.inference.StateNode;
import beast.base.inference.distribution.Normal;
import beast.base.util.Randomizer;

@Description("Transforms distance matrix to/from tree and capture correlations through a spanning tree")
public class SimpleMatrixTransform extends Transform {
	final public Input<TreeInterface> treeInput = new Input<>("tree", "phylogenetic beast.tree with sequence data in the leafs", Validate.REQUIRED);

	protected TreeInterface tree;
	protected int [][] pairs;

	private List<double[]> trace;
	private double [] mean;
	private double [] stdev;
	private double[][] covar0;
	private double [] threshold;
	private int [] from, to;

	@Override
	public void initAndValidate() {
		super.initAndValidate();
		tree = treeInput.get();
		parameter.setDimension(tree.getLeafNodeCount()*2-3);
	}
	
	@Override
	public double forward(double x) {
		return forward()[0];
	}

	
	
	@Override
	public double backward(double z) {
		if (mean == null) {
			if (Boolean.getBoolean("beast.debug")) {
				for (int i = 0; i < pairs.length; i++) {
					beast.base.core.Log.warning(Arrays.toString(pairs[i]));
				}
				beast.base.core.Log.warning(Arrays.toString(tree.getTaxonset().asStringList().toArray()));
				// store trace in /tmp/beast.log
				try {
					PrintStream out;
					out = new PrintStream("/tmp/beast.log");
					out.print("Sample\t");
					for (int i = 0; i < trace.get(0).length; i++) {
						out.print("trace" + i + "\t");
					}
					out.println();
					for (int i = 0; i < trace.size(); i++) {
						out.print(i+"\t");
						double [] v = trace.get(i);
						for (double d : v) {
							out.print(d+"\t");
						}
						out.println();
					}
					out.close();
				} catch (FileNotFoundException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
			}
			
			initMeanStdevCovar();
			calcProbOffDiagonal();
			calcSpanningTree();
		}
		
		int n = tree.getInternalNodeCount();
		int [][] myPairs = new int [n][];
		int [] src = new int[n];
		double [] rnd = new double[n];
		for (int i = 0; i < n; i++) {
			double r = Randomizer.nextGaussian();
			rnd[i] = r;
			if (r < threshold[i]) {
				myPairs[i] = pairs[n + i];
				src[i] = n + i;
			} else {
				myPairs[i] = pairs[i];
				src[i] = i;
			}
		}
		
		double [] h = new double[tree.getInternalNodeCount()];
		double r = rnd[0];
		int start = src[to[0]];
		h[start % n] = mean[start] + r * stdev[start];
		
		for (int i = 1; i < myPairs.length; i++) {
			r =rnd[i];
			int k = src[from[i]];
			int j = src[to[i]];
			double rho = covar0[j][k] / (stdev[j] * stdev[k]);
			h[j % n] = mean[j] + rho * (stdev[j]/stdev[k])*(h[k % n] - mean[k]);
			h[j % n] += r * (1-rho * rho) * stdev[j];
		}
		
		for (int i = 0; i < myPairs.length; i++) {
			h[i] = Math.exp(h[i]);
		}
		matrix2Tree(myPairs, h, tree);
		//System.out.println(tree.getRoot().getHeight());
		return 0;
	}
	
	
	private void calcSpanningTree() {
		int n = tree.getInternalNodeCount();
		to = new int[n * 2];
		from = new int[n * 2];
		// find max mean as starting point
		int iMax = 0;
		double max = mean[0];
		for (int i = 1; i < n; i++) {
			if (mean[i] > max) {
				max = mean[i];
				iMax = i;
			}
		}
		to[0] = iMax;
		
		Tuple [] tuples = new Tuple[n * (n-1)/2];
		int k = 0;
		for (int i = 0; i < n; i++) {
			for (int j = i + 1; j < n; j++) {
				tuples[k++] = new Tuple(i, j, covar0[i][j]);
			}
		}
		
		Arrays.sort(tuples);
		
		
		// find largest covar to already used items
		boolean [] done = new boolean[n];
		done[iMax] = true;
		for (int i = 1; i < n; i++) {
			k = 0;
			while ((!done[tuples[k].i] && !done[tuples[k].j])||
					(done[tuples[k].i] &&  done[tuples[k].j])){
				k++;
			}
			if (done[tuples[k].i]) {
				from[i] = tuples[k].i;
				to[i] = tuples[k].j;
				from[i+n] = tuples[k].i;
				to[i+n] = tuples[k].j + n;
				done[tuples[k].j] = true;
			} else {
				from[i] = tuples[k].j;
				to[i] = tuples[k].i;
				from[i+n] = tuples[k].j;
				to[i+n] = tuples[k].i + n;
				done[tuples[k].i] = true;
			}
		}
	}
	
	private void calcProbOffDiagonal() {
		int n = tree.getInternalNodeCount();
		double [] probOffDiagonal = new double[n];
		
		for (double [] v : trace) {
			for (int i = 0; i < n; i++) {
				if (v[i] > v[n + i]) {
					probOffDiagonal[i]++;
				}
			}
		}
		for (int i = 0; i < n; i++) {
			probOffDiagonal[i] /= trace.size();
		}

		if (Boolean.getBoolean("beast.debug")) {
			double [] probBelowDiagonal = new double[n];
			
			for (double [] v : trace) {
				for (int i = 0; i < n; i++) {
					if (v[i] < v[n + i]) {
						probBelowDiagonal[i]++;
					}
				}
			}
			for (int i = 0; i < n; i++) {
				probBelowDiagonal[i] /= trace.size();
			}
			beast.base.core.Log.warning("probOffDiagonal   = " +Arrays.toString(probOffDiagonal));
			beast.base.core.Log.warning("probBelowDiagonal = " +Arrays.toString(probBelowDiagonal));
		}
		
		threshold = new double[n];
		Normal normal = new Normal();
		normal.initAndValidate();
		for (int i = 0; i < n; i++) {
			try {
				threshold[i] = normal.inverseCumulativeProbability(probOffDiagonal[i]);
			} catch (MathException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
		
	}

	private void initMeanStdevCovar() {
		int n = pairs.length;
		mean = new double[n];
		stdev = new double[n];
		covar0 = new double[n][n];
		TransformUtils.initMeanStdevCovar(trace, n, mean, stdev, covar0);
	}



	@Override
	public double jacobian_det(double x) {
		return 0;
	}

	@Override
	public double backward(int dim, double oldValue, double newValue) {
		// TODO Auto-generated method stub
		throw new RuntimeException("not implemented yet");
		// return 0;
	}
	
	@Override
	public double[] forward() {
		// tree to matrix
		if (pairs == null) {
			pairs = findpairs();
			trace = new ArrayList<>();
		}
		double [] heights = new double[pairs.length];
		
		for (int i = 0; i < pairs.length; i++) {
			int [] pair = pairs[i];
			Node node1 = tree.getNode(pair[0]);
			Node node2 = tree.getNode(pair[1]);
			Set<String> leafNodes = new HashSet<>();
			leafNodes.add(node1.getID());
			leafNodes.add(node2.getID());
			Node mrca = TreeUtils.getCommonAncestorNode((Tree) tree, leafNodes);
			heights[i] = mrca.getHeight();
		}
		
		for (int i = 0; i < heights.length; i++) {
			heights[i] = Math.log(heights[i]);
		}
		trace.add(heights);
		
		heights = new double[getDimension()];
		heights[0] = Randomizer.nextGaussian();
		return heights;
	}

	protected int[][] findpairs() {
		int [] order = new int[tree.getNodeCount()];
		List<int[]> pairs = new ArrayList<>();
		// get a cube = diagonal
		getCube(tree.getRoot(), pairs, order, new int[1]);
		
		// get 2-off-diagonal
		int [] pair = new int[2];
		pair[0] = order[tree.getLeafNodeCount() - 1];
		pair[1] = order[1];
		pairs.add(pair);
		for (int i = 2; i < tree.getLeafNodeCount(); i++) {
			pair = new int[2];
			pair[0] = order[i-2];
			pair[1] = order[i];
			pairs.add(pair);
		}
		Object [] os = pairs.toArray();
		int [][] pairs_ = new int[os.length][];
		for (int i = 0; i < os.length; i++) {
			pairs_[i] = (int []) os[i];
		}
		return pairs_;
	}
	
	
	private void getCube(Node node, List<int[]> pairs,
			int [] order,
			int [] index) {
		if (node.isLeaf()) {
			order[index[0]] = node.getNr();
			if (index[0] > 0) {
				int [] pair = new int[2];
				pair[0] = order[index[0] - 1];
				pair[1] = order[index[0]];
				pairs.add(pair);
			}
		} else {
			Node first = null, second = null;
			if (Randomizer.nextBoolean()) {
				first = node.getChild(0);
				second = node.getChild(1);
			} else {
				first = node.getChild(1);
				second = node.getChild(0);
			}
			getCube(first, pairs, order, index);
			order[index[0]+tree.getLeafNodeCount()] = node.getNr();
			index[0]++;
			getCube(second, pairs, order, index);
		}		
	}
	
	@Override
	public double[] backward(double[] z) {
		double [] values = new double[z.length];
		for (int i = 0; i < z.length; i++) {
			values[i] = Math.exp(z[i]);
		}
		matrix2Tree(pairs, values, tree);
		return values;		
	}

	private void matrix2Tree(int[][] pairs, double[] values, TreeInterface tree) {
		int n = tree.getTaxonset().getTaxonCount();
		double [] distanceMatrix = new double[n * n];
		Arrays.fill(distanceMatrix, Double.MAX_VALUE);
		for (int i = 0; i < values.length; i++) {
			int [] pair = pairs[i];
			distanceMatrix[pair[0] * n + pair[1]] = values[i] * 2;
			distanceMatrix[pair[1] * n + pair[0]] = values[i] * 2;
		}
		
		Distance distance = new Distance() {
			@Override
			public double pairwiseDistance(int taxon1, int taxon2) {
				return distanceMatrix[taxon1 * n + taxon2];
			}
		};
		
		ClusterTree clusterer = new ClusterTree();
		clusterer.initByName("clusterType", "single",
				"distance", distance,
				"taxonset", tree.getTaxonset(),
				"initial", tree);
	}

	@Override
	public double jacobian_det(int dim, double x) {
		return 0;
	}

	@Override
	public double jacobian_det(double[] x) {
		return 0;
	}
	
	@Override
	public StateNode getTransformNode() {
		return parameter;
	}
	
	@Override
	public StateNode getStateNode() {
		return (StateNode) tree;
	}
	
	@Override
	public int getDimension() {		
		return 1;
	}

}
