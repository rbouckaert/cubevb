package cubevb.transform;


import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

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
import beast.base.util.Randomizer;

@Description("Transforms distance matrix to/from tree")
public class MatrixTransform extends Transform {
	final public Input<TreeInterface> treeInput = new Input<>("tree", "phylogenetic beast.tree with sequence data in the leafs", Validate.REQUIRED);

	protected TreeInterface tree;
	protected int [][] pairs;
	
	@Override
	public void initAndValidate() {
		super.initAndValidate();
		tree = treeInput.get();
		parameter.setDimension(tree.getLeafNodeCount()*2-3);
	}
	
	@Override
	public double forward(double x) {
		throw new RuntimeException("Cannot use single dimension call on " + getClass().getName());
	}

	@Override
	public double backward(double z) {
		throw new RuntimeException("Cannot use single dimension call on " + getClass().getName());
	}

	@Override
	public double jacobian_det(double x) {
		throw new RuntimeException("Cannot use single dimension call on " + getClass().getName());
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
		double [] heights = new double[parameter.getDimension()];
		if (pairs == null) {
			pairs = findpairs();
		}
		
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
		return heights;
	}

	protected int[][] findpairs() {
		int [] order = new int[tree.getNodeCount()];
		List<int[]> pairs = new ArrayList<>();
		// get a cube = diagonal
		getCube(tree.getRoot(), pairs, order, new int[1]);
		
		// get 2-off-diagonal
		for (int i = 2; i < tree.getLeafNodeCount(); i++) {
			int [] pair = new int[2];
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
		return x;
	}

	@Override
	public double jacobian_det(double[] x) {
		double sum = 0;
		for (double d : x) {
			sum += d;
		}
		return sum;
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
		return parameter.getDimension();
	}

}
