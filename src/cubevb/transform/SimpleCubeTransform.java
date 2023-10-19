package cubevb.transform;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.inference.StateNode;
import beast.base.core.Input.Validate;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.evolution.tree.TreeInterface;
import beast.base.evolution.tree.TreeUtils;
import beast.base.util.Randomizer;

@Description("Represents a tree as a cube through a RealParameter and capture correlations through a spanning tree")
public class SimpleCubeTransform extends Transform {
	final public Input<TreeInterface> treeInput = new Input<>("tree", "phylogenetic beast.tree with sequence data in the leafs", Validate.REQUIRED);
	final public Input<String> orderInput = new Input<>("order", "comma separated list of taxa names in tree in desirded cube order. If not specified, a location "
			+ "provider is used, and if that is not specified either, a random order is used.");
//    final public Input<LocationProvider> locationInput = new Input<>("location", "location provider for informing cube order. Should not be used when order is "
//    		+ "already specified.");
 	
	protected TreeInterface tree;
	protected int [] order;
	protected int [] inverseorder;
	
	
	private List<double[]> trace;
	private double [] mean;
	private double [] stdev;
	private double[][] covar0;
	private int [] to, from;

	public int [] getOrder() {
		if (order == null) {
			forward();
		}
		return order;
	}
	
	public int [] getInverseOrder() {
		if (order == null) {
			forward();
		}
		return inverseorder;
	}
	
	
	@Override
	public void initAndValidate() {
		super.initAndValidate();
		tree = treeInput.get();
		parameter.setDimension(tree.getLeafNodeCount()-1);
		Node[] nodes = tree.getNodesAsArray();
		
		Map<String, Integer> map = new LinkedHashMap<>();
		for (int i = 0; i < tree.getLeafNodeCount(); i++) {
			map.put(nodes[i].getID(), i);
		}

		if (orderInput.get() != null) {
//			if(	locationInput.get()!= null) {
//				throw new IllegalArgumentException("Specify at most one of order or location, but not both.");
//			}
			
			order = new int[nodes.length];
			String str = orderInput.get().trim();
			String [] strs = str.split(",");
			if (strs.length != tree.getLeafNodeCount()) {
				throw new IllegalArgumentException("number of taxa in tree does not match number of taxa in order");
			}
			for (int i = 0; i < strs.length; i++) {
				if (!map.containsKey(strs[i])) {
					throw new IllegalArgumentException("could not find taxon " + strs[i] + " in tree (typo?)");
				}
				order[i] = map.get(strs[i]);
			}
			for (int i = strs.length; i < order.length; i++) {
				order[i] = i;
			}
		}
	}
	
	
	@Override
	public double[] forward() {
		if (order == null) {
			Node [] internalNodes = new Node[tree.getLeafNodeCount() - 1];
			order = new int[tree.getNodeCount()];
			getCube(tree.getRoot(), order, new int[]{0}, internalNodes);
			inverseorder = new int[tree.getNodeCount()];
			for (int i = 0; i < order.length; i++) {
				inverseorder[order[i]] = i;
			}
			
			trace = new ArrayList<>();
		}

		double [] heights = new double[tree.getInternalNodeCount()];
		tree2LogCube(order, tree, heights);
		trace.add(heights);
		
		heights = new double[getDimension()];
		heights[0] = Randomizer.nextGaussian();
		return heights;
	}

	
	protected void tree2LogCube(int[] order, TreeInterface tree, double[] heights) {
		for (int j = 0; j < heights.length; j++) {
			Node node1 = tree.getNode(order[j]);
			Node node2 = tree.getNode(order[j+1]);
			Set<String> leafNodes = new HashSet<>();
			leafNodes.add(node1.getID());
			leafNodes.add(node2.getID());
			Node mrca = TreeUtils.getCommonAncestorNode((Tree) tree, leafNodes);
			heights[j] = mrca.getHeight();
			if (heights[j] < 1e-5) {
				heights[j] = 1e-5;
			}
			heights[j] = Math.log(heights[j]);
		}
	}

	/** recursively determines order of tree compatible with (random) planar drawing of the tree 
	 * @param node current node in tree
	 * @param order is populated with ordering of leaf nodes
	 * @param heights of the cube
	 * @param index keeps track of next order index to add
	 * @param internalNodes nodes for which heights are recorded
	 */
	protected void getCube(Node node, int [] order, int [] index, Node [] internalNodes) {
		if (node.isLeaf()) {
			order[index[0]] = node.getNr();
		} else {
			Node first = null, second = null;
			if (Randomizer.nextBoolean()) {
				first = node.getChild(0);
				second = node.getChild(1);
			} else {
				first = node.getChild(1);
				second = node.getChild(0);
			}
			getCube(first, order, index, internalNodes);
			internalNodes[index[0]] = node;
			order[index[0]+tree.getLeafNodeCount()] = node.getNr();
			index[0]++;
			getCube(second, order, index, internalNodes);
		}		
	}
	
	
	@Override
	public double[] forward(double[] x) {
		double [] values = new double[x.length];
		values[0] = Randomizer.nextGaussian();
		return values;
	}

	@Override
	public double[] backward(double[] z) {
		double [] values = new double[z.length];
		for (int i = 0; i < z.length; i++) {
			values[i] = Math.exp(z[i]);
		}
		cube2Tree(order, values, tree);
		return values;
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

	
	protected Node cube2Tree(int [] order, double [] values, TreeInterface tree) {
		int N = tree.getLeafNodeCount();
		
		Node [] nodes = tree.getNodesAsArray();
		
		// TODO: reserve memory only once, and check dimension if not already reserved
		double [] heights = new double[tree.getNodeCount() + 1];
		int [] parent_ = new int[tree.getNodeCount()];
		int [] left = new int[tree.getNodeCount() + 1];
		int [] right = new int[tree.getNodeCount() + 1];
		Arrays.fill(parent_, -1);
		Arrays.fill(left, -1);
		Arrays.fill(right, -1);
		
		int origin_ = tree.getNodeCount();
		heights[origin_] = Double.MAX_VALUE;
		left[origin_] = order[0];
		parent_[order[0]] = origin_;
		int root_ = order[0]; // should be updated later on
		int current_;
		int next_ = N;
		for (int i = 0; i < N - 1; i++) {
			current_ = order[i];
			double target = values[i];
			while (target > heights[parent_[current_]]) {
				current_ = parent_[current_];
			}
			int parent = parent_[current_];
			int internal = order[next_];
			if (left[parent] == current_) {
				left[parent] = internal;
			} else {
				// note, then right[parent] == current_
				right[parent] = internal;
			}
			parent_[internal] = parent;
			left[internal] = current_;
			parent_[current_] = internal;
			if (internal == origin_) {
				root_ = current_;
			}
			right[internal] = order[i+1];
			parent_[order[i+1]] = internal;
			heights[internal] = target;
			next_++;
		}
		
		List<Node> dirtyNodes = new ArrayList<>();
		List<Node> heightChanged = new ArrayList<>();
		for (int i = N; i < nodes.length; i++) {
			Node node = nodes[i];
			if (Math.abs(node.getHeight() - heights[i]) > 1e-10) {
				node.setHeight(heights[i]);
				heightChanged.add(node);
			}
			if (node.getLeft().getNr() != left[i]) {
				node.setLeft(nodes[left[i]]);
				nodes[left[i]].setParent(node);
				dirtyNodes.add(node);
			} 
			if (node.getRight().getNr() != right[i]) {
				node.setRight(nodes[right[i]]);
				nodes[right[i]].setParent(node);
				dirtyNodes.add(node);
			}
			if (parent_[i] == nodes.length) {
				node.setParent(null);
			}
		}
		if (left[nodes.length] != nodes.length-1) {
			int treeroot = left[nodes.length];
			((Tree) tree).setRoot(nodes[treeroot]);
		}
		
		// TODO: the following crudely enforces mapping onto positive branch lengths
		// this should happen at node.setHeight above -- only a subset of nodes should be checked
		for (Node node : heightChanged) {
			if (!node.isRoot()&&node.getLength() <= 1e-16) {
				node.setHeight(node.getParent().getHeight() - 1e-15);
			}
			final Node nodeLeft = node.getLeft();
			if (nodeLeft.getLength() <= 1e-16) {
				nodeLeft.setHeight(node.getHeight() - 1e-15);
			}
			
			final Node nodeRight = node.getRight();
			if (nodeRight.getLength() <= 1e-16) {
				nodeRight.setHeight(node.getHeight() - 1e-15);
			}
		}
		boolean foundWrongHeight = false;
		do {
			foundWrongHeight = false;
			for (int i = N; i < nodes.length-1; i++) {
				if (nodes[i].getLength() <= 1e-16) {
					nodes[i].setHeight(nodes[i].getParent().getHeight() - 1e-15);
					foundWrongHeight = true;
				}
			}
		} while (foundWrongHeight);
		
		int max = 0;
		for (Node node : dirtyNodes) {
			int i = 0;
			while (node.getParent() != null) {
				node = node.getParent();
				i++;
			}
			max = Math.max(i,  max);
		}
		
//		if (left[nodes.length] != nodes.length-1) {
//			String newick2 = tree.toString();
//			for (int i = N; i < nodes.length; i++) {
//				nodes[i].removeAllChildren(false);
//			}
//			Node root = new Node();
//			root.setHeight(Double.MAX_VALUE);
//			root.addChild(nodes[order[0]]);
//	
//			Node current;
//			int next = N;
//			for (int i = 0; i < N - 1; i++) {
//				current = nodes[order[i]];
//				double target = values[i];
//				while (target > current.getParent().getHeight()) {
//					current = current.getParent();
//				}
//				Node parent = current.getParent();
//				parent.removeChild(current);
//				Node internal = nodes[order[next]];
//				parent.addChild(internal);
//				internal.addChild(current);
//				internal.addChild(nodes[order[i+1]]);
//				internal.setHeight(target);
//				next++;
//			}
//			
//			root = root.getChild(0);
//			root.setParent(null);
//			((Tree)tree).setRoot(root);
//			//return root;
//
//			String newick = tree.toString();
//			if (!newick.equals(newick2)) {
////				System.err.println("Something went wrong: \n" + newick + "\n" + newick2);
//				cube2Tree(order, values, tree);
//			}
//		}
		return tree.getRoot();
	}

	
	@Override
	public double forward(double x) {
		forward();
		return 0;
	}

	@Override
	public double backward(double z) {
		if (mean == null) {
			initMeanStdevCovar();
			calcSpanningTree();
		}
		
		double [] h = new double[tree.getInternalNodeCount()];
		double r = Randomizer.nextGaussian();
		h[to[0]] = mean[to[0]] + r * stdev[to[0]];
		
		for (int i = 1; i < mean.length; i++) {
			r = Randomizer.nextGaussian();
			int j = to[i];
			int k = from[i];
			double rho = covar0[j][k] / (stdev[j] * stdev[k]);
			h[j] = mean[j] + rho * (stdev[j]/stdev[k])*(h[k] - mean[k]);
			h[j] += r * (1-rho * rho) * stdev[j];
		}
		
		for (int i = 0; i < mean.length; i++) {
			h[i] = Math.exp(h[i]);
		}
		cube2Tree(order, h, tree);
		return 0;
	}

	private void calcSpanningTree() {
		int n = tree.getInternalNodeCount();
		to = new int[n];
		from = new int[n];
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
				done[tuples[k].j] = true;
			} else {
				from[i] = tuples[k].j;
				to[i] = tuples[k].i;
				done[tuples[k].i] = true;
			}
		}
	}

	private void initMeanStdevCovar() {
		int n = tree.getInternalNodeCount();
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
		double logHR = propose((Tree)tree, this, false, dim, Math.exp(newValue) - Math.exp(oldValue));
		if (Double.isInfinite(logHR)) {
			return Double.NEGATIVE_INFINITY;
		}
		return newValue;
	}

	@Override
	public double jacobian_det(int dim, double x) {
		return x;
	}


    static public double propose(Tree tree, SimpleCubeTransform transform, boolean markClades, int indexInOrder, double delta) {
        Set<Node> path = new HashSet<>(); 
        Set<Node> leftpath = new HashSet<>(); 
        Node parent = getPath(tree, transform, indexInOrder, path, leftpath);
        final double newHeight = parent.getHeight() + delta;
    	
        Node targetChild = leftpath.contains(parent.getLeft()) ? parent.getLeft() : parent.getRight();
        Node otherChild = getOtherChildS(parent, targetChild);
        Node grantParent = parent.getParent();

        final double oldHeight = parent.getHeight();
        // 3. if the move is up
        if (delta > 0) {

            // 3.1 if the topology will change
            if (grantParent != null && grantParent.getHeight() < newHeight) {

                // slide up in the topology
                do {
                	if (isLeftOf(transform, indexInOrder, getOtherChildS(grantParent, parent))) {
                		if (!grantParent.isRoot()) {
                			Node ggParent = grantParent.getParent();
                			ggParent.removeChild(grantParent);
                			ggParent.addChild(parent);
            			
                			grantParent.removeChild(parent);
                			grantParent.addChild(targetChild);
                			grantParent.makeDirty(Tree.IS_FILTHY);
            		
                			parent.removeChild(targetChild);
                			parent.addChild(grantParent);
                			targetChild = grantParent;
                		} else {
                			// create a new root
                			grantParent.removeChild(parent);
                			grantParent.addChild(targetChild);
                			grantParent.makeDirty(Tree.IS_FILTHY);
                			
                			parent.removeChild(targetChild);
                			parent.addChild(grantParent);

                			targetChild = grantParent;
                            parent.setParent(null);
                            tree.setRoot(parent);
                		}
                	} else {
                		if (!grantParent.isRoot()) {
                            Node ggParent = grantParent.getParent();
                            ggParent.removeChild(grantParent);
                            ggParent.addChild(parent);
                			grantParent.makeDirty(Tree.IS_FILTHY);
                			
                			grantParent.removeChild(parent);
                			grantParent.addChild(otherChild);
                		
                			parent.removeChild(otherChild);
                			parent.addChild(grantParent);
                		} else {
                			// create a new root
                			grantParent.removeChild(parent);
                			grantParent.addChild(otherChild);
                			grantParent.makeDirty(Tree.IS_FILTHY);

                			parent.removeChild(otherChild);
                			parent.addChild(grantParent);
                			
                            parent.setParent(null);
                            tree.setRoot(parent);
                		}
                	}
            		// parent.setHeight((Math.max(parent.getLeft().getHeight(), parent.getRight().getHeight()) + grantParent.getHeight())/2);
            		
                    otherChild = getOtherChildS(parent, targetChild);
                    grantParent = parent.getParent();
                } while (grantParent != null && grantParent.getHeight() < newHeight);
                	
            }
        }
        // 4 if we are sliding the subtree down.
        else {

            // 4.0 is it a valid move?
            if (targetChild.getHeight() > newHeight) {
        		return Double.NEGATIVE_INFINITY;    		
            }

            // 4.1 will the move change the topology
            if (otherChild.getHeight() > newHeight) {
                final Node newChild = getNodeIntersectingHeightOnPath(otherChild, newHeight, path);

                // if no valid destinations then return a failure
                if (newChild == null) {
            		return Double.NEGATIVE_INFINITY;    		
                }

                final Node newParent = newChild.getParent();

                // 4.1.1 if p was root
                if (parent.isRoot()) {
                    // new root is CiP
                    replaceS(parent, otherChild, newChild);
                    replaceS(newParent, newChild, parent);

                    otherChild.setParent(null);
                    tree.setRoot(otherChild);

                } else {
                    replaceS(parent, otherChild, newChild);
                    replaceS(grantParent, parent, otherChild);
                    replaceS(newParent, newChild, parent);
                }

                if( markClades ) {
                    // make dirty the path from the (down) moved node back up to former parent.
                    Node n = parent;
                    while( n != otherChild ) {
                        n.makeDirty(Tree.IS_FILTHY); // JH
                        n = n.getParent();
                    }
                }
            }
        }
        
        parent.setHeight(newHeight);


    	return 0.0;    	
    }
    
    
    public static boolean isLeftOf(SimpleCubeTransform transform, int indexInOrder, Node node) {
    	while (!node.isLeaf()) {
    		node = node.getLeft();
    	}
    	int [] inverseOrder = transform.getInverseOrder();
    	int nodeInOrder = inverseOrder[node.getNr()];
    	if (nodeInOrder <= indexInOrder) {
    		return true;
    	}
		return false;
	}
    

	static private Node getNodeIntersectingHeightOnPath(Node node, double height, Set<Node> path) {
        final Node parent = node.getParent();
        
        if (parent == null) {
        	// can happen with non-standard non-mutable trees
        	return null;
        }

        if (parent.getHeight() < height) { 
        	return null;
        }

        if (node.getHeight() < height) {
        	if (path.contains(node)) {
        		return node;
        	} else {
        		return null;
        	}
        }

        if (node.isLeaf()) {
            // makes sense, since node has no children
        	// if this leaf had a branch above crossing `height`
        	// it would return in the previous `if` statement
            return null;
        } else {
        	if (path.contains(node.getLeft())) {
        		return getNodeIntersectingHeightOnPath(node.getLeft(), height, path);
        	}
        	if (path.contains(node.getRight())) {
                return getNodeIntersectingHeightOnPath(node.getRight(), height, path);
        	}
        }
        return null;
    }    
    
    static public Node getPath(TreeInterface tree, SimpleCubeTransform transform, int indexInOrder, Set<Node> nodesTraversed, Set<Node> leftPath) {
		Node [] nodes = tree.getNodesAsArray();
		int [] order = transform.getOrder();
		Node n1 = nodes[order[indexInOrder]];
		Node n2 = nodes[order[indexInOrder + 1]];
    	
    	nodesTraversed.clear();
    	leftPath.clear();
    	nodesTraversed.add(n1);
        leftPath.add(n1);
    	nodesTraversed.add(n2);
        while (n1 != n2) {
	        double h1 = n1.getHeight();
	        double h2 = n2.getHeight();
	        if ( h1 < h2 ) {
	            n1 = n1.getParent();
	            nodesTraversed.add(n1);
	            leftPath.add(n1);
	        } else if( h2 < h1 ) {
	            n2 = n2.getParent();
	            nodesTraversed.add(n2);
	        } else {
	            //zero length branches hell
	            Node n;
	            double b1 = n1.getLength();
	            double b2 = n2.getLength();
	            if( b1 > 0 ) {
	                n = n2;
	            } else { // b1 == 0
	                if( b2 > 0 ) {
	                    n = n1;
	                } else {
	                    // both 0
	                    n = n1;
	                    while( n != null && n != n2 ) {
	                        n = n.getParent();
	                    }
	                    if( n == n2 ) {
	                        // n2 is an ancestor of n1
	                        n = n1;
	                    } else {
	                        // always safe to advance n2
	                        n = n2;
	                    }
	                }
	            }
	            if( n == n1 ) {
                    n = n1 = n.getParent();
    	            leftPath.add(n1);
                } else {
                    n = n2 = n.getParent();
                }
	            nodesTraversed.add(n);
	        }
        }
        return n1;
    }		
    /**
     * @param parent the parent
     * @param child  the child that you want the sister of
     * @return the other child of the given parent.
     */
    static protected Node getOtherChildS(final Node parent, final Node child) {
        if (parent.getLeft().getNr() == child.getNr()) {
            return parent.getRight();
        } else {
            return parent.getLeft();
        }
    }

    /**
     * replace child with another node
     *
     * @param node
     * @param child
     * @param replacement
     */
    static public void replaceS(final Node node, final Node child, final Node replacement) {
    	node.removeChild(child);
    	node.addChild(replacement);
        node.makeDirty(Tree.IS_FILTHY);
        replacement.makeDirty(Tree.IS_FILTHY);
    }
}
