package cubevb.transform;

import java.util.Arrays;
import java.util.List;

import beast.base.core.BEASTObject;
import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.evolution.alignment.Alignment;
import beast.base.evolution.alignment.Sequence;
import beast.base.evolution.distance.Distance;
import beast.base.evolution.distance.JukesCantorDistance;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.TreeInterface;

@Description("Order cube by recursively minimizing distance between outside of already ordered sub trees")
public class DistanceCubeOrderer extends BEASTObject implements CubeOrderer {
	final public Input<Alignment> alignmentInput = new Input<>("data", "alignment used for pairwise distances", Validate.REQUIRED);
	enum ordermode {outside, average, minimum};
	final public Input<ordermode> modeInput = new Input<>("mode", "mode with which to calculate distance between clades. One of " + Arrays.toString(ordermode.values()), ordermode.minimum, ordermode.values());
	
		
	@Override
	public void initAndValidate() {
	}

	private int [] map;
	private TreeInterface tree;
	private ordermode mode;
	
	@Override
	public int[] order(TreeInterface tree, double [] heights) {
		this.tree = tree;
		Alignment data = alignmentInput.get();
		JukesCantorDistance d = new JukesCantorDistance();
		d.setPatterns(data);
		
		map = new int[tree.getNodeCount()];
		Node [] nodes = tree.getNodesAsArray();
		for (int i = 0; i < tree.getLeafNodeCount(); i++) {
			Node node = nodes[i];
			map[node.getNr()] = nodenr(node.getID(), data);
		}
		System.err.println(Arrays.toString(map));
		
		boolean [] rotate = new boolean[tree.getNodeCount()];
		mode = modeInput.get();
		switch (mode) {
		case outside : traverse(tree.getRoot(), d, rotate); break;
		case average :
		case minimum :
			c = new int[tree.getNodeCount()][];
			traverseAverage(tree.getRoot(), d, rotate); break;
		}
		
		int [] order = new int[tree.getNodeCount()]; 
		setOrder(tree.getLeafNodeCount(), tree.getRoot(), order, new int[] {0}, heights, rotate, false);
		System.err.println(Arrays.toString(order));
		return order;
	}

	private int nodenr(String taxon, Alignment data) {
		List<Sequence> seqs = data.sequenceInput.get();
		for (int i = 0; i < seqs.size(); i++) {
			if (seqs.get(i).taxonInput.get().equals(taxon)) {
				return i;
			}
		}
		throw new IllegalArgumentException("Cannot find taxon " + taxon + " in alignment '" + data.getID() + "'");
	}

	private void setOrder(int offset, Node node, int[] order, int [] index,double[] heights, boolean[] rotate, boolean rotated) {
		if (node.isLeaf()) {
			order[index[0]] = node.getNr();
		} else {
			Node first = null, second = null;
			if ((!rotated && !rotate[node.getNr()]) || (rotated && rotate[node.getNr()])) {
				first = node.getChild(0);
				second = node.getChild(1);
			} else {
				first = node.getChild(1);
				second = node.getChild(0);
			}
			if (rotate[node.getNr()]) {
				rotated = !rotated;
			}
			setOrder(offset, first, order, index, heights, rotate, rotated);
			heights[index[0]] = node.getHeight();
			order[index[0]+offset] = node.getNr();
			index[0]++;
			setOrder(offset, second, order, index, heights, rotate, rotated);
		}		
	}

	int [][] c;
	private void traverseAverage(Node node, Distance d, boolean [] rotate) {
		if (node.isLeaf()) {
			c[node.getNr()] = new int[] {node.getNr()}; 
		} else {
			Node left = node.getLeft();
			Node right = node.getRight();
			
			traverseAverage(left, d, rotate);
			traverseAverage(right, d, rotate);
			
			if (left.isLeaf()) {
				if (right.isLeaf()) {
					c[node.getNr()] = new int[] {left.getNr(), right.getNr()};
					return;
				}
			}
			if (!right.isLeaf()) {
				int [] clade1 = c[left.getNr()];
				int [] clade2 = c[right.getLeft().getNr()];
				int [] clade3 = c[right.getRight().getNr()];
				double dist12 = averagePairwiseDistance(d, clade1, clade2);	
				double dist13 = averagePairwiseDistance(d, clade1, clade3);	
				if (dist12 <= dist13) {
					rotate[right.getNr()] = false;
				} else {
					rotate[right.getNr()] = true;
				}
			}
				
			//} else if (right.isLeaf()) {
			if (!left.isLeaf()) {
				int [] clade1 = c[right.getNr()];
				int [] clade2 = c[left.getLeft().getNr()];
				int [] clade3 = c[left.getRight().getNr()];
				double dist21 = averagePairwiseDistance(d, clade2, clade1);	
				double dist31 = averagePairwiseDistance(d, clade3, clade1);	
				if (dist21 <= dist31) {
					rotate[left.getNr()] = true;
				} else {
					rotate[left.getNr()] = false;
				}
			}
//			} else {
//				int [] clade1 = c[left.getLeft().getNr()];
//				int [] clade2 = c[left.getRight().getNr()];
//				int [] clade3 = c[right.getLeft().getNr()];
//				int [] clade4 = c[right.getRight().getNr()];
//				double dist13 = averagePairwiseDistance(d, clade1, clade3);	
//				double dist14 = averagePairwiseDistance(d, clade1, clade4);	
//				double dist23 = averagePairwiseDistance(d, clade2, clade3);	
//				double dist24 = averagePairwiseDistance(d, clade2, clade4);	
//				if (dist13 <= dist14 && dist13 <= dist23 && dist13 <= dist24) {
//					System.err.println("d11");
//					rotate[left.getNr()] = true;
//					rotate[right.getNr()] = false;
//				} else if (dist14 <= dist13 && dist14 <= dist23 && dist14 <= dist24) {
//					System.err.println("d12");
//					rotate[left.getNr()] = true;
//					rotate[right.getNr()] = true;
//				} else if (dist23 <= dist13 && dist23 <= dist14 && dist23 <= dist24) {
//					System.err.println("d21");
//					rotate[left.getNr()] = false;
//					rotate[right.getNr()] = false;
//				} else {
//					System.err.println("d22");
//					rotate[left.getNr()] = false;
//					rotate[right.getNr()] = true;
//				}
//			}
			
			int [] clade1 = c[left.getNr()];
			int [] clade2 = c[left.getNr()];
			int [] clade = new int[clade1.length + clade2.length];
			c[node.getNr()] = clade;
			System.arraycopy(clade1, 0, clade, 0, clade1.length);
			System.arraycopy(clade2, 0, clade, clade1.length, clade2.length);
			
//			Node upper, lower;
//			// TODO deal with case where leaf is higher than other child
//			// this only happens with tip dated trees
//			if (right.isLeaf() || (!left.isLeaf() && left.getHeight() > right.getHeight())) {
//				upper = left;
//				lower = right;
//			} else {
//				lower = left;
//				upper = right;
//			}
//			
//			int [] clade1 = c[lower.getNr()];
//			int [] clade2 = c[upper.getLeft().getNr()];
//			int [] clade3 = c[upper.getRight().getNr()];
//			double dist12 = averagePairwiseDistance(d, clade1, clade2);	
//			double dist13 = averagePairwiseDistance(d, clade1, clade3);	
//			if (dist12 <= dist13) {
//				rotate[upper.getLeft().getNr()] = true;
//			} else {
//				rotate[upper.getLeft().getNr()] = false;
//			}
//			int [] clade = new int[clade1.length + clade2.length + clade3.length];
//			c[node.getNr()] = clade;
//			System.arraycopy(clade1, 0, clade, 0, clade1.length);
//			System.arraycopy(c[upper.getNr()], 0, clade, clade1.length, clade2.length +  + clade3.length);
			
			
		}		
	}
	
	private double averagePairwiseDistance(Distance d,int [] clade1, int [] clade2) {
		double distance = 0;
		switch (mode) {
		case average:
			for (int i : clade1) {
				for (int j : clade2) {
					distance += d.pairwiseDistance(map[i], map[j]);
				}
			}
			distance /= (clade1.length * clade2.length);
			return distance;
		case minimum:
			distance = Double.POSITIVE_INFINITY;
			for (int i : clade1) {
				for (int j : clade2) {
					distance = Math.min(d.pairwiseDistance(map[i], map[j]), distance);
				}
			}
			return distance;
		default:
			return 0;
		}
	}

	private int [] traverse(Node node, Distance d, boolean [] rotate) {
		if (node.isLeaf()) {
			return new int[] {node.getNr(), node.getNr()};
		} else {
			int [] clade1 = traverse(node.getLeft(), d, rotate);
			int [] clade2 = traverse(node.getRight(), d, rotate);
			double dist11 = pairwiseDistance(d, map[clade1[0]], map[clade2[0]]);	
			double dist12 = pairwiseDistance(d, map[clade1[0]], map[clade2[1]]);	
			double dist21 = pairwiseDistance(d, map[clade1[1]], map[clade2[0]]);	
			double dist22 = pairwiseDistance(d, map[clade1[1]], map[clade2[1]]);
			rotate[node.getNr()] = true;
			if (dist11 <= dist12 && dist11 <= dist21 && dist11 <= dist22) {
				System.err.println("d11");
				rotate[node.getLeft().getNr()] = true;
				rotate[node.getRight().getNr()] = false;
				return new int[] {clade1[1], clade2[1]};
			} else if (dist12 <= dist11 && dist12 <= dist21 && dist12 <= dist22) {
				System.err.println("d12");
				rotate[node.getLeft().getNr()] = true;
				rotate[node.getRight().getNr()] = true;
				return new int[] {clade1[1], clade2[0]};
			} else if (dist21 <= dist11 && dist21 <= dist12 && dist21 <= dist22) {
				System.err.println("d21");
				rotate[node.getLeft().getNr()] = false;
				rotate[node.getRight().getNr()] = false;
				return new int[] {clade1[0], clade2[1]};
			} else {
				System.err.println("d22");
				rotate[node.getLeft().getNr()] = false;
				rotate[node.getRight().getNr()] = true;
				return new int[] {clade1[0], clade2[0]};
			}
		}
	}

	private double pairwiseDistance(Distance d, int taxon1, int taxon2) {
		double distance = d.pairwiseDistance(taxon1, taxon2);
		System.err.println(tree.getNode(taxon1).getID() + " " + tree.getNode(taxon2).getID() + " " + distance);
		return distance;
	}

}
