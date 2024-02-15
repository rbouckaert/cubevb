package cubevb.transform;

import beast.base.core.BEASTObject;
import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.evolution.alignment.Alignment;
import beast.base.evolution.distance.Distance;
import beast.base.evolution.distance.JukesCantorDistance;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.TreeInterface;
import beast.base.util.Randomizer;

@Description("Order cube by recursively minimizing distance between outside of already ordered sub trees")
public class DistanceCubeOrderer extends BEASTObject implements CubeOrderer {
	final public Input<Alignment> alignmentInput = new Input<>("data", "alignment used for pairwise distances", Validate.REQUIRED);
		
	@Override
	public void initAndValidate() {
	}

	@Override
	public int[] order(TreeInterface tree, double [] heights) {
		Alignment data = alignmentInput.get();
		JukesCantorDistance d = new JukesCantorDistance();
		d.setPatterns(data);
		boolean [] rightFirst = new boolean[tree.getNodeCount()];
		traverse(tree.getRoot(), d, rightFirst);
		
		int [] order = new int[tree.getNodeCount()]; 
		setOrder(tree.getLeafNodeCount(), tree.getRoot(), order, new int[] {0}, heights, rightFirst, false);
		return order;
	}

	private void setOrder(int offset, Node node, int[] order, int [] index,double[] heights, boolean[] rightFirst, boolean rotated) {
		if (node.isLeaf()) {
			order[index[0]] = node.getNr();
		} else {
			Node first = null, second = null;
			if ((!rotated && rightFirst[node.getNr()])||(rotated && !rightFirst[node.getNr()])) {
				first = node.getChild(0);
				second = node.getChild(1);
			} else {
				first = node.getChild(1);
				second = node.getChild(0);
				rotated = !rotated;
			}
			setOrder(offset, first, order, index, heights, rightFirst, rotated);
			heights[index[0]] = node.getHeight();
			order[index[0]+offset] = node.getNr();
			index[0]++;
			setOrder(offset, second, order, index, heights, rightFirst, rotated);
		}		
	}

	private int [] traverse(Node node, Distance d, boolean [] rightFirst) {
		if (node.isLeaf()) {
			return new int[] {node.getNr(), node.getNr()};
		} else {
			int [] clade1 = traverse(node.getLeft(), d, rightFirst);
			int [] clade2 = traverse(node.getLeft(), d, rightFirst);
			double dist11 = d.pairwiseDistance(clade1[0], clade2[0]);	
			double dist12 = d.pairwiseDistance(clade1[0], clade2[clade2.length-1]);	
			double dist21 = d.pairwiseDistance(clade1[clade1.length-1], clade2[0]);	
			double dist22 = d.pairwiseDistance(clade1[clade1.length-1], clade2[clade2.length-1]);
			if (dist11 < dist12 && dist11 < dist21 && dist11 < dist22) {
				rightFirst[node.getLeft().getNr()] = false;
				rightFirst[node.getRight().getNr()] = true;
				return new int[] {clade1[1], clade2[1]};
			} else if (dist12 < dist11 && dist12 < dist21 && dist12 < dist22) {
				rightFirst[node.getLeft().getNr()] = false;
				rightFirst[node.getRight().getNr()] = false;
				return new int[] {clade1[1], clade2[0]};
			} else if (dist21 < dist11 && dist21 < dist21 && dist21 < dist22) {
				rightFirst[node.getLeft().getNr()] = true;
				rightFirst[node.getRight().getNr()] = true;
				return new int[] {clade1[0], clade2[1]};
			} else {
				rightFirst[node.getLeft().getNr()] = true;
				rightFirst[node.getRight().getNr()] = false;
				return new int[] {clade1[0], clade2[0]};
			}
		}
	}

}
