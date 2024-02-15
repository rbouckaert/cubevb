package cubevb.transform;

import beast.base.core.BEASTObject;
import beast.base.core.Description;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.TreeInterface;
import beast.base.util.Randomizer;

@Description("Find random order consistent with given tree")
public class RandomCuberOrder extends BEASTObject implements CubeOrderer {

	@Override
	public void initAndValidate() {
	}

	@Override
	public int[] order(TreeInterface tree, double [] heights) {
		Node [] internalNodes = new Node[tree.getLeafNodeCount() - 1];
		int [] order = new int[tree.getNodeCount()];
		getRandomCube(tree.getLeafNodeCount(), tree.getRoot(), order, new int[]{0}, internalNodes, heights);
		return order;
	}
	
	/** recursively determines order of tree compatible with (random) planar drawing of the tree 
	 * @param node current node in tree
	 * @param order is populated with ordering of leaf nodes
	 * @param heights of the cube
	 * @param index keeps track of next order index to add
	 * @param internalNodes nodes for which heights are recorded
	 */
	protected void getRandomCube(int offset, Node node, int [] order,  int [] index, Node [] internalNodes, double [] heights) {
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
			getRandomCube(offset, first, order, index, internalNodes, heights);
			heights[index[0]] = node.getHeight();
			internalNodes[index[0]] = node;
			order[index[0]+offset] = node.getNr();
			index[0]++;
			getRandomCube(offset, second, order, index, internalNodes, heights);
		}		
	}

}
