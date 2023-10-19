package cubevb.operators;


import java.util.HashSet;
import java.util.Set;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Log;
import beast.base.evolution.operator.Exchange;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.util.Randomizer;
import cubevb.transform.CubeTransform;

@Description(
        "The narrow exchange is very similar to a rooted-beast.tree nearest-neighbour " +
        "interchange but with the restriction that node height must remain consistent " + 
        "and the cube order is maintained.")
public class CubeExchange extends Exchange {

	final public Input<CubeTransform> cubeTransformInput = new Input<>("cube", "cube transform. If specified uses order as specified by the transform, otherwise use a random order");

	CubeOperatorHelper cubeOperatorHelper;
	CubeTransform transform;
	
    @Override
    public void initAndValidate() {
    	super.initAndValidate();
    	
    	transform = cubeTransformInput.get();
    }

    
    @Override
    public double proposal() {    	
    	Node n = narrow_(treeInput.get(), transform);

    	if (n == null) {
    		return Double.NEGATIVE_INFINITY;
    	}
    	
        if (System.getProperty("beast.debug") != null && System.getProperty("beast.debug").equals("true")) {
        	if (cubeOperatorHelper == null) {
        		cubeOperatorHelper = new CubeOperatorHelper(cubeTransformInput.get(), treeInput.get());
        	}
        	if (!cubeOperatorHelper.treeIsCompatibleWithCube()) {
        		Log.warning("CubeExchange.proposal programmer error. Should not get here");
        		return Double.NEGATIVE_INFINITY;
    		}    		
    	}
    	return 0.0;
    }
    
    /**
     * WARNING: Assumes strictly bifurcating beast.tree.
     */
    private Set<Node> path = new HashSet<>();
    private Set<Node> leftpath = new HashSet<>();
    public Node narrow_(final Tree tree, CubeTransform transform) {

        final int internalNodes = tree.getInternalNodeCount();
        if (internalNodes <= 1) {
            return null;
        }

        int k = Randomizer.nextInt(internalNodes);
        
        Node grandParent = CubeTransform.getPath(tree, transform, k, path, leftpath);
        while (grandParent.getLeft().isLeaf() && grandParent.getRight().isLeaf()) {
        	k = Randomizer.nextInt(internalNodes);
        	grandParent = CubeTransform.getPath(tree, transform, k, path, leftpath);
        }

        Node parentIndex = grandParent.getLeft();
        Node uncle = grandParent.getRight();
        if (parentIndex.getHeight() < uncle.getHeight()) {
            parentIndex = grandParent.getRight();
            uncle = grandParent.getLeft();
        }

        if( parentIndex.isLeaf() ) {
            // tree with dated tips
            return null;
        }

    
        Node i = null;
        if (!path.contains(parentIndex.getLeft())) {
        	i = parentIndex.getLeft();
        }
        if (!path.contains(parentIndex.getRight())) {
        	i = parentIndex.getRight();
        }
        if (i == null) {
        	// should not get here
        	Log.warning("Programmer error in CubeExchange.narrow()");
        	return null;
        }
        exchangeNodes(i, uncle, parentIndex, grandParent);

        return grandParent;
    }
}
