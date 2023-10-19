package cubevb.operators;



import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.evolution.operator.SubtreeSlide;
import beast.base.evolution.tree.Tree;
import beast.base.inference.util.InputUtil;
import beast.base.util.Randomizer;
import cubevb.transform.CubeTransform;

@Description("Moves the height of an internal node along the branch. " +
        "If it moves up, it can exceed the root and become a new root. " +
        "If it moves down, it chooses the branch to slide down into so " +
        "it remains compatible with a cube ordering.")
public class CubeSubTreeSlide extends SubtreeSlide {
	final public Input<CubeTransform> cubeTransformInput = new Input<>("cube", "cube transform. If specified uses order as specified by the transform, otherwise use a random order");

	CubeOperatorHelper cubeOperatorHelper;
	
    @Override
    public void initAndValidate() {
    	super.initAndValidate();
    }

    @Override
    public double proposal() {
        final Tree tree = (Tree) InputUtil.get(treeInput, this);

        final boolean markClades = markCladesInput.get();
        // 1. choose a random node avoiding root
        final int nodeCount = tree.getNodeCount();
        if (nodeCount == 1) {
        	// test for degenerate case (https://github.com/CompEvol/beast2/issues/887)
    		return Double.NEGATIVE_INFINITY;    		
        }

        
        int indexInOrder = -1;
        indexInOrder = Randomizer.nextInt(tree.getLeafNodeCount() - 1);


        // 2. choose a delta to move
        final double delta = getDelta();
        
        double logHR = CubeTransform.propose(tree, cubeTransformInput.get(), markClades, indexInOrder, delta);
    	
        if (System.getProperty("beast.debug") != null && System.getProperty("beast.debug").equals("true")) {
        	if (cubeOperatorHelper == null) {
        		cubeOperatorHelper = new CubeOperatorHelper(cubeTransformInput.get(), treeInput.get());
        	}
        	if (!cubeOperatorHelper.treeIsCompatibleWithCube()) {
        		return Double.NEGATIVE_INFINITY;
    		}    		
    	}

    	return logHR;
    }

    private double getDelta() {
        if (!gaussianInput.get()) {
            return (Randomizer.nextDouble() * size) - (size / 2.0);
        } else {
            return Randomizer.nextGaussian() * size;
        }
    }
}
