package cubevb.transform;

import beast.base.core.Description;
import beast.base.evolution.tree.TreeInterface;

@Description("Method for finding taxon order of a cube")
public interface CubeOrderer {
	
	public int [] order(TreeInterface tree, double [] heights);

}
