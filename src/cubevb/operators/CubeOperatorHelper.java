package cubevb.operators;

import java.util.Map;
import java.util.HashMap;

import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.evolution.tree.TreeInterface;
import cubevb.transform.CubeTransform;

public class CubeOperatorHelper {
	
	private CubeTransform transform;
	private TreeInterface tree;
	private Map<Integer, Integer> node2order;
	private Node lastRootSeen;
	
	public CubeOperatorHelper(CubeTransform transform, TreeInterface tree) {
		this.transform = transform;
		this.tree = tree;
		node2order = initMap();
	}

	public boolean treeIsCompatibleWithCube() {
    	int [] order = transform.getOrder();

    	boolean [] done = new boolean[tree.getNodeCount()];
    	Node current = tree.getNode(order[0]);
    	done[current.getNr()] = true;
    	for (int i = 1; i < tree.getLeafNodeCount()-1; i++) {
    		Node n = tree.getNode(order[i]);
        	done[n.getNr()] = true;
    		Node p = n.getParent();
    		Node other = otherChild(p, n);
        	do {
        		p = n.getParent();
        		if (p != null) {
        			other = otherChild(p, n);
        			n = p;
        		}
        	} while (!done[other.getNr()] && p != null);
        	if (!done[other.getNr()] || p == null) {
        		return false;
        	}
        	done[n.getNr()] = true;
    	}
    	
    	// for debugging:
    	if (System.getProperty("beast.debug") != null && System.getProperty("beast.debug").equals("true") && initMap() == null) {
    		initMap();
    		throw new RuntimeException("Programmer error: CubeOperatorHelper.treeIsCompatibleWithCube() passed while it should not");
    	}
    	
    	return true;
	}
	
	private Map<Integer, Integer> initMap() {
		int [] order = transform.getOrder();
		Node [] nodes = tree.getNodesAsArray();

		Map<Integer, Integer> node2order = new HashMap<>();
		for (int i = 0; i < tree.getLeafNodeCount() - 1; i++) {
			Node n = getCommonAncestorNode(nodes[order[i]], nodes[order[i+1]]);
			if (node2order.containsKey(n.getNr())) {
				// node found that is MRCA of 2 pairs in order
				// this means that the tree cannot be compatible with the ordering
				return null;
				// throw new IllegalArgumentException("Tree not compatible with node order");
			}
			node2order.put(n.getNr(), i);
		}

		lastRootSeen = tree.getRoot();
		return node2order;
	}
	
    protected Node getCommonAncestorNode(Node n1, Node n2) {
    	boolean [] nodesTraversed = new boolean[tree.getNodeCount()];
    	// assert n1.getTree() == n2.getTree();
        if( ! nodesTraversed[n1.getNr()] ) {
            nodesTraversed[n1.getNr()] = true;
        }
        if( ! nodesTraversed[n2.getNr()] ) {
            nodesTraversed[n2.getNr()] = true;
        }
        while (n1 != n2) {
	        double h1 = n1.getHeight();
	        double h2 = n2.getHeight();
	        if ( h1 < h2 ) {
	            n1 = n1.getParent();
	            if( ! nodesTraversed[n1.getNr()] ) {
	                nodesTraversed[n1.getNr()] = true;
	            }
	        } else if( h2 < h1 ) {
	            n2 = n2.getParent();
	            if( ! nodesTraversed[n2.getNr()] ) {
	                nodesTraversed[n2.getNr()] = true;
	            }
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
                } else {
                    n = n2 = n.getParent();
                }
	            if( ! nodesTraversed[n.getNr()] ) {
	                nodesTraversed[n.getNr()] = true;
	            } 
	        }
        }
        return n1;
    }	
	
    private Node otherChild(final Node parent, final Node child) {
        if (parent.getLeft().getNr() == child.getNr()) {
            return parent.getRight();
        } else {
            return parent.getLeft();
        }
    }

    @Override
    public String toString() {
    	String [] taxa = ((Tree)tree).getTaxaNames();
    	StringBuilder b = new StringBuilder();
    	int [] order = transform.getOrder();
    	for (int i = 0; i < taxa.length; i++) {
    		b.append(taxa[order[i]]);
    		if (i < taxa.length - 1) {
    			b.append(", ");
    		}
    	} 
    	return b.toString();
    }
}
