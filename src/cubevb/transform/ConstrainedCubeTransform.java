package cubevb.transform;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import beast.base.core.BEASTInterface;
import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Log;
import beast.base.inference.State;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.MRCAPrior;
//import beastlabs.math.distributions.MultiMRCAPriors;

@Description("Represents a tree as a cube through a RealParameter that has monophyletic constraints")
public class ConstrainedCubeTransform extends CubeTransform {
//	final public Input<MultiMRCAPriors> constraintInput = new Input<>("constraints", "set of monophyletic constraints on the tree");
// 	
//
//	/** 
//	 * relHeights are heights relative to the highest of left and right 
//	 * clade that are constrained to be monophyletic
//	 * If both sides are not constrained, relHeights equals the height of 
//	 * the gap in a cube 
//	 **/
//	private double [] relheights;
//	
//	/**
//	 * from and to are arrays indicating the start of gaps to the left till 
//	 * end of gaps to the right of largest monophyletic constraint that
//	 * is next to the current gap
//	 **/
//	private int [] from, to;
//	
//	
//	@Override
//	public double[] forward() {
//		double [] heights = new double[parameter.getDimension()];
//		if (order == null) {
//			// order will be set by initialiseOrder, so we should only get here once
//			initialiseOrder(heights);
//		} else {		
//			tree2Cube(order, tree, heights);
//		}
//		calcRelHeights(heights);
//		for (int i = 0; i < relheights.length; i++) {
//			relheights[i] = Math.log(relheights[i]);
//		}
//		return relheights;
//	}
//	
//	
//	
//	private void initialiseOrder(double [] heights) {
//		Node [] internalNodes = new Node[tree.getLeafNodeCount() - 1];
//		order = new int[tree.getNodeCount()];
////		if (locationInput.get() != null) {
////			LocationProvider d = locationInput.get();
////			
////			// mark all state nodes dirty to guarantee recaclculation of posterior  
////			for (BEASTInterface o : ((BEASTInterface)tree).getOutputs()) {
////				if (o instanceof State) {
////					State state = (State) o;
////					state.setEverythingDirty(true);
////					state.checkCalculationNodesDirtiness();
////				}
////			}
////			
////			orderByDistance(tree.getRoot(), order, heights, new int[]{0}, internalNodes, d);
////		} else {
//			// random order consistent with tree
//			getCube(tree.getRoot(), order, heights, new int[]{0}, internalNodes);
////		}
//		
//		
//		// allocate memory
//		relheights = new double[tree.getLeafNodeCount()-1];
//		from = new int[relheights.length];
//		to = new int[relheights.length];
//		Arrays.fill(from, tree.getNodeCount());
//		Arrays.fill(to, -1);
//		
//		calcBoundaries();
//	}
//
//
//
//	// set up boundaries in `from` and `to` arrays
//	// it goes through all monophyletic constraints
//	// every such constraint implies that both left and right gap next to the set of constrained nodes in the cube
//	// has to be higher than the highest gap inside the set of constrained nodes (if either left of right gap would
//	// be lower, the clade would not be monophyletic any more).
//	// from[i] and to[i]  represent the range of gaps that must be lower than gap i
//	// by going through all the constraints, we can expand the range till its final range covers all constraints
//	private void calcBoundaries() {
//		// taxaList is list of taxa names in order of tree taxon list, not in order of cube 
//		String [] taxaList = new String[tree.getLeafNodeCount()];
//        for(int k = 0; k < taxaList.length; ++k) {
//            Node n = tree.getNode(k);
//            taxaList[k] = n.getID();
//        }
//
//        // reverorder = reverse of order[] array, so revorder[order[i]]==i
//        int [] revorder = new int[order.length];
//		for (int i = 0; i < order.length; i++) {
//			revorder[order[i]] = i;
//		}
//		final int n = tree.getLeafNodeCount();
//		
//		// deal with the monophyletic constraints specified by Newick tree first
//		MultiMRCAPriors priors = constraintInput.get();
//		List<List<Integer>> taxonIDList = parse(priors.newickInput.get(), taxaList);
//		for (List<Integer> IDList : taxonIDList) {
//			int min = order.length;
//			int max = -1;
//			for (int i : IDList) {
//				int j = revorder[i];
//				min = Math.min(min, j);
//				max = Math.max(max, j);
//			}
//			
////			// debug code			
////			if (max - min +1 != IDList.size()) {
////				int h = 3;
////				h++;
////			}
//			updateRange(min, max, n);
//		}
//
//		// deal with the monophyletic constraints specified by MRCA priors
//		List<MRCAPrior> mrcaPriors = priors.calibrationsInput.get();
//		for (MRCAPrior prior : mrcaPriors) {
//			if (prior.isMonophyleticInput.get()) {
//				int min = order.length;
//				int max = -1;
//				List<String> taxa = prior.taxonsetInput.get().asStringList();
//				for (String taxon : taxa) {
//					taxon = taxon.trim();
//					int i = indexOf(taxon, taxaList);
//					if (i == -1 && (taxon.startsWith("'") || taxon.startsWith("\""))) {
//						i = indexOf(taxon.substring(1, taxon.length() - 1), taxaList);
//					}
//		            if (i == -1) {
//		                throw new RuntimeException("Cannot find taxon " + taxon + "  in taxon list");
//		            }
//					int j = revorder[i];
//					min = Math.min(min, j);
//					max = Math.max(max, j);
//				}
//				updateRange(min, max, n);
//			}
//		}
//	}
//	
//	/* expand from-to range for gaps left and right of clade from min to max if necessary 
//	 * n = nr of leafs in tree = size of cube + 1 */
//    private void updateRange(int min, int max, int n) {
//		if (min > 0) {
//			// update gap on the left of the monophyletic constraint
//			// if min == 0 there is no gap to update since we are at the boundary of the cube
//			from[min-1] = Math.min(min, from[min-1]);
//			to[min-1] = Math.max(max, to[min-1]);
//		}
//		if (max < n-1) {
//			// update gap on the right of the monophyletic constraint
//			// if max == n-1 there is no gap to update since we are at the boundary of the cube
//			from[max] = Math.min(min, from[max]);
//			to[max] = Math.max(max, to[max]);
//		}
//	}
//
//
//
//	/** extract clades from Newick string, and add constraints for all internal nodes
//     ** (except the root if it contains all taxa). This code populates taxonIDList.
//     **
//	 **/
//	private List<List<Integer>> parse(String newick, String [] taxaList) {
//		List<List<Integer>> taxonIDList = new ArrayList<>();
//        assert taxonIDList.size() == 0;
//
//		// get rid of initial and trailing spaces
//		newick = newick.trim();
//		// remove comments
//		newick = newick.replaceAll(".*\\[[^\\]]*\\].*", "");
//		// remove branch lengths
//		newick = newick.replaceAll(":[^,\\(\\)]*", "");
//		Pattern pattern = Pattern.compile("\\(([^\\(\\)]*)\\)");
//		Matcher m = pattern.matcher(newick);
//		
//		String prev = "";
//		while (m.find()) {
//			String group = m.group();
//			String [] taxa = group.substring(1,group.length()-1).split(",");
//			List<Integer> list = new ArrayList<Integer>();
//			for (String taxon : taxa) {
//				taxon = taxon.trim();
//				if (taxon.length() > 0) {
//					int i = indexOf(taxon, taxaList);
//					if (i == -1 && (taxon.startsWith("'") || taxon.startsWith("\""))) {
//						i = indexOf(taxon.substring(1, taxon.length() - 1), taxaList);
//					}
//                    if (i == -1) {
//                        throw new RuntimeException("Cannot find taxon " + taxon + "  in taxon list");
//                    }
//					list.add(i);
//				}
//			}
//			// 1. only add when it is not the complete taxonset
//			// 2. make sure it is not equal to previous set -- happens with one-node-branches
//			if (list.size() < tree.getLeafNodeCount() && !group.equals(prev)) {
//				taxonIDList.add(list);
//				//Log.trace.println("Constraining " + group);// + " " + Arrays.toString(list.toArray()));			
//			}
//			newick = newick.replaceFirst("\\(([^\\(\\)]*)\\)", ",$1,");
//			newick = newick.replaceAll("([\\(,]),", "$1");
//			newick = newick.replaceAll(",\\)", ")");
//			m = pattern.matcher(newick);
//			prev = group;
//		}
//		return taxonIDList;
//	}
//
//
//    private int indexOf(String taxon, String [] taxaList) {
//		for (int k = 0; k < taxaList.length; k++) {
//			if (taxon.equals(taxaList[k])) {
//				return k;
//			}
//		}
//		Log.warning("Could not find taxon " + taxon + "\nPerhaps a typo in the taxon name?");
//		return -1;
//	}
//
//	private void calcRelHeights(double [] heights) {
//		for (int i = 0; i < heights.length; i++) {
//			double max = 0;
//			for (int j = from[i]; j < to[i]; j++) {
//				if (heights[j] > max) {
//					max = heights[j];
//				}
//			}
//			if (max > 0) {
//				relheights[i] = heights[i] / max;
//			} else {
//				relheights[i] = heights[i];
//			}
//		}			
//	}
//
//	
//	@Override
//	public double[] forward(double[] x) {
//		double [] values = new double[x.length];
//		for (int i = 0; i < x.length; i++) {
//			values[i] = Math.log(x[i]);
//		}
//		cube2Tree(order, values, tree);
//		calcRelHeights(values);
//		System.arraycopy(relheights, 0, values, 0, values.length);
//		return values;
//	}
//
//	@Override
//	public double[] backward(double[] z) {
//		double [] values = new double[z.length];
//		for (int i = 0; i < z.length; i++) {
//			values[i] = Math.exp(z[i]);
//		}
//		calcHeights(values);
//		cube2Tree(order, values, tree);
//		return values;
//	}
//
//
//	
//	private void calcHeights(double[] relheights) {
//		final int n = relheights.length;
//		double [] heights = new double[n];
//		for (int i = 0; i < n; i++) {
//			setHeight(heights, relheights, i);
//		}
//		System.arraycopy(heights, 0, relheights, 0, n);
//		
//	}
//
//	//int depth = 0;
//	// set height of gap k based on relative heights and from-to constraints
//	// if not all necessary gap heights are available, set those first recursively
//	// until we have all gap heights to determine the height of the current gap
//	private void setHeight(double[] heights, double[] relheights2, int k) {
//		//depth++;
//		if (heights[k] > 0) {
//			// height is already set
//			//depth--;
//			return;
//		}
//		// first get heights of the highest clades left and right and store in max
//		// if any of these are monophyletic, 
//		// -1 < from[k] < i for left clade and 
//		//  i < to[k] < n for right clade
//		double max = 0;
//		for (int i = from[k]; i < to[k]; i++) {
//			if (i != k) {
//				if (heights[i] == 0) {
//					setHeight(heights, relheights2, i);
//				}
//				max = Math.max(heights[i], max);
//			}
//		}
//		if (max == 0) {
//			// no clade constrain this gap in the cube, so just use relative height as height
//			heights[k] = relheights2[k];
//		} else {
//			// one or two clades are constrained, so interpret relative height wrt to the highest clade height
//			heights[k] = relheights2[k] * max;
//		}
//		//depth--;
//	}
//
//	@Override
//	public double forward(double x) {
//		throw new RuntimeException("Cannot use single dimension call on StickBreaking transform");
//	}
//
//	@Override
//	public double backward(double z) {
//		throw new RuntimeException("Cannot use single dimension call on StickBreaking transform");
//	}
//
//	@Override
//	public double jacobian_det(double x) {
//		throw new RuntimeException("Cannot use single dimension call on StickBreaking transform");
//	}
//
}
