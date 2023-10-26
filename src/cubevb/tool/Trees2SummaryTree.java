package cubevb.tool;


import java.io.IOException;
import java.io.PrintStream;
import java.util.*;


import beastfx.app.tools.Application;
import beastfx.app.treeannotator.TreeAnnotator;
import beastfx.app.treeannotator.TreeAnnotator.MemoryFriendlyTreeSet;
import beastfx.app.treeannotator.TreeAnnotator.TreeSet;
import beastfx.app.treeannotator.services.TopologySettingService;
import beastfx.app.util.OutFile;
import beastfx.app.util.TreeFile;
import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.core.Log;
import beast.base.evolution.alignment.TaxonSet;
import beast.base.evolution.distance.Distance;
import beast.base.evolution.alignment.Alignment;
import beast.base.evolution.alignment.Sequence;
import beast.base.evolution.alignment.Taxon;
import beast.base.evolution.tree.ClusterTree;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;

@Description("Convert posterior tree set to distance matrix and return UPGMA tree on the matrix")
public class Trees2SummaryTree extends beast.base.inference.Runnable implements TopologySettingService {
	final public Input<TreeFile> treesInput = new Input<>("trees", "beast.trees on which this operation is performed",
			new TreeFile("[[none]]"));
	final public Input<Integer> burnInPercentageInput = new Input<>("burnin",
			"percentage of trees to used as burn-in (and will be ignored). Set >= 100 if only a cube is required", 10);
	final public Input<OutFile> outputInput = new Input<>("out", "output file. Fails if not specified", Validate.REQUIRED);
	
	public Trees2SummaryTree() {
	}

	@Override
	public void initAndValidate() {
		// nothing to do
	}

	
	@Override
	public void run() throws Exception {
		MemoryFriendlyTreeSet treeSet = new TreeAnnotator().new MemoryFriendlyTreeSet(treesInput.get().getAbsolutePath(), burnInPercentageInput.get());

		Tree clusterTree = setTopology(treeSet, Log.warning, null);
		
		PrintStream out = System.out;
		if (outputInput.get() != null && !outputInput.get().getName().equals("[[none]]")) {
			out = new PrintStream(outputInput.get());
		}
		clusterTree.init(out);
		out.println();
		clusterTree.log(0, out);
		out.println();
		clusterTree.close(out);
		if (outputInput.get() != null && !outputInput.get().getName().equals("[[none]]")) {
			out.close();
		}

		
		Log.warning("Done");
	}
		
	

	/**
	 * calculate the distance between leafs in a consensus tree and update the
	 * distance matrix weighted with the relative frequency of the tree
	 * 
	 * @param node:
	 *            current node
	 * @param nOrder:
	 *            mapping of node label to [0,...,NrOfLeafs-1]
	 * @param fDistMatrix:
	 *            distance matrix to be updated
	 * @param iLabel:
	 *            used to report set of leafs in sub tree below node
	 * @param fLength:
	 *            used to report set of lengths from current node to leafs in
	 *            iLabel
	 */
	private void calcDistance(Node node, int index, double[][] fDistMatrix, List<Integer> iLabel,
			List<Double> fLength) {
		if (node == null) {
			return;
		}
		if (node.isLeaf()) {
			// iLabel.add(nOrder[node.m_iLabel]);
			iLabel.add(node.getNr());
			// fLength.add(node.m_fLength);
			fLength.add(node.getLength() + node.getHeight());
		} else {
			iLabel.clear();
			List<Integer> iLeft = new ArrayList<>();
			List<Integer> iRight = new ArrayList<>();
			List<Double> fLeft = new ArrayList<>();
			List<Double> fRight = new ArrayList<>();
			calcDistance(node.getChild(0), index, fDistMatrix, iLeft, fLeft);
			calcDistance(node.getChild(1), index, fDistMatrix, iRight, fRight);
			for (int i = 0; i < iLeft.size(); i++) {
				int i1 = iLeft.get(i);
				double f1 = fLeft.get(i);
				for (int j = 0; j < iRight.size(); j++) {
					int i2 = iRight.get(j);
					double f2 = f1 + fRight.get(j);
					f2 = Math.log(f2);
					fDistMatrix[i1][i2] += f2;
					fDistMatrix[i2][i1] += f2;
				}
			}
			for (int i = 0; i < fLeft.size(); i++) {
				iLabel.add(iLeft.get(i));
				// fLength.add(fLeft.get(i) + node.m_fLength);
				fLength.add(fLeft.get(i) + node.getLength());
				// fLength.add(fLeft.get(i) + 1.0f);
			}
			for (int i = 0; i < fRight.size(); i++) {
				iLabel.add(iRight.get(i));
				// fLength.add(fRight.get(i) + node.m_fLength);
				fLength.add(fRight.get(i) + node.getLength());
				// fLength.add(fRight.get(i) + 1.0f);
			}
		}
	} // calcDistance




	@Override
	public Tree setTopology(TreeSet treeSet, PrintStream progressStream, TreeAnnotator annotator) throws IOException {
        treeSet.reset();
        Tree tree = treeSet.next();
        treeSet.reset();
    	String[] taxa = tree.getTaxaNames();
        int n = taxa.length;
        
        progressStream.println("0              25             50             75            100");
        progressStream.println("|--------------|--------------|--------------|--------------|");
        
		double[][] distance = new double[n][n];
		int k = 0;
		int numTrees = treeSet.totalTrees - treeSet.burninCount;
		int percentageDone = 0;
        while (treeSet.hasNext()) {
        	tree = treeSet.next();
    		calcDistance(tree.getRoot(), k, distance, new ArrayList<Integer>(), new ArrayList<Double>());
			while ((62*k) /numTrees > percentageDone) {
				progressStream.print("*");
				progressStream.flush();
				percentageDone++;
			}
    		k++;
        }
        progressStream.print("\n" + k + " trees fetched ");
		
		TaxonSet taxonset = new TaxonSet();
		List<Sequence> sequences = new ArrayList<>();
		for (String s : taxa) {
			taxonset.taxonsetInput.get().add(new Taxon(s));
			sequences.add(new Sequence(s, "?"));
		}
		Alignment data = new Alignment(sequences, "nucleotide");

		
		// normalise matrix
		for (int i = 0; i < n; i++) {
			double [] d = distance[i];
			for (int j = 0; j < n; j++) {
				d[j] /= k;
				d[j] = Math.exp(d[j]);
			}
		}
		
		Distance distance_ = new Distance() {
			@Override
			public double pairwiseDistance(int taxon1, int taxon2) {
				return distance[taxon1][taxon2];
			}
		};
		
		
		Tree clusterTree = new ClusterTree();
		clusterTree.initByName("clusterType", "single", "taxa", data, "distance", distance_);
		for (int i = 0; i < clusterTree.getLeafNodeCount(); i++) {
			clusterTree.getNode(i).setHeight(tree.getNode(i).getHeight());
		}
		return clusterTree;
	}

	@Override
	public String getDescription() {
		return "Matrix on average distance matrix from trees";
	}

	@Override
	public String getServiceName() {
		return "matrix";
	}

	public static void main(String[] args) throws Exception {
		new Application(new Trees2SummaryTree(), "Trees2SummaryTree", args);
	}


}
