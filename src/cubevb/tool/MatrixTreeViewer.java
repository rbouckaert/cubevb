package cubevb.tool;




import java.awt.BorderLayout;
import java.awt.GridLayout;
import java.util.*;

import javax.swing.JCheckBox;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JSlider;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;

import beast.base.core.Description;
import beast.base.evolution.alignment.Alignment;
import beast.base.evolution.alignment.Sequence;
import beast.base.evolution.alignment.Taxon;
import beast.base.evolution.alignment.TaxonSet;
import beast.base.evolution.distance.Distance;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;

@Description("Represent tree as pairwise symmetric distance matrix, use single link clustering to resolve tree")
public class MatrixTreeViewer extends CubeTreeViewer {
	private static final long serialVersionUID = 1L;
	ClusterTreeViewerControlPanel controlPanel;
	
	public MatrixTreeViewer() {
		super(null);
	}
	
	@Override
	protected void initialise(String startValues) {
		N = 3;
		setLayout(new BorderLayout());
		
		controlPanel = new ClusterTreeViewerControlPanel();
		add(controlPanel, BorderLayout.WEST);
		treePanel = new TreePanel();
		add(treePanel, BorderLayout.CENTER);
	}
	

	boolean linkDistances = false;
	
	public class ClusterTreeViewerControlPanel extends JPanel {
		private static final long serialVersionUID = 1L;
		
		JSlider [][] slider;
		
		public ClusterTreeViewerControlPanel() {
			setLayout(new GridLayout(N+1,N));
			
			slider = new JSlider[N+1][N+1];
			for (int i = 1; i <= N; i++) {
				for (int j = 0; j < i; j++) {
					slider[i][j] = new JSlider(0, 100, 100); // Randomizer.nextInt(100));
					slider[j][i] = slider[i][j];
					slider[i][j].setOrientation(JSlider.VERTICAL);
					slider[i][j].setName(i + "," + j);
					slider[i][j].addChangeListener(new ChangeListener() {				
						@Override
						public void stateChanged(ChangeEvent e) {
							String name = ((JSlider) e.getSource()).getName();
							String [] str = name.split(",");
							int i = Integer.parseInt(str[0]);
							int j = Integer.parseInt(str[1]);
							updateClusterTree(i, j);
						}
					});
					add(slider[i][j]);
				}
				for (int j = i; j < N; j++) {
					add(new JLabel("-"));
				}
			}
			JCheckBox linkbox = new JCheckBox("<html>link<br>distances</html>");
			linkbox.setSelected(linkDistances);
			linkbox.addActionListener(e -> {
				updateClusterTree(0, 1);				
				linkDistances = linkbox.isSelected();
				setSliderValues();
				updateClusterTree(0, 1);
			});
			add(linkbox);
		}

		
		private void updateClusterTree(int x, int y) {
			if (linkDistances) {
				double target = slider[x][y].getValue();
				target *= 0.01;
				Set<Node> path = new HashSet<>();
				Node node = commonAncestor(nodes[x], nodes[y]);
				if (node.getHeight() == target) {
				} else if (node.getHeight() > target) {
					// sliding down
					if (target < node.getLeft().getHeight() || 
							target < node.getRight().getHeight()) {
						Node c = node.getLeft().getAllLeafNodes().contains(nodes[x]) ? 
								node.getLeft() : node.getRight();
						if (!node.isRoot()) {
							Node p = node.getParent();
							p.removeChild(node);
							p.addChild(c);
						}
						Node c2 = c.getLeft().getAllLeafNodes().contains(nodes[x]) ? 
								c.getLeft() : c.getRight();
						node.removeChild(c);
						node.addChild(c2);
					}
				} else {
					// sliding up
					if (!node.isRoot() && target > node.getParent().getHeight()) {
						Node p = node.getParent();
						p.removeChild(node);
						node.setParent(null);
						p.addChild(node.getLeft());
						node.removeChild(node.getLeft());
						node.addChild(p);
						if (node.isRoot()) {
							root = node;
						}
					}
				}
				node.setHeight(target);

				setSliderValues();
				
			} else {
			
				double [][] matrix = new double [N+1][N+1];
				for (int i = 0; i < N+1; i++) {
					for (int j = 0; j < i; j++) {
						matrix[i][j] = slider[i][j].getValue();
						matrix[j][i] = matrix[i][j];
					}
				}
	
				
				Tree clusterTree = matrix2Tree(matrix, "single");
				
				root = clusterTree.getRoot();
				nodes = clusterTree.getNodesAsArray();
				clusterTree.scale(0.02);
				
				System.out.print(root.toNewick());
			}
			treePanel.repaint();
			
		}


		private void setSliderValues() {
			// set slider values
			for (int i = 1; i <= N; i++) {
				for (int j = 0; j < i; j++) {
					Node node = commonAncestor(nodes[i], nodes[j]);
					slider[i][j].setValue((int) (100 * node.getHeight()));
				}
			}
		}


		private Node commonAncestor(Node node, Node node2) {
			if (node == node2) {
				return node;
			}
			if (node.isRoot()) {
				return node;
			}
			if (node2.isRoot()) {
				return node2;
			}
			if (node == node2.getParent()) {
				return node;
			}
			if (node2 == node.getParent()) {
				return node2;
			}
			if (node.getParent().getHeight() > node2.getParent().getHeight()) {
				return commonAncestor(node, node2.getParent());
			}
			return commonAncestor(node.getParent(), node2);
		}

	}

	public static Tree matrix2Tree(double[][] matrix, String clustertype) {
		int N = matrix.length;
		Distance distance = new Distance() {			
			@Override
			public double pairwiseDistance(int taxon1, int taxon2) {
				return matrix[taxon1][taxon2];
			}
		};
		
		TaxonSet taxonset = new TaxonSet();
		List<Sequence> sequences = new ArrayList<>();
		for (int i = 0; i < N; i++) {
			taxonset.taxonsetInput.get().add(new Taxon("taxon" + i));
			sequences.add(new Sequence("taxon" + i, "?"));
		}
		taxonset.initAndValidate();
		Alignment data = new Alignment(sequences, "nucleotide");
		
//		Matrix2Tree matrixTree = new Matrix2Tree();
//		matrixTree.initByName( 
//				"taxonset", taxonset,
//				"distance", distance);
//		return matrixTree;

		
		ClusterTree clusterTree = new ClusterTree();
		clusterTree.initByName("clusterType", clustertype, 
				"taxa", data,
				"distance", distance);
		return clusterTree;
	}
		
	
	
	public static void main(String[] args) {
		JFrame frame = new JFrame();
		frame.setSize(1024, 768);
		MatrixTreeViewer viewer = new MatrixTreeViewer();
		frame.add(viewer);
		frame.setVisible(true);
		frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
	}

}
