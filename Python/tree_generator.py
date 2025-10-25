import sys
import os
import toytree
import toyplot.svg

def generate_tree(tree_file):
    """
    Summary:
        Generates a .svg file from a newick string phylogenetic tree.
    Parameters:
        tree_file: Path to the file containing the Newick tree string.
    """
    
    # Ensure the input file path exists
    if not os.path.isfile(tree_file):
        raise FileNotFoundError(f"{tree_file} not found!!!")
    
    # Read the input file
    with open(tree_file, "r") as file:
        newick_file = file.read().strip()
        
    # Create a toytree object in order to manipulate the newick tree
    tree = toytree.tree(newick_file)
    
    # Set plot dimensions based on the number of taxa to help tree visualization
    taxa = len(tree.get_tip_labels())
    width = max(2000, taxa * 50)
    heigth = max(1000, taxa * 30)
    
    # Get support values
    support = tree.get_node_data("support")
    labels = [f"{round(float(s), 2)}" if s else "" for s in support]
    
    # Define the output file path
    pwd = os.getcwd()
    outfile = os.path.join(pwd, "phylo_tree.svg") 
    
    # Create the .svg file and save it
    canvas, axes, mark = tree.draw(node_labels=labels, node_sizes=30, width=width, height=heigth, tree_style="n")
    toyplot.svg.render(canvas, outfile)

# Just ensuring the code is only executed when the script is run as a standalone program.
if __name__ == "__main__":
    if len(sys.argv)!=2:
        print("Usage: python3 {path/to/script.py} {path/to/support/tree/file}")
    else:
        # Generate SVG file
        generate_tree(sys.argv[1])
        
    
        
    