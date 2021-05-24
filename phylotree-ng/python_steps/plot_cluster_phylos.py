import argparse
import glob

import matplotlib
import toyplot
import toyplot.pdf
import toytree

matplotlib.use('Agg')

def main(outgroup: str):
    tree_file_paths = glob.glob('ska_outputs/*/*.treefile')

    for tree_file in tree_file_paths:
        file = open(tree_file, 'r')
        lines = file.readlines()
        newick = lines[0].strip()

        # set the outgroup to first sample unless otherwise specified
        # TODO: create a more robust method for picking outgroup...
        # ... the current method of specifying it will break / error in cases where there
        # ... is > 1 cluster
        wildcard_value = newick.split('.')[0].split('(')[1]
        if(outgroup == ""):
            wildcard_value = newick.split('.')[0].split('(')[1]
        else:
            wildcard_value = outgroup

        tre0 = toytree.tree(newick, tree_format=0)
        rtre = tre0.root(wildcard=wildcard_value)

        style = {
            "tip_labels_align": False,
            "tip_labels_style": {
                "font-size": "9px"
            },
        }
        canvas, axes, _ = rtre.draw(tip_labels_colors='indigo', **style);
        axes.show = True
        axes.x.ticks.show = True
        axes.y.ticks.show = False
        print("about to create plot")
        print('phylo_plot_outputs/' + tree_file.split('/')[1] + '.treeplot.pdf')
        toyplot.pdf.render(canvas, 'phylo_plot_outputs/' + tree_file.split('/')[1] + '.treeplot.pdf')
        #plt.savefig('phylo_plot_outputs/' + tree_file.split('/')[1] + '.treeplot.png', bbox_inches='tight')

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--outgroup")
    args = parser.parse_args()
    main(args.outgroup)
