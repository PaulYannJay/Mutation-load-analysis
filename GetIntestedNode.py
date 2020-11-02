## From a newick tree and a list of sample, find the TMRCA of these samples
from ete3 import PhyloTree
from optparse import OptionParser

parser = OptionParser()
parser.add_option("-f", "--file", dest="filename",
	help="File containing newick tree", metavar="FILE")
parser.add_option("-s", "--species", dest="speciesList",
	help="file containing list of wanted species, one per line", metavar="FILE")
(options, args) = parser.parse_args()

t= PhyloTree(options.filename, format=1)

with open(options.speciesList) as f:
	liste = f.read().splitlines()

pp=t.get_common_ancestor(liste)
print pp.name
