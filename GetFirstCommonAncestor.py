#From a Newick tree and two list of samples, extract the first coalescent event between samples of this two groups
from ete3 import PhyloTree
from optparse import OptionParser

parser = OptionParser()
parser.add_option("-f", "--file", dest="filename",
	help="File containing newick tree", metavar="FILE")
parser.add_option("-s", "--species1", dest="species1List",
	help="file containing list of samples, group1, one per line", metavar="FILE")
parser.add_option("-p", "--species2", dest="species2List",
	help="file containing list of samples, group2, one per line", metavar="FILE")
(options, args) = parser.parse_args()

t= PhyloTree(options.filename, format=1)

with open(options.species1List) as f:
	liste1 = f.read().splitlines()

with open(options.species2List) as f:
	liste2 = f.read().splitlines()

for s1 in liste1:
	for s2 in liste2:
		pp=t.get_common_ancestor(s1,s2)
		print pp.name
