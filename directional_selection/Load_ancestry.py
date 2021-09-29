
import msprime
import pyslim
import numpy as np


# Load the resulting .trees file

ts = pyslim.load("Slim_output.tree")

### Only keeping females


female_nodes = []
for ind in ts.individuals():
  if ind.metadata['sex'] == pyslim.INDIVIDUAL_TYPE_FEMALE:
    female_nodes.extend(ind.nodes)

female_ts = ts.simplify(female_nodes, keep_input_roots=True)


# Assess true local ancestry
breaks = np.zeros(female_ts.num_trees + 1)
ancestry = np.zeros(female_ts.num_trees + 1)
    
for tree in female_ts.trees(sample_counts=True):
	subpop_sum, subpop_weights = 0, 0
	for root in tree.roots:
		leaves_count = tree.num_samples(root) - 1  # subtract one for the root, which is a sample
		subpop_sum += tree.population(root) * leaves_count
		subpop_weights += leaves_count
	breaks[tree.index] = tree.interval[0]
	ancestry[tree.index] = subpop_sum / subpop_weights
	
breaks[-1] = ts.sequence_length
ancestry[-1] = ancestry[-2]



anc = np.matrix(ancestry)
with open("ancestry.txt", "wb") as f:
   for line in anc:
      np.savetxt(f, line, fmt="%.4f", delimiter="\n")
    
    
mat = np.matrix(breaks)
with open("breaks.txt", "wb") as f:
    for line in mat:
        np.savetxt(f, line, delimiter="\n", fmt="%e")
