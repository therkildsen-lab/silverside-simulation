import pyslim, msprime
import numpy as np
import random
import time
import re
import ast
import sys
import os

#dat = pyslim.load("1672548051455_SilverSide_Inversion.tree")
#dat = pyslim.load("1612388845805_SilverSide_Inversion.tree")
dat = pyslim.load("1614460993787_SilverSide_Inversion.tree")

print(f"The tree sequence has {dat.num_trees} trees on a genome of length {dat.sequence_length},"
      f" {dat.num_individuals} individuals, {dat.num_samples} 'sample' genomes,"
      f" and {dat.num_mutations} mutations.")

mig_mat = [[0.0, 0.0, 0.0], [0.0, 0.0, 0.001], [0.0, 0.001, 0.0]]

positions = []
rates = []
with open('recomb_rates2.tsv', 'r') as file:
   header = file.readline().strip().split("\t")
   assert(header[0] == "end_position" and header[1] == "rate")
   for line in file:
      components = line.split("\t")
      positions.append(float(components[0]))
      rates.append(float(components[1]))

# step 1
positions.insert(0, 0)
# step 2
rates.append(0.0)
# step 3
#positions[-1] += 1

recomb_map = msprime.RecombinationMap(positions, rates)

start = time.time()
recap = dat.recapitate(recombination_map=recomb_map, migration_matrix=mig_mat, Ne=5000)
end = time.time()
print("Time it took to recapitate:", end - start)

print(f"The tree sequence has {dat.num_trees} trees,"
      f" and {dat.num_mutations} mutations.")

#recap.site(recap.mutation(0).site)
start = time.time()
mutated = pyslim.SlimTreeSequence(msprime.mutate(recap, rate=1e-7, keep=False))
end = time.time()
print("Time it took to add mutations:", end - start)

print(f"The tree sequence now has {mutated.num_trees} trees,"
      f" and {mutated.num_mutations} mutations.")

mutated.dump("1614460993787_SilverSide_Inversion_recap_mut2.tree")

indA = mutated.individuals_alive_at(0)

indivlist = []
indivnames = []

with open("spatial_sim_individuals_1614460993787.txt", "w") as indfile:
	indfile.writelines("\t".join(["vcf_label", "tskit_id", "slim_id"] + ["birth_time_ago", "age", "x", "y"]) + "\n")
	for i in indA:
		indivlist.append(i)
		ind = mutated.individual(i)
		vcf_label = f"tsk_{ind.id}"
		indivnames.append(vcf_label)
		data = [vcf_label, str(ind.id), str(ind.metadata["pedigree_id"]), str(ind.time), str(ind.metadata["age"]), str(ind.location[0]), str(ind.location[1])]
		indfile.writelines("\t".join(data) + "\n")

with open("silverside_sim_1614460993787_2.vcf", "w") as vcffile:
  mutated.write_vcf(vcffile, individuals=indivlist, individual_names=indivnames)