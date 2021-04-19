import os
import pickle
import timeit

import numpy as np
import pandas as pd

from pyscenic.utils import load_motifs
from pyscenic.transform import df2regulons

if __name__ == "__main__":
	## set up basic parameters
	repeats = 10 # how often to repeat the GRN inference procedure?
	no_workers_grn = 40 
	no_workers_pruning = 6 # there might be a different value required here as 
						# as this seems to crash with large no. of workers
			# (e.g. 32 workers results in "'std::system_error', what(): resource unavailable"
	np.random.seed(42)
	seeds = np.random.randint(low=1, high=1000, size=repeats)
	
	## set up dataframe for basic profiling
	#times = pd.DataFrame(np.zeros((repeats, 2)), 
	#                     index=range(1,(repeats+1)), columns=["grn-inference", "regulon-pruning"])
	#times = pd.read_csv("data/day61/results/repeats-whole-data/times.csv")	
	
	## prepare filenames
	cm_fname = "Table_HIO_mes_count_data_t.csv"
	tf_fname = "/nas/groups/treutlein/USERS/Qianhui_Yu/Annotation/pySCENIC/2020_Nat_protocol/hs_hgnc_tfs.txt"
	ranking_dbs = list(map(lambda fn: os.path.join("/nas/groups/treutlein/USERS/Qianhui_Yu/Annotation/pySCENIC/2020_Nat_protocol/", fn),
						   ['hg19-tss-centered-5kb-7species.mc9nr.feather']))
	dbs_param = ' '.join(ranking_dbs)
	mtf_annotations = "/nas/groups/treutlein/USERS/Qianhui_Yu/Annotation/pySCENIC/2020_Nat_protocol/motifs-v9-nr.hgnc-m0.001-o0.0.tbl"
	for repeat in range(1,(repeats+1)):
		print("\n##### repeat {} #####".format(repeat))
		
		results_dir = "repeat0{}".format(repeat)
		if os.path.exists(results_dir):
			pass
		else: 
			os.makedirs(results_dir)
		
		adj_fname = os.path.join(results_dir, "adjacency.tsv")
		mtf_fname = os.path.join(results_dir, "motifs.csv")
		reg_fname = os.path.join(results_dir, "regulons.p")
		#aucell_mtx_fname = os.path.join(results_dir, 'aucell-mtx.csv')
		#ranking_dbs = list(map(lambda fn: os.path.join("/nas/groups/treutlein/USERS/Qianhui_Yu/Annotation/pySCENIC/2020_Nat_protocol/", fn),
		#					   ['hg19-500bp-upstream-7species.mc9nr.feather',
		#						'hg19-tss-centered-10kb-7species.mc9nr.feather',
		#						'hg19-tss-centered-5kb-7species.mc9nr.feather']))
		seed = seeds[repeat-1]
		 
		## step 1: GRN inference
		print("### GRN inference ##\n")
		#start_time = timeit.default_timer()
		os.system("python arboreto_with_multiprocessing.py {cm} {tf_file}".format(cm=cm_fname, 
																				  tf_file=tf_fname) +
			" -o {adj} --num_workers {nw} --seed {s} --m grnboost2".format(adj=adj_fname, 
																		   nw=no_workers_grn,
																		   s=seed))
		#end_time = timeit.default_timer()
		#times.loc[repeat, "grn-inference"] = round((end_time - start_time)/60, 2)
		
		## step 2 & 3: create regulons
		print("### regulon pruning #\n")
		#start_time = timeit.default_timer()
		os.system("pyscenic ctx {adj} {dbp} --annotations_fname {mtf_an}".format(adj=adj_fname,
																				 dbp=dbs_param,
																				 mtf_an=mtf_annotations) +
			" --expression_mtx_fname {cm} --output {mtf} --num_workers {nw} --mode {modep}".format(cm=cm_fname,
																				   mtf=mtf_fname,
																				   nw=no_workers_pruning,
																				   modep="custom_multiprocessing"))
		df_motifs = load_motifs(mtf_fname)
		regulons = df2regulons(df_motifs)
		with open(reg_fname, 'wb') as f:
			pickle.dump(regulons, f)
		
		name = np.concatenate([[r.name]*len(r.gene2weight) for r in regulons])
		score = np.concatenate([[r.score]*len(r.gene2weight) for r in regulons])
		weight = np.concatenate([list(dict(r.gene2weight).values())
		    for r in regulons])
		genes = np.concatenate([list(dict(r.gene2weight).keys())
			for r in regulons])
		regulon_df = pd.DataFrame({
			'name':name,
			'weight':weight,
			'score':score,
			'genes':genes
		})
		
		# Save the enriched motifs to disk.
		regulons_df_out = os.path.join(results_dir, 'regulons.csv')
		regulon_df.to_csv(regulons_df_out, index=False)
	
		## step 4: cellular enrichment
		#os.system("pyscenic aucell "+cm_fname+" "+mtf_fname+" --output "+aucell_mtx_fname+" --num_workers "+str(no_workers_pruning))

		#end_time = timeit.default_timer()
		#times.loc[repeat, "regulon-pruning"] = round((end_time - start_time)/60, 2)
	
	#times.to_csv("data/day61/results/repeats-whole-data/times.csv")
