#!/usr/bin/env python

import crag
import argparse

parser = argparse.ArgumentParser()

# set argparse
parser.add_argument('--gfile', type=str, default="",
		    help='gen file needed')
parser.add_argument('--sfile', type=str, default="",
		    help='sample file needed')
parser.add_argument('--ofile', type=str, default="",
		    help='output file stem, e.g. <ofile>.out')
parser.add_argument('--chr', type=str,
		    help='chromosome number to use in output file')
parser.add_argument('--obs', type=int, default=1,
		    help='which obs indicator is the competing risk (1 is default, do not specify 0 as it is reserved for censored observations)')
parser.add_argument('--t_pheno', type=str, default="",
		    help='time to event in SAMPLE file')
parser.add_argument('--et_pheno', type=str, default="",
		    help='event type indicator in SAMPLE file')
parser.add_argument('--covs', type=str, default="",
		    help='comma delimited list of covariates in SAMPLE file')
parser.add_argument('--int', type=str, default="",
		    help='covariate to adjust for SNP covariate interaction, must be from SAMPLE file')
parser.add_argument('--crtype', type=int, default=0,
		    help='use 0 for subdistribution competing risks (default), 1 for cause-specific competing risks and 2 for conventional Cox PH survival model')
parser.add_argument('--r', type=int, default=0,
		    help='use 1 to perform analysis with R package cmprsk or 0 to perform analysis in Python package lifelines')
parser.add_argument('--verbose', type=int, default=0,
		    help='use 1 to print results of covariates tested at each SNP, includes interaction terms of specified')
parser.add_argument('--thread', type=int, default=1,
		    help='specify number of threads to use (please consider the hardware limitations of your computer/cluster)')
args=parser.parse_args()

if args.r == 1:
	import rpy2.robjects as robjects
	r = robjects.r
	from rpy2.robjects.packages import importr
	from rpy2.robjects import r, pandas2ri
	pandas2ri.activate()
	#import R's "base" package
	base = importr('base')
	#import R's "utils" package
	utils = importr('utils')
	cmprsk = importr('cmprsk')
	robjects.r['options'](warn=-1)
	from rpy2.rinterface import RRuntimeError
	#warnings.filterwarnings("ignore", category=RRuntimeError) # Sort out why this generates error!
	robjects.r('sink("/dev/null")')

if __name__ == "__main__":
	print("##############################")
	print("#      WELCOME TO CRAG       #")
	print("##############################\n")

	# Reduce sample to needed covs

	sample = pd.read_csv(args.sfile,sep=" ") 
	list=[args.t_pheno,args.et_pheno]
	if args.covs!="":
		list.extend(args.covs.split(','))
	if args.int!="":
		list.extend(args.int.split(','))
	sub=sample[list]
	sub=sub.drop([0])
	sub[list] = sub[list].apply(pd.to_numeric)
	sub['bin'] = np.where(sub[args.et_pheno] == args.obs, 1, 0)

	if args.crtype==0:

		tmax=sub.loc[sub['bin']==1][args.t_pheno].max()
		sub[args.t_pheno] = np.where((sub[args.et_pheno] != args.obs) & (sub[args.et_pheno] != 0), tmax, sub[args.t_pheno])

	if args.crtype==2:

		sub['bin'] = np.where(sub[args.et_pheno] > 0, 1, 0)

	sub=sub.drop(columns=[args.et_pheno])



	# Check for low variance
	if args.covs:

		print("##############################")
		print("#    CLINICAL ANALYSIS...    #")
		print("##############################\n")
		low_var = pd.DataFrame(sub.var(0) < 10e-5)
		print("\nResults of the low variance (<10e-5) test... (if True these are removed to avoid convergence issues)")
		if args.int:

			print(low_var.loc[args.covs.split(',') + args.int.split(',')])

		if not args.int:

			print(low_var.loc[args.covs.split(',')])

		lvf=low_var[low_var.iloc[:,0] == True]

		if not lvf.empty:

			sub = sub.drop(lvf.index,axis = 1)

		print("\nResults of pre-GWAS multivariable analysis...")

		cph.fit(sub,duration_col=args.t_pheno,event_col='bin')
		cph.print_summary()

	# Open gen file and perform calculation every line
	colNames = ('chr','snp','test','bp','ea','nea','p','info','maf','hwe','beta','se')
	#CHR SNP BP P INFO MAF HWE BETA SE CONV
	print("##############################")
	print("#   BEGIN GWAS ANALYSIS...   #")
	print("##############################\n")
	num_workers = args.thread

	manager = Manager()
	results = manager.list()
	work = manager.Queue(num_workers)

	#start for workers
	pool = []
	for i in range(num_workers):
		p = Process(target=do_work, args=(work, results,sub))
		p.start()
		pool.append(p)
	i=1
	if args.gfile.endswith('gz'):

		with gzip.open(args.gfile,'rt') as f:
			iters = itertools.chain(f, (None,)*num_workers)
			for num_and_line in enumerate(iters):
				gwasprogress()
				i+=1

	else:
		with open(args.gfile,'rt') as f:
			iters = itertools.chain(f, (None,)*num_workers)
			for num_and_line in enumerate(iters):
				gwasprogress()
				i+=1
	for p in pool:
		p.join()

	f_results=pd.concat(results)
	final=f_results.ix[:,colNames]
	final.to_csv(args.ofile + ".out",sep=" ",index=False)

	print("\n\n##############################")
	print("#      END OF ANALYSIS       #")
	print("##############################\n")
