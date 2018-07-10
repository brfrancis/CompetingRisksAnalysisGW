#!/usr/bin/env python

from crag import sample_read, crtype_sample, low_var_test, QC, pycrgwas, output, oresults
import argparse
import threading
from multiprocessing import Process, Manager, Pool
import itertools
import time
start=time.time()
import warnings
import gzip
import sys

def do_work(in_queue, out_list, sub):
	while True:
		item = in_queue.get()
		line_no, line = item
	# exit signal
		if line == None:
			return

	# Using lifelines work for cause-specific
		gp, info, maf, hwe, gl = QC(line, line_no, nchr)
		df3_1 = pycrgwas(gp, sub, x, t_pheno)
		result = output(df3_1, info, maf, hwe, gl, verbose)
		# output
		out_list.append(result)

def gwasprogress():
	sys.stdout.write('\r')
	progress = "#" + str(i-thread) + " SNPS ANALYSED @ {0:0.01f} SECS PER SNP #".format((time.time()-start)/i)
	#sys.stdout.write(progress.center(os.get_terminal_size().columns))
	sys.stdout.write(progress)
	sys.stdout.flush()
	return

def get_args():
	parser = argparse.ArgumentParser()
	# set argparse
	parser.add_argument('--gfile', type=str, default="", 
	help='gen file needed')
	parser.add_argument('--sfile', type=str, default="", 
	help='sample file needed')
	parser.add_argument('--ofile', type=str, default="crag", 
	help='output file stem, e.g. <ofile>.out')
	parser.add_argument('--chr', type=str, default="NA",
	help='chromosome number to use in output file')
	parser.add_argument('--obs', type=int, default=1, 
	help='which obs indicator is the competing risk (1 is default, do not specify 0 as it is reserved for censored observations)')
	parser.add_argument('--t_pheno', type=str, default="", 
	help='time to event in SAMPLE file')
	parser.add_argument('--et_pheno', type=str, default="", 
	help='event type indicator in SAMPLE file')
	parser.add_argument('--covs', type=str, default="", 
	help='comma delimited list of covariates in SAMPLE file')
	parser.add_argument('--x', type=str, default="",
	help='covariate to adjust for SNP covariate interaction, must be from SAMPLE file')
	parser.add_argument('--crtype', type=int, default=0, 
	help='use 0 for subdistribution competing risks (default), 1 for cause-specific competing risks and 2 for conventional Cox PH survival model')
	parser.add_argument('--verbose', type=int, default=0, 
	help='use 1 to print results of covariates tested at each SNP, includes interaction terms of specified')
	parser.add_argument('--thread', type=int, default=1, 
	help='specify number of threads to use (does not work on Windows)')
	args=parser.parse_args()
	gfile = args.gfile
	sfile = args.sfile
	ofile = args.ofile
	nchr = args.chr
	obs = args.obs
	t_pheno = args.t_pheno
	et_pheno = args.et_pheno
	covs = args.covs.split(',')
	x = args.x.split(',')
	crtype = args.crtype
	verbose = args.verbose
	thread = args.thread
	return gfile, sfile, ofile, nchr, obs, t_pheno, et_pheno, covs, x, crtype, verbose, thread


if __name__ == "__main__":
	gfile, sfile, ofile, nchr, obs, t_pheno, et_pheno, covs, x, crtype, verbose, thread = get_args()
	print("##############################")
	print("#      WELCOME TO CRAG       #")
	print("##############################\n")
	# Reduce sample to needed covs
	sub=sample_read(sfile,covs,x,t_pheno,et_pheno,obs)
	crtype_sample(crtype, sub, t_pheno, et_pheno,obs)

	# Check for low variance
	if not '' in covs:

		print("##############################")
		print("#    CLINICAL ANALYSIS...    #")
		print("##############################\n")
		sub = low_var_test(sub,x,covs,t_pheno)
		
	# Open gen file and perform calculation every line
	print("##############################")
	print("#   BEGIN GWAS ANALYSIS...   #")
	print("##############################\n")
	
	colNames = ('chr','snp','test','bp','ea','nea','p','info','maf','hwe','beta','se')
	num_workers = thread

	manager = Manager()
	results = manager.list()
	work = manager.Queue(num_workers)

	#start for workers
	pool = []
	for i in range(num_workers):
		p = Process(target=do_work, args=(work, results, sub))
		p.start()
		pool.append(p)

	i=1
	if gfile.endswith('gz'):

		with gzip.open(gfile,'rt') as f:
			iters = itertools.chain(f, (None,)*num_workers)
			for num_and_line in enumerate(iters):
				gwasprogress()
				work.put(num_and_line)
				i+=1

	else:
		with open(gfile,'rt') as f:
			iters = itertools.chain(f, (None,)*num_workers)
			for num_and_line in enumerate(iters):
				gwasprogress()
				work.put(num_and_line)
				i+=1

	for p in pool:
		p.join()

	oresults(results,colNames,ofile)


	print("\n\n##############################")
	print("#      END OF ANALYSIS       #")
	print("##############################\n")
