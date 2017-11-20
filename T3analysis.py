
# run_max.py
import subprocess
import sys
import csv
import numpy as np
import pandas as pd
import gzip
import time
#import StringIO
import argparse
#import pyper as pr
#r = pr.R(use_pandas = True)
#r("library(cmprsk)")
import rpy2.robjects as robjects
r = robjects.r
#print(rpy2.__version__)
from rpy2.robjects.packages import importr
from rpy2.robjects import r, pandas2ri
pandas2ri.activate()
# import R's "base" package
base = importr('base')
#import R's "utils" package
utils = importr('utils')
cmprsk = importr('cmprsk')
import threading
from multiprocessing import Process, Manager
import itertools
import time
#from multiprocessing.dummy import Pool as ThreadPool 
	
def do_work(in_queue, out_list, sub):
	while True:
		item = in_queue.get()
		line_no, line = item

        # exit signal 
		if line == None:
			return

		# work
		#result = crgwas(line, line_no )
		# fake work
		i=line_no
		ssline=ss.iloc[i,:];
		
		gl=line.split(' ');	
		df2 = pd.DataFrame({'snp':[ssline['RSID']],'bp':[ssline['position']],'info':[ssline['information']],'maf':[ssline['MAF']],'hwe':[ssline['HWE']]})
		gp0=gl[5::3];
		gp1=gl[6::3]; 
		gp2=gl[7::3];
		gp = [];a_l=len(gp1)
		for item in range(a_l):
			if (float(gp0[item])+float(gp1[item])+float(gp2[item]))>0:
				gp.append(float(gp1[item])+(2*float(gp2[item])))
			else:
				gp.append(np.nan)
		
		
		# sub = pd.DataFrame({'time_to': [150,200,300,140,500],
								# 'Reason_stop': [1,1,2,2,0],
								# 'S200_NEUROL_EXAM_RESULT': [1,0,1,1,0], 
								# 'S383_LEARNING_DISABILITY': [0,1,1,1,1],
								# 'S188_POSITIVE_FAMHX': [1,1,0,0,0], 
								# 'Age': [21,22,40,50,10]})		
		sub['gz']=gp
		#print(ft)
		#sub[args.covs.split(',')]
		sub[[args.t_pheno, args.et_pheno]] = sub[[args.t_pheno, args.et_pheno]].astype(int)
		robjects.globalenv['sub'] = sub
		#print(sub[args.t_pheno])
		#print(sub[args.et_pheno])
		#print(sub['gz'])
		test=cmprsk.crr(ftime=sub[args.t_pheno],fstatus=sub[args.et_pheno],cov1=sub.iloc[:,2:],failcode=args.obs,cencode=0)
		res=base.summary(test).rx2('coef')
		#resconv=base.summary(test).rx2('conv')
		#print(resconv)
		df3 = pd.DataFrame({'beta':[res.rx('gz',1)],'se':[res.rx('gz',3)],'p':[res.rx('gz',5)],'conv':[1]})
		df3['beta'] = df3['beta'].str[0]
		df3['se'] = df3['se'].str[0]
		df3['p'] = df3['p'].str[0]
		df1 = pd.DataFrame({'chr':[args.chr]})
		df2=df1.join(df2)
		result=df2.join(df3)
		#print(pydata)
		
		#result = rdata
		#result = (line_no, line)

		out_list.append(result)
	
def crgwas(line, line_no ):
		i=line_no
		ssline=ss.iloc[i,:];
		gl=line.split(' ');	
		df2 = pd.DataFrame({'snp':[ssline['RSID']],'bp':[ssline['position']],'info':[ssline['information']],'maf':[ssline['MAF']],'hwe':[ssline['HWE']]})
		gp0=gl[5::3];
		gp1=gl[6::3]; 
		gp2=gl[7::3];
		gp = [];a_l=len(gp1)
		for item in range(a_l):
			if (float(gp0[item])+float(gp1[item])+float(gp2[item]))>0:
				gp.append(float(gp1[item])+(2*float(gp2[item])))
			else:
				gp.append(-9)
		#print(' '.join(map(str,gp)));
		# PASS DATA FROM PYTHON TO R
		rdata=robjects.DataFrame(sub)
		gz=robjects.DataFrame(gp)
		# SHOW DATA SUMMARY
		#print(r("summary(rdata)"))

		# r("library(cmprsk)")
		r("a=rdata[-1,]")
		r("gz[gz==-9]<-NA")
		r("test<-crr(ftime=as.numeric(a[,1]),fstatus=a[,2],cov1=cbind(gz,a[,c(-1,-2)]),failcode=" + args.obs  + ",cencode=0)")
		r("beta=summary(test)$coef[1,1]")
		r("se=summary(test)$coef[1,3]")
		r("pval=summary(test)$coef[1,5]")
		r("conv=summary(test)$conv")
		df3 = pd.DataFrame(np.transpose(r.get("rbind(beta,se,pval,conv)")),columns=('beta','se','p','conv'))
		df1 = pd.DataFrame({'chr':[args.chr]})
		df2=df1.join(df2)
		pydata=df2.join(df3)
		print(pydata)
		return pydata;
		#print(pydata)
		
		#print(pydata.iloc[:,0])


parser = argparse.ArgumentParser()

# set argparse
parser.add_argument('--file', type=str, default="",
                    help='file stem needed')
parser.add_argument('--chr', type=str,
                    help='which chromosome')
parser.add_argument('--chunk', type=str,
                    help='which 10k chunk is required')
parser.add_argument('--obs', type=str, default=1,
                    help='which obs indicator is the reference group (0 is ALWAYS censored, 1 is default)')
parser.add_argument('--t_pheno', type=str, default="",
                    help='time to event in SAMPLE file')
parser.add_argument('--et_pheno', type=str, default="",
                    help='event type indicator in SAMPLE file')
parser.add_argument('--covs', type=str, default="",
                    help='comma delimited list of covariates in SAMPLE file')
args=parser.parse_args()



# Reduce sample to needed covs
sample = pd.read_csv(args.file + ".sample",sep=" ") 
list=args.covs.split(',')
list.insert(0,args.et_pheno)
list.insert(0,args.t_pheno)
sub=sample[list]
sub=sub.drop([0])

# Open gen file and perform calculation every line
colNames = ('chr','snp','bp','p','info','maf','hwe','beta','se','conv')
ss=pd.read_csv(args.file + "." + args.chr + ".snp-stats",sep="\t")
#CHR SNP BP P INFO MAF HWE BETA SE CONV
if __name__ == "__main__":	
	num_workers = 1
	
	manager = Manager()
	results = manager.list()
	work = manager.Queue(num_workers)
	
	#start for workers
	pool = []
	for i in range(num_workers):
		p = Process(target=do_work, args=(work, results, sub))
		p.start()
		pool.append(p)
	
	with gzip.open(args.file + "_" + args.chr + "_impute2c.gz",'rt') as f:
		iters = itertools.chain(f, (None,)*num_workers)
		for num_and_line in enumerate(iters):
			work.put(num_and_line)

	for p in pool:
		p.join()
		
	f_results=pd.concat(results)
	final=f_results.ix[:,colNames]
	print(final)
	#for i in range(1, len(results)):
	#	print(results[i])
	

