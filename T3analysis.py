
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
import pyper as pr
r = pr.R(use_pandas = True)
r("library(cmprsk)")
#import rpy2.robjects as robjects
#print(rpy2.__version__)
#from rpy2.robjects.packages import rpackages
# import R's "base" package
#base = rpackages.importr('base')
# import R's "utils" package
#utils = rpackages.importr('utils')
#cmprsk = rpackages.importr('cmprsk')

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
#print(sub.head(10))
# Open gen file and perform calculation every line
import gzip
colNames = ('chr','snp','bp','p','info','maf','hwe','beta','se','conv')
ss=pd.read_csv(args.file + "." + args.chr + ".snp-stats",sep="\t")
#CHR SNP BP P INFO MAF HWE BETA SE CONV
with gzip.open(args.file + "_" + args.chr + "_impute2c.gz",'rt') as f:
	masterDF = pd.DataFrame(columns = colNames)
	df1 = pd.DataFrame({'chr':[args.chr]})
	
	for i, line in enumerate(f, 0):
		ssline=ss.iloc[i,:];
		gl=line.split(' ');	#df2=pd.DataFrame(ssline[['RSID','position','information','MAF','HWE']],columns=('snp','bp','info','maf','hwe'))
		df2 = pd.DataFrame({'snp':[ssline['RSID']],'bp':[ssline['position']],'info':[ssline['information']],'maf':[ssline['MAF']],'hwe':[ssline['HWE']]})
		#print(df2)
		#print(gl[1]) 
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
		r.assign("rdata", sub)
		r.assign("gz", gp)
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
		#df1['chr'] = args.chr
		df2=df1.join(df2)
		pydata=df2.join(df3)
		masterDF=masterDF.append(pydata,ignore_index=True)
		#print(r("summary(test)"))
		#print(pydata)
		
		#print(pydata.iloc[:,0])
	mdf=masterDF.ix[:,colNames]
	print(mdf)
	
#cbind(as.numeric(stat_pre[i,1,with=FALSE]),as.character(stat_pre[i,2,with=FALSE]),as.numeric(stat_pre[i,3,with=FALSE]),summary(ressnps1)$coef[1,5],as.numeric(stat_pre[i,4,with=FALSE]),as.numeric(stat_pre[i,5,with=FALSE]),as.numeric(stat_pre[i,6,with=FALSE]),summary(ressnps1)$coef[1,1],summary(ressnps1)$coef[1,3],as.numeric(ressnps1$converged))

#In [15]: r("m <- betareg(LEV_LT3 ~ SIZE1 + PROF2 + GROWTH2 + AGE + IND3A, data = rdata, subset = LEV_LT3 > 0)")
#Out[15]: 'try({m <- betareg(LEV_LT3 ~ SIZE1 + PROF2 + GROWTH2 + AGE + IND3A, data = rdata, subset = LEV_LT3 > 0)})\n'
 
#In [16]: # OUTPUT MODEL SUMMARY
 
#In [17]: print r("summary(m)")
#try({summary(m)})

#path = '/work/projects/epipgx/users/LIV/Task3Analysis/EpiPGX/temp.sample'
#sub.to_csv(path,sep=" ",index=False)
#samR=StringIO.StringIO(sub)

# Iteration over all arguments:
#for eachArg in sys.argv:   
#        print eachArg

# Define command and arguments
#command = 'Rscript'
#path2script = '/work/projects/epipgx/users/LIV/Task3Analysis/EpiPGX/Competing2Risks.R'

# Variable number of args in a list
#argsR = [args.file, args.chr, args.chunk, args.obs]


# Build subprocess command
#cmd = [command, path2script] + argsR

# check_output will run the command and store to result
#x = subprocess.check_output(cmd, universal_newlines=True)

#print('The maximum of the numbers is:', x)
