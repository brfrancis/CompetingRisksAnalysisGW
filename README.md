
![](https://i.imgur.com/wWIRzgA.png)

[![Coverage Status](https://coveralls.io/repos/github/brfrancis/CompetingRisksAnalysisGW/badge.svg?branch=master)](https://coveralls.io/github/brfrancis/CompetingRisksAnalysisGW?branch=master)
[![Join the chat at https://gitter.im/CompetingRisksAnalysisGW/Lobby](https://badges.gitter.im/CompetingRisksAnalysisGW/Lobby.svg)](https://gitter.im/CompetingRisksAnalysisGW/Lobby)
[![DOI](https://zenodo.org/badge/110241944.svg)](https://zenodo.org/badge/latestdoi/110241944)

### What can *CRAG* do for you? 

Survival analysis is the study of "time-to-event" (TTE) outcomes, taking the focus off *if* an event will occur and instead looking at *when* an event will occur. Often, the when question is answered by a combination of factors, such as age or sex, survival analysis tests whether these associations are clear from data collected. So, if I'm asking *when will I respond to a medication?*, survival analysis can help to answer this question. Further information can be found on this great page by Cam Davidson-Pilon (author of lifelines) [here](http://lifelines.readthedocs.org/en/latest/Survival%20Analysis%20intro.html).

A genome-wide association study (GWAS), is a study of a genetic variants across the genome. Whether we are faced with a *when* or *if*, GWAS allows us to consider whether genetic factors are part of the combination of factors that answer these questions. 

*CRAG* is a novel software that performs survival analysis genome-wide and goes one step further: handling a particular type of TTE outcomes that occur for different reasons, such as drug withdrawal occuring due to either the drug not working or causing a bad reaction, appropriately modelling via “competing risks”. 

*CRAG* is an easy to use software implemented using Python, and is compatible with Linux, Mac OS X & Windows operating systems. *CRAG* can handle large-scale genome-wide data including imputed genotypes. Analyses can be undertaken with a cause-specific or sub-distribution model. The software can adjust for multiple covariates and incorporate SNP-covariate interaction effects.

The main use of CRAG is expected to be in pharmacogenetic studies to identify genetic biomarkers of patient response to treatment, with the ultimate goal of personalising therapeutic intervention for an array of diseases.

### Installation:
#### Dependencies:

Please install the following  Python packages prior to use of *CRAG* (some may be installed with Python): subprocess, sys, math, lifelines, csv, numpy, pandas, gzip, time, argparse and warnings.

#### Installing

Type the following command in the terminal/command line to intall *CRAG*

       pip install --upgrade --no-deps git+https://github.com/brfrancis/CompetingRisksAnalysisGW.git

from the command line.

#### Testing CRAG


### CRAG guide

If you are new to survival analysis, wondering why it is useful, or are interested in *lifelines* examples and syntax,
please check out the [Documentation and Tutorials page](http://lifelines.readthedocs.org/en/latest/index.html)


### Contacting
Please chat with author Ben Francis @ [Gitter](https://gitter.im/CompetingRisksAnalysisGW/Lobby). Email is also possible but likely to not be as timely!


### CRAG citation

If you find CRAG useful, the best way to show your appreciation is through citing, click the badge below to generate a DOI that can be put into most citation managers. A manuscript for CRAG is currently in production...

 [![DOI](https://zenodo.org/badge/110241944.svg)](https://zenodo.org/badge/latestdoi/110241944)
