## A meta-analysis of epigenome-wide association study of ultra-processed food consumption with DNA methylation in European children 

Code used in the following study: a meta-analysis of epigenome-wide association study of ultra-processed food consumption with DNA methylation in European children.  

## Description
List of scripts used to perform the different steps in the EWAS. 

1 - Descriptive Analysis: Given a dataframe with nrow/samples x ncol/variables get summary tables and plots.
2 - Association Analysis: 
3 - Run EWAS: Run robsut linear models in parallel given a model to assess the associaiton between UPF consumption (exposure) and DNAm (outcome)
4 - preQC EASIER: From the output (a dataframe) after running the rlm arrange the names to match the input speciffications to run the Quality Control from EASIER package. 
5 - Run QC EASIER: Run the Quality Control from EASIER package to remove problematic CpGs and check the DNAm data.
6 - Run GWAMA EASIER: Given different dataframes with the output from QC, use GWAMA to perfrom a meta analaysis. 
7 - Read GWAMA EASIER: Read the ouptut from GWAMA and extract top and significant CpGs.
A - Functional Analysis - preEnrichment: Prepare input to perfrom Enrichment Analysis using Enrichment module from EASIER package
A - Functional Analysis - Enrichment: Run statistical analysis to search biological terms and patwhays related to top CpGs. 
B - DMRs - CorrMatrix for DMRs - Given a dataframe with methyaltion values compute a correaltion matrix to meta analyze algonside other data and search for significant DMRs using dmrff package. 
B -  DMRs -  Meta DMRs - Given a groups of correlation matrix search DMRs using dmrff package 

## Considerations

To run GWAMA software it must be installed where the code is running. 
The code msolty prepare inputs-ouputs to use EASIER package, speifically designed to run EWAS. We suggest to also look at its documentation. 
