# T-REx-CNV-scripts
This python script extracts, summarizes, and visualize the CNVs called by CODEX2 [https://github.com/yuchaojiang/CODEX2] into 
* Distribution of CNV sizes 
* Distribution of number of CNVs found accross samples
* Top known genes from overlapped with called CNVs
* Top novel genes 
* Tables such as T-REx overlapped with DGV, T-REx non-overlapped with DGV, etc. 

## Notes
1. A CNV is annotated as known if it is 50% overlapping with a record in DGV [http://dgv.tcag.ca/dgv/app/home].
2. A CNV is annotated with a gene name if it is 50% overlapping with a record in GENCODE [https://www.gencodegenes.org/human/].
3. The visualization package is Bokeh [https://bokeh.org/].

# How to run 
> python summarizeCNVsFromCODEX2.py frac_files.txt results/CODEX2/CODEX2_combined/ dvgfile kitname gencodefile sampleSourceFiles.txt combined_sampname

where,
* frac_files.txt contains the list of result files (files endining with \_frac.txt) from CODEX2
* workingDir is the directory containing the \_frac files
* dgvfile is directory containing the DGV file e.g., DGV/GRCh37_hg19_variants_2016-05-15.txt
* kitname specifies a capture kit, e.g., SS5, or just a tag for describing the analyzed data
* gencodefile is directory containing the GENCODE anntoation file e.g., annotation/gencode.v30lift37.annotation.gff3
* sampleSourceFile list of sample sources, which we used to tag individual samples when building a table of CNV summary for individual samples.  
* sampNameFile contains the list of bam files, which we used only the file name as the representation of sample name for building a table of CNV summary for individual samples.
