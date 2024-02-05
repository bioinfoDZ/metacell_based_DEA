Scripts for COVID-19 time series scRNA seq analysis using metacells.

A collection of scripts and associated metadata for regression-based analysis of PBMC expression changes through time. 
4 scripts are related to maSigPro analysis by cell type for SEACells metacells, random metacells, all cells, and pseudobulked samples.
1 script is related to analysis of SEACells metacell purity.
1 Python script for creating metacells by time point.

Note that the SEACells algorithm will not produce identical results with each run. Please email kevin.oleary@einsteinmed.edu if you are interested in obtaining the same SEACell metacell assignments that we use in our analysis.

Updates to packages, such as Seurat, may mean that some code needs slight modification.

We describe the contents of each file below:

**all.cells.masigpro.R**

We performed quadratic regression using maSigPro using all single cells. We subset by cell type and grouped cells by time since symptom onset before running the analysis to find differentially expressed gene through time.

**all.cells.masigpro.design.2.csv**

This was metadata for SEACells metacells that we used for the maSigPro design matrix.

**combined.metacells.design.xlsx**

This excel file has a variety of sheets, each of which represnt the design matrix used for random metacell maSigPro analysis. Each sheet is specific to a cell type.

**metacells.random.R**

In this script, we create random metacells after subsetting all cells by time and cell type. We create metacells by averaging the normalized expression of 20 random cells at each time point for each cell type. We then do regression analysis with maSigPro to find differentially expressed genes through time.

**pseudobulk.R**

We pseudobulked by cell type, sample, and time point before using maSigPro to find differentially expressed genes through time.

**puirty.tests.R**

In this script, we subset SEACElls metacells based on their purity score and compare the effect this has on variance across genes.

**seacells.allcovid.ipynb**

This Python script shows the creation of SEACElls metacells for each time point. Note that the output (metacell assignemnts) will differ slightly with each run due to the SEACell algorithm's random initialization step.

**seacells.masigpro.allcovid.1.2024.R**

This outlines the regression-based analysis we did for each cell type using SEACells generated metacells. Prior to this, we transformed Seurat objects to annData objects for export to Python so that we could use the SEACells algorithm, which was not developed for use in R. This script also includes a variety of plots we generated when exploring individual genes we found to be differentially expressed through time. While we find DEGs for all cell types, we did not further explore some of them due to low overall number of cells or osberved model overfitting.



