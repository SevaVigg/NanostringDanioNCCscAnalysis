# NanostringDanioNCCscAnalysis

This repo contains two datasets of expression profiling of gene panels of Danio rerio developing neural crest cells
Expression profiling was performed in the lab of Prof. Robert Kelsh (Univ. of Bath, https://researchportal.bath.ac.uk/en/persons/robert-kelsh)
and conducted by Nanostring nCounter (R) for a gene panel (main poriton of files in SourceData folder) and with TaqMan assay
(tables in the Taqman subfolder of the SourceData folder)

The code in the R directory performes analysis and obtains Figures for the manuscript intended for submission in Nature

Title: Zebrafish pigment cells develop directly from highly multipotent progenitors.
Authors: Masataka Nikaido^*,1,a, Tatiana Subkhankulova^*,1,b, Leonid Uroshlev^2, Artem Kasianov^2,3, Karen Camargo Sosa^1, 
Gemma Bavister^1, Xueyan Yang^1,,c, Frederico S. L. M. Rodrigues^1, Thomas J. Carney^1,d, Jonathan H.P. Dawes^4, 
Andrea Rocco^5, Vsevelod Makeev^2,3 and Robert N. Kelsh^1

*These individuals contributed equally to this work

1 Department of Biology & Biochemistry, University of Bath, Claverton Down, Bath BA2 7AY, UK
2 Vavilov Institute of General Genetics, Russian Academy of Sciences, Ul. Gubkina 3, Moscow, 119991, Russian Federation.
3 Department of Medical and Biological Physics, Moscow Institute of Physics and Technology, 9 Institutskiy per., Dolgoprudny, Moscow Region, 141701, Russian Federation
4 Department of Mathematical Sciences, University of Bath, Claverton Down, Bath BA2 7AY, UK
5 Department of Microbial Sciences, FHMS, University of Surrey, GU2 7XH Guildford, UK

BA2 7AY, UK
Current addresses:
a Current address: Graduate School of Life Science, University of Hyogo, Ako-gun, Hyogo Pref., 678-1297, Japan
b Current address: Department of Genetics, Evolution and Environment, Division of Biosciences
Medical Sciences Building, University College London, Gower Street, London WC1E 6BT
c Current address: The MOE Key Laboratory of Contemporary Anthropology, School of Life Sciences, Fudan University, Shanghai 200438, PR China
d Current address: Lee Kong Chian School of Medicine, Experimental Medicine Building, Yunnan Garden Campus, 59 Nanyang Drive, Nanyang Technological University, 636921, Singapore


The code was tested at my rather old MacBook Pro  (13-inch, Mid 2012, 2.5 GHz Intel Core i5, 16 GB 1600 MHz DDR3) running MacOS 10.12.6.
I used R console ver. 3.5.2, using Seurat 2.3.4, slingshot 1.0.0, and tidyverse 1.2.1 with ggplot2 3.3.3

To start the program pull the repo into a folder, launch R console, set the folder as working by typing in R console setwd( "RunFolder") and type source("R/ProcessData.r")

The repository contains the demo version of the code, which reads the optimal coarse grain cluster structure from the file
/Res/ClusterData/bestUmap_regular_renamed.rObj
and the vizualization perspective from the file
/Res/ClusterData/visualisation2Dumap_renamed.rObj

This demo version must create a number of subfolders of the Res folder containing intermediary files (e.g. scTables, InitialTables, QualityControl) and plots
(in the subfolder Plots and its subfolders). Note that the main quality control plot is found in /Res/QualityCotrol/Plot/
Demo version running takes about 15mins on my rather old MacBookPro. Normally the version must reproduce all plots found in the main Figures, Extended Data Figures and Supplementary Materials. Some small deviations may appear from the published figures as the algorithms are stochastic. The most non-robust figure is the UMAP plot with mutants ( initCellPlotMut.png, aka Extended Figure 1b).

To run the complete version of the code please go to the script ProcessData.r and uncomment execution of findBestUmapCluster (line 197)
and comment reading of the saved cluster structure (line 201). It may take up to 40h to complete optimization

To optimize for visualization perspective please go to the script ProcessData.r and uncomment execution of map2umap (lines 234)
and comment reading of the saved 2D umap map (line 238). Visulization runs quickly, but about 80% of maps have strange layout of clusters with
ltHMPs outlying above eHMPs, and trajectories starting from eHMP and bypassing back eHMP cluster ih hiD, which is visually misleading (recall that trajectories are calculated in 42D probe count space).

The data are distributed under  CC BY licence, the software under MIT licence.


