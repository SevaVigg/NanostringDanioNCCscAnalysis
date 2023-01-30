# NanostringDanioNCCscAnalysis
 
Analysis of Nanostring and Taqman data of developing neural crest single cell RNA profiling in Danio rerio. 
Code written by Uroshlev L.A. [1], Kasyanov A.S.[1,2], Makeev V.J.[1,3,4]

This repository contains software and two datasets of developing neural crest single cell expression profiling of Danio rerio. 
Expression profiling was performed in the lab of Prof. Robert Kelsh 
(Univ. of Bath, https://researchportal.bath.ac.uk/en/persons/robert-kelsh) 
with the help of the Nanostring nCounter (R) for the inhouse panel of 45 genes
(main portion of files in the SourceData folder) 
and with the TaqMan assay (tables in the Taqman subfolder of the SourceData folder)

The code in the R directory conducts the analysis and plots Figures for the manuscript accepted to Nature Communications.  

Title: Zebrafish pigment cells develop directly from highly multipotent progenitors.
Authors: Tatiana Subkhankulova[1], Karen Camargo Sosa[1], Leonid A. Uroshlev[2], Masataka Nikaido[1,a], 
Noah Shriever[1], Artem S. Kasianov[2,3,4], Xueyan Yang[1,b], Frederico S. L. M. Rodrigues[1], 
Thomas J. Carney[1,c], Gemma Bavister[1], Hartmut Schwetlick[5], Jonathan H.P. Dawes[5], Andrea Rocco[6,7], 
Vsevolod Makeev[2,3,8] and Robert N. Kelsh[*,1]

These authors contributed equally: Tatiana Subkhankulova and Karen Camargo Sosa

1 Department of Life Sciences, University of Bath, Claverton Down, Bath BA2 7AY, UK
2 Vavilov Institute of General Genetics, Russian Academy of Sciences, Ul. Gubkina 3, Moscow, 119991, Russian Federation.
3 Department of Medical and Biological Physics, Moscow Institute of Physics and Technology, 9 Institutskiy per., Dolgoprudny, Moscow Region, 141701, Russian Federation
4 A.A. Kharkevich Institute for Information Transmission Problems (IITP), RAS Bolshoy Karetny per. 19, build.1, Moscow 127051, Russian Federation
5 Department of Mathematical Sciences, University of Bath, Claverton Down, Bath BA2 7AY, UK
6 Department of Microbial Sciences, FHMS, University of Surrey, GU2 7XH Guildford, UK
7 Department of Physics, FEPS, University of Surrey, GU2 7XH Guildford, UK
8 Laboratory ‘Regulatory Genomics’, Institute of Fundamental Medicine and Biology, Kazan Federal University, 18 Kremlyovskaya street, Kazan 420008, Russian Federation

Current addresses:
a Current address: Graduate School of Science, University of Hyogo, Ako-gun, Hyogo Pref., 678-1297, Japan 
b Current address: The MOE Key Laboratory of Contemporary Anthropology, School of Life Sciences, Fudan University, Shanghai 200438, PR China
c Current address: Lee Kong Chian School of Medicine, Experimental Medicine Building, Yunnan Garden Campus, 59 Nanyang Drive, Nanyang Technological University, 636921, Singapore


The code was tested on my rather old MacBook Pro  (13-inch, Mid 2012, 2.5 GHz Intel Core i5, 16 GB 1600 MHz DDR3) running MacOS 10.12.6.
I used R console ver. 3.5.3, using Seurat 2.3.4, slingshot 1.0.0, and tidyverse 1.2.1 with ggplot2 3.3.3

To start the program pull this repo into a folder, launch R console, set the  working folder by typing in R console setwd( "RunFolder") and type source("R/ProcessData.r")

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

Published single cell data integration was conducted on the server with 128 cores (Intel Core I5 processors) and 512 Gb of RAM running 18.04.1-Ubuntu Linux 5.4.0-122-generic kernel. Batch correction and data integration scripts was developed on R 4.1.1 version using batchelor 1.14.0 package

The data are distributed under  CC BY licence, the software under MIT licence.


