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

The code was tested under R ver. 3.5.2, using Seurat 2.3.4 and slingshot 1.0.0

to start the program write source("R/ProcessData.r") from the R console. 
