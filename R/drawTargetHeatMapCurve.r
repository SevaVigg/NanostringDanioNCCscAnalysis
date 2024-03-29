drawTargetHeatMapCurve	<- function( seuratObj,  targetCurve, heatMapHeight, heatMapWidth){

# This snippet makes a heatmap of gene expression of cells along the target curve
# Curve gives smoothed expression values, use drawTargetHeatMapCells for original gene expression values for cells
#
# written by Seva Makeev, 2017 - 2021

require( ComplexHeatmap)
require( circlize)
require( viridis)
require( proxy)
library( dendsort)


source("R/setClusterColors.r")

heatMapGenes 	<- c( 	"tfap2e"	, "tfap2a"	, "her9"	, "sox9b"	, "foxg1b" 	,
			"snai1b"	, "alx4b"	, "hmx1"	, "otx2b"	, "sox10"	, 
			"impdh1b"	, "foxo1b"	, "tyr"		, "pax7b"	, "mc1r"	, 
			"id2a"		, "hmx4"	, "foxd3"	, "ednrba"	, "kita"	, 
			"mbpa"		, "phox2bb"	, "tfec"	, "mitfa"	, "foxp4" 	, 
			"foxo1a"	, "hbp1"	, "ltk"		, "mycla"	, "pax7a"	, 
			"tyrp1b"	, "slc24a5"	, "oca2"	, "mlphb"	, "pmela"	, 
			"myo5aa"	, "pnp4a"	, "ets1"	, "fgfr3"	, "pax3a"	, 
			"dpf3"		, "smad9")

dataMatrix 		<- seuratObj@data

row_dend 		<- dendsort(hclust( dist(dataMatrix, method = "cosine", pairwise = TRUE)))

nCol			<- 1024
colPanelFun	 	<- colorRamp2( quantile( dataMatrix, seq(0, 1, by = 1/(nCol - 1))), viridis( nCol))

 
lineageCells		<- rownames(targetCurve$curve)
lineageType		<- targetCurve$target


clusterColors		<- setClusterColors( seuratObj)
names(clusterColors)	<- levels( seuratObj@ident)


curveClust 		<- seuratObj@ident[ lineageCells]

annotbar		<- HeatmapAnnotation( 	Cluster = seuratObj@ident[lineageCells],
				    		col 	= list( Cluster = clusterColors), 
				    		height = unit(30, "points")
				    #annotation_legend_param = list( nrow = 10)
					)

hMap <- Heatmap( 	t(targetCurve$curve), 
			name = "gene log expression", 
			width 	= unit( heatMapWidth,  "inches"), 
			height 	= unit( heatMapHeight, "inches"),
			heatmap_legend_param = list(
				title = "log10(Exp)",
				legend_height = unit( 2, "inches"),
				grid_width = unit( 0.2, "inches"),
				title_gp = gpar( fontsize = 12, fontface = "plain")		
						),
			col	=  colPanelFun,
			cluster_columns = FALSE,
			cluster_rows = row_dend,
			row_dend_reorder = TRUE,
			show_column_names = FALSE, 
			row_names_gp = gpar(fontsize = 12), 
#			column_names_gp = gpar(fontsize = 10), 
			bottom_annotation = annotbar,
			column_title = paste0("Expression of ", lineageType, " lineage"), 
			clustering_distance_rows = function(x) dist( x, method = "cosine", pairwise = TRUE) , 
			use_raster = FALSE
#			raster_device = "png" 
			)

return( hMap)
}
