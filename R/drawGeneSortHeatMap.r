drawGeneSortHeatMap	<- function( seuratObj,  heatMapHeight, heatMapWidth, showCellNames = FALSE){

# This snippet makes a heatmap of cells sorted with decreasing ltk expression
# the lines are set explicitely in geneOrder
# the parameters control page dimensions and if to show the cell names at the X axis
#
# written by Seva Makeev, 2017 - 2021

require( ComplexHeatmap)
require( circlize)
require( viridis)
require( proxy)
library( dendsort)


heatMapGenes	<- rownames( seuratObj@data) 

nCol		<- 1024


geneOrder	<- 	c( 	"ltk"		,"sox9b"	, "foxo1a"	, "tfap2e"	, "tfap2a"	, "her9"	, 
				"foxg1b" 	, "snai1b"	, "alx4b"	, "hmx1"	, 
				"otx2b"		, "sox10"	, "impdh1b"	, "foxo1b"	, "tyr"		, 
				"pax7b"		, "mc1r"	, "id2a"	, "hmx4"	, "foxd3"	, 
				"ednrba"	, "kita"	, "mbpa"	, "phox2bb"	, "tfec"	, 
				"mitfa"		, "foxp4" 	, "hbp1"	,  "mycla"	, "pax7a"	, 
				"tyrp1b"	, "slc24a5"	, "oca2"	, "mlphb"	, "pmela"	, 
				"myo5aa"	, "pnp4a"	, "ets1"	, "fgfr3"	, "pax3a"	, 
				"dpf3"		, "smad9")


dataMatrix 	<- seuratObj@data[ geneOrder , order( seuratObj@data[ "ltk", ],  decreasing = TRUE)]

colPanelFun	 = colorRamp2( quantile( dataMatrix, seq(0, 1, by = 1/(nCol - 1))), viridis( nCol))

row_dend 		<- dendsort(hclust( dist(dataMatrix, method = "cosine", pairwise = TRUE)))

hMap		<- Heatmap( 	dataMatrix, 
			name 	= "gene expression",
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
			cluster_rows = TRUE,
			row_dend_reorder = TRUE,
			show_column_names = showCellNames, 
			row_names_gp = gpar(fontsize = 10), 
			column_names_gp = gpar(fontsize = 10),
			clustering_distance_rows = function(x) dist( x, method = "cosine", pairwise = TRUE), 
			use_raster = FALSE
#			 raster_device = "png"
			)


return( hMap)

}
