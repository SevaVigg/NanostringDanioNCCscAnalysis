makeVlnPlots	<- function( seuratObj, name = "vlnPlots", orientation = "portrait", plotDPI = "100"){

# this functions makes all VlnPlots either for taqman or nanostring genes
# the genes are listed in explicit lists, controlling their position in the plot
# paramteres control resolution (plotDPI) and the figure name
#
# written by Seva Makeev, 2017 - 2021

source("R/setClusterColors.r")

resDir		<- file.path(getwd(), "Res")

plotDir		<- file.path(resDir, "Plots")
dir.create(plotDir, showWarnings = FALSE)

vlnPlotDir <- file.path( plotDir, "vlnPlots")
dir.create( vlnPlotDir, showWarnings = FALSE)


if (orientation == "landscape") { pageWidth = 18; pageHeight = 13 }
if (orientation == "portrait") { pageWidth = 18; pageHeight = 13 }

if( seuratObj@project.name == "taqman"){ 
	
vlnPlot <- VlnPlot( seuratObj, c( 	"sox9b"		, "snai1b"	, "sox10"	, "pax7b"	, "mbpa"	,
					"phox2bb"	, "mitfa"	, "ltk"		, "neurog1"	, "pnp4a"	,
					 "tyrp1b"	, "xdh"		, "elavl3"),
			nCol = 5, cols.use = setClusterColors( seuratObj), do.return = TRUE) 

}else{

vlnPlot <- VlnPlot( seuratObj, c( 	"tfap2e"	, "tfap2a"	, "her9"	, "sox9b"	, "foxg1b" 	,
					"snai1b"	, "alx4b"	, "hmx1"	, "otx2b"	, "sox10"	, 
					"impdh1b"	, "foxo1b"	, "tyr"		, "pax7b"	, "mc1r"	, 
					"id2a"		, "hmx4"	, "foxd3"	, "ednrba"	, "kita"	, 
					"mbpa"		, "phox2bb"	, "tfec"	, "mitfa"	, "foxp4" 	, 
					"foxo1a"	, "hbp1"	, "ltk"		, "mycla"	, "pax7a"	, 
					"tyrp1b"	, "slc24a5"	, "oca2"	, "mlphb"	, "pmela"	, 
					"myo5aa"	, "pnp4a"	, "ets1"	, "fgfr3"	, "pax3a"	, 
					"dpf3"		, "smad9"),

	nCol = 5,  cols.use = setClusterColors(seuratObj), do.return = TRUE)
}

if( plotDPI == 600){
	
  vlnPlotDir600dpi <- file.path( vlnPlotDir, "600dpi")
  dir.create( vlnPlotDir600dpi, showWarnings = FALSE)

  ggsave( paste0( name, ".png"), path = vlnPlotDir600dpi, device = "png" , plot = vlnPlot, width = pageWidth, height = pageHeight, units = "cm", dpi = 600, scale = 4)

}else if( plotDPI == 100){

  vlnPlotDir100dpi <- file.path( vlnPlotDir, "100dpi")
  dir.create( vlnPlotDir100dpi, showWarnings = FALSE)

  ggsave( paste0( name, ".png"), path = vlnPlotDir100dpi, device = "png" , plot = vlnPlot, width = pageWidth, height = pageHeight, units = "cm", dpi = 100, scale = 4)

}else{ cat("Select plotDPI 600 or plotDPI 100")}

}

