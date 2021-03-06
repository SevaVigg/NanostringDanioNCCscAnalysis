getTargetCurve <- function(

# extracts the target curve (ends in special cell type) out of the slingshot object
#
# written by Seva Makeev, 2017 - 2021

				slingObjs, 
				target = "I", 
				dimRed = "umap", 
				dimsUseHD = 1:2,
				cellsKeepThresh = 0.95){

require( "slingshot")

slingLins 	<- slingObjs$DimH@lineages

LineageId <- which(unlist(lapply( slingLins, function(x) tail(x, 1) == target)))  

curve 		<- slingObjs$DimH@curves[[ LineageId]]
curveOrd 	<- curve$s[ curve$ord, ]
curveWt		<- curve$w[ curve$ord ]
curveData	<- curveOrd[ curveWt > cellsKeepThresh, ]
targetCurve 	<- list( target = target, curve = curveData)

return( targetCurve)
}
 
