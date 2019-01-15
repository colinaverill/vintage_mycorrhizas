#' Extract PRISM climate data for a set of latitiude longitude.
#' depends on raster, and sp packages.
#'
#' @param longitude vector of longitude
#' @param latitude  vector of latitude
#' @param prism.dir path to directory that contains prism files. Currently defaults to directory in Colin's data folder.
#'
#' @return returns of a data frame of 30-year mat and map normals as well as location-month specific tmean, tmin, tmax, and ppt.
#' @export
#'
#' @examples
prism_query <- function(latitude,longitude,prism.dir = '/fs/data3/caverill/PRISM/'){
  #grab longitude and latitude from data frame.
  points <- cbind(longitude,latitude)
  
  #extract 30-year products, create a data frame.
  map30 <- extract(raster(paste0(prism.dir,'PRISM_ppt_30yr_normal_800mM2_annual_bil/PRISM_ppt_30yr_normal_800mM2_annual_bil.bil')),points)
  mat30 <- extract(raster(paste0(prism.dir,'PRISM_tmean_30yr_normal_800mM2_annual_bil/PRISM_tmean_30yr_normal_800mM2_annual_bil.bil')),points)
  
  to_return <- data.frame(map30, mat30)
  
  #return a dataframe with the climate variables
  return(to_return)
}
