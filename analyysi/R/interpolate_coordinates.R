# Given floor polygon in data frame (vertex) format
# converts angles dataset y, inclusion / exclusion polygon data frame and
# path data frame from WGS84 to meters
interpolate_coordinates <- function(floor_polygon, y,
                                    exclusion_polygons=NA,
                                    inclusion_polygons=NA,
                                    path=NA) {
  
  # Set distance division factor, if linear interpolation is used this is disabled
  if("data.frame" %in% class(floor_polygon)) {
    conversion_factor <- 1
    
    # Do the linear interpolation to all of the coordinates
    bounding_box <- sf::st_bbox(sfheaders::sf_polygon(floor_polygon, x = "x", y = "y", polygon_id = "id"), crs = sf::st_crs(4326))
    bounding_box_min_lon <- bounding_box[1]
    bounding_box_min_lat <- bounding_box[2]
    bounding_box_max_lon <- bounding_box[3]
    bounding_box_max_lat <- bounding_box[4]
    bounding_box_lon_distance <- raster::pointDistance(c(bounding_box_min_lon, bounding_box_min_lat), c(bounding_box_max_lon, bounding_box_min_lat), lonlat=T)
    bounding_box_lat_distance <- raster::pointDistance(c(bounding_box_min_lon, bounding_box_min_lat), c(bounding_box_min_lon, bounding_box_max_lat), lonlat=T)
    
    # Then start converting
    lon_conversion <- function(x) ((x - bounding_box_min_lon) / (bounding_box_max_lon - bounding_box_min_lon)) * bounding_box_lon_distance 
    lat_conversion <- function(x) ((x - bounding_box_min_lat) / (bounding_box_max_lat - bounding_box_min_lat)) * bounding_box_lat_distance 
    
    floor_polygon$x <- vapply(floor_polygon$x, lon_conversion, FUN.VALUE=0)
    floor_polygon$y <- vapply(floor_polygon$y, lat_conversion, FUN.VALUE=0)
    floor_polygon$label_position_x <- vapply(floor_polygon$label_position_x, lon_conversion, FUN.VALUE=0)
    floor_polygon$label_position_y <- vapply(floor_polygon$label_position_y, lat_conversion, FUN.VALUE=0)
    
    y$lon <- vapply(y$lon, lon_conversion, FUN.VALUE=0)
    y$lat <- vapply(y$lat, lat_conversion, FUN.VALUE=0)
    
    if("data.frame" %in% class(exclusion_polygons)) {
      exclusion_polygons$x <- vapply(exclusion_polygons$x, lon_conversion, FUN.VALUE=0)
      exclusion_polygons$y <- vapply(exclusion_polygons$y, lat_conversion, FUN.VALUE=0)
      exclusion_polygons$label_position_x <- vapply(exclusion_polygons$label_position_x, lon_conversion, FUN.VALUE=0)
      exclusion_polygons$label_position_y <- vapply(exclusion_polygons$label_position_y, lat_conversion, FUN.VALUE=0)
    }
    if("data.frame" %in% class(inclusion_polygons)) {
      inclusion_polygons$x <- vapply(inclusion_polygons$x, lon_conversion, FUN.VALUE=0)
      inclusion_polygons$y <- vapply(inclusion_polygons$y, lat_conversion, FUN.VALUE=0)
      inclusion_polygons$label_position_x <- vapply(inclusion_polygons$label_position_x, lon_conversion, FUN.VALUE=0)
      inclusion_polygons$label_position_y <- vapply(inclusion_polygons$label_position_y, lat_conversion, FUN.VALUE=0)
    }
    if("data.frame" %in% class(path)) {
      path$x <- vapply(path$x, lon_conversion, FUN.VALUE=0)
      path$y <- vapply(path$y, lat_conversion, FUN.VALUE=0)
    }
  }
  return(list(y, exclusion_polygons, inclusion_polygons, path))
}

