# Positioning for a single time-step using ToTaL
total <- function(y) {
  # Check we can do triangelation
  if(nrow(y)<3) return("Not enough angles") 
  
  # 0. Prepare data
  # Select the three best angles
  y_subset <- y %>% tidyr::drop_na(azimuth_location_mdf) %>%
    dplyr::group_by(locator_mac) %>%
    dplyr::slice_max(rssi) %>%
    dplyr::arrange(desc(rssi))
  y_subset <- y_subset[1:3,] 
  y_subset <- y_subset %>% dplyr::rename(x_coordinate=x_m, 
                                         y_coordinate=y_m)
  
  
  # Get locator coordinates
  l1 <- c(y_subset$x_coordinate[1], y_subset$y_coordinate[1])
  l2 <- c(y_subset$x_coordinate[2], y_subset$y_coordinate[2])
  l3 <- c(y_subset$x_coordinate[3], y_subset$y_coordinate[3])
  
  # Get angles
  a1 <- as.numeric(y_subset[1,"azimuth_location_mdf"]) * 180/pi
  a2 <- as.numeric(y_subset[2,"azimuth_location_mdf"]) * 180/pi
  a3 <- as.numeric(y_subset[3,"azimuth_location_mdf"]) * 180/pi
  
  # Skip if two angles the same: Note: This could be avoided by just going down
  # further down the list of angles, however this happens so rarely that we'll just skip it here
  if((a2==a1)|(a2==a3)|(a1==a3)) return("Identical angles")
           
  # Get complement of angle y
  a1 <- 360-a1
  a2 <- 360-a2
  a3 <- 360-a3
  
  # Get the opposite angle (foculator to locator), convert to rads
  a1 <- ((a1+180)%%360)*(pi/180)
  a2 <- ((a2+180)%%360)*(pi/180)
  a3 <- ((a3+180)%%360)*(pi/180)
  
  # Preparations done,
  # Implement the ToTal algorithm (http://www.telecom.ulg.ac.be/publi/publications/pierlot/Pierlot2014ANewThree)
  # 1. Compute modified locator coordinates
  x_prime_1 <- l1[1]-l2[1]
  y_prime_1 <- l1[2]-l2[2]
  x_prime_3 <- l3[1]-l2[1]
  y_prime_3 <- l3[2]-l2[2]
  
  # 2. Compute cotangents
  T12 <- pracma::cot(a2-a1)
  T23 <- pracma::cot(a3-a2)
  T31 <- (1-(T12*T23))/(T12+T23)
  
  # 3. Compute the modified circle center coordinates:
  x_prime_12 <- x_prime_1 + T12*y_prime_1
  y_prime_12 <- y_prime_1 - T12*x_prime_1
  x_prime_23 <- x_prime_3 - T23*y_prime_3
  y_prime_23 <- y_prime_3 + T23*x_prime_3
  x_prime_31 <- (x_prime_3+x_prime_1)+T31*(y_prime_3-y_prime_1)
  y_prime_31 <- (y_prime_3+y_prime_1)-T31*(x_prime_3-x_prime_1)
  
  # 4. Compute k_prime_31
  k_prime_31 <- (x_prime_1*x_prime_3)+(y_prime_1*y_prime_3)+T31*(x_prime_1*y_prime_3-x_prime_3*y_prime_1)
  
  # 5 Compute D (if D = 0 return error)
  D <- ((x_prime_12-x_prime_23)*(y_prime_23-y_prime_31))-((y_prime_12-y_prime_23)*(x_prime_23-x_prime_31))
  if(D==0 | is.na(D)) return("D equals zero")
  
  # Return position + positioning error
  x_R <- l2[1]+(k_prime_31*(y_prime_12-y_prime_23))/D
  y_R <- l2[2]+(k_prime_31*(x_prime_23-x_prime_12))/D
  
  return(c(x_R,y_R,D))
}