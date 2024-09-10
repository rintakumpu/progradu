#' @title
#' Particle filter positioning

#' @description 
#' `pf` is a simple positioning particle filter for research purposes. 
#' @param y A data frame of angle observations
#' @param N Number of particles. Larger the number of particles the more likely 
#' the algorithm is to estimate the hidden positioning state effectively. The smaller the number of used particles the faster the algorithm runs.
#' @param q Random noise (in meters) that is added to the movement model. Standard deviation of a truncated Normal distribution.
#' @param q_min Left truncation of the above truncated Normal distribution.
#' @param q_max Right truncation of the above truncated Normal distribution.
#' @param resampling The value \eqn{\tau} in adaptive resampling.
#' @param proposal Proposal distribution \eqn{q(\cdot)} to be used. One of `("likelihood", "prior")`
#' where `"likelihood"` is \eqn{q(x_k|x_{1:k-1}, y_k) = p(y_k|x_k)} and `"prior"` is 
#' \eqn{q(x_k|x_{1:k-1}, y_k) = p(x_k|x_k-1)}. **Note!** Only `"likelihood"` is implemented.
#' @param likelihood How likelihood is computed. 
#' One of `("multiplicative", "additive")`. If set to `"multiplicative"` 
#' the likelihoods from angles are multiplied, otherwise they are added to each other.
#' @param likelihood_threshold If a particle likelihood for azimuth or elevation is below this value and 
#' `likelihood_function_p1` or `likelihood_elevation_function_p1` is set to `FALSE` the weight of the particles
#' below this linear scale value is set to 1. This is done after likelihood evaluation but before factoring in elevation. Defaults
#' to zero, set to negative value to disable. 
#' @param likelihood_use_elevation If set to `TRUE` elevation azimuth is added to the likelihood model. This
#' is done either by multiplying the azimuth and elevation likelihoods (if `likelihood_elevation=="multiplicative"`) or
#' by adding them together (if `likelihood_elevation=="additive"`).
#' @param likelihood_elevation How elevation is factored into the likelihood.
#' One of `("multiplicative", "additive")`. If set to `"multiplicative"` 
#' the likelihoods from azimuth and elevation are multiplied, otherwise they are added to each other.
#' @param likelihood_function One of `c("vonmises_truncnorm", "vonmises", "normal")`. 
#' 
#' If set to `"vonmises_truncnorm"`, the PDF of von Mises distribution is used for unconstrainted azimuth
#' angles and truncated Normal distribution is used for constrainted angle.
#' If set to `"normal"` the likelihood 
#' function used for each locator is the likelihood function of Normal distribution.
#' 
#' @param likelihood_elevation_function Which likelihood function to use with elevation? One of `c("normal", "truncnorm")`.
#' `"truncnorm"` as well as their `"_angle"` counterparts are as above except parameters a and b are always 0 and 90 degrees, respectively. 
#' Likewise `"full"` and `"exponential` 
#' `"normal"` is the same as `"full"` but [stats::dnorm()] is used instead of the explicit likelihood function. 
#' @param likelihood_function_p1 Add 1 to azimuth likelihood?
#' @param likelihood_elevation_function_p1 Add 1 to azimuth likelihood?
#' @param likelihood_power_weight Power weight applied to total elevation. 
#' @param likelihood_elevation_power_weight Power weight applied to likelihood elevation.
#' @param likelihood_apply_constraints If azimuth angle constraints exists, sets the weights of the particles outside these constraints to zero when evaluating the likelihood.
#' @param rssi_threshold Azimuth/elevation angles that have RSSI value below this threshold
#' value are dropped from likelihood evaluation.
#' @param asset_tag_height If elevation is used when computing the likelihood, the assumed asset tag height
#' for computing elevation angles between locators and particles.
#' @param true_path A path data frame from the function [get_path_recording_data()]. A data frame with only one row is expected,
#' if larger data frame is provided only the first row is used.
#' @param verbose Print output?
#' @param verbose Print diagnostic output?
#' @param initial_positions How to generate initial particle positions. 
#' One of `("triangulation", "polygon")`. If set to `polygon`, `grid_polygon`
#' is used (see below). If set to `triangulation`, the triangulated position based on the first time step
#' data is used as the initial position of the particles.
#' @param position_estimate One of `("max", "weighted_mean")`. If set to `"max"`, the particle
#' with maximum weight is used as the position estimate of each iteration. If set to `"weighted_mean"`
#' a weighted mean of all of the particles is used.
#' @param zones_for_plotting Data frame of zone polygons in the polygon-to-plot format provided 
#' by [galoRe::create_polygons_for_plotting()]. Used in creating particle maps.
#' @param floor_for_plotting Data frame of floor polygons in the polygon-to-plot format provided 
#' by [galoRe::create_polygons_for_plotting()]. Used in creating particle maps.
#' @param disable_movement_model_on_sleep If set to `TRUE` the movement model is dampened by a factor
#' of 10 whenever tag is determined to be in a sleep mode. Algorithm described in this
#' document is used for sleep determination 
#' @param stationary_attenuation If tag is determined to be stationary (either because of `disable_movement_model_on_sleep`
#' above or because of interpolation) a factor with which to reduce the `q` noise value.
#' @param stationary_attenuation_mode How to apply `stationary_attenuation` when applicable. One of `("constant", "progressive")`.
#' If set to `constant`, constant `stationary_attenuation` factor is applied. If set to `progressive` the noise is
#' further attenuated by applying an additional factor `c` which equals the ordinal of the time step for which attenuation
#' is applied. So for example first attenuated time step the factor is `1/stationary_attenuation` then `1/2*stationary_attenuation`
#' and so on.
#' @param inclusion_polygons A data frame in the polygon-to-plot format provided by [galoRe::create_polygons_for_plotting()].
#' If provided, the weights of particles outside these polygons will be set to zero.
#' @param session_polygons a data frame in the polygon-to-plot format provided by [galoRe::create_polygons_for_plotting()].
#' If provided, [galoRe::pf.apply_session_polygon()] is used to make transitions between these polygon areas harder.
#' @param session_polygon_penalty when session polygon penalty is applied when using [galoRe::pf.apply_session_polygon()]
#' a factor by which to weight down particles transitioning to a different session polygon. Defaults to 1000.
#' @param starting_area_id If set to a session area polygon ID the penalty will only be applied
#' to non-starting areas. I.e. it's harder to get to the non-starting areas but easy to get back
#' to starting area again.
#' @param exclusion_polygons A data frame in the polygon-to-plot format provided by [galoRe::create_polygons_for_plotting()].
#' If provided, particles falling on these polygons will be removed by resampling whereas particles passing
#' through these polygons accrue a penalty set by `P` (see below).
#' @param exclusion_raycasting_threshold When applying raycasting to
#' exclusion polygons, allow intersections of length `exclusion_raycasting_threshold`
#' and below (in meters).
#' @param P The multiplicative factor applied to particles that have crossed an exclusion
#' polygon when raycasting. Defaults to zero, i.e. sets particles crossing exclusion polygons to zero.
#' @param interpolation Interpolate missing positions. This adds missing time steps
#' into the original data frame. That is, if interpolation is used, the value `t` is increased 
#' to cover all the missing time steps, for `f=1` this will be 
#' `t_new = max(timestamp, t) - min(timestamp)`. When interpolated time steps are encountered during the run of
#' the algorithm, no observation data will be found. When this happens the the
#' value of `q` will be set to `q_interpolation` (see below) and the particle weights
#' will be set to last valid values. 
#' @param q_interpolation `q` value to use for interpolated time steps.
#' @param q_min_interpolation `q_min` value to use for interpolated time steps.
#' @param q_max_interpolation `q_max` value to use for interpolated time steps.
#' @param burn_in An integer denoting how many time steps should be used as burn-in. For
#' this many steps the q value is increased by `burn_in_multiplier` 
#' and movement restrictions are removed. This allows the algorithm to converge near
#' the true position before the dynamic model kicks in. Burn-in steps 
#' are not discarded from the initial results but should be discarded later on. If
#' set to `NA` no burn-in is used.
#' @param burn_in_multiplier For `burn_in` time steps the `q`values are
#' multiplied by `burn_in_multiplier`.
#' @param plotting If set to `TRUE`, plots and corresponding animations are returned.
#' @param TODO: ADD SMOOTHING, SEQUENCE NUMBER STUFF.

#' @return A list with the following elements: 

#' @export
pf_positioning <- function(y, 
                N = 1024,
                q = 3, 
                q_min = 0, 
                q_max = 10, 
                resampling=2/3,
                rssi_threshold=-120,
                asset_tag_height = 1.2,
                test_path=NA,
                map_matching=F,
                P=0,
                exclusion_polygons=NA, 
                inclusion_polygons=NA,
                disable_movement_model_on_sleep = T,
                stationary_attenuation=10,
                burn_in=0, 
                burn_in_multiplier=10,
                smoothing=F,
                verbose=T) {
  
  ##########################
  # STEP 0. Initialization #
  ##########################
  
  # Set initial movement model disablement value
  disable_movement <- 0
  
  # Drop all NA angles if exist
  y <- tidyr::drop_na(y, "azimuth_location_mdf")
  y <- tidyr::drop_na(y, "azimuth_scale")
  
  # Convert input to data.table for faster processing
  y <- data.table::as.data.table(y)
  data.table::setkey(y, ts)
  
  # If exclusion zones provided, create a polygon set out of the data frame
  if("data.frame" %in% class(exclusion_polygons)) {
    exclusion_polygons_list <- list()
    for(id in unique(exclusion_polygons$id)) {
      exclusion_polygons_list[[length(exclusion_polygons_list)+1]] <- sf::st_polygon(list(cbind(exclusion_polygons[exclusion_polygons$id==id,]$x, 
                                                                                                exclusion_polygons[exclusion_polygons$id==id,]$y)))
    }
    exclusion_polygon_set <- sf::st_sfc(exclusion_polygons_list)
  }
  
  # If inclusion polygons provided, create a polygon set out of the data frame
  if("data.frame" %in% class(inclusion_polygons)) {
    inclusion_polygons_list <- list()
    for(id in unique(inclusion_polygons$id)) {
      inclusion_polygons_list[[length(inclusion_polygons_list)+1]] <- sf::st_polygon(list(cbind(inclusion_polygons[inclusion_polygons$id==id,]$x, 
                                                                                                inclusion_polygons[inclusion_polygons$id==id,]$y)))
    }
    inclusion_polygon_set <- sf::st_sfc(inclusion_polygons_list)
    inclusion_polygon_ids <- unique(inclusion_polygons$id)
  }
  
  t <- length(unique(y$ts))
  if(verbose) {
    cat(paste("Running particle filter with N: ", N, "; t: ", t, sep=""))
  }
  
  # Create weight, model and observation matrices based
  # on the N and t values, this creates all combinations of 0:t, 1:N
  w <- data.table::as.data.table(expand.grid(0:t,1:N))
  # k is the current time step, i the particle
  # in-line with notation used in the literature
  colnames(w) <- c("k","i")
  data.table::setkey(w, k)
  x <- w
  
  # Add matrix-specific variable "w" to matrix w
  w[, w:=rep(0, nrow(w))]
  
  # Initial (uniform) weights for time step k, set to -log(N)
  # i.e. log(1/N)
  w[.(0), w:=rep(-log(N), N)]
  
  # Copy w to w_original which will be used to store
  # non-resampled weights 
  w_original <- data.table::copy(w)
  # Same for x
  x_original <- data.table::copy(x)
  
  # Creates movement length noise from truncated normal distribution
  replicated_truncnorm <- replicate(t, truncnorm::rtruncnorm(N, a=q_min, 
                                                             b=q_max, 
                                                             mean=0, sd=q)) 
  # and direction from uniform
  psi <- replicate(t, runif(N, min=0, max=2*pi)) 
  lat_psi <- cos(psi)
  lon_psi <- sin(psi)
  # then store noise in the lists
  lat_noise_list <- split(replicated_truncnorm * lat_psi, rep(1:t, each=N))
  lon_noise_list <- split(replicated_truncnorm * lon_psi, rep(1:t, each=N))
  
  # Creates various vectors and lists used
  # to store results and diagnostics
  v <- rep(NA, t) # Variance
  movement_disabled <- c()
  v_lon <- rep(NA, t) # Weighted longitude variance
  v_lat <- rep(NA, t) # Weighted latitude variance
  qv <- rep(q, t) # q vector
  N_eff_vector <- rep(N, t) # Effective sample sizes
  resample_vector <- rep(FALSE, t) # Resampling performed?
  #L_vector <- rep(0, t) # List of used smoothing lags
  #L_var_vector <- rep(-Inf, t) # List of variances used with adaptive-lag smoothing
  alpha_vector <- rep(1, t) # Square of mean squared errors of likelihood function
  locator_n_vector <- rep(0, t) # How many locators were used each time step
  likelihood_sum_vector <- rep(0, t) # Sum of likelihoods
  likelihood_max_vector <- rep(0, t) # Max of likelihoods
  likelihood_mean_vector_normalized <- rep(0, t) # Normalized mean of likelihoods
  likelihood_list <- list() # All likelihoods
  azimuth_ch_list <- list() # List used to track convex hull standard deviations
  ch_df <- data.frame() # A data frame of convex hulls used in angle selection
  angles_df <- data.frame() # A data frame of angles used in the likelihood evaluation
  locators_used_list <- split(rep("", t), rep(1:t, each=1)) # Mac addresses of locators used
  interpolated <- rep(FALSE, t) # Interpolation indicator
  interpolation_intervals <- 0 # Interpolation counter(s)
  
  # If auxiliary particle filter is used to estimate variance, store resampling indicies (possible with smoothing
  # too)
  I <- list() # Resampling indices.
  I[[1]] <- 1:N
  E <- list() # List of lists. Enoch indices for each time step as presented in the paper on ALvar variance estimator.
  # Alternatively fixed-lag indices if SFLS is used
  E[[1]] <- I
  E_lambda <- rep(50, t) # Lambda value corresponding to each Enoch index run, size of the list stored in E. Let's start with
  # fixed lambda.
  r <- rep(1, t) # How many times we have resampled (APF/ALvar)
  
  # List to store the position estimates & candidates
  position_estimates <- list()
  
  # Or use a random grid created by sampling lat and lon axes
  # from uniform distribution
  lat <- runif(n = N, min = min(y$y_m, na.rm = T), max = max(y$y_m, na.rm = T))
  lon <- runif(n = N, min = min(y$x_m, na.rm = T), max = max(y$x_m, na.rm = T))
  
  # Save initial distribution to the x df
  x_0 <- as.data.frame(matrix(ncol=2,nrow=length(lat)))
  lat_lon_names <- c("lat", "lon")
  colnames(x_0) <- lat_lon_names
  x_0$lat <- lat
  x_0$lon <- lon
  x[.(0), lat := x_0$lat]
  x[.(0), lon := x_0$lon]
  
  # Get all the unique timestamps 
  all_ts <- sort(unique(y$ts))  
  # and set other initial values
  k <- 1
  ts_previous <- NA
  
  
  #######################
  # START THE ALGORITHM #
  #######################
  
  for(k_ts in all_ts) {
    
    # Output time
    if(verbose) print(paste("Time step ", k, ": ", k_ts, sep=""))
    
    # Compute time diff (in seconds)
    if(!is.na(ts_previous)) {
      ts_diff <- abs(k_ts - ts_previous)
    } else {
      ts_diff <- 1
    }
    
    # Get correct time-slice of observations / particles
    y_subset <- y[.(k_ts)]
    
    # Get previous lon/lat
    x_lon_previous <- x[eval(.(k-1)), "lon"][[1]]
    x_lat_previous <- x[eval(.(k-1)), "lat"][[1]]
    
    
    #################################
    # STEP 1: Apply movement model. #
    #################################
    
    # In which movement/dynamic model is applied to the particles and weights
    # are re-adjusted from proposal distribution.
    # Note on naming: Dynamic model is the entire model describing the particle
    # movement. Movement model is the dead-reckoning part of the dynamic model.
    
    # Apply sleep detection
    if(disable_movement_model_on_sleep) {
      # Movement isn't disabled, check for sleep
      if(disable_movement!=0) {
        disable_movement <- disable_movement - 1
      }
      # Check forward for ten seconds (or t-k if that's smaller)
      disable_movement_lag <- min(10, t-k)
      for(l in 1:disable_movement_lag) {
        # Check how many rows of data we have for that second
        if(nrow(y[.(k_ts + l)]) == 1) { # Missing time stap has one row of NAs. 
          # Movement is not disabled, let's use this lag as an initial estimate
          disable_movement = l
          disable_movement_max = l
        } 
      } 
    }
    movement_disabled[k] <- FALSE
    if(disable_movement_model_on_sleep && disable_movement > 0 && k > burn_in) {
      movement_disabled[k] <- TRUE
      generated_truncnorm <- truncnorm::rtruncnorm(N, a=q_min, 
                                                   b=q_max, 
                                                   mean=0, sd=q / stationary_attenuation)
      psi <- runif(N, min=0, max=2*pi)
      lat_psi <- cos(psi)
      lon_psi <- sin(psi)
      lat_noise_list[[k]] <- generated_truncnorm * lat_psi
      lon_noise_list[[k]] <- generated_truncnorm * lon_psi 
    }
    
    if(is.na(burn_in)) burn_in <- 0
    
    # Adjust the noise term further, if required
    # i.e. if we are missing seconds
    if(ts_diff != 1) {
      lat_noise_list[[k]] <- lat_noise_list[[k]] * ts_diff
      lon_noise_list[[k]] <- lon_noise_list[[k]] * ts_diff  
    }
    
    # Apply model
    x_lon <- x_lon_previous + lon_noise_list[[k]]
    x_lat <- x_lat_previous + lat_noise_list[[k]]
    x[eval(.(k)), lon := x_lon]
    x[eval(.(k)), lat := x_lat]
    
    ######################
    # STEP 2: Re-weight. #
    ######################
    
    # Use prior distribution
    w_k <- w[eval(.(k-1)), w]
    
    #################################################
    # STEP 3: Handle inclusion and exclusion zones. #
    #################################################
    
    # Remove particles outside inclusion_polygons
    # (set their weight to zero)
    if("data.frame" %in% class(inclusion_polygons)) {
      
      within_polygon_particles <- sp::point.in.polygon(point.x=x_lon,
                                                       point.y=x_lat,
                                                       pol.x=inclusion_polygons$x,
                                                       pol.y=inclusion_polygons$y,
                                                       mode.checked=F)
      # 0 is outside, as the weights are e^w
      # we set weights to -Inf which equals zero
      w_k[within_polygon_particles==0] <- -Inf
    }
    
    # Apply exclusion polygons if polygons provided, not in a burn-in period
    # and exclusion polygons are applied to particles 
    if("data.frame" %in% class(exclusion_polygons) && (is.na(burn_in) | (k > burn_in))) {
      # First we simply remove the particles that are within an exclusion polygon
      # First create a set from data frame
      particles_set <- sf::st_as_sf(x[eval(.(k)), ], coords=c("lon", "lat"))
      
      # And this line does a check for each particle if that particle is within an exclusion polygon
      # using sf::st_within
      # this results in a matrix with each particle/polygon combibnation
      # From which we can easily get the intersecting particles
      intersecting_particles <- rowSums(sf::st_within(particles_set, exclusion_polygon_set, sparse=F))
      
      # These are within exclcusion polygons, so just set to zero
      w_k[intersecting_particles!=0] <- -Inf
      
      # Next check the intersections from current position
      if(k>1 && P) {
        
        # Create lines
        lines_matrix_list <- lapply(apply(x[eval(.(k)), c("lon", "lat")], MARGIN=1, FUN=base::c, simplify=FALSE), 
                                    FUN = rbind, c(position_estimates[[(k-1)]]))
        lines_st <- sf::st_sfc(lapply(lines_matrix_list, FUN=sf::st_linestring))
        
        # Check whether each of the matrix's lines intersects with the exclusion polygon
        intersecting_particles <- rowSums(as.data.frame(sf::st_intersects(lines_st, 
                                                                          exclusion_polygon_set, 
                                                                          sparse = F, prepared = T)))
        
        # If intersections found, check if they are below threshold
        if(sum(intersecting_particles) > 0) {
          
          intersection_length <- 0
          # If threshold set, check if intersections less than that.
          if(exclusion_raycasting_threshold > 0) {
            for(p_index in which(intersecting_particles>0)) {
              
              # Then compute the number of intersections
              no_intersecting_points <- length(sf::st_intersection(lines_st[p_index], 
                                                                   sf::st_cast(exclusion_polygon_set, 
                                                                               "MULTILINESTRING"))) / 2
              
              # If (mod 2) isn't zero this particle cannot be included
              # because of the threshold
              if(no_intersecting_points %% 2 != 0) next()
              
              # Check if below threshold
              intersection_length <- sf::st_length(sf::st_intersection(lines_st[p_index], 
                                                                       exclusion_polygon_set), 
                                                   sparse = F, prepared = T) * conversion_factor  
              
              # Let the jumps of less than the threshold pass
              intersecting_particles[intersection_length<=exclusion_raycasting_threshold] <- 0
            }
          } else {
            
          }
        }  
        
        # For those that are >0 we downweight set the corresponding weights
      }
      
      if(is.na(P)) {
        w_k[intersecting_particles!=0] <- -Inf
      } else {
        w_k[intersecting_particles!=0] <- log(exp(w_k[intersecting_particles!=0])*P)
      }
      
      
      # Sum exp to see if we have eliminated all the particles
      # which is the case if sum is zero
      w_sum_exp <- sum(exp(w_k))
      
      # However now these particles don't estimate the position anymore
      # This can be side-stepped by re-sampling (which cannot be done
      # if all probabilities -Inf)
      
      if(w_sum_exp!=0) {
        # We draw from the original particles a sample with the weights
        # where intersecting_particles are removed, the sample size is
        # sum(intersecting_particles) i.e. the amount of TRUEs in that vector
        resample_index <- sample(x=1:N, size=sum(intersecting_particles!=0), 
                                 replace=T, prob=exp(w_k))
        # And then set the current particles + weights to the ones resampled
        x_resample <- data.table::copy(x[eval(.(k))][resample_index, ..lat_lon_names])
        x[eval(.(k))][intersecting_particles!=0, (lat_lon_names) := x_resample]
        w_k[intersecting_particles!=0] <- w_k[resample_index]
      }
    } 
    
    # If all the weights are -Inf i.e. all the particles are excluded,
    # skip the movement, exclusion, inclusion and session polygons
    # (reset to previous)
    w_k <- tidyr::replace_na(w_k, -Inf)
    if(all(w_k==-Inf) && k > 1) {
      w_k <- w_k_previous
      x[eval(.(k)), lon := x_lon_previous]
      x[eval(.(k)), lat := x_lat_previous]
    }
    
    ################################################
    # STEP 3b: Auxiliary step                      #
    # TODO: Store for Eve/Enoch index computation, #
    # "commented out"                              #
    ################################################
    
    if(1==2) {
      
      # Use these weights to sample parent trajectories, "presampling"
      if(((resampling=="adaptive") && (N_eff < (resampling_threshold*N))) || (resampling=="each-step"))  {
        resample_index <- sample(x=1:N, size=N, replace=T, prob=exp(v_k_normalized))
        r[k] <- r[(k-1)] + 1 # Increment resampling counter
      } else {
        resample_index <- 1:N # This is how we handle Eve/Enoch indices
        r[k] <- r[(k-1)]
      }
      
      ######################################################
      # Store resample index and compute Eve/Enoch indices #
      ######################################################
      
      I[[k]] <- resample_index
      
      # Get current lambda
      current_lambda <- E_lambda[k]
      if((k-current_lambda)<1) current_lambda <- (k - 1)
      if(current_lambda > 0 || eve_index) {
        
        # This list is created always
        E_current <- list()
        
        # Implementation of Eve index
        # By going through resample indices recursively, starting from the k=1 + 1
        if(eve_index) {
          for(i in 1:k) {
            
            # For i==1 the resample index is the original 1:N
            if(i==1) {
              E_current[[i]] <- I[[i]]
            } else {
              E_current[[i]] <- E_current[[(i-1)]][I[[i]]]
            }
          }  
        } else {
          # Otherwise use Enoch index, starting from k-current_lambda
          j <- 1
          for(i in unique(r[(k-current_lambda):k])) {
            if(i==unique(r[(k-current_lambda):k])[1]) {
              E_current[[j]] <- I[[i]]
            } else {
              E_current[[j]] <- E_current[[(j-1)]][I[[i]]]
            }
            j <- j + 1
          }
        }
        # Store index list (of lists)
        E[[k]] <- E_current
      } else {
        # If lambda zero, we don't store the eve indices at all
        E[[k]] <- list(I[[k]])
      }
      
      
    } else {
      # Store Eve/Enoch index, it's just i here, and it will be i always
      # if APF is not used
      current_E <- list()
      current_E[[1]] <- 1:N 
      E_list <- current_E
    }
    
    ###############################
    # STEP 4: Measurement update. #
    ###############################
    
    # Compute normalizing constant, note that this is a vector, the 
    # normalization is done in the next step, this evaluates the
    # likelihood p(y_k|x_k) for current set of observations y_k
    
    # Create a copy of angle data for likelihood computation, will be
    # subset according to several criteria
    y_subset_likelihood <- data.table::copy(y_subset)
    
    # First RSSI
    rssi_cutoff_count <- nrow(y_subset_likelihood[rssi >= rssi_threshold])
    if(rssi_cutoff_count>0) {
      if(verbose) print(paste("Dropping ", rssi_cutoff_count,
                                    " azimuths below RSSI threshold.", sep=""))
      y_subset_likelihood <- y_subset_likelihood[rssi >= rssi_threshold]
      if(verbose) print(paste("Angle observations remaining:", nrow(y_subset_likelihood)))  
    }
    
    # If less or just one row remaining
    # this step needs to be interpolated, TODO check how to handle
    if(nrow(y_subset_likelihood)>1) {
      
      # Data for likelihood evaluation selected, next, create
      # a dataframe of locator/angle-particle pairs
      
      # Creates a new data frame for particle positions x
      x_likelihood <- data.table::as.data.table(matrix(ncol=2, nrow=N, dimnames = list(c(), c("lat_b", "lon_b"))))
      # And stores the input x latitude and longitude into this data frame
      x_likelihood$lat_b <- x_lat
      x_likelihood$lon_b <- x_lon
      
      # If we want to avoid joining by lat/lon values, let's add an id here
      x_likelihood$id <- 1:nrow(x_likelihood)
      y_subset_likelihood$id <- 1:nrow(y_subset_likelihood)
      # Do some renaming
      data.table::setnames(y_subset_likelihood, old = "y_m", new = "lat_a")
      data.table::setnames(y_subset_likelihood, old = "x_m", new = "lon_a")
      
      # Now we have all the locator lat-lon pairs in y_input and
      # particle lat-lon pairs in x_input
      # So combine them to xy
      xy <- data.table::CJ(x_id=x_likelihood$id, y_id=y_subset_likelihood$id)
      xy <- merge(xy, x_likelihood, by.x="x_id", by.y="id")
      xy <- merge(xy, y_subset_likelihood, by.x="y_id", by.y="id")
      
      # Now we have a dataset where we have each locator lat/lon
      # paired with each particle lat/lon
      # Compute angles between all the combinations
      # NOTE! Uses atan2(x,y) not atan2(y,x)
      xy <- xy[, angle := terra::atan2((as.numeric(lon_b)-as.numeric(lon_a)), 
                                       (as.numeric(lat_b)-as.numeric(lat_a)))]
      
      # Elevation
      # Create its own data frame for elevation
      xy_elevation <- data.table::copy(xy)
      # Compute angle between each particle and locator
      xy_elevation <- xy_elevation %>%
        dplyr::mutate(d = raster::pointDistance(matrix(c(lon_a, lat_a), ncol=2),
                                                matrix(c(as.numeric(lon_b), as.numeric(lat_b)), ncol=2), 
                                                lonlat = F, allpairs = F),
                      h = height - asset_tag_height,
                      elevation_angle_complement_rads = atan(d/h),
                      elevation_angle_computed = pi/2-elevation_angle_complement_rads)
      
      # Compute squared error(s)
      xy <- xy %>% dplyr::mutate(diff = pmin(abs(angle-azimuth_location_mdf), 
                                             2*pi-abs(angle-azimuth_location_mdf)), se = diff^2)
      xy_elevation <- xy_elevation %>% dplyr::mutate(diff_elevation = abs(elevation_angle_computed - elevation_location),
                                                       se_elevation = diff_elevation^2)
      
      xy[, L := .(galoRe::pf.dvonmises(angle, azimuth_location_mdf, azimuth_scale)), by = y_id]
      xy_elevation <- xy_elevation[, L := .(truncnorm::dtruncnorm(elevation_angle_computed, 
                                                                  a=0, b=pi/2, 
                                                                  mean=data.table::first(elevation_location), 
                                                                  sd=sqrt(data.table::first(elevation_scale)))), by=y_id]
      
      
      # Convert to log scale, compute MSE
      xy <- xy[, `:=`(l=log(L), mse=mean(se))]
      xy_elevation <- xy_elevation[, `:=`(l=log(L), mse=mean(se_elevation))]
      xy$L <- xy$L * xy_elevation$L 
      # Compute likelihood product.
      l_full <- xy[order(x_id),.(logSum=sum(log(L))), by=x_id]  
      log_p <- l_full$logSum
    }
    
    likelihood_list[[k]] <- log_p
    c_k_vector <- log_p + w_k
    
    # c_k is a normalizing constant, logSumExp simply
    # takes the exponents of c_k_vector then sums them and takes the
    # logarithm again. 
    # This function is to avoid numerical underflow which
    # would happen if using just log(sum(exp(x)))
    c_k <- matrixStats::logSumExp(c_k_vector)
    if(verbose) print(paste("Normalising constant: ", c_k, sep=""))
    
    #########################
    # STEP 5: Normalization #
    #########################
    
    # Make sure that c_k_vector doesn't have NaNs
    c_k_vector[is.nan(c_k_vector)] <- -Inf
    c_k <- matrixStats::logSumExp(c_k_vector)
    c_k_vector_normalized <- c_k_vector-c_k
    
    # Store to weights, however use c_k_vector_normalized in the future
    # to speed up rest of the algorithm
    w_k <- c_k_vector_normalized
    
    # Store previous weights,
    # this is used by auxiliary PF and the regular one in the case
    # that movement/exclusion pushes everything to zero (-Inf on log scale)
    w_k_previous <- w_k
    
    ###############################
    # STEP 6 Position estimation. #
    ###############################
    
    # Set weights to data frame (pre resample, done also here because of
    # position estimation function requiring the possibility to access earlier
    # weights)
    w[eval(.(k)), w := w_k]
    
    # Get current x_lon, x_lat vectors
    x_lon <- x[eval(.(k)), "lon"][[1]]
    x_lat <- x[eval(.(k)), "lat"][[1]]
    
    # Remove points with zero likelihood (i.e. -Inf on exp scale)
    non_zero_weights <- which(w_k!=-Inf)
    
    # Data frame of current particle positions with their weights, 
    position_df <- as.data.frame(cbind(x_lat[non_zero_weights],
                                       x_lon[non_zero_weights],
                                       w_k[non_zero_weights]))
    colnames(position_df) <- c("lat", "lon", "w")
    
    # If weighted mean is used
    position_estimate_k <- c(stats::weighted.mean(x=position_df$lat, w=exp(position_df$w)), 
                             stats::weighted.mean(x=position_df$lon, w=exp(position_df$w)))
    
    # Check if weighted mean falls within an exclusion polygon
    if("data.frame" %in% class(exclusion_polygons) && (is.na(burn_in) || (k > burn_in))) {
      
      # Convert position estimate to point
      position_estimate_point <- lapply(mapply(base::c, position_estimate_k[2], position_estimate_k[1], SIMPLIFY = F), sf::st_point)
      # Check for each particle if the position is within an exclusion polygon
      # using sf::st_within
      intersecting_position <- unlist(lapply(position_estimate_point, sf::st_within, y = exclusion_polygon_set, sparse = F))
      # this results in a vector of boolean values where FALSE is a valid position (i.e. it's not within
      # the exclusion polygon).
      # If it does, use max weight instead
      if(sum(intersecting_position)>0) {
        use_max_estimate_for_this_iteration <- T
        if(verbose) {
          print("Weighted mean position within exclusion polygon. Using position with max weight instead.")
        }
      } else {
        # Otherwise we're good.
        use_max_estimate_for_this_iteration <- F
      }
    } else {
      use_max_estimate_for_this_iteration <- F
    }
  
    
    # If max is requested or required
    if(use_max_estimate_for_this_iteration) {
      position_estimate_k <- as.numeric(position_df[which.max(position_df$w), 
                                                    c("lat","lon")])
      
      # If used as a backup for weighted mean, set indicator back to false
      use_max_estimate_for_this_iteration <- F 
    }
    
    # Store the estimate in a list
    position_estimates[[k]] <- position_estimate_k
    names(position_estimates)[[k]] <- k_ts
    
    ################################
    # STEP 7: Variance estimation. #
    ################################
    
    if(1==2) {
      variance_estimate_k <- pf.estimate_variance(x, w_k, k_current, E,
                                                  auxiliary = auxiliary,
                                                  variance_estimate = variance_estimate,
                                                  position_estimates = position_estimates,
                                                  conversion_factor = conversion_factor)
      v[k_current] <- variance_estimate_k[[1]]
      v_lon[k_current] <- variance_estimate_k[[2]]
      v_lat[k_current] <- variance_estimate_k[[3]]
      E_lambda[(k_current+1)] <- variance_estimate_k[[4]]   
    }
    
    ###################################
    # STEP 8: Resampling (optional). #
    ###################################
    
    # Compute the effective sample size
    N_eff_vector[k] <- round(1/sum(exp(w_k+w_k)))
    if(is.na(N_eff_vector[k])) {
      N_eff_vector[k] <- 0  
    }
    N_eff <- N_eff_vector[k]
    
    # Before resampling, store the original weight values
    w_original[eval(.(k)),w := w_k]
    x_original[eval(.(k)),lon := x_lon]
    x_original[eval(.(k)),lat := x_lat]
    
    # First check if resampling is used and that the step is not interpolated
    if(resampling>0) {
      # Then check what resampling mode is used and whether a possible adaptive
      # threshold is exceed
      if(N_eff < (resampling*N)) {
        # Check if all zeros, in that case can't resample, just set all weights to uniform
        if(sum(exp(w_k))==0) {
          if(verbose) print(paste("Effective sample size: ", N_eff, ". Weight vector equals zero, skipping resampling.", sep=""))
          w_k <- -log(N)
        } else {
          if(verbose) print(paste("Effective sample size: ", N_eff, ". Resampling.", sep=""))
          resample_vector[k] <- T
          # Resampling. Get particles at random using exp-weights as probabilities.
          resample_index <- sample(x=1:N, size=N, replace=T, prob=exp(w_k))
          
          # Store descendants for online smoothing
          # TODO: This might still be required, commenting out.
          #if(smoothing %in% c("SFLS")) {
          #  I[[k]] <- resample_index
          #} 
          
          # And then set the current particles to the ones resampled
          x_resample <- data.table::copy(x[eval(.(k))][resample_index, ..lat_lon_names])
          x[eval(.(k)), (lat_lon_names) := x_resample]
          
          # And reset weights to uniform
          w_k <- -log(N) 
        }
      } else {
        if(verbose) print(paste("Effective sample size ", N_eff, " larger than threshold value. Skipping resampling.", sep=""))
        
        # Store descendants for online smoothing
        # TODO: Again, indices need to be computed, leaving this in.
        #if(smoothing %in% c("SFLS")) {
        #  I[[k]] <- 1:N # This is how we handle Eve/Enoch indices
        #}
      }
    } else {
      # TODO ???
      I[[k]] <- 1:N # Interpolation on the first time step
    }
    
    # No resampling, just reset weights
    if(resampling==0) {
      w_k <- -log(N) 
    }
    
    # Compute particle ancestors (cf. Enoch index) if required
    # TODO, "commented out" for now
    if(1==2) {
      if(smoothing %in% c("SFLS")) {
        
        E_current <- list()
        
        # Fixed lag smoothing descendants
        if(k < smoothing_lag) {
          # Start from the first time step if not over
          # the lag horizon yet
          lookback_k <- 1
        } else {
          # Otherwise start lookback from whatever the lag is
          lookback_k <- k - smoothing_lag
          #print(length(lookback_k:k))
        }
        
        j <- 1
        for(i in lookback_k:k) {
          if(i==lookback_k) {
            E_current[[j]] <- 1:N
          } else {
            E_current[[j]] <- E_current[[(j-1)]][I[[i]]]
          }
          j <- j + 1
        }
        # Store index list (of lists)
        E[[k]] <- E_current
        
      }
    }
    
    ########################
    # STEP 9: Time update. #
    ########################
    
    # If at the end of burn-in period, replace burn-in period estimates 
    # with the first proper estimates
    if(!is.na(burn_in) && burn_in > 0 && (k-1) == burn_in) {
      for(b in 1:(k-1)) {
        position_estimates[[b]] <- position_estimates[[k]]
      }
    }
    
    # Set weights to data frame (post resample)
    w[eval(.(k)), w := w_k]
    if(k >= t) { break } else { k <- k + 1 }
    
    # Store current time as previous for computing time difference
    ts_previous <- k_ts
  }
  
  ###########################
  # ALL DONE, WRAP-UP STEPS #
  ###########################
  
  if(verbose) print("Positioning done. Wrapping up.")
  
  return(positioning_estimates)
  
  ###############################
  # A. Create a path data frame #
  ###############################
  
  # Also used to return various state values for each time step
  path_df <- as.data.frame(matrix(ncol=5, nrow=length(position_estimates)))  
  colnames(path_df) <- c("timestamp", "x", "y", "variance", "interpolated")
  path_df$timestamp <- as.numeric(names(position_estimates))
  path_df$y <- unlist(position_estimates)[seq(from=1, to=length(unlist(position_estimates)), by=2)]
  path_df$x <- unlist(position_estimates)[seq(from=2, to=length(unlist(position_estimates)), by=2)]
  path_df$y_prev <- c(NA, path_df$y[1:(nrow(path_df)-1)])
  path_df$x_prev <- c(NA, path_df$x[1:(nrow(path_df)-1)])
  path_df$variance <- v[1:nrow(path_df)]
  path_df$variance_lon <- v_lon[1:nrow(path_df)]
  path_df$variance_lat <- v_lat[1:nrow(path_df)]
  path_df$sd_m <-sqrt(path_df$variance) * conversion_factor
  path_df$asset_tag_id <- 0
  path_df$interpolated <- interpolated[1:nrow(path_df)]
  path_df$movement_disabled <- movement_disabled[1:nrow(path_df)]
  path_df$resampling_used <- resample_vector[1:nrow(path_df)]
  path_df$locator_n <- locator_n_vector[1:nrow(path_df)]
  path_df$likelihood_sum <- likelihood_sum_vector[1:nrow(path_df)]
  path_df$likelihood_max <- likelihood_max_vector[1:nrow(path_df)]
  path_df$likelihood_mean_normalized <- likelihood_mean_vector_normalized[1:nrow(path_df)]
  path_df$q <- qv[1:nrow(path_df)]
  path_df$q_min <- q_min
  path_df$q_max <- q_max
  path_df$group <- 1
  # TODO: Store smoothing once re-implemented.
  #path_df$lag <- L_vector[1:nrow(path_df)]
  #if(!is.na(smoothing_e)) path_df$lag_var <- L_var_vector[1:nrow(path_df)]
  
  # Set k values
  path_df$k <- 1:nrow(path_df)
  
  # Add some statistics from lists
  #path_df$mean_sd_m <- unlist(lapply(sd_list, mean))[1:nrow(path_df)]
  #path_df$median_sd_m <- unlist(lapply(sd_list, median))[1:nrow(path_df)]
  #path_df$mean_sd_m <- round(path_df$mean_sd_m, 4)
  #path_df$median_sd_m <- round(path_df$median_sd_m, 4)
  #path_df$mean_sd_m_silabs <- round(path_df$mean_sd_m_silabs, 4)
  #path_df$median_sd_m_silabs <- round(path_df$median_sd_m_silabs, 4)
 
  # Add true path and compute error statistics
  if("data.frame" %in% class(true_path)) { 
    # If multiple paths are provided just use the first in data frame (as only one
    # row is supposed to be provided)
    true_path_string <- true_path$path[1]  
    # Convert db coordinate string to vectors
    true_path_string <- stringr::str_replace_all(true_path_string, "\\[", "")
    true_path_string <- stringr::str_replace_all(true_path_string, "\\]", "")
    true_path_string <- stringr::str_replace_all(true_path_string, "\\)\\,\\(", "\\,")
    true_path_string <- stringr::str_replace_all(true_path_string, "\\(", "c\\(")
    true_path_string <- eval(parse(text=true_path_string))
    # Convert to df
    true_path_df <- as.data.frame(matrix(true_path_string, ncol=2, dimnames=list(c(), c("lat","lon")), byrow=T))
    # Assume constant speed and add time stamp
    true_path_df$timestamp <- as.numeric(seq(from=true_path$start_time[1], 
                                             to=true_path$end_time[1], 
                                             length.out=nrow(true_path_df)))
    
    # Densify the true path for error computation
    true_path_smooth <- as.data.frame(smoothr::smooth_densify(as.matrix(true_path_df[,1:2]), n=1000))
    colnames(true_path_smooth) <- c("lat", "lon")
    # Do a loop of path_df and compute the errors
    path_df$positioning_error <- rep(NA, nrow(path_df))
    if(apply_kalman_filter) path_df_kf$positioning_error <- rep(NA, nrow(path_df_kf))
    for(path_df_i in 1:nrow(path_df)) {
      # Compute distance to the true path i.e. positioning error for this time step
      path_df[path_df_i, "positioning_error"] <- min(raster::pointDistance(t(matrix(c(path_df[path_df_i, "x"], 
                                                                                      path_df[path_df_i, "y"]))), 
                                                                           cbind(true_path_smooth$lon, true_path_smooth$lat), lonlat=lonlat))
      
      
      # Kalman filtered df has the same number of rows
      # so we can do that in this same loop
      if(apply_kalman_filter) {
        path_df_kf[path_df_i, "positioning_error"] <- min(raster::pointDistance(t(matrix(c(path_df_kf[path_df_i, "x"], 
                                                                                           path_df_kf[path_df_i, "y"]))), 
                                                                                cbind(true_path_smooth$lon, true_path_smooth$lat), lonlat=lonlat))
      }
    }
  }
  
  # Add names/timestamps to locators used list
  names(locators_used_list) <- path_df$timestamp
  
  # Add row for k=0
  path_df[nrow(path_df)+1, ] <- rep(NA, ncol(path_df))
  path_df[nrow(path_df), "k"] <- 0
  path_df$N <- N
  path_df$N_eff <- c(N, N_eff_vector)[1:nrow(path_df)]
  
  #######################
  # F. Summarize errors #
  #######################
  
  errors_data_frame <- as.data.frame(matrix(ncol=0, nrow=10))
  rownames(errors_data_frame) <- c("min", "max", "mean", "rmse", "q50", "q80", "q90", "q95", "q99", "sd")
  if("data.frame" %in% class(true_path)) {
    errors_functions <- function(x) c(min(x, na.rm=T), max(x, na.rm=T), mean(x, na.rm=T), 
                                      sqrt(mean(x^2, na.rm=T)), quantile(x, 0.5, na.rm=T), 
                                      quantile(x, 0.8, na.rm=T), quantile(x, 0.9, na.rm=T), 
                                      quantile(x, 0.95, na.rm=T), quantile(x,0.99, na.rm=T),
                                      sd(x, na.rm=T))
    errors_data_frame$positions <- errors_functions(path_df$positioning_error)
  }
  
  ##########################################
  # G. Return variables, plots, parameters #
  ##########################################
  
  
  #variance_return_list <- list(variance_estimate, E_lambda, v, sqrt(v)*conversion_factor, v_lon, v_lat)
  #names(variance_return_list) <- c("Used variance estimate", "Lag values", "Variance", "Standard deviation (in meters)", "Longitudinal variance", "Latitudinal variance")
  return_list <- list(x,w,p_maps,#variance_return_list,
                      qv,path_df,path_df_kf,
                      locators_used_list,w_original,p_path,p_path_kf,
                      sitemaps,
                      list(lat_noise_list, lon_noise_list),
                      errors_data_frame, list(I, E), pf_animation_f)
  names(return_list) <- c("Particle positions",
                          "Particle weights", "Particle and angle maps",
                          "Vector of used q values",
                          "Position estimates", 
                          "Position estimates (Kalman filtered)",
                          "Locators used","Particle weights before resampling","Path plot", 
                          "Path plot (Kalman filtered)","Sitemaps", "List of noise terms (lat, lon)",
                          "Positioning errors", "Genealogy", "Animation save function")
  if(verbose) print("All done.")
  return(return_list)
}
