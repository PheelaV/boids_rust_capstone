# Generated by extendr: Do not edit by hand

# nolint start

#
# This file was created with the following call:
#   .Call("wrap__make_boidr_wrappers", use_symbols = TRUE, package_name = "boidr")

#' @docType package
#' @usage NULL
#' @useDynLib boidr, .registration = TRUE
NULL

flock <- function(no_iter, init_boids, save_locations_path) .Call(wrap__flock, no_iter, init_boids, save_locations_path)

#' executes flocking and returns a dataframe with the location data
flock_return <- function(no_iter, init_boids, save_locations_path, sample_rate, init_width, init_height) .Call(wrap__flock_return, no_iter, init_boids, save_locations_path, sample_rate, init_width, init_height)

#' executes flocking and returns a dataframe with the location data
flock_detailed <- function(no_iter, init_boids, save_locations_path, sample_rate, init_width, init_height, sensory_distance, alignment_coef, cohesion_coef, separation_coef, alignment_trs_coef, cohesion_trs_coef, separation_trs_coef, min_speed, agent_steering, max_speed, max_steering, dbscan_clustering, boundary_config, distance_config, field_of_vision, rules_impl, wander_on, wander_coef, wander_rate, wander_radius, wander_distance, baseline_speed, wander_random) .Call(wrap__flock_detailed, no_iter, init_boids, save_locations_path, sample_rate, init_width, init_height, sensory_distance, alignment_coef, cohesion_coef, separation_coef, alignment_trs_coef, cohesion_trs_coef, separation_trs_coef, min_speed, agent_steering, max_speed, max_steering, dbscan_clustering, boundary_config, distance_config, field_of_vision, rules_impl, wander_on, wander_coef, wander_rate, wander_radius, wander_distance, baseline_speed, wander_random)

#' executes flocking and returns a dataframe with the location data
flock_detailed_no_return <- function(no_iter, init_boids, save_locations_path, sample_rate, init_width, init_height, sensory_distance, alignment_coef, cohesion_coef, separation_coef, alignment_trs_coef, cohesion_trs_coef, separation_trs_coef, min_speed, agent_steering, max_speed, max_steering, dbscan_clustering, boundary_config, distance_config, field_of_vision, rules_impl, wander_on, wander_coef, wander_rate, wander_radius, wander_distance, baseline_speed, wander_random) .Call(wrap__flock_detailed_no_return, no_iter, init_boids, save_locations_path, sample_rate, init_width, init_height, sensory_distance, alignment_coef, cohesion_coef, separation_coef, alignment_trs_coef, cohesion_trs_coef, separation_trs_coef, min_speed, agent_steering, max_speed, max_steering, dbscan_clustering, boundary_config, distance_config, field_of_vision, rules_impl, wander_on, wander_coef, wander_rate, wander_radius, wander_distance, baseline_speed, wander_random)

force_recompile <- function() .Call(wrap__force_recompile)

get_convex_hull <- function(x, y) .Call(wrap__get_convex_hull, x, y)

get_voronoi_areas <- function(x, y, width, height) .Call(wrap__get_voronoi_areas, x, y, width, height)

preprocess_file <- function(file_path, init_width, init_height, no_boids) .Call(wrap__preprocess_file, file_path, init_width, init_height, no_boids)

test_execution <- function() .Call(wrap__test_execution)


# nolint end
