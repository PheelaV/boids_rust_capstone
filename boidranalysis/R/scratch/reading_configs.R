list.files("../boids_app/configs")
library(RcppTOML)
config <- parseTOML("../boids_app/configs/string_ball_hybrid.toml")

test_config <- list(
no_iter = 8000,
init_boids = 512,
save_locations_path = "", # data_folder,
sample_rate = 32,
init_width = 4000,
init_height = 4000,
sensory_distance = 50,
alignment_coef = .02,
alignment_trs_coef = 1.15,
cohesion_coef = 0.002,
cohesion_trs_coef = .95,
separation_coef = 4.1,
separation_trs_coef = .3,
min_speed = .5,
max_speed = 2.,
max_steering = .65,
dbscan_clustering = T,
boundary_config = "{\"type\": \"Toroidal\"}"
)

test_config$init_boids = config$no_boids
test_config$save_locations_path = "" #config$save # config$save_timestamp
test_config$sample_rate = config$sample_rate
test_config$init_width = config$init_width
test_config$init_height = config$init_height
test_config$sensory_distance = config$sensory_distance
test_config$alignment_coef = config$alignment_coefficient
test_config$alignment_trs_coef = config$alignment_treshold_coefficient
test_config$cohesion_coef = config$cohesion_coefficient
test_config$cohesion_trs_coef = config$cohesion_treshold_coefficient
test_config$separation_coef = config$separation_coefficient
test_config$separation_trs_coef = config$separation_treshold_coefficient
test_config$min_speed = config$min_speed
test_config$max_speed = config$max_speed
test_config$max_steering = config$max_steering
test_config$dbscan_clustering = TRUE
test_config$boundary_config = "{\"type\": \"Toroidal\"}"


rlang::exec(flock_detailed, !!!test_config)


config_name <- "string_ball_hybrid.toml"

get_config <- function(config_name) {
  config_path <- paste0("../boids_app/configs/", config_name)
  if (!file.exists(config_path)) {
    stop("The config file does not exist")
  }

  config <- parseTOML(config_path)

  converted_config <- list()

  converted_config$init_boids = config$no_boids
  converted_config$save_locations_path = "" #config$save # config$save_timestamp
  converted_config$sample_rate = config$sample_rate
  converted_config$init_width = config$init_width
  converted_config$init_height = config$init_height
  converted_config$sensory_distance = config$sensory_distance
  converted_config$alignment_coef = config$alignment_coefficient
  converted_config$alignment_trs_coef = config$alignment_treshold_coefficient
  converted_config$cohesion_coef = config$cohesion_coefficient
  converted_config$cohesion_trs_coef = config$cohesion_treshold_coefficient
  converted_config$separation_coef = config$separation_coefficient
  converted_config$separation_trs_coef = config$separation_treshold_coefficient
  converted_config$min_speed = config$min_speed
  converted_config$max_speed = config$max_speed
  converted_config$max_steering = config$max_steering
  converted_config$dbscan_clustering = TRUE
  converted_config$boundary_config = "{\"type\": \"Toroidal\"}"

  return(converted_config)
}

