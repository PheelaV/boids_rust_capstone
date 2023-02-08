for (i in 1:10) {
  flock_return(ino_iter = 2400, init_boids = 128, save_locations_path = data_folder, sample_rate = 4, init_width = 5000, init_height = 5000)
}
installed.packages("microbenchmark")


# Unit: milliseconds
# expr
# flock_return(no_iter = 100, init_boids = 128, save_locations_path = data_folder,      sample_rate = 4, init_width = 5000, init_height = 5000)
# flock_return(no_iter = 100, init_boids = 256, save_locations_path = data_folder,      sample_rate = 4, init_width = 5000, init_height = 5000)
# flock_return(no_iter = 100, init_boids = 512, save_locations_path = data_folder,      sample_rate = 4, init_width = 5000, init_height = 5000)
# flock_return(no_iter = 1000, init_boids = 128, save_locations_path = data_folder,      sample_rate = 4, init_width = 5000, init_height = 5000)
# flock_return(no_iter = 1000, init_boids = 256, save_locations_path = data_folder,      sample_rate = 4, init_width = 5000, init_height = 5000)
# flock_return(no_iter = 1000, init_boids = 512, save_locations_path = data_folder,      sample_rate = 4, init_width = 5000, init_height = 5000)
# min         lq      mean     median        uq       max neval
# 2.428963   2.534928   3.36073   2.617686   2.74946  38.14984   100
# 8.074048   8.355574  10.98037   8.549300   8.70471 101.44913   100
# 29.635456  30.187316  38.50114  30.427351  30.98438 109.39657   100
# 29.921390  32.099289  40.01446  33.116377  34.69734  93.49902   100
# 116.945940 130.657304 146.85568 140.534244 165.69687 189.44472   100
# 470.473401 524.442029 555.60955 544.520406 579.11922 683.97610   100
library(microbenchmark)
data_folder <- ""
microbenchmark(
  flock_return(no_iter = 100, init_boids = 128, save_locations_path = data_folder, sample_rate = 4, init_width = 5000, init_height = 5000),
  flock_return(no_iter = 100, init_boids = 256, save_locations_path = data_folder, sample_rate = 4, init_width = 5000, init_height = 5000),
  flock_return(no_iter = 100, init_boids = 512, save_locations_path = data_folder, sample_rate = 4, init_width = 5000, init_height = 5000),
  flock_return(no_iter = 1000, init_boids = 128, save_locations_path = data_folder, sample_rate = 4, init_width = 5000, init_height = 5000),
  flock_return(no_iter = 1000, init_boids = 256, save_locations_path = data_folder, sample_rate = 4, init_width = 5000, init_height = 5000),
  flock_return(no_iter = 1000, init_boids = 512, save_locations_path = data_folder, sample_rate = 4, init_width = 5000, init_height = 5000),
  flock_return(no_iter = 1000, init_boids = 1024, save_locations_path = data_folder, sample_rate = 4, init_width = 5000, init_height = 5000),
  flock_return(no_iter = 1000, init_boids = 4096, save_locations_path = data_folder, sample_rate = 4, init_width = 5000, init_height = 5000)
)
