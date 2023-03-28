if (FALSE) {
  for (f in grep("./Data/0327_experiment2_s_noags[0-9]1_noise[0-9][0-9]*", list.dirs(), value = TRUE)) {
    for (ff in list.files(f, "boids")) {
      print(file.path(f, ff))
      # file.remove(file.path(f, ff))
    }
  }
}
