# source("./R/helpers.R")

# returns (nrow - unique(id)) records as it does a diff to normalize the coordinates
# results range [0, 2pi)
get_headings <- function(x, y, radians = T) {
  if(radians) {
    mapply(function(x, y) (atan2(x, y) + pi) %% (2 * pi), diff(x, 1), diff(y,1))
  } else {
    mapply(function(x, y) rad2deg((atan2(x, y) + pi) %% (2 * pi)), diff(x, 1), diff(y,1))
  }
}

get_headings2 <- function(dx, dy, radians = T) {
  if(radians) {
    mapply(function(dx, dy) (atan2(dx, dy) + pi) %% (2 * pi), dx, dy)
  } else {
    mapply(function(dx, dy) rad2deg((atan2(dx, dy) + pi) %% (2 * pi)), dx, dy)
  }
}

# returns (nrow - unique(id)) records as it does a diffs headings to get a bearing
#  results range [-pi, pi)
get_bearings <- function(headings) {
  # because the previous heading could have been close to but larger than -pi
  # and the next one close to but smaller than pi

  # that would translate to >0 and <2pi
  # now subtracting those two, yields one end of the range, -2pi
  # and considering the switched case yields the second end of the range, 2pi

  # thus a correction is needed as we are interested in unique values, not wrap-arounds
  bearings_wrapped <- diff(headings + pi) * -1
  (if_else(bearings_wrapped < 0, (2 * pi) + bearings_wrapped, bearings_wrapped + 0) + pi) %% (2 * pi) - pi
}
