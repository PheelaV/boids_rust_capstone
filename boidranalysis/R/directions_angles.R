# source("./R/helpers.R")

get_headings <- function(x, y, radians = F) {
  if(radians) {
    mapply(function(x, y) atan2(x, y), diff(x, 1), diff(y,1))
  } else {
    mapply(function(x, y) rad2deg(atan2(x, y)), diff(x, 1), diff(y,1))
  }
}

get_bearings <- function(headings) {
  # because the previous heading could have been close to but larger than -180
  # and the next one close to but smaller than 180

  # that would translate to >0 and <360
  # now subtracting those two, yields one end of the range, -360
  # and considering the switched case yields the second end of the range, 360

  # thus a correction is needed as we are interested in unique values, not wrap-arounds
  bearings_wrapped <- diff(headings + 180) * -1
  (if_else(bearings_wrapped < 0, 360 + bearings_wrapped, bearings_wrapped + 0) + 180) %% 360 - 180
}
