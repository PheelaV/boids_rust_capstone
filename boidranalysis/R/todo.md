# TODO
1. Fristly, the heading are -180 to 180
  - make them positive only
  - restrict them to [0,360) - 360 and not 361 distinct values
2. fix bearings


either way, bearings do not look so interesting, maybe sampling, messing with counts will help?

actually, sampling should healp make the graph more intereseting by a lot - if you think about it, I now
track every single step, whose change in bearing is possibly very small, add those changes over a couple of steps and they will increase in the extremes - spread the distribution a bit wider


# buff out the bearings
#
# does this distribution of the turning angles settle down, how fast, which parameters, ammount of sampling?
#
# decide on ammount of sampling first
#
# then try some extreme values of the forces
#
# add in some clustering measures
#
# number of individuals of the clusters, how long do they last,
# [convex hull, minimum spanning tree = massive space, minimal number of boids] calculate sizes for different parameters

