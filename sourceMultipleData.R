####
# Sourcing multiple data
###

require(mp)
require(ggplot2)
require(ZedXL)
set.seed(102528)
computeproq3 <- F

get_density <- function(x, y, n = 1000) {
  dens <- MASS::kde2d(x = x, y = y, n = n)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}

compute.density.centromer <- function(funnelTable) {

  centromerX <- sum(funnelTable$projectionX * funnelTable$density)/sum(funnelTable$density)
  centromerY <- sum(funnelTable$projectionY * funnelTable$density)/sum(funnelTable$density)

  return(c("centromerX" = centromerX, "centromerY" = centromerY))
}

# SALB3.INIT.001
init_gscoreLogPath = '/home/guilherme/datasets/SALB3.NOCST.001/gscore/gscore-TMscore-050.dat'
init_loglistLocation = '/home/guilherme/datasets/SALB3.NOCST.001/loglist.txt'
init_topolinkLogsDirectory = '/home/guilherme/datasets/SALB3.NOCST.001/topolink_observed/'
init_alignlistPath = '/home/guilherme/datasets/SALB3.NOCST.001/alignlist.txt'
init_compactlogPath = '/home/guilherme/datasets/SALB3.NOCST.001/gscore/compactlog-TMscore.dat'
init_lovoalignLogPath = '/home/guilherme/datasets/SALB3.NOCST.001/lovoalign.log'
init_proq3listLocation = '/home/guilherme/datasets/SALB3.NOCST.001/proq3list.txt'
cryslistLocation = '/home/guilherme/datasets/SALB3.NOCST.001/utility/cryslist.txt'
optlistLocation = '/home/guilherme/datasets/SALB3.NOCST.001/utility/optlist.txt'
init_distanceTableLocation = '/home/guilherme/datasets/SALB3.NOCST.001/utility/distance_table.log'

init_modelScores <- compute.model.scores(type = "gscore",
                                         init_gscoreLogPath,
                                         init_lovoalignLogPath,
                                         init_proq3listLocation,
                                         computeproq3)

init_optimumXlinkMirttable <- create.xlink.mirttable(
  prepare.topolink.logs(
    read.topolink.output(mode,
                         init_loglistLocation,
                         init_topolinkLogsDirectory)
  )
)

init_optimumSimilarityTable <- create.dissimilarity.matrix(mode = "similarity",
                                                           diagonal = 1,
                                                           init_alignlistPath,
                                                           init_compactlogPath,
                                                           nmodels = nrow(init_optimumXlinkMirttable))

init_optimumXlinkMirttable <- digest.xlink.mirttable(init_optimumXlinkMirttable)

init_optimumDescript <- ltm::descript(init_optimumXlinkMirttable)

init_globalCronbachAlpha <- init_optimumDescript$alpha[1]

init_restrictionScores <- data.frame(bis = init_optimumDescript$bisCorr,
                                right = init_optimumDescript$perc[,2],
                                logit = init_optimumDescript$perc[,3],
                                alpha = init_optimumDescript$alpha[-1],
                                deltaAlpha = init_optimumDescript$alpha[2:length(init_optimumDescript$alpha)] - init_globalCronbachAlpha)

init_restrictionScores <- attribute.crys.and.opt(init_restrictionScores, cryslistLocation, optlistLocation)

init_modelDistanceTable <- 1/(init_optimumSimilarityTable)

init_projection.coordinates <- forceScheme(init_modelDistanceTable,
                                      Y = NULL,
                                      #max.iter = ,
                                      tol = 1e-04,
                                      #fraction = 8,
                                      eps = 1e-05)

init_projection.coordinates <- scale(init_projection.coordinates, center = T, scale = T)

init_funnelTable <- data.frame("projectionX" = init_projection.coordinates[,1],
                          "projectionY" = init_projection.coordinates[,2],
                          "tmscore" = init_modelScores$`TM-Score`,
                          "density" = get_density(x = init_projection.coordinates[,1], y = init_projection.coordinates[,2]))

init_densityCentromer <- compute.density.centromer(init_funnelTable)

# First iteration

first_gscoreLogPath = '/home/guilherme/datasets/SALB3.NOCST.BISCORE_PROQ3_TM.002/gscore/gscore-TMscore-050.dat'
first_loglistLocation = '/home/guilherme/datasets/SALB3.NOCST.BISCORE_PROQ3_TM.002/loglist.txt'
first_topolinkLogsDirectory = '/home/guilherme/datasets/SALB3.NOCST.BISCORE_PROQ3_TM.002/topolink_observed/'
first_alignlistPath = '/home/guilherme/datasets/SALB3.NOCST.BISCORE_PROQ3_TM.002/alignlist.txt'
first_compactlogPath = '/home/guilherme/datasets/SALB3.NOCST.BISCORE_PROQ3_TM.002/gscore/compactlog-TMscore.dat'
first_lovoalignLogPath = '/home/guilherme/datasets/SALB3.NOCST.BISCORE_PROQ3_TM.002/lovoalign.log'
first_proq3listLocation = '/home/guilherme/datasets/SALB3.NOCST.BISCORE_PROQ3_TM.002/proq3list.txt'
first_cryslistLocation = '/home/guilherme/datasets/SALB3.NOCST.BISCORE_PROQ3_TM.002/utility/cryslist.txt'
first_optlistLocation = '/home/guilherme/datasets/SALB3.NOCST.BISCORE_PROQ3_TM.002/utility/optlist.txt'
first_distanceTableLocation = '/home/guilherme/datasets/SALB3.NOCST.BISCORE_PROQ3_TM.002/utility/distance_table.log'

first_modelScores <- compute.model.scores(type = "gscore",
                                         first_gscoreLogPath,
                                         first_lovoalignLogPath,
                                         first_proq3listLocation,
                                         computeproq3)

first_optimumXlinkMirttable <- create.xlink.mirttable(
  prepare.topolink.logs(
    read.topolink.output(mode,
                         first_loglistLocation,
                         first_topolinkLogsDirectory)
  )
)

first_optimumSimilarityTable <- create.dissimilarity.matrix(mode = "similarity",
                                                           diagonal = 1,
                                                           first_alignlistPath,
                                                           first_compactlogPath,
                                                           nmodels = nrow(first_optimumXlinkMirttable))

first_optimumXlinkMirttable <- digest.xlink.mirttable(first_optimumXlinkMirttable)

first_optimumDescript <- ltm::descript(first_optimumXlinkMirttable)

first_globalCronbachAlpha <- first_optimumDescript$alpha[1]

first_restrictionScores <- data.frame(bis = first_optimumDescript$bisCorr,
                                     right = first_optimumDescript$perc[,2],
                                     logit = first_optimumDescript$perc[,3],
                                     alpha = first_optimumDescript$alpha[-1],
                                     deltaAlpha = first_optimumDescript$alpha[2:length(first_optimumDescript$alpha)] - first_globalCronbachAlpha)

first_restrictionScores <- attribute.crys.and.opt(first_restrictionScores, cryslistLocation, optlistLocation)

first_modelDistanceTable <- 1/(first_optimumSimilarityTable)

first_projection.coordinates <- forceScheme(first_modelDistanceTable,
                                           Y = NULL,
                                           #max.iter = ,
                                           tol = 1e-04,
                                           #fraction = 8,
                                           eps = 1e-05)

first_projection.coordinates <- scale(first_projection.coordinates, center = T, scale = T)

first_funnelTable <- data.frame("projectionX" = first_projection.coordinates[,1],
                               "projectionY" = first_projection.coordinates[,2],
                               "tmscore" = first_modelScores$`TM-Score`,
                               "density" = get_density(x = first_projection.coordinates[,1], y = first_projection.coordinates[,2]))

first_densityCentromer <- compute.density.centromer(first_funnelTable)

# Second Iteration

second_gscoreLogPath = '/home/guilherme/datasets/SALB3.NOCST.BISCORE_PROQ3_TM.003/gscore/gscore-TMscore-050.dat'
second_loglistLocation = '/home/guilherme/datasets/SALB3.NOCST.BISCORE_PROQ3_TM.003/loglist.txt'
second_topolinkLogsDirectory = '/home/guilherme/datasets/SALB3.NOCST.BISCORE_PROQ3_TM.003/topolink_observed/'
second_alignlistPath = '/home/guilherme/datasets/SALB3.NOCST.BISCORE_PROQ3_TM.003/alignlist.txt'
second_compactlogPath = '/home/guilherme/datasets/SALB3.NOCST.BISCORE_PROQ3_TM.003/gscore/compactlog-TMscore.dat'
second_lovoalignLogPath = '/home/guilherme/datasets/SALB3.NOCST.BISCORE_PROQ3_TM.003/lovoalign.log'
second_proq3listLocation = '/home/guilherme/datasets/SALB3.NOCST.BISCORE_PROQ3_TM.003/proq3list.txt'
second_cryslistLocation = '/home/guilherme/datasets/SALB3.NOCST.BISCORE_PROQ3_TM.003/utility/cryslist.txt'
second_optlistLocation = '/home/guilherme/datasets/SALB3.NOCST.BISCORE_PROQ3_TM.003/utility/optlist.txt'
second_distanceTableLocation = '/home/guilherme/datasets/SALB3.NOCST.BISCORE_PROQ3_TM.003/utility/distance_table.log'

second_modelScores <- compute.model.scores(type = "gscore",
                                          second_gscoreLogPath,
                                          second_lovoalignLogPath,
                                          second_proq3listLocation,
                                          computeproq3)

second_optimumXlinkMirttable <- create.xlink.mirttable(
  prepare.topolink.logs(
    read.topolink.output(mode,
                         second_loglistLocation,
                         second_topolinkLogsDirectory)
  )
)

second_optimumSimilarityTable <- create.dissimilarity.matrix(mode = "similarity",
                                                            diagonal = 1,
                                                            second_alignlistPath,
                                                            second_compactlogPath,
                                                            nmodels = nrow(second_optimumXlinkMirttable))

second_optimumXlinkMirttable <- digest.xlink.mirttable(second_optimumXlinkMirttable)

second_optimumDescript <- ltm::descript(second_optimumXlinkMirttable)

second_globalCronbachAlpha <- second_optimumDescript$alpha[1]

second_restrictionScores <- data.frame(bis = second_optimumDescript$bisCorr,
                                      right = second_optimumDescript$perc[,2],
                                      logit = second_optimumDescript$perc[,3],
                                      alpha = second_optimumDescript$alpha[-1],
                                      deltaAlpha = second_optimumDescript$alpha[2:length(second_optimumDescript$alpha)] - second_globalCronbachAlpha)

second_restrictionScores <- attribute.crys.and.opt(second_restrictionScores, cryslistLocation, optlistLocation)

second_modelDistanceTable <- 1/(second_optimumSimilarityTable)

second_projection.coordinates <- forceScheme(second_modelDistanceTable,
                                            Y = NULL,
                                            #max.iter = ,
                                            tol = 1e-04,
                                            #fraction = 8,
                                            eps = 1e-05)

second_projection.coordinates <- scale(second_projection.coordinates, center = T, scale = T)

second_funnelTable <- data.frame("projectionX" = second_projection.coordinates[,1],
                                "projectionY" = second_projection.coordinates[,2],
                                "tmscore" = second_modelScores$`TM-Score`,
                                "density" = get_density(x = second_projection.coordinates[,1], y = second_projection.coordinates[,2]))

second_densityCentromer <- compute.density.centromer(second_funnelTable)

# Third Iteration

third_gscoreLogPath = '/home/guilherme/datasets/SALB3.NOCST.BISCORE_PROQ3_TM.004/gscore/gscore-TMscore-050.dat'
third_loglistLocation = '/home/guilherme/datasets/SALB3.NOCST.BISCORE_PROQ3_TM.004/loglist.txt'
third_topolinkLogsDirectory = '/home/guilherme/datasets/SALB3.NOCST.BISCORE_PROQ3_TM.004/topolink_observed/'
third_alignlistPath = '/home/guilherme/datasets/SALB3.NOCST.BISCORE_PROQ3_TM.004/alignlist.txt'
third_compactlogPath = '/home/guilherme/datasets/SALB3.NOCST.BISCORE_PROQ3_TM.004/gscore/compactlog-TMscore.dat'
third_lovoalignLogPath = '/home/guilherme/datasets/SALB3.NOCST.BISCORE_PROQ3_TM.004/lovoalign.log'
third_proq3listLocation = '/home/guilherme/datasets/SALB3.NOCST.BISCORE_PROQ3_TM.004/proq3list.txt'
third_cryslistLocation = '/home/guilherme/datasets/SALB3.NOCST.BISCORE_PROQ3_TM.004/utility/cryslist.txt'
third_optlistLocation = '/home/guilherme/datasets/SALB3.NOCST.BISCORE_PROQ3_TM.004/utility/optlist.txt'
third_distanceTableLocation = '/home/guilherme/datasets/SALB3.NOCST.BISCORE_PROQ3_TM.004/utility/distance_table.log'

third_modelScores <- compute.model.scores(type = "gscore",
                                           third_gscoreLogPath,
                                           third_lovoalignLogPath,
                                           third_proq3listLocation,
                                           computeproq3)

third_optimumXlinkMirttable <- create.xlink.mirttable(
  prepare.topolink.logs(
    read.topolink.output(mode,
                         third_loglistLocation,
                         third_topolinkLogsDirectory)
  )
)

third_optimumSimilarityTable <- create.dissimilarity.matrix(mode = "similarity",
                                                             diagonal = 1,
                                                             third_alignlistPath,
                                                             third_compactlogPath,
                                                             nmodels = nrow(third_optimumXlinkMirttable))

third_optimumXlinkMirttable <- digest.xlink.mirttable(third_optimumXlinkMirttable)

third_optimumDescript <- ltm::descript(third_optimumXlinkMirttable)

third_globalCronbachAlpha <- third_optimumDescript$alpha[1]

third_restrictionScores <- data.frame(bis = third_optimumDescript$bisCorr,
                                       right = third_optimumDescript$perc[,2],
                                       logit = third_optimumDescript$perc[,3],
                                       alpha = third_optimumDescript$alpha[-1],
                                       deltaAlpha = third_optimumDescript$alpha[2:length(third_optimumDescript$alpha)] - third_globalCronbachAlpha)

third_restrictionScores <- attribute.crys.and.opt(third_restrictionScores, cryslistLocation, optlistLocation)

third_modelDistanceTable <- 1/(third_optimumSimilarityTable)

third_projection.coordinates <- forceScheme(third_modelDistanceTable,
                                             Y = NULL,
                                             #max.iter = ,
                                             tol = 1e-04,
                                             #fraction = 8,
                                             eps = 1e-05)

third_projection.coordinates <- scale(third_projection.coordinates, center = T, scale = T)

third_funnelTable <- data.frame("projectionX" = third_projection.coordinates[,1],
                                 "projectionY" = third_projection.coordinates[,2],
                                 "tmscore" = third_modelScores$`TM-Score`,
                                 "density" = get_density(x = third_projection.coordinates[,1], y = third_projection.coordinates[,2]))

third_densityCentromer <- compute.density.centromer(third_funnelTable)

# Fourth Iteration

fourth_gscoreLogPath = '/home/guilherme/datasets/SALB3.NOCST.BISCORE_PROQ3_TM.005/gscore/gscore-TMscore-050.dat'
fourth_loglistLocation = '/home/guilherme/datasets/SALB3.NOCST.BISCORE_PROQ3_TM.005/loglist.txt'
fourth_topolinkLogsDirectory = '/home/guilherme/datasets/SALB3.NOCST.BISCORE_PROQ3_TM.005/topolink_observed/'
fourth_alignlistPath = '/home/guilherme/datasets/SALB3.NOCST.BISCORE_PROQ3_TM.005/alignlist.txt'
fourth_compactlogPath = '/home/guilherme/datasets/SALB3.NOCST.BISCORE_PROQ3_TM.005/gscore/compactlog-TMscore.dat'
fourth_lovoalignLogPath = '/home/guilherme/datasets/SALB3.NOCST.BISCORE_PROQ3_TM.005/lovoalign.log'
fourth_proq3listLocation = '/home/guilherme/datasets/SALB3.NOCST.BISCORE_PROQ3_TM.005/proq3list.txt'

fourth_modelScores <- compute.model.scores(type = "gscore",
                                          fourth_gscoreLogPath,
                                          fourth_lovoalignLogPath,
                                          fourth_proq3listLocation,
                                          computeproq3)

fourth_optimumXlinkMirttable <- create.xlink.mirttable(
  prepare.topolink.logs(
    read.topolink.output(mode,
                         fourth_loglistLocation,
                         fourth_topolinkLogsDirectory)
  )
)

fourth_optimumSimilarityTable <- create.dissimilarity.matrix(mode = "similarity",
                                                            diagonal = 1,
                                                            fourth_alignlistPath,
                                                            fourth_compactlogPath,
                                                            nmodels = nrow(fourth_optimumXlinkMirttable))

fourth_optimumXlinkMirttable <- digest.xlink.mirttable(fourth_optimumXlinkMirttable)

fourth_optimumDescript <- ltm::descript(fourth_optimumXlinkMirttable)

fourth_globalCronbachAlpha <- fourth_optimumDescript$alpha[1]

fourth_restrictionScores <- data.frame(bis = fourth_optimumDescript$bisCorr,
                                      right = fourth_optimumDescript$perc[,2],
                                      logit = fourth_optimumDescript$perc[,3],
                                      alpha = fourth_optimumDescript$alpha[-1],
                                      deltaAlpha = fourth_optimumDescript$alpha[2:length(fourth_optimumDescript$alpha)] - fourth_globalCronbachAlpha)

fourth_restrictionScores <- attribute.crys.and.opt(fourth_restrictionScores, cryslistLocation, optlistLocation)

fourth_modelDistanceTable <- 1/(fourth_optimumSimilarityTable)

fourth_projection.coordinates <- forceScheme(fourth_modelDistanceTable,
                                            Y = NULL,
                                            #max.iter = ,
                                            tol = 1e-04,
                                            #fraction = 8,
                                            eps = 1e-05)

fourth_projection.coordinates <- scale(fourth_projection.coordinates, center = T, scale = T)

fourth_funnelTable <- data.frame("projectionX" = fourth_projection.coordinates[,1],
                                "projectionY" = fourth_projection.coordinates[,2],
                                "tmscore" = fourth_modelScores$`TM-Score`,
                                "density" = get_density(x = fourth_projection.coordinates[,1], y = fourth_projection.coordinates[,2]))

fourth_densityCentromer <- compute.density.centromer(fourth_funnelTable)

# Fifth iteration

fifth_gscoreLogPath = '/home/guilherme/datasets/SALB3.NOCST.BISCORE_PROQ3_TM.006/gscore/gscore-TMscore-050.dat'
fifth_loglistLocation = '/home/guilherme/datasets/SALB3.NOCST.BISCORE_PROQ3_TM.006/loglist.txt'
fifth_topolinkLogsDirectory = '/home/guilherme/datasets/SALB3.NOCST.BISCORE_PROQ3_TM.006/topolink_observed/'
fifth_alignlistPath = '/home/guilherme/datasets/SALB3.NOCST.BISCORE_PROQ3_TM.006/alignlist.txt'
fifth_compactlogPath = '/home/guilherme/datasets/SALB3.NOCST.BISCORE_PROQ3_TM.006/gscore/compactlog-TMscore.dat'
fifth_lovoalignLogPath = '/home/guilherme/datasets/SALB3.NOCST.BISCORE_PROQ3_TM.006/lovoalign.log'
fifth_proq3listLocation = '/home/guilherme/datasets/SALB3.NOCST.BISCORE_PROQ3_TM.006/proq3list.txt'

fifth_modelScores <- compute.model.scores(type = "gscore",
                                           fifth_gscoreLogPath,
                                           fifth_lovoalignLogPath,
                                           fifth_proq3listLocation,
                                           computeproq3)

fifth_optimumXlinkMirttable <- create.xlink.mirttable(
  prepare.topolink.logs(
    read.topolink.output(mode,
                         fifth_loglistLocation,
                         fifth_topolinkLogsDirectory)
  )
)

fifth_optimumSimilarityTable <- create.dissimilarity.matrix(mode = "similarity",
                                                             diagonal = 1,
                                                             fifth_alignlistPath,
                                                             fifth_compactlogPath,
                                                             nmodels = nrow(fifth_optimumXlinkMirttable))

fifth_optimumXlinkMirttable <- digest.xlink.mirttable(fifth_optimumXlinkMirttable)

fifth_optimumDescript <- ltm::descript(fifth_optimumXlinkMirttable)

fifth_globalCronbachAlpha <- fifth_optimumDescript$alpha[1]

fifth_restrictionScores <- data.frame(bis = fifth_optimumDescript$bisCorr,
                                       right = fifth_optimumDescript$perc[,2],
                                       logit = fifth_optimumDescript$perc[,3],
                                       alpha = fifth_optimumDescript$alpha[-1],
                                       deltaAlpha = fifth_optimumDescript$alpha[2:length(fifth_optimumDescript$alpha)] - fifth_globalCronbachAlpha)

fifth_restrictionScores <- attribute.crys.and.opt(fifth_restrictionScores, cryslistLocation, optlistLocation)

fifth_modelDistanceTable <- 1/(fifth_optimumSimilarityTable)

fifth_projection.coordinates <- forceScheme(fifth_modelDistanceTable,
                                             Y = NULL,
                                             #max.iter = ,
                                             tol = 1e-04,
                                             #fraction = 8,
                                             eps = 1e-05)

fifth_projection.coordinates <- scale(fifth_projection.coordinates, center = T, scale = T)

fifth_funnelTable <- data.frame("projectionX" = fifth_projection.coordinates[,1],
                                 "projectionY" = fifth_projection.coordinates[,2],
                                 "tmscore" = fifth_modelScores$`TM-Score`,
                                 "density" = get_density(x = fifth_projection.coordinates[,1], y = fifth_projection.coordinates[,2]))

fifth_densityCentromer <- compute.density.centromer(fifth_funnelTable)

## Collective Plots

### Tm_Score

# Initial Modeling:

ggplot(init_funnelTable[order(init_funnelTable$tmscore),], aes(x = projectionX, y = projectionY)) +
  geom_point(aes(color = tmscore), alpha = 0.5, size = 2.5) +
  xlim(-3, 3) + ylim(-3,3) +
  annotate("point", x = init_densityCentromer["centromerX"], y = init_densityCentromer["centromerY"], size = 10,color = "orange", alpha = 0.5) +
  scale_colour_viridis_c(limits = c(0.1, 0.8), direction = -1) + theme_bw()

# First Iteration:

ggplot(first_funnelTable[order(first_funnelTable$tmscore),], aes(x = projectionX, y = projectionY)) +
  geom_point(aes(color = tmscore), alpha = 0.5, size = 2.5) +
  xlim(-3, 3) + ylim(-3,3) +
  annotate("point", x = first_densityCentromer["centromerX"], y = first_densityCentromer["centromerY"], size = 10,color = "orange", alpha = 0.5) +
  scale_colour_viridis_c(limits = c(0.1, 0.8), direction = -1) + theme_bw()

# Second Iteration:

ggplot(second_funnelTable[order(second_funnelTable$tmscore),], aes(x = projectionX, y = projectionY)) +
  geom_point(aes(color = tmscore), alpha = 0.5, size = 2.5) +
  xlim(-3, 3) + ylim(-3,3) +
  annotate("point", x = second_densityCentromer["centromerX"], y = second_densityCentromer["centromerY"], size = 10,color = "orange", alpha = 0.5) +
  scale_colour_viridis_c(limits = c(0.1, 0.8), direction = -1) + theme_bw()

# Third Iteration:

ggplot(third_funnelTable[order(third_funnelTable$tmscore),], aes(x = projectionX, y = projectionY)) +
  geom_point(aes(color = tmscore), alpha = 0.5, size = 2.5) +
  xlim(-3, 3) + ylim(-3,3) +
  annotate("point", x = third_densityCentromer["centromerX"], y = third_densityCentromer["centromerY"], size = 10,color = "orange", alpha = 0.5) +
  scale_colour_viridis_c(limits = c(0.1, 0.8), direction = -1) + theme_bw()

# Fourth Iteration:

ggplot(fourth_funnelTable[order(fourth_funnelTable$tmscore),], aes(x = projectionX, y = projectionY)) +
  geom_point(aes(color = tmscore), alpha = 0.5, size = 2.5) +
  xlim(-3, 3) + ylim(-3,3) +
  annotate("point", x = fourth_densityCentromer["centromerX"], y = fourth_densityCentromer["centromerY"], size = 10,color = "orange", alpha = 0.5) +
  scale_colour_viridis_c(limits = c(0.1, 0.8), direction = -1) + theme_bw()

# Fifth Iteration:

ggplot(fifth_funnelTable[order(fifth_funnelTable$tmscore),], aes(x = projectionX, y = projectionY)) +
  geom_point(aes(color = tmscore), alpha = 0.5, size = 2.5) +
  xlim(-3, 3) + ylim(-3,3) +
  annotate("point", x = fifth_densityCentromer["centromerX"], y = fifth_densityCentromer["centromerY"], size = 10,color = "orange", alpha = 0.5) +
  scale_colour_viridis_c(limits = c(0.1, 0.8), direction = -1) + theme_bw()

## Combination of all modelings:

combined_funnelTable <- rbind(init_funnelTable,
                              first_funnelTable,
                              second_funnelTable,
                              third_funnelTable,
                              fourth_funnelTable,
                              fifth_funnelTable)

combined_funnelTable$density <- get_density(x = combined_funnelTable$projectionX, y = combined_funnelTable$projectionY)

ggplot(combined_funnelTable[order(combined_funnelTable$tmscore),], aes(x = projectionX, y = projectionY)) +
  geom_point(aes(color = tmscore), alpha = 0.2, size = 2.5) +
  xlim(-3, 3) + ylim(-3,3) +
  annotate("point", x = fifth_densityCentromer["centromerX"], y = fifth_densityCentromer["centromerY"], size = 10,color = "orange", alpha = 0.5) +
  scale_colour_viridis_c(limits = c(0.1, 0.8), direction = -1) + theme_bw()

ggplot(combined_funnelTable, aes(x = projectionX, y = projectionY)) +
  geom_point(aes(color = density), alpha = 0.2, size = 2.5) +
  xlim(-3, 3) + ylim(-3,3) +
  annotate("point", x = fifth_densityCentromer["centromerX"], y = fifth_densityCentromer["centromerY"], size = 10,color = "orange", alpha = 0.5) +
  scale_colour_viridis_c(direction = -1) + theme_bw()

### Density

# Initial Modeling:

ggplot(init_funnelTable, aes(x = projectionX, y = projectionY)) +
  geom_point(aes(color = density), alpha = 0.5, size = 2.5) +
  xlim(-3, 3) + ylim(-3,3) +
  annotate("point", x = init_densityCentromer["centromerX"], y = init_densityCentromer["centromerY"], size = 10,color = "orange", alpha = 0.5) +
  scale_colour_viridis_c(direction = -1) + theme_bw()

# First Iteration:
ggplot(first_funnelTable, aes(x = projectionX, y = projectionY)) +
  geom_point(aes(color = density), alpha = 0.5, size = 2.5) +
  xlim(-3, 3) + ylim(-3,3) +
  annotate("point", x = first_densityCentromer["centromerX"], y = first_densityCentromer["centromerY"], size = 10,color = "orange", alpha = 0.5) +
  scale_colour_viridis_c(direction = -1) + theme_bw()

# First Iteration:
ggplot(second_funnelTable, aes(x = projectionX, y = projectionY)) +
  geom_point(aes(color = density), alpha = 0.5, size = 2.5) +
  xlim(-3, 3) + ylim(-3,3) +
  annotate("point", x = second_densityCentromer["centromerX"], y = second_densityCentromer["centromerY"], size = 10,color = "orange", alpha = 0.5) +
  scale_colour_viridis_c(direction = -1) + theme_bw()

ggplot(third_funnelTable, aes(x = projectionX, y = projectionY)) +
  geom_point(aes(color = density), alpha = 0.9, size = 3) +
  xlim(-3, 3) + ylim(-3,3) +
  annotate("point", x = third_densityCentromer["centromerX"], y = third_densityCentromer["centromerY"], size = 10,color = "orange", alpha = 0.5) +
  scale_colour_viridis_c(direction = -1) + theme_bw()

ggplot(fourth_funnelTable, aes(x = projectionX, y = projectionY)) +
  geom_point(aes(color = density), alpha = 0.5, size = 2.5) +
  xlim(-3, 3) + ylim(-3,3) +
  annotate("point", x = fourth_densityCentromer["centromerX"], y = fourth_densityCentromer["centromerY"], size = 10,color = "orange", alpha = 0.5) +
  scale_colour_viridis_c(direction = -1) + theme_bw()

ggplot(fifth_funnelTable, aes(x = projectionX, y = projectionY)) +
  geom_point(aes(color = density), alpha = 0.5, size = 2.5) +
  xlim(-3, 3) + ylim(-3,3) +
  annotate("point", x = fifth_densityCentromer["centromerX"], y = fifth_densityCentromer["centromerY"], size = 10,color = "orange", alpha = 0.5) +
  scale_colour_viridis_c(direction = -1) + theme_bw()


### Baripapa

bigfunnel <- read.table('/home/guilherme/logs/2018-10-07/coordinates.txt')
bigfunnelScores <- read.table('/home/guilherme/logs/2018-10-07/lovoalign.log')
bigfunnel <- scale(bigfunnel, center = T, scale = T)

bigfunnelTable <- data.frame("projectionX" = bigfunnel[,1],
                                "projectionY" = bigfunnel[,2],
                                "tmscore" = bigfunnelScores$V3,
                                "density" = get_density(x = bigfunnel[,1], y = bigfunnel[,2]))

bigfunnelCentromer <- compute.density.centromer(bigfunnelTable)

ggplot(bigfunnelTable[order(bigfunnelTable$tmscore),], aes(x = projectionX, y = projectionY)) +
  geom_point(aes(color = tmscore), alpha = 0.5, size = 2.5) +
  xlim(-3, 3) + ylim(-3,3) +
  annotate("point", x = bigfunnelCentromer["centromerX"], y = bigfunnelCentromer["centromerY"], size = 10,color = "orange", alpha = 0.5) +
  scale_colour_viridis_c(limits = c(0.1, 0.85), direction = -1) + theme_bw()

ggplot(bigfunnelTable[order(bigfunnelTable$density),], aes(x = projectionX, y = projectionY)) +
  geom_point(aes(color = density), alpha = 0.5, size = 2.5) +
  xlim(-3, 3) + ylim(-3,3) +
  annotate("point", x = bigfunnelCentromer["centromerX"], y = bigfunnelCentromer["centromerY"], size = 10,color = "orange", alpha = 0.5) +
  scale_colour_viridis_c(direction = -1) + theme_bw()

library(plot3D)
library(plotly)

scatter3D(bigfunnelTable$projectionX,
          bigfunnelTable$projectionY,
          -bigfunnelTable$tmscore, phi = 0, bty ="b2")

plot_ly(bigfunnelTable,
        x = ~projectionX,
        y = ~projectionY,
        z = ~-tmscore,
        marker=list(
          color=~-tmscore,
          colorbar=list(
            title='Colorbar'
          ),
          colorscale='viridis',
          reversescale = F
        )) %>%
  add_markers()

# Animation

animationTable <- bigfunnelTable
animationTable$frame <- c(rep(1, 5000), rep(2, 5000), rep(3, 5000), rep(4, 5000), rep(5, 5000), rep(6, 5000))

scatter3D(animationTable[animationTable$frame == 1, ]$projectionX,
          animationTable[animationTable$frame == 1, ]$projectionY,
          -animationTable[animationTable$frame == 1, ]$tmscore, phi = 0, bty ="b2",
          xlab = "Projection X", ylab = "Projection Y", zlab = "-(TM-Score)",
          col = viridis(100))

scatter3D(animationTable[animationTable$frame == 2, ]$projectionX,
          animationTable[animationTable$frame == 2, ]$projectionY,
          -animationTable[animationTable$frame == 2, ]$tmscore, phi = 0, bty ="b2",
          xlab = "Projection X", ylab = "Projection Y", zlab = "-(TM-Score)",
          col = viridis(100))

scatter3D(animationTable[animationTable$frame == 3, ]$projectionX,
          animationTable[animationTable$frame == 3, ]$projectionY,
          -animationTable[animationTable$frame == 3, ]$tmscore, phi = 0, bty ="b2",
          xlab = "Projection X", ylab = "Projection Y", zlab = "-(TM-Score)",
          col = viridis(100))

scatter3D(animationTable[animationTable$frame == 4, ]$projectionX,
          animationTable[animationTable$frame == 4, ]$projectionY,
          -animationTable[animationTable$frame == 4, ]$tmscore, phi = 0, bty ="b2",
          xlab = "Projection X", ylab = "Projection Y", zlab = "-(TM-Score)",
          col = viridis(100))

scatter3D(animationTable[animationTable$frame == 5, ]$projectionX,
          animationTable[animationTable$frame == 5, ]$projectionY,
          -animationTable[animationTable$frame == 5, ]$tmscore, phi = 0, bty ="b2",
          xlab = "Projection X", ylab = "Projection Y", zlab = "-(TM-Score)",
          col = viridis(100))

scatter3D(animationTable[animationTable$frame == 6, ]$projectionX,
          animationTable[animationTable$frame == 6, ]$projectionY,
          -animationTable[animationTable$frame == 6, ]$tmscore, phi = 0, bty ="b2",
          xlab = "Projection X", ylab = "Projection Y", zlab = "-(TM-Score)",
          col = viridis(100))

scatter3D(bigfunnelTable$projectionX,
          bigfunnelTable$projectionY,
          -bigfunnelTable$tmscore, phi = 0, bty ="b2",
          xlab = "Projection X", ylab = "Projection Y", zlab = "-(TM-Score)",
          col = viridis(100))

animationTable <- do.call(rbind, animationTable)

funnelPlot <- plot_ly(animationTable,
                      x = ~projectionX,
                      y = ~projectionY,
                      z = ~-tmscore,
                      frame = ~frame,
                      color = ~-tmscore,
                      colors = viridis(10)) %>%
  add_markers() %>%
  animation_opts(frame = 5000,
                 transition = 2000,
                 easing = "elastic",
                 redraw = F)

funnelPlot

Sys.setenv("plotly_username"="GuilhermeFahurBottino")
Sys.setenv("plotly_api_key"="le7AUsxZNFDYFKUPEQ0Y")

chart_link = api_create(funnelPlot, filename="animations-animation-options")

htmlwidgets::saveWidget(widget=funnelPlot,"index.html")

library(plotly)
library(htmlwidgets)
