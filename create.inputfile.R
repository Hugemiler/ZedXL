#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

create.inputfile <- function(modelingDirectory, nconst, computeproq3, indicator) {

  inputFileBody <- NULL

  inputFileBody[1]  <- paste0('nconst = ', nconst)
  inputFileBody[2]  <- paste0('computeproq3 = ', computeproq3)
  inputFileBody[3]  <- paste0('indicator = \'', indicator, '\' \n')

  inputFileBody[4]  <- paste0('gscoreLogPath = \'', modelingDirectory, '/gscore/gscore-TMscore-050.dat\'')
  inputFileBody[5]  <- paste0('loglistLocation = \'', modelingDirectory, '/loglist.txt\'')
  inputFileBody[6]  <- paste0('topolinkLogsDirectory = \'', modelingDirectory, '/topolink_observed/\'')

  inputFileBody[7]  <- paste0('alignlistPath = \'', modelingDirectory, '/alignlist.txt\'')
  inputFileBody[8]  <- paste0('compactlogPath = \'', modelingDirectory, '/gscore/compactlog-TMscore.dat\'')
  inputFileBody[9]  <- paste0('lovoalignLogPath = \'', modelingDirectory, '/lovoalign.log\'')

  inputFileBody[10] <- paste0('proq3listLocation = \'', modelingDirectory, '/proq3list.txt\'')

  inputFileBody[11] <- paste0('cryslistLocation = \'', modelingDirectory, '/utility/cryslist.txt\'')
  inputFileBody[12] <- paste0('optlistLocation = \'', modelingDirectory, '/utility/optlist.txt\'')
  inputFileBody[13] <- paste0('distanceTableLocation = \'', modelingDirectory, '/utility/distance_table.log\'')

  return(inputFileBody)

}

write(create.inputfile(modelingDirectory = args[1],
                       nconst = args[2],
                       computeproq3 = args[3],
                       indicator = args[4]),
      file = args[5])
