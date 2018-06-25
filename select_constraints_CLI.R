#!/usr/bin/env Rscript

## Loading Package
library(ZedXL)

## Read Argument passed via Command-Line
args <- commandArgs(trailingOnly=TRUE)

print(args[1])
print(args[2])

## Run Analysis
select.constraints.CLI(inputfile = args[1], outputfile = as.character(args[2]))
