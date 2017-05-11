#print-args.R this script 

args <- commandArgs(trailingOnly = TRUE) #return command line arguments as a vector
cat(args, sep = "\n") #send output to standard output, print vector element on own line
