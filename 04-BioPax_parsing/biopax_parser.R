
library("rBiopaxParser")
library("argparse")

parser <- ArgumentParser()

parser$add_argument("-i", "--input",  
    help="Number of random normals to generate [default %(default)s]")

parser$add_argument("-o", "--output",  
    help="Number of random normals to generate [default %(default)s]")
# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults, 
args = parser$parse_args()

biopax = readBiopax(args$input)
biopax$dt$property_value = gsub("[\r\n]", " ", biopax$dt$property_value)
write.table(biopax$dt, file = args$output, quote = F, row.names = F, sep = '\t')