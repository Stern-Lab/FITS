# toBi.R
# Written by Tal Zinger
# 
# This script converts dataset compatible with FITS to biallelic.
# The WT allele is defined as the most common at the first generation.
# the frequency of the WT allele is then taken for all other generations,
# and the mutant allele frequency complements the WT to 1.
#
# in order to install the package - run the command:
# install.packages("plyr")

library(plyr)

args.vec <- commandArgs(trailingOnly = T)

if ( length(args.vec) != 2 ) {
	print("Syntax is Rscript toBi.R input_file.txt output_file.txt")
	quit()
}

in.file.name <- args.vec[1]
out.file.name <- args.vec[2]


quad.df <- read.table( in.file.name, header = T )

quad.df <- quad.df[ order(quad.df$gen), ]


# search for the WT allele
first.gen <- quad.df$gen[1]
first.gen.df <- quad.df[ quad.df$gen==first.gen, ]
first.gen.df <- first.gen.df[ order(-first.gen.df$freq), ]
wt.allele <- first.gen.df$allele[1]

wf.allele.df <- quad.df[ quad.df$allele==wt.allele, ]
wf.allele.df$allele <- 0

bi.df <- ddply( wf.allele.df, "gen", function(current.generation.df) {
	
	mutant.freq <- 1 - current.generation.df$freq
	
	tmp.df <- data.frame( gen=current.generation.df$gen, allele=1, freq=mutant.freq)
	
	rbind( current.generation.df, tmp.df)
})


write.table( bi.df, file = out.file.name, quote = F, row.names = F, sep = "\t" )