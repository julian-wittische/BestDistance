# I strongly suspect that there is an issue stemming from the code from Winiarski et al. 
# The df2genind(cd_df) part of their code, where they load the cdpop output grid 
# does not put out an error but incorrectly interprets allele counts as allele IDs. 
# This creates a functional object with an incorrect number of loci and weird allelic
# frequencies. If they used 2 alleles, they may have miraculously avoided the issue.

alleles <- 10
loci <- 16

# First let's get rid of the comma at the last column (ISSUE: forces character)
cd_df[,ncol(cd_df)] <- gsub(",","",cd_df[,ncol(cd_df)])
cd_df <- apply(as.matrix(cd_df),2,as.numeric)


fakedf <- data.frame(matrix(rep(paste(sample(1:alleles,nrow(cd_df),replace=TRUE),
                                      sample(1:alleles,nrow(cd_df),replace=TRUE),
                                      sep="/"), loci), ncol=loci)) 

colnames(fakedf) <- .genlab("Locus",loci)

lool <- df2genind(fakedf, ploidy=2, sep="/")
lool@tab <- cd_df




