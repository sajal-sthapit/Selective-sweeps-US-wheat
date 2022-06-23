#Spillover code, in case I need it later. Delete after checking completed

# snp.types and snp.order ----
# numeric codes may not be needed for unpolarized analysis.
# might want to know snp.types for plotting post analysis.
snp.types <- read_tsv("data/90k_SNP_type.txt") %>% dplyr::select(Name, SNPTYPE)
snp.order <- left_join(geno$nuc[1:6], snp.types, by = "Name")
#subset of snp.types for 10391 snps in the same order as in g.poly

# create numeric codes for AA and BB based on snp types, as required for rehh.
AA <- substr(snp.order$SNPTYPE, 1, 1)
AA[AA == "A"] <- 1
AA[AA == "C"] <- 2
AA[AA == "G"] <- 3
AA[AA == "T"] <- 4
AA <- as.numeric(AA)
BB <- substr(snp.order$SNPTYPE, 2, 2)
BB[BB == "A"] <- 1
BB[BB == "C"] <- 2
BB[BB == "G"] <- 3
BB[BB == "T"] <- 4
BB <- as.numeric(BB)

snp.order <- mutate(snp.order, AA = AA, BB = BB)
head(snp.order)
rm(AA, BB)