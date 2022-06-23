# function definitions for all analysis for selective sweep in the US regional wheat populations.

calc_major_minor_alleles_df <- function(lines, geno, snp.order){
  # lines = vector of variety names in the population to be assessed
  # geno = genotype calls with variety names as columns, markers as rows
  # snp.order = data.frame with SNPs in the same order as geno. The function only uses columns SNPTYPE, AA and BB for calculations. 
  # But all the columns are used to save output with marker details.
  g.inp <- dplyr::select(geno, all_of(lines))     # genotype of specified lines
  alleles <- tibble(num.A = rowSums(g.inp == "A"),  # allele counts in specified lines
                    num.C = rowSums(g.inp == "C"),
                    num.G = rowSums(g.inp == "G"),
                    num.T = rowSums(g.inp == "T"))
  
  input <- cbind(snp.order, alleles)
  
  allele.freq <- list(0)
  for(i in 1:nrow(alleles)){
    allele.freq[i] <- list(c(alleles$num.A[i], alleles$num.C[i], alleles$num.G[i], alleles$num.T[i]) )
  }
  
  # function definition
  determine_major_allele <- function(inp = c(35,0,5,0), AA = 1, BB = 3){
    # function to determine major allele. 
    # Input (inp) is a vector of 4 allele frequencies in the order of ACGT
    # AA and BB specify the two alleles for the SNP, where A = 1, C = 2, G = 3, T = 4. 
    # Therefore, C/G SNP is defined as AA = 2, BB =3.
    # The output is going to be 0 for indeterminate (both alleles in equal frequency) or 2 or 3.
    if_else(inp[AA] == inp[BB], 0, # 0 means both alleles are in equal frequency
            if_else(inp[AA] > inp[BB], AA, BB) )
  }
  
  # identify major allele  
  major.allele <- rep(NA, nrow(alleles)) # create empty vector
  for(i in 1:nrow(alleles)){
    major.allele[i] <- determine_major_allele(inp = allele.freq[[i]], input$AA[i], input$BB[i])
  }
  
  # identify minor allele as the complementary allele to major allele from the SNPTYPE information
  A <- as.character(major.allele)
  A[A == "1"] <- "A"
  A[A == "2"] <- "C"
  A[A == "3"] <- "G"
  A[A == "4"] <- "T"
  snp.list <- map(.x = snp.order$SNPTYPE, .f = ~ c(substr(.x, 1, 1), substr(.x, 2, 2)))
  B <- map2_chr(.x = A, .y = snp.list, .f = ~ .y[!.y %in% .x]) # identify minor allele based on major allele
  output <- mutate(input, major = A, minor = B)
  
  return(output)
}

first_half <- function(inp){
  mid <- round(length(inp)/2, digits = 0)
  return(inp[1:mid])
}

second_half <- function(inp){
  mid <- round(length(inp)/2, digits = 0)
  return(inp[(mid+1):length(inp)])
}

scan_populations <- function(inp, pop.name, chrm, param = param){
  # Input descriptions
  # inp = a tibble of genotype data with varieties starting on the 7th column
  # pop.name = Name of the population so that we can track the state of operations
  # chrm = chromosome labels used, e.g. 1A to 7D
  # param = parameters for the scans
  
  # WRITE
  # one thap genotype file per chromosome is written (n = 21)
  # output directory is assumed to be "rehh_files/genotype/"
  message(str_c("Writing .thap genotype files for ", pop.name))
  map(.x = chrm,
      .f = ~ filter(inp, Chrom == .x) %>% dplyr::select(7:ncol(inp)) %>% 
        write_delim(file = str_c("rehh_files/genotype/", .x, ".thap"), 
                    col_names = FALSE, delim = " "))    
  
  # READ
  # one thap file per chromosome is loaded (n = 21)
  # genotype directory assumed to be "rehh_files/genotype/"
  # map info directory assumed to be "rehh_files/map/"
  message(str_c("Reading .thap genotype and map files for ", pop.name))
  hh <- map(.x = chrm,
            .f = ~ data2haplohh(hap_file = str_c("rehh_files/genotype/", .x, ".thap"),
                                map_file = str_c("rehh_files/map/", .x, ".inp"),
                                min_maf = 0.00, # use all markers
                                chr.name = .x,
                                allele_coding = "none", # allele_coding to map to account for ancestral allele
                                haplotype.in.columns = TRUE)
  ) # close the map function
  names(hh) <- chrm
  
  # SCAN
  # one scan output per chromosome is generated (n = 21)
  message(str_c("Running rehh::scan_hh on ", pop.name))
  scan <- map_df(.x = chrm,                    
                 # map_df to concatenate genomewide scan results
                 # use map if you want to do chrom by chrom scan results
                 .f = ~ hh[[.x]] %>% 
                   scan_hh(., polarized = param$polarized.value,      # True means use ancestral allele information
                           scalegap = param$scale.gap.value,
                           maxgap = param$max.gap.value,
                           discard_integration_at_border = param$discard.integration.at.border) %>% 
                   cbind(tibble(SNPid = row.names(.)), .))
  #  names(scan) <- chrm
  
# OUTPUT ----
  output <- list(hh = hh, scan = scan)
  return(output)
}

calc_Fst <- function(inp, Total.pop, Sub.pop, save.cols = 6){
  # inp = genotype file with markers as rows, varieties as columns, and genotypes as AA/BB
  # Total.pop = vector of variety names to consider as members of the whole population
  # Sub.pop = list of named vectors of variety names for different sub.populations to calculate Fst for
  # save.cols = number of first columns from the inp to include in the output
  # names(Sub.pop) <- str_c("Fst.", names(Sub.pop)) 
  # add informative prefixe to sub-population names, which will get inherited as column names in the Fst tibble
  # Do not need the prefix if I keep Fst and PIC data separate.
  p <- rowSums(inp[Total.pop] == "AA")/ncol(inp[Total.pop])
  p.q <- p * (1 - p)
  inp.sub <- map(.x = Sub.pop, .f = ~ inp[.x]) # subset tibbles by each sub population
  p.sub <- map(.x = inp.sub, .f = ~ rowSums(.x == "AA")/ncol(.x))
  # Fst.sub <- map(.x = p.sub, .f = ~ (.x - p)^2 / p.q ) # outputs as a vector for each sub.pop. useful for testing.
  Fst.sub <- map_df(.x = p.sub, .f = ~ (.x - p)^2 / p.q ) # output in dataframe
  return(cbind(inp[1:save.cols], Fst.sub))
}

calc_PIC <- function(inp, Sub.pop, save.cols = 6){
  # inp = genotype file with markers as rows, varieties as columns, and genotypes as AA/BB
  # Sub.pop = list of names vectors of variety names for different sub.populations to calcuate PIC for
  # save.cols = # save.cols = number of first columns from the inp to include in the output
  # names(Sub.pop) <- str_c("PIC.", names(Sub.pop)) 
  # add informative prefix to sub-population names, which will get inherited as column names for the pic tibble.
  # Do not need the prefix if I keep Fst and PIC data separate.
  g.pop <- map(.x = Sub.pop, .f = ~ inp[.x])
  # pic <- map(.x = g.pop, .f = ~ 1 - ((rowSums(.x == "AA")/ncol(.x))^2 + (rowSums(.x == "BB")/ncol(.x))^2) )
  # outputs as a vector for each sub.pop. Useful for testing
  pic.sub <- map_df(.x = g.pop, 
                    .f = ~ 1 - ((rowSums(.x == "AA")/ncol(.x))^2 + (rowSums(.x == "BB")/ncol(.x))^2) )
  return(cbind(inp[1:save.cols], pic.sub))
}

conv_names_to_labels <- function(inp = inp, sep = "\\."){
  # inp is a list of data frames with names
  # sep specifies the character to split the names by. Default is "\\."
  # By default the labels are assigned to two columns: Habit and Pops
  rows <- map_dbl(.x = inp, .f = ~nrow(.x))
  named.rows <- map(.x = 1:length(rows),
                    .f = ~ rep(names(rows)[.x], times = rows[.x]) )
  df <- map2_df(.x = named.rows, .y = inp,
                .f = ~ mutate(.y, Category = .x))
  df <- separate(df, col = Category, into = c("Habit", "Pops"), sep = sep, remove = FALSE, 
                 extra = "merge", fill = "right")
  return(df)  
}
