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
        write_delim(path = str_c("rehh_files/genotype/", .x, ".thap"), 
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

calc_fst <- function(pop1, pop2){
  # calculates Fst using pegas::fst implementation of Wc 1984.
  # pop1 and pop2 are used to select the vectors of variety names in res$pop, which are selected out of the geno$diploid dataframe to calculate pairwise Fst.
  pop1.df <- select(geno$diploid, res[["pop"]][[pop1]]) %>% t() %>% as.data.frame(.) %>% mutate(., population = "pop1")
  pop2.df <- select(geno$diploid, res[["pop"]][[pop2]]) %>% t() %>% as.data.frame(.) %>% mutate(., population = "pop2")
  input <- rbind(pop1.df, pop2.df)
  input <- as.data.frame(unclass(input), stringsAsFactors = TRUE)
  fst <- Fst(as.loci(input))
  out <- tibble(CHR = geno$nuc$Chrom, POSITION = geno$nuc$pos, FST = fst[,2])
  return(out)
}

calc_allele_freq <- function(pop.name){
  inp <- select(geno$nuc, res[["pop"]][[pop.name]])
  out <- select(geno$nuc, SNPid, Chrom, pos) %>% 
    mutate(N = ncol(inp), # better than length(lines) in cases where not all the lines are genotyped
           frq.A = rowSums(inp == "A")/N,
           frq.C = rowSums(inp == "C")/N,
           frq.G = rowSums(inp == "G")/N,
           frq.T = rowSums(inp == "T")/N)
  return(out)
}

consolidate_results <- function(pop1, pop2){
  matrix.column.names <- c(str_c("N_", pop1), str_c("N_", pop2),
                           str_c("A_", pop1), str_c("C_", pop1), str_c("G_", pop1), str_c("T_", pop1),
                           str_c("A_", pop2), str_c("C_", pop2), str_c("G_", pop2), str_c("T_", pop2),
                           "FST", "RSB", "XPEHH" )
  out <- matrix(data = NA, nrow = nrow(geno$nuc), ncol = length(matrix.column.names))
  out[,1] <- allele.freq[[pop1]][["N"]]
  out[,2] <- allele.freq[[pop2]][["N"]]
  out[,3] <- allele.freq[[pop1]][["frq.A"]]
  out[,4] <- allele.freq[[pop1]][["frq.C"]]
  out[,5] <- allele.freq[[pop1]][["frq.G"]]
  out[,6] <- allele.freq[[pop1]][["frq.T"]]
  out[,7] <- allele.freq[[pop2]][["frq.A"]]
  out[,8] <- allele.freq[[pop2]][["frq.C"]]
  out[,9] <- allele.freq[[pop2]][["frq.G"]]
  out[,10] <- allele.freq[[pop2]][["frq.T"]]
  out[,11] <- fst.id[[pop1]][["FST"]]
  out[,12] <- rsb.id[[pop1]][["RSB"]]
  out[,13] <- xpe.id[[pop1]][["XPEHH"]]
  colnames(out) <- matrix.column.names
  out <- as_tibble(out)
  # add SNPid, cM, Chrm, Position
  out <- cbind(tibble(SNPid = geno$nuc$SNPid,
                      cM    = geno$nuc$cM,
                      CHR   = geno$nuc$Chrom,
                      POSITION = geno$nuc$pos),
               out)
  return(out)
}

extract_most_extreme_fst_snp <- function(pop1 = NULL, cr.list = cr.fst, scan.list = fst.id){
  # expects fst, rsb, and xpe results in lists named fst.id, rsb.id, and xpe.id, which also have SNPids
  num.cr <- nrow(cr.list[[pop1]])
  out <- matrix(nrow = num.cr, ncol = 3)
  for(i in 1:num.cr){
    extremes <- filter(scan.list[[pop1]], 
                       CHR == cr.list[[pop1]][["CHR"]][[i]],
                       POSITION >= cr.list[[pop1]][["START"]][[i]],
                       POSITION <= cr.list[[pop1]][["END"]][[i]] ) %>% arrange(., desc(FST)) %>% .[1,]
    out[i, ] <- c(extremes$CHR, extremes$POSITION, extremes$SNPid)
  }
  out.df <- data.frame(CHR = out[,1], POSITION = as.numeric(out[,2]))
  row.names(out.df) <- out[,3]
  return(out.df)
}

extract_most_extreme_rsb_snp <- function(pop1 = NULL, cr.list = cr.rsb, scan.list = rsb.id){
  # expects fst, rsb, and xpe results in lists named fst.id, rsb.id, and xpe.id, which also have SNPids
  num.cr <- nrow(cr.list[[pop1]])
  out <- matrix(nrow = num.cr, ncol = 3)
  for(i in 1:num.cr){
    extremes <- filter(scan.list[[pop1]], 
                       CHR == cr.list[[pop1]][["CHR"]][[i]],
                       POSITION >= cr.list[[pop1]][["START"]][[i]],
                       POSITION <= cr.list[[pop1]][["END"]][[i]] ) %>% arrange(., desc(RSB)) %>% .[1,]
    out[i, ] <- c(extremes$CHR, extremes$POSITION, extremes$SNPid)
  }
  out.df <- data.frame(CHR = out[,1], POSITION = as.numeric(out[,2]))
  row.names(out.df) <- out[,3]
  return(out.df)
}

extract_most_extreme_xpehh_snp <- function(pop1 = NULL, cr.list = cr.xpe, scan.list = xpe.id){
  num.cr <- nrow(cr.list[[pop1]])
  out <- matrix(nrow = num.cr, ncol = 3)
  for(i in 1:num.cr){
    extremes <- filter(scan.list[[pop1]], 
                       CHR == cr.list[[pop1]][["CHR"]][[i]],
                       POSITION >= cr.list[[pop1]][["START"]][[i]],
                       POSITION <= cr.list[[pop1]][["END"]][[i]] ) %>% arrange(., desc(XPEHH)) %>% .[1,]
    out[i, ] <- c(extremes$CHR, extremes$POSITION, extremes$SNPid)
  }
  out.df <- data.frame(CHR = out[,1], POSITION = as.numeric(out[,2]))
  row.names(out.df) <- out[,3]
  return(out.df)
}

custom_manhattan <- function(pop1 = NULL){
  # expects that fst, rsb, and xpe results are in lists named fst, rsb, and xpe.
  # expects that candidate regions are in cr.fst, cr.rsb, and cr.xpe.
  # expects the fst thresholds are in quantile.threshold named vector
  # pop1 name provide here is used to extract the right data from the above lists.
  par(mfrow = c(3,1))
  manhattanplot(data = fst.id[[pop1]], cr = cr.fst[[pop1]], threshold = quantile.threshold[[pop1]], mrk = fst.mrk[[pop1]])
  mtext(pop1) # add a label above the first plot
  manhattanplot(data = rsb.id[[pop1]], cr = cr.rsb[[pop1]], mrk = rsb.mrk[[pop1]])
  manhattanplot(data = xpe.id[[pop1]], cr = cr.xpe[[pop1]], mrk = xpe.mrk[[pop1]])
  #out <- plot_grid(p.fst, p.rsb, p.xpe, nrow = 3, ncol = 1)
  #return(out)
}

create_map_file_by_category <- function(cr, kim = kim, chr.beg = chr.beg, chr.end = chr.end){
  # cr is the input candidate region data frame
  # Font size based on no. of extremal markers
  #cr <- mutate(cr, Font.size = if_else(N_EXTR_MRK < 17, "S14",
  #                                     if_else(N_EXTR_MRK >= 17 & N_EXTR_MRK < 39, "S18", # from 75-93rd percentile
  #                                     if_else(N_EXTR_MRK >= 39 & N_EXTR_MRK < 56, "S22", # from 93-98 percentile
  #                                             if_else(N_EXTR_MRK >= 56, "S26", "")))))
  
  # Font size based on KB size
  q.75 <- quantile(cr$KB, 0.75)
  q.9375 <- quantile(cr$KB, 0.9375)
  q.99 <- quantile(cr$KB, 0.99)
  map <- mutate(cr, Font.size = if_else(KB >= q.75 & KB < q.9375, "S12", # from 75-93rd percentile
                                        if_else(KB >= q.9375 & KB < q.99, "S16", # from 93-99 percentile
                                                if_else(KB >= q.99, "S20", ""))))
  map <- tibble(Chrom = map$CHR, Tag = map$label, # Use Pop.H if plotting both habit together
                Pos = map$START, Color = NA, 
                Habit = map$Habit, Stat = map$Stat, Size = map$Font.size)
  # Update color based on stat
  map <- mutate(map, Color = if_else(Stat == "Fst", "C2",
                                     if_else(Stat == "Rsb", "C3", 
                                             if_else(Stat == "xpEHH", "C4", ""))))
  map <- rbind(map, kim, chr.beg, chr.end) %>% arrange(., Chrom, Pos)
  map <- mutate(map, Genome = substr(Chrom, 2, 2))
  return(map)
}

write_mapchart_files <- function(inp, category.name = NULL){
  grp <- list(one   = c("1A", "1B", "1D"), two  = c("2A", "2B", "2D"), three = c("3A", "3B", "3D"),
              four  = c("4A", "4B", "4D"), five = c("5A", "5B", "5D"), six   = c("6A", "6B", "6D"),
              seven = c("7A", "7B", "7D"))
  for(i in 1:7){
    # loops through the 7 chromosome groups to generate the output files. I am using loops because it is only 7 iterations.
    df <- map(.x = grp[[i]], .f = ~filter(inp, Chrom == .x))
    df.labels <- map(.x = grp[[i]], 
                     .f = ~ tibble(Chrom = "", 
                                   # columns exported in the mapchart input
                                   Tag = "group", Pos = str_c(.x), Size = "", Color = "",
                                   Habit = "", Stat = "", Genome = ""))
    # combine the filtered tibbles based on matching column names.
    out <- map2_df(.x = df.labels, .y = df, .f = ~ rbind(.x, .y))
    out[2:5] %>% write_delim(str_c("output/mapchart/", category.name, "_grp_", i, ".mct"), delim = " ", col_names = FALSE)
    rm(df, df.labels, out)
  }
}
