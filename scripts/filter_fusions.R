
################################################################################
####### -------                   Goals                          ------- #######

# I am using this script to employ additional filtering to extract potential fusions

# filter_vcfs.sh will produce many files and some of those may still be of 
# interest for confirmation


################################################################################
####### -------               R packages and paths               ------- #######

## --- setwd()
setwd("/Volumes/Partition 2/Core/hayashi_twist_Dec2025/")

## --- Read in R packages
library(plyr)
library(dplyr)
library(openxlsx)
library(tidyverse)
library(rtracklayer)

## --- Paths to input files
shared_files <- "/Volumes/Partition 2/Core/hayashi_sarcoma_fusion_files/Data/"
path_fusions <- paste0(shared_files, "sarcoma_genes/canonical_fusions.xlsx") # fusions with CANONICAL gene order
path_gene_info <- paste0(shared_files, "sarcoma_genes/gene_info.xlsx") # Info for each gene partner
path_fusions_special <- paste0(shared_files, "sarcoma_genes/special_fusions.xlsx") # Intrachromosomal and outside of gene body
path_twist_canonical_region_bed <- paste0(shared_files, "targets_CU-Hayashi_PedSolidTumorTransloc_TE-92680772_hg38.bed")
# Note that `path_twist_canonical_region_bed` contains the desired canonical regions
# we wished were covered - not regions specifically having probes covered

# This file will be used to filter canonical translocations:
path_unfiltered_no_header_vcf <- "Analysis/vcfs_loose_FUSIONS/merged_vcfs/canonical_fusions_no_header.vcf"
# Note that pass_BND would often have the same effect but I want to make sure
# I am not missing anything, even poorly covered (but canonical) translocations

## --- IGV parameters
igv_padding <- 250 # bases to add and subtract to each breakpoint
igv_xml_template <- "Data/igv_session_template.xml"

## --- primer3 parameters
# v1 primer pairs will work if a deletion from a breakpoint is <50 bp
primer_fragment_from_breakpoint_v1 <- 250 # Meaning the full fragment size **could** be 500bp
primer_internal_to_include_v1 <- 50 # Meaning in the chimeric sequence, bases 200-300 need to be included but not part of the primer pair
# v2 primer pairs will work if a deletion is up to 250 bp
primer_fragment_from_breakpoint_v2 <- 500 # Meaning the full fragment size **could** be 1000bp
primer_internal_to_include_v2 <- 250 # Meaning in the chimeric sequence, bases 250-750 need to be included but not part of the primer pair
# Files
path_fasta <- paste0(shared_files, "WholeGenomeFasta/Homo_sapiens_assembly38.fasta")
primer3_template <- paste0(shared_files, "primer3_template")
# Clean up files because shared directory has a space
path_fasta <- gsub(" ", "\\\\ ", path_fasta)
primer3_template <- gsub(" ", "\\\\ ", primer3_template)
output_primer3 <- "Analysis/primer3/"
dir.create(output_primer3)

## --- Paths to output files and directories 
output_directory <- "Analysis/main_filtered_output/" 
dir.create(output_directory)


################################################################################
####### -------        Extract canonical fusions from VCF        ------- #######


## --- Filter `path_unfiltered_no_header_vcf` by canonical translocations

fusion_partners <- read.xlsx(path_fusions)

# Note that alterations failing TIDDIT FILTERS will be here as well but that seems like
# an advantage because any evidence for a canonical translocation is worthy of 
# validating with Sanger sequencing

# I found that SNPEFF annotates the EWSR1-ERG fusion as 'ERG&EWSR1' so I know 
# for sure that I need to look for both canonical and noncanonical orders 

# Duplicate fusion_partners by using fusions in both gene orders
fusion_partners_reversed <- fusion_partners
fusion_partners_reversed$fusion <- paste0(fusion_partners_reversed$gene2,
                                          "&",
                                          fusion_partners_reversed$gene1)
fusion_partners <- rbind(fusion_partners, fusion_partners_reversed)

VCF <- read.table(file = path_unfiltered_no_header_vcf,
                  sep = "\t", header = T, stringsAsFactors = F,
                  comment.char = "", check.names = F)
saveRDS(VCF, paste0(output_directory, "rds_unfiltered_vcf_no_header_loose_v3.rds")) # Increases speed if needs to be run again
VCF <- readRDS(paste0(output_directory, "rds_unfiltered_vcf_no_header_loose_v3.rds"))

n_fusions <- c()
fusion_names <- c()
# Loop through each fusion and if found, record the row # and fusion
for(fusion in fusion_partners$fusion){
  
  print(paste0("searching for: ", fusion))
  
  n_temp <- grep(fusion, VCF$INFO)
  
  if(length(n_temp) > 0){
    n_fusions <- c(n_fusions,
                   n_temp)
    fusion_names <- c(fusion_names,
                      rep(fusion, length(n_temp)))
  }
  
}

n_fusions
length(n_fusions)

# # of n_fusions == 56

# Subset VCF to the alterations with the fusion
VCF_fusions <- cbind.data.frame(fusion = "",
                                VCF)
VCF_fusions$fusion[n_fusions] <- fusion_names
VCF_fusions <- VCF_fusions[n_fusions,]

# Bind gene1 and gene2 and add gene info
gene_info <- read.xlsx(path_gene_info)
gene1_info <- gene_info
colnames(gene1_info) <- c("gene1", paste0("gene1_", colnames(gene1_info)[2:ncol(gene1_info)]))
gene2_info <- gene_info
colnames(gene2_info) <- c("gene2", paste0("gene2_", colnames(gene2_info)[2:ncol(gene2_info)]))
fusion_partners <- join(fusion_partners, gene1_info)
fusion_partners <- join(fusion_partners, gene2_info)
VCF_fusions <- join(fusion_partners,
                    VCF_fusions)
VCF_fusions <- VCF_fusions[!(is.na(VCF_fusions$POS)),] # remove fusions not present
# VCF_fusions should now be the same length as length(n_fusions)
length(n_fusions) == nrow(VCF_fusions)

saveRDS(VCF_fusions, paste0(output_directory, "rds_filtered_vcf_loose.rds")) # Just saving for future speed
VCF_fusions <- readRDS(paste0(output_directory, "rds_filtered_vcf_loose.rds"))


## --- Need to extract the FORMAT to the sample level
# Hold sample-level information within `VCF_fusions_sample_level`
VCF_fusions_sample_level <- VCF_fusions
VCF_fusions_sample_level <- cbind.data.frame(Sample = NA,
                                             Sample_FORMAT=NA,
                                             VCF_fusions_sample_level)
VCF_fusions_sample_level <- VCF_fusions_sample_level[1,]
VCF_fusions_sample_level[1,] <- NA
# Define where the samples are within `VCF_fusions`
sample_start <- which(colnames(VCF_fusions) == "FORMAT") + 1
#sample_end <- which(colnames(VCF_fusions) == "SVTYPE") - 1
sample_end <- ncol(VCF_fusions)
samples_sub <- VCF_fusions[,sample_start:sample_end]
# Loop through each alteration and sample-level info to `VCF_fusions_sample_level`
for(i in 1:nrow(VCF_fusions)){
  variant_vector <- as.character(as.vector(samples_sub[i,]))
  n_sample <- which(!(variant_vector == "./.:.:.:.:.:.:.:."))
  
  samples <- colnames(samples_sub)[n_sample]
  formats <- as.character(as.vector(samples_sub[i,n_sample]))
  temp_df <- VCF_fusions[i,]
  
  temp_df <- cbind.data.frame(Sample = samples,
                              Sample_FORMAT = formats,
                              temp_df)
  
  VCF_fusions_sample_level <- rbind(VCF_fusions_sample_level,
                                    temp_df)
}
VCF_fusions_sample_level <- VCF_fusions_sample_level[!(is.na(VCF_fusions_sample_level$Sample)),] # remove fusions without samples
#VCF_fusions_sample_level <- na.omit(VCF_fusions_sample_level)
VCF_fusions_sample_level <- VCF_fusions_sample_level[,1:(which(colnames(VCF_fusions_sample_level) == "FORMAT"))]

# Extract sample_FORMAT column and add back to the DF
# These statistics seem to have more accurate coverage than INFO
# GT:CN:COV:DV:RV:LQ:RR:DR
FORMAT_vector <- strsplit(VCF_fusions_sample_level$FORMAT[1], split = ":")[[1]]
for(metric in FORMAT_vector){
  n_metrics_in_FORMAT <- which(FORMAT_vector == metric)
  VCF_fusions_sample_level[,metric] <- sapply(strsplit(VCF_fusions_sample_level$Sample_FORMAT,
                                                       split = ":"),
                                              function(x){x[[n_metrics_in_FORMAT]]})
}
VCF_fusions_sample_level$RR_A <- sapply(strsplit(VCF_fusions_sample_level$RR,
                                                 split = ","), function(x){x[[1]]})
VCF_fusions_sample_level$RR_B <- sapply(strsplit(VCF_fusions_sample_level$RR,
                                                 split = ","), function(x){x[[2]]})
VCF_fusions_sample_level$DR_A <- sapply(strsplit(VCF_fusions_sample_level$DR,
                                                 split = ","), function(x){x[[1]]})
VCF_fusions_sample_level$DR_B <- sapply(strsplit(VCF_fusions_sample_level$DR,
                                                 split = ","), function(x){x[[2]]})


## --- Extract INFO column and add back to the DF
# I don't think the INFO-level statistics could ever have been trusted at the 
# sample-level outside of the positions, but we need those at least
# Adopted from: https://stackoverflow.com/questions/60523778/parse-vcf-files-info-to-an-r-dataframe
VCF_fusions_sample_level <- cbind.data.frame(row = 1:nrow(VCF_fusions_sample_level),
                                             VCF_fusions_sample_level)
VCF_fusions_sample_level_INFO_df <- VCF_fusions_sample_level %>% 
                        mutate(across(everything(), as.character)) %>% 
                        select(c(row, INFO)) %>%
                        separate(INFO, sep = ";", into = paste0("cols", as.character(1:8))) %>%
                        #mutate(sample = row_number()) %>% # I added row instead
                        pivot_longer(-row) %>%
                        # here you have to unselect the name column, since this is actually wrong
                        select(-name) %>% 
                        # separate the strings that contain the data
                        separate(value, into = c("group", "value"), sep = "=") %>% 
                        # put it back to wide format
                        pivot_wider( names_from = group, 
                                     values_from = value)
VCF_fusions_sample_level <- join(VCF_fusions_sample_level,
                                 VCF_fusions_sample_level_INFO_df)

# Write file containing all rows
write.table(VCF_fusions_sample_level,
            file = paste0(output_directory, "temp_all_rows_unfiltered_vcf_subset_to_canonical_fusions.txt"),
            sep = "\t", row.names = F, quote = F)
VCF_fusions_sample_level <- read.table(file = paste0(output_directory, "temp_all_rows_unfiltered_vcf_subset_to_canonical_fusions.txt"),
                          sep = "\t", header = T, stringsAsFactors = F,
                          comment.char = "", check.names = F)


## --- Subset to only a single alteration per fusion
# TIDDIT returns the break point and the reverse on two different rows but they
# essentially contain the same information
# I will keep only the row where gene1_chr == #CHROM
table(VCF_fusions_sample_level$fusion)
VCF_fusions_sample_level <- VCF_fusions_sample_level[VCF_fusions_sample_level$gene1_chr == VCF_fusions_sample_level$'#CHROM',]
table(VCF_fusions_sample_level$fusion)
# This does not cut intrachromosomal fusions in half but I think that is ok
# because they don't seem to have the same mirroring issue as the interchromosomal fusions


## --- Extract gene1 and gene2 breakpoints
VCF_fusions_sample_level$gene1_chr_breakpoint <- VCF_fusions_sample_level$POS
VCF_fusions_sample_level$gene2_chr_breakpoint <- gsub("\\]|N|\\[", "", VCF_fusions_sample_level$ALT)
VCF_fusions_sample_level$gene2_chr_breakpoint <- gsub("CHR", "chr", VCF_fusions_sample_level$gene2_chr_breakpoint)
check_gene2_chr <- sapply(strsplit(VCF_fusions_sample_level$gene2_chr_breakpoint,
                                   split = ":"),
                          function(x){x[[1]]})
check_gene2_pos <- check_gene2_chr == VCF_fusions_sample_level$gene2_chr # Should all be true
check_gene2_pos # Should be true for interchromosomal fusions
if(FALSE %in% check_gene2_pos){
  
  likely_intrachromosomal <- which(FALSE == check_gene2_pos)
  print(paste0("Intrachromosomal fusions expected in ", paste0(likely_intrachromosomal,
                                                        collapse = " ")))
  
  # Replacing `gene2_chr_breakpoint` with `#CHROM` and `END`
  VCF_fusions_sample_level$gene2_chr_breakpoint[likely_intrachromosomal] <- paste0(
    VCF_fusions_sample_level$gene2_chr[likely_intrachromosomal],
    ":",
    VCF_fusions_sample_level$END[likely_intrachromosomal])
  
}

VCF_fusions_sample_level$gene2_chr_breakpoint <- as.numeric(sapply(strsplit(VCF_fusions_sample_level$gene2_chr_breakpoint,
                                                                 split = ":"),
                                                        function(x){x[[2]]}))


################################################################################
##### ----- For IGV, generate gene1 and gene2 ranges and a BED file ----- ######

VCF_fusions_sample_level$gene1_igv_range <- paste0(VCF_fusions_sample_level$gene1_chr,
                                                   ":", VCF_fusions_sample_level$gene1_chr_breakpoint-igv_padding,
                                                   "-", VCF_fusions_sample_level$gene1_chr_breakpoint+igv_padding)
VCF_fusions_sample_level$gene2_igv_range <- paste0(VCF_fusions_sample_level$gene2_chr,
                                                   ":", VCF_fusions_sample_level$gene2_chr_breakpoint-igv_padding,
                                                   "-", VCF_fusions_sample_level$gene2_chr_breakpoint+igv_padding)
bed_file_igv <- cbind.data.frame(chr = VCF_fusions_sample_level$gene1_chr,
                                 chromStart = VCF_fusions_sample_level$gene1_chr_breakpoint-igv_padding,
                                 chromEnd = VCF_fusions_sample_level$gene1_chr_breakpoint+igv_padding,
                                 name = paste0(VCF_fusions_sample_level$Sample, ".",
                                               VCF_fusions_sample_level$fusion, ".gene1.",
                                               VCF_fusions_sample_level$row)) 
# row is added to make sure they are unique even if the same sample has two fusions of the same pairing
bed_file_igv <- rbind(bed_file_igv,
                      cbind.data.frame(chr = VCF_fusions_sample_level$gene2_chr,
                                       chromStart = VCF_fusions_sample_level$gene2_chr_breakpoint-igv_padding,
                                       chromEnd = VCF_fusions_sample_level$gene2_chr_breakpoint+igv_padding,
                                       name = paste0(VCF_fusions_sample_level$Sample, ".",
                                                     VCF_fusions_sample_level$fusion, ".gene2.",
                                                     VCF_fusions_sample_level$row)))
# Write bed file
write.table(bed_file_igv,
            file = paste0(output_directory, "igv_ranges.bed"),
            sep = "\t", col.names = F, row.names = F, quote = F)

## --- For primer3, generate gene1 and gene2 PRIMER ranges 

# To get the chimeric sequence:
#   If positive strand, subtract from breakpoint and extract the forward strand
#   If negative strand, add from the breakpoint and extract the reverse complement
# Also indicate the excluded middle region of the primer ranges and the full
# primer ranges as a bed file

bed_primer_regions <- cbind.data.frame(chr = NA,
                                       chromStart = NA,
                                       chromEnd = NA,
                                       name = NA) 
for(i in 1:nrow(VCF_fusions_sample_level)){
  
  print(i)
  gene1 <- VCF_fusions_sample_level$gene1[i]
  gene2 <- VCF_fusions_sample_level$gene2[i]
  sample <- VCF_fusions_sample_level$Sample[i]
  
  gene1_excluded_region_v1 <- paste0(VCF_fusions_sample_level$gene1_chr[i], ":",
                               (VCF_fusions_sample_level$gene1_chr_breakpoint[i]-primer_internal_to_include_v1),
                               "-", (VCF_fusions_sample_level$gene1_chr_breakpoint[i]+primer_internal_to_include_v1))
  gene1_excluded_region_v2 <- paste0(VCF_fusions_sample_level$gene1_chr[i], ":",
                               (VCF_fusions_sample_level$gene1_chr_breakpoint[i]-primer_internal_to_include_v2),
                               "-", (VCF_fusions_sample_level$gene1_chr_breakpoint[i]+primer_internal_to_include_v2))
  gene2_excluded_region_v1 <- paste0(VCF_fusions_sample_level$gene2_chr[i], ":",
                                     (VCF_fusions_sample_level$gene2_chr_breakpoint[i]-primer_internal_to_include_v1),
                                     "-", (VCF_fusions_sample_level$gene2_chr_breakpoint[i]+primer_internal_to_include_v1))
  gene2_excluded_region_v2 <- paste0(VCF_fusions_sample_level$gene2_chr[i], ":",
                                     (VCF_fusions_sample_level$gene2_chr_breakpoint[i]-primer_internal_to_include_v2),
                                     "-", (VCF_fusions_sample_level$gene2_chr_breakpoint[i]+primer_internal_to_include_v2))
  
  # Define gene1 primer regions
  if(VCF_fusions_sample_level$gene1_strand[i] == "positive"){
    gene1_v1 <- paste0(VCF_fusions_sample_level$gene1_chr[i], ":",
                       (VCF_fusions_sample_level$gene1_chr_breakpoint[i] - primer_fragment_from_breakpoint_v1),
                       "-", VCF_fusions_sample_level$gene1_chr_breakpoint[i])
    gene1_v2 <- paste0(VCF_fusions_sample_level$gene1_chr[i], ":",
                       (VCF_fusions_sample_level$gene1_chr_breakpoint[i] - primer_fragment_from_breakpoint_v2),
                       "-", VCF_fusions_sample_level$gene1_chr_breakpoint[i])
    
    seq_output <- system(paste0("samtools faidx ", path_fasta, " ",
                                gene1_v1), intern = TRUE) # -i means rev complement
    VCF_fusions_sample_level$gene1_v1_primer_region[i] <- paste0(seq_output[2:length(seq_output)], collapse = "")
    seq_output <- system(paste0("samtools faidx ", path_fasta, " ",
                                gene1_v2), intern = TRUE) # -i means rev complement
    VCF_fusions_sample_level$gene1_v2_primer_region[i] <- paste0(seq_output[2:length(seq_output)], collapse = "")
    
  }else{
    gene1_v1 <- paste0(VCF_fusions_sample_level$gene1_chr[i], ":",
                       VCF_fusions_sample_level$gene1_chr_breakpoint[i], "-",
                       (VCF_fusions_sample_level$gene1_chr_breakpoint[i] + primer_fragment_from_breakpoint_v1))
    gene1_v2 <- paste0(VCF_fusions_sample_level$gene1_chr[i], ":",
                       VCF_fusions_sample_level$gene1_chr_breakpoint[i], "-",
                       (VCF_fusions_sample_level$gene1_chr_breakpoint[i] + primer_fragment_from_breakpoint_v2))
    
    seq_output <- system(paste0("samtools faidx ", path_fasta, " ",
                                gene1_v1, " -i"), intern = TRUE) # -i means rev complement
    VCF_fusions_sample_level$gene1_v1_primer_region[i] <- paste0(seq_output[2:length(seq_output)], collapse = "")
    seq_output <- system(paste0("samtools faidx ", path_fasta, " ",
                                gene1_v2, " -i"), intern = TRUE) # -i means rev complement
    VCF_fusions_sample_level$gene1_v2_primer_region[i] <- paste0(seq_output[2:length(seq_output)], collapse = "")
    
  }
  
  # Define gene2 primer regions (flip sign compared to gene1)
  # Also need to flip the 'start'/'end'
  # P2 always reads from the breakpoint. P1 always reads to it
  if(VCF_fusions_sample_level$gene2_strand[i] == "positive"){
    gene2_v1 <- paste0(VCF_fusions_sample_level$gene2_chr[i], ":",
                       VCF_fusions_sample_level$gene2_chr_breakpoint[i], "-",
                       (VCF_fusions_sample_level$gene2_chr_breakpoint[i] + primer_fragment_from_breakpoint_v1))
    gene2_v2 <- paste0(VCF_fusions_sample_level$gene2_chr[i], ":",
                       VCF_fusions_sample_level$gene2_chr_breakpoint[i], "-",
                       (VCF_fusions_sample_level$gene2_chr_breakpoint[i] + primer_fragment_from_breakpoint_v2))
    
    seq_output <- system(paste0("samtools faidx ", path_fasta, " ",
                                gene2_v1), intern = TRUE) # -i means rev complement
    VCF_fusions_sample_level$gene2_v1_primer_region[i] <- paste0(seq_output[2:length(seq_output)], collapse = "")
    seq_output <- system(paste0("samtools faidx ", path_fasta, " ",
                                gene2_v2), intern = TRUE) # -i means rev complement
    VCF_fusions_sample_level$gene2_v2_primer_region[i] <- paste0(seq_output[2:length(seq_output)], collapse = "")
    
  }else{
    gene2_v1 <- paste0(VCF_fusions_sample_level$gene2_chr[i], ":",
                       (VCF_fusions_sample_level$gene2_chr_breakpoint[i] - primer_fragment_from_breakpoint_v1),
                       "-",
                       VCF_fusions_sample_level$gene2_chr_breakpoint[i])
    gene2_v2 <- paste0(VCF_fusions_sample_level$gene2_chr[i], ":",
                       (VCF_fusions_sample_level$gene2_chr_breakpoint[i] - primer_fragment_from_breakpoint_v2),
                       "-",
                       VCF_fusions_sample_level$gene2_chr_breakpoint[i])
    
    seq_output <- system(paste0("samtools faidx ", path_fasta, " ",
                                gene2_v1, " -i"), intern = TRUE) # -i means rev complement
    VCF_fusions_sample_level$gene2_v1_primer_region[i] <- paste0(seq_output[2:length(seq_output)], collapse = "")
    seq_output <- system(paste0("samtools faidx ", path_fasta, " ",
                                gene2_v2, " -i"), intern = TRUE) # -i means rev complement
    VCF_fusions_sample_level$gene2_v2_primer_region[i] <- paste0(seq_output[2:length(seq_output)], collapse = "")
    
  }
  
  regions <- c(gene1_v1, gene1_v2, gene2_v1, gene2_v2,
               gene1_excluded_region_v1, gene1_excluded_region_v2, 
               gene2_excluded_region_v1, gene2_excluded_region_v2)
  names(regions) <- c(paste0(sample, "_", gene1, "_v1_primer_region"),
                      paste0(sample, "_", gene1, "_v2_primer_region"),
                      paste0(sample, "_", gene2, "_v1_primer_region"),
                      paste0(sample, "_", gene2, "_v2_primer_region"),
                      paste0(sample, "_", gene1, "_v1_excluded_primer_region"),
                      paste0(sample, "_", gene1, "_v2_excluded_primer_region"),
                      paste0(sample, "_", gene2, "_v1_excluded_primer_region"),
                      paste0(sample, "_", gene2, "_v2_excluded_primer_region"))
  for(ii in 1:length(regions)){
    region <- regions[ii]
    split_region <- strsplit(region, split = ":|-")[[1]]
    bed_primer_regions <- rbind(bed_primer_regions,
                                cbind.data.frame(chr = split_region[1],
                                                 chromStart = as.numeric(split_region[2]),
                                                 chromEnd = as.numeric(split_region[3]),
                                                 name = names(regions)[ii]))
  }
   
}

VCF_fusions_sample_level$chimeric_v1 <- paste0(VCF_fusions_sample_level$gene1_v1_primer_region,
                                               VCF_fusions_sample_level$gene2_v1_primer_region)
VCF_fusions_sample_level$chimeric_v2 <- paste0(VCF_fusions_sample_level$gene1_v2_primer_region,
                                               VCF_fusions_sample_level$gene2_v2_primer_region)

bed_primer_regions <- na.omit(bed_primer_regions)
# Write bed file
write.table(bed_primer_regions,
            file = paste0(output_directory, "primer_ranges.bed"),
            sep = "\t", col.names = F, row.names = F, quote = F)


################################################################################
####### -------       Generate primer pairs using primer3        ------- #######

# Using primer3_template, replace example, chimeric_sequence, and replace_internal_exclusion

for(i in 1:nrow(VCF_fusions_sample_level)){
  print(i) # If some fusions fail, params can be relaxed in `primer3_template`
  
  # V1 primer set
  primer_set_name <- paste0(VCF_fusions_sample_level$Sample[i], ".",
                            VCF_fusions_sample_level$fusion[i], ".v1.",
                            VCF_fusions_sample_level$row[i])
  dir.create("Data")
  file.copy(gsub("\\\\ ", " ", primer3_template), "Data/temp_primer3", overwrite = T, copy.mode = F)
  system(paste0("sed -i '' 's/chimeric_sequence/", VCF_fusions_sample_level$chimeric_v1[i], "/g' Data/temp_primer3"))
  system(paste0("sed -i '' 's/example/", primer_set_name, "/g' Data/temp_primer3"))
  system(paste0("sed -i '' 's/replace_internal_exclusion/", primer_fragment_from_breakpoint_v1-primer_internal_to_include_v1,
                ",", primer_internal_to_include_v1*2,
                "/g' Data/temp_primer3"))
  primer3_output <- system("primer3_core < Data/temp_primer3", intern = TRUE)
  n_p1_seq <- grep("PRIMER_LEFT_0_SEQUENCE=", primer3_output)
  n_p2_seq <- grep("PRIMER_RIGHT_0_SEQUENCE=", primer3_output)
  n_p1_TM <- grep("PRIMER_LEFT_0_TM=", primer3_output)
  n_p2_TM <- grep("PRIMER_RIGHT_0_TM=", primer3_output)
  n_p1_location <- grep("PRIMER_LEFT_0=", primer3_output)
  n_p2_location <- grep("PRIMER_RIGHT_0=", primer3_output)
  # Write to table
  VCF_fusions_sample_level$P1_sequence_v1[i] <- strsplit(primer3_output[n_p1_seq], split = "=")[[1]][2]
  VCF_fusions_sample_level$P2_sequence_v1[i] <- strsplit(primer3_output[n_p2_seq], split = "=")[[1]][2]
  VCF_fusions_sample_level$P1_TM_v1[i] <- strsplit(primer3_output[n_p1_TM], split = "=")[[1]][2]
  VCF_fusions_sample_level$P2_TM_v1[i] <- strsplit(primer3_output[n_p2_TM], split = "=")[[1]][2]
  # Parse location of primers
  # P1
  if(VCF_fusions_sample_level$gene1_strand[i] == "positive"){
    p1_start <- VCF_fusions_sample_level$gene1_chr_breakpoint[i] - primer_fragment_from_breakpoint_v1
    p1_start <- p1_start + as.numeric(strsplit(primer3_output[n_p1_location], split = "=|,")[[1]][2])
    p1_start <- p1_start - 1 # Included to fit IGV
    p1_end <- p1_start + as.numeric(strsplit(primer3_output[n_p1_location], split = "=|,")[[1]][3])
  }else{
    p1_start <- VCF_fusions_sample_level$gene1_chr_breakpoint[i] + primer_fragment_from_breakpoint_v1
    p1_start <- p1_start - as.numeric(strsplit(primer3_output[n_p1_location], split = "=|,")[[1]][2])
    p1_end <- p1_start - as.numeric(strsplit(primer3_output[n_p1_location], split = "=|,")[[1]][3])
  }
  # P2 - more complicated because the primer3 positions mean that the p1 (gene1) is already included
  # Need to subtract the p1 half from the distance reported in primer3
  if(VCF_fusions_sample_level$gene2_strand[i] == "positive"){
    p2_start <- VCF_fusions_sample_level$gene2_chr_breakpoint[i] - primer_fragment_from_breakpoint_v1
    p2_start <- p2_start + as.numeric(strsplit(primer3_output[n_p2_location], split = "=|,")[[1]][2])
    p2_start <- p2_start - 1 # Included to fit IGV
    p2_end <- p2_start - as.numeric(strsplit(primer3_output[n_p2_location], split = "=|,")[[1]][3])
  }else{
    p2_start <- VCF_fusions_sample_level$gene2_chr_breakpoint[i] + primer_fragment_from_breakpoint_v1
    p2_start <- p2_start - as.numeric(strsplit(primer3_output[n_p2_location], split = "=|,")[[1]][2])
    #p2_start <- p2_start + 1 # Included to fit IGV
    p2_end <- p2_start + as.numeric(strsplit(primer3_output[n_p2_location], split = "=|,")[[1]][3])
  }
  # Add primer locations
  VCF_fusions_sample_level$P1_location_v1[i] <- paste0(VCF_fusions_sample_level$gene1_chr[i], ":",
                                                    p1_start, "-", p1_end)
  VCF_fusions_sample_level$P2_location_v1[i] <- paste0(VCF_fusions_sample_level$gene2_chr[i], ":",
                                                       p2_start, "-", p2_end)
  VCF_fusions_sample_level$primer_set_name_v1[i] <- primer_set_name
  write.table(primer3_output,
              file = paste0(output_primer3, primer_set_name, ".txt"),
              sep = "\t", row.names = F, col.names = F, quote = F)
  
  
  # V2 primer set
  primer_set_name <- paste0(VCF_fusions_sample_level$Sample[i], ".",
                            VCF_fusions_sample_level$fusion[i], ".v2.",
                            VCF_fusions_sample_level$row[i])
  dir.create("Data")
  file.copy(gsub("\\\\ ", " ", primer3_template), "Data/temp_primer3", overwrite = T, copy.mode = F)
  system(paste0("sed -i '' 's/chimeric_sequence/", VCF_fusions_sample_level$chimeric_v2[i], "/g' Data/temp_primer3"))
  system(paste0("sed -i '' 's/example/", primer_set_name, "/g' Data/temp_primer3"))
  system(paste0("sed -i '' 's/replace_internal_exclusion/", primer_fragment_from_breakpoint_v2-primer_internal_to_include_v2,
                ",", primer_internal_to_include_v2*2,
                "/g' Data/temp_primer3"))
  # V2 allows for a large deletion so I can add to the expected product size
  system(paste0("sed -i '' 's/100-500/600-900/g' Data/temp_primer3"))
  primer3_output <- system("primer3_core < Data/temp_primer3", intern = TRUE)
  n_p1_seq <- grep("PRIMER_LEFT_0_SEQUENCE=", primer3_output)
  n_p2_seq <- grep("PRIMER_RIGHT_0_SEQUENCE=", primer3_output)
  n_p1_TM <- grep("PRIMER_LEFT_0_TM=", primer3_output)
  n_p2_TM <- grep("PRIMER_RIGHT_0_TM=", primer3_output)
  n_p1_location <- grep("PRIMER_LEFT_0=", primer3_output)
  n_p2_location <- grep("PRIMER_RIGHT_0=", primer3_output)
  # Write to table
  VCF_fusions_sample_level$P1_sequence_v2[i] <- strsplit(primer3_output[n_p1_seq], split = "=")[[1]][2]
  VCF_fusions_sample_level$P2_sequence_v2[i] <- strsplit(primer3_output[n_p2_seq], split = "=")[[1]][2]
  VCF_fusions_sample_level$P1_TM_v2[i] <- strsplit(primer3_output[n_p1_TM], split = "=")[[1]][2]
  VCF_fusions_sample_level$P2_TM_v2[i] <- strsplit(primer3_output[n_p2_TM], split = "=")[[1]][2]
  # Parse location of primers
  # P1
  if(VCF_fusions_sample_level$gene1_strand[i] == "positive"){
    p1_start <- VCF_fusions_sample_level$gene1_chr_breakpoint[i] - primer_fragment_from_breakpoint_v2
    p1_start <- p1_start + as.numeric(strsplit(primer3_output[n_p1_location], split = "=|,")[[1]][2])
    p1_start <- p1_start - 1 # Included to fit IGV
    p1_end <- p1_start + as.numeric(strsplit(primer3_output[n_p1_location], split = "=|,")[[1]][3])
  }else{
    p1_start <- VCF_fusions_sample_level$gene1_chr_breakpoint[i] + primer_fragment_from_breakpoint_v2
    p1_start <- p1_start - as.numeric(strsplit(primer3_output[n_p1_location], split = "=|,")[[1]][2])
    p1_end <- p1_start - as.numeric(strsplit(primer3_output[n_p1_location], split = "=|,")[[1]][3])
  }
  # P2 - more complicated because the primer3 positions mean that the p1 (gene1) is already included
  # Need to subtract the p1 half from the distance reported in primer3
  if(VCF_fusions_sample_level$gene2_strand[i] == "positive"){
    p2_start <- VCF_fusions_sample_level$gene2_chr_breakpoint[i] - primer_fragment_from_breakpoint_v2
    p2_start <- p2_start + as.numeric(strsplit(primer3_output[n_p2_location], split = "=|,")[[1]][2])
    p2_start <- p2_start - 1 # Included to fit IGV
    p2_end <- p2_start - as.numeric(strsplit(primer3_output[n_p2_location], split = "=|,")[[1]][3])
  }else{
    p2_start <- VCF_fusions_sample_level$gene2_chr_breakpoint[i] + primer_fragment_from_breakpoint_v2
    p2_start <- p2_start - as.numeric(strsplit(primer3_output[n_p2_location], split = "=|,")[[1]][2])
    #p2_start <- p2_start + 1 # Included to fit IGV
    p2_end <- p2_start + as.numeric(strsplit(primer3_output[n_p2_location], split = "=|,")[[1]][3])
  }
  # Add primer locations
  VCF_fusions_sample_level$P1_location_v2[i] <- paste0(VCF_fusions_sample_level$gene1_chr[i], ":",
                                                       p1_start, "-", p1_end)
  VCF_fusions_sample_level$P2_location_v2[i] <- paste0(VCF_fusions_sample_level$gene2_chr[i], ":",
                                                       p2_start, "-", p2_end)
  VCF_fusions_sample_level$primer_set_name_v2[i] <- primer_set_name
  write.table(primer3_output,
              file = paste0(output_primer3, primer_set_name, ".txt"),
              sep = "\t", row.names = F, col.names = F, quote = F)
  
}

# Extract primer bed file
bed_primers <- cbind.data.frame(chr = c(sapply(strsplit(VCF_fusions_sample_level$P1_location_v1, split = ":|-"),
                                               function(x){x[1]}),
                                        sapply(strsplit(VCF_fusions_sample_level$P2_location_v1, split = ":|-"),
                                               function(x){x[1]}),
                                        sapply(strsplit(VCF_fusions_sample_level$P1_location_v2, split = ":|-"),
                                               function(x){x[1]}),
                                        sapply(strsplit(VCF_fusions_sample_level$P2_location_v2, split = ":|-"),
                                               function(x){x[1]})
                                        ),
                                chromStart = c(sapply(strsplit(VCF_fusions_sample_level$P1_location_v1, split = ":|-"),
                                                      function(x){x[2]}),
                                               sapply(strsplit(VCF_fusions_sample_level$P2_location_v1, split = ":|-"),
                                                      function(x){x[2]}),
                                               sapply(strsplit(VCF_fusions_sample_level$P1_location_v2, split = ":|-"),
                                                      function(x){x[2]}),
                                               sapply(strsplit(VCF_fusions_sample_level$P2_location_v2, split = ":|-"),
                                                      function(x){x[2]})
                                ),
                                chromEnd = c(sapply(strsplit(VCF_fusions_sample_level$P1_location_v1, split = ":|-"),
                                                    function(x){x[3]}),
                                             sapply(strsplit(VCF_fusions_sample_level$P2_location_v1, split = ":|-"),
                                                    function(x){x[3]}),
                                             sapply(strsplit(VCF_fusions_sample_level$P1_location_v2, split = ":|-"),
                                                    function(x){x[3]}),
                                             sapply(strsplit(VCF_fusions_sample_level$P2_location_v2, split = ":|-"),
                                                    function(x){x[3]})
                                ),
                                name = c(paste0("p1_", VCF_fusions_sample_level$primer_set_name_v1),
                                         paste0("p2_", VCF_fusions_sample_level$primer_set_name_v1),
                                         paste0("p1_", VCF_fusions_sample_level$primer_set_name_v2),
                                         paste0("p2_", VCF_fusions_sample_level$primer_set_name_v2))
                                ) 
# Write bed file
write.table(bed_primers,
            file = paste0(output_directory, "primer_locations.bed"),
            sep = "\t", col.names = F, row.names = F, quote = F)


################################################################################
####### -------       Flag fusions based on TIDDIT metrics       ------- #######

# TIDDIT pass/fail does not work well for enriched datasets (e.g. twist) and 
# furthermore, twist enrichment seems to identify more false positive gene fusions
# than previously detected in WGS - perhaps because of the very high coverage or
# generation of chimeric regions during the PCR step following target enrichment. 

####### ------- Flags I will use
#post_TIDDIT_flag_gene1_outside_twist # Is the gene1 breakpoint outside of canonical fusion region?
#post_TIDDIT_flag_gene2_outside_twist # Is the gene2 breakpoint outside of canonical fusion region?
#post_TIDDIT_flag_low_split_read_coverage # Is the split read coverage < 1?
#post_TIDDIT_flag_low_paired_read_coverage # Is the paired read coverage < 1?
#post_TIDDIT_flag_low_variant_coverage # Is the aggregate fusion read evidence < 5?
#post_TIDDIT_flag_low_MAPQ # Is the MAPQ < 20?
#post_TIDDIT_flag_identical_fusion_in_multiple_samples # Is the identical fusion found in multiple samples?

# post_TIDDIT_flag_summary # Summary of the above

####### ------- Flag fusions when one of the breakpoints is outside of the Twist target regions

# VCF_fusions_sample_level <- openxlsx::read.xlsx("Analysis/main_filtered_output/final_summary_table.xlsx")

# Import twist BED file
bed <- import(path_twist_canonical_region_bed)

# GRanges and bed overlap for gene1
query_gene1 <- GRanges(seqnames = VCF_fusions_sample_level$gene1_chr,
                       ranges = IRanges(start = VCF_fusions_sample_level$gene1_chr_breakpoint, 
                                        end = VCF_fusions_sample_level$gene1_chr_breakpoint))
#hits_gene1 <- findOverlaps(query_gene1, bed) # Find overlap with bed file
overlap_bool_gene1 <- countOverlaps(query_gene1, bed) > 0 # Boolean: does each query overlap?

# GRanges and bed overlap for gene2
query_gene2 <- GRanges(seqnames = VCF_fusions_sample_level$gene2_chr,
                       ranges = IRanges(start = VCF_fusions_sample_level$gene2_chr_breakpoint, 
                                        end = VCF_fusions_sample_level$gene2_chr_breakpoint))
#hits_gene2 <- findOverlaps(query_gene2, bed) # Find overlap with bed file
overlap_bool_gene2 <- countOverlaps(query_gene2, bed) > 0 # Boolean: does each query overlap?

# Add flags to table
VCF_fusions_sample_level$post_TIDDIT_flag_gene1_outside_twist <- ifelse(overlap_bool_gene1, "", "gene1_outside_twist")
VCF_fusions_sample_level$post_TIDDIT_flag_gene2_outside_twist <- ifelse(overlap_bool_gene2, "", "gene2_outside_twist")

####### ------- Flag split read evidence < 1

VCF_fusions_sample_level$post_TIDDIT_flag_low_split_read_coverage <- ifelse(VCF_fusions_sample_level$RV >= 1, 
                                                                            "", "low_split_read_coverage")

####### ------- Flag paired read evidence < 1
VCF_fusions_sample_level$post_TIDDIT_flag_low_paired_read_coverage <- ifelse(VCF_fusions_sample_level$DV >= 1, 
                                                                             "", "low_paired_read_coverage")

####### ------- Flag aggregate fusion read evidence < 5
VCF_fusions_sample_level$post_TIDDIT_flag_low_variant_coverage <- ifelse(VCF_fusions_sample_level$DV + VCF_fusions_sample_level$RV >= 5, 
                                                                         "", "low_variant_coverage")

####### ------- Flag MAPQ < 20
VCF_fusions_sample_level$post_TIDDIT_flag_low_MAPQ <- ifelse(VCF_fusions_sample_level$QUAL >= 20, 
                                                             "", "low_MAPQ")

####### ------- Flag identical fusions found in multiple samples
aggregate_fusions <- paste0(VCF_fusions_sample_level$fusion,
                            VCF_fusions_sample_level$gene1_chr_breakpoint,
                            VCF_fusions_sample_level$gene2_chr_breakpoint)
duplicated_fusions <- unique(aggregate_fusions[duplicated(aggregate_fusions)])
if(length(duplicated_fusions) > 0){
  VCF_fusions_sample_level$post_TIDDIT_flag_identical_fusion_in_multiple_samples <- ifelse(aggregate_fusions %in% duplicated_fusions,
                                                                                           "identical_fusion_in_multiple_samples", "")
}

####### ------- Collapse flags into one column
VCF_fusions_sample_level$post_TIDDIT_flag_summary <- apply(VCF_fusions_sample_level[,grep("post_TIDDIT_flag_", colnames(VCF_fusions_sample_level))],
                                                           1, function(x){
                                                             paste(x[x != ""], collapse = ";")
                                                           })
VCF_fusions_sample_level$post_TIDDIT_flag_summary[VCF_fusions_sample_level$post_TIDDIT_flag_summary == ""] <- "PASS"
VCF_fusions_sample_level$post_TIDDIT_flag_summary == "PASS"
sum(VCF_fusions_sample_level$post_TIDDIT_flag_summary == "PASS")


################################################################################
####### -------                 Write final table                ------- #######

# Note that additional sheets will be added in later R scripts to this final summary table

write.table(VCF_fusions_sample_level,
            file = paste0(output_directory, "final_summary_table.txt"),
            sep = "\t", col.names = T, row.names = F, quote = F)
write.xlsx(VCF_fusions_sample_level,
           file = paste0(output_directory, "final_summary_table.xlsx"),
           overwrite = T, colNames=T, rowNames=F,
           firstActiveCol = 3, firstActiveRow = 1)

# Write a table to store information about primers after samples are inspected in IGV
# Note that overwrite is false to limit accidental overwriting

write.xlsx(cbind.data.frame(row = VCF_fusions_sample_level$row,
                            Sample = VCF_fusions_sample_level$Sample,
                            fusion = VCF_fusions_sample_level$fusion,
                            p1_v1_status = "skip", # Will delete if the primers will work fine
                            p2_v1_status = "skip",
                            p1_v2_status = "skip",
                            p2_v2_status = "skip",
                            primer_notes = "",
                            non_sample_specific_primer_notes = ""),
           file = paste0(output_directory, "primer_summary_post_inspection.xlsx"),
           overwrite = F,
           colNames=T, rowNames=F,
           firstActiveCol = 3, firstActiveRow = 1)
# Looking for deletions, low coverage (mostly concern of deletion), or consistent
# mutations that make some primers suboptimal 


################################################################################
# -- Filter possible intrachromosomal fusions and those outside of gene body - #

stop("I stopped here for this dataset because I already identified some intrachromosomal gene fusions") 

fusions_outside_gene_body <- read.xlsx(path_fusions_special)

candidate_fusions_outside_gene_body <- list()

for(i in 1:nrow(fusions_outside_gene_body)){
  
  print(i)
  
  # Filter by gene1 breakpoint
  VCF_acceptable_breakpoint1 <- VCF[VCF$'#CHROM' == fusions_outside_gene_body$gene1_chr[i] & 
                                    VCF$POS > fusions_outside_gene_body$gene1_ideal_start[i] & 
                                    VCF$POS < fusions_outside_gene_body$gene1_ideal_end[i],]
  
  # Extract INFO
  VCF_acceptable_breakpoint1 <- cbind.data.frame(row = 1:nrow(VCF_acceptable_breakpoint1),
                                                 VCF_acceptable_breakpoint1)
  VCF_INFO <- VCF_acceptable_breakpoint1 %>% mutate(across(everything(), as.character)) %>% 
    select(c(row, INFO)) %>%
    separate(INFO, sep = ";", into = paste0("cols", as.character(1:8))) %>%
    #mutate(sample = row_number()) %>% # I added row instead
    pivot_longer(-row) %>%
    # here you have to unselect the name column, since this is actually wrong
    select(-name) %>% 
    # separate the strings that contain the data
    separate(value, into = c("group", "value"), sep = "=") %>% 
    # put it back to wide format
    pivot_wider( names_from = group, 
                 values_from = value)
  VCF_acceptable_breakpoint1_w_info <- join(VCF_acceptable_breakpoint1, VCF_INFO)
  
  # Extract gene2 info
  ## --- Extract gene1 and gene2 breakpoints
  VCF_acceptable_breakpoint1_w_info$gene2_chr_breakpoint <- gsub("\\]|N|\\[", "", VCF_acceptable_breakpoint1_w_info$ALT)
  VCF_acceptable_breakpoint1_w_info$gene2_chr_breakpoint <- gsub("CHR", "chr", VCF_acceptable_breakpoint1_w_info$gene2_chr_breakpoint)
  VCF_acceptable_breakpoint1_w_info$gene2_chr <- sapply(strsplit(VCF_acceptable_breakpoint1_w_info$gene2_chr_breakpoint,
                                     split = ":"),
                            function(x){x[[1]]})
  VCF_acceptable_breakpoint1_w_info$gene2_chr_breakpoint <- as.numeric(sapply(strsplit(VCF_acceptable_breakpoint1_w_info$gene2_chr_breakpoint,
                                                                              split = ":"),
                                                                     function(x){x[[2]]}))
  
  
  
  ## Filter by gene2 breakpoint
  #print(sum(VCF_acceptable_breakpoint1_w_info$gene2_chr == fusions_outside_gene_body$gene2_chr[i]))
  #chromosomes <- unique(VCF_acceptable_breakpoint1_w_info$gene2_chr)
  #print(chromosomes[-grep("chr", chromosomes)])
  if(fusions_outside_gene_body$fusion_type[i] == "intrachromosomal"){
    dup_n <- grep("<DUP",VCF_acceptable_breakpoint1_w_info$gene2_chr)
    if(length(dup_n) > 0){
      VCF_acceptable_breakpoint1_w_info$gene2_chr_breakpoint[dup_n] <- as.numeric(VCF_acceptable_breakpoint1_w_info$END[dup_n])
    }
    VCF_acceptable <- VCF_acceptable_breakpoint1_w_info[VCF_acceptable_breakpoint1_w_info$gene2_chr %in% c("<DUP", fusions_outside_gene_body$gene2_chr[i]) & 
                                                          VCF_acceptable_breakpoint1_w_info$gene2_chr_breakpoint > fusions_outside_gene_body$gene2_ideal_start[i] & 
                                                          VCF_acceptable_breakpoint1_w_info$gene2_chr_breakpoint < fusions_outside_gene_body$gene2_ideal_end[i],]
  } else{
    VCF_acceptable <- VCF_acceptable_breakpoint1_w_info[VCF_acceptable_breakpoint1_w_info$gene2_chr == fusions_outside_gene_body$gene2_chr[i] & 
                                                          VCF_acceptable_breakpoint1_w_info$gene2_chr_breakpoint > fusions_outside_gene_body$gene2_ideal_start[i] & 
                                                          VCF_acceptable_breakpoint1_w_info$gene2_chr_breakpoint < fusions_outside_gene_body$gene2_ideal_end[i],]
  }
  
  
  if(nrow(VCF_acceptable) > 0){
    candidate_fusions_outside_gene_body[[length(candidate_fusions_outside_gene_body) + 1]] <- VCF_acceptable
    names(candidate_fusions_outside_gene_body) <- fusions_outside_gene_body$fusion[i]
  }
  
  #print(nrow(VCF_acceptable))

}

length(candidate_fusions_outside_gene_body)

# IF > 0, INSPECT THE CONTENTS FOR CANDIDATE FUSIONS

names(candidate_fusions_outside_gene_body)
View(candidate_fusions_outside_gene_body[1])




