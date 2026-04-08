
################################################################################
####### -------                   Goals                          ------- #######

# After running filter_fusions.R, make xml files to ease the opening of each
# fusion in IGV

# I will use paths relative to the xml files so I should be able to share the
# parent directory to others

################################################################################
####### -------               R packages and paths               ------- #######

## --- setwd()
main_dir <- "/Volumes/Partition 2/Core/hayashi_twist_Dec2025//"
setwd(main_dir)

## --- Read in R packages
library(plyr)
library(dplyr)

## --- Paths to input files
filter_fusion_R_output <- paste0(main_dir, "/Analysis/main_filtered_output/")
shared_files <- "/Volumes/Partition 2/Core/hayashi_sarcoma_fusion_files/Data/"
path_igv_template <- paste0(shared_files, "igv_session_template.xml")
path_bed_primer_regions <- paste0(filter_fusion_R_output, "primer_ranges.bed")
path_bed_primer_sequences <- paste0(filter_fusion_R_output, "primer_locations.bed")
path_fusion_table <- paste0(filter_fusion_R_output, "final_summary_table.txt")
bam_subsets <- paste0(main_dir, "/Analysis/igv_subset_bams/")

## --- Paths to output files and directories 
output_directory <- paste0(bam_subsets, "IGV_files/")
dir.create(output_directory)
bed_file_dir_name <- "bed_files/"
bed_output_directory <- paste0(bam_subsets, bed_file_dir_name)
dir.create(bed_output_directory)


################################################################################
####### -------                   Read in data                   ------- #######

fusion_table <- read.table(file = path_fusion_table,
                           sep = "\t", header = T, stringsAsFactors = F,
                           check.names = F, comment.char = "")

# Adding a column that matches samples - these steps might be specific to the sample format
fusion_table$Sample
# For this data, the sample name is always printed twice, separate by "_"
dedup <- vapply(strsplit(fusion_table$Sample, "_"), function(parts) {
  half <- length(parts) / 2
  paste(parts[1:half], collapse = "_")
}, character(1))
dedup
fusion_table$Sample_similar <- dedup

#fusion_table$Sample_similar <- fusion_table$Sample
#fusion_table$Sample_similar <- gsub("_dup", "dup", fusion_table$Sample_similar)
#fusion_table$Sample_similar <- sapply(strsplit(fusion_table$Sample_similar, split = "_"),
#                                      function(x){x[[1]]})
#fusion_table$Sample_similar <- gsub("dup", "_dup", fusion_table$Sample_similar)

# Using the sample list from the **bam** files because not all samples are in
# the fusion table (an alternative could be to only use the samples in the table)
samples <- list.files(bam_subsets)
samples <- samples[grep(".bam$", samples)]
samples <- gsub(".bam", "", samples)

all(fusion_table$Sample_similar %in% samples) # should be TRUE


################################################################################
####### -------         Copy non-sample specific bed files       ------- #######

name_path_bed_primer_regions <- "primer_ranges.bed"
file.copy(path_bed_primer_regions, 
          paste0(bed_output_directory, name_path_bed_primer_regions), 
          overwrite = T)
name_path_bed_primer_sequences <- "primer_locations.bed"
file.copy(path_bed_primer_sequences, 
          paste0(bed_output_directory, name_path_bed_primer_sequences), 
          overwrite = T)


################################################################################
####### -------        Generate IGV files for each sample        ------- #######

# I reworked this section from previous to allow for multiple fusions per sample,
# which I didn't need to have before
# This rework basically means that samples with multiple fusions will have 
# '_fusionname_rowX' appended to their file name. Samples with one fusion will
# just have the fusion name appended (without _rowX) and samples without any fusions
# will just be their sample name. 'rowX' corresponds to the row column in the
# fusion table. 
# It also means that the loop is messier now

# Using the sample list from the **bam** files because not all samples are in
# the fusion table

for(sample in samples){
  
  print(sample)
  
  ## If sample is in the fusion table, I will make a more specific xml file
  if(sample %in% fusion_table$Sample_similar){ 
    
    n_sample_fusion_table <- grep(sample, fusion_table$Sample_similar) 
    # In instances where multiple canonical fusions are found, there would be multiple
    # rows matching the sample
    
    for(n_sample in n_sample_fusion_table){
      
      # File names will be different based on if there are multiple fusions
      if(length(n_sample_fusion_table) == 1){
        sample_xml <- paste0(output_directory, sample, "_", 
                             fusion_table$fusion[n_sample], ".xml")
      }else{
        sample_xml <- paste0(output_directory, sample, "_", 
                             fusion_table$fusion[n_sample], "_row",
                             fusion_table$row[n_sample], ".xml")
      }
      
      # Write and edit xml file
      file.copy(path_igv_template, sample_xml, overwrite = T)
      sample_xml <- gsub(" ", "\\\\ ", sample_xml)
      # Replace holder variables
      # <!--  replace `bam_name` with bam file name (ex. 12.bam)   -->
      # <!--  replace `bam_path` with full (OR RELATIVE) bam path (ex. /Andrew/directory/12.bam)   -->
      # <!--  replace `gene1_range` with chromosome of gene1 (ex. chr1:1200901-1201001)   -->
      # <!--  replace `gene1_range` with chromosome of gene2 (ex. chr1:1200901-1201001)   -->
      # <!--  replace `path_primer_loc_bed` with full (OR RELATIVE) path   -->
      # <!--  replace `path_primer_range_bed` with full (OR RELATIVE) path   -->
      # maximum="419.0" minimum="0.0"
      # might need to be replaced
      # Escape the ampersand for the sed commands
      sample_xml <- gsub("&", "\\\\&", sample_xml)
      system(paste0("sed -i '' 's/bam_name/", paste0(sample, ".bam"), "/g' ", sample_xml))
      # These are for building an xml with the full paths
      #system(paste0("sed -i '' 's#bam_path#", paste0(bam_subsets, sample, ".bam"), "#g' ", sample_xml)) # This is the full path
      #system(paste0("sed -i '' 's#path_primer_loc_bed#", path_bed_primer_sequences, "#g' ", sample_xml)) 
      #system(paste0("sed -i '' 's#path_primer_range_bed#", path_bed_primer_regions, "#g' ", sample_xml))
      # These are xml files with relative paths
      system(paste0("sed -i '' 's#bam_path#", paste0("../", sample, ".bam"), "#g' ", sample_xml))
      system(paste0("sed -i '' 's#path_primer_loc_bed#", 
                    paste0("../", bed_file_dir_name, name_path_bed_primer_sequences), 
                    "#g' ", sample_xml)) 
      system(paste0("sed -i '' 's#path_primer_range_bed#", 
                    paste0("../", bed_file_dir_name, name_path_bed_primer_regions), 
                    "#g' ", sample_xml))
      
      # Specify fusion region
      system(paste0("sed -i '' 's/gene1_range/", fusion_table$gene1_igv_range[n_sample], "/g' ", sample_xml))
      system(paste0("sed -i '' 's/gene2_range/", fusion_table$gene2_igv_range[n_sample], "/g' ", sample_xml))
      
    }
    
    ## If samples are not in the fusion table, a generic xml file will be made  
  }else{
    
    # Write and edit xml file
    sample_xml <- paste0(output_directory, sample, ".xml")
    file.copy(path_igv_template, sample_xml, overwrite = T)
    sample_xml <- gsub(" ", "\\\\ ", sample_xml)
    # Replace holder variables
    # <!--  replace `bam_name` with bam file name (ex. 12.bam)   -->
    # <!--  replace `bam_path` with full (OR RELATIVE) bam path (ex. /Andrew/directory/12.bam)   -->
    # <!--  replace `gene1_range` with chromosome of gene1 (ex. chr1:1200901-1201001)   -->
    # <!--  replace `gene1_range` with chromosome of gene2 (ex. chr1:1200901-1201001)   -->
    # <!--  replace `path_primer_loc_bed` with full (OR RELATIVE) path   -->
    # <!--  replace `path_primer_range_bed` with full (OR RELATIVE) path   -->
    # maximum="419.0" minimum="0.0"
    # might need to be replaced
    system(paste0("sed -i '' 's/bam_name/", paste0(sample, ".bam"), "/g' ", sample_xml))
    # These are for building an xml with the full paths
    #system(paste0("sed -i '' 's#bam_path#", paste0(bam_subsets, sample, ".bam"), "#g' ", sample_xml)) # This is the full path
    #system(paste0("sed -i '' 's#path_primer_loc_bed#", path_bed_primer_sequences, "#g' ", sample_xml)) 
    #system(paste0("sed -i '' 's#path_primer_range_bed#", path_bed_primer_regions, "#g' ", sample_xml))
    # These are xml files with relative paths
    system(paste0("sed -i '' 's#bam_path#", paste0("../", sample, ".bam"), "#g' ", sample_xml))
    system(paste0("sed -i '' 's#path_primer_loc_bed#", 
                  paste0("../", bed_file_dir_name, name_path_bed_primer_sequences), 
                  "#g' ", sample_xml)) 
    system(paste0("sed -i '' 's#path_primer_range_bed#", 
                  paste0("../", bed_file_dir_name, name_path_bed_primer_regions), 
                  "#g' ", sample_xml))
    
    # No fusion detected so will default to chr1 and chr2
    system(paste0("sed -i '' 's/gene1_range/", "chr1:1-1000", "/g' ", sample_xml))
    system(paste0("sed -i '' 's/gene2_range/", "chr1:1-1000", "/g' ", sample_xml))
  }
  
}




