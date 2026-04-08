
################################################################################
####### -------                   Goals                          ------- #######

# Make a sample summary table for the analysis 
#   Table will contain sample info, QC, and current translocations
#   Changing the underlying data would impact this summary table

################################################################################
####### -------               Read in data files                 ------- #######

setwd("/Volumes/Partition 2/Core/hayashi_twist_Dec2025/")

library(openxlsx)
library(plyr)
library(dplyr)
library(xlsx)

#multiqc_table <- openxlsx::read.xlsx("Analysis/QC_reports/multiqc_table.xlsx")
multiqc_table <- read.table(file = "Analysis/QC_reports/multiqc_general_stats.txt",
                            sep = "\t", header = T, stringsAsFactors = F)
sample_table <- openxlsx::read.xlsx("samples.xlsx")
#fusion_summary <- read.xlsx("Analysis/pass_sarcoma.xlsx")
fusion_table <- openxlsx::read.xlsx("Analysis/main_filtered_output/final_summary_table.xlsx")
primer_summary <- openxlsx::read.xlsx("Analysis/main_filtered_output/primer_summary_post_inspection.xlsx")
twist_coverage <- read.table(file = "Analysis/QC_reports/twist_coverage_summary_results.txt",
                             sep = "\t", header = T, stringsAsFactors = F)

shared_files <- "/Volumes/Partition 2/Core/hayashi_sarcoma_fusion_files/Data/"
column_descriptions <- paste0(shared_files, "summary_table_column_descriptions_twist.xlsx") 

output <- "Analysis/main_filtered_output/"


################################################################################
####### -------                 Process data                     ------- #######

####### ------- sample_table

sample_table <- sample_table[,c("patient", "Study.ID", "Fusion")]
colnames(sample_table) <- c("Sample", "Study.ID", "Cohort")
sample_table <- distinct(sample_table)

####### ------- primer_summary

# The `row` column will match up with fusion_table
#   This is being somewhat as a fusion-level ID (there can be multiple fusions per sample)
# Some of the other columns are already in fusion_table

primer_summary <- primer_summary[,!(colnames(primer_summary) %in% c("Sample", "fusion"))]


####### ------- fusion_table

# Output of sarek has sample_sample. Editing to match the other samples
fusion_table$Sample <- vapply(strsplit(fusion_table$Sample, "_"), function(parts) {
  half <- length(parts) / 2
  paste(parts[1:half], collapse = "_")
}, character(1))
#fusion_table$Sample <- gsub("_dup", "dup", fusion_table$Sample)
#fusion_table$Sample <- sapply(strsplit(fusion_table$Sample, split = "_"),
#                                      function(x){x[[1]]})
#fusion_table$Sample <- gsub("dup", "_dup", fusion_table$Sample)


####### ------- fusion_summary

#fusion_summary <- fusion_summary[fusion_summary$keep_discard == "keep",]
#fusion_summary <- fusion_summary[!(is.na(fusion_summary$keep_discard)),]
#n_first_sample <- which(colnames(fusion_summary) == "FORMAT") + 1
#fusion <- cbind.data.frame(Sample = NA,
#                           Identified_fusion = NA)
#for(i in n_first_sample:ncol(fusion_summary)){
#  fusion_temp <- fusion_summary[!(fusion_summary[,i] == "./.:.:.:.:.:.:.:."),]
#  if(nrow(fusion_temp) > 0){
#    fusion_temp <- fusion_temp[,c("AG_note", colnames(fusion_temp)[i])]
#    fusion_temp <- distinct(fusion_temp)
#    fusion <- rbind(fusion,
#                    cbind.data.frame(Sample = colnames(fusion_temp)[2],
#                                     Identified_fusion = fusion_temp$AG_note))
#  }
#}

#fusion <- na.omit(fusion)
#fusion <- distinct(fusion)

#n_dup <- grep("dup", fusion$Sample)
#fusion$Sample
#fusion$Sample[-n_dup] <- sapply(strsplit(fusion$Sample[-n_dup], split = "_"),
#                                function(x){paste0(x[[1]])})
#fusion$Sample[n_dup] <- sapply(strsplit(fusion$Sample[n_dup], split = "_"),
#                               function(x){paste0(x[[1]], "_", x[[2]])})
#fusion$Sample



####### ------- multiqc_table

#colnames(multiqc_table)

# In later versions where I pulled the multiQC table directly from the nf-core pipeline,
# the names changed a bit so I will convert them back
colnames(multiqc_table) <- gsub("FastQC..raw._mqc.generalstats.fastqc_raw", 
                                "FastQC_mqc.generalstats.fastqc", 
                                colnames(multiqc_table))
colnames(multiqc_table) <- gsub("FastP..Read.preprocessing._mqc.generalstats.fastp_read_preprocessing", 
                                "fastp_mqc.generalstats.fastp", 
                                colnames(multiqc_table))
colnames(multiqc_table) <- gsub("Samtools.Flagstat_mqc.generalstats.samtools_flagstat", 
                                "Samtools_mqc.generalstats.samtools", 
                                colnames(multiqc_table))
colnames(multiqc_table) <- gsub("Mosdepth_mqc.generalstats", 
                                "mosdepth_mqc.generalstats", 
                                colnames(multiqc_table))
colnames(multiqc_table) <- gsub("SNPeff_mqc.generalstats.snpeff", 
                                "SnpEff_mqc.generalstats.snpeff", 
                                colnames(multiqc_table))

# With the pipeline table, I am missing Picard_mqc.generalstats.picard.PERCENT_DUPLICATION


# Lengths should be equal to or double the number of samples
# Since different files are used, some names only match up to one type of metric
#length(grep("recal", multiqc_table$Sample.Name))
#length(grep(".recal.cram", multiqc_table$Sample.Name))
#length(grep("-L[0-9].0", multiqc_table$Sample.Name))
length(grep("recal", multiqc_table$Sample))
length(grep(".recal.cram", multiqc_table$Sample))
length(grep("-L[0-9].0", multiqc_table$Sample))

#qc_coverage <- cbind.data.frame(Sample = multiqc_table$Sample.Name,
#                                Fraction_genome_30x = multiqc_table$'≥.30X',
#                                Median_coverage_X = multiqc_table$Median)
qc_coverage <- cbind.data.frame(Sample = multiqc_table$Sample,
                                Fraction_genome_30x = multiqc_table$mosdepth_mqc.generalstats.mosdepth.30_x_pc,
                                Median_coverage_X = multiqc_table$mosdepth_mqc.generalstats.mosdepth.median_coverage)
qc_coverage <- na.omit(qc_coverage)
qc_coverage <- qc_coverage[grep("recal", qc_coverage$Sample),]
qc_coverage$Median_coverage_X <- gsub("X", "", qc_coverage$Median_coverage_X)
qc_coverage$Sample <- gsub(".recal", "", qc_coverage$Sample)

#qc_reads <- cbind.data.frame(Sample = multiqc_table$Sample.Name,
#                             Total_seqs_M = multiqc_table$M.Total.seqs,
#                             Fraction_reads_mapped = multiqc_table$'%.Mapped')
qc_reads <- cbind.data.frame(Sample = multiqc_table$Sample,
                             Total_seqs_M = multiqc_table$Samtools_mqc.generalstats.samtools.raw_total_sequences,
                              # Total_seqs_M = READS (so 1 pair is counted as two here)
                             Fraction_reads_mapped = multiqc_table$Samtools_mqc.generalstats.samtools.reads_mapped_percent)
qc_reads <- na.omit(qc_reads)
qc_reads <- qc_reads[grep("recal", qc_reads$Sample),]
qc_reads$Sample <- gsub(".recal.cram", "", qc_reads$Sample)

# Begin join
qc_summary <- join(qc_coverage,
                   qc_reads)

if("Picard_mqc.generalstats.picard.PERCENT_DUPLICATION" %in% colnames(multiqc_table)){
  message("Using Picard duplication metrics from multiQC")

  qc_dups <- cbind.data.frame(Sample = multiqc_table$Sample,
                              Fraction_duplicated_reads = multiqc_table$Picard_mqc.generalstats.picard.PERCENT_DUPLICATION)
  qc_dups <- na.omit(qc_dups)
  qc_dups$Sample <- sapply(strsplit(qc_dups$Sample, split = "-"),
                           function(x){x[[1]]})
  # Taking the average if there are multiple lanes:
  qc_dups <- ddply(qc_dups, "Sample", numcolwise(mean))

  # Join
  qc_summary <- join(qc_summary,
                     qc_dups)
} 

qc_summary <- join(qc_summary,
                   twist_coverage)


################################################################################
####### -------                 Join tables                      ------- #######

# Join
summary_table <- join(sample_table,
                      fusion_table)
if(any(is.na(summary_table$fusion))){
  warning("Confirm the sample names to merge are correct")
}
# A warning above is normal if not all samples have fusions

summary_table <- join(summary_table,
                      primer_summary,
                      by = "row", type = "left")
summary_table <- join(summary_table,
                      qc_summary)

# Replace NA with ""
summary_table2 <- sapply(summary_table, as.character)
summary_table2[is.na(summary_table2)] <- ""

# Ordering by post-TIDDIT filtering summary
summary_table2 <- as.data.frame(summary_table2)
summary_table2 <- summary_table2[order(summary_table2$post_TIDDIT_flag_summary,
                                       decreasing = T),]


################################################################################
####### -------                     Write tables                 ------- #######

# Write tables
write.table(summary_table2,
            file = paste0(output, "final_merged_table.txt"),
            sep = "\t", row.names = F, quote = F)

# Write as xlsx tables with an existing column description sheet
file.copy(column_descriptions,
          paste0(output, "final_merged_table.xlsx"),
          overwrite = T)
# A TRIMMED RESULTS TABLE
select_columns <- c("Sample",
                    "Study.ID",
                    "Cohort",
                    "row",
                    "fusion",
                    "gene1",
                    "gene2",
                    "post_TIDDIT_flag_summary",
                    "FILTER",
                    "P1_sequence_v1",
                    "P2_sequence_v1",
                    "P1_TM_v1",
                    "P2_TM_v1",
                    "P1_sequence_v2",
                    "P2_sequence_v2",
                    "P1_TM_v2",
                    "P2_TM_v2",
                    "p1_v1_status",
                    "p2_v1_status",
                    "p1_v2_status",
                    "p2_v2_status",
                    "primer_notes",
                    "Fraction_genome_30x",
                    "Median_coverage_X",
                    "Total_seqs_M",
                    "Fraction_reads_mapped",
                    "Fraction_duplicated_reads",
                    "MedianInside",
                    "PctAt10XInside",
                    "PctAt200XInside",
                    "MedianOutside")
if(!("Picard_mqc.generalstats.picard.PERCENT_DUPLICATION" %in% colnames(multiqc_table))){

  select_columns <- select_columns[select_columns != "Fraction_duplicated_reads"]
    
}
xlsx::write.xlsx(summary_table2[,select_columns], 
                 paste0(output, "final_merged_table.xlsx"),
                 sheetName="Select results", append=TRUE, row.names = F)
# FULL
xlsx::write.xlsx(summary_table2, 
                 paste0(output, "final_merged_table.xlsx"),
                 sheetName="Full results", append=TRUE, row.names = F)





