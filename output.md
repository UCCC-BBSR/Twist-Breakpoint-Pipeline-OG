**Goals**

Identify if canonical sarcoma fusions are present in these
twist-enriched samples

**Background**

I am following work performed in \`hayashi_cfDNA_twist_Jun2025\` but I
am reorganizing the directory on petalibrary to make it a bit more
organized

Alpine directory:
/pl/active/bbsrinformatics/analysts/andrew/core/hayashi_twist_Dec2025

Local analysis directory: /Volumes/Partition\\
2/Core/hayashi_twist_Dec2025

This directory should contain:

./outline.docx (this file)

./samples.xlsx

> (sample information, see /Volumes/Partition
> 2/Core/hayashi_twist_Dec2025/samples.xlsx)

./Analysis (empty folder)

./Scripts

Location of local shared files:

\<shared_files\> = /Volumes/Partition
2/Core/hayashi_sarcoma_fusion_files/Data

**nf-core/sarek**

Edit scripts/config_loose.sh (sourced by several shell scripts on
Alpine)

Run script run_sarek_loose.sh

[Location of important files]{.underline}

QC reports = \<outdir\>/reports/\*

Ideally, multiqc was run within here but I had to run separate

Annotated TIDDIT VCFs = \<outdir\>/annotation/tiddit/\*

Sample bam files = \<outdir\>/preprocessing/recalibrated/\*

**Sample-level quality control**

[Transfer nf-core-generated multiqc]{.underline}

mkdir -p \"/Volumes/Partition
2/Core/hayashi_twist_Dec2025/Analysis/QC_reports\" \\

&& cd \"/Volumes/Partition
2/Core/hayashi_twist_Dec2025/Analysis/QC_reports\"

rsync -aP
aeg22@xsede.org@login-ci.rc.colorado.edu:/pl/active/bbsrinformatics/analysts/andrew/core/hayashi_twist_Dec2025/nf-core_sarek_loose/multiqc/multiqc_report.html
.

rsync -aP
aeg22@xsede.org@login-ci.rc.colorado.edu:/pl/active/bbsrinformatics/analysts/andrew/core/hayashi_twist_Dec2025/nf-core_sarek_loose/multiqc/multiqc_data/multiqc_general_stats.txt
.

\# In older examples, I needed to re-run multiqc manually

\# I will need to change some column names in this file

While I can look at QC here using the html file or excel table, I
extract important metrics and combine with fusions later.

[Notes on QC from multiqc at this stage:]{.underline}

The QC metrics when using the Twist panel will look pretty poor because
the whole genome is considered -- I calculate coverage for just the
Twist region later

All samples are 2.2 -- 6.2M reads except for 039 and 158, which are 0.2M
and 0.7M, respectively. For all samples, there are far fewer reads than
I have previously seen.

Adapter sequences are a bit higher than I have previously seen. By
basepair 140, about 20% of reads have adapter where as it was \<12% for
all samples in a previous experiment

[Calculate coverage based on Twist target regions]{.underline}

twist_coverage.sh

Copy to local

rsync -aP
aeg22@xsede.org@login-ci.rc.colorado.edu:/pl/active/bbsrinformatics/analysts/andrew/core/hayashi_twist_Dec2025/twist_coverage/summary_results.txt
/Volumes/Partition\\
2/Core/hayashi_twist_Dec2025/Analysis/QC_reports/twist_coverage_summary_results.txt

[Notes on twist coverage]{.underline}

Median coverage inside twist region is \> 200 except for two samples:
039 and 158 which only have a median coverage of 0X and 44X,
respectively

Note that in a previous twist experiment, median coverage was \> 5000
for all samples so all samples in this dataset is way lower

**Filter VCF files**

The following shell script employs many types of filtering. In the end,
I will likely only use one or two files but they are helpful to keep
around to be able to inspect specific samples under different filtering
conditions

[filter_vcfs_loose_FUSIONS.sh (run on alpine)]{.underline}

This script is much different than I have run for non-Twist data. Here,
I will mimic my R script by focusing on alterations with canonically
annotated fusions (e.g. PAX3&NCOA1) I am still keeping the other
filtered versions of the files but there will not be a merged unfiltered
as I normally was able to do

Copy vcfs to local hard drive:

rsync -aP
aeg22@xsede.org@login-ci.rc.colorado.edu:/pl/active/bbsrinformatics/analysts/andrew/core/hayashi_twist_Dec2025/vcfs_loose_FUSIONS/merged_vcfs
/Volumes/Partition\\
2/Core/hayashi_twist_Dec2025/Analysis/vcfs_loose_FUSIONS/

**~~Filter VCF files for non-canonical gene fusions~~**

~~\# I am interested in fusions we would never expect to see because
these might~~

~~\# make ideal false positive thresholds for each sample\'s canonical
fusions to~~

~~\# pass over~~

~~I technically ran this in
hayashi_sarcoma_twist_enriched_removeduplicates_Aug2025 but that is fine
because the TIDDIT results were the same~~

**~~[filter_vcfs_loose_NONCANONICALFUSIONS.sh]{.underline}~~**

~~This script is much different than I have run for non-Twist data.
Here, I will mimic my R script by focusing on alterations with
**NON-**canonically annotated fusions.~~

**~~Copy local:~~**

~~rsync -aP
aeg22@xsede.org@login-ci.rc.colorado.edu:/pl/active/bbsrinformatics/analysts/andrew/core/hayashi_twist_Dec2025/vcfs_loose_NONCANONICALFUSIONS/merged_vcfs
/Volumes/Partition\\
2/Core/hayashi_twist_Dec2025/Analysis/vcfs_loose_NONCANONICALFUSIONS/~~

~~Process locally in R~~

~~[filter_noncanonical_fusions.R]{.underline}~~

~~The most important output is printed to console before the final lines
of the script~~

~~For summary statistics, see: /Volumes/Partition
2/Core/hayashi_sarcoma_twist_enriched_Aug2025/noncanonical_fusions_summary_statistics.xlsx~~

**Generate a table of sample fusions**

This R script uses filtered files from filter_vcf.sh to extract
canonical inter- and intrachromosome gene fusions expected in these
sarcoma samples. I am discarded alterations that are not these canonical
fusions.

[filter_fusions.R]{.underline}

Located here: \<local_directory\>/scripts/filter_fusions.R

Reads in several files placed in \<shared_files\> as well as the sample
info file

Incorporating post TIDDIT filters:

#post_TIDDIT_flag_gene1_outside_twist \# Is the gene1 breakpoint outside
of canonical fusion region?

#post_TIDDIT_flag_gene2_outside_twist \# Is the gene2 breakpoint outside
of canonical fusion region?

#post_TIDDIT_flag_low_split_read_coverage \# Is the split read coverage
\< 1?

#post_TIDDIT_flag_low_paired_read_coverage \# Is the paired read
coverage \< 1?

#post_TIDDIT_flag_low_variant_coverage \# Is the aggregate fusion read
evidence \< 5?

#post_TIDDIT_flag_low_MAPQ \# Is the MAPQ \< 20?

#post_TIDDIT_flag_identical_fusion_in_multiple_samples \# Is the
identical fusion found in multiple samples?

[Things the tail end of this script does:]{.underline}

1.  I think fusion_partners should be altered to have gene1 and gene2 in
    the same/expected order

2.  Bind gene1 and gene2 to filtered table

3.  Keep only the VCF row where gene1 is first

4.  Construct IGV ranges

5.  Bed file of IGV ranges

6.  Bed files of possible primer regions for IGV

    a.  To deal with some deletions I have already seen, I will include
        2 primer set options (v1 and v2) based on how close they are to
        the breakpoint

    b.  I need to extract sequences from a local fasta file (which I
        don't have yet) and some need to be the reverse complement (any
        gene where that is on the negative strand)

    c.  Could this be a single bed file that spans all samples?

        i.  Might also consider adding the primer positions to this file
            as well (next step, below)

7.  Generate inputs for primer3 and calculate R1 and R2 primers

    a.  Looks like I need to extract a chimeric sequence

    b.  Also generate a bed file of the possible primer regions I chose
        (will be included in igv xml file below)

8.  Run primer3

    a.  In addition to the primer sequences, can I also extract the
        chromosomal positions and add these to a bed file as well?

**Subset bam files by IGV regions**

After filter_fusion.R, I have suspected IGV regions for the aggregate of
fusions. To minimize the pain and unnecessary reads of large bam files,
I will subset the bam files only the IGV regions.

Move main_filtered_output/primer_ranges.bed to Alpine

rsync -aP /Volumes/Partition\\
2/Core/hayashi_twist_Dec2025/Analysis/main_filtered_output/igv_ranges.bed
aeg22@xsede.org@login-ci.rc.colorado.edu:/pl/active/bbsrinformatics/analysts/andrew/core/hayashi_twist_Dec2025/

**Run a script to subset bam files on Alpine:**

**subset_bam.sh**

Move the subset bam files local:

rsync -aP
aeg22@xsede.org@login-ci.rc.colorado.edu:/pl/active/bbsrinformatics/analysts/andrew/core/hayashi_twist_Dec2025/bam_subset_igv_discordant/\*
/Volumes/Partition\\
2/Core/hayashi_twist_Dec2025/Analysis/igv_subset_bams/

**Generate IGV xml files**

[make_IGV_xml.R]{.underline}

Located here: /Volumes/Partition
2/Core/hayashi_twist_Dec2025/Scripts/make_IGV_xml.R

\# Samples with multiple fusions will have

\# \'\_fusionname_rowX\' appended to their file name. Samples with one
fusion will

\# just have the fusion name appended (without \_rowX) and samples
without any fusions

\# will just be their sample name. \'rowX\' corresponds to the row
column in the

\# fusion table.

**Inspect fusions in IGV**

Open the xml files within /Volumes/Partition\\
2/Core/hayashi_twist_Dec2025/Analysis/igv_subset_bams/

using IGV_2.16.2.app

Record status of primers in
\<output\>/primer_summary_post_inspection.xlsx

- This file was created in a previous R script

- **I skipped this step for now**

The goal is to identify any primers (or pairs) that should not be used
because of a likely translocation deletion, too close to the breakpoint,
consistent mutations, or low coverage that could also indicate the first
category.

The xml files are written with relative paths so the direction
'igv_subset_bams' can be shared and the xml files will work properly.

**Summarize QC and identified fusions**

Table results from several of the scripts above

[summary_table.R]{.underline}

Located here: /Volumes/Partition
2/Core/hayashi_twist_Dec2025/Scripts/summary_table.R

Takes in the QC table, sample info (format was not described in this doc
-- only needs 3 columns: patient, Study.ID, Fusion), primer summary, and
fusion output from filter_fusions.R and combined them into a summary
table.

**Share results**

The WGS data was processed using the nf-core/sarek pipeline and VCF
files after TIDDIT and snpEff were filtered and summarized to canonical
sarcoma fusions.

final_merged_table.xlsx

- This file summarizes sample- and fusion-level information across the
  categories described in the first sheet

- The results table is displayed as a subset of the most useful data as
  well as the full table

- Sample-level information includes cohorts and WGS quality control

- The fusion-level information often corresponds to a single sample
  although some samples may have multiple fusions. The 'row' column
  identifies unique fusion events.

Within the igv_subset_bams/ folder, there is a subfolder named
IGV_files/. Within that folder, there are xml files that can be opened
as [IGV](https://igv.org/doc/desktop/#DownloadPage/) sessions (might
need specifically version 2.16.2). These files describe the paths to
each sample's bam files and amplicon and primer regions. The xml files
will streamline the process of reading in these files and zooming into
each sample's breakpoint regions.

**Subsequent emails**

**Sent on 12/11**

Hi Mas,

Attached are the canonical fusions results for the newest Twist dataset.

Most of the samples with "PASS" in the column "post_TIDDIT_flag_summary"
should be considered good candidates.

There are some exceptions:

- 062-UAB has two EWSR1-FLI1 breakpoints that pass filters. Both look
  decent but I think the row 42 fusion has better support.

- 010-MCC-20339 also has two FOXO1-PAX3 breakpoints that pass filters
  but they are right next to each other so I don't think that is a large
  concern. Both primer sets (row 48 or 49) should cover the same regions
  to test

- The 136-UM EWSR1-FLI1 row 28 fusion did not pass filters because it
  was an identical breakpoint as another sample (010-MCC-20339 row 27).
  However, that second sample had a much better FOXO1-PAX3 fusion and
  the EWSR1-FLI1 only had 3 reads of paired support and 0 split reads.
  So, this isn't like the older dataset where we saw the same fusion in
  all samples so I would consider the 136-UM EWSR1-FLI1 row 28 fusion
  good to include.

- There is one fusion in sample 145-UM (EWSR1&FLI1). It did not pass
  filters because of low coverage but there is both split and paired
  read evidence and the reference coverage at this region is also low. I
  think this fusion could be included in a test.

If 136-UM and 145-UM are included, 11/16 samples had a canonical fusion
and 2/5 that did not had very low total reads.

Additional notes after I went through the inconsistent fusion data from
previous datasets:

- 010-MCC-20339 has the plateau structure with seemingly 2 breakpoints
  in each gene

- I would place the 078-CCM breakpoint over by about 100bp but the
  primer sets should still work

Let me know if you want me to send the IGV-related data so you can view
each of the breakpoints yourself,
