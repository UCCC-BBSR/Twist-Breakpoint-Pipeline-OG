## Twist

The repository contains documentation for multiple projects related to **Twist-enriched sequencing** for gene fusion identification in sarcoma research, primarily for the Hayashi lab:

### Main Twist Projects:

1. **hayashi_twist_Dec2025** - The cleanest processing and analysis of Twist-enriched data, recommended for future similar work
2. **hayashi_cfDNA_twist_Jun2025** - cfDNA samples analyzed using the Twist assay
3. **hayashi_twist_design_Feb2025** - Custom Twist fusion panel design covering canonical sarcoma-related translocations
4. **hayashi_twist_inconsistent_samples_Dec2025** - Investigation of inconsistencies between bioinformatics-predicted fusions and lab results
5. **hayashi_sarcoma_twist_enriched_Aug2025** - Sarcoma samples with Twist enrichment
6. **hayashi_sarcoma_twist_enriched_removeduplicates_Aug2025** - Testing duplicate removal impact

### Key Details:

- The Twist assay is a **custom panel for detecting gene fusions** in sarcomas
- Projects include WGS and Twist-enriched sequencing data
- Analysis involves TIDDIT fusion caller with modified parameters for Twist coverage
- Twist-enriched QC includes additional coverage metrics in expected regions
- Part of the "Identification of gene fusions" project collection

## Twist Pipeline Overview

The **hayashi_twist_Dec2025** project is documented as "**the cleanest processing and analysis of Twist-enriched data**" and should be used as the reference for future work.

### Pipeline Steps:

#### 1. **Initial Processing with nf-core/sarek**

- Run the **nf-core/sarek pipeline** on raw sequencing data (example script: `run_sarek_loose.sh`)
- Uses **TIDDIT** as the main fusion caller with **SNPEff** for annotation
- For Twist-enriched data, TIDDIT parameters need to be **loosened** to prevent enriched coverage from interfering with fusion detection

#### 2. **Quality Control**

- Transfer nf-core multiqc report locally
- For Twist-enriched samples, run **`twist_coverage. sh`** to check coverage in expected regions (specific to Twist assay)

#### 3. **VCF Filtering (on Alpine/HPC)**

- Filter VCF files using script like **`filter_vcfs_loose_FUSIONS.sh`**
- Filters for annotated canonical fusions

#### 4. **Downstream Analysis (Local)**

Scripts are typically in `/Volumes/Partition 2/Core/<project_name>/`:

- **`filter_fusions. R`**:
  - Filters to canonical fusions
  - Constructs IGV ranges and bed files of breakpoints
  - Runs primer3 for PCR-testable breakpoints
  - For Twist data: provides post-TIDDIT thresholds to fusion calls
- **`subset_bam.sh`** (on Alpine): Subsets BAM files to discordant reads for easier breakpoint viewing
- **`make_IGV_xml.R`**: Generates IGV session files for visualizing breakpoints
- **`summary_table.R`**: Creates final summary table
  - Output: `<project>/Analysis/main_filtered_output/final_merged_table.xlsx`
  - This is shared with the Hayashi lab

#### 5. **Reference Files**

Located in `hayashi_sarcoma_fusion_files/`:

- Twist bed files (probe regions)
- Primer3 templates
- IGV session templates
- Column descriptions for summary tables

### Key Documentation

The most current workflow is documented in:

Code

```
/Volumes/Partition 2/Core/hayashi_twist_Dec2025/outline.docx
```

For complete methodology, see:

Code

```
overview_documents/masanori_hayashi.md
```

### Data Locations

- **Raw data**: `/pl/active/bbsrinformaticsdata/bbsr/Hayashi/<run_folder>/`
- **Analysis scripts**: `/pl/active/bbsrinformatics/analysts/andrew/core/<project_name>/` (Alpine)
- **Local analysis**: `/Volumes/Partition 2/Core/<project_name>/`

---

2026-01-28