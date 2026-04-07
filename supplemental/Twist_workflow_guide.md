The Andrew Goodspeed Twist analysis pipeline is a specialized bioinformatics workflow designed for the identification of gene fusions from Twist-enriched DNA sequencing data (sarcoma-related). It builds upon standard WGS (Whole Genome Sequencing) methods but is optimized for the high-coverage, targeted nature of Twist assays.

### **1. Goals and Basic Steps**

**Goals:**

- **Identification of Gene Fusions:** Detect canonical and novel gene fusions (translocations) specifically in sarcoma and cfDNA (cell-free DNA) samples.
- **High Sensitivity:** Utilize "loosened" calling parameters to ensure canonical fusions aren't missed due to the specific coverage patterns of Twist-enrichment.
- **Validation Readiness:** Provide actionable data (IGV sessions, PCR primers) for lab-based verification of predicted breakpoints.

**Basic Steps:**

1. **Alignment & Quality Control:** Standard processing via `nf-core/sarek` with additional Twist-specific metrics.
2. **Structural Variant Calling:** Running `TIDDIT` with relaxed parameters.
3. **Targeted Filtering:** Filtering VCFs for annotated canonical fusions and known hotspots.
4. **Downstream Annotation:** Constructing summary tables, IGV viewing sessions, and generating PCR primer sequences for lab validation.

------

### **2. Workflow Diagram**

Mermaid diagram code:

```
graph TD
    A[Raw FASTQ Data] --> B[nf-core/sarek Pipeline]
    B --> C{QC Metrics}
    C -->|Twist Specific| D[Coverage in Expected Regions]
    B --> E[BAM Alignment]
    E --> F[Subset BAMs to Discordant Reads]
    F --> G[TIDDIT SV Calling]
    G --> H[VCF Filtering]
    H --> I[Annotation with SNPEff]
    I --> J[Local R Analysis]
    J --> K[Final Summary Table]
    J --> L[IGV Session XML]
    J --> M[Primer3 PCR Designs]
```

Diagram:

![twist-workflow](/Users/mkaufman/BBSR_Projects/Hayashi_Twist_Jan2026/supplemental/twist-workflow.png)

### **3. Detailed Process Steps**

- **nf-core/sarek Alignment:**
  - Initial processing of raw data using the `nf-core/sarek` pipeline.
  - Includes duplicate removal (though tests in `hayashi_sarcoma_twist_enriched_removeduplicates_Aug2025` suggested this has a minor impact on results).
- **Twist-Enriched QC:**
  - Running `twist_coverage.sh` to assess coverage specifically within the target regions of the custom Twist panel.
  - Standard MultiQC reports are generated to monitor alignment quality.
- **Relaxed TIDDIT Calling:**
  - TIDDIT is used as the primary structural variant (SV) caller.
  - Parameters are explicitly "loosened" (relaxed) to prevent the Twist-enrichment coverage spikes from interfering with the variant calling algorithm.
- **BAM Subsetting:**
  - Running `subset_bam.sh` on the cluster (Alpine) to extract discordant reads. This creates smaller BAM files that allow for faster and easier manual inspection of breakpoints.
- **VCF Filtering & Annotation:**
  - Filtering is performed using scripts like `filter_vcfs_loose_FUSIONS.sh`.
  - Variants are annotated using **SNPEff** to identify which genes are involved in potential fusions.
- **Downstream R Analysis (`filter_fusions.R`):**
  - Filters results to focus on canonical fusions.
  - Calculates post-TIDDIT thresholds to rank the probability of fusions.
  - Constructs BED files and IGV ranges for the breakpoints.
  - **Primer3 Integration:** Automatically generates primer sequences for the lab to test the predicted breakpoints via PCR.
- **Reporting and Visualization:**
  - **IGV Sessions:** `make_IGV_xml.R` generates XML files so researchers can quickly view the breakpoints in the IGV browser.
  - **Summary Table:** `summary_table.R` produces the `final_merged_table.xlsx`, which is the primary deliverable shared with the lab.

*Note: For the most current implementation of this pipeline, the repository recommends referencing the `hayashi_twist_Dec2025` project as it contains the "cleanest" version of the processing logic.*

