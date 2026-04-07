# Twist Fusion Pipeline: Getting Started Guide

This guide details the steps required to get the **Twist-enriched sequencing pipeline** running on new data. The pipeline is designed to detect gene fusions (specifically sarcoma translocations) using a combination of `nf-core/sarek` (TIDDIT caller) and custom R filtering scripts.

---

## 🏗️ Prerequisites

### Infrastructure

* **HPC Access**: Alpine/CU Research Computing (CURC) environment.
* **Software Modules**:
  * `nextflow/23.04`
  * `singularity/3.7.4`
  * `bcftools/1.16`, `htslib/1.16`, `bedtools/2.29.1`, `samtools`
* **Singularity Containers**: Pre-built containers for VCTools, R, and MultiQC (located in project reference paths).

### Reference Files

* **Genome**: GRCh38 (GATK version recommended).
* **Twist Bed File**: `Probes_merged_ok_CU-Hayashi_PedSolidTumorTransloc_TE-92680772_hg38_250404204242.bed` (Probe regions for target enrichment).

---

## 📂 Data Setup

### 1. Project Directory Structure

Create a new directory for your project and copy the pipeline scripts:

```bash
mkdir -p /pl/active/bbsrinformatics/analysts/$USER/projects/<new_project>
cd /pl/active/bbsrinformatics/analysts/$USER/projects/<new_project>
cp -r /Users/mkaufman/BBSR_Projects/Hayashi_Twist_Jan2026/pipeline/scripts .
```

### 2. Create `samples.csv`

Prepare a CSV file mapping your samples to their FastQ files.
Header: `patient,status,sample,lane,fastq_1,fastq_2`
> [!NOTE]
> `status=1` indicates tumor/affected samples.

### 3. Configure `config_loose.sh`

Update the following variables in `pipeline/scripts/config_loose.sh`:

* `samplefile`: Absolute path to your new `samples.csv`.
* `projbn`: Your new project base name.
* `plhome`: Your project home directory on Alpine.

---

## 🚀 Pipeline Execution

### Step 1: Run nf-core/sarek

Execute the main pipeline using `run_sarek_loose.sh`. This runs TIDDIT with "loose" parameters optimized for Twist coverage.

```bash
sbatch run_sarek_loose.sh
```

* **Output**: `nf-core_sarek_loose/`

### Step 2: Quality Control (Twist Specific)

Check the coverage specifically within the Twist probe regions.

```bash
sbatch twist_coverage.sh
```

* **Output**: `twist_coverage/summary_results.txt` (Contains MedianInside, PctAt10X, etc.)

### Step 3: VCF Filtering

Extract canonical fusions from the massive TIDDIT output.

```bash
sbatch filter_vcfs_loose_FUSIONS.sh
```

* **Output**: `vcfs_loose_FUSIONS/` (Contains filtered VCFs and a merged `canonical_fusions_no_header.vcf`)

---

## 📊 Downstream Analysis (R)

After the HPC steps, transfer the results to your local environment for final analysis.

### Required R Scripts

Located in `hayashi_twist_Dec2025/scripts/`:

1. **`filter_fusions.R`**:
    * **Action**: Update `setwd()` and `shared_files` path (likely `/Volumes/Partition 2/Core/hayashi_sarcoma_fusion_files/Data/`).
    * **Result**: Generates `final_summary_table.xlsx` with fusion calls and flags.
2. **`summary_table.R`**:
    * **Action**: Merges MultiQC results, Twist coverage, and fusion calls.
    * **Result**: Produces the final `final_merged_table.xlsx` for collaborator sharing.

### ⚠️ Common Pitfalls

* **Hardcoded Paths**: The R scripts contain hardcoded volume paths (`/Volumes/Partition 2/...`). These **MUST** be updated to match your local setup.
* **Sample Naming**: Ensure sample IDs in `samples.csv` match the naming convention used in the R join steps (often involves stripping suffixes like `_L1`).

---

## ✅ Expected Results

The final output is `Analysis/main_filtered_output/final_merged_table.xlsx`. This table includes:

* Sample metadata.
* Fusion identification (Gene 1, Gene 2, Breakpoints).
* Evidence metrics (Split reads, Paired reads).
* Quality flags (e.g., `PASS`, `low_variant_coverage`, `gene1_outside_twist`).
* Primer sequences for validation.
