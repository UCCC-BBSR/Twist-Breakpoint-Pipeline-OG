#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=72G
#SBATCH --time=72:00:00
#SBATCH --partition=amilan
#SBATCH --error=job_outs/filter_vcfs-%j.err  
#SBATCH --output=job_outs/filter_vcfs-%j.out
#SBATCH --qos long


### Differences from previous work
# Very different from previous work
# I have been having trouble merging the loosened version of these files because of their size 
# Here, I will mimic my R script by focusing on alterations with canonically annotated fusions (e.g. PAX3&NCOA1)
# I am still keeping the other filtered versions of the files but there will not be a merged unfiltered as I normally was able to do

### Load module
module load bcftools/1.16
module load htslib/1.16

### Source config to get paths
. config_loose.sh # Using $project, $pipeoutdir
path_nfcore_vcfs=${pipeoutdir}/annotation/tiddit

### Additional paths
nfcore_suffix=tiddit_snpEff.ann.vcf.gz
twist_bed=${project}/"Probes_merged_ok_CU-Hayashi_PedSolidTumorTransloc_TE-92680772_hg38_250404204242.bed"
# Major change to filtering strategy because the vcf files are too large with loose tiddit params to merge                                              
# Adding these new types of files as sarek_tiddit_snpeff_twist_region    

### Output paths
main_output=${project}/vcfs_loose_FUSIONS
mkdir ${main_output}
input=${main_output}/sarek_tiddit_snpeff # This will be the location for copied vcf files
mkdir ${input}

# **Description of output folders**
# sarek_tiddit_snpeff                                                Unfiltered VCF
# sarek_tiddit_snpeff_no_header                                      Unfiltered VCF without the header
# sarek_tiddit_snpeff_no_header_INV_chrX                             Unfiltered VCF without the header subset to ChrX inversions (including inverted duplications)
# sarek_tiddit_snpeff_PASS                                           Variants that pass TIDDIT's guidelines. Header contains meta information 
# sarek_tiddit_snpeff_PASS_AllSarcomaGenes                           sarek_tiddit_snpeff_PASS & contains a sarcoma-related fusion gene (all variants, not limited to fusions)
# sarek_tiddit_snpeff_PASS_BND                                       sarek_tiddit_snpeff_PASS & contains a break end alteration (BND, not limited to gene fusions)
# sarek_tiddit_snpeff_PASS_genefusion                                sarek_tiddit_snpeff_PASS & contains a gene_fusion
# sarek_tiddit_snpeff_PASS_BND_InterChrSarcomaGenes                  sarek_tiddit_snpeff_PASS_BND & contains an INTERCHROMOSOME sarcoma-related fusion gene
# sarek_tiddit_snpeff_PASS_genefusion_AllSarcomaGenes                sarek_tiddit_snpeff_PASS_genefusion & contains a sarcoma-related fusion gene
# sarek_tiddit_snpeff_PASS_IntraChrSarcomaGenes                      sarek_tiddit_snpeff_PASS & contains BOTH BCOR & CCNB3 
# sarek_tiddit_snpeff_PASS_BND_InterChrSarcomaGenes_ProteinCoding    sarek_tiddit_snpeff_PASS_BND_InterChrSarcomaGenes & contains protein_coding
# sarek_tiddit_snpeff_PASS_INV_chrX                                  sarek_tiddit_snpeff_PASS & contains an inversion on chromosome X
# merged_vcf                                                         Contains versions of the above in a merged VCF format

snpeff_unfiltered_no_header=${main_output}/sarek_tiddit_snpeff_unfiltered_no_header
snpeff_unfiltered_no_header_INV_chrX=${main_output}/sarek_tiddit_snpeff_unfiltered_no_header_INV_chrX
snpeff_PASS=${main_output}/sarek_tiddit_snpeff_PASS
snpeff_PASS_AllSarcomaGenes=${main_output}/sarek_tiddit_snpeff_PASS_SarcomaGenes
snpeff_PASS_BND=${main_output}/sarek_tiddit_snpeff_PASS_BND
snpeff_PASS_genefusion=${main_output}/sarek_tiddit_snpeff_PASS_genefusion
snpeff_PASS_BND_InterChrSarcomaGenes=${main_output}/sarek_tiddit_snpeff_PASS_BND_InterChrSarcomaGenes
snpeff_PASS_genefusion_AllSarcomaGenes=${main_output}/sarek_tiddit_snpeff_PASS_genefusion_SarcomaGenes
snpeff_PASS_IntraChrSarcomaGenes=${main_output}/sarek_tiddit_snpeff_PASS_IntraChrSarcomaGenes
snpeff_PASS_BND_InterChrSarcomaGenes_ProteinCoding=${main_output}/sarek_tiddit_snpeff_PASS_BND_InterChrSarcomaGenes_ProteinCoding
snpeff_PASS_INV_chrX=${main_output}/sarek_tiddit_snpeff_PASS_INV_chrX
merged_vcf=${main_output}/merged_vcfs
sarek_tiddit_snpeff_twist_region=${main_output}/sarek_tiddit_snpeff_twist_region
snpeff_unfiltered_no_header_FUSIONS=${main_output}/sarek_tiddit_snpeff_unfiltered_FUSIONS

mkdir ${snpeff_unfiltered_no_header} ${snpeff_unfiltered_no_header_INV_chrX} ${snpeff_PASS} ${snpeff_PASS_AllSarcomaGenes} \
      ${snpeff_PASS_BND} ${snpeff_PASS_genefusion} ${snpeff_PASS_BND_InterChrSarcomaGenes} ${snpeff_PASS_genefusion_AllSarcomaGenes} \
      ${snpeff_PASS_IntraChrSarcomaGenes} ${snpeff_PASS_BND_InterChrSarcomaGenes_ProteinCoding} ${snpeff_PASS_INV_chrX} \
      ${sarek_tiddit_snpeff_twist_region} ${snpeff_unfiltered_no_header_FUSIONS} ${merged_vcf}


### Sarcoma fusions
fusions_canonical_order="BCOR&CCNB3|CIC&DUX4|CIC&FOXO4|EWSR1&FLI1|EWSR1&ERG|EWSR1&WT1|PAX3&FOXO1|PAX7&FOXO1|SS18&SSX1|SS18&SSX2|SS18&SSX4B|ASPSCR1&TFE3|PAX3&MAML3|ETV6&NTRK3|TPM3&NTRK1|TEAD1&NCOA2|VGLL2&NCOA2|VGLL2&CITED2|EWSR1&NR4A3|TAF15&NR4A3|TCF12&NR4A3|TFG&NR4A3|HEY1&NCOA2|EWSR1&ETV1|EWSR1&ETV4|EWSR1&FEV|EWSR1&ATF1|EWSR1&DDIT3|EWSR1&TFCP2|EWSR1&CREB1|EWSR1&NFATC2|FUS&TFCP2|FUS&ERG|FUS&FEV|FUS&DDIT3|EML4&NTRK3|LMNA&NTRK1|PAX3&NCOA1"
fusions_reversed_order="CCNB3&BCOR|DUX4&CIC|FOXO4&CIC|FLI1&EWSR1|ERG&EWSR1|WT1&EWSR1|FOXO1&PAX3|FOXO1&PAX7|SSX1&SS18|SSX2&SS18|SSX4B&SS18|TFE3&ASPSCR1|MAML3&PAX3|NTRK3&ETV6|NTRK1&TPM3|NCOA2&TEAD1|NCOA2&VGLL2|CITED2&VGLL2|NR4A3&EWSR1|NR4A3&TAF15|NR4A3&TCF12|NR4A3&TFG|NCOA2&HEY1|ETV1&EWSR1|ETV4&EWSR1|FEV&EWSR1|ATF1&EWSR1|DDIT3&EWSR1|TFCP2&EWSR1|CREB1&EWSR1|NFATC2&EWSR1|TFCP2&FUS|ERG&FUS|FEV&FUS|DDIT3&FUS|NTRK3&EML4|NTRK1&LMNA|NCOA1&PAX3"


### Sarcoma genes
sarcoma_genes="BCOR|CIC|EWSR1|PAX3|PAX7|SS18|ASPSCR1|ETV6|TPM3|TEAD1|VGLL2|TAF15|TCF12|TFG|HEY1|FUS|EML4|LMNA|CCNB3|DUX4|FOXO4|FLI1|ERG|WT1|FOXO1|SSX1|SSX2|SSX4B|TFE3|MAML3|NTRK3|NTRK1|NCOA2|CITED2|NR4A3|ETV1|ETV4|FEV|ATF1|DDIT3|TFCP2|CREB1|NFATC2|NCOA1"
sarcoma_ensembl="ENSG00000183337|ENSG00000079432|ENSG00000182944|ENSG00000135903|ENSG00000009709|ENSG00000141380|ENSG00000169696|ENSG00000139083|ENSG00000143549|ENSG00000187079|ENSG00000170162|ENSG00000276833|ENSG00000140262|ENSG00000114354|ENSG00000164683|ENSG00000089280|ENSG00000143924|ENSG00000160789|ENSG00000147082|ENSG00000283949|ENSG00000184481|ENSG00000151702|ENSG00000157554|ENSG00000184937|ENSG00000150907|ENSG00000126752|ENSG00000241476|ENSG00000269791|ENSG00000068323|ENSG00000196782|ENSG00000140538|ENSG00000198400|ENSG00000140396|ENSG00000164442|ENSG00000119508|ENSG00000006468|ENSG00000175832|ENSG00000163497|ENSG00000123268|ENSG00000175197|ENSG00000135457|ENSG00000118260|ENSG00000101096|ENSG00000084676"
sarcoma_genes_interchromosomal="CIC|EWSR1|PAX3|PAX7|SS18|ASPSCR1|ETV6|TEAD1|VGLL2|TAF15|TCF12|TFG|HEY1|FUS|EML4|DUX4|FOXO4|FLI1|ERG|WT1|FOXO1|SSX1|SSX2|SSX4B|TFE3|MAML3|NTRK3|NCOA2|NR4A3|ETV1|ETV4|FEV|ATF1|DDIT3|TFCP2|CREB1|NFATC2"
sarcoma_ensembl_interchromosomal="ENSG00000079432|ENSG00000182944|ENSG00000135903|ENSG00000009709|ENSG00000141380|ENSG00000169696|ENSG00000139083|ENSG00000187079|ENSG00000170162|ENSG00000276833|ENSG00000140262|ENSG00000114354|ENSG00000164683|ENSG00000089280|ENSG00000143924|ENSG00000283949|ENSG00000184481|ENSG00000151702|ENSG00000157554|ENSG00000184937|ENSG00000150907|ENSG00000126752|ENSG00000241476|ENSG00000269791|ENSG00000068323|ENSG00000196782|ENSG00000140538|ENSG00000140396|ENSG00000119508|ENSG00000006468|ENSG00000175832|ENSG00000163497|ENSG00000123268|ENSG00000175197|ENSG00000135457|ENSG00000118260|ENSG00000101096"
# Not using below but keeping for clarity
sarcoma_genes_intrachromosomal="BCOR|TPM3|VGLL2|LMNA|PAX3|CCNB3|NTRK1|CITED2|NCOA1"
sarcoma_ensembl_intrachromosomal="ENSG00000183337|ENSG00000143549|ENSG00000170162|ENSG00000160789|ENSG00000135903|ENSG00000147082|ENSG00000198400|ENSG00000164442|ENSG00000084676"


### Copy vcfs from nf-core pipeline
rsync -aP ${path_nfcore_vcfs}/*/*${nfcore_suffix} ${input}/


### Get base names of vcd files from $input
vcfs=( "$input"/*.vcf.gz )  # load literal filenames
vcfs=( "${vcfs[@]##*/}" )   # strip off directory names
vcfs=( "${vcfs[@]%.gz}" )   # strip off .gz extension


### Filter sample-specific VCFs
for vcf in "${vcfs[@]}"
do
   zgrep -E "#ALT|#INFO|#TIDDIT|#FILTER|#FORMAT|#LibraryStats|#SnpEff|#fileformat|#source|PASS" "$input"/${vcf}.gz > "$snpeff_PASS"/${vcf}
   zgrep -v "##" "$input"/${vcf}.gz > "$snpeff_unfiltered_no_header"/${vcf}
   grep -E "#|SVTYPE=INV|SVTYPE=DUP:INV" "$snpeff_unfiltered_no_header"/${vcf} | grep -E "#|chrX" > "$snpeff_unfiltered_no_header_INV_chrX"/${vcf}
   grep -E "#|""$sarcoma_ensembl" "$snpeff_PASS"/${vcf} > "$snpeff_PASS_AllSarcomaGenes"/${vcf}
   grep -E "#|SVTYPE=BND" "$snpeff_PASS"/${vcf} > "$snpeff_PASS_BND"/${vcf}
   grep -E "#|SVTYPE=INV" "$snpeff_PASS"/${vcf} | grep -E "#|chrX" > "$snpeff_PASS_INV_chrX"/${vcf}
   grep -E "#|gene_fusion" "$snpeff_PASS"/${vcf} > "$snpeff_PASS_genefusion"/${vcf}
   grep -E "#|""$sarcoma_genes_interchromosomal" "$snpeff_PASS_BND"/${vcf} > "$snpeff_PASS_BND_InterChrSarcomaGenes"/${vcf}
   grep -E "#|""$sarcoma_genes" "$snpeff_PASS_genefusion"/${vcf} > "$snpeff_PASS_genefusion_AllSarcomaGenes"/${vcf}
   grep -E "#|ENSG00000183337" "$snpeff_PASS"/${vcf} | grep -E "#|ENSG00000147082" > "$snpeff_PASS_IntraChrSarcomaGenes"/${vcf}
   grep -E "#|protein_coding" "$snpeff_PASS_BND_InterChrSarcomaGenes"/${vcf} > "$snpeff_PASS_BND_InterChrSarcomaGenes_ProteinCoding"/${vcf}
   #bedtools intersect -a "$input"/${vcf}.gz -b ${twist_bed} -header -wa | gunzip -c > "$sarek_tiddit_snpeff_twist_region"/${vcf}
   #bedtools intersect -a "$input"/${vcf}.gz -b ${twist_bed} > "$sarek_tiddit_snpeff_twist_region"/${vcf}
   #cp "$sarek_tiddit_snpeff_twist_region"/${vcf} "$sarek_tiddit_snpeff_twist_region"/SAVE_${vcf}

   #grep -E "$fusions_canonical_order"|"$fusions_reversed_order" "$snpeff_unfiltered_no_header"/${vcf} > ${snpeff_unfiltered_no_header_FUSIONS}/${vcf}

#   grep -E "($fusions_canonical_order|$fusions_reversed_order)" \
#  "$input"/${vcf}.gz \
#  > "$snpeff_unfiltered_no_header_FUSIONS/${vcf}"

   zgrep -E "^(##|#CHROM)|($fusions_canonical_order|$fusions_reversed_order)" \
  "$input/${vcf}.gz" \
  > "$snpeff_unfiltered_no_header_FUSIONS/${vcf}"
   
done

### Create a merged VCF to perform the same filtering
cd "${snpeff_unfiltered_no_header_FUSIONS}"
#gunzip *.vcf.gz
for F in *.vcf ; do bgzip ${F} ; done
for i in *.vcf.gz ; do bcftools index -t -f $i ; done
bcftools merge *.vcf.gz > "$merged_vcf"/"canonical_fusions_w_header.vcf"
grep -v "##" "$merged_vcf"/"canonical_fusions_w_header.vcf" > "$merged_vcf"/"canonical_fusions_no_header.vcf"
