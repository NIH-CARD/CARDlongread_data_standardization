#!/bin/bash
# initial script for cleaning up input variants, whether SV or SNV
# load necessary modules
module load plink bcftools
# print help statements
helpFunction()
{
   echo "Script for initial variant cleanup pre-allele specific QTL, using plink and bcftools. Splits multiallelic variants by default."
   echo "Usage: $0 -i input_prefix -s sample_list -a sample_list_action -m maf_cutoff -g missing_genotype_rate -h hwe_pvalue -p indep_pairwise_ld_pruning_values -t threads"
   echo -e "\t-i Input VCF file prefix"
   echo -e "\t-s Path to list of samples to include or exclude"
   echo -e "\t-a Action to take on sample list (either include or exclude)"
   echo -e "\t-m Minor allele frequency (MAF) cutoff (default 0.05)"
   echo -e "\t-g Missing genotype rate cutoff (default 0.05)"
   echo -e "\t-h Hardy-Weinberg equilibrium p-value cutoff (default 0.001)"
   echo -e "\t-p Independent pairwise LD pruning settings (default \"1000 50 0.3\"). Numbers indicate window, step size, and r2 value for indep-pairwise pruning."
   echo -e "\t-t Threads to use for bcftools decompression (default 0; plink uses all available CPUs.)"
   exit 1 # Exit script after printing help
}

# set default values
maf_cutoff=0.05
missing_genotype_rate=0.05
hwe_pvalue=0.001
indep_pairwise_ld_pruning_values="1000 50 0.3"
sample_list_action="exclude"
threads=0

# get options listed above
# defaults overridden if necessary
while getopts "i:s:a:m:g:h:p:t:" opt
do
   case "$opt" in
      i ) input_prefix="$OPTARG" ;;
      s ) sample_list="$OPTARG" ;;
      a ) sample_list_action="$OPTARG" ;;
      m ) maf_cutoff="$OPTARG" ;;
      g ) missing_genotype_rate="$OPTARG" ;;
      h ) hwe_pvalue="$OPTARG" ;;
      p ) indep_pairwise_ld_pruning_values="$OPTARG" ;;
      t ) threads="$OPTARG" ;;
      ? ) helpFunction ;; # Print helpFunction in case parameter is non-existent
   esac
done

# Print helpFunction in case required parameters are empty
if [ -z "$input_prefix" ]
then
   echo "The required parameter (-i input_prefix) is empty.";
   helpFunction
fi

# exit if sample list action is set to neither exclude nor include
if [[ $sample_list_action != "exclude" && $sample_list_action != "include" ]]
then
  echo "Sample list action not specified as 'include' or 'exclude', exiting."
  exit
fi

# different chained commands depending on whether excluded sample list provided or not
# in both cases, split all multiallelic variants, filter based on MAF/genotype/HWE, and then LD prune with indep-pairwise argument
# note that variant IDs over 16000 characters should be rewritten as chr_pos_alt to escape plink errors with bcftools annotate before plink
# bcftools annotate -I '%CHROM\_%POS\_%ALT' $input_prefix.vcf.gz -O z -o $input_prefix"_"relabeled.vcf.gz --threads $threads
if [ -z "$sample_list" ]
then
    	bcftools norm -m -both $input_prefix.vcf.gz -O z -o $input_prefix"_"multsplit.vcf.gz --threads $threads
        plink2 --vcf $input_prefix"_"multsplit.vcf.gz --split-par hg38 --autosome-par --output-chr chrMT --maf $maf_cutoff --geno $missing_genotype_rate --hwe $hwe_pvalue --vcf-half-call m --rm-dup force-first --export vcf bgz --out $input_prefix"_"multsplit_filtered
        plink2 --vcf $input_prefix"_"multsplit_filtered.vcf.gz --output-chr chrMT --indep-pairwise $indep_pairwise_ld_pruning_values --out $input_prefix"_"multsplit_filtered_pruned
	# extract variants listed as pruned 
	plink2 --extract $input_prefix"_"multsplit_filtered_pruned.prune.in --output-chr chrMT --vcf $input_prefix"_"multsplit_filtered.vcf.gz --export vcf bgz --out $input_prefix"_"multsplit_filtered_pruned
else
	# first check sample list action
	if [[ $sample_list_action == "exclude" ]]
	then
		# first exclude samples with bcftools view -S ^exclude_list_filename
		bcftools view -S "^"$sample_list $input_prefix.vcf.gz -O z -o $input_prefix"_"sample_excluded.vcf.gz --threads $threads
		bcftools norm -m -both $input_prefix"_"sample_excluded.vcf.gz -O z -o $input_prefix"_"sample_excluded_multsplit.vcf.gz --threads $threads
		plink2 --vcf $input_prefix"_"sample_excluded_multsplit.vcf.gz --split-par hg38 --autosome-par --output-chr chrMT --maf $maf_cutoff --geno $missing_genotype_rate --hwe $hwe_pvalue --vcf-half-call m --rm-dup force-first --export vcf bgz --out $input_prefix"_"sample_excluded_multsplit_filtered
		plink2 --vcf $input_prefix"_"sample_excluded_multsplit_filtered.vcf.gz --output-chr chrMT --indep-pairwise $indep_pairwise_ld_pruning_values --out $input_prefix"_"sample_excluded_multsplit_filtered_pruned
		plink2 --extract $input_prefix"_"sample_excluded_multsplit_filtered_pruned.prune.in --vcf $input_prefix"_"sample_excluded_multsplit_filtered.vcf.gz --output-chr chrMT --export vcf bgz --out $input_prefix"_"sample_excluded_multsplit_filtered_pruned
	else
		bcftools view -S $sample_list $input_prefix.vcf.gz -O z -o $input_prefix"_"sample_included.vcf.gz --threads $threads
		bcftools norm -m -both $input_prefix"_"sample_included.vcf.gz -O z -o $input_prefix"_"sample_included_multsplit.vcf.gz --threads $threads
		plink2 --vcf $input_prefix"_"sample_included_multsplit.vcf.gz --split-par hg38 --autosome-par --output-chr chrMT --maf $maf_cutoff --geno $missing_genotype_rate --hwe $hwe_pvalue --vcf-half-call m --rm-dup force-first --export vcf bgz --out $input_prefix"_"sample_included_multsplit_filtered
		plink2 --vcf $input_prefix"_"sample_included_multsplit_filtered.vcf.gz --output-chr chrMT --indep-pairwise $indep_pairwise_ld_pruning_values --out $input_prefix"_"sample_included_multsplit_filtered_pruned
		plink2 --extract $input_prefix"_"sample_included_multsplit_filtered_pruned.prune.in --vcf $input_prefix"_"sample_included_multsplit_filtered.vcf.gz --output-chr chrMT --export vcf bgz --out $input_prefix"_"sample_included_multsplit_filtered_pruned
	fi
fi
