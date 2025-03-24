#!/usr/bin/bash
# script to generate initial CSV from bcftools query for genetic data and genetic map generation in python
# script depends on bcftools
module load bcftools
helpFunction()
{
   echo ""
   echo "Usage: $0 -v variant_type -c chromosome -i input.vcf(.gz) -o output.csv"
   echo -e "\t-v Variant type - structural variants (SV) or single nucleotide variants (SNV)"
   echo -e "\t-c Chromosome to subset (optional)"
   echo -e "\t-i Input VCF file"
   echo -e "\t-o Output CSV file for genetic data/map generation with helper Python script"
   exit 1 # Exit script after printing help
}

# get options listed above
while getopts "v:c:i:o:" opt
do
   case "$opt" in
      v ) variant_type="$OPTARG" ;;
      c ) chromosome="$OPTARG" ;;
      i ) input="$OPTARG" ;;
      o ) output="$OPTARG" ;;
      ? ) helpFunction ;; # Print helpFunction in case parameter is non-existent
   esac
done

# Print helpFunction in case required parameters are empty
if [ -z "$variant_type" ] || [ -z "$input" ] || [ -z "$output" ]
then
   echo "Some or all of the required parameters are empty.";
   helpFunction
fi

# Begin script in case all parameters are correct
echo "Variant type: ${variant_type}"
echo "Chromosome: ${chromosome}"
echo "Input: ${input}"
echo "Output: ${output}"

# Call bcftools in different ways depending on parameters
if [ ${variant_type} == "SNV" ]; then
	# SNV conversion
	# use bcftools norm -m -both to collapse all multiallelic sites
	# use bcftools query -HH to make initial preprocessed CSV with header
	if [ -z "$chromosome" ]; then
		bcftools norm -m -both ${input} | bcftools query -HH -f "%ID,%CHROM,%POS,%END,%REF,%ALT[,%GT]" | sed 's/\#//g' | sed 's/:GT//g' > ${output}
	else
		bcftools norm -m -both -r ${chromosome} ${input} | bcftools query -HH -f "%ID,%CHROM,%POS,%END,%REF,%ALT[,%GT]" | sed 's/\#//g' | sed 's/:GT//g' > ${output}
	fi
fi
if [ ${variant_type} == "SV" ]; then
	# SV conversion
	# use bcftools norm -m -both to collapse all multiallelic sites
	if [ -z "$chromosome" ]; then
		bcftools norm -m -both ${input} | bcftools query -HH -f "%ID,%CHROM,%POS,%END,%INFO/SVTYPE,%INFO/SVLEN,%REF,%ALT[,%GT]" | sed 's/\#//g' | sed 's/:GT//g' > ${output}
	else
		bcftools norm -m -both -r ${chromosome} ${input} | bcftools query -HH -f "%ID,%CHROM,%POS,%END,%INFO/SVTYPE,%INFO/SVLEN,%REF,%ALT[,%GT]" | sed 's/\#//g' | sed 's/:GT//g' > ${output}
	fi
fi
