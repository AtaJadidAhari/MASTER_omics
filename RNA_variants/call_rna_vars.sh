#!/bin/bash

input_bam = $1
input_bai = $2

ref = 

gatk --java-options -Djava.io.tmpdir=${tmpdir} HaplotypeCaller \
  -I $input_bam \
  -R $ref \
  -L my_mapped_genes.bed \
  --dont-use-soft-clipped-bases \
  -stand-call-conf 20.0 \
  --dbsnp "${tmp_vcf}" \
  --output-mode EMIT_ALL_CONFIDENT_SITES \
  -ERC GVCF \
  $hcArgs \
  -O $output_gVCF 2>&1 | tee -a $log