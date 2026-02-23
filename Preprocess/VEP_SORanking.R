VEP_SORanking <- c(
  # HIGH
  "transcript_ablation",
  "splice_acceptor_variant", 
  "splice_donor_variant",
  "stop_gained",
  "frameshift_variant",
  "stop_lost",
  "start_lost",
  "transcript_amplification",
  "feature_elongation",
  "feature_truncation",
  
  # MODERAT
  "inframe_insertion",
  "inframe_deletion",
  "missense_variant",
  "protein_altering_variant",
  
  # LOW
  "splice_donor_5th_base_variant",
  "splice_region_variant",
  "splice_donor_region_variant",
  "splice_polypyrimidine_tract_variant",
  "incomplete_terminal_codon_variant",
  "start_retained_variant",
  "stop_retained_variant",
  "synonymous_variant",
  
  # MODIFIER
  "coding_sequence_variant",
  "mature_miRNA_variant",
  "5_prime_UTR_variant",
  "3_prime_UTR_variant",
  "non_coding_transcript_exon_variant",
  "intron_variant",
  "NMD_transcript_variant",
  "non_coding_transcript_variant",
  "coding_transcript_variant",
  "upstream_gene_variant",
  "downstream_gene_variant",
  "TFBS_ablation",
  "TFBS_amplification",
  "TF_binding_site_variant",
  "regulatory_region_ablation",
  "regulatory_region_amplification",
  "regulatory_region_variant",
  "intergenic_variant",
  "sequence_variant"
)

severity_rank <- setNames(seq_along(VEP_SORanking), VEP_SORanking)


splice_region_variants <- c("splice_region_variant", "splice_donor_region_variant", "splice_polypyrimidine_tract_variant", "splice_donor_5th_base_variant")
splice_site_variants <- c("splice_donor_variant", "splice_acceptor_variant")
stop_variants <- c("stop_gained", "stop_lost") # stop_retained_variant
inframe_variants <- c("inframe_insertion", "inframe_deletion")
up_down_stream_vars <- c("downstream_gene_variant", "upstream_gene_variant")


