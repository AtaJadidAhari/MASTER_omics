# VEP Sequence Ontology severity ordering (most severe first)
VEP_SORanking = [
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

    # MODERATE
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
    "sequence_variant",
]

# Severity rank: lower = more severe (same as seq_along in R)
severity_rank = {term: i + 1 for i, term in enumerate(VEP_SORanking)}

# Variant category helpers (same semantics as R vectors)
splice_region_variants = {
    "splice_region_variant",
    "splice_donor_region_variant",
    "splice_polypyrimidine_tract_variant",
    "splice_donor_5th_base_variant",
}

splice_site_variants = {
    "splice_donor_variant",
    "splice_acceptor_variant",
}

stop_variants = {
    "stop_gained",
    "stop_lost",   # stop_retained_variant intentionally excluded
}

inframe_variants = {
    "inframe_insertion",
    "inframe_deletion",
}

up_down_stream_vars = {
    "downstream_gene_variant",
    "upstream_gene_variant",
}
