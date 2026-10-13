# Data dictionary

Columns of every archived derived table. Values are as written by the scripts under
`pipeline/scripts/`; the `.prov` file next to each `results/figures/*.csv` records the
UTC timestamp and the source files of that particular build.

## `results/dedup/sites_tierA.tsv`
union variant sites for the Tier A matrix, with the column index in the genotype matrix
- rows: 12926
- columns: `chrom`, `pos`, `column`

## `results/dedup/sites_tierB.tsv`
union variant sites for the Tier B matrix
- rows: 67711
- columns: `chrom`, `pos`, `column`

## `results/dedup/tierA_clusters.conservative.tsv`
- rows: 209
- columns: `cluster_id`, `size`, `representative`, `diameter`, `members`

## `results/dedup/tierA_clusters.full.tsv`
- rows: 206
- columns: `cluster_id`, `size`, `representative`, `diameter`, `members`

## `results/dedup/tierA_clusters.off.tsv`
- rows: 175
- columns: `cluster_id`, `size`, `representative`, `diameter`, `members`

## `results/dedup/tierA_clusters.threestate.tsv`
- rows: 205
- columns: `cluster_id`, `size`, `representative`, `diameter`, `members`

## `results/dedup/tierB_backbone.tsv`
- rows: 621
- columns: `accession`

## `results/dedup/tierB_clusters.tsv`
- rows: 621
- columns: `cluster_id`, `size`, `representative`, `diameter`, `members`

## `results/figures/fig1a_called_frac.csv`
- rows: 3697
- columns: `tier`, `sample`, `called_frac`

## `results/figures/fig1b_distance_dist.csv`
- rows: 12
- columns: `percentile`, `raw_distance`

## `results/figures/fig1b_distance_hist.csv`
- rows: 60
- columns: `bin_lo`, `bin_hi`, `n_pairs`

## `results/figures/fig1b_same_strain.csv`
- rows: 5
- columns: `strain`, `acc_a`, `acc_b`, `raw_distance_tierA`, `raw_distance_tierB`

## `results/figures/fig1b_scalars.csv`
- rows: 6
- columns: `statistic`, `value`

## `results/figures/fig1c_self_blocks.csv`
- rows: 216
- columns: `block_index`, `length_bp`, `replicon`, `start`, `end`

## `results/figures/fig1d_cluster_internal.csv`
- rows: 20
- columns: `distance_snp`, `n_pairs`

## `results/figures/fig1d_cluster_years.csv`
- rows: 7
- columns: `statistic`, `value`

## `results/figures/fig1d_qinghai.csv`
- rows: 419
- columns: `sample`, `year`, `country`

## `results/figures/fig2a_metric_flags.csv`
- rows: 765
- columns: `accession`, `region`, `metric`, `value`, `flagged`

## `results/figures/fig2b_regression.csv`
- rows: 255
- columns: `accession`, `distance_to_CO92`, `called_frac`, `resid_z`, `flag`, `anchor_label`

## `results/figures/fig2c_dnadiff.csv`
- rows: 256
- columns: `accession`, `ref_pct`, `qry_pct`, `plasmid_submitted`, `strain_label`

## `results/figures/fig2d_checkm.csv`
- rows: 255
- columns: `accession`, `completeness`, `anchor_label`

## `results/figures/fig3a_block_lengths.csv`
- rows: 8227
- columns: `accession`, `replicon`, `start`, `length_bp`

## `results/figures/fig3b_level_frac.csv`
- rows: 4
- columns: `level`, `n_total`, `n_with_block`, `frac`

## `results/figures/fig3b_recurrence.csv`
- rows: 161
- columns: `replicon`, `start`, `length_bp`, `n_genomes`

## `results/figures/fig3c_pgm_align.csv`
- rows: 33
- columns: `strain`, `accession`, `start`, `end`, `length_bp`

## `results/figures/fig3c_pgm_genes.csv`
- rows: 15
- columns: `gene`, `start`, `end`, `group`

## `results/figures/fig4a_mechanism.csv`
- rows: 4
- columns: `site_class`, `handling`, `affects_distance`

## `results/figures/fig4b_co92_pair.csv`
- rows: 2
- columns: `accession`, `label`, `raw_distance`, `replicons_absent`, `big_blocks`

## `results/figures/fig4b_co92_tracks.csv`
- rows: 1
- columns: `label`, `accession`, `replicon`, `start`, `end`, `length_bp`

## `results/figures/fig4c_constraint.csv`
- rows: 2
- columns: `tier`, `n_before_constraint`, `n_after`, `mixed_profile_clusters_before`, `mixed_profile_clusters_after`

## `results/figures/fig4c_constraints.csv`
- rows: 4
- columns: `constraint`, `retained`, `clusters`, `multi_member_clusters`

## `results/figures/fig4d_sweep.csv`
- rows: 11
- columns: `threshold`, `retained`, `max_diameter`, `max_size`, `dup_pairs_merged`

## `results/figures/fig5a_region_polarity.csv`
- rows: 13
- columns: `region_id`, `replicon`, `start`, `end`, `length_bp`, `n_genomes`, `og1_callable`, `og2_callable`, `polarity`, `n_loss_votes`, `highlight`

## `results/figures/fig5b_background.csv`
- rows: 1000
- columns: `window_id`, `delta`

## `results/figures/fig5b_block_delta.csv`
- rows: 18
- columns: `region_id`, `delta`, `percentile`, `within_background`

## `results/figures/fig_m1_edge_labels.csv`
- rows: 1
- columns: `x;y;label;hjust`

## `results/figures/fig_m1_edges.csv`
- rows: 27
- columns: `edge_id;seq;x;y;arrow_end;dashed`

## `results/figures/fig_m1_nodes.csv`
- rows: 8
- columns: `node_id;xmin;xmax;ymin;ymax;label;sublabel;style`

## `results/figures/fig_m1_sites.csv`
- rows: 32
- columns: `site_index`, `genome`, `state`

## `results/qc/archive_bias.tsv`
- rows: 20
- columns: `group`, `n`, `indicator`, `n_positive`, `pct_positive`, `baseline_n`, `baseline_n_positive`, `baseline_pct`, `fisher_p`

## `results/qc/assembly_metadata.tsv`
NCBI metadata for every assembly considered, after GCA/GCF mirror merge
- rows: 4585
- columns: `accession`, `assembly_level`, `assembly_status`, `total_sequence_length`, `number_of_contigs`, `contig_n50`, `gc_percent`, `ambiguous_fraction`, `checkm_completeness`, `checkm_contamination`, `checkm_version`, `source_database`, `paired_accession`, `reference_genome`, `strain`, `collection_date`, `geo_loc_name`, `host`, `isolation_source`, `bioproject`, `release_date`

## `results/qc/c4_decomposition.tsv`
dnadiff unaligned-reference decomposition into the three mutually exclusive categories of Methods 2.8
- rows: 256
- columns: `accession`, `unaligned_mcoords`, `unaligned_report`, `diff`, `del_replicon`, `del_block`, `residual`, `absent_replicons`, `n_carried_replicons`

## `results/qc/cluster349_internal.tsv`
- rows: 16
- columns: `statistic`, `value`

## `results/qc/contamination_scale.txt`
- rows: 3

## `results/qc/dedup_mixed_classification.tsv`
classification of the multi-member clusters formed without a constraint
- rows: 36
- columns: `cluster_id`, `size`, `representative`, `category`, `members`

## `results/qc/dedup_modes.txt`
- rows: 45

## `results/qc/dedup_modes_summary.txt`
cluster counts under the four deletion-profile constraints
- rows: 11

## `results/qc/deletion_profile_all.tsv`
one row per profiled assembly: which reference replicons are absent and which >=10 kb chromosomal blocks are uncalled (5 kb bins)
- rows: 3441
- columns: `accession`, `replicons_absent`, `big_blocks_5kb`, `n_big`, `n_absent`

## `results/qc/excluded_with_reason.tsv`
- rows: 210
- columns: `accession`, `tier`, `assembly_level`, `total_bp`, `contigs`, `contig_n50`, `checkm_contamination`, `checkm_completeness_COVARIATE_ONLY`, `dnadiff_aligned_pct`, `anchor`, `collection_date`, `geo_loc_name`, `host`, `strain`, `notes`

## `results/qc/lab_strain_inventory.tsv`
genomes whose strain field matches a laboratory naming pattern, with their absent replicons and blocks
- rows: 25
- columns: `query`, `accession`, `strain_field`, `geo_loc_name`, `replicons_absent`, `big_blocks_5kb`

## `results/qc/lineage_stratified_qc.tsv`
- rows: 256
- columns: `accession`, `group`, `group_size`, `called_frac`, `lineage_median_cf`, `within_lineage_z`, `flag`, `deletion_verdict`, `strain`, `geo_loc_name`

## `results/qc/metadata_completeness.tsv`
- rows: 8
- columns: `dataset`, `field`, `nonnull`, `total`, `fraction`

## `results/qc/outgroup_concordance.tsv`
per-block callable fractions of the two outgroups
- rows: 18
- columns: `replicon`, `start_bin`, `length_bin`, `og1`, `og2`, `delta`, `background_percentile`, `polarity`, `n_loss_votes`, `verdict`

## `results/qc/pseudogene_entanglement.tsv`
curated CO92 pseudogenes against the automated annotation (Table S6)
- rows: 217
- columns: `replicon`, `ncbi_start`, `ncbi_end`, `ncbi_span`, `ncbi_gene`, `ncbi_product`, `bakta_span`, `bakta_frac_of_ncbi`, `bakta_gene`, `category`

## `results/qc/quality_metric_geography.tsv`
- rows: 29
- columns: `metric`, `accession`, `value`, `unit`, `region`, `geo_loc_name`, `strain`

## `results/qc/recurrent_block_polarity.singleOG.tsv`
- rows: 18
- columns: `replicon`, `start_bin`, `length_bin`, `n_genomes_absent`, `pct_genomes_absent`, `outgroup_called_frac`, `polarity`, `n_genes`, `gene_names`, `products`, `note`

## `results/qc/recurrent_block_polarity.tsv`
one row per recurrent >=10 kb block with its outgroup callable fractions and assigned polarity
- rows: 18
- columns: `replicon`, `start_bin`, `length_bin`, `n_genomes_absent`, `pct_genomes_absent`, `outgroup_called_frac`, `outgroup_fracs_all`, `n_outgroups_informative`, `n_loss_votes`, `n_gain_votes`, `polarity`, `n_genes`, `gene_names`, `products`, `note`

## `results/qc/reference_bias.tsv`
per-genome reference-relative metrics: distance to CO92, called fraction, aligned bases, CheckM completeness, residual z
- rows: 255
- columns: `accession`, `dist_to_reference`, `called_frac`, `dnadiff_aligned_pct`, `checkm_completeness`, `strain`, `geo_loc_name`, `predicted_called_frac`, `residual`, `residual_z`, `flag`

## `results/qc/replicon_migration_P258.tsv`
- rows: 258
- columns: `accession`, `n_seq`, `n_plasmid_seq`, `old_replicons_absent`, `pCD1_state`, `pPCP1_state`, `pMT1_state`

## `results/qc/replicon_three_state.tsv`
- rows: 3443
- columns: `accession`, `n_seq`, `n_plasmid_seq`, `total_bp`, `pCD1_state`, `pPCP1_state`, `pMT1_state`, `old_replicons_absent`, `evidence`

## `results/qc/replicon_three_state.txt`
- rows: 0

## `results/qc/replicon_three_state_seqreport.tsv`
per-genome plasmid state from the submitted sequence inventory: CARRIED / LOST / UNRESOLVED
- rows: 258
- columns: `accession`, `n_seq`, `n_plasmid_seq`, `total_bp`, `pCD1_state`, `pPCP1_state`, `pMT1_state`, `old_replicons_absent`, `evidence`

## `results/qc/replicon_three_state_seqreport.txt`
- rows: 30

## `results/qc/s3_pgm_only.tsv`
the pgm-only directional test and its refined sensitivity analysis
- rows: 2
- columns: `test`, `arm1_labnamed`, `arm1_n`, `arm2_labnamed`, `arm2_n`, `fisher_p`, `odds_ratio`

## `results/qc/table2_regions.tsv`
the 13 merged recurrent regions (Table 2): coordinates, carrier counts, polarity, background percentile
- rows: 12
- columns: `replicon`, `start`, `end`, `length`, `n_blocks_merged`, `n_genomes_absent`, `pct_genomes_absent`, `outgroup_fracs`, `polarity`, `polarities_in_merge`, `n_loss_votes`, `background_pctile`, `support`, `gene_names`

## `results/qc/tierA_candidates.tsv`
- rows: 256
- columns: `accession`, `tier`, `assembly_level`, `total_bp`, `contigs`, `contig_n50`, `checkm_contamination`, `checkm_completeness_COVARIATE_ONLY`, `dnadiff_aligned_pct`, `anchor`, `collection_date`, `geo_loc_name`, `host`, `strain`, `notes`

## `results/qc/tierA_refined.tsv`
the 256-genome analysis set (A256) with its per-genome QC covariates
- rows: 256
- columns: `accession`, `tier`, `assembly_level`, `total_bp`, `contigs`, `contig_n50`, `checkm_contamination`, `checkm_completeness_COVARIATE_ONLY`, `dnadiff_aligned_pct`, `anchor`, `collection_date`, `geo_loc_name`, `host`, `strain`, `notes`, `refine_note`

## `results/qc/tierB_candidates.tsv`
the 3,441-genome Tier B set
- rows: 3441
- columns: `accession`, `tier`, `assembly_level`, `total_bp`, `contigs`, `contig_n50`, `checkm_contamination`, `checkm_completeness_COVARIATE_ONLY`, `dnadiff_aligned_pct`, `anchor`, `collection_date`, `geo_loc_name`, `host`, `strain`, `notes`

## `results/qc/tier_report.tsv`
set reduction at each stage of Methods 2.1 with the reason for each exclusion
- rows: 3651
- columns: `accession`, `tier`, `assembly_level`, `total_bp`, `contigs`, `contig_n50`, `checkm_contamination`, `checkm_completeness_COVARIATE_ONLY`, `dnadiff_aligned_pct`, `anchor`, `collection_date`, `geo_loc_name`, `host`, `strain`, `notes`

## `results/qc/uncalled_blocks.tsv`
the same blocks restricted to the closed genomes used for deletion-profile analysis
- rows: 8227
- columns: `accession`, `replicon`, `start`, `end`, `length`

## `results/qc/uncalled_blocks_all.tsv`
every >=10 kb uncalled block in every profiled assembly (5 kb-binned start)
- rows: 3682
- columns: `accession`, `replicon`, `start`, `end`, `length`

## `results/qc/uncalled_profile.tsv`
- rows: 50
- columns: `accession`, `total_bp`, `uncalled_bp`, `called_frac`, `n_blocks_ge500bp`, `n_IS_like_500_3000bp`, `bp_in_IS_like`, `n_blocks_ge10kb`, `bp_in_blocks_ge10kb`, `n_pgm_sized_80_130kb`, `pgm_sized_coords`, `replicons_absent`, `largest_block_bp`, `verdict`

## `results/qc/variant_freq_strata.tsv`
loss frequencies for pgm, pCD1, pPCP1 and pMT1 in three nested sets
- rows: 3
- columns: `set`, `n`, `pgm_pct`, `pPCP1_pct`, `pCD1_pct`, `pMT1_pct`
