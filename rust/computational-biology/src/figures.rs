//! Computational Biology's own live figures (ax_engine::figkit), drawn for the visitor's inputs from what the model
//! computes: sequences, alignment tables, trees, crosses, population trajectories, kinetics and fits
//! (docs/FIGURE_PARITY.md). The figure families live in figures_seq.rs, figures_align.rs, figures_trees.rs,
//! figures_pop.rs and figures_kin.rs.

use crate::figures_align::{dp_table, Mode};
use crate::{figures_kin as kin, figures_pop as pop, figures_seq as seq, figures_trees as tr};
use ax_engine::figkit::Solved;
use ax_engine::figures::Figure;

/// The figures for one solved request (the dispatcher ax_engine calls: Pack::figures).
pub fn live(s: &Solved) -> Vec<Figure> {
    let v: Vec<Option<Figure>> = match s.key {
        // sequences
        "bio_aromaticity" => vec![seq::aromatic_example(s), seq::aromatic_sequence(s)],
        "bio_gc_content" => vec![seq::gc_profile(s), seq::gc_sequence(s)],
        "bio_gravy" => vec![seq::gravy_profile(s)],
        "bio_instability_index" => vec![seq::instability_profile(s)],
        "bio_kmer_counts" => vec![seq::kmer_sequence(s)],
        "bio_longest_orf" => vec![seq::longest_orf(s)],
        "bio_reverse_complement" => vec![seq::duplex(s)],
        "bio_transcription" => vec![seq::transcription(s)],
        "bio_translation" => vec![seq::codons(s)],
        "bio_hamming_distance" => vec![seq::hamming(s)],
        "bio_restriction_sites" => vec![seq::digest(s), seq::gel(s)],
        "bio_isoelectric_point" => vec![seq::titration(s)],
        "bio_protein_extinction" => vec![seq::extinction_heat(s)],
        "bio_protein_mass" => vec![seq::mass_bars(s)],
        "bio_amino_acid_composition" => vec![seq::composition_bars(s)],
        "bio_secondary_structure_fraction" => vec![seq::structure_bars(s)],
        "bio_nucleotide_diversity" => vec![seq::nd_sample(s)],
        "bio_watterson_theta" => vec![seq::watterson_sample(s)],
        "bio_tajima_d" => vec![seq::tajima_bars(s), seq::tajima_sample(s)],
        // alignment tables
        "bio_edit_distance" => vec![dp_table(s, Mode::Edit, "Edit distance table with one cheapest series of edits",
                                             "Edit distance (unit costs)")],
        "bio_global_alignment_score" => vec![dp_table(s, Mode::GlobalLinear, "Global alignment table with an optimal path",
                                                      "Needleman-Wunsch global alignment score")],
        "bio_global_alignment_affine" => vec![dp_table(s, Mode::GlobalAffine, "Global alignment score of every pair of prefixes",
                                                       "Global alignment score with affine gaps")],
        "bio_protein_global_blosum62" => vec![dp_table(s, Mode::GlobalAffine, "Global alignment score of every pair of prefixes",
                                                       "Global protein alignment score (BLOSUM62, affine gaps)")],
        "bio_semiglobal_alignment_score" => vec![dp_table(s, Mode::SemiGlobal,
                                                          "Semi-global alignment score of every pair of prefixes",
                                                          "Semi-global (end-gap-free) alignment score")],
        "bio_local_alignment_score" => vec![dp_table(s, Mode::Local, "Local alignment: best score within growing prefixes",
                                                     "Local alignment score (linear gaps)")],
        "bio_local_alignment_affine" => vec![dp_table(s, Mode::Local, "Local alignment: best score within growing prefixes",
                                                      "Local alignment score with affine gaps")],
        "bio_protein_local_blosum62" => vec![dp_table(s, Mode::Local, "Local alignment: best score within growing prefixes",
                                                      "Local protein alignment score (BLOSUM62, affine gaps)")],
        // trees, pedigrees, networks
        "bio_faith_pd" => vec![tr::faith_bars(s), tr::faith_tree(s)],
        "bio_fitch_parsimony" => vec![tr::fitch_tree(s)],
        "bio_neighbor_joining" => vec![tr::nj_matrix(s), tr::nj_tree(s)],
        "bio_upgma" => vec![tr::upgma_matrix(s), tr::upgma_tree(s)],
        "bio_robinson_foulds" => vec![tr::rf_trees(s)],
        "bio_pedigree_inbreeding" => vec![tr::pedigree_bars(s), tr::pedigree(s)],
        "bio_coalescent_times" => vec![tr::genealogy(s)],
        "bio_flux_balance_analysis" => vec![tr::flux_network(s)],
        // populations, genetics, diversity
        "bio_ace" => vec![pop::ace_bars(s)],
        "bio_chao1" => vec![pop::chao1_bars(s)],
        "bio_beta_diversity" => vec![pop::beta_bars(s)],
        "bio_ne_fluctuating_size" => vec![pop::ne_fluct(s)],
        "bio_ne_sex_ratio" => vec![pop::ne_sex_field(s)],
        "bio_cross_genotype_probability" => vec![pop::genotype_bars(s), pop::genotype_punnett(s)],
        "bio_cross_phenotype_ratio" => vec![pop::phenotype_punnett(s)],
        "bio_drift_heterozygosity" => vec![pop::drift_heat(s)],
        "bio_fst" => vec![pop::fst_heat(s)],
        "bio_heterozygote_advantage" => vec![pop::het_adv_heat(s)],
        "bio_haploid_selection" => vec![pop::haploid_traj(s)],
        "bio_recessive_lethal" => vec![pop::lethal_traj(s)],
        "bio_hardy_weinberg_test" => vec![pop::hw_bars(s), pop::hw_triangle(s)],
        "bio_hardy_weinberg_exact_test" => vec![pop::hw_exact_bars(s)],
        "bio_leslie_projection" => vec![pop::leslie_population(s), pop::leslie_trajectory(s)],
        "bio_linkage_disequilibrium" => vec![pop::ld_table(s)],
        // kinetics and fits
        "bio_michaelis_menten_fit" => vec![kin::mm_fit(s), kin::mm_direct(s)],
        "bio_competitive_inhibition_fit" => vec![kin::competitive_fit(s), kin::competitive_linear(s)],
        "bio_eadie_hofstee" => vec![kin::eadie(s)],
        "bio_hanes_woolf" => vec![kin::hanes(s)],
        "bio_lineweaver_burk" => vec![kin::lineweaver(s)],
        "bio_hill_fit" => vec![kin::hill_fit(s)],
        "bio_hill_plot" => vec![kin::hill_plot(s)],
        "bio_schumacher_eschmeyer" => vec![kin::schumacher(s)],
        "bio_species_area" => vec![kin::species_area(s)],
        "bio_michaelis_menten_progress" => vec![kin::mm_progress(s), kin::mm_progress_family(s)],
        "bio_chemostat_steady_state" => vec![kin::chemostat(s)],
        _ => vec![],
    };
    v.into_iter().flatten().collect()
}
