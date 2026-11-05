//! computational_biology: procedures translated from the package's Python by tools/py2rust.py (do not edit; regenerate).
//! Source: models/computational_biology/capabilities/procedures.py. Each keeps Python's number rules through
//! ax_engine::py and is checked against the package by tools/rust_accuracy.py before the build is listed.
#![allow(clippy::all, non_snake_case, unused_parens, unused_braces, unused_mut, unused_variables, unreachable_code)]

use ax_engine::py::{self, Py, PyErr, R};
use ax_engine::Val;
use serde_json::Value;
use std::collections::HashMap;

/// The procedures this crate implements (tools/export_rust_pack.py reads this list).
pub const NATIVE: &[&str] = &[
    "bio_ace",
    "bio_amino_acid_composition",
    "bio_aromaticity",
    "bio_beta_diversity",
    "bio_brillouin_diversity",
    "bio_chao1",
    "bio_chemostat_steady_state",
    "bio_coalescent_times",
    "bio_competitive_inhibition_fit",
    "bio_cross_genotype_probability",
    "bio_cross_phenotype_ratio",
    "bio_dominance_indices",
    "bio_drift_heterozygosity",
    "bio_eadie_hofstee",
    "bio_edit_distance",
    "bio_faith_pd",
    "bio_fitch_parsimony",
    "bio_flux_balance_analysis",
    "bio_fst",
    "bio_gc_content",
    "bio_global_alignment_affine",
    "bio_global_alignment_score",
    "bio_goods_coverage",
    "bio_gravy",
    "bio_hamming_distance",
    "bio_hanes_woolf",
    "bio_haploid_selection",
    "bio_hardy_weinberg_exact_test",
    "bio_hardy_weinberg_test",
    "bio_heterozygote_advantage",
    "bio_hill_fit",
    "bio_hill_plot",
    "bio_instability_index",
    "bio_isoelectric_point",
    "bio_jukes_cantor_distance",
    "bio_kimura_distance",
    "bio_kmer_counts",
    "bio_leslie_growth_rate",
    "bio_leslie_projection",
    "bio_lineweaver_burk",
    "bio_linkage_disequilibrium",
    "bio_local_alignment_affine",
    "bio_local_alignment_score",
    "bio_longest_orf",
    "bio_map_distance",
    "bio_michaelis_menten_fit",
    "bio_michaelis_menten_progress",
    "bio_michaelis_menten_time",
    "bio_moran_fixation",
    "bio_mutation_selection_balance",
    "bio_ne_fluctuating_size",
    "bio_ne_sex_ratio",
    "bio_neighbor_joining",
    "bio_nucleic_acid_mass",
    "bio_nucleotide_diversity",
    "bio_pedigree_inbreeding",
    "bio_petersen_estimate",
    "bio_protein_extinction",
    "bio_protein_global_blosum62",
    "bio_protein_local_blosum62",
    "bio_protein_mass",
    "bio_protein_net_charge",
    "bio_rarefaction",
    "bio_recessive_lethal",
    "bio_recombination_fraction",
    "bio_restriction_sites",
    "bio_reverse_complement",
    "bio_richness_indices",
    "bio_robinson_foulds",
    "bio_schumacher_eschmeyer",
    "bio_secondary_structure_fraction",
    "bio_semiglobal_alignment_score",
    "bio_shannon_diversity",
    "bio_simpson_diversity",
    "bio_species_area",
    "bio_tajima_d",
    "bio_tm_nearest_neighbor",
    "bio_tm_wallace",
    "bio_transcription",
    "bio_translation",
    "bio_tree_count",
    "bio_upgma",
    "bio_watterson_theta",
];

type Out = Result<(Vec<(String, Value)>, f64), String>;

/// The run's result as the engine's answer: a dict of outputs and the residual of the package's own check; a
/// ValueError is the refusal's reason; any other error is an internal error (None).
fn finish(r: R, residual: impl Fn(&Py, &Py) -> f64, v: &HashMap<String, Val>) -> Option<Out> {
    let o = match r {
        Ok(o) => o,
        Err(PyErr::Value(m)) => return Some(Err(m)),
        Err(e) => {
            if std::env::var_os("AX_DEBUG").is_some() {        // the error behind an internal error, for diagnosis
                eprintln!("internal error: {e:?}");
            }
            return None;
        }
    };
    let Py::Dict(d) = &o else { return None };
    let outs: Vec<(String, Value)> = d.borrow().iter().map(|(k, x)| (py::to_str(k), py::to_json(x))).collect();
    let i = Py::dict(v.iter().map(|(k, x)| (Py::str(k), py::from_val(x))).collect());
    let res = residual(&i, &o);
    if std::env::var_os("AX_DEBUG").is_some() && !(res <= 1e-9) {
        eprintln!("check failed: residual {res}");
    }
    Some(Ok((outs, res)))
}

/// The engine's entry point (ax_engine::pack::NativeFn): None for a key that is not one of these procedures.
pub fn native(key: &str, v: &HashMap<String, Val>) -> Option<Out> {
    match key {
        "bio_ace" => {
            let r = f_ace(match v.get("counts") { Some(x) => py::from_val(x), None => { return None } }, match v.get("rare_threshold") { Some(x) => py::from_val(x), None => { return None } }, (|| -> R { Ok(Py::Bool(false)) })().ok()?);
            finish(r, |i: &Py, o: &Py| -> f64 { match s_36(i.clone(), o.clone()).and_then(|r| r.f()) { Ok(x) => x, Err(e) => { if std::env::var_os("AX_DEBUG").is_some() { eprintln!("check error: {e:?}"); } f64::NAN } } }, v)
        }
        "bio_amino_acid_composition" => {
            let r = f_aa_composition(match v.get("seq") { Some(x) => py::from_val(x), None => { return None } });
            finish(r, |i: &Py, o: &Py| -> f64 { match s_118(i.clone(), o.clone()).and_then(|r| r.f()) { Ok(x) => x, Err(e) => { if std::env::var_os("AX_DEBUG").is_some() { eprintln!("check error: {e:?}"); } f64::NAN } } }, v)
        }
        "bio_aromaticity" => {
            let r = f_aromaticity(match v.get("seq") { Some(x) => py::from_val(x), None => { return None } });
            finish(r, |i: &Py, o: &Py| -> f64 { match s_182(i.clone(), o.clone()).and_then(|r| r.f()) { Ok(x) => x, Err(e) => { if std::env::var_os("AX_DEBUG").is_some() { eprintln!("check error: {e:?}"); } f64::NAN } } }, v)
        }
        "bio_beta_diversity" => {
            let r = f_beta_diversity(match v.get("x") { Some(x) => py::from_val(x), None => { return None } }, match v.get("y") { Some(x) => py::from_val(x), None => { return None } }, (|| -> R { Ok(Py::Bool(false)) })().ok()?);
            finish(r, |i: &Py, o: &Py| -> f64 { match s_256(i.clone(), o.clone()).and_then(|r| r.f()) { Ok(x) => x, Err(e) => { if std::env::var_os("AX_DEBUG").is_some() { eprintln!("check error: {e:?}"); } f64::NAN } } }, v)
        }
        "bio_brillouin_diversity" => {
            let r = f_brillouin(match v.get("counts") { Some(x) => py::from_val(x), None => { return None } }, (|| -> R { Ok(Py::Bool(false)) })().ok()?);
            finish(r, |i: &Py, o: &Py| -> f64 { match s_330(i.clone(), o.clone()).and_then(|r| r.f()) { Ok(x) => x, Err(e) => { if std::env::var_os("AX_DEBUG").is_some() { eprintln!("check error: {e:?}"); } f64::NAN } } }, v)
        }
        "bio_chao1" => {
            let r = f_chao1(match v.get("counts") { Some(x) => py::from_val(x), None => { return None } }, match v.get("form") { Some(x) => py::from_val(x), None => { return None } }, (|| -> R { Ok(Py::Bool(false)) })().ok()?);
            finish(r, |i: &Py, o: &Py| -> f64 { match s_404(i.clone(), o.clone()).and_then(|r| r.f()) { Ok(x) => x, Err(e) => { if std::env::var_os("AX_DEBUG").is_some() { eprintln!("check error: {e:?}"); } f64::NAN } } }, v)
        }
        "bio_chemostat_steady_state" => {
            let r = f_chemostat(match v.get("mu_max") { Some(x) => py::from_val(x), None => { return None } }, match v.get("Ks") { Some(x) => py::from_val(x), None => { return None } }, match v.get("Y") { Some(x) => py::from_val(x), None => { return None } }, match v.get("D") { Some(x) => py::from_val(x), None => { return None } }, match v.get("S0") { Some(x) => py::from_val(x), None => { return None } }, (|| -> R { Ok(Py::Bool(false)) })().ok()?);
            finish(r, |i: &Py, o: &Py| -> f64 { match s_469(i.clone(), o.clone()).and_then(|r| r.f()) { Ok(x) => x, Err(e) => { if std::env::var_os("AX_DEBUG").is_some() { eprintln!("check error: {e:?}"); } f64::NAN } } }, v)
        }
        "bio_coalescent_times" => {
            let r = f_coalescent_times(match v.get("n") { Some(x) => py::from_val(x), None => { return None } }, match v.get("Ne") { Some(x) => py::from_val(x), None => { return None } });
            finish(r, |i: &Py, o: &Py| -> f64 { match s_535(i.clone(), o.clone()).and_then(|r| r.f()) { Ok(x) => x, Err(e) => { if std::env::var_os("AX_DEBUG").is_some() { eprintln!("check error: {e:?}"); } f64::NAN } } }, v)
        }
        "bio_competitive_inhibition_fit" => {
            let r = f_competitive_fit(match v.get("S") { Some(x) => py::from_val(x), None => { return None } }, match v.get("I") { Some(x) => py::from_val(x), None => { return None } }, match v.get("v") { Some(x) => py::from_val(x), None => { return None } }, (|| -> R { Ok(Py::str("competitive")) })().ok()?);
            finish(r, |i: &Py, o: &Py| -> f64 { match s_776(i.clone(), o.clone()).and_then(|r| r.f()) { Ok(x) => x, Err(e) => { if std::env::var_os("AX_DEBUG").is_some() { eprintln!("check error: {e:?}"); } f64::NAN } } }, v)
        }
        "bio_cross_genotype_probability" => {
            let r = f_cross_genotype(match v.get("parent1") { Some(x) => py::from_val(x), None => { return None } }, match v.get("parent2") { Some(x) => py::from_val(x), None => { return None } }, match v.get("offspring") { Some(x) => py::from_val(x), None => { return None } });
            finish(r, |i: &Py, o: &Py| -> f64 { match s_847(i.clone(), o.clone()).and_then(|r| r.f()) { Ok(x) => x, Err(e) => { if std::env::var_os("AX_DEBUG").is_some() { eprintln!("check error: {e:?}"); } f64::NAN } } }, v)
        }
        "bio_cross_phenotype_ratio" => {
            let r = f_cross_phenotypes(match v.get("parent1") { Some(x) => py::from_val(x), None => { return None } }, match v.get("parent2") { Some(x) => py::from_val(x), None => { return None } }, (|| -> R { Ok(Py::Bool(true)) })().ok()?);
            finish(r, |i: &Py, o: &Py| -> f64 { match s_938(i.clone(), o.clone()).and_then(|r| r.f()) { Ok(x) => x, Err(e) => { if std::env::var_os("AX_DEBUG").is_some() { eprintln!("check error: {e:?}"); } f64::NAN } } }, v)
        }
        "bio_dominance_indices" => {
            let r = f_dominance_indices(match v.get("counts") { Some(x) => py::from_val(x), None => { return None } }, (|| -> R { Ok(Py::Bool(false)) })().ok()?);
            finish(r, |i: &Py, o: &Py| -> f64 { match s_1014(i.clone(), o.clone()).and_then(|r| r.f()) { Ok(x) => x, Err(e) => { if std::env::var_os("AX_DEBUG").is_some() { eprintln!("check error: {e:?}"); } f64::NAN } } }, v)
        }
        "bio_drift_heterozygosity" => {
            let r = f_drift_heterozygosity(match v.get("N") { Some(x) => py::from_val(x), None => { return None } }, match v.get("i") { Some(x) => py::from_val(x), None => { return None } }, match v.get("t") { Some(x) => py::from_val(x), None => { return None } });
            finish(r, |i: &Py, o: &Py| -> f64 { match s_1081(i.clone(), o.clone()).and_then(|r| r.f()) { Ok(x) => x, Err(e) => { if std::env::var_os("AX_DEBUG").is_some() { eprintln!("check error: {e:?}"); } f64::NAN } } }, v)
        }
        "bio_eadie_hofstee" => {
            let r = f_eadie_hofstee(match v.get("S") { Some(x) => py::from_val(x), None => { return None } }, match v.get("v") { Some(x) => py::from_val(x), None => { return None } });
            finish(r, |i: &Py, o: &Py| -> f64 { match s_1170(i.clone(), o.clone()).and_then(|r| r.f()) { Ok(x) => x, Err(e) => { if std::env::var_os("AX_DEBUG").is_some() { eprintln!("check error: {e:?}"); } f64::NAN } } }, v)
        }
        "bio_edit_distance" => {
            let r = f_edit_distance(match v.get("seq1") { Some(x) => py::from_val(x), None => { return None } }, match v.get("seq2") { Some(x) => py::from_val(x), None => { return None } });
            finish(r, |i: &Py, o: &Py| -> f64 { match s_1241(i.clone(), o.clone()).and_then(|r| r.f()) { Ok(x) => x, Err(e) => { if std::env::var_os("AX_DEBUG").is_some() { eprintln!("check error: {e:?}"); } f64::NAN } } }, v)
        }
        "bio_faith_pd" => {
            let r = f_faith_pd(match v.get("tree") { Some(x) => py::from_val(x), None => { return None } }, match v.get("taxa") { Some(x) => py::from_val(x), None => { return None } }, (|| -> R { Ok(Py::Bool(false)) })().ok()?);
            finish(r, |i: &Py, o: &Py| -> f64 { match s_1349(i.clone(), o.clone()).and_then(|r| r.f()) { Ok(x) => x, Err(e) => { if std::env::var_os("AX_DEBUG").is_some() { eprintln!("check error: {e:?}"); } f64::NAN } } }, v)
        }
        "bio_fitch_parsimony" => {
            let r = f_fitch_parsimony(match v.get("tree") { Some(x) => py::from_val(x), None => { return None } }, match v.get("states") { Some(x) => py::from_val(x), None => { return None } }, (|| -> R { Ok(Py::Bool(false)) })().ok()?);
            finish(r, |i: &Py, o: &Py| -> f64 { match s_1462(i.clone(), o.clone()).and_then(|r| r.f()) { Ok(x) => x, Err(e) => { if std::env::var_os("AX_DEBUG").is_some() { eprintln!("check error: {e:?}"); } f64::NAN } } }, v)
        }
        "bio_flux_balance_analysis" => {
            let r = f_fba(match v.get("S") { Some(x) => py::from_val(x), None => { return None } }, match v.get("lb") { Some(x) => py::from_val(x), None => { return None } }, match v.get("ub") { Some(x) => py::from_val(x), None => { return None } }, match v.get("c") { Some(x) => py::from_val(x), None => { return None } }, (|| -> R { Ok(Py::Bool(false)) })().ok()?);
            finish(r, |i: &Py, o: &Py| -> f64 { match s_1602(i.clone(), o.clone()).and_then(|r| r.f()) { Ok(x) => x, Err(e) => { if std::env::var_os("AX_DEBUG").is_some() { eprintln!("check error: {e:?}"); } f64::NAN } } }, v)
        }
        "bio_fst" => {
            let r = f_fst(match v.get("p") { Some(x) => py::from_val(x), None => { return None } });
            finish(r, |i: &Py, o: &Py| -> f64 { match s_1677(i.clone(), o.clone()).and_then(|r| r.f()) { Ok(x) => x, Err(e) => { if std::env::var_os("AX_DEBUG").is_some() { eprintln!("check error: {e:?}"); } f64::NAN } } }, v)
        }
        "bio_gc_content" => {
            let r = f_gc_content(match v.get("seq") { Some(x) => py::from_val(x), None => { return None } });
            finish(r, |i: &Py, o: &Py| -> f64 { match s_1741(i.clone(), o.clone()).and_then(|r| r.f()) { Ok(x) => x, Err(e) => { if std::env::var_os("AX_DEBUG").is_some() { eprintln!("check error: {e:?}"); } f64::NAN } } }, v)
        }
        "bio_global_alignment_affine" => {
            let r = f_global_affine(match v.get("seq1") { Some(x) => py::from_val(x), None => { return None } }, match v.get("seq2") { Some(x) => py::from_val(x), None => { return None } }, match v.get("match") { Some(x) => py::from_val(x), None => { return None } }, match v.get("mismatch") { Some(x) => py::from_val(x), None => { return None } }, match v.get("open") { Some(x) => py::from_val(x), None => { return None } }, match v.get("extend") { Some(x) => py::from_val(x), None => { return None } });
            finish(r, |i: &Py, o: &Py| -> f64 { match s_1863(i.clone(), o.clone()).and_then(|r| r.f()) { Ok(x) => x, Err(e) => { if std::env::var_os("AX_DEBUG").is_some() { eprintln!("check error: {e:?}"); } f64::NAN } } }, v)
        }
        "bio_global_alignment_score" => {
            let r = f_global_linear(match v.get("seq1") { Some(x) => py::from_val(x), None => { return None } }, match v.get("seq2") { Some(x) => py::from_val(x), None => { return None } }, match v.get("match") { Some(x) => py::from_val(x), None => { return None } }, match v.get("mismatch") { Some(x) => py::from_val(x), None => { return None } }, match v.get("gap") { Some(x) => py::from_val(x), None => { return None } });
            finish(r, |i: &Py, o: &Py| -> f64 { match s_1982(i.clone(), o.clone()).and_then(|r| r.f()) { Ok(x) => x, Err(e) => { if std::env::var_os("AX_DEBUG").is_some() { eprintln!("check error: {e:?}"); } f64::NAN } } }, v)
        }
        "bio_goods_coverage" => {
            let r = f_goods_coverage(match v.get("counts") { Some(x) => py::from_val(x), None => { return None } });
            finish(r, |i: &Py, o: &Py| -> f64 { match s_2044(i.clone(), o.clone()).and_then(|r| r.f()) { Ok(x) => x, Err(e) => { if std::env::var_os("AX_DEBUG").is_some() { eprintln!("check error: {e:?}"); } f64::NAN } } }, v)
        }
        "bio_gravy" => {
            let r = f_gravy(match v.get("seq") { Some(x) => py::from_val(x), None => { return None } });
            finish(r, |i: &Py, o: &Py| -> f64 { match s_2123(i.clone(), o.clone()).and_then(|r| r.f()) { Ok(x) => x, Err(e) => { if std::env::var_os("AX_DEBUG").is_some() { eprintln!("check error: {e:?}"); } f64::NAN } } }, v)
        }
        "bio_hamming_distance" => {
            let r = f_hamming(match v.get("seq1") { Some(x) => py::from_val(x), None => { return None } }, match v.get("seq2") { Some(x) => py::from_val(x), None => { return None } });
            finish(r, |i: &Py, o: &Py| -> f64 { match s_2194(i.clone(), o.clone()).and_then(|r| r.f()) { Ok(x) => x, Err(e) => { if std::env::var_os("AX_DEBUG").is_some() { eprintln!("check error: {e:?}"); } f64::NAN } } }, v)
        }
        "bio_hanes_woolf" => {
            let r = f_hanes_woolf(match v.get("S") { Some(x) => py::from_val(x), None => { return None } }, match v.get("v") { Some(x) => py::from_val(x), None => { return None } });
            finish(r, |i: &Py, o: &Py| -> f64 { match s_2283(i.clone(), o.clone()).and_then(|r| r.f()) { Ok(x) => x, Err(e) => { if std::env::var_os("AX_DEBUG").is_some() { eprintln!("check error: {e:?}"); } f64::NAN } } }, v)
        }
        "bio_haploid_selection" => {
            let r = f_haploid_selection(match v.get("p0") { Some(x) => py::from_val(x), None => { return None } }, match v.get("s") { Some(x) => py::from_val(x), None => { return None } }, match v.get("t") { Some(x) => py::from_val(x), None => { return None } });
            finish(r, |i: &Py, o: &Py| -> f64 { match s_2357(i.clone(), o.clone()).and_then(|r| r.f()) { Ok(x) => x, Err(e) => { if std::env::var_os("AX_DEBUG").is_some() { eprintln!("check error: {e:?}"); } f64::NAN } } }, v)
        }
        "bio_hardy_weinberg_exact_test" => {
            let r = f_hwe_exact(match v.get("n_AA") { Some(x) => py::from_val(x), None => { return None } }, match v.get("n_Aa") { Some(x) => py::from_val(x), None => { return None } }, match v.get("n_aa") { Some(x) => py::from_val(x), None => { return None } });
            finish(r, |i: &Py, o: &Py| -> f64 { match s_2414(i.clone(), o.clone()).and_then(|r| r.f()) { Ok(x) => x, Err(e) => { if std::env::var_os("AX_DEBUG").is_some() { eprintln!("check error: {e:?}"); } f64::NAN } } }, v)
        }
        "bio_hardy_weinberg_test" => {
            let r = f_hardy_weinberg_test(match v.get("n_AA") { Some(x) => py::from_val(x), None => { return None } }, match v.get("n_Aa") { Some(x) => py::from_val(x), None => { return None } }, match v.get("n_aa") { Some(x) => py::from_val(x), None => { return None } });
            finish(r, |i: &Py, o: &Py| -> f64 { match s_2482(i.clone(), o.clone()).and_then(|r| r.f()) { Ok(x) => x, Err(e) => { if std::env::var_os("AX_DEBUG").is_some() { eprintln!("check error: {e:?}"); } f64::NAN } } }, v)
        }
        "bio_heterozygote_advantage" => {
            let r = f_het_advantage(match v.get("w_AA") { Some(x) => py::from_val(x), None => { return None } }, match v.get("w_Aa") { Some(x) => py::from_val(x), None => { return None } }, match v.get("w_aa") { Some(x) => py::from_val(x), None => { return None } });
            finish(r, |i: &Py, o: &Py| -> f64 { match s_2548(i.clone(), o.clone()).and_then(|r| r.f()) { Ok(x) => x, Err(e) => { if std::env::var_os("AX_DEBUG").is_some() { eprintln!("check error: {e:?}"); } f64::NAN } } }, v)
        }
        "bio_hill_fit" => {
            let r = f_hill_fit(match v.get("S") { Some(x) => py::from_val(x), None => { return None } }, match v.get("v") { Some(x) => py::from_val(x), None => { return None } }, (|| -> R { Ok(Py::Bool(false)) })().ok()?);
            finish(r, |i: &Py, o: &Py| -> f64 { match s_2764(i.clone(), o.clone()).and_then(|r| r.f()) { Ok(x) => x, Err(e) => { if std::env::var_os("AX_DEBUG").is_some() { eprintln!("check error: {e:?}"); } f64::NAN } } }, v)
        }
        "bio_hill_plot" => {
            let r = f_hill_plot(match v.get("S") { Some(x) => py::from_val(x), None => { return None } }, match v.get("Y") { Some(x) => py::from_val(x), None => { return None } });
            finish(r, |i: &Py, o: &Py| -> f64 { match s_2860(i.clone(), o.clone()).and_then(|r| r.f()) { Ok(x) => x, Err(e) => { if std::env::var_os("AX_DEBUG").is_some() { eprintln!("check error: {e:?}"); } f64::NAN } } }, v)
        }
        "bio_instability_index" => {
            let r = f_instability(match v.get("seq") { Some(x) => py::from_val(x), None => { return None } });
            finish(r, |i: &Py, o: &Py| -> f64 { match s_2941(i.clone(), o.clone()).and_then(|r| r.f()) { Ok(x) => x, Err(e) => { if std::env::var_os("AX_DEBUG").is_some() { eprintln!("check error: {e:?}"); } f64::NAN } } }, v)
        }
        "bio_isoelectric_point" => {
            let r = f_isoelectric_point(match v.get("seq") { Some(x) => py::from_val(x), None => { return None } });
            finish(r, |i: &Py, o: &Py| -> f64 { match s_3056(i.clone(), o.clone()).and_then(|r| r.f()) { Ok(x) => x, Err(e) => { if std::env::var_os("AX_DEBUG").is_some() { eprintln!("check error: {e:?}"); } f64::NAN } } }, v)
        }
        "bio_jukes_cantor_distance" => {
            let r = f_jukes_cantor(match v.get("seq1") { Some(x) => py::from_val(x), None => { return None } }, match v.get("seq2") { Some(x) => py::from_val(x), None => { return None } });
            finish(r, |i: &Py, o: &Py| -> f64 { match s_3144(i.clone(), o.clone()).and_then(|r| r.f()) { Ok(x) => x, Err(e) => { if std::env::var_os("AX_DEBUG").is_some() { eprintln!("check error: {e:?}"); } f64::NAN } } }, v)
        }
        "bio_kimura_distance" => {
            let r = f_kimura(match v.get("seq1") { Some(x) => py::from_val(x), None => { return None } }, match v.get("seq2") { Some(x) => py::from_val(x), None => { return None } });
            finish(r, |i: &Py, o: &Py| -> f64 { match s_3251(i.clone(), o.clone()).and_then(|r| r.f()) { Ok(x) => x, Err(e) => { if std::env::var_os("AX_DEBUG").is_some() { eprintln!("check error: {e:?}"); } f64::NAN } } }, v)
        }
        "bio_kmer_counts" => {
            let r = f_kmer_counts(match v.get("seq") { Some(x) => py::from_val(x), None => { return None } }, match v.get("k") { Some(x) => py::from_val(x), None => { return None } });
            finish(r, |i: &Py, o: &Py| -> f64 { match s_3324(i.clone(), o.clone()).and_then(|r| r.f()) { Ok(x) => x, Err(e) => { if std::env::var_os("AX_DEBUG").is_some() { eprintln!("check error: {e:?}"); } f64::NAN } } }, v)
        }
        "bio_leslie_growth_rate" => {
            let r = f_leslie_growth(match v.get("F") { Some(x) => py::from_val(x), None => { return None } }, match v.get("S") { Some(x) => py::from_val(x), None => { return None } }, (|| -> R { Ok(Py::Bool(false)) })().ok()?);
            finish(r, |i: &Py, o: &Py| -> f64 { match s_3438(i.clone(), o.clone()).and_then(|r| r.f()) { Ok(x) => x, Err(e) => { if std::env::var_os("AX_DEBUG").is_some() { eprintln!("check error: {e:?}"); } f64::NAN } } }, v)
        }
        "bio_leslie_projection" => {
            let r = f_leslie_projection(match v.get("F") { Some(x) => py::from_val(x), None => { return None } }, match v.get("S") { Some(x) => py::from_val(x), None => { return None } }, match v.get("n0") { Some(x) => py::from_val(x), None => { return None } }, match v.get("t") { Some(x) => py::from_val(x), None => { return None } }, (|| -> R { Ok(Py::Bool(false)) })().ok()?);
            finish(r, |i: &Py, o: &Py| -> f64 { match s_3532(i.clone(), o.clone()).and_then(|r| r.f()) { Ok(x) => x, Err(e) => { if std::env::var_os("AX_DEBUG").is_some() { eprintln!("check error: {e:?}"); } f64::NAN } } }, v)
        }
        "bio_lineweaver_burk" => {
            let r = f_lineweaver_burk(match v.get("S") { Some(x) => py::from_val(x), None => { return None } }, match v.get("v") { Some(x) => py::from_val(x), None => { return None } });
            finish(r, |i: &Py, o: &Py| -> f64 { match s_3621(i.clone(), o.clone()).and_then(|r| r.f()) { Ok(x) => x, Err(e) => { if std::env::var_os("AX_DEBUG").is_some() { eprintln!("check error: {e:?}"); } f64::NAN } } }, v)
        }
        "bio_linkage_disequilibrium" => {
            let r = f_linkage_disequilibrium(match v.get("n_AB") { Some(x) => py::from_val(x), None => { return None } }, match v.get("n_Ab") { Some(x) => py::from_val(x), None => { return None } }, match v.get("n_aB") { Some(x) => py::from_val(x), None => { return None } }, match v.get("n_ab") { Some(x) => py::from_val(x), None => { return None } });
            finish(r, |i: &Py, o: &Py| -> f64 { match s_3680(i.clone(), o.clone()).and_then(|r| r.f()) { Ok(x) => x, Err(e) => { if std::env::var_os("AX_DEBUG").is_some() { eprintln!("check error: {e:?}"); } f64::NAN } } }, v)
        }
        "bio_local_alignment_affine" => {
            let r = f_local_affine(match v.get("seq1") { Some(x) => py::from_val(x), None => { return None } }, match v.get("seq2") { Some(x) => py::from_val(x), None => { return None } }, match v.get("match") { Some(x) => py::from_val(x), None => { return None } }, match v.get("mismatch") { Some(x) => py::from_val(x), None => { return None } }, match v.get("open") { Some(x) => py::from_val(x), None => { return None } }, match v.get("extend") { Some(x) => py::from_val(x), None => { return None } });
            finish(r, |i: &Py, o: &Py| -> f64 { match s_3804(i.clone(), o.clone()).and_then(|r| r.f()) { Ok(x) => x, Err(e) => { if std::env::var_os("AX_DEBUG").is_some() { eprintln!("check error: {e:?}"); } f64::NAN } } }, v)
        }
        "bio_local_alignment_score" => {
            let r = f_local_linear(match v.get("seq1") { Some(x) => py::from_val(x), None => { return None } }, match v.get("seq2") { Some(x) => py::from_val(x), None => { return None } }, match v.get("match") { Some(x) => py::from_val(x), None => { return None } }, match v.get("mismatch") { Some(x) => py::from_val(x), None => { return None } }, match v.get("gap") { Some(x) => py::from_val(x), None => { return None } });
            finish(r, |i: &Py, o: &Py| -> f64 { match s_3925(i.clone(), o.clone()).and_then(|r| r.f()) { Ok(x) => x, Err(e) => { if std::env::var_os("AX_DEBUG").is_some() { eprintln!("check error: {e:?}"); } f64::NAN } } }, v)
        }
        "bio_longest_orf" => {
            let r = f_longest_orf(match v.get("seq") { Some(x) => py::from_val(x), None => { return None } });
            finish(r, |i: &Py, o: &Py| -> f64 { match s_3997(i.clone(), o.clone()).and_then(|r| r.f()) { Ok(x) => x, Err(e) => { if std::env::var_os("AX_DEBUG").is_some() { eprintln!("check error: {e:?}"); } f64::NAN } } }, v)
        }
        "bio_map_distance" => {
            let r = f_map_distance(match v.get("r") { Some(x) => py::from_val(x), None => { return None } }, match v.get("map") { Some(x) => py::from_val(x), None => { return None } });
            finish(r, |i: &Py, o: &Py| -> f64 { match s_4062(i.clone(), o.clone()).and_then(|r| r.f()) { Ok(x) => x, Err(e) => { if std::env::var_os("AX_DEBUG").is_some() { eprintln!("check error: {e:?}"); } f64::NAN } } }, v)
        }
        "bio_michaelis_menten_fit" => {
            let r = f_mm_fit(match v.get("S") { Some(x) => py::from_val(x), None => { return None } }, match v.get("v") { Some(x) => py::from_val(x), None => { return None } });
            finish(r, |i: &Py, o: &Py| -> f64 { match s_4179(i.clone(), o.clone()).and_then(|r| r.f()) { Ok(x) => x, Err(e) => { if std::env::var_os("AX_DEBUG").is_some() { eprintln!("check error: {e:?}"); } f64::NAN } } }, v)
        }
        "bio_michaelis_menten_progress" => {
            let r = f_mm_progress(match v.get("S0") { Some(x) => py::from_val(x), None => { return None } }, match v.get("Vmax") { Some(x) => py::from_val(x), None => { return None } }, match v.get("Km") { Some(x) => py::from_val(x), None => { return None } }, match v.get("t") { Some(x) => py::from_val(x), None => { return None } }, (|| -> R { Ok(Py::Bool(false)) })().ok()?);
            finish(r, |i: &Py, o: &Py| -> f64 { match s_4277(i.clone(), o.clone()).and_then(|r| r.f()) { Ok(x) => x, Err(e) => { if std::env::var_os("AX_DEBUG").is_some() { eprintln!("check error: {e:?}"); } f64::NAN } } }, v)
        }
        "bio_michaelis_menten_time" => {
            let r = f_mm_time(match v.get("S0") { Some(x) => py::from_val(x), None => { return None } }, match v.get("Vmax") { Some(x) => py::from_val(x), None => { return None } }, match v.get("Km") { Some(x) => py::from_val(x), None => { return None } }, match v.get("X") { Some(x) => py::from_val(x), None => { return None } }, (|| -> R { Ok(Py::Bool(false)) })().ok()?);
            finish(r, |i: &Py, o: &Py| -> f64 { match s_4344(i.clone(), o.clone()).and_then(|r| r.f()) { Ok(x) => x, Err(e) => { if std::env::var_os("AX_DEBUG").is_some() { eprintln!("check error: {e:?}"); } f64::NAN } } }, v)
        }
        "bio_moran_fixation" => {
            let r = f_moran_fixation(match v.get("N") { Some(x) => py::from_val(x), None => { return None } }, match v.get("i") { Some(x) => py::from_val(x), None => { return None } }, match v.get("r") { Some(x) => py::from_val(x), None => { return None } });
            finish(r, |i: &Py, o: &Py| -> f64 { match s_4419(i.clone(), o.clone()).and_then(|r| r.f()) { Ok(x) => x, Err(e) => { if std::env::var_os("AX_DEBUG").is_some() { eprintln!("check error: {e:?}"); } f64::NAN } } }, v)
        }
        "bio_mutation_selection_balance" => {
            let r = f_mutation_selection(match v.get("mu") { Some(x) => py::from_val(x), None => { return None } }, match v.get("s") { Some(x) => py::from_val(x), None => { return None } }, match v.get("h") { Some(x) => py::from_val(x), None => { return None } });
            finish(r, |i: &Py, o: &Py| -> f64 { match s_4516(i.clone(), o.clone()).and_then(|r| r.f()) { Ok(x) => x, Err(e) => { if std::env::var_os("AX_DEBUG").is_some() { eprintln!("check error: {e:?}"); } f64::NAN } } }, v)
        }
        "bio_ne_fluctuating_size" => {
            let r = f_ne_fluctuating(match v.get("sizes") { Some(x) => py::from_val(x), None => { return None } });
            finish(r, |i: &Py, o: &Py| -> f64 { match s_4585(i.clone(), o.clone()).and_then(|r| r.f()) { Ok(x) => x, Err(e) => { if std::env::var_os("AX_DEBUG").is_some() { eprintln!("check error: {e:?}"); } f64::NAN } } }, v)
        }
        "bio_ne_sex_ratio" => {
            let r = f_ne_sex_ratio(match v.get("N_m") { Some(x) => py::from_val(x), None => { return None } }, match v.get("N_f") { Some(x) => py::from_val(x), None => { return None } });
            finish(r, |i: &Py, o: &Py| -> f64 { match s_4633(i.clone(), o.clone()).and_then(|r| r.f()) { Ok(x) => x, Err(e) => { if std::env::var_os("AX_DEBUG").is_some() { eprintln!("check error: {e:?}"); } f64::NAN } } }, v)
        }
        "bio_neighbor_joining" => {
            let r = f_neighbor_joining(match v.get("D") { Some(x) => py::from_val(x), None => { return None } }, (|| -> R { Ok(Py::Bool(true)) })().ok()?);
            finish(r, |i: &Py, o: &Py| -> f64 { match s_4753(i.clone(), o.clone()).and_then(|r| r.f()) { Ok(x) => x, Err(e) => { if std::env::var_os("AX_DEBUG").is_some() { eprintln!("check error: {e:?}"); } f64::NAN } } }, v)
        }
        "bio_nucleic_acid_mass" => {
            let r = f_nucleic_mass(match v.get("seq") { Some(x) => py::from_val(x), None => { return None } }, match v.get("molecule") { Some(x) => py::from_val(x), None => { return None } }, match v.get("strands") { Some(x) => py::from_val(x), None => { return None } });
            finish(r, |i: &Py, o: &Py| -> f64 { match s_4842(i.clone(), o.clone()).and_then(|r| r.f()) { Ok(x) => x, Err(e) => { if std::env::var_os("AX_DEBUG").is_some() { eprintln!("check error: {e:?}"); } f64::NAN } } }, v)
        }
        "bio_nucleotide_diversity" => {
            let r = f_nucleotide_diversity(match v.get("seqs") { Some(x) => py::from_val(x), None => { return None } });
            finish(r, |i: &Py, o: &Py| -> f64 { match s_4898(i.clone(), o.clone()).and_then(|r| r.f()) { Ok(x) => x, Err(e) => { if std::env::var_os("AX_DEBUG").is_some() { eprintln!("check error: {e:?}"); } f64::NAN } } }, v)
        }
        "bio_pedigree_inbreeding" => {
            let r = f_pedigree_inbreeding(match v.get("pedigree") { Some(x) => py::from_val(x), None => { return None } }, match v.get("individual") { Some(x) => py::from_val(x), None => { return None } });
            finish(r, |i: &Py, o: &Py| -> f64 { match s_4980(i.clone(), o.clone()).and_then(|r| r.f()) { Ok(x) => x, Err(e) => { if std::env::var_os("AX_DEBUG").is_some() { eprintln!("check error: {e:?}"); } f64::NAN } } }, v)
        }
        "bio_petersen_estimate" => {
            let r = f_petersen_mle(match v.get("M") { Some(x) => py::from_val(x), None => { return None } }, match v.get("C") { Some(x) => py::from_val(x), None => { return None } }, match v.get("R") { Some(x) => py::from_val(x), None => { return None } }, (|| -> R { Ok(Py::Bool(false)) })().ok()?);
            finish(r, |i: &Py, o: &Py| -> f64 { match s_5028(i.clone(), o.clone()).and_then(|r| r.f()) { Ok(x) => x, Err(e) => { if std::env::var_os("AX_DEBUG").is_some() { eprintln!("check error: {e:?}"); } f64::NAN } } }, v)
        }
        "bio_protein_extinction" => {
            let r = f_extinction(match v.get("seq") { Some(x) => py::from_val(x), None => { return None } });
            finish(r, |i: &Py, o: &Py| -> f64 { match s_5107(i.clone(), o.clone()).and_then(|r| r.f()) { Ok(x) => x, Err(e) => { if std::env::var_os("AX_DEBUG").is_some() { eprintln!("check error: {e:?}"); } f64::NAN } } }, v)
        }
        "bio_protein_global_blosum62" => {
            let r = f_protein_global(match v.get("seq1") { Some(x) => py::from_val(x), None => { return None } }, match v.get("seq2") { Some(x) => py::from_val(x), None => { return None } }, match v.get("open") { Some(x) => py::from_val(x), None => { return None } }, match v.get("extend") { Some(x) => py::from_val(x), None => { return None } });
            finish(r, |i: &Py, o: &Py| -> f64 { match s_5227(i.clone(), o.clone()).and_then(|r| r.f()) { Ok(x) => x, Err(e) => { if std::env::var_os("AX_DEBUG").is_some() { eprintln!("check error: {e:?}"); } f64::NAN } } }, v)
        }
        "bio_protein_local_blosum62" => {
            let r = f_protein_local(match v.get("seq1") { Some(x) => py::from_val(x), None => { return None } }, match v.get("seq2") { Some(x) => py::from_val(x), None => { return None } }, match v.get("open") { Some(x) => py::from_val(x), None => { return None } }, match v.get("extend") { Some(x) => py::from_val(x), None => { return None } });
            finish(r, |i: &Py, o: &Py| -> f64 { match s_5347(i.clone(), o.clone()).and_then(|r| r.f()) { Ok(x) => x, Err(e) => { if std::env::var_os("AX_DEBUG").is_some() { eprintln!("check error: {e:?}"); } f64::NAN } } }, v)
        }
        "bio_protein_mass" => {
            let r = f_protein_mass(match v.get("seq") { Some(x) => py::from_val(x), None => { return None } });
            finish(r, |i: &Py, o: &Py| -> f64 { match s_5426(i.clone(), o.clone()).and_then(|r| r.f()) { Ok(x) => x, Err(e) => { if std::env::var_os("AX_DEBUG").is_some() { eprintln!("check error: {e:?}"); } f64::NAN } } }, v)
        }
        "bio_protein_net_charge" => {
            let r = f_net_charge(match v.get("seq") { Some(x) => py::from_val(x), None => { return None } }, match v.get("pH") { Some(x) => py::from_val(x), None => { return None } });
            finish(r, |i: &Py, o: &Py| -> f64 { match s_5521(i.clone(), o.clone()).and_then(|r| r.f()) { Ok(x) => x, Err(e) => { if std::env::var_os("AX_DEBUG").is_some() { eprintln!("check error: {e:?}"); } f64::NAN } } }, v)
        }
        "bio_rarefaction" => {
            let r = f_rarefaction(match v.get("counts") { Some(x) => py::from_val(x), None => { return None } }, match v.get("n") { Some(x) => py::from_val(x), None => { return None } }, (|| -> R { Ok(Py::Bool(false)) })().ok()?);
            finish(r, |i: &Py, o: &Py| -> f64 { match s_5588(i.clone(), o.clone()).and_then(|r| r.f()) { Ok(x) => x, Err(e) => { if std::env::var_os("AX_DEBUG").is_some() { eprintln!("check error: {e:?}"); } f64::NAN } } }, v)
        }
        "bio_recessive_lethal" => {
            let r = f_recessive_lethal(match v.get("q0") { Some(x) => py::from_val(x), None => { return None } }, match v.get("t") { Some(x) => py::from_val(x), None => { return None } });
            finish(r, |i: &Py, o: &Py| -> f64 { match s_5660(i.clone(), o.clone()).and_then(|r| r.f()) { Ok(x) => x, Err(e) => { if std::env::var_os("AX_DEBUG").is_some() { eprintln!("check error: {e:?}"); } f64::NAN } } }, v)
        }
        "bio_recombination_fraction" => {
            let r = f_recombination_fraction(match v.get("d_cM") { Some(x) => py::from_val(x), None => { return None } }, match v.get("map") { Some(x) => py::from_val(x), None => { return None } });
            finish(r, |i: &Py, o: &Py| -> f64 { match s_5724(i.clone(), o.clone()).and_then(|r| r.f()) { Ok(x) => x, Err(e) => { if std::env::var_os("AX_DEBUG").is_some() { eprintln!("check error: {e:?}"); } f64::NAN } } }, v)
        }
        "bio_restriction_sites" => {
            let r = f_restriction_sites(match v.get("seq") { Some(x) => py::from_val(x), None => { return None } }, match v.get("enzyme") { Some(x) => py::from_val(x), None => { return None } });
            finish(r, |i: &Py, o: &Py| -> f64 { match s_5793(i.clone(), o.clone()).and_then(|r| r.f()) { Ok(x) => x, Err(e) => { if std::env::var_os("AX_DEBUG").is_some() { eprintln!("check error: {e:?}"); } f64::NAN } } }, v)
        }
        "bio_reverse_complement" => {
            let r = f_reverse_complement(match v.get("seq") { Some(x) => py::from_val(x), None => { return None } });
            finish(r, |i: &Py, o: &Py| -> f64 { match s_5857(i.clone(), o.clone()).and_then(|r| r.f()) { Ok(x) => x, Err(e) => { if std::env::var_os("AX_DEBUG").is_some() { eprintln!("check error: {e:?}"); } f64::NAN } } }, v)
        }
        "bio_richness_indices" => {
            let r = f_richness_indices(match v.get("counts") { Some(x) => py::from_val(x), None => { return None } }, (|| -> R { Ok(Py::Bool(false)) })().ok()?);
            finish(r, |i: &Py, o: &Py| -> f64 { match s_5933(i.clone(), o.clone()).and_then(|r| r.f()) { Ok(x) => x, Err(e) => { if std::env::var_os("AX_DEBUG").is_some() { eprintln!("check error: {e:?}"); } f64::NAN } } }, v)
        }
        "bio_robinson_foulds" => {
            let r = f_robinson_foulds(match v.get("tree1") { Some(x) => py::from_val(x), None => { return None } }, match v.get("tree2") { Some(x) => py::from_val(x), None => { return None } }, (|| -> R { Ok(Py::Bool(false)) })().ok()?);
            finish(r, |i: &Py, o: &Py| -> f64 { match s_6033(i.clone(), o.clone()).and_then(|r| r.f()) { Ok(x) => x, Err(e) => { if std::env::var_os("AX_DEBUG").is_some() { eprintln!("check error: {e:?}"); } f64::NAN } } }, v)
        }
        "bio_schumacher_eschmeyer" => {
            let r = f_schumacher(match v.get("C") { Some(x) => py::from_val(x), None => { return None } }, match v.get("R") { Some(x) => py::from_val(x), None => { return None } }, match v.get("M") { Some(x) => py::from_val(x), None => { return None } }, (|| -> R { Ok(Py::Bool(false)) })().ok()?);
            finish(r, |i: &Py, o: &Py| -> f64 { match s_6108(i.clone(), o.clone()).and_then(|r| r.f()) { Ok(x) => x, Err(e) => { if std::env::var_os("AX_DEBUG").is_some() { eprintln!("check error: {e:?}"); } f64::NAN } } }, v)
        }
        "bio_secondary_structure_fraction" => {
            let r = f_ss_fraction(match v.get("seq") { Some(x) => py::from_val(x), None => { return None } }, (|| -> R { Ok(Py::tuple(vec![Py::str("helix"), Py::str("turn"), Py::str("sheet")])) })().ok()?);
            finish(r, |i: &Py, o: &Py| -> f64 { match s_6175(i.clone(), o.clone()).and_then(|r| r.f()) { Ok(x) => x, Err(e) => { if std::env::var_os("AX_DEBUG").is_some() { eprintln!("check error: {e:?}"); } f64::NAN } } }, v)
        }
        "bio_semiglobal_alignment_score" => {
            let r = f_semiglobal(match v.get("seq1") { Some(x) => py::from_val(x), None => { return None } }, match v.get("seq2") { Some(x) => py::from_val(x), None => { return None } }, match v.get("match") { Some(x) => py::from_val(x), None => { return None } }, match v.get("mismatch") { Some(x) => py::from_val(x), None => { return None } }, match v.get("gap") { Some(x) => py::from_val(x), None => { return None } });
            finish(r, |i: &Py, o: &Py| -> f64 { match s_6294(i.clone(), o.clone()).and_then(|r| r.f()) { Ok(x) => x, Err(e) => { if std::env::var_os("AX_DEBUG").is_some() { eprintln!("check error: {e:?}"); } f64::NAN } } }, v)
        }
        "bio_shannon_diversity" => {
            let r = f_shannon(match v.get("counts") { Some(x) => py::from_val(x), None => { return None } }, match v.get("base") { Some(x) => py::from_val(x), None => { return None } }, (|| -> R { Ok(Py::Bool(false)) })().ok()?);
            finish(r, |i: &Py, o: &Py| -> f64 { match s_6372(i.clone(), o.clone()).and_then(|r| r.f()) { Ok(x) => x, Err(e) => { if std::env::var_os("AX_DEBUG").is_some() { eprintln!("check error: {e:?}"); } f64::NAN } } }, v)
        }
        "bio_simpson_diversity" => {
            let r = f_simpson(match v.get("counts") { Some(x) => py::from_val(x), None => { return None } }, match v.get("estimator") { Some(x) => py::from_val(x), None => { return None } }, (|| -> R { Ok(Py::Bool(false)) })().ok()?);
            finish(r, |i: &Py, o: &Py| -> f64 { match s_6441(i.clone(), o.clone()).and_then(|r| r.f()) { Ok(x) => x, Err(e) => { if std::env::var_os("AX_DEBUG").is_some() { eprintln!("check error: {e:?}"); } f64::NAN } } }, v)
        }
        "bio_species_area" => {
            let r = f_species_area(match v.get("A") { Some(x) => py::from_val(x), None => { return None } }, match v.get("S") { Some(x) => py::from_val(x), None => { return None } }, (|| -> R { Ok(Py::Bool(false)) })().ok()?);
            finish(r, |i: &Py, o: &Py| -> f64 { match s_6520(i.clone(), o.clone()).and_then(|r| r.f()) { Ok(x) => x, Err(e) => { if std::env::var_os("AX_DEBUG").is_some() { eprintln!("check error: {e:?}"); } f64::NAN } } }, v)
        }
        "bio_tajima_d" => {
            let r = f_tajima_d(match v.get("seqs") { Some(x) => py::from_val(x), None => { return None } }, (|| -> R { Ok(Py::Bool(true)) })().ok()?);
            finish(r, |i: &Py, o: &Py| -> f64 { match s_6598(i.clone(), o.clone()).and_then(|r| r.f()) { Ok(x) => x, Err(e) => { if std::env::var_os("AX_DEBUG").is_some() { eprintln!("check error: {e:?}"); } f64::NAN } } }, v)
        }
        "bio_tm_nearest_neighbor" => {
            let r = f_tm_nn(match v.get("seq") { Some(x) => py::from_val(x), None => { return None } }, match v.get("Na") { Some(x) => py::from_val(x), None => { return None } }, match v.get("dnac1") { Some(x) => py::from_val(x), None => { return None } }, match v.get("dnac2") { Some(x) => py::from_val(x), None => { return None } });
            finish(r, |i: &Py, o: &Py| -> f64 { match s_6705(i.clone(), o.clone()).and_then(|r| r.f()) { Ok(x) => x, Err(e) => { if std::env::var_os("AX_DEBUG").is_some() { eprintln!("check error: {e:?}"); } f64::NAN } } }, v)
        }
        "bio_tm_wallace" => {
            let r = f_tm_wallace(match v.get("seq") { Some(x) => py::from_val(x), None => { return None } });
            finish(r, |i: &Py, o: &Py| -> f64 { match s_6769(i.clone(), o.clone()).and_then(|r| r.f()) { Ok(x) => x, Err(e) => { if std::env::var_os("AX_DEBUG").is_some() { eprintln!("check error: {e:?}"); } f64::NAN } } }, v)
        }
        "bio_transcription" => {
            let r = f_transcription(match v.get("seq") { Some(x) => py::from_val(x), None => { return None } }, match v.get("strand") { Some(x) => py::from_val(x), None => { return None } });
            finish(r, |i: &Py, o: &Py| -> f64 { match s_6837(i.clone(), o.clone()).and_then(|r| r.f()) { Ok(x) => x, Err(e) => { if std::env::var_os("AX_DEBUG").is_some() { eprintln!("check error: {e:?}"); } f64::NAN } } }, v)
        }
        "bio_translation" => {
            let r = f_translation(match v.get("seq") { Some(x) => py::from_val(x), None => { return None } }, match v.get("stops") { Some(x) => py::from_val(x), None => { return None } });
            finish(r, |i: &Py, o: &Py| -> f64 { match s_6911(i.clone(), o.clone()).and_then(|r| r.f()) { Ok(x) => x, Err(e) => { if std::env::var_os("AX_DEBUG").is_some() { eprintln!("check error: {e:?}"); } f64::NAN } } }, v)
        }
        "bio_tree_count" => {
            let r = f_tree_count(match v.get("n") { Some(x) => py::from_val(x), None => { return None } }, match v.get("rooted") { Some(x) => py::from_val(x), None => { return None } });
            finish(r, |i: &Py, o: &Py| -> f64 { match s_6964(i.clone(), o.clone()).and_then(|r| r.f()) { Ok(x) => x, Err(e) => { if std::env::var_os("AX_DEBUG").is_some() { eprintln!("check error: {e:?}"); } f64::NAN } } }, v)
        }
        "bio_upgma" => {
            let r = f_upgma(match v.get("D") { Some(x) => py::from_val(x), None => { return None } }, (|| -> R { Ok(Py::str("average")) })().ok()?);
            finish(r, |i: &Py, o: &Py| -> f64 { match s_7065(i.clone(), o.clone()).and_then(|r| r.f()) { Ok(x) => x, Err(e) => { if std::env::var_os("AX_DEBUG").is_some() { eprintln!("check error: {e:?}"); } f64::NAN } } }, v)
        }
        "bio_watterson_theta" => {
            let r = f_watterson(match v.get("seqs") { Some(x) => py::from_val(x), None => { return None } }, (|| -> R { Ok(Py::Bool(false)) })().ok()?);
            finish(r, |i: &Py, o: &Py| -> f64 { match s_7134(i.clone(), o.clone()).and_then(|r| r.f()) { Ok(x) => x, Err(e) => { if std::env::var_os("AX_DEBUG").is_some() { eprintln!("check error: {e:?}"); } f64::NAN } } }, v)
        }
        _ => None,
    }
}

#[allow(non_snake_case)]
fn c_AMINO() -> R {
    thread_local!(static V: R = (|| -> R { Ok(Py::str("ACDEFGHIKLMNPQRSTVWY")) })());
    V.with(|v| v.clone())
}

#[allow(non_snake_case)]
fn c_BLOSUM62() -> R {
    thread_local!(static V: R = (|| -> R { Ok(Py::dict(vec![(Py::tuple(vec![Py::str("A"), Py::str("A")]), Py::Float(4.0)), (Py::tuple(vec![Py::str("A"), Py::str("R")]), Py::Float(-1.0)), (Py::tuple(vec![Py::str("A"), Py::str("N")]), Py::Float(-2.0)), (Py::tuple(vec![Py::str("A"), Py::str("D")]), Py::Float(-2.0)), (Py::tuple(vec![Py::str("A"), Py::str("C")]), Py::Float(0.0)), (Py::tuple(vec![Py::str("A"), Py::str("Q")]), Py::Float(-1.0)), (Py::tuple(vec![Py::str("A"), Py::str("E")]), Py::Float(-1.0)), (Py::tuple(vec![Py::str("A"), Py::str("G")]), Py::Float(0.0)), (Py::tuple(vec![Py::str("A"), Py::str("H")]), Py::Float(-2.0)), (Py::tuple(vec![Py::str("A"), Py::str("I")]), Py::Float(-1.0)), (Py::tuple(vec![Py::str("A"), Py::str("L")]), Py::Float(-1.0)), (Py::tuple(vec![Py::str("A"), Py::str("K")]), Py::Float(-1.0)), (Py::tuple(vec![Py::str("A"), Py::str("M")]), Py::Float(-1.0)), (Py::tuple(vec![Py::str("A"), Py::str("F")]), Py::Float(-2.0)), (Py::tuple(vec![Py::str("A"), Py::str("P")]), Py::Float(-1.0)), (Py::tuple(vec![Py::str("A"), Py::str("S")]), Py::Float(1.0)), (Py::tuple(vec![Py::str("A"), Py::str("T")]), Py::Float(0.0)), (Py::tuple(vec![Py::str("A"), Py::str("W")]), Py::Float(-3.0)), (Py::tuple(vec![Py::str("A"), Py::str("Y")]), Py::Float(-2.0)), (Py::tuple(vec![Py::str("A"), Py::str("V")]), Py::Float(0.0)), (Py::tuple(vec![Py::str("R"), Py::str("A")]), Py::Float(-1.0)), (Py::tuple(vec![Py::str("R"), Py::str("R")]), Py::Float(5.0)), (Py::tuple(vec![Py::str("R"), Py::str("N")]), Py::Float(0.0)), (Py::tuple(vec![Py::str("R"), Py::str("D")]), Py::Float(-2.0)), (Py::tuple(vec![Py::str("R"), Py::str("C")]), Py::Float(-3.0)), (Py::tuple(vec![Py::str("R"), Py::str("Q")]), Py::Float(1.0)), (Py::tuple(vec![Py::str("R"), Py::str("E")]), Py::Float(0.0)), (Py::tuple(vec![Py::str("R"), Py::str("G")]), Py::Float(-2.0)), (Py::tuple(vec![Py::str("R"), Py::str("H")]), Py::Float(0.0)), (Py::tuple(vec![Py::str("R"), Py::str("I")]), Py::Float(-3.0)), (Py::tuple(vec![Py::str("R"), Py::str("L")]), Py::Float(-2.0)), (Py::tuple(vec![Py::str("R"), Py::str("K")]), Py::Float(2.0)), (Py::tuple(vec![Py::str("R"), Py::str("M")]), Py::Float(-1.0)), (Py::tuple(vec![Py::str("R"), Py::str("F")]), Py::Float(-3.0)), (Py::tuple(vec![Py::str("R"), Py::str("P")]), Py::Float(-2.0)), (Py::tuple(vec![Py::str("R"), Py::str("S")]), Py::Float(-1.0)), (Py::tuple(vec![Py::str("R"), Py::str("T")]), Py::Float(-1.0)), (Py::tuple(vec![Py::str("R"), Py::str("W")]), Py::Float(-3.0)), (Py::tuple(vec![Py::str("R"), Py::str("Y")]), Py::Float(-2.0)), (Py::tuple(vec![Py::str("R"), Py::str("V")]), Py::Float(-3.0)), (Py::tuple(vec![Py::str("N"), Py::str("A")]), Py::Float(-2.0)), (Py::tuple(vec![Py::str("N"), Py::str("R")]), Py::Float(0.0)), (Py::tuple(vec![Py::str("N"), Py::str("N")]), Py::Float(6.0)), (Py::tuple(vec![Py::str("N"), Py::str("D")]), Py::Float(1.0)), (Py::tuple(vec![Py::str("N"), Py::str("C")]), Py::Float(-3.0)), (Py::tuple(vec![Py::str("N"), Py::str("Q")]), Py::Float(0.0)), (Py::tuple(vec![Py::str("N"), Py::str("E")]), Py::Float(0.0)), (Py::tuple(vec![Py::str("N"), Py::str("G")]), Py::Float(0.0)), (Py::tuple(vec![Py::str("N"), Py::str("H")]), Py::Float(1.0)), (Py::tuple(vec![Py::str("N"), Py::str("I")]), Py::Float(-3.0)), (Py::tuple(vec![Py::str("N"), Py::str("L")]), Py::Float(-3.0)), (Py::tuple(vec![Py::str("N"), Py::str("K")]), Py::Float(0.0)), (Py::tuple(vec![Py::str("N"), Py::str("M")]), Py::Float(-2.0)), (Py::tuple(vec![Py::str("N"), Py::str("F")]), Py::Float(-3.0)), (Py::tuple(vec![Py::str("N"), Py::str("P")]), Py::Float(-2.0)), (Py::tuple(vec![Py::str("N"), Py::str("S")]), Py::Float(1.0)), (Py::tuple(vec![Py::str("N"), Py::str("T")]), Py::Float(0.0)), (Py::tuple(vec![Py::str("N"), Py::str("W")]), Py::Float(-4.0)), (Py::tuple(vec![Py::str("N"), Py::str("Y")]), Py::Float(-2.0)), (Py::tuple(vec![Py::str("N"), Py::str("V")]), Py::Float(-3.0)), (Py::tuple(vec![Py::str("D"), Py::str("A")]), Py::Float(-2.0)), (Py::tuple(vec![Py::str("D"), Py::str("R")]), Py::Float(-2.0)), (Py::tuple(vec![Py::str("D"), Py::str("N")]), Py::Float(1.0)), (Py::tuple(vec![Py::str("D"), Py::str("D")]), Py::Float(6.0)), (Py::tuple(vec![Py::str("D"), Py::str("C")]), Py::Float(-3.0)), (Py::tuple(vec![Py::str("D"), Py::str("Q")]), Py::Float(0.0)), (Py::tuple(vec![Py::str("D"), Py::str("E")]), Py::Float(2.0)), (Py::tuple(vec![Py::str("D"), Py::str("G")]), Py::Float(-1.0)), (Py::tuple(vec![Py::str("D"), Py::str("H")]), Py::Float(-1.0)), (Py::tuple(vec![Py::str("D"), Py::str("I")]), Py::Float(-3.0)), (Py::tuple(vec![Py::str("D"), Py::str("L")]), Py::Float(-4.0)), (Py::tuple(vec![Py::str("D"), Py::str("K")]), Py::Float(-1.0)), (Py::tuple(vec![Py::str("D"), Py::str("M")]), Py::Float(-3.0)), (Py::tuple(vec![Py::str("D"), Py::str("F")]), Py::Float(-3.0)), (Py::tuple(vec![Py::str("D"), Py::str("P")]), Py::Float(-1.0)), (Py::tuple(vec![Py::str("D"), Py::str("S")]), Py::Float(0.0)), (Py::tuple(vec![Py::str("D"), Py::str("T")]), Py::Float(-1.0)), (Py::tuple(vec![Py::str("D"), Py::str("W")]), Py::Float(-4.0)), (Py::tuple(vec![Py::str("D"), Py::str("Y")]), Py::Float(-3.0)), (Py::tuple(vec![Py::str("D"), Py::str("V")]), Py::Float(-3.0)), (Py::tuple(vec![Py::str("C"), Py::str("A")]), Py::Float(0.0)), (Py::tuple(vec![Py::str("C"), Py::str("R")]), Py::Float(-3.0)), (Py::tuple(vec![Py::str("C"), Py::str("N")]), Py::Float(-3.0)), (Py::tuple(vec![Py::str("C"), Py::str("D")]), Py::Float(-3.0)), (Py::tuple(vec![Py::str("C"), Py::str("C")]), Py::Float(9.0)), (Py::tuple(vec![Py::str("C"), Py::str("Q")]), Py::Float(-3.0)), (Py::tuple(vec![Py::str("C"), Py::str("E")]), Py::Float(-4.0)), (Py::tuple(vec![Py::str("C"), Py::str("G")]), Py::Float(-3.0)), (Py::tuple(vec![Py::str("C"), Py::str("H")]), Py::Float(-3.0)), (Py::tuple(vec![Py::str("C"), Py::str("I")]), Py::Float(-1.0)), (Py::tuple(vec![Py::str("C"), Py::str("L")]), Py::Float(-1.0)), (Py::tuple(vec![Py::str("C"), Py::str("K")]), Py::Float(-3.0)), (Py::tuple(vec![Py::str("C"), Py::str("M")]), Py::Float(-1.0)), (Py::tuple(vec![Py::str("C"), Py::str("F")]), Py::Float(-2.0)), (Py::tuple(vec![Py::str("C"), Py::str("P")]), Py::Float(-3.0)), (Py::tuple(vec![Py::str("C"), Py::str("S")]), Py::Float(-1.0)), (Py::tuple(vec![Py::str("C"), Py::str("T")]), Py::Float(-1.0)), (Py::tuple(vec![Py::str("C"), Py::str("W")]), Py::Float(-2.0)), (Py::tuple(vec![Py::str("C"), Py::str("Y")]), Py::Float(-2.0)), (Py::tuple(vec![Py::str("C"), Py::str("V")]), Py::Float(-1.0)), (Py::tuple(vec![Py::str("Q"), Py::str("A")]), Py::Float(-1.0)), (Py::tuple(vec![Py::str("Q"), Py::str("R")]), Py::Float(1.0)), (Py::tuple(vec![Py::str("Q"), Py::str("N")]), Py::Float(0.0)), (Py::tuple(vec![Py::str("Q"), Py::str("D")]), Py::Float(0.0)), (Py::tuple(vec![Py::str("Q"), Py::str("C")]), Py::Float(-3.0)), (Py::tuple(vec![Py::str("Q"), Py::str("Q")]), Py::Float(5.0)), (Py::tuple(vec![Py::str("Q"), Py::str("E")]), Py::Float(2.0)), (Py::tuple(vec![Py::str("Q"), Py::str("G")]), Py::Float(-2.0)), (Py::tuple(vec![Py::str("Q"), Py::str("H")]), Py::Float(0.0)), (Py::tuple(vec![Py::str("Q"), Py::str("I")]), Py::Float(-3.0)), (Py::tuple(vec![Py::str("Q"), Py::str("L")]), Py::Float(-2.0)), (Py::tuple(vec![Py::str("Q"), Py::str("K")]), Py::Float(1.0)), (Py::tuple(vec![Py::str("Q"), Py::str("M")]), Py::Float(0.0)), (Py::tuple(vec![Py::str("Q"), Py::str("F")]), Py::Float(-3.0)), (Py::tuple(vec![Py::str("Q"), Py::str("P")]), Py::Float(-1.0)), (Py::tuple(vec![Py::str("Q"), Py::str("S")]), Py::Float(0.0)), (Py::tuple(vec![Py::str("Q"), Py::str("T")]), Py::Float(-1.0)), (Py::tuple(vec![Py::str("Q"), Py::str("W")]), Py::Float(-2.0)), (Py::tuple(vec![Py::str("Q"), Py::str("Y")]), Py::Float(-1.0)), (Py::tuple(vec![Py::str("Q"), Py::str("V")]), Py::Float(-2.0)), (Py::tuple(vec![Py::str("E"), Py::str("A")]), Py::Float(-1.0)), (Py::tuple(vec![Py::str("E"), Py::str("R")]), Py::Float(0.0)), (Py::tuple(vec![Py::str("E"), Py::str("N")]), Py::Float(0.0)), (Py::tuple(vec![Py::str("E"), Py::str("D")]), Py::Float(2.0)), (Py::tuple(vec![Py::str("E"), Py::str("C")]), Py::Float(-4.0)), (Py::tuple(vec![Py::str("E"), Py::str("Q")]), Py::Float(2.0)), (Py::tuple(vec![Py::str("E"), Py::str("E")]), Py::Float(5.0)), (Py::tuple(vec![Py::str("E"), Py::str("G")]), Py::Float(-2.0)), (Py::tuple(vec![Py::str("E"), Py::str("H")]), Py::Float(0.0)), (Py::tuple(vec![Py::str("E"), Py::str("I")]), Py::Float(-3.0)), (Py::tuple(vec![Py::str("E"), Py::str("L")]), Py::Float(-3.0)), (Py::tuple(vec![Py::str("E"), Py::str("K")]), Py::Float(1.0)), (Py::tuple(vec![Py::str("E"), Py::str("M")]), Py::Float(-2.0)), (Py::tuple(vec![Py::str("E"), Py::str("F")]), Py::Float(-3.0)), (Py::tuple(vec![Py::str("E"), Py::str("P")]), Py::Float(-1.0)), (Py::tuple(vec![Py::str("E"), Py::str("S")]), Py::Float(0.0)), (Py::tuple(vec![Py::str("E"), Py::str("T")]), Py::Float(-1.0)), (Py::tuple(vec![Py::str("E"), Py::str("W")]), Py::Float(-3.0)), (Py::tuple(vec![Py::str("E"), Py::str("Y")]), Py::Float(-2.0)), (Py::tuple(vec![Py::str("E"), Py::str("V")]), Py::Float(-2.0)), (Py::tuple(vec![Py::str("G"), Py::str("A")]), Py::Float(0.0)), (Py::tuple(vec![Py::str("G"), Py::str("R")]), Py::Float(-2.0)), (Py::tuple(vec![Py::str("G"), Py::str("N")]), Py::Float(0.0)), (Py::tuple(vec![Py::str("G"), Py::str("D")]), Py::Float(-1.0)), (Py::tuple(vec![Py::str("G"), Py::str("C")]), Py::Float(-3.0)), (Py::tuple(vec![Py::str("G"), Py::str("Q")]), Py::Float(-2.0)), (Py::tuple(vec![Py::str("G"), Py::str("E")]), Py::Float(-2.0)), (Py::tuple(vec![Py::str("G"), Py::str("G")]), Py::Float(6.0)), (Py::tuple(vec![Py::str("G"), Py::str("H")]), Py::Float(-2.0)), (Py::tuple(vec![Py::str("G"), Py::str("I")]), Py::Float(-4.0)), (Py::tuple(vec![Py::str("G"), Py::str("L")]), Py::Float(-4.0)), (Py::tuple(vec![Py::str("G"), Py::str("K")]), Py::Float(-2.0)), (Py::tuple(vec![Py::str("G"), Py::str("M")]), Py::Float(-3.0)), (Py::tuple(vec![Py::str("G"), Py::str("F")]), Py::Float(-3.0)), (Py::tuple(vec![Py::str("G"), Py::str("P")]), Py::Float(-2.0)), (Py::tuple(vec![Py::str("G"), Py::str("S")]), Py::Float(0.0)), (Py::tuple(vec![Py::str("G"), Py::str("T")]), Py::Float(-2.0)), (Py::tuple(vec![Py::str("G"), Py::str("W")]), Py::Float(-2.0)), (Py::tuple(vec![Py::str("G"), Py::str("Y")]), Py::Float(-3.0)), (Py::tuple(vec![Py::str("G"), Py::str("V")]), Py::Float(-3.0)), (Py::tuple(vec![Py::str("H"), Py::str("A")]), Py::Float(-2.0)), (Py::tuple(vec![Py::str("H"), Py::str("R")]), Py::Float(0.0)), (Py::tuple(vec![Py::str("H"), Py::str("N")]), Py::Float(1.0)), (Py::tuple(vec![Py::str("H"), Py::str("D")]), Py::Float(-1.0)), (Py::tuple(vec![Py::str("H"), Py::str("C")]), Py::Float(-3.0)), (Py::tuple(vec![Py::str("H"), Py::str("Q")]), Py::Float(0.0)), (Py::tuple(vec![Py::str("H"), Py::str("E")]), Py::Float(0.0)), (Py::tuple(vec![Py::str("H"), Py::str("G")]), Py::Float(-2.0)), (Py::tuple(vec![Py::str("H"), Py::str("H")]), Py::Float(8.0)), (Py::tuple(vec![Py::str("H"), Py::str("I")]), Py::Float(-3.0)), (Py::tuple(vec![Py::str("H"), Py::str("L")]), Py::Float(-3.0)), (Py::tuple(vec![Py::str("H"), Py::str("K")]), Py::Float(-1.0)), (Py::tuple(vec![Py::str("H"), Py::str("M")]), Py::Float(-2.0)), (Py::tuple(vec![Py::str("H"), Py::str("F")]), Py::Float(-1.0)), (Py::tuple(vec![Py::str("H"), Py::str("P")]), Py::Float(-2.0)), (Py::tuple(vec![Py::str("H"), Py::str("S")]), Py::Float(-1.0)), (Py::tuple(vec![Py::str("H"), Py::str("T")]), Py::Float(-2.0)), (Py::tuple(vec![Py::str("H"), Py::str("W")]), Py::Float(-2.0)), (Py::tuple(vec![Py::str("H"), Py::str("Y")]), Py::Float(2.0)), (Py::tuple(vec![Py::str("H"), Py::str("V")]), Py::Float(-3.0)), (Py::tuple(vec![Py::str("I"), Py::str("A")]), Py::Float(-1.0)), (Py::tuple(vec![Py::str("I"), Py::str("R")]), Py::Float(-3.0)), (Py::tuple(vec![Py::str("I"), Py::str("N")]), Py::Float(-3.0)), (Py::tuple(vec![Py::str("I"), Py::str("D")]), Py::Float(-3.0)), (Py::tuple(vec![Py::str("I"), Py::str("C")]), Py::Float(-1.0)), (Py::tuple(vec![Py::str("I"), Py::str("Q")]), Py::Float(-3.0)), (Py::tuple(vec![Py::str("I"), Py::str("E")]), Py::Float(-3.0)), (Py::tuple(vec![Py::str("I"), Py::str("G")]), Py::Float(-4.0)), (Py::tuple(vec![Py::str("I"), Py::str("H")]), Py::Float(-3.0)), (Py::tuple(vec![Py::str("I"), Py::str("I")]), Py::Float(4.0)), (Py::tuple(vec![Py::str("I"), Py::str("L")]), Py::Float(2.0)), (Py::tuple(vec![Py::str("I"), Py::str("K")]), Py::Float(-3.0)), (Py::tuple(vec![Py::str("I"), Py::str("M")]), Py::Float(1.0)), (Py::tuple(vec![Py::str("I"), Py::str("F")]), Py::Float(0.0)), (Py::tuple(vec![Py::str("I"), Py::str("P")]), Py::Float(-3.0)), (Py::tuple(vec![Py::str("I"), Py::str("S")]), Py::Float(-2.0)), (Py::tuple(vec![Py::str("I"), Py::str("T")]), Py::Float(-1.0)), (Py::tuple(vec![Py::str("I"), Py::str("W")]), Py::Float(-3.0)), (Py::tuple(vec![Py::str("I"), Py::str("Y")]), Py::Float(-1.0)), (Py::tuple(vec![Py::str("I"), Py::str("V")]), Py::Float(3.0)), (Py::tuple(vec![Py::str("L"), Py::str("A")]), Py::Float(-1.0)), (Py::tuple(vec![Py::str("L"), Py::str("R")]), Py::Float(-2.0)), (Py::tuple(vec![Py::str("L"), Py::str("N")]), Py::Float(-3.0)), (Py::tuple(vec![Py::str("L"), Py::str("D")]), Py::Float(-4.0)), (Py::tuple(vec![Py::str("L"), Py::str("C")]), Py::Float(-1.0)), (Py::tuple(vec![Py::str("L"), Py::str("Q")]), Py::Float(-2.0)), (Py::tuple(vec![Py::str("L"), Py::str("E")]), Py::Float(-3.0)), (Py::tuple(vec![Py::str("L"), Py::str("G")]), Py::Float(-4.0)), (Py::tuple(vec![Py::str("L"), Py::str("H")]), Py::Float(-3.0)), (Py::tuple(vec![Py::str("L"), Py::str("I")]), Py::Float(2.0)), (Py::tuple(vec![Py::str("L"), Py::str("L")]), Py::Float(4.0)), (Py::tuple(vec![Py::str("L"), Py::str("K")]), Py::Float(-2.0)), (Py::tuple(vec![Py::str("L"), Py::str("M")]), Py::Float(2.0)), (Py::tuple(vec![Py::str("L"), Py::str("F")]), Py::Float(0.0)), (Py::tuple(vec![Py::str("L"), Py::str("P")]), Py::Float(-3.0)), (Py::tuple(vec![Py::str("L"), Py::str("S")]), Py::Float(-2.0)), (Py::tuple(vec![Py::str("L"), Py::str("T")]), Py::Float(-1.0)), (Py::tuple(vec![Py::str("L"), Py::str("W")]), Py::Float(-2.0)), (Py::tuple(vec![Py::str("L"), Py::str("Y")]), Py::Float(-1.0)), (Py::tuple(vec![Py::str("L"), Py::str("V")]), Py::Float(1.0)), (Py::tuple(vec![Py::str("K"), Py::str("A")]), Py::Float(-1.0)), (Py::tuple(vec![Py::str("K"), Py::str("R")]), Py::Float(2.0)), (Py::tuple(vec![Py::str("K"), Py::str("N")]), Py::Float(0.0)), (Py::tuple(vec![Py::str("K"), Py::str("D")]), Py::Float(-1.0)), (Py::tuple(vec![Py::str("K"), Py::str("C")]), Py::Float(-3.0)), (Py::tuple(vec![Py::str("K"), Py::str("Q")]), Py::Float(1.0)), (Py::tuple(vec![Py::str("K"), Py::str("E")]), Py::Float(1.0)), (Py::tuple(vec![Py::str("K"), Py::str("G")]), Py::Float(-2.0)), (Py::tuple(vec![Py::str("K"), Py::str("H")]), Py::Float(-1.0)), (Py::tuple(vec![Py::str("K"), Py::str("I")]), Py::Float(-3.0)), (Py::tuple(vec![Py::str("K"), Py::str("L")]), Py::Float(-2.0)), (Py::tuple(vec![Py::str("K"), Py::str("K")]), Py::Float(5.0)), (Py::tuple(vec![Py::str("K"), Py::str("M")]), Py::Float(-1.0)), (Py::tuple(vec![Py::str("K"), Py::str("F")]), Py::Float(-3.0)), (Py::tuple(vec![Py::str("K"), Py::str("P")]), Py::Float(-1.0)), (Py::tuple(vec![Py::str("K"), Py::str("S")]), Py::Float(0.0)), (Py::tuple(vec![Py::str("K"), Py::str("T")]), Py::Float(-1.0)), (Py::tuple(vec![Py::str("K"), Py::str("W")]), Py::Float(-3.0)), (Py::tuple(vec![Py::str("K"), Py::str("Y")]), Py::Float(-2.0)), (Py::tuple(vec![Py::str("K"), Py::str("V")]), Py::Float(-2.0)), (Py::tuple(vec![Py::str("M"), Py::str("A")]), Py::Float(-1.0)), (Py::tuple(vec![Py::str("M"), Py::str("R")]), Py::Float(-1.0)), (Py::tuple(vec![Py::str("M"), Py::str("N")]), Py::Float(-2.0)), (Py::tuple(vec![Py::str("M"), Py::str("D")]), Py::Float(-3.0)), (Py::tuple(vec![Py::str("M"), Py::str("C")]), Py::Float(-1.0)), (Py::tuple(vec![Py::str("M"), Py::str("Q")]), Py::Float(0.0)), (Py::tuple(vec![Py::str("M"), Py::str("E")]), Py::Float(-2.0)), (Py::tuple(vec![Py::str("M"), Py::str("G")]), Py::Float(-3.0)), (Py::tuple(vec![Py::str("M"), Py::str("H")]), Py::Float(-2.0)), (Py::tuple(vec![Py::str("M"), Py::str("I")]), Py::Float(1.0)), (Py::tuple(vec![Py::str("M"), Py::str("L")]), Py::Float(2.0)), (Py::tuple(vec![Py::str("M"), Py::str("K")]), Py::Float(-1.0)), (Py::tuple(vec![Py::str("M"), Py::str("M")]), Py::Float(5.0)), (Py::tuple(vec![Py::str("M"), Py::str("F")]), Py::Float(0.0)), (Py::tuple(vec![Py::str("M"), Py::str("P")]), Py::Float(-2.0)), (Py::tuple(vec![Py::str("M"), Py::str("S")]), Py::Float(-1.0)), (Py::tuple(vec![Py::str("M"), Py::str("T")]), Py::Float(-1.0)), (Py::tuple(vec![Py::str("M"), Py::str("W")]), Py::Float(-1.0)), (Py::tuple(vec![Py::str("M"), Py::str("Y")]), Py::Float(-1.0)), (Py::tuple(vec![Py::str("M"), Py::str("V")]), Py::Float(1.0)), (Py::tuple(vec![Py::str("F"), Py::str("A")]), Py::Float(-2.0)), (Py::tuple(vec![Py::str("F"), Py::str("R")]), Py::Float(-3.0)), (Py::tuple(vec![Py::str("F"), Py::str("N")]), Py::Float(-3.0)), (Py::tuple(vec![Py::str("F"), Py::str("D")]), Py::Float(-3.0)), (Py::tuple(vec![Py::str("F"), Py::str("C")]), Py::Float(-2.0)), (Py::tuple(vec![Py::str("F"), Py::str("Q")]), Py::Float(-3.0)), (Py::tuple(vec![Py::str("F"), Py::str("E")]), Py::Float(-3.0)), (Py::tuple(vec![Py::str("F"), Py::str("G")]), Py::Float(-3.0)), (Py::tuple(vec![Py::str("F"), Py::str("H")]), Py::Float(-1.0)), (Py::tuple(vec![Py::str("F"), Py::str("I")]), Py::Float(0.0)), (Py::tuple(vec![Py::str("F"), Py::str("L")]), Py::Float(0.0)), (Py::tuple(vec![Py::str("F"), Py::str("K")]), Py::Float(-3.0)), (Py::tuple(vec![Py::str("F"), Py::str("M")]), Py::Float(0.0)), (Py::tuple(vec![Py::str("F"), Py::str("F")]), Py::Float(6.0)), (Py::tuple(vec![Py::str("F"), Py::str("P")]), Py::Float(-4.0)), (Py::tuple(vec![Py::str("F"), Py::str("S")]), Py::Float(-2.0)), (Py::tuple(vec![Py::str("F"), Py::str("T")]), Py::Float(-2.0)), (Py::tuple(vec![Py::str("F"), Py::str("W")]), Py::Float(1.0)), (Py::tuple(vec![Py::str("F"), Py::str("Y")]), Py::Float(3.0)), (Py::tuple(vec![Py::str("F"), Py::str("V")]), Py::Float(-1.0)), (Py::tuple(vec![Py::str("P"), Py::str("A")]), Py::Float(-1.0)), (Py::tuple(vec![Py::str("P"), Py::str("R")]), Py::Float(-2.0)), (Py::tuple(vec![Py::str("P"), Py::str("N")]), Py::Float(-2.0)), (Py::tuple(vec![Py::str("P"), Py::str("D")]), Py::Float(-1.0)), (Py::tuple(vec![Py::str("P"), Py::str("C")]), Py::Float(-3.0)), (Py::tuple(vec![Py::str("P"), Py::str("Q")]), Py::Float(-1.0)), (Py::tuple(vec![Py::str("P"), Py::str("E")]), Py::Float(-1.0)), (Py::tuple(vec![Py::str("P"), Py::str("G")]), Py::Float(-2.0)), (Py::tuple(vec![Py::str("P"), Py::str("H")]), Py::Float(-2.0)), (Py::tuple(vec![Py::str("P"), Py::str("I")]), Py::Float(-3.0)), (Py::tuple(vec![Py::str("P"), Py::str("L")]), Py::Float(-3.0)), (Py::tuple(vec![Py::str("P"), Py::str("K")]), Py::Float(-1.0)), (Py::tuple(vec![Py::str("P"), Py::str("M")]), Py::Float(-2.0)), (Py::tuple(vec![Py::str("P"), Py::str("F")]), Py::Float(-4.0)), (Py::tuple(vec![Py::str("P"), Py::str("P")]), Py::Float(7.0)), (Py::tuple(vec![Py::str("P"), Py::str("S")]), Py::Float(-1.0)), (Py::tuple(vec![Py::str("P"), Py::str("T")]), Py::Float(-1.0)), (Py::tuple(vec![Py::str("P"), Py::str("W")]), Py::Float(-4.0)), (Py::tuple(vec![Py::str("P"), Py::str("Y")]), Py::Float(-3.0)), (Py::tuple(vec![Py::str("P"), Py::str("V")]), Py::Float(-2.0)), (Py::tuple(vec![Py::str("S"), Py::str("A")]), Py::Float(1.0)), (Py::tuple(vec![Py::str("S"), Py::str("R")]), Py::Float(-1.0)), (Py::tuple(vec![Py::str("S"), Py::str("N")]), Py::Float(1.0)), (Py::tuple(vec![Py::str("S"), Py::str("D")]), Py::Float(0.0)), (Py::tuple(vec![Py::str("S"), Py::str("C")]), Py::Float(-1.0)), (Py::tuple(vec![Py::str("S"), Py::str("Q")]), Py::Float(0.0)), (Py::tuple(vec![Py::str("S"), Py::str("E")]), Py::Float(0.0)), (Py::tuple(vec![Py::str("S"), Py::str("G")]), Py::Float(0.0)), (Py::tuple(vec![Py::str("S"), Py::str("H")]), Py::Float(-1.0)), (Py::tuple(vec![Py::str("S"), Py::str("I")]), Py::Float(-2.0)), (Py::tuple(vec![Py::str("S"), Py::str("L")]), Py::Float(-2.0)), (Py::tuple(vec![Py::str("S"), Py::str("K")]), Py::Float(0.0)), (Py::tuple(vec![Py::str("S"), Py::str("M")]), Py::Float(-1.0)), (Py::tuple(vec![Py::str("S"), Py::str("F")]), Py::Float(-2.0)), (Py::tuple(vec![Py::str("S"), Py::str("P")]), Py::Float(-1.0)), (Py::tuple(vec![Py::str("S"), Py::str("S")]), Py::Float(4.0)), (Py::tuple(vec![Py::str("S"), Py::str("T")]), Py::Float(1.0)), (Py::tuple(vec![Py::str("S"), Py::str("W")]), Py::Float(-3.0)), (Py::tuple(vec![Py::str("S"), Py::str("Y")]), Py::Float(-2.0)), (Py::tuple(vec![Py::str("S"), Py::str("V")]), Py::Float(-2.0)), (Py::tuple(vec![Py::str("T"), Py::str("A")]), Py::Float(0.0)), (Py::tuple(vec![Py::str("T"), Py::str("R")]), Py::Float(-1.0)), (Py::tuple(vec![Py::str("T"), Py::str("N")]), Py::Float(0.0)), (Py::tuple(vec![Py::str("T"), Py::str("D")]), Py::Float(-1.0)), (Py::tuple(vec![Py::str("T"), Py::str("C")]), Py::Float(-1.0)), (Py::tuple(vec![Py::str("T"), Py::str("Q")]), Py::Float(-1.0)), (Py::tuple(vec![Py::str("T"), Py::str("E")]), Py::Float(-1.0)), (Py::tuple(vec![Py::str("T"), Py::str("G")]), Py::Float(-2.0)), (Py::tuple(vec![Py::str("T"), Py::str("H")]), Py::Float(-2.0)), (Py::tuple(vec![Py::str("T"), Py::str("I")]), Py::Float(-1.0)), (Py::tuple(vec![Py::str("T"), Py::str("L")]), Py::Float(-1.0)), (Py::tuple(vec![Py::str("T"), Py::str("K")]), Py::Float(-1.0)), (Py::tuple(vec![Py::str("T"), Py::str("M")]), Py::Float(-1.0)), (Py::tuple(vec![Py::str("T"), Py::str("F")]), Py::Float(-2.0)), (Py::tuple(vec![Py::str("T"), Py::str("P")]), Py::Float(-1.0)), (Py::tuple(vec![Py::str("T"), Py::str("S")]), Py::Float(1.0)), (Py::tuple(vec![Py::str("T"), Py::str("T")]), Py::Float(5.0)), (Py::tuple(vec![Py::str("T"), Py::str("W")]), Py::Float(-2.0)), (Py::tuple(vec![Py::str("T"), Py::str("Y")]), Py::Float(-2.0)), (Py::tuple(vec![Py::str("T"), Py::str("V")]), Py::Float(0.0)), (Py::tuple(vec![Py::str("W"), Py::str("A")]), Py::Float(-3.0)), (Py::tuple(vec![Py::str("W"), Py::str("R")]), Py::Float(-3.0)), (Py::tuple(vec![Py::str("W"), Py::str("N")]), Py::Float(-4.0)), (Py::tuple(vec![Py::str("W"), Py::str("D")]), Py::Float(-4.0)), (Py::tuple(vec![Py::str("W"), Py::str("C")]), Py::Float(-2.0)), (Py::tuple(vec![Py::str("W"), Py::str("Q")]), Py::Float(-2.0)), (Py::tuple(vec![Py::str("W"), Py::str("E")]), Py::Float(-3.0)), (Py::tuple(vec![Py::str("W"), Py::str("G")]), Py::Float(-2.0)), (Py::tuple(vec![Py::str("W"), Py::str("H")]), Py::Float(-2.0)), (Py::tuple(vec![Py::str("W"), Py::str("I")]), Py::Float(-3.0)), (Py::tuple(vec![Py::str("W"), Py::str("L")]), Py::Float(-2.0)), (Py::tuple(vec![Py::str("W"), Py::str("K")]), Py::Float(-3.0)), (Py::tuple(vec![Py::str("W"), Py::str("M")]), Py::Float(-1.0)), (Py::tuple(vec![Py::str("W"), Py::str("F")]), Py::Float(1.0)), (Py::tuple(vec![Py::str("W"), Py::str("P")]), Py::Float(-4.0)), (Py::tuple(vec![Py::str("W"), Py::str("S")]), Py::Float(-3.0)), (Py::tuple(vec![Py::str("W"), Py::str("T")]), Py::Float(-2.0)), (Py::tuple(vec![Py::str("W"), Py::str("W")]), Py::Float(11.0)), (Py::tuple(vec![Py::str("W"), Py::str("Y")]), Py::Float(2.0)), (Py::tuple(vec![Py::str("W"), Py::str("V")]), Py::Float(-3.0)), (Py::tuple(vec![Py::str("Y"), Py::str("A")]), Py::Float(-2.0)), (Py::tuple(vec![Py::str("Y"), Py::str("R")]), Py::Float(-2.0)), (Py::tuple(vec![Py::str("Y"), Py::str("N")]), Py::Float(-2.0)), (Py::tuple(vec![Py::str("Y"), Py::str("D")]), Py::Float(-3.0)), (Py::tuple(vec![Py::str("Y"), Py::str("C")]), Py::Float(-2.0)), (Py::tuple(vec![Py::str("Y"), Py::str("Q")]), Py::Float(-1.0)), (Py::tuple(vec![Py::str("Y"), Py::str("E")]), Py::Float(-2.0)), (Py::tuple(vec![Py::str("Y"), Py::str("G")]), Py::Float(-3.0)), (Py::tuple(vec![Py::str("Y"), Py::str("H")]), Py::Float(2.0)), (Py::tuple(vec![Py::str("Y"), Py::str("I")]), Py::Float(-1.0)), (Py::tuple(vec![Py::str("Y"), Py::str("L")]), Py::Float(-1.0)), (Py::tuple(vec![Py::str("Y"), Py::str("K")]), Py::Float(-2.0)), (Py::tuple(vec![Py::str("Y"), Py::str("M")]), Py::Float(-1.0)), (Py::tuple(vec![Py::str("Y"), Py::str("F")]), Py::Float(3.0)), (Py::tuple(vec![Py::str("Y"), Py::str("P")]), Py::Float(-3.0)), (Py::tuple(vec![Py::str("Y"), Py::str("S")]), Py::Float(-2.0)), (Py::tuple(vec![Py::str("Y"), Py::str("T")]), Py::Float(-2.0)), (Py::tuple(vec![Py::str("Y"), Py::str("W")]), Py::Float(2.0)), (Py::tuple(vec![Py::str("Y"), Py::str("Y")]), Py::Float(7.0)), (Py::tuple(vec![Py::str("Y"), Py::str("V")]), Py::Float(-1.0)), (Py::tuple(vec![Py::str("V"), Py::str("A")]), Py::Float(0.0)), (Py::tuple(vec![Py::str("V"), Py::str("R")]), Py::Float(-3.0)), (Py::tuple(vec![Py::str("V"), Py::str("N")]), Py::Float(-3.0)), (Py::tuple(vec![Py::str("V"), Py::str("D")]), Py::Float(-3.0)), (Py::tuple(vec![Py::str("V"), Py::str("C")]), Py::Float(-1.0)), (Py::tuple(vec![Py::str("V"), Py::str("Q")]), Py::Float(-2.0)), (Py::tuple(vec![Py::str("V"), Py::str("E")]), Py::Float(-2.0)), (Py::tuple(vec![Py::str("V"), Py::str("G")]), Py::Float(-3.0)), (Py::tuple(vec![Py::str("V"), Py::str("H")]), Py::Float(-3.0)), (Py::tuple(vec![Py::str("V"), Py::str("I")]), Py::Float(3.0)), (Py::tuple(vec![Py::str("V"), Py::str("L")]), Py::Float(1.0)), (Py::tuple(vec![Py::str("V"), Py::str("K")]), Py::Float(-2.0)), (Py::tuple(vec![Py::str("V"), Py::str("M")]), Py::Float(1.0)), (Py::tuple(vec![Py::str("V"), Py::str("F")]), Py::Float(-1.0)), (Py::tuple(vec![Py::str("V"), Py::str("P")]), Py::Float(-2.0)), (Py::tuple(vec![Py::str("V"), Py::str("S")]), Py::Float(-2.0)), (Py::tuple(vec![Py::str("V"), Py::str("T")]), Py::Float(0.0)), (Py::tuple(vec![Py::str("V"), Py::str("W")]), Py::Float(-3.0)), (Py::tuple(vec![Py::str("V"), Py::str("Y")]), Py::Float(-1.0)), (Py::tuple(vec![Py::str("V"), Py::str("V")]), Py::Float(4.0))])) })());
    V.with(|v| v.clone())
}

#[allow(non_snake_case)]
fn c_CODON() -> R {
    thread_local!(static V: R = (|| -> R { Ok(Py::dict(vec![(Py::str("TTT"), Py::str("F")), (Py::str("TTC"), Py::str("F")), (Py::str("TTA"), Py::str("L")), (Py::str("TTG"), Py::str("L")), (Py::str("TCT"), Py::str("S")), (Py::str("TCC"), Py::str("S")), (Py::str("TCA"), Py::str("S")), (Py::str("TCG"), Py::str("S")), (Py::str("TAT"), Py::str("Y")), (Py::str("TAC"), Py::str("Y")), (Py::str("TAA"), Py::str("*")), (Py::str("TAG"), Py::str("*")), (Py::str("TGT"), Py::str("C")), (Py::str("TGC"), Py::str("C")), (Py::str("TGA"), Py::str("*")), (Py::str("TGG"), Py::str("W")), (Py::str("CTT"), Py::str("L")), (Py::str("CTC"), Py::str("L")), (Py::str("CTA"), Py::str("L")), (Py::str("CTG"), Py::str("L")), (Py::str("CCT"), Py::str("P")), (Py::str("CCC"), Py::str("P")), (Py::str("CCA"), Py::str("P")), (Py::str("CCG"), Py::str("P")), (Py::str("CAT"), Py::str("H")), (Py::str("CAC"), Py::str("H")), (Py::str("CAA"), Py::str("Q")), (Py::str("CAG"), Py::str("Q")), (Py::str("CGT"), Py::str("R")), (Py::str("CGC"), Py::str("R")), (Py::str("CGA"), Py::str("R")), (Py::str("CGG"), Py::str("R")), (Py::str("ATT"), Py::str("I")), (Py::str("ATC"), Py::str("I")), (Py::str("ATA"), Py::str("I")), (Py::str("ATG"), Py::str("M")), (Py::str("ACT"), Py::str("T")), (Py::str("ACC"), Py::str("T")), (Py::str("ACA"), Py::str("T")), (Py::str("ACG"), Py::str("T")), (Py::str("AAT"), Py::str("N")), (Py::str("AAC"), Py::str("N")), (Py::str("AAA"), Py::str("K")), (Py::str("AAG"), Py::str("K")), (Py::str("AGT"), Py::str("S")), (Py::str("AGC"), Py::str("S")), (Py::str("AGA"), Py::str("R")), (Py::str("AGG"), Py::str("R")), (Py::str("GTT"), Py::str("V")), (Py::str("GTC"), Py::str("V")), (Py::str("GTA"), Py::str("V")), (Py::str("GTG"), Py::str("V")), (Py::str("GCT"), Py::str("A")), (Py::str("GCC"), Py::str("A")), (Py::str("GCA"), Py::str("A")), (Py::str("GCG"), Py::str("A")), (Py::str("GAT"), Py::str("D")), (Py::str("GAC"), Py::str("D")), (Py::str("GAA"), Py::str("E")), (Py::str("GAG"), Py::str("E")), (Py::str("GGT"), Py::str("G")), (Py::str("GGC"), Py::str("G")), (Py::str("GGA"), Py::str("G")), (Py::str("GGG"), Py::str("G"))])) })());
    V.with(|v| v.clone())
}

#[allow(non_snake_case)]
fn c_DIWV() -> R {
    thread_local!(static V: R = (|| -> R { Ok(Py::dict(vec![(Py::str("A"), Py::dict(vec![(Py::str("A"), Py::Float(1.0)), (Py::str("C"), Py::Float(44.94)), (Py::str("E"), Py::Float(1.0)), (Py::str("D"), Py::Float(-7.49)), (Py::str("G"), Py::Float(1.0)), (Py::str("F"), Py::Float(1.0)), (Py::str("I"), Py::Float(1.0)), (Py::str("H"), Py::Float(-7.49)), (Py::str("K"), Py::Float(1.0)), (Py::str("M"), Py::Float(1.0)), (Py::str("L"), Py::Float(1.0)), (Py::str("N"), Py::Float(1.0)), (Py::str("Q"), Py::Float(1.0)), (Py::str("P"), Py::Float(20.26)), (Py::str("S"), Py::Float(1.0)), (Py::str("R"), Py::Float(1.0)), (Py::str("T"), Py::Float(1.0)), (Py::str("W"), Py::Float(1.0)), (Py::str("V"), Py::Float(1.0)), (Py::str("Y"), Py::Float(1.0))])), (Py::str("C"), Py::dict(vec![(Py::str("A"), Py::Float(1.0)), (Py::str("C"), Py::Float(1.0)), (Py::str("E"), Py::Float(1.0)), (Py::str("D"), Py::Float(20.26)), (Py::str("G"), Py::Float(1.0)), (Py::str("F"), Py::Float(1.0)), (Py::str("I"), Py::Float(1.0)), (Py::str("H"), Py::Float(33.6)), (Py::str("K"), Py::Float(1.0)), (Py::str("M"), Py::Float(33.6)), (Py::str("L"), Py::Float(20.26)), (Py::str("N"), Py::Float(1.0)), (Py::str("Q"), Py::Float(-6.54)), (Py::str("P"), Py::Float(20.26)), (Py::str("S"), Py::Float(1.0)), (Py::str("R"), Py::Float(1.0)), (Py::str("T"), Py::Float(33.6)), (Py::str("W"), Py::Float(24.68)), (Py::str("V"), Py::Float(-6.54)), (Py::str("Y"), Py::Float(1.0))])), (Py::str("E"), Py::dict(vec![(Py::str("A"), Py::Float(1.0)), (Py::str("C"), Py::Float(44.94)), (Py::str("E"), Py::Float(33.6)), (Py::str("D"), Py::Float(20.26)), (Py::str("G"), Py::Float(1.0)), (Py::str("F"), Py::Float(1.0)), (Py::str("I"), Py::Float(20.26)), (Py::str("H"), Py::Float(-6.54)), (Py::str("K"), Py::Float(1.0)), (Py::str("M"), Py::Float(1.0)), (Py::str("L"), Py::Float(1.0)), (Py::str("N"), Py::Float(1.0)), (Py::str("Q"), Py::Float(20.26)), (Py::str("P"), Py::Float(20.26)), (Py::str("S"), Py::Float(20.26)), (Py::str("R"), Py::Float(1.0)), (Py::str("T"), Py::Float(1.0)), (Py::str("W"), Py::Float(-14.03)), (Py::str("V"), Py::Float(1.0)), (Py::str("Y"), Py::Float(1.0))])), (Py::str("D"), Py::dict(vec![(Py::str("A"), Py::Float(1.0)), (Py::str("C"), Py::Float(1.0)), (Py::str("E"), Py::Float(1.0)), (Py::str("D"), Py::Float(1.0)), (Py::str("G"), Py::Float(1.0)), (Py::str("F"), Py::Float(-6.54)), (Py::str("I"), Py::Float(1.0)), (Py::str("H"), Py::Float(1.0)), (Py::str("K"), Py::Float(-7.49)), (Py::str("M"), Py::Float(1.0)), (Py::str("L"), Py::Float(1.0)), (Py::str("N"), Py::Float(1.0)), (Py::str("Q"), Py::Float(1.0)), (Py::str("P"), Py::Float(1.0)), (Py::str("S"), Py::Float(20.26)), (Py::str("R"), Py::Float(-6.54)), (Py::str("T"), Py::Float(-14.03)), (Py::str("W"), Py::Float(1.0)), (Py::str("V"), Py::Float(1.0)), (Py::str("Y"), Py::Float(1.0))])), (Py::str("G"), Py::dict(vec![(Py::str("A"), Py::Float(-7.49)), (Py::str("C"), Py::Float(1.0)), (Py::str("E"), Py::Float(-6.54)), (Py::str("D"), Py::Float(1.0)), (Py::str("G"), Py::Float(13.34)), (Py::str("F"), Py::Float(1.0)), (Py::str("I"), Py::Float(-7.49)), (Py::str("H"), Py::Float(1.0)), (Py::str("K"), Py::Float(-7.49)), (Py::str("M"), Py::Float(1.0)), (Py::str("L"), Py::Float(1.0)), (Py::str("N"), Py::Float(-7.49)), (Py::str("Q"), Py::Float(1.0)), (Py::str("P"), Py::Float(1.0)), (Py::str("S"), Py::Float(1.0)), (Py::str("R"), Py::Float(1.0)), (Py::str("T"), Py::Float(-7.49)), (Py::str("W"), Py::Float(13.34)), (Py::str("V"), Py::Float(1.0)), (Py::str("Y"), Py::Float(-7.49))])), (Py::str("F"), Py::dict(vec![(Py::str("A"), Py::Float(1.0)), (Py::str("C"), Py::Float(1.0)), (Py::str("E"), Py::Float(1.0)), (Py::str("D"), Py::Float(13.34)), (Py::str("G"), Py::Float(1.0)), (Py::str("F"), Py::Float(1.0)), (Py::str("I"), Py::Float(1.0)), (Py::str("H"), Py::Float(1.0)), (Py::str("K"), Py::Float(-14.03)), (Py::str("M"), Py::Float(1.0)), (Py::str("L"), Py::Float(1.0)), (Py::str("N"), Py::Float(1.0)), (Py::str("Q"), Py::Float(1.0)), (Py::str("P"), Py::Float(20.26)), (Py::str("S"), Py::Float(1.0)), (Py::str("R"), Py::Float(1.0)), (Py::str("T"), Py::Float(1.0)), (Py::str("W"), Py::Float(1.0)), (Py::str("V"), Py::Float(1.0)), (Py::str("Y"), Py::Float(33.601))])), (Py::str("I"), Py::dict(vec![(Py::str("A"), Py::Float(1.0)), (Py::str("C"), Py::Float(1.0)), (Py::str("E"), Py::Float(44.94)), (Py::str("D"), Py::Float(1.0)), (Py::str("G"), Py::Float(1.0)), (Py::str("F"), Py::Float(1.0)), (Py::str("I"), Py::Float(1.0)), (Py::str("H"), Py::Float(13.34)), (Py::str("K"), Py::Float(-7.49)), (Py::str("M"), Py::Float(1.0)), (Py::str("L"), Py::Float(20.26)), (Py::str("N"), Py::Float(1.0)), (Py::str("Q"), Py::Float(1.0)), (Py::str("P"), Py::Float(-1.88)), (Py::str("S"), Py::Float(1.0)), (Py::str("R"), Py::Float(1.0)), (Py::str("T"), Py::Float(1.0)), (Py::str("W"), Py::Float(1.0)), (Py::str("V"), Py::Float(-7.49)), (Py::str("Y"), Py::Float(1.0))])), (Py::str("H"), Py::dict(vec![(Py::str("A"), Py::Float(1.0)), (Py::str("C"), Py::Float(1.0)), (Py::str("E"), Py::Float(1.0)), (Py::str("D"), Py::Float(1.0)), (Py::str("G"), Py::Float(-9.37)), (Py::str("F"), Py::Float(-9.37)), (Py::str("I"), Py::Float(44.94)), (Py::str("H"), Py::Float(1.0)), (Py::str("K"), Py::Float(24.68)), (Py::str("M"), Py::Float(1.0)), (Py::str("L"), Py::Float(1.0)), (Py::str("N"), Py::Float(24.68)), (Py::str("Q"), Py::Float(1.0)), (Py::str("P"), Py::Float(-1.88)), (Py::str("S"), Py::Float(1.0)), (Py::str("R"), Py::Float(1.0)), (Py::str("T"), Py::Float(-6.54)), (Py::str("W"), Py::Float(-1.88)), (Py::str("V"), Py::Float(1.0)), (Py::str("Y"), Py::Float(44.94))])), (Py::str("K"), Py::dict(vec![(Py::str("A"), Py::Float(1.0)), (Py::str("C"), Py::Float(1.0)), (Py::str("E"), Py::Float(1.0)), (Py::str("D"), Py::Float(1.0)), (Py::str("G"), Py::Float(-7.49)), (Py::str("F"), Py::Float(1.0)), (Py::str("I"), Py::Float(-7.49)), (Py::str("H"), Py::Float(1.0)), (Py::str("K"), Py::Float(1.0)), (Py::str("M"), Py::Float(33.6)), (Py::str("L"), Py::Float(-7.49)), (Py::str("N"), Py::Float(1.0)), (Py::str("Q"), Py::Float(24.64)), (Py::str("P"), Py::Float(-6.54)), (Py::str("S"), Py::Float(1.0)), (Py::str("R"), Py::Float(33.6)), (Py::str("T"), Py::Float(1.0)), (Py::str("W"), Py::Float(1.0)), (Py::str("V"), Py::Float(-7.49)), (Py::str("Y"), Py::Float(1.0))])), (Py::str("M"), Py::dict(vec![(Py::str("A"), Py::Float(13.34)), (Py::str("C"), Py::Float(1.0)), (Py::str("E"), Py::Float(1.0)), (Py::str("D"), Py::Float(1.0)), (Py::str("G"), Py::Float(1.0)), (Py::str("F"), Py::Float(1.0)), (Py::str("I"), Py::Float(1.0)), (Py::str("H"), Py::Float(58.28)), (Py::str("K"), Py::Float(1.0)), (Py::str("M"), Py::Float(-1.88)), (Py::str("L"), Py::Float(1.0)), (Py::str("N"), Py::Float(1.0)), (Py::str("Q"), Py::Float(-6.54)), (Py::str("P"), Py::Float(44.94)), (Py::str("S"), Py::Float(44.94)), (Py::str("R"), Py::Float(-6.54)), (Py::str("T"), Py::Float(-1.88)), (Py::str("W"), Py::Float(1.0)), (Py::str("V"), Py::Float(1.0)), (Py::str("Y"), Py::Float(24.68))])), (Py::str("L"), Py::dict(vec![(Py::str("A"), Py::Float(1.0)), (Py::str("C"), Py::Float(1.0)), (Py::str("E"), Py::Float(1.0)), (Py::str("D"), Py::Float(1.0)), (Py::str("G"), Py::Float(1.0)), (Py::str("F"), Py::Float(1.0)), (Py::str("I"), Py::Float(1.0)), (Py::str("H"), Py::Float(1.0)), (Py::str("K"), Py::Float(-7.49)), (Py::str("M"), Py::Float(1.0)), (Py::str("L"), Py::Float(1.0)), (Py::str("N"), Py::Float(1.0)), (Py::str("Q"), Py::Float(33.6)), (Py::str("P"), Py::Float(20.26)), (Py::str("S"), Py::Float(1.0)), (Py::str("R"), Py::Float(20.26)), (Py::str("T"), Py::Float(1.0)), (Py::str("W"), Py::Float(24.68)), (Py::str("V"), Py::Float(1.0)), (Py::str("Y"), Py::Float(1.0))])), (Py::str("N"), Py::dict(vec![(Py::str("A"), Py::Float(1.0)), (Py::str("C"), Py::Float(-1.88)), (Py::str("E"), Py::Float(1.0)), (Py::str("D"), Py::Float(1.0)), (Py::str("G"), Py::Float(-14.03)), (Py::str("F"), Py::Float(-14.03)), (Py::str("I"), Py::Float(44.94)), (Py::str("H"), Py::Float(1.0)), (Py::str("K"), Py::Float(24.68)), (Py::str("M"), Py::Float(1.0)), (Py::str("L"), Py::Float(1.0)), (Py::str("N"), Py::Float(1.0)), (Py::str("Q"), Py::Float(-6.54)), (Py::str("P"), Py::Float(-1.88)), (Py::str("S"), Py::Float(1.0)), (Py::str("R"), Py::Float(1.0)), (Py::str("T"), Py::Float(-7.49)), (Py::str("W"), Py::Float(-9.37)), (Py::str("V"), Py::Float(1.0)), (Py::str("Y"), Py::Float(1.0))])), (Py::str("Q"), Py::dict(vec![(Py::str("A"), Py::Float(1.0)), (Py::str("C"), Py::Float(-6.54)), (Py::str("E"), Py::Float(20.26)), (Py::str("D"), Py::Float(20.26)), (Py::str("G"), Py::Float(1.0)), (Py::str("F"), Py::Float(-6.54)), (Py::str("I"), Py::Float(1.0)), (Py::str("H"), Py::Float(1.0)), (Py::str("K"), Py::Float(1.0)), (Py::str("M"), Py::Float(1.0)), (Py::str("L"), Py::Float(1.0)), (Py::str("N"), Py::Float(1.0)), (Py::str("Q"), Py::Float(20.26)), (Py::str("P"), Py::Float(20.26)), (Py::str("S"), Py::Float(44.94)), (Py::str("R"), Py::Float(1.0)), (Py::str("T"), Py::Float(1.0)), (Py::str("W"), Py::Float(1.0)), (Py::str("V"), Py::Float(-6.54)), (Py::str("Y"), Py::Float(-6.54))])), (Py::str("P"), Py::dict(vec![(Py::str("A"), Py::Float(20.26)), (Py::str("C"), Py::Float(-6.54)), (Py::str("E"), Py::Float(18.38)), (Py::str("D"), Py::Float(-6.54)), (Py::str("G"), Py::Float(1.0)), (Py::str("F"), Py::Float(20.26)), (Py::str("I"), Py::Float(1.0)), (Py::str("H"), Py::Float(1.0)), (Py::str("K"), Py::Float(1.0)), (Py::str("M"), Py::Float(-6.54)), (Py::str("L"), Py::Float(1.0)), (Py::str("N"), Py::Float(1.0)), (Py::str("Q"), Py::Float(20.26)), (Py::str("P"), Py::Float(20.26)), (Py::str("S"), Py::Float(20.26)), (Py::str("R"), Py::Float(-6.54)), (Py::str("T"), Py::Float(1.0)), (Py::str("W"), Py::Float(-1.88)), (Py::str("V"), Py::Float(20.26)), (Py::str("Y"), Py::Float(1.0))])), (Py::str("S"), Py::dict(vec![(Py::str("A"), Py::Float(1.0)), (Py::str("C"), Py::Float(33.6)), (Py::str("E"), Py::Float(20.26)), (Py::str("D"), Py::Float(1.0)), (Py::str("G"), Py::Float(1.0)), (Py::str("F"), Py::Float(1.0)), (Py::str("I"), Py::Float(1.0)), (Py::str("H"), Py::Float(1.0)), (Py::str("K"), Py::Float(1.0)), (Py::str("M"), Py::Float(1.0)), (Py::str("L"), Py::Float(1.0)), (Py::str("N"), Py::Float(1.0)), (Py::str("Q"), Py::Float(20.26)), (Py::str("P"), Py::Float(44.94)), (Py::str("S"), Py::Float(20.26)), (Py::str("R"), Py::Float(20.26)), (Py::str("T"), Py::Float(1.0)), (Py::str("W"), Py::Float(1.0)), (Py::str("V"), Py::Float(1.0)), (Py::str("Y"), Py::Float(1.0))])), (Py::str("R"), Py::dict(vec![(Py::str("A"), Py::Float(1.0)), (Py::str("C"), Py::Float(1.0)), (Py::str("E"), Py::Float(1.0)), (Py::str("D"), Py::Float(1.0)), (Py::str("G"), Py::Float(-7.49)), (Py::str("F"), Py::Float(1.0)), (Py::str("I"), Py::Float(1.0)), (Py::str("H"), Py::Float(20.26)), (Py::str("K"), Py::Float(1.0)), (Py::str("M"), Py::Float(1.0)), (Py::str("L"), Py::Float(1.0)), (Py::str("N"), Py::Float(13.34)), (Py::str("Q"), Py::Float(20.26)), (Py::str("P"), Py::Float(20.26)), (Py::str("S"), Py::Float(44.94)), (Py::str("R"), Py::Float(58.28)), (Py::str("T"), Py::Float(1.0)), (Py::str("W"), Py::Float(58.28)), (Py::str("V"), Py::Float(1.0)), (Py::str("Y"), Py::Float(-6.54))])), (Py::str("T"), Py::dict(vec![(Py::str("A"), Py::Float(1.0)), (Py::str("C"), Py::Float(1.0)), (Py::str("E"), Py::Float(20.26)), (Py::str("D"), Py::Float(1.0)), (Py::str("G"), Py::Float(-7.49)), (Py::str("F"), Py::Float(13.34)), (Py::str("I"), Py::Float(1.0)), (Py::str("H"), Py::Float(1.0)), (Py::str("K"), Py::Float(1.0)), (Py::str("M"), Py::Float(1.0)), (Py::str("L"), Py::Float(1.0)), (Py::str("N"), Py::Float(-14.03)), (Py::str("Q"), Py::Float(-6.54)), (Py::str("P"), Py::Float(1.0)), (Py::str("S"), Py::Float(1.0)), (Py::str("R"), Py::Float(1.0)), (Py::str("T"), Py::Float(1.0)), (Py::str("W"), Py::Float(-14.03)), (Py::str("V"), Py::Float(1.0)), (Py::str("Y"), Py::Float(1.0))])), (Py::str("W"), Py::dict(vec![(Py::str("A"), Py::Float(-14.03)), (Py::str("C"), Py::Float(1.0)), (Py::str("E"), Py::Float(1.0)), (Py::str("D"), Py::Float(1.0)), (Py::str("G"), Py::Float(-9.37)), (Py::str("F"), Py::Float(1.0)), (Py::str("I"), Py::Float(1.0)), (Py::str("H"), Py::Float(24.68)), (Py::str("K"), Py::Float(1.0)), (Py::str("M"), Py::Float(24.68)), (Py::str("L"), Py::Float(13.34)), (Py::str("N"), Py::Float(13.34)), (Py::str("Q"), Py::Float(1.0)), (Py::str("P"), Py::Float(1.0)), (Py::str("S"), Py::Float(1.0)), (Py::str("R"), Py::Float(1.0)), (Py::str("T"), Py::Float(-14.03)), (Py::str("W"), Py::Float(1.0)), (Py::str("V"), Py::Float(-7.49)), (Py::str("Y"), Py::Float(1.0))])), (Py::str("V"), Py::dict(vec![(Py::str("A"), Py::Float(1.0)), (Py::str("C"), Py::Float(1.0)), (Py::str("E"), Py::Float(1.0)), (Py::str("D"), Py::Float(-14.03)), (Py::str("G"), Py::Float(-7.49)), (Py::str("F"), Py::Float(1.0)), (Py::str("I"), Py::Float(1.0)), (Py::str("H"), Py::Float(1.0)), (Py::str("K"), Py::Float(-1.88)), (Py::str("M"), Py::Float(1.0)), (Py::str("L"), Py::Float(1.0)), (Py::str("N"), Py::Float(1.0)), (Py::str("Q"), Py::Float(1.0)), (Py::str("P"), Py::Float(20.26)), (Py::str("S"), Py::Float(1.0)), (Py::str("R"), Py::Float(1.0)), (Py::str("T"), Py::Float(-7.49)), (Py::str("W"), Py::Float(1.0)), (Py::str("V"), Py::Float(1.0)), (Py::str("Y"), Py::Float(-6.54))])), (Py::str("Y"), Py::dict(vec![(Py::str("A"), Py::Float(24.68)), (Py::str("C"), Py::Float(1.0)), (Py::str("E"), Py::Float(-6.54)), (Py::str("D"), Py::Float(24.68)), (Py::str("G"), Py::Float(-7.49)), (Py::str("F"), Py::Float(1.0)), (Py::str("I"), Py::Float(1.0)), (Py::str("H"), Py::Float(13.34)), (Py::str("K"), Py::Float(1.0)), (Py::str("M"), Py::Float(44.94)), (Py::str("L"), Py::Float(1.0)), (Py::str("N"), Py::Float(1.0)), (Py::str("Q"), Py::Float(1.0)), (Py::str("P"), Py::Float(13.34)), (Py::str("S"), Py::Float(1.0)), (Py::str("R"), Py::Float(-15.91)), (Py::str("T"), Py::Float(-7.49)), (Py::str("W"), Py::Float(-9.37)), (Py::str("V"), Py::Float(1.0)), (Py::str("Y"), Py::Float(13.34))]))])) })());
    V.with(|v| v.clone())
}

#[allow(non_snake_case)]
fn c_DNA_NN3() -> R {
    thread_local!(static V: R = (|| -> R { Ok(Py::dict(vec![(Py::str("init"), Py::tuple(vec![Py::Int(0i64), Py::Int(0i64)])), (Py::str("init_A/T"), Py::tuple(vec![Py::Float(2.3), Py::Float(4.1)])), (Py::str("init_G/C"), Py::tuple(vec![Py::Float(0.1), Py::Float(-2.8)])), (Py::str("init_oneG/C"), Py::tuple(vec![Py::Int(0i64), Py::Int(0i64)])), (Py::str("init_allA/T"), Py::tuple(vec![Py::Int(0i64), Py::Int(0i64)])), (Py::str("init_5T/A"), Py::tuple(vec![Py::Int(0i64), Py::Int(0i64)])), (Py::str("sym"), Py::tuple(vec![Py::Int(0i64), Py::Float(-1.4)])), (Py::str("AA/TT"), Py::tuple(vec![Py::Float(-7.9), Py::Float(-22.2)])), (Py::str("AT/TA"), Py::tuple(vec![Py::Float(-7.2), Py::Float(-20.4)])), (Py::str("TA/AT"), Py::tuple(vec![Py::Float(-7.2), Py::Float(-21.3)])), (Py::str("CA/GT"), Py::tuple(vec![Py::Float(-8.5), Py::Float(-22.7)])), (Py::str("GT/CA"), Py::tuple(vec![Py::Float(-8.4), Py::Float(-22.4)])), (Py::str("CT/GA"), Py::tuple(vec![Py::Float(-7.8), Py::Float(-21.0)])), (Py::str("GA/CT"), Py::tuple(vec![Py::Float(-8.2), Py::Float(-22.2)])), (Py::str("CG/GC"), Py::tuple(vec![Py::Float(-10.6), Py::Float(-27.2)])), (Py::str("GC/CG"), Py::tuple(vec![Py::Float(-9.8), Py::Float(-24.4)])), (Py::str("GG/CC"), Py::tuple(vec![Py::Float(-8.0), Py::Float(-19.9)]))])) })());
    V.with(|v| v.clone())
}

#[allow(non_snake_case)]
fn c_DNA_W() -> R {
    thread_local!(static V: R = (|| -> R { Ok(Py::dict(vec![(Py::str("A"), Py::Float(331.2218)), (Py::str("C"), Py::Float(307.1971)), (Py::str("G"), Py::Float(347.2212)), (Py::str("T"), Py::Float(322.2085))])) })());
    V.with(|v| v.clone())
}

#[allow(non_snake_case)]
fn c_ENZYMES() -> R {
    thread_local!(static V: R = (|| -> R { Ok(Py::dict(vec![(Py::str("EcoRI"), Py::tuple(vec![Py::str("GAATTC"), Py::Int(1i64)])), (Py::str("BamHI"), Py::tuple(vec![Py::str("GGATCC"), Py::Int(1i64)])), (Py::str("HindIII"), Py::tuple(vec![Py::str("AAGCTT"), Py::Int(1i64)])), (Py::str("PstI"), Py::tuple(vec![Py::str("CTGCAG"), Py::Int(5i64)])), (Py::str("SmaI"), Py::tuple(vec![Py::str("CCCGGG"), Py::Int(3i64)])), (Py::str("KpnI"), Py::tuple(vec![Py::str("GGTACC"), Py::Int(5i64)])), (Py::str("SacI"), Py::tuple(vec![Py::str("GAGCTC"), Py::Int(5i64)])), (Py::str("XbaI"), Py::tuple(vec![Py::str("TCTAGA"), Py::Int(1i64)])), (Py::str("SalI"), Py::tuple(vec![Py::str("GTCGAC"), Py::Int(1i64)])), (Py::str("XhoI"), Py::tuple(vec![Py::str("CTCGAG"), Py::Int(1i64)])), (Py::str("NotI"), Py::tuple(vec![Py::str("GCGGCCGC"), Py::Int(2i64)])), (Py::str("EcoRV"), Py::tuple(vec![Py::str("GATATC"), Py::Int(3i64)])), (Py::str("NcoI"), Py::tuple(vec![Py::str("CCATGG"), Py::Int(1i64)])), (Py::str("NdeI"), Py::tuple(vec![Py::str("CATATG"), Py::Int(2i64)])), (Py::str("HaeIII"), Py::tuple(vec![Py::str("GGCC"), Py::Int(2i64)])), (Py::str("AluI"), Py::tuple(vec![Py::str("AGCT"), Py::Int(2i64)])), (Py::str("TaqI"), Py::tuple(vec![Py::str("TCGA"), Py::Int(1i64)])), (Py::str("MspI"), Py::tuple(vec![Py::str("CCGG"), Py::Int(1i64)])), (Py::str("ClaI"), Py::tuple(vec![Py::str("ATCGAT"), Py::Int(2i64)])), (Py::str("BglII"), Py::tuple(vec![Py::str("AGATCT"), Py::Int(1i64)])), (Py::str("SpeI"), Py::tuple(vec![Py::str("ACTAGT"), Py::Int(1i64)])), (Py::str("NheI"), Py::tuple(vec![Py::str("GCTAGC"), Py::Int(1i64)])), (Py::str("ApaI"), Py::tuple(vec![Py::str("GGGCCC"), Py::Int(5i64)])), (Py::str("SphI"), Py::tuple(vec![Py::str("GCATGC"), Py::Int(5i64)])), (Py::str("MluI"), Py::tuple(vec![Py::str("ACGCGT"), Py::Int(1i64)])), (Py::str("ScaI"), Py::tuple(vec![Py::str("AGTACT"), Py::Int(3i64)])), (Py::str("PvuII"), Py::tuple(vec![Py::str("CAGCTG"), Py::Int(3i64)]))])) })());
    V.with(|v| v.clone())
}

#[allow(non_snake_case)]
fn c_KD() -> R {
    thread_local!(static V: R = (|| -> R { Ok(Py::dict(vec![(Py::str("A"), Py::Float(1.8)), (Py::str("R"), Py::Float(-4.5)), (Py::str("N"), Py::Float(-3.5)), (Py::str("D"), Py::Float(-3.5)), (Py::str("C"), Py::Float(2.5)), (Py::str("Q"), Py::Float(-3.5)), (Py::str("E"), Py::Float(-3.5)), (Py::str("G"), Py::Float(-0.4)), (Py::str("H"), Py::Float(-3.2)), (Py::str("I"), Py::Float(4.5)), (Py::str("L"), Py::Float(3.8)), (Py::str("K"), Py::Float(-3.9)), (Py::str("M"), Py::Float(1.9)), (Py::str("F"), Py::Float(2.8)), (Py::str("P"), Py::Float(-1.6)), (Py::str("S"), Py::Float(-0.8)), (Py::str("T"), Py::Float(-0.7)), (Py::str("W"), Py::Float(-0.9)), (Py::str("Y"), Py::Float(-1.3)), (Py::str("V"), Py::Float(4.2))])) })());
    V.with(|v| v.clone())
}

#[allow(non_snake_case)]
fn c_PK_CTERM() -> R {
    thread_local!(static V: R = (|| -> R { Ok(Py::dict(vec![(Py::str("D"), Py::Float(4.55)), (Py::str("E"), Py::Float(4.75))])) })());
    V.with(|v| v.clone())
}

#[allow(non_snake_case)]
fn c_PK_NEG() -> R {
    thread_local!(static V: R = (|| -> R { Ok(Py::dict(vec![(Py::str("Cterm"), Py::Float(3.55)), (Py::str("D"), Py::Float(4.05)), (Py::str("E"), Py::Float(4.45)), (Py::str("C"), Py::Float(9.0)), (Py::str("Y"), Py::Float(10.0))])) })());
    V.with(|v| v.clone())
}

#[allow(non_snake_case)]
fn c_PK_NTERM() -> R {
    thread_local!(static V: R = (|| -> R { Ok(Py::dict(vec![(Py::str("A"), Py::Float(7.59)), (Py::str("M"), Py::Float(7.0)), (Py::str("S"), Py::Float(6.93)), (Py::str("P"), Py::Float(8.36)), (Py::str("T"), Py::Float(6.82)), (Py::str("V"), Py::Float(7.44)), (Py::str("E"), Py::Float(7.7))])) })());
    V.with(|v| v.clone())
}

#[allow(non_snake_case)]
fn c_PK_POS() -> R {
    thread_local!(static V: R = (|| -> R { Ok(Py::dict(vec![(Py::str("Nterm"), Py::Float(7.5)), (Py::str("K"), Py::Float(10.0)), (Py::str("R"), Py::Float(12.0)), (Py::str("H"), Py::Float(5.98))])) })());
    V.with(|v| v.clone())
}

#[allow(non_snake_case)]
fn c_PROTEIN_W() -> R {
    thread_local!(static V: R = (|| -> R { Ok(Py::dict(vec![(Py::str("A"), Py::Float(89.0932)), (Py::str("C"), Py::Float(121.1582)), (Py::str("D"), Py::Float(133.1027)), (Py::str("E"), Py::Float(147.1293)), (Py::str("F"), Py::Float(165.1891)), (Py::str("G"), Py::Float(75.0666)), (Py::str("H"), Py::Float(155.1546)), (Py::str("I"), Py::Float(131.1729)), (Py::str("K"), Py::Float(146.1876)), (Py::str("L"), Py::Float(131.1729)), (Py::str("M"), Py::Float(149.2113)), (Py::str("N"), Py::Float(132.1179)), (Py::str("P"), Py::Float(115.1305)), (Py::str("Q"), Py::Float(146.1445)), (Py::str("R"), Py::Float(174.201)), (Py::str("S"), Py::Float(105.0926)), (Py::str("T"), Py::Float(119.1192)), (Py::str("V"), Py::Float(117.1463)), (Py::str("W"), Py::Float(204.2252)), (Py::str("Y"), Py::Float(181.1885))])) })());
    V.with(|v| v.clone())
}

#[allow(non_snake_case)]
fn c_RNA_W() -> R {
    thread_local!(static V: R = (|| -> R { Ok(Py::dict(vec![(Py::str("A"), Py::Float(347.2212)), (Py::str("C"), Py::Float(323.1965)), (Py::str("G"), Py::Float(363.2206)), (Py::str("U"), Py::Float(324.1813))])) })());
    V.with(|v| v.clone())
}

#[allow(non_snake_case)]
fn c_STOPS() -> R {
    thread_local!(static V: R = (|| -> R { Ok(Py::tuple(vec![Py::str("TAA"), Py::str("TAG"), Py::str("TGA")])) })());
    V.with(|v| v.clone())
}

#[allow(non_snake_case)]
fn c_WATER() -> R {
    thread_local!(static V: R = (|| -> R { Ok(Py::Float(18.0153)) })());
    V.with(|v| v.clone())
}

#[allow(non_snake_case)]
fn c__COMP() -> R {
    thread_local!(static V: R = (|| -> R { Ok(Py::dict(vec![(Py::Int(65i64), Py::Int(84i64)), (Py::Int(67i64), Py::Int(71i64)), (Py::Int(71i64), Py::Int(67i64)), (Py::Int(84i64), Py::Int(65i64)), (Py::Int(85i64), Py::Int(65i64))])) })());
    V.with(|v| v.clone())
}

#[allow(non_snake_case)]
fn c__COMP_RNA() -> R {
    thread_local!(static V: R = (|| -> R { Ok(Py::dict(vec![(Py::Int(65i64), Py::Int(85i64)), (Py::Int(67i64), Py::Int(71i64)), (Py::Int(71i64), Py::Int(67i64)), (Py::Int(85i64), Py::Int(65i64))])) })());
    V.with(|v| v.clone())
}

#[allow(non_snake_case)]
fn c__HELIX() -> R {
    thread_local!(static V: R = (|| -> R { Ok(Py::str("EMALK")) })());
    V.with(|v| v.clone())
}

#[allow(non_snake_case)]
fn c__LOGS() -> R {
    thread_local!(static V: R = (|| -> R { Ok(Py::dict(vec![(Py::str("e"), Py::Float(2.718281828459045)), (Py::str("2"), Py::Float(2.0)), (Py::str("10"), Py::Float(10.0))])) })());
    V.with(|v| v.clone())
}

#[allow(non_snake_case)]
fn c__NEG() -> R {
    thread_local!(static V: R = (|| -> R { Ok(Py::Float(f64::NEG_INFINITY)) })());
    V.with(|v| v.clone())
}

#[allow(non_snake_case)]
fn c__PUR() -> R {
    thread_local!(static V: R = (|| -> R { Ok(py::set(&Py::str("AG"))?) })());
    V.with(|v| v.clone())
}

#[allow(non_snake_case)]
fn c__SHEET() -> R {
    thread_local!(static V: R = (|| -> R { Ok(Py::str("VIYFWLT")) })());
    V.with(|v| v.clone())
}

#[allow(non_snake_case)]
fn c__TURN() -> R {
    thread_local!(static V: R = (|| -> R { Ok(Py::str("NPGSD")) })());
    V.with(|v| v.clone())
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn f__abund(mut v_counts: Py, mut v_name: Py) -> R {
    let mut v_c: Py = Py::None;
    let mut v_pos: Py = Py::None;
    v_c = f__counts(v_name.clone(), v_counts.clone(), Py::Int(0i64), Py::Bool(false))?;
    v_pos = { let mut out = vec![]; for t6428 in py::iter(&(v_c.clone()))? { let mut v_x: Py; 
v_x = t6428;
if ({ let t6429 = v_x.clone(); Py::Bool({ let t6430 = Py::Int(0i64); if py::gt(&t6429, &t6430)? { true } else { false } }) }).truthy() { out.push(v_x.clone()); }} Py::list(out) };
    let _ = f__need(v_pos.clone(), Py::str(&[py::format(&v_name.clone(), "")?, ": at least one species must have a positive count".to_string()].concat()))?;
    return Ok(v_pos.clone());
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn f__align(mut v_a: Py, mut v_b: Py, mut v_sc: Py, mut v_op: Py, mut v_ex: Py, mut v_mode: Py, mut v_wrong_gaps: Py) -> R {
    let mut v_E: Py = Py::None;
    let mut v_Ei: Py = Py::None;
    let mut v_F: Py = Py::None;
    let mut v_Fi: Py = Py::None;
    let mut v_Fp: Py = Py::None;
    let mut v_H: Py = Py::None;
    let mut v_Hi: Py = Py::None;
    let mut v_Hp: Py = Py::None;
    let mut v_ai: Py = Py::None;
    let mut v_best: Py = Py::None;
    let mut v_e: Py = Py::None;
    let mut v_f: Py = Py::None;
    let mut v_h: Py = Py::None;
    let mut v_i: Py = Py::None;
    let mut v_j: Py = Py::None;
    let mut v_m: Py = Py::None;
    let mut v_n: Py = Py::None;
    let mut v_o2: Py = Py::None;
    let t6263 = (Py::tuple(vec![py::len(&v_a.clone())?, py::len(&v_b.clone())?])).items()?;
    if t6263.len() != 2 { return Err(PyErr::Value("unpack".into())); }
    v_n = t6263[0].clone();
    v_m = t6263[1].clone();
    v_o2 = (if (v_wrong_gaps.clone()).truthy() { py::add(&v_op.clone(), &v_ex.clone())? } else { v_op.clone() });
    v_H = { let mut out = vec![]; for t6264 in py::range_iter(&[(py::add(&v_n.clone(), &Py::Int(1i64))?).i()?]) { let mut v__: Py; 
v__ = t6264;
out.push(py::mul(&Py::list(vec![Py::Float(0.0)]), &py::add(&v_m.clone(), &Py::Int(1i64))?)?); } Py::list(out) };
    v_E = { let mut out = vec![]; for t6265 in py::range_iter(&[(py::add(&v_n.clone(), &Py::Int(1i64))?).i()?]) { let mut v__: Py; 
v__ = t6265;
out.push(py::mul(&Py::list(vec![c__NEG()?]), &py::add(&v_m.clone(), &Py::Int(1i64))?)?); } Py::list(out) };
    v_F = { let mut out = vec![]; for t6266 in py::range_iter(&[(py::add(&v_n.clone(), &Py::Int(1i64))?).i()?]) { let mut v__: Py; 
v__ = t6266;
out.push(py::mul(&Py::list(vec![c__NEG()?]), &py::add(&v_m.clone(), &Py::Int(1i64))?)?); } Py::list(out) };
    for t6267 in py::range_iter(&[(Py::Int(1i64)).i()?, (py::add(&v_n.clone(), &Py::Int(1i64))?).i()?]) {
        v_i = t6267;
        if ({ let t6268 = v_mode.clone(); Py::Bool({ let t6269 = Py::str("global"); if py::eq(&t6268, &t6269) { true } else { false } }) }).truthy() {
            let t6270 = py::add(&v_o2.clone(), &py::mul(&py::sub(&v_i.clone(), &Py::Int(1i64))?, &v_ex.clone())?)?;
            py::setitem(&py::getitem(&v_H.clone(), &v_i.clone())?, &Py::Int(0i64), t6270.clone())?;
            py::setitem(&py::getitem(&v_F.clone(), &v_i.clone())?, &Py::Int(0i64), t6270.clone())?;
        }
    }
    for t6271 in py::range_iter(&[(Py::Int(1i64)).i()?, (py::add(&v_m.clone(), &Py::Int(1i64))?).i()?]) {
        v_j = t6271;
        if ({ let t6272 = v_mode.clone(); Py::Bool({ let t6273 = Py::str("global"); if py::eq(&t6272, &t6273) { true } else { false } }) }).truthy() {
            let t6274 = py::add(&v_o2.clone(), &py::mul(&py::sub(&v_j.clone(), &Py::Int(1i64))?, &v_ex.clone())?)?;
            py::setitem(&py::getitem(&v_H.clone(), &Py::Int(0i64))?, &v_j.clone(), t6274.clone())?;
            py::setitem(&py::getitem(&v_E.clone(), &Py::Int(0i64))?, &v_j.clone(), t6274.clone())?;
        }
    }
    v_best = Py::Float(0.0);
    for t6275 in py::range_iter(&[(Py::Int(1i64)).i()?, (py::add(&v_n.clone(), &Py::Int(1i64))?).i()?]) {
        v_i = t6275;
        let t6276 = (Py::tuple(vec![py::getitem(&v_H.clone(), &v_i.clone())?, py::getitem(&v_H.clone(), &py::sub(&v_i.clone(), &Py::Int(1i64))?)?, py::getitem(&v_E.clone(), &v_i.clone())?, py::getitem(&v_F.clone(), &v_i.clone())?, py::getitem(&v_F.clone(), &py::sub(&v_i.clone(), &Py::Int(1i64))?)?])).items()?;
        if t6276.len() != 5 { return Err(PyErr::Value("unpack".into())); }
        v_Hi = t6276[0].clone();
        v_Hp = t6276[1].clone();
        v_Ei = t6276[2].clone();
        v_Fi = t6276[3].clone();
        v_Fp = t6276[4].clone();
        v_ai = py::getitem(&v_a.clone(), &py::sub(&v_i.clone(), &Py::Int(1i64))?)?;
        for t6277 in py::range_iter(&[(Py::Int(1i64)).i()?, (py::add(&v_m.clone(), &Py::Int(1i64))?).i()?]) {
            v_j = t6277;
            v_e = py::minmax(&[py::add(&py::getitem(&v_Ei.clone(), &py::sub(&v_j.clone(), &Py::Int(1i64))?)?, &v_ex.clone())?, py::add(&py::getitem(&v_Hi.clone(), &py::sub(&v_j.clone(), &Py::Int(1i64))?)?, &v_o2.clone())?], true)?;
            v_f = py::minmax(&[py::add(&py::getitem(&v_Fp.clone(), &v_j.clone())?, &v_ex.clone())?, py::add(&py::getitem(&v_Hp.clone(), &v_j.clone())?, &v_o2.clone())?], true)?;
            v_h = py::minmax(&[py::add(&py::getitem(&v_Hp.clone(), &py::sub(&v_j.clone(), &Py::Int(1i64))?)?, &py::call(&v_sc.clone(), &[v_ai.clone(), py::getitem(&v_b.clone(), &py::sub(&v_j.clone(), &Py::Int(1i64))?)?])?)?, v_e.clone(), v_f.clone()], true)?;
            if ({ let t6282 = { let t6278 = v_mode.clone(); Py::Bool({ let t6279 = Py::str("local"); if py::eq(&t6278, &t6279) { true } else { false } }) }; if !t6282.truthy() { t6282 } else { { let t6280 = v_h.clone(); Py::Bool({ let t6281 = Py::Int(0i64); if py::lt(&t6280, &t6281)? { true } else { false } }) } } }).truthy() {
                v_h = Py::Float(0.0);
            }
            let t6283 = (Py::tuple(vec![v_e.clone(), v_f.clone(), v_h.clone()])).items()?;
            if t6283.len() != 3 { return Err(PyErr::Value("unpack".into())); }
            py::setitem(&v_Ei.clone(), &v_j.clone(), t6283[0].clone())?;
            py::setitem(&v_Fi.clone(), &v_j.clone(), t6283[1].clone())?;
            py::setitem(&v_Hi.clone(), &v_j.clone(), t6283[2].clone())?;
            if ({ let t6284 = v_h.clone(); Py::Bool({ let t6285 = v_best.clone(); if py::gt(&t6284, &t6285)? { true } else { false } }) }).truthy() {
                v_best = v_h.clone();
            }
        }
    }
    if ({ let t6286 = v_mode.clone(); Py::Bool({ let t6287 = Py::str("global"); if py::eq(&t6286, &t6287) { true } else { false } }) }).truthy() {
        return Ok(py::getitem(&py::getitem(&v_H.clone(), &v_n.clone())?, &v_m.clone())?);
    }
    if ({ let t6288 = v_mode.clone(); Py::Bool({ let t6289 = Py::str("local"); if py::eq(&t6288, &t6289) { true } else { false } }) }).truthy() {
        return Ok(v_best.clone());
    }
    return Ok(py::minmax(&[py::add(&{ let mut out = vec![]; for t6290 in py::range_iter(&[(py::add(&v_m.clone(), &Py::Int(1i64))?).i()?]) { let mut v_j: Py; 
v_j = t6290;
out.push(py::getitem(&py::getitem(&v_H.clone(), &v_n.clone())?, &v_j.clone())?); } Py::list(out) }, &{ let mut out = vec![]; for t6291 in py::range_iter(&[(py::add(&v_n.clone(), &Py::Int(1i64))?).i()?]) { let mut v_i: Py; 
v_i = t6291;
out.push(py::getitem(&py::getitem(&v_H.clone(), &v_i.clone())?, &v_m.clone())?); } Py::list(out) })?], true)?);
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn f__aligned_dna(mut v_seq1: Py, mut v_seq2: Py) -> R {
    let mut v_a: Py = Py::None;
    let mut v_b: Py = Py::None;
    let t3208 = (Py::tuple(vec![f__seq(v_seq1.clone(), Py::str("seq1"), Py::str("DNA"))?, f__seq(v_seq2.clone(), Py::str("seq2"), Py::str("DNA"))?])).items()?;
    if t3208.len() != 2 { return Err(PyErr::Value("unpack".into())); }
    v_a = t3208[0].clone();
    v_b = t3208[1].clone();
    let _ = f__need({ let t3209 = py::len(&v_a.clone())?; Py::Bool({ let t3210 = py::len(&v_b.clone())?; if py::eq(&t3209, &t3210) { true } else { false } }) }, Py::str("the two aligned sequences must have the same length (no gaps)"))?;
    return Ok(Py::tuple(vec![v_a.clone(), v_b.clone()]));
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn f__alignment(mut v_seqs: Py, mut v_lo: Py) -> R {
    let mut v_r: Py = Py::None;
    let mut v_rows: Py = Py::None;
    let _ = f__need(Py::Bool(matches!(v_seqs.clone(), Py::Str(_))), Py::str("give the aligned sequences separated by spaces"))?;
    v_rows = { let mut out = vec![]; for t7106 in py::iter(&(py::method(&v_seqs.clone(), "split", &[])?))? { let mut v_s: Py; 
v_s = t7106;
out.push(py::method(&v_s.clone(), "upper", &[])?); } Py::list(out) };
    let _ = f__need({ let t7107 = py::len(&v_rows.clone())?; Py::Bool({ let t7108 = v_lo.clone(); if py::ge(&t7107, &t7108)? { true } else { false } }) }, Py::str(&["at least ".to_string(), py::format(&v_lo.clone(), "")?, " sequences are needed".to_string()].concat()))?;
    for t7109 in py::iter(&(v_rows.clone()))? {
        v_r = t7109;
        let _ = f__need({ let t7110 = ax_engine::pyre::module_call("fullmatch", &Py::str("[ACGT]+"), &[v_r.clone()], 0)?; Py::Bool({ let t7111 = Py::None; if !py::is(&t7110, &t7111) { true } else { false } }) }, Py::str("sequences must be DNA (A, C, G, T) without gaps"))?;
    }
    let _ = f__need({ let t7112 = py::len(&py::set(&{ let mut out = vec![]; for t7114 in py::iter(&(v_rows.clone()))? { let mut v_r: Py; 
v_r = t7114;
out.push(py::len(&v_r.clone())?); } Py::list(out) })?)?; Py::Bool({ let t7113 = Py::Int(1i64); if py::eq(&t7112, &t7113) { true } else { false } }) }, Py::str("the sequences must be aligned (equal lengths)"))?;
    return Ok(v_rows.clone());
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn f__bisect(mut v_f: Py, mut v_a: Py, mut v_b: Py, mut v_n: Py) -> R {
    let mut v__: Py = Py::None;
    let mut v_fa: Py = Py::None;
    let mut v_fb: Py = Py::None;
    let mut v_fm: Py = Py::None;
    let mut v_m: Py = Py::None;
    let t4490 = (Py::tuple(vec![py::call(&v_f.clone(), &[v_a.clone()])?, py::call(&v_f.clone(), &[v_b.clone()])?])).items()?;
    if t4490.len() != 2 { return Err(PyErr::Value("unpack".into())); }
    v_fa = t4490[0].clone();
    v_fb = t4490[1].clone();
    if ({ let t4491 = v_fa.clone(); Py::Bool({ let t4492 = Py::Int(0i64); if py::eq(&t4491, &t4492) { true } else { false } }) }).truthy() {
        return Ok(v_a.clone());
    }
    if ({ let t4493 = v_fb.clone(); Py::Bool({ let t4494 = Py::Int(0i64); if py::eq(&t4493, &t4494) { true } else { false } }) }).truthy() {
        return Ok(v_b.clone());
    }
    let _ = f__need({ let t4495 = { let t4497 = v_fa.clone(); Py::Bool({ let t4498 = Py::Int(0i64); if py::gt(&t4497, &t4498)? { true } else { false } }) }; Py::Bool({ let t4496 = { let t4499 = v_fb.clone(); Py::Bool({ let t4500 = Py::Int(0i64); if py::gt(&t4499, &t4500)? { true } else { false } }) }; if py::ne(&t4495, &t4496) { true } else { false } }) }, Py::str("no root in the interval"))?;
    for t4501 in py::range_iter(&[(v_n.clone()).i()?]) {
        v__ = t4501;
        v_m = py::mul(&Py::Float(0.5), &py::add(&v_a.clone(), &v_b.clone())?)?;
        if ({ let t4506 = { let t4502 = v_m.clone(); Py::Bool({ let t4503 = v_a.clone(); if py::eq(&t4502, &t4503) { true } else { false } }) }; if t4506.truthy() { t4506 } else { { let t4504 = v_m.clone(); Py::Bool({ let t4505 = v_b.clone(); if py::eq(&t4504, &t4505) { true } else { false } }) } } }).truthy() {
            break;
        }
        v_fm = py::call(&v_f.clone(), &[v_m.clone()])?;
        if ({ let t4507 = v_fm.clone(); Py::Bool({ let t4508 = Py::Int(0i64); if py::eq(&t4507, &t4508) { true } else { false } }) }).truthy() {
            return Ok(v_m.clone());
        }
        if ({ let t4509 = { let t4511 = v_fm.clone(); Py::Bool({ let t4512 = Py::Int(0i64); if py::gt(&t4511, &t4512)? { true } else { false } }) }; Py::Bool({ let t4510 = { let t4513 = v_fa.clone(); Py::Bool({ let t4514 = Py::Int(0i64); if py::gt(&t4513, &t4514)? { true } else { false } }) }; if py::eq(&t4509, &t4510) { true } else { false } }) }).truthy() {
            let t4515 = (Py::tuple(vec![v_m.clone(), v_fm.clone()])).items()?;
            if t4515.len() != 2 { return Err(PyErr::Value("unpack".into())); }
            v_a = t4515[0].clone();
            v_fa = t4515[1].clone();
        } else {
            v_b = v_m.clone();
        }
    }
    return Ok(py::mul(&Py::Float(0.5), &py::add(&v_a.clone(), &v_b.clone())?)?);
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn f__blo(mut v_x: Py, mut v_y: Py) -> R {
    return Ok(py::getitem(&c_BLOSUM62()?, &Py::tuple(vec![v_x.clone(), v_y.clone()]))?);
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn f__chain_mass(mut v_s: Py, mut v_table: Py) -> R {
    return Ok(py::sub(&py::sum(&{ let mut out = vec![]; for t5425 in py::iter(&(v_s.clone()))? { let mut v_x: Py; 
v_x = t5425;
out.push(py::getitem(&v_table.clone(), &v_x.clone())?); } Py::list(out) }, None)?, &py::mul(&py::sub(&py::len(&v_s.clone())?, &Py::Int(1i64))?, &c_WATER()?)?)?);
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn f__charge(mut v_s: Py, mut v_pH: Py, mut v_termini: Py) -> R {
    let mut v_aa: Py = Py::None;
    let mut v_minus: Py = Py::None;
    let mut v_n: Py = Py::None;
    let mut v_neg: Py = Py::None;
    let mut v_pk: Py = Py::None;
    let mut v_plus: Py = Py::None;
    let mut v_pos: Py = Py::None;
    let t5510 = (Py::tuple(vec![py::dict_from(&c_PK_POS()?, vec![])?, py::dict_from(&c_PK_NEG()?, vec![])?])).items()?;
    if t5510.len() != 2 { return Err(PyErr::Value("unpack".into())); }
    v_pos = t5510[0].clone();
    v_neg = t5510[1].clone();
    if ({ let t5511 = py::getitem(&v_s.clone(), &Py::Int(0i64))?; Py::Bool({ let t5512 = c_PK_NTERM()?; if py::contains(&t5512, &t5511)? { true } else { false } }) }).truthy() {
        py::setitem(&v_pos.clone(), &Py::str("Nterm"), py::getitem(&c_PK_NTERM()?, &py::getitem(&v_s.clone(), &Py::Int(0i64))?)?)?;
    }
    if ({ let t5513 = py::getitem(&v_s.clone(), &py::neg(&Py::Int(1i64))?)?; Py::Bool({ let t5514 = c_PK_CTERM()?; if py::contains(&t5514, &t5513)? { true } else { false } }) }).truthy() {
        py::setitem(&v_neg.clone(), &Py::str("Cterm"), py::getitem(&c_PK_CTERM()?, &py::getitem(&v_s.clone(), &py::neg(&Py::Int(1i64))?)?)?)?;
    }
    v_n = { let mut out = vec![]; for t5515 in py::iter(&(Py::str("KRHDECY")))? { let mut v_aa: Py; 
v_aa = t5515;
out.push((v_aa.clone(), py::float(&py::method(&v_s.clone(), "count", &[v_aa.clone()])?)?)); } Py::dict(out) };
    let t5516 = (if (v_termini.clone()).truthy() { Py::Float(1.0) } else { Py::Float(0.0) });
    py::setitem(&v_n.clone(), &Py::str("Nterm"), t5516.clone())?;
    py::setitem(&v_n.clone(), &Py::str("Cterm"), t5516.clone())?;
    v_plus = Py::Float(0.0);
    for t5517 in py::iter(&(py::method(&v_pos.clone(), "items", &[])?))? {
        let t5518 = (t5517).items()?;
        if t5518.len() != 2 { return Err(PyErr::Value("unpack".into())); }
        v_aa = t5518[0].clone();
        v_pk = t5518[1].clone();
        v_plus = py::iadd(&v_plus, &py::mul(&py::getitem(&v_n.clone(), &v_aa.clone())?, &py::truediv(&Py::Float(1.0), &py::add(&py::pow(&Py::Int(10i64), &py::sub(&v_pH.clone(), &v_pk.clone())?)?, &Py::Float(1.0))?)?)?)?;
    }
    v_minus = Py::Float(0.0);
    for t5519 in py::iter(&(py::method(&v_neg.clone(), "items", &[])?))? {
        let t5520 = (t5519).items()?;
        if t5520.len() != 2 { return Err(PyErr::Value("unpack".into())); }
        v_aa = t5520[0].clone();
        v_pk = t5520[1].clone();
        v_minus = py::iadd(&v_minus, &py::mul(&py::getitem(&v_n.clone(), &v_aa.clone())?, &py::truediv(&Py::Float(1.0), &py::add(&py::pow(&Py::Int(10i64), &py::sub(&v_pk.clone(), &v_pH.clone())?)?, &Py::Float(1.0))?)?)?)?;
    }
    return Ok(py::sub(&v_plus.clone(), &v_minus.clone())?);
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn f__check_align(mut v_a: Py, mut v_b: Py, mut v_scores: Py) -> R {
    let mut v_name: Py = Py::None;
    let mut v_v: Py = Py::None;
    let _ = f__need({ let t6240 = py::mul(&py::len(&v_a.clone())?, &py::len(&v_b.clone())?)?; Py::Bool({ let t6241 = Py::Int(1000000i64); if py::le(&t6240, &t6241)? { true } else { false } }) }, Py::str("sequences too long (length product at most 1 000 000)"))?;
    for t6242 in py::iter(&(py::method(&v_scores.clone(), "items", &[])?))? {
        let t6246 = (t6242).items()?;
        if t6246.len() != 2 { return Err(PyErr::Value("unpack".into())); }
        v_name = t6246[0].clone();
        v_v = t6246[1].clone();
        let _ = f__fin(v_name.clone(), Py::tuple(vec![v_v.clone()]))?;
    }
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn f__ci_model(mut v_S: Py, mut v_I: Py) -> R {
    let mut v_m: Py = Py::None;
    v_m = { let v_I = v_I.clone(); let v_S = v_S.clone(); let v_m = v_m.clone(); py::func(move |a: &[Py]| -> R {
            #[allow(unused_mut, unused_variables, unused_assignments)]
{
            let mut v_p: Py = a.get(0).cloned().unwrap_or(Py::None);
            let mut v_J: Py = Py::None;
            let mut v_Ki: Py = Py::None;
            let mut v_Km: Py = Py::None;
            let mut v_Vm: Py = Py::None;
            let mut v_den: Py = Py::None;
            let mut v_f: Py = Py::None;
            let mut v_i: Py = Py::None;
            let mut v_s: Py = Py::None;
            let t691 = (v_p.clone()).items()?;
            if t691.len() != 3 { return Err(PyErr::Value("unpack".into())); }
            v_Vm = t691[0].clone();
            v_Km = t691[1].clone();
            v_Ki = t691[2].clone();
            let t692 = (Py::tuple(vec![Py::list(vec![]), Py::list(vec![])])).items()?;
            if t692.len() != 2 { return Err(PyErr::Value("unpack".into())); }
            v_f = t692[0].clone();
            v_J = t692[1].clone();
            for t693 in py::iter(&(py::zip(&[v_S.clone(), v_I.clone()])?))? {
                let t694 = (t693).items()?;
                if t694.len() != 2 { return Err(PyErr::Value("unpack".into())); }
                v_s = t694[0].clone();
                v_i = t694[1].clone();
                v_den = py::add(&py::mul(&v_Km.clone(), &py::add(&Py::Int(1i64), &py::truediv(&v_i.clone(), &v_Ki.clone())?)?)?, &v_s.clone())?;
                let _ = py::method(&v_f.clone(), "append", &[py::truediv(&py::mul(&v_Vm.clone(), &v_s.clone())?, &v_den.clone())?])?;
                let _ = py::method(&v_J.clone(), "append", &[Py::list(vec![py::truediv(&v_s.clone(), &v_den.clone())?, py::truediv(&py::mul(&py::mul(&py::neg(&v_Vm.clone())?, &v_s.clone())?, &py::add(&Py::Int(1i64), &py::truediv(&v_i.clone(), &v_Ki.clone())?)?)?, &py::pow(&v_den.clone(), &Py::Int(2i64))?)?, py::truediv(&py::truediv(&py::mul(&py::mul(&py::mul(&v_Vm.clone(), &v_s.clone())?, &v_Km.clone())?, &v_i.clone())?, &py::mul(&v_Ki.clone(), &v_Ki.clone())?)?, &py::pow(&v_den.clone(), &Py::Int(2i64))?)?])])?;
            }
            return Ok(Py::tuple(vec![v_f.clone(), v_J.clone()]));
            #[allow(unreachable_code)]
            Ok(Py::None)
        }
        }) };
    return Ok(v_m.clone());
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn f__counts(mut v_name: Py, mut v_v: Py, mut v_lo: Py, mut v_positive: Py) -> R {
    let mut v_out: Py = Py::None;
    let _ = f__need({ let t6415 = Py::Bool(matches!(v_v.clone(), Py::List(_) | Py::Tuple(_))); if !t6415.truthy() { t6415 } else { { let t6413 = py::len(&v_v.clone())?; Py::Bool({ let t6414 = Py::Int(1i64); if py::ge(&t6413, &t6414)? { true } else { false } }) } } }, Py::str(&[py::format(&v_name.clone(), "")?, " must be a list of counts".to_string()].concat()))?;
    v_out = { let mut out = vec![]; for t6416 in py::iter(&(v_v.clone()))? { let mut v_x: Py; 
v_x = t6416;
out.push(f__whole(v_name.clone(), v_x.clone(), v_lo.clone())?); } Py::list(out) };
    if (v_positive.clone()).truthy() {
        let _ = f__need(Py::Bool(({ let mut out = vec![]; for t6425 in py::iter(&(v_out.clone()))? { let mut v_x: Py; 
v_x = t6425;
out.push({ let t6426 = v_x.clone(); Py::Bool({ let t6427 = Py::Int(0i64); if py::gt(&t6426, &t6427)? { true } else { false } }) }); } Py::list(out) }).items()?.iter().all(|x| x.truthy())), Py::str(&[py::format(&v_name.clone(), "")?, ": every count must be positive".to_string()].concat()))?;
    }
    return Ok(v_out.clone());
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn f__dfact(mut v_k: Py) -> R {
    let mut v_out: Py = Py::None;
    v_out = Py::Int(1i64);
    while ({ let t6962 = v_k.clone(); Py::Bool({ let t6963 = Py::Int(1i64); if py::gt(&t6962, &t6963)? { true } else { false } }) }).truthy() {
        v_out = py::mul(&v_out, &v_k.clone())?;
        v_k = py::sub(&v_k, &Py::Int(2i64))?;
    }
    return Ok(v_out.clone());
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn f__distance_matrix(mut v_D: Py, mut v_name: Py, mut v_lo: Py, mut v_hi: Py) -> R {
    let mut v_i: Py = Py::None;
    let mut v_j: Py = Py::None;
    let mut v_n: Py = Py::None;
    v_n = py::len(&v_D.clone())?;
    let _ = f__need({ let t7011 = { let t7005 = v_lo.clone(); Py::Bool({ let t7006 = v_n.clone(); if py::le(&t7005, &t7006)? { { let t7007 = v_hi.clone(); if py::le(&t7006, &t7007)? { true } else { false } } } else { false } }) }; if !t7011.truthy() { t7011 } else { Py::Bool(({ let mut out = vec![]; for t7008 in py::iter(&(v_D.clone()))? { let mut v_r: Py; 
v_r = t7008;
out.push({ let t7009 = py::len(&v_r.clone())?; Py::Bool({ let t7010 = v_n.clone(); if py::eq(&t7009, &t7010) { true } else { false } }) }); } Py::list(out) }).items()?.iter().all(|x| x.truthy())) } }, Py::str(&[py::format(&v_name.clone(), "")?, " must be a square matrix of ".to_string(), py::format(&v_lo.clone(), "")?, " to ".to_string(), py::format(&v_hi.clone(), "")?, " taxa".to_string()].concat()))?;
    for t7012 in py::range_iter(&[(v_n.clone()).i()?]) {
        v_i = t7012;
        let _ = f__need({ let t7013 = py::getitem(&py::getitem(&v_D.clone(), &v_i.clone())?, &v_i.clone())?; Py::Bool({ let t7014 = Py::Int(0i64); if py::eq(&t7013, &t7014) { true } else { false } }) }, Py::str(&[py::format(&v_name.clone(), "")?, ": the diagonal must be 0".to_string()].concat()))?;
        for t7015 in py::range_iter(&[(py::add(&v_i.clone(), &Py::Int(1i64))?).i()?, (v_n.clone()).i()?]) {
            v_j = t7015;
            let _ = f__need({ let t7016 = py::getitem(&py::getitem(&v_D.clone(), &v_i.clone())?, &v_j.clone())?; Py::Bool({ let t7017 = py::getitem(&py::getitem(&v_D.clone(), &v_j.clone())?, &v_i.clone())?; if py::eq(&t7016, &t7017) { true } else { false } }) }, Py::str(&[py::format(&v_name.clone(), "")?, " must be symmetric".to_string()].concat()))?;
            let _ = f__need({ let t7018 = py::getitem(&py::getitem(&v_D.clone(), &v_i.clone())?, &v_j.clone())?; Py::Bool({ let t7019 = Py::Int(0i64); if py::gt(&t7018, &t7019)? { true } else { false } }) }, Py::str(&[py::format(&v_name.clone(), "")?, ": distances between different taxa must be positive".to_string()].concat()))?;
        }
    }
    return Ok({ let mut out = vec![]; for t7020 in py::iter(&(v_D.clone()))? { let mut v_r: Py; 
v_r = t7020;
out.push({ let mut out = vec![]; for t7021 in py::iter(&(v_r.clone()))? { let mut v_x: Py; 
v_x = t7021;
out.push(py::float(&v_x.clone())?); } Py::list(out) }); } Py::list(out) });
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn f__euler_lotka(mut v_F: Py, mut v_S: Py) -> R {
    let mut v_l: Py = Py::None;
    let mut v_s: Py = Py::None;
    v_l = Py::list(vec![Py::Float(1.0)]);
    for t3385 in py::iter(&(v_S.clone()))? {
        v_s = t3385;
        let _ = py::method(&v_l.clone(), "append", &[py::mul(&py::getitem(&v_l.clone(), &py::neg(&Py::Int(1i64))?)?, &v_s.clone())?])?;
    }
    return Ok(v_l.clone());
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn f__fin(mut v_name: Py, mut v_vals: Py) -> R {
    let mut v_v: Py = Py::None;
    for t6670 in py::iter(&(v_vals.clone()))? {
        v_v = t6670;
        let _ = f__need({ let t6672 = Py::Bool(matches!(v_v.clone(), Py::Int(_) | Py::Big(_) | Py::Bool(_) | Py::Float(_))); if !t6672.truthy() { t6672 } else { { let t6671 = py::not(&Py::Bool(matches!(v_v.clone(), Py::Bool(_)))); if !t6671.truthy() { t6671 } else { py::math1("isfinite", &v_v.clone())? } } } }, Py::str(&[py::format(&v_name.clone(), "")?, " must be a finite number".to_string()].concat()))?;
    }
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn f__fr_text(mut v_f: Py) -> R {
    return Ok(Py::str(&[py::format(&py::numerator(&v_f.clone())?, "")?, "/".to_string(), py::format(&py::denominator(&v_f.clone())?, "")?].concat()));
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn f__frac(mut v_name: Py, mut v_v: Py, mut v_lo_open: Py, mut v_hi_open: Py) -> R {
    let _ = f__fin(v_name.clone(), Py::tuple(vec![v_v.clone()]))?;
    let _ = f__need({ let t5640 = (if (v_lo_open.clone()).truthy() { { let t5632 = v_v.clone(); Py::Bool({ let t5633 = Py::Int(0i64); if py::gt(&t5632, &t5633)? { true } else { false } }) } } else { { let t5634 = v_v.clone(); Py::Bool({ let t5635 = Py::Int(0i64); if py::ge(&t5634, &t5635)? { true } else { false } }) } }); if !t5640.truthy() { t5640 } else { (if (v_hi_open.clone()).truthy() { { let t5636 = v_v.clone(); Py::Bool({ let t5637 = Py::Int(1i64); if py::lt(&t5636, &t5637)? { true } else { false } }) } } else { { let t5638 = v_v.clone(); Py::Bool({ let t5639 = Py::Int(1i64); if py::le(&t5638, &t5639)? { true } else { false } }) } }) } }, Py::str(&[py::format(&v_name.clone(), "")?, " must be a fraction between 0 and 1".to_string()].concat()))?;
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn f__genotype(mut v_g: Py, mut v_name: Py) -> R {
    let mut v_i: Py = Py::None;
    let mut v_letters: Py = Py::None;
    let mut v_loci: Py = Py::None;
    let mut v_x: Py = Py::None;
    let mut v_y: Py = Py::None;
    let _ = f__need({ let t894 = Py::Bool(matches!(v_g.clone(), Py::Str(_))); if !t894.truthy() { t894 } else { { let t893 = { let t888 = ax_engine::pyre::module_call("fullmatch", &Py::str("[A-Za-z]+"), &[{ let t890 = v_g.clone(); if t890.truthy() { t890 } else { Py::str("") } }], 0)?; Py::Bool({ let t889 = Py::None; if !py::is(&t888, &t889) { true } else { false } }) }; if !t893.truthy() { t893 } else { { let t891 = py::modulo(&py::len(&v_g.clone())?, &Py::Int(2i64))?; Py::Bool({ let t892 = Py::Int(0i64); if py::eq(&t891, &t892) { true } else { false } }) } } } } }, Py::str(&[py::format(&v_name.clone(), "")?, " must be letter pairs, one pair per locus (e.g. AaBb)".to_string()].concat()))?;
    v_loci = Py::list(vec![]);
    for t895 in py::range_iter(&[(Py::Int(0i64)).i()?, (py::len(&v_g.clone())?).i()?, (Py::Int(2i64)).i()?]) {
        v_i = t895;
        let t896 = (Py::tuple(vec![py::getitem(&v_g.clone(), &v_i.clone())?, py::getitem(&v_g.clone(), &py::add(&v_i.clone(), &Py::Int(1i64))?)?])).items()?;
        if t896.len() != 2 { return Err(PyErr::Value("unpack".into())); }
        v_x = t896[0].clone();
        v_y = t896[1].clone();
        let _ = f__need({ let t897 = py::method(&v_x.clone(), "upper", &[])?; Py::Bool({ let t898 = py::method(&v_y.clone(), "upper", &[])?; if py::eq(&t897, &t898) { true } else { false } }) }, Py::str(&[py::format(&v_name.clone(), "")?, ": '".to_string(), py::format(&v_x.clone(), "")?, py::format(&v_y.clone(), "")?, "' is not one locus (both letters must name the same gene)".to_string()].concat()))?;
        let _ = py::method(&v_loci.clone(), "append", &[Py::tuple((py::sorted_key(&Py::tuple(vec![v_x.clone(), v_y.clone()]), &{ let v_g = v_g.clone(); let v_i = v_i.clone(); let v_loci = v_loci.clone(); let v_name = v_name.clone(); let v_x = v_x.clone(); let v_y = v_y.clone(); py::func(move |a: &[Py]| -> R {
            #[allow(unused_mut, unused_variables, unused_assignments)]
{
            let mut v_c: Py = a.get(0).cloned().unwrap_or(Py::None);
            return Ok(Py::tuple(vec![py::method(&v_c.clone(), "islower", &[])?, v_c.clone()]));
            #[allow(unreachable_code)]
            Ok(Py::None)
        }
        }) }, (Py::Bool(false)).truthy())?).items()?)])?;
    }
    v_letters = { let mut out = vec![]; for t899 in py::iter(&(v_loci.clone()))? { let mut v_p: Py; 
v_p = t899;
out.push(py::method(&py::getitem(&v_p.clone(), &Py::Int(0i64))?, "upper", &[])?); } Py::list(out) };
    let _ = f__need({ let t900 = py::len(&py::set(&v_letters.clone())?)?; Py::Bool({ let t901 = py::len(&v_letters.clone())?; if py::eq(&t900, &t901) { true } else { false } }) }, Py::str(&[py::format(&v_name.clone(), "")?, ": each locus must use a different letter".to_string()].concat()))?;
    let _ = f__need({ let t902 = py::len(&v_loci.clone())?; Py::Bool({ let t903 = Py::Int(6i64); if py::le(&t902, &t903)? { true } else { false } }) }, Py::str(&[py::format(&v_name.clone(), "")?, ": at most 6 loci".to_string()].concat()))?;
    return Ok(v_loci.clone());
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn f__gn(mut v_model: Py, mut v_p: Py, mut v_ys: Py, mut v_ok: Py) -> R {
    let mut v_A: Py = Py::None;
    let mut v_J: Py = Py::None;
    let mut v__: Py = Py::None;
    let mut v_d: Py = Py::None;
    let mut v_f: Py = Py::None;
    let mut v_g: Py = Py::None;
    let mut v_last: Py = Py::None;
    let mut v_q: Py = Py::None;
    let mut v_r: Py = Py::None;
    let mut v_size: Py = Py::None;
    v_last = Py::Float(f64::INFINITY);
    for t2721 in py::range_iter(&[(Py::Int(60i64)).i()?]) {
        v__ = t2721;
        let t2722 = (py::call(&v_model.clone(), &[v_p.clone()])?).items()?;
        if t2722.len() != 2 { return Err(PyErr::Value("unpack".into())); }
        v_f = t2722[0].clone();
        v_J = t2722[1].clone();
        v_r = { let mut out = vec![]; for t2723 in py::iter(&(py::zip(&[v_ys.clone(), v_f.clone()])?))? { let mut v_a: Py; let mut v_y: Py; 
let t2724 = (t2723).items()?;
if t2724.len() != 2 { return Err(PyErr::Value("unpack".into())); }
v_y = t2724[0].clone();
v_a = t2724[1].clone();
out.push(py::sub(&v_y.clone(), &v_a.clone())?); } Py::list(out) };
        v_A = { let mut out = vec![]; for t2725 in py::range_iter(&[(Py::Int(3i64)).i()?]) { let mut v_i: Py; 
v_i = t2725;
out.push({ let mut out = vec![]; for t2726 in py::range_iter(&[(Py::Int(3i64)).i()?]) { let mut v_j: Py; 
v_j = t2726;
out.push(py::fsum(&{ let mut out = vec![]; for t2727 in py::iter(&(v_J.clone()))? { let mut v_row: Py; 
v_row = t2727;
out.push(py::mul(&py::getitem(&v_row.clone(), &v_i.clone())?, &py::getitem(&v_row.clone(), &v_j.clone())?)?); } Py::list(out) })?); } Py::list(out) }); } Py::list(out) };
        v_g = { let mut out = vec![]; for t2728 in py::range_iter(&[(Py::Int(3i64)).i()?]) { let mut v_i: Py; 
v_i = t2728;
out.push(py::fsum(&{ let mut out = vec![]; for t2729 in py::iter(&(py::zip(&[v_J.clone(), v_r.clone()])?))? { let mut v_ri: Py; let mut v_row: Py; 
let t2730 = (t2729).items()?;
if t2730.len() != 2 { return Err(PyErr::Value("unpack".into())); }
v_row = t2730[0].clone();
v_ri = t2730[1].clone();
out.push(py::mul(&py::getitem(&v_row.clone(), &v_i.clone())?, &v_ri.clone())?); } Py::list(out) })?); } Py::list(out) };
        let t2731: Result<Option<Py>, PyErr> = (|| -> Result<Option<Py>, PyErr> {
                v_d = f__solve3(v_A.clone(), v_g.clone())?;
                #[allow(unreachable_code)]
                Ok(None)
            })();
        match t2731 {
            Ok(Some(r)) => { return Ok(r); }
            Ok(None) => {}
            Err(e) if (matches!(e, PyErr::Value(_))) => {
                return Ok(v_p.clone());
            }
            Err(e) => return Err(e),
        }
        v_q = { let mut out = vec![]; for t2732 in py::range_iter(&[(Py::Int(3i64)).i()?]) { let mut v_i: Py; 
v_i = t2732;
out.push(py::add(&py::getitem(&v_p.clone(), &v_i.clone())?, &py::getitem(&v_d.clone(), &v_i.clone())?)?); } Py::list(out) };
        if (py::not(&py::call(&v_ok.clone(), &[v_q.clone()])?)).truthy() {
            return Ok(v_p.clone());
        }
        v_size = py::minmax(&[{ let mut out = vec![]; for t2733 in py::range_iter(&[(Py::Int(3i64)).i()?]) { let mut v_i: Py; 
v_i = t2733;
out.push(py::truediv(&py::abs(&py::getitem(&v_d.clone(), &v_i.clone())?)?, &py::minmax(&[py::abs(&py::getitem(&v_q.clone(), &v_i.clone())?)?, Py::Float(1e-300)], true)?)?); } Py::list(out) }], true)?;
        v_p = v_q.clone();
        if ({ let t2738 = { let t2734 = v_size.clone(); Py::Bool({ let t2735 = Py::Float(1e-16); if py::lt(&t2734, &t2735)? { true } else { false } }) }; if t2738.truthy() { t2738 } else { { let t2736 = v_size.clone(); Py::Bool({ let t2737 = v_last.clone(); if py::ge(&t2736, &t2737)? { true } else { false } }) } } }).truthy() {
            return Ok(v_p.clone());
        }
        v_last = v_size.clone();
    }
    return Ok(v_p.clone());
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn f__has_any(mut v_t: Py, mut v_want: Py) -> R {
    if (py::not(&py::getitem(&v_t.clone(), &Py::Int(2i64))?)).truthy() {
        return Ok({ let t1341 = py::getitem(&v_t.clone(), &Py::Int(0i64))?; Py::Bool({ let t1342 = v_want.clone(); if py::contains(&t1342, &t1341)? { true } else { false } }) });
    }
    return Ok(Py::Bool(({ let mut out = vec![]; for t1343 in py::iter(&(py::getitem(&v_t.clone(), &Py::Int(2i64))?))? { let mut v_k: Py; 
v_k = t1343;
out.push(f__has_any(v_k.clone(), v_want.clone())?); } Py::list(out) }).items()?.iter().any(|x| x.truthy())));
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn f__hill_model(mut v_S: Py) -> R {
    let mut v_m: Py = Py::None;
    v_m = { let v_S = v_S.clone(); let v_m = v_m.clone(); py::func(move |a: &[Py]| -> R {
            #[allow(unused_mut, unused_variables, unused_assignments)]
{
            let mut v_p: Py = a.get(0).cloned().unwrap_or(Py::None);
            let mut v_J: Py = Py::None;
            let mut v_K: Py = Py::None;
            let mut v_Vm: Py = Py::None;
            let mut v_dx: Py = Py::None;
            let mut v_f: Py = Py::None;
            let mut v_n: Py = Py::None;
            let mut v_s: Py = Py::None;
            let mut v_u: Py = Py::None;
            let mut v_x: Py = Py::None;
            let t2741 = (v_p.clone()).items()?;
            if t2741.len() != 3 { return Err(PyErr::Value("unpack".into())); }
            v_Vm = t2741[0].clone();
            v_K = t2741[1].clone();
            v_n = t2741[2].clone();
            let t2742 = (Py::tuple(vec![Py::list(vec![]), Py::list(vec![])])).items()?;
            if t2742.len() != 2 { return Err(PyErr::Value("unpack".into())); }
            v_f = t2742[0].clone();
            v_J = t2742[1].clone();
            for t2743 in py::iter(&(v_S.clone()))? {
                v_s = t2743;
                v_u = py::pow(&py::truediv(&v_s.clone(), &v_K.clone())?, &v_n.clone())?;
                v_x = py::truediv(&v_u.clone(), &py::add(&Py::Int(1i64), &v_u.clone())?)?;
                let _ = py::method(&v_f.clone(), "append", &[py::mul(&v_Vm.clone(), &v_x.clone())?])?;
                v_dx = py::mul(&v_x.clone(), &py::sub(&Py::Int(1i64), &v_x.clone())?)?;
                let _ = py::method(&v_J.clone(), "append", &[Py::list(vec![v_x.clone(), py::truediv(&py::mul(&py::mul(&py::neg(&v_Vm.clone())?, &v_dx.clone())?, &v_n.clone())?, &v_K.clone())?, py::mul(&py::mul(&v_Vm.clone(), &v_dx.clone())?, &py::math1("log", &py::truediv(&v_s.clone(), &v_K.clone())?)?)?])])?;
            }
            return Ok(Py::tuple(vec![v_f.clone(), v_J.clone()]));
            #[allow(unreachable_code)]
            Ok(Py::None)
        }
        }) };
    return Ok(v_m.clone());
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn f__hwe_probs(mut v_nA: Py, mut v_N: Py) -> R {
    let mut v_k: Py = Py::None;
    let mut v_nAA: Py = Py::None;
    let mut v_nB: Py = Py::None;
    let mut v_nBB: Py = Py::None;
    let mut v_out: Py = Py::None;
    let mut v_total: Py = Py::None;
    let mut v_ways: Py = Py::None;
    v_nB = py::sub(&py::mul(&Py::Int(2i64), &v_N.clone())?, &v_nA.clone())?;
    v_total = py::comb(&py::mul(&Py::Int(2i64), &v_N.clone())?, &v_nA.clone())?;
    v_out = Py::dict(vec![]);
    for t2409 in py::range_iter(&[(py::modulo(&v_nA.clone(), &Py::Int(2i64))?).i()?, (py::add(&py::minmax(&[v_nA.clone(), v_nB.clone()], false)?, &Py::Int(1i64))?).i()?, (Py::Int(2i64)).i()?]) {
        v_k = t2409;
        let t2410 = (Py::tuple(vec![py::floordiv(&py::sub(&v_nA.clone(), &v_k.clone())?, &Py::Int(2i64))?, py::floordiv(&py::sub(&v_nB.clone(), &v_k.clone())?, &Py::Int(2i64))?])).items()?;
        if t2410.len() != 2 { return Err(PyErr::Value("unpack".into())); }
        v_nAA = t2410[0].clone();
        v_nBB = t2410[1].clone();
        v_ways = py::mul(&py::floordiv(&py::factorial(&v_N.clone())?, &py::mul(&py::mul(&py::factorial(&v_nAA.clone())?, &py::factorial(&v_k.clone())?)?, &py::factorial(&v_nBB.clone())?)?)?, &py::pow(&Py::Int(2i64), &v_k.clone())?)?;
        py::setitem(&v_out.clone(), &v_k.clone(), py::fraction(&[v_ways.clone(), v_total.clone()])?)?;
    }
    return Ok(v_out.clone());
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn f__kin_data(mut v_S: Py, mut v_v: Py, mut v_lo: Py) -> R {
    let mut v_x: Py = Py::None;
    let _ = f__need({ let t4109 = Py::Bool(matches!(v_S.clone(), Py::List(_))); if !t4109.truthy() { t4109 } else { { let t4108 = Py::Bool(matches!(v_v.clone(), Py::List(_))); if !t4108.truthy() { t4108 } else { { let t4107 = { let t4103 = py::len(&v_S.clone())?; Py::Bool({ let t4104 = py::len(&v_v.clone())?; if py::eq(&t4103, &t4104) { true } else { false } }) }; if !t4107.truthy() { t4107 } else { { let t4105 = py::len(&v_S.clone())?; Py::Bool({ let t4106 = v_lo.clone(); if py::ge(&t4105, &t4106)? { true } else { false } }) } } } } } } }, Py::str(&["S and v must be lists of equal length (at least ".to_string(), py::format(&v_lo.clone(), "")?, " points)".to_string()].concat()))?;
    for t4110 in py::iter(&(py::add(&v_S.clone(), &v_v.clone())?))? {
        v_x = t4110;
        let _ = f__pos(Py::str("every concentration and rate"), Py::tuple(vec![v_x.clone()]))?;
    }
    return Ok(Py::tuple(vec![{ let mut out = vec![]; for t4117 in py::iter(&(v_S.clone()))? { let mut v_x: Py; 
v_x = t4117;
out.push(py::float(&v_x.clone())?); } Py::list(out) }, { let mut out = vec![]; for t4118 in py::iter(&(v_v.clone()))? { let mut v_x: Py; 
v_x = t4118;
out.push(py::float(&v_x.clone())?); } Py::list(out) }]));
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn f__leslie_in(mut v_F: Py, mut v_S: Py) -> R {
    let mut v_f: Py = Py::None;
    let mut v_s: Py = Py::None;
    let _ = f__need({ let t3485 = Py::Bool(matches!(v_F.clone(), Py::List(_))); if !t3485.truthy() { t3485 } else { { let t3484 = Py::Bool(matches!(v_S.clone(), Py::List(_))); if !t3484.truthy() { t3484 } else { { let t3483 = { let t3479 = py::len(&v_F.clone())?; Py::Bool({ let t3480 = Py::Int(2i64); if py::ge(&t3479, &t3480)? { true } else { false } }) }; if !t3483.truthy() { t3483 } else { { let t3481 = py::len(&v_S.clone())?; Py::Bool({ let t3482 = py::sub(&py::len(&v_F.clone())?, &Py::Int(1i64))?; if py::eq(&t3481, &t3482) { true } else { false } }) } } } } } } }, Py::str("F needs one fecundity per age class (at least 2) and S one survival fewer"))?;
    for t3486 in py::iter(&(v_F.clone()))? {
        v_f = t3486;
        let _ = f__fin(Py::str("F"), Py::tuple(vec![v_f.clone()]))?;
        let _ = f__need({ let t3490 = v_f.clone(); Py::Bool({ let t3491 = Py::Int(0i64); if py::ge(&t3490, &t3491)? { true } else { false } }) }, Py::str("fecundities must be zero or positive"))?;
    }
    for t3492 in py::iter(&(v_S.clone()))? {
        v_s = t3492;
        let _ = f__fin(Py::str("S"), Py::tuple(vec![v_s.clone()]))?;
        let _ = f__need({ let t3493 = Py::Int(0i64); Py::Bool({ let t3494 = v_s.clone(); if py::lt(&t3493, &t3494)? { { let t3495 = Py::Int(1i64); if py::le(&t3494, &t3495)? { true } else { false } } } else { false } }) }, Py::str("survival probabilities must lie in (0, 1]"))?;
    }
    let _ = f__need(Py::Bool(({ let mut out = vec![]; for t3496 in py::iter(&(v_F.clone()))? { let mut v_f: Py; 
v_f = t3496;
out.push({ let t3497 = v_f.clone(); Py::Bool({ let t3498 = Py::Int(0i64); if py::gt(&t3497, &t3498)? { true } else { false } }) }); } Py::list(out) }).items()?.iter().any(|x| x.truthy())), Py::str("at least one age class must reproduce"))?;
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn f__levenshtein(mut v_a: Py, mut v_b: Py) -> R {
    let mut v_cur: Py = Py::None;
    let mut v_i: Py = Py::None;
    let mut v_j: Py = Py::None;
    let mut v_prev: Py = Py::None;
    v_prev = Py::list((Py::list((py::range_iter(&[(py::add(&py::len(&v_b.clone())?, &Py::Int(1i64))?).i()?])).into_iter().collect::<Vec<Py>>())).items()?);
    for t1237 in py::range_iter(&[(Py::Int(1i64)).i()?, (py::add(&py::len(&v_a.clone())?, &Py::Int(1i64))?).i()?]) {
        v_i = t1237;
        v_cur = py::add(&Py::list(vec![v_i.clone()]), &py::mul(&Py::list(vec![Py::Int(0i64)]), &py::len(&v_b.clone())?)?)?;
        for t1238 in py::range_iter(&[(Py::Int(1i64)).i()?, (py::add(&py::len(&v_b.clone())?, &Py::Int(1i64))?).i()?]) {
            v_j = t1238;
            py::setitem(&v_cur.clone(), &v_j.clone(), py::minmax(&[py::add(&py::getitem(&v_prev.clone(), &v_j.clone())?, &Py::Int(1i64))?, py::add(&py::getitem(&v_cur.clone(), &py::sub(&v_j.clone(), &Py::Int(1i64))?)?, &Py::Int(1i64))?, py::add(&py::getitem(&v_prev.clone(), &py::sub(&v_j.clone(), &Py::Int(1i64))?)?, &{ let t1239 = py::getitem(&v_a.clone(), &py::sub(&v_i.clone(), &Py::Int(1i64))?)?; Py::Bool({ let t1240 = py::getitem(&v_b.clone(), &py::sub(&v_j.clone(), &Py::Int(1i64))?)?; if py::ne(&t1239, &t1240) { true } else { false } }) })?], false)?)?;
        }
        v_prev = v_cur.clone();
    }
    return Ok(py::getitem(&v_prev.clone(), &py::neg(&Py::Int(1i64))?)?);
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn f__lm(mut v_model: Py, mut v_p: Py, mut v_ys: Py, mut v_ok: Py) -> R {
    let mut v_A: Py = Py::None;
    let mut v_Ad: Py = Py::None;
    let mut v_J: Py = Py::None;
    let mut v__: Py = Py::None;
    let mut v_c: Py = Py::None;
    let mut v_cost: Py = Py::None;
    let mut v_cq: Py = Py::None;
    let mut v_d: Py = Py::None;
    let mut v_f: Py = Py::None;
    let mut v_g: Py = Py::None;
    let mut v_lam: Py = Py::None;
    let mut v_q: Py = Py::None;
    let mut v_r: Py = Py::None;
    let mut v_small: Py = Py::None;
    v_cost = { let v_cost = v_cost.clone(); let v_model = v_model.clone(); let v_ok = v_ok.clone(); let v_p = v_p.clone(); let v_ys = v_ys.clone(); py::func(move |a: &[Py]| -> R {
            #[allow(unused_mut, unused_variables, unused_assignments)]
{
            let mut v_pp: Py = a.get(0).cloned().unwrap_or(Py::None);
            let mut v__: Py = Py::None;
            let mut v_f: Py = Py::None;
            let t2683 = (py::call(&v_model.clone(), &[v_pp.clone()])?).items()?;
            if t2683.len() != 2 { return Err(PyErr::Value("unpack".into())); }
            v_f = t2683[0].clone();
            v__ = t2683[1].clone();
            return Ok(py::fsum(&{ let mut out = vec![]; for t2684 in py::iter(&(py::zip(&[v_ys.clone(), v_f.clone()])?))? { let mut v_a: Py; let mut v_y: Py; 
let t2685 = (t2684).items()?;
if t2685.len() != 2 { return Err(PyErr::Value("unpack".into())); }
v_y = t2685[0].clone();
v_a = t2685[1].clone();
out.push(py::pow(&py::sub(&v_y.clone(), &v_a.clone())?, &Py::Int(2i64))?); } Py::list(out) })?);
            #[allow(unreachable_code)]
            Ok(Py::None)
        }
        }) };
    let t2686 = (Py::tuple(vec![Py::Float(0.001), py::call(&v_cost.clone(), &[v_p.clone()])?])).items()?;
    if t2686.len() != 2 { return Err(PyErr::Value("unpack".into())); }
    v_lam = t2686[0].clone();
    v_c = t2686[1].clone();
    for t2687 in py::range_iter(&[(Py::Int(2000i64)).i()?]) {
        v__ = t2687;
        let t2688 = (py::call(&v_model.clone(), &[v_p.clone()])?).items()?;
        if t2688.len() != 2 { return Err(PyErr::Value("unpack".into())); }
        v_f = t2688[0].clone();
        v_J = t2688[1].clone();
        v_r = { let mut out = vec![]; for t2689 in py::iter(&(py::zip(&[v_ys.clone(), v_f.clone()])?))? { let mut v_a: Py; let mut v_y: Py; 
let t2690 = (t2689).items()?;
if t2690.len() != 2 { return Err(PyErr::Value("unpack".into())); }
v_y = t2690[0].clone();
v_a = t2690[1].clone();
out.push(py::sub(&v_y.clone(), &v_a.clone())?); } Py::list(out) };
        v_A = { let mut out = vec![]; for t2691 in py::range_iter(&[(Py::Int(3i64)).i()?]) { let mut v_i: Py; 
v_i = t2691;
out.push({ let mut out = vec![]; for t2692 in py::range_iter(&[(Py::Int(3i64)).i()?]) { let mut v_j: Py; 
v_j = t2692;
out.push(py::fsum(&{ let mut out = vec![]; for t2693 in py::iter(&(v_J.clone()))? { let mut v_row: Py; 
v_row = t2693;
out.push(py::mul(&py::getitem(&v_row.clone(), &v_i.clone())?, &py::getitem(&v_row.clone(), &v_j.clone())?)?); } Py::list(out) })?); } Py::list(out) }); } Py::list(out) };
        v_g = { let mut out = vec![]; for t2694 in py::range_iter(&[(Py::Int(3i64)).i()?]) { let mut v_i: Py; 
v_i = t2694;
out.push(py::fsum(&{ let mut out = vec![]; for t2695 in py::iter(&(py::zip(&[v_J.clone(), v_r.clone()])?))? { let mut v_ri: Py; let mut v_row: Py; 
let t2696 = (t2695).items()?;
if t2696.len() != 2 { return Err(PyErr::Value("unpack".into())); }
v_row = t2696[0].clone();
v_ri = t2696[1].clone();
out.push(py::mul(&py::getitem(&v_row.clone(), &v_i.clone())?, &v_ri.clone())?); } Py::list(out) })?); } Py::list(out) };
        while (Py::Bool(true)).truthy() {
            v_Ad = { let mut out = vec![]; for t2697 in py::range_iter(&[(Py::Int(3i64)).i()?]) { let mut v_i: Py; 
v_i = t2697;
out.push({ let mut out = vec![]; for t2698 in py::range_iter(&[(Py::Int(3i64)).i()?]) { let mut v_j: Py; 
v_j = t2698;
out.push(py::add(&py::getitem(&py::getitem(&v_A.clone(), &v_i.clone())?, &v_j.clone())?, &(if ({ let t2699 = v_i.clone(); Py::Bool({ let t2700 = v_j.clone(); if py::eq(&t2699, &t2700) { true } else { false } }) }).truthy() { py::mul(&v_lam.clone(), &py::getitem(&py::getitem(&v_A.clone(), &v_i.clone())?, &v_i.clone())?)? } else { Py::Float(0.0) }))?); } Py::list(out) }); } Py::list(out) };
            let t2712: Result<Option<Py>, PyErr> = (|| -> Result<Option<Py>, PyErr> {
                    v_d = f__solve3(v_Ad.clone(), v_g.clone())?;
                    #[allow(unreachable_code)]
                    Ok(None)
                })();
            match t2712 {
                Ok(Some(r)) => { return Ok(r); }
                Ok(None) => {}
                Err(e) if (matches!(e, PyErr::Value(_))) => {
                    v_d = Py::None;
                }
                Err(e) => return Err(e),
            }
            v_q = (if (v_d.clone()).truthy() { { let mut out = vec![]; for t2713 in py::range_iter(&[(Py::Int(3i64)).i()?]) { let mut v_i: Py; 
v_i = t2713;
out.push(py::add(&py::getitem(&v_p.clone(), &v_i.clone())?, &py::getitem(&v_d.clone(), &v_i.clone())?)?); } Py::list(out) } } else { Py::None });
            if ({ let t2714 = v_q.clone(); if !t2714.truthy() { t2714 } else { py::call(&v_ok.clone(), &[v_q.clone()])? } }).truthy() {
                v_cq = py::call(&v_cost.clone(), &[v_q.clone()])?;
                if ({ let t2715 = v_cq.clone(); Py::Bool({ let t2716 = v_c.clone(); if py::lt(&t2715, &t2716)? { true } else { false } }) }).truthy() {
                    v_small = { let t2717 = py::minmax(&[{ let mut out = vec![]; for t2719 in py::range_iter(&[(Py::Int(3i64)).i()?]) { let mut v_i: Py; 
v_i = t2719;
out.push(py::truediv(&py::abs(&py::getitem(&v_d.clone(), &v_i.clone())?)?, &py::minmax(&[py::abs(&py::getitem(&v_p.clone(), &v_i.clone())?)?, Py::Float(1e-300)], true)?)?); } Py::list(out) }], true)?; Py::Bool({ let t2718 = Py::Float(1e-15); if py::lt(&t2717, &t2718)? { true } else { false } }) };
                    let t2720 = (Py::tuple(vec![v_q.clone(), v_cq.clone()])).items()?;
                    if t2720.len() != 2 { return Err(PyErr::Value("unpack".into())); }
                    v_p = t2720[0].clone();
                    v_c = t2720[1].clone();
                    v_lam = py::minmax(&[py::truediv(&v_lam.clone(), &Py::Int(10i64))?, Py::Float(1e-15)], true)?;
                    if (v_small.clone()).truthy() {
                        return Ok(f__gn(v_model.clone(), v_p.clone(), v_ys.clone(), v_ok.clone())?);
                    }
                    break;
                }
            }
            v_lam = py::mul(&v_lam, &Py::Int(10i64))?;
            if ({ let t2739 = v_lam.clone(); Py::Bool({ let t2740 = Py::Float(1e+20); if py::gt(&t2739, &t2740)? { true } else { false } }) }).truthy() {
                return Ok(f__gn(v_model.clone(), v_p.clone(), v_ys.clone(), v_ok.clone())?);
            }
        }
    }
    return Ok(f__gn(v_model.clone(), v_p.clone(), v_ys.clone(), v_ok.clone())?);
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn f__locus_dist(mut v_x: Py, mut v_y: Py) -> R {
    let mut v_d: Py = Py::None;
    let mut v_k: Py = Py::None;
    let mut v_u: Py = Py::None;
    let mut v_v: Py = Py::None;
    v_d = Py::dict(vec![]);
    for t911 in py::iter(&(v_x.clone()))? {
        v_u = t911;
        for t912 in py::iter(&(v_y.clone()))? {
            v_v = t912;
            v_k = Py::tuple((py::sorted_key(&Py::tuple(vec![v_u.clone(), v_v.clone()]), &{ let v_d = v_d.clone(); let v_u = v_u.clone(); let v_v = v_v.clone(); let v_x = v_x.clone(); let v_y = v_y.clone(); py::func(move |a: &[Py]| -> R {
            #[allow(unused_mut, unused_variables, unused_assignments)]
{
            let mut v_c: Py = a.get(0).cloned().unwrap_or(Py::None);
            return Ok(Py::tuple(vec![py::method(&v_c.clone(), "islower", &[])?, v_c.clone()]));
            #[allow(unreachable_code)]
            Ok(Py::None)
        }
        }) }, (Py::Bool(false)).truthy())?).items()?);
            py::setitem(&v_d.clone(), &v_k.clone(), py::add(&py::method(&v_d.clone(), "get", &[v_k.clone(), Py::Int(0i64)])?, &py::fraction(&[Py::Int(1i64), Py::Int(4i64)])?)?)?;
        }
    }
    return Ok(v_d.clone());
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn f__log_abs_expm1(mut v_z: Py) -> R {
    if ({ let t4415 = v_z.clone(); Py::Bool({ let t4416 = Py::Int(0i64); if py::gt(&t4415, &t4416)? { true } else { false } }) }).truthy() {
        return Ok((if ({ let t4417 = v_z.clone(); Py::Bool({ let t4418 = Py::Int(30i64); if py::gt(&t4417, &t4418)? { true } else { false } }) }).truthy() { py::add(&v_z.clone(), &py::math1("log1p", &py::neg(&py::math1("exp", &py::neg(&v_z.clone())?)?)?)?)? } else { py::math1("log", &py::math1("expm1", &v_z.clone())?)? }));
    }
    return Ok(py::math1("log", &py::neg(&py::math1("expm1", &v_z.clone())?)?)?);
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn f__mm(mut v_match: Py, mut v_mismatch: Py) -> R {
    return Ok({ let v_match = v_match.clone(); let v_mismatch = v_mismatch.clone(); py::func(move |a: &[Py]| -> R {
            #[allow(unused_mut, unused_variables, unused_assignments)]
{
            let mut v_x: Py = a.get(0).cloned().unwrap_or(Py::None);
            let mut v_y: Py = a.get(1).cloned().unwrap_or(Py::None);
            return Ok((if ({ let t6292 = v_x.clone(); Py::Bool({ let t6293 = v_y.clone(); if py::eq(&t6292, &t6293) { true } else { false } }) }).truthy() { v_match.clone() } else { v_mismatch.clone() }));
            #[allow(unreachable_code)]
            Ok(Py::None)
        }
        }) });
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn f__mm_check(mut v_S0: Py, mut v_Vmax: Py, mut v_Km: Py) -> R {
    let _ = f__pos(Py::str("S0"), Py::tuple(vec![v_S0.clone()]))?;
    let _ = f__pos(Py::str("Vmax"), Py::tuple(vec![v_Vmax.clone()]))?;
    let _ = f__pos(Py::str("Km"), Py::tuple(vec![v_Km.clone()]))?;
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn f__mm_varpro(mut v_S: Py, mut v_v: Py, mut v_K: Py) -> R {
    let mut v_svx: Py = Py::None;
    let mut v_svxd: Py = Py::None;
    let mut v_sxx: Py = Py::None;
    let mut v_sxxd: Py = Py::None;
    let mut v_x: Py = Py::None;
    let mut v_xd: Py = Py::None;
    v_x = { let mut out = vec![]; for t4125 in py::iter(&(v_S.clone()))? { let mut v_s: Py; 
v_s = t4125;
out.push(py::truediv(&v_s.clone(), &py::add(&v_K.clone(), &v_s.clone())?)?); } Py::list(out) };
    v_xd = { let mut out = vec![]; for t4126 in py::iter(&(v_S.clone()))? { let mut v_s: Py; 
v_s = t4126;
out.push(py::truediv(&py::neg(&v_s.clone())?, &py::pow(&py::add(&v_K.clone(), &v_s.clone())?, &Py::Int(2i64))?)?); } Py::list(out) };
    v_svx = py::fsum(&{ let mut out = vec![]; for t4127 in py::iter(&(py::zip(&[v_v.clone(), v_x.clone()])?))? { let mut v_a: Py; let mut v_b: Py; 
let t4128 = (t4127).items()?;
if t4128.len() != 2 { return Err(PyErr::Value("unpack".into())); }
v_a = t4128[0].clone();
v_b = t4128[1].clone();
out.push(py::mul(&v_a.clone(), &v_b.clone())?); } Py::list(out) })?;
    v_sxx = py::fsum(&{ let mut out = vec![]; for t4129 in py::iter(&(v_x.clone()))? { let mut v_a: Py; 
v_a = t4129;
out.push(py::mul(&v_a.clone(), &v_a.clone())?); } Py::list(out) })?;
    v_svxd = py::fsum(&{ let mut out = vec![]; for t4130 in py::iter(&(py::zip(&[v_v.clone(), v_xd.clone()])?))? { let mut v_a: Py; let mut v_b: Py; 
let t4131 = (t4130).items()?;
if t4131.len() != 2 { return Err(PyErr::Value("unpack".into())); }
v_a = t4131[0].clone();
v_b = t4131[1].clone();
out.push(py::mul(&v_a.clone(), &v_b.clone())?); } Py::list(out) })?;
    v_sxxd = py::fsum(&{ let mut out = vec![]; for t4132 in py::iter(&(py::zip(&[v_x.clone(), v_xd.clone()])?))? { let mut v_a: Py; let mut v_b: Py; 
let t4133 = (t4132).items()?;
if t4133.len() != 2 { return Err(PyErr::Value("unpack".into())); }
v_a = t4133[0].clone();
v_b = t4133[1].clone();
out.push(py::mul(&v_a.clone(), &v_b.clone())?); } Py::list(out) })?;
    return Ok(Py::tuple(vec![v_svx.clone(), v_sxx.clone(), py::sub(&py::mul(&v_svxd.clone(), &v_sxx.clone())?, &py::mul(&v_svx.clone(), &v_sxxd.clone())?)?]));
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn f__mutsel_map(mut v_q: Py, mut v_mu: Py, mut v_s: Py, mut v_h: Py) -> R {
    let mut v_p: Py = Py::None;
    let mut v_qs: Py = Py::None;
    let mut v_wbar: Py = Py::None;
    v_p = py::sub(&Py::Int(1i64), &v_q.clone())?;
    v_wbar = py::add(&py::add(&py::mul(&v_p.clone(), &v_p.clone())?, &py::mul(&py::mul(&py::mul(&Py::Int(2i64), &v_p.clone())?, &v_q.clone())?, &py::sub(&Py::Int(1i64), &py::mul(&v_h.clone(), &v_s.clone())?)?)?)?, &py::mul(&py::mul(&v_q.clone(), &v_q.clone())?, &py::sub(&Py::Int(1i64), &v_s.clone())?)?)?;
    v_qs = py::truediv(&py::add(&py::mul(&py::mul(&v_p.clone(), &v_q.clone())?, &py::sub(&Py::Int(1i64), &py::mul(&v_h.clone(), &v_s.clone())?)?)?, &py::mul(&py::mul(&v_q.clone(), &v_q.clone())?, &py::sub(&Py::Int(1i64), &v_s.clone())?)?)?, &v_wbar.clone())?;
    return Ok(py::add(&v_qs.clone(), &py::mul(&v_mu.clone(), &py::sub(&Py::Int(1i64), &v_qs.clone())?)?)?);
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn f__nc_model(mut v_S: Py, mut v_I: Py) -> R {
    let mut v_m: Py = Py::None;
    v_m = { let v_I = v_I.clone(); let v_S = v_S.clone(); let v_m = v_m.clone(); py::func(move |a: &[Py]| -> R {
            #[allow(unused_mut, unused_variables, unused_assignments)]
{
            let mut v_p: Py = a.get(0).cloned().unwrap_or(Py::None);
            let mut v_J: Py = Py::None;
            let mut v_Ki: Py = Py::None;
            let mut v_Km: Py = Py::None;
            let mut v_Vm: Py = Py::None;
            let mut v_a: Py = Py::None;
            let mut v_f: Py = Py::None;
            let mut v_i: Py = Py::None;
            let mut v_s: Py = Py::None;
            let t695 = (v_p.clone()).items()?;
            if t695.len() != 3 { return Err(PyErr::Value("unpack".into())); }
            v_Vm = t695[0].clone();
            v_Km = t695[1].clone();
            v_Ki = t695[2].clone();
            let t696 = (Py::tuple(vec![Py::list(vec![]), Py::list(vec![])])).items()?;
            if t696.len() != 2 { return Err(PyErr::Value("unpack".into())); }
            v_f = t696[0].clone();
            v_J = t696[1].clone();
            for t697 in py::iter(&(py::zip(&[v_S.clone(), v_I.clone()])?))? {
                let t698 = (t697).items()?;
                if t698.len() != 2 { return Err(PyErr::Value("unpack".into())); }
                v_s = t698[0].clone();
                v_i = t698[1].clone();
                v_a = py::add(&Py::Int(1i64), &py::truediv(&v_i.clone(), &v_Ki.clone())?)?;
                let _ = py::method(&v_f.clone(), "append", &[py::truediv(&py::mul(&v_Vm.clone(), &v_s.clone())?, &py::mul(&py::add(&v_Km.clone(), &v_s.clone())?, &v_a.clone())?)?])?;
                let _ = py::method(&v_J.clone(), "append", &[Py::list(vec![py::truediv(&v_s.clone(), &py::mul(&py::add(&v_Km.clone(), &v_s.clone())?, &v_a.clone())?)?, py::truediv(&py::mul(&py::neg(&v_Vm.clone())?, &v_s.clone())?, &py::mul(&py::pow(&py::add(&v_Km.clone(), &v_s.clone())?, &Py::Int(2i64))?, &v_a.clone())?)?, py::truediv(&py::truediv(&py::mul(&py::mul(&v_Vm.clone(), &v_s.clone())?, &v_i.clone())?, &py::mul(&v_Ki.clone(), &v_Ki.clone())?)?, &py::mul(&py::mul(&py::add(&v_Km.clone(), &v_s.clone())?, &v_a.clone())?, &v_a.clone())?)?])])?;
            }
            return Ok(Py::tuple(vec![v_f.clone(), v_J.clone()]));
            #[allow(unreachable_code)]
            Ok(Py::None)
        }
        }) };
    return Ok(v_m.clone());
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn f__need(mut v_ok: Py, mut v_msg: Py) -> R {
    if (py::not(&v_ok.clone())).truthy() {
        return Err(py::value_error(py::to_str(&v_msg.clone())));
    }
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn f__ols(mut v_x: Py, mut v_y: Py) -> R {
    let mut v_b: Py = Py::None;
    let mut v_mx: Py = Py::None;
    let mut v_my: Py = Py::None;
    let mut v_n: Py = Py::None;
    let mut v_r2: Py = Py::None;
    let mut v_sxx: Py = Py::None;
    let mut v_sxy: Py = Py::None;
    let mut v_syy: Py = Py::None;
    v_n = py::len(&v_x.clone())?;
    let t6494 = (Py::tuple(vec![py::truediv(&py::fsum(&v_x.clone())?, &v_n.clone())?, py::truediv(&py::fsum(&v_y.clone())?, &v_n.clone())?])).items()?;
    if t6494.len() != 2 { return Err(PyErr::Value("unpack".into())); }
    v_mx = t6494[0].clone();
    v_my = t6494[1].clone();
    v_sxx = py::fsum(&{ let mut out = vec![]; for t6495 in py::iter(&(v_x.clone()))? { let mut v_a: Py; 
v_a = t6495;
out.push(py::pow(&py::sub(&v_a.clone(), &v_mx.clone())?, &Py::Int(2i64))?); } Py::list(out) })?;
    v_syy = py::fsum(&{ let mut out = vec![]; for t6496 in py::iter(&(v_y.clone()))? { let mut v_b: Py; 
v_b = t6496;
out.push(py::pow(&py::sub(&v_b.clone(), &v_my.clone())?, &Py::Int(2i64))?); } Py::list(out) })?;
    v_sxy = py::fsum(&{ let mut out = vec![]; for t6497 in py::iter(&(py::zip(&[v_x.clone(), v_y.clone()])?))? { let mut v_a: Py; let mut v_b: Py; 
let t6498 = (t6497).items()?;
if t6498.len() != 2 { return Err(PyErr::Value("unpack".into())); }
v_a = t6498[0].clone();
v_b = t6498[1].clone();
out.push(py::mul(&py::sub(&v_a.clone(), &v_mx.clone())?, &py::sub(&v_b.clone(), &v_my.clone())?)?); } Py::list(out) })?;
    let _ = f__need({ let t6499 = v_sxx.clone(); Py::Bool({ let t6500 = Py::Int(0i64); if py::gt(&t6499, &t6500)? { true } else { false } }) }, Py::str("the x values must not all be equal"))?;
    v_b = py::truediv(&v_sxy.clone(), &v_sxx.clone())?;
    v_r2 = (if ({ let t6501 = v_syy.clone(); Py::Bool({ let t6502 = Py::Int(0i64); if py::gt(&t6501, &t6502)? { true } else { false } }) }).truthy() { py::truediv(&py::mul(&v_sxy.clone(), &v_sxy.clone())?, &py::mul(&v_sxx.clone(), &v_syy.clone())?)? } else { Py::Float(1.0) });
    return Ok(Py::tuple(vec![v_b.clone(), py::sub(&v_my.clone(), &py::mul(&v_b.clone(), &v_mx.clone())?)?, v_r2.clone()]));
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn f__orf_lengths(mut v_s: Py) -> R {
    let mut v_c: Py = Py::None;
    let mut v_f: Py = Py::None;
    let mut v_i: Py = Py::None;
    let mut v_res: Py = Py::None;
    let mut v_start: Py = Py::None;
    v_res = Py::list(vec![]);
    for t3989 in py::range_iter(&[(Py::Int(3i64)).i()?]) {
        v_f = t3989;
        v_start = Py::None;
        for t3990 in py::range_iter(&[(v_f.clone()).i()?, (py::sub(&py::len(&v_s.clone())?, &Py::Int(2i64))?).i()?, (Py::Int(3i64)).i()?]) {
            v_i = t3990;
            v_c = py::slice(&v_s.clone(), Some((v_i.clone()).i()?), Some((py::add(&v_i.clone(), &Py::Int(3i64))?).i()?), None)?;
            if ({ let t3991 = v_start.clone(); Py::Bool({ let t3992 = Py::None; if py::is(&t3991, &t3992) { true } else { false } }) }).truthy() {
                if ({ let t3993 = v_c.clone(); Py::Bool({ let t3994 = Py::str("ATG"); if py::eq(&t3993, &t3994) { true } else { false } }) }).truthy() {
                    v_start = v_i.clone();
                }
            } else {
                if ({ let t3995 = v_c.clone(); Py::Bool({ let t3996 = c_STOPS()?; if py::contains(&t3996, &t3995)? { true } else { false } }) }).truthy() {
                    let _ = py::method(&v_res.clone(), "append", &[py::sub(&py::add(&v_i.clone(), &Py::Int(3i64))?, &v_start.clone())?])?;
                    v_start = Py::None;
                }
            }
        }
    }
    return Ok(v_res.clone());
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn f__osd(mut v_c: Py) -> R {
    return Ok(Py::tuple(vec![py::len(&v_c.clone())?, py::sum(&{ let mut out = vec![]; for t393 in py::iter(&(v_c.clone()))? { let mut v_x: Py; 
v_x = t393;
if ({ let t394 = v_x.clone(); Py::Bool({ let t395 = Py::Int(1i64); if py::eq(&t394, &t395) { true } else { false } }) }).truthy() { out.push(Py::Int(1i64)); }} Py::list(out) }, None)?, py::sum(&{ let mut out = vec![]; for t396 in py::iter(&(v_c.clone()))? { let mut v_x: Py; 
v_x = t396;
if ({ let t397 = v_x.clone(); Py::Bool({ let t398 = Py::Int(2i64); if py::eq(&t397, &t398) { true } else { false } }) }).truthy() { out.push(Py::Int(1i64)); }} Py::list(out) }, None)?]));
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn f__out(mut v_kw: Py) -> R {
    return Ok({ let mut out = vec![]; for t7120 in py::iter(&(py::method(&v_kw.clone(), "items", &[])?))? { let mut v_k: Py; let mut v_v: Py; 
let t7121 = (t7120).items()?;
if t7121.len() != 2 { return Err(PyErr::Value("unpack".into())); }
v_k = t7121[0].clone();
v_v = t7121[1].clone();
out.push((v_k.clone(), f__r(v_v.clone())?)); } Py::dict(out) });
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn f__pairwise_k(mut v_rows: Py) -> R {
    let mut v_n: Py = Py::None;
    let mut v_tot: Py = Py::None;
    v_n = py::len(&v_rows.clone())?;
    v_tot = py::sum(&{ let mut out = vec![]; for t6575 in py::range_iter(&[(v_n.clone()).i()?]) { let mut v_i: Py; 
v_i = t6575;
for t6576 in py::range_iter(&[(py::add(&v_i.clone(), &Py::Int(1i64))?).i()?, (v_n.clone()).i()?]) { let mut v_j: Py; 
v_j = t6576;
out.push(py::sum(&{ let mut out = vec![]; for t6577 in py::iter(&(py::zip(&[py::getitem(&v_rows.clone(), &v_i.clone())?, py::getitem(&v_rows.clone(), &v_j.clone())?])?))? { let mut v_x: Py; let mut v_y: Py; 
let t6578 = (t6577).items()?;
if t6578.len() != 2 { return Err(PyErr::Value("unpack".into())); }
v_x = t6578[0].clone();
v_y = t6578[1].clone();
if ({ let t6579 = v_x.clone(); Py::Bool({ let t6580 = v_y.clone(); if py::ne(&t6579, &t6580) { true } else { false } }) }).truthy() { out.push(Py::Int(1i64)); }} Py::list(out) }, None)?); }} Py::list(out) }, None)?;
    return Ok(py::fraction(&[v_tot.clone(), py::floordiv(&py::mul(&v_n.clone(), &py::sub(&v_n.clone(), &Py::Int(1i64))?)?, &Py::Int(2i64))?])?);
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn f__parse_newick(mut v_text: Py) -> R {
    let mut v_leaves: Py = Py::None;
    let mut v_node: Py = Py::None;
    let mut v_num: Py = Py::None;
    let mut v_pos: Py = Py::None;
    let mut v_root: Py = Py::None;
    let mut v_s: Py = Py::None;
    let mut v_walk: Py = Py::None;
    let _ = f__need(Py::Bool(matches!(v_text.clone(), Py::Str(_))), Py::str("the tree must be a Newick string"))?;
    v_s = py::method(&ax_engine::pyre::module_call("sub", &Py::str("\\s+"), &[Py::str(""), v_text.clone(), Py::Int(0)], 0)?, "rstrip", &[Py::str(";")])?;
    v_pos = Py::list(vec![Py::Int(0i64)]);
    v_num = { let v_num = v_num.clone(); let v_pos = v_pos.clone(); let v_s = v_s.clone(); let v_text = v_text.clone(); py::func(move |a: &[Py]| -> R {
            #[allow(unused_mut, unused_variables, unused_assignments)]
{
            let mut v_m: Py = Py::None;
            v_m = ax_engine::pyre::module_call("match", &Py::str("[-+]?(?:\\d+(?:\\.\\d*)?|\\.\\d+)(?:[eE][-+]?\\d+)?"), &[py::slice(&v_s.clone(), Some((py::getitem(&v_pos.clone(), &Py::Int(0i64))?).i()?), None, None)?], 0)?;
            let _ = f__need({ let t5974 = v_m.clone(); Py::Bool({ let t5975 = Py::None; if !py::is(&t5974, &t5975) { true } else { false } }) }, Py::str("a branch length must be a number"))?;
            let t5976 = v_pos.clone(); let t5977 = Py::Int(0i64);
            py::setitem(&t5976, &t5977, py::iadd(&py::getitem(&t5976, &t5977)?, &py::method(&v_m.clone(), "end", &[])?)?)?;
            return Ok(py::float(&py::method(&v_m.clone(), "group", &[Py::Int(0i64)])?)?);
            #[allow(unreachable_code)]
            Ok(Py::None)
        }
        }) };
    let rec5978 = py::rec_cell();
    v_node = { let v_num = v_num.clone(); let v_pos = v_pos.clone(); let v_s = v_s.clone(); let v_text = v_text.clone(); let rec5978 = rec5978.clone(); py::func(move |a: &[Py]| -> R {
            #[allow(unused_mut, unused_variables, unused_assignments)]
{
            let v_node: Py = py::rec_get(&rec5978)?;
            let mut v_kids: Py = Py::None;
            let mut v_length: Py = Py::None;
            let mut v_m: Py = Py::None;
            let mut v_name: Py = Py::None;
            v_kids = Py::list(vec![]);
            if ({ let t5983 = { let t5979 = py::getitem(&v_pos.clone(), &Py::Int(0i64))?; Py::Bool({ let t5980 = py::len(&v_s.clone())?; if py::lt(&t5979, &t5980)? { true } else { false } }) }; if !t5983.truthy() { t5983 } else { { let t5981 = py::getitem(&v_s.clone(), &py::getitem(&v_pos.clone(), &Py::Int(0i64))?)?; Py::Bool({ let t5982 = Py::str("("); if py::eq(&t5981, &t5982) { true } else { false } }) } } }).truthy() {
                let t5984 = v_pos.clone(); let t5985 = Py::Int(0i64);
                py::setitem(&t5984, &t5985, py::iadd(&py::getitem(&t5984, &t5985)?, &Py::Int(1i64))?)?;
                let _ = py::method(&v_kids.clone(), "append", &[py::call(&v_node.clone(), &[])?])?;
                while ({ let t5992 = { let t5988 = py::getitem(&v_pos.clone(), &Py::Int(0i64))?; Py::Bool({ let t5989 = py::len(&v_s.clone())?; if py::lt(&t5988, &t5989)? { true } else { false } }) }; if !t5992.truthy() { t5992 } else { { let t5990 = py::getitem(&v_s.clone(), &py::getitem(&v_pos.clone(), &Py::Int(0i64))?)?; Py::Bool({ let t5991 = Py::str(","); if py::eq(&t5990, &t5991) { true } else { false } }) } } }).truthy() {
                    let t5986 = v_pos.clone(); let t5987 = Py::Int(0i64);
                    py::setitem(&t5986, &t5987, py::iadd(&py::getitem(&t5986, &t5987)?, &Py::Int(1i64))?)?;
                    let _ = py::method(&v_kids.clone(), "append", &[py::call(&v_node.clone(), &[])?])?;
                }
                let _ = f__need({ let t5997 = { let t5993 = py::getitem(&v_pos.clone(), &Py::Int(0i64))?; Py::Bool({ let t5994 = py::len(&v_s.clone())?; if py::lt(&t5993, &t5994)? { true } else { false } }) }; if !t5997.truthy() { t5997 } else { { let t5995 = py::getitem(&v_s.clone(), &py::getitem(&v_pos.clone(), &Py::Int(0i64))?)?; Py::Bool({ let t5996 = Py::str(")"); if py::eq(&t5995, &t5996) { true } else { false } }) } } }, Py::str("unbalanced parentheses in the tree"))?;
                let t5998 = v_pos.clone(); let t5999 = Py::Int(0i64);
                py::setitem(&t5998, &t5999, py::iadd(&py::getitem(&t5998, &t5999)?, &Py::Int(1i64))?)?;
            }
            v_m = ax_engine::pyre::module_call("match", &Py::str("[A-Za-z0-9_.]*"), &[py::slice(&v_s.clone(), Some((py::getitem(&v_pos.clone(), &Py::Int(0i64))?).i()?), None, None)?], 0)?;
            v_name = py::method(&v_m.clone(), "group", &[Py::Int(0i64)])?;
            let t6000 = v_pos.clone(); let t6001 = Py::Int(0i64);
            py::setitem(&t6000, &t6001, py::iadd(&py::getitem(&t6000, &t6001)?, &py::method(&v_m.clone(), "end", &[])?)?)?;
            v_length = Py::None;
            if ({ let t6006 = { let t6002 = py::getitem(&v_pos.clone(), &Py::Int(0i64))?; Py::Bool({ let t6003 = py::len(&v_s.clone())?; if py::lt(&t6002, &t6003)? { true } else { false } }) }; if !t6006.truthy() { t6006 } else { { let t6004 = py::getitem(&v_s.clone(), &py::getitem(&v_pos.clone(), &Py::Int(0i64))?)?; Py::Bool({ let t6005 = Py::str(":"); if py::eq(&t6004, &t6005) { true } else { false } }) } } }).truthy() {
                let t6007 = v_pos.clone(); let t6008 = Py::Int(0i64);
                py::setitem(&t6007, &t6008, py::iadd(&py::getitem(&t6007, &t6008)?, &Py::Int(1i64))?)?;
                v_length = py::call(&v_num.clone(), &[])?;
            }
            if (py::not(&v_kids.clone())).truthy() {
                let _ = f__need({ let t6009 = ax_engine::pyre::module_call("fullmatch", &Py::str("[A-Za-z][A-Za-z0-9_]*"), &[{ let t6011 = v_name.clone(); if t6011.truthy() { t6011 } else { Py::str("") } }], 0)?; Py::Bool({ let t6010 = Py::None; if !py::is(&t6009, &t6010) { true } else { false } }) }, Py::str("every leaf needs a name"))?;
            }
            return Ok(Py::tuple(vec![v_name.clone(), v_length.clone(), v_kids.clone()]));
            #[allow(unreachable_code)]
            Ok(Py::None)
        }
        }) };
    py::rec_set(&rec5978, &v_node);
    v_root = py::call(&v_node.clone(), &[])?;
    let _ = f__need({ let t6012 = py::getitem(&v_pos.clone(), &Py::Int(0i64))?; Py::Bool({ let t6013 = py::len(&v_s.clone())?; if py::eq(&t6012, &t6013) { true } else { false } }) }, Py::str("unreadable text after the tree"))?;
    v_leaves = Py::list(vec![]);
    let rec6014 = py::rec_cell();
    v_walk = { let v_leaves = v_leaves.clone(); let v_node = v_node.clone(); let v_num = v_num.clone(); let v_pos = v_pos.clone(); let v_root = v_root.clone(); let v_s = v_s.clone(); let v_text = v_text.clone(); let rec6014 = rec6014.clone(); py::func(move |a: &[Py]| -> R {
            #[allow(unused_mut, unused_variables, unused_assignments)]
{
            let mut v_t: Py = a.get(0).cloned().unwrap_or(Py::None);
            let v_walk: Py = py::rec_get(&rec6014)?;
            let mut v_k: Py = Py::None;
            if (py::not(&py::getitem(&v_t.clone(), &Py::Int(2i64))?)).truthy() {
                let _ = py::method(&v_leaves.clone(), "append", &[py::getitem(&v_t.clone(), &Py::Int(0i64))?])?;
            }
            for t6015 in py::iter(&(py::getitem(&v_t.clone(), &Py::Int(2i64))?))? {
                v_k = t6015;
                let _ = py::call(&v_walk.clone(), &[v_k.clone()])?;
            }
            #[allow(unreachable_code)]
            Ok(Py::None)
        }
        }) };
    py::rec_set(&rec6014, &v_walk);
    let _ = py::call(&v_walk.clone(), &[v_root.clone()])?;
    let _ = f__need({ let t6016 = py::len(&py::set(&v_leaves.clone())?)?; Py::Bool({ let t6017 = py::len(&v_leaves.clone())?; if py::eq(&t6016, &t6017) { true } else { false } }) }, Py::str("leaf names must be unique"))?;
    return Ok(Py::tuple(vec![v_root.clone(), v_leaves.clone()]));
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn f__pedigree(mut v_pedigree: Py) -> R {
    let mut v_c: Py = Py::None;
    let mut v_done: Py = Py::None;
    let mut v_e: Py = Py::None;
    let mut v_p: Py = Py::None;
    let mut v_parents: Py = Py::None;
    let mut v_ps: Py = Py::None;
    let mut v_seen: Py = Py::None;
    let mut v_visit: Py = Py::None;
    v_parents = Py::dict(vec![]);
    for t4939 in py::iter(&(v_pedigree.clone()))? {
        v_e = t4939;
        let _ = f__need({ let t4940 = py::len(&v_e.clone())?; Py::Bool({ let t4941 = Py::Int(2i64); if py::eq(&t4940, &t4941) { true } else { false } }) }, Py::str("pedigree edges are parent-child pairs without weights"))?;
        let t4942 = (v_e.clone()).items()?;
        if t4942.len() != 2 { return Err(PyErr::Value("unpack".into())); }
        v_p = t4942[0].clone();
        v_c = t4942[1].clone();
        let _ = f__need({ let t4943 = v_p.clone(); Py::Bool({ let t4944 = v_c.clone(); if py::ne(&t4943, &t4944) { true } else { false } }) }, Py::str("an individual cannot be its own parent"))?;
        let _ = py::method(&v_parents.clone(), "setdefault", &[v_c.clone(), Py::list(vec![])])?;
        let _ = f__need({ let t4945 = v_p.clone(); Py::Bool({ let t4946 = py::getitem(&v_parents.clone(), &v_c.clone())?; if !py::contains(&t4946, &t4945)? { true } else { false } }) }, Py::str(&[py::format(&v_p.clone(), "")?, " is listed twice as a parent of ".to_string(), py::format(&v_c.clone(), "")?].concat()))?;
        let _ = py::method(&py::getitem(&v_parents.clone(), &v_c.clone())?, "append", &[v_p.clone()])?;
    }
    for t4947 in py::iter(&(py::method(&v_parents.clone(), "items", &[])?))? {
        let t4950 = (t4947).items()?;
        if t4950.len() != 2 { return Err(PyErr::Value("unpack".into())); }
        v_c = t4950[0].clone();
        v_ps = t4950[1].clone();
        let _ = f__need({ let t4948 = py::len(&v_ps.clone())?; Py::Bool({ let t4949 = Py::Int(2i64); if py::eq(&t4948, &t4949) { true } else { false } }) }, Py::str(&[py::format(&v_c.clone(), "")?, " must have exactly two parents listed (founders have none)".to_string()].concat()))?;
    }
    let t4951 = (Py::tuple(vec![py::set(&Py::list(vec![]))?, py::set(&Py::list(vec![]))?])).items()?;
    if t4951.len() != 2 { return Err(PyErr::Value("unpack".into())); }
    v_seen = t4951[0].clone();
    v_done = t4951[1].clone();
    let rec4952 = py::rec_cell();
    v_visit = { let v_c = v_c.clone(); let v_done = v_done.clone(); let v_e = v_e.clone(); let v_p = v_p.clone(); let v_parents = v_parents.clone(); let v_pedigree = v_pedigree.clone(); let v_ps = v_ps.clone(); let v_seen = v_seen.clone(); let rec4952 = rec4952.clone(); py::func(move |a: &[Py]| -> R {
            #[allow(unused_mut, unused_variables, unused_assignments)]
{
            let mut v_x: Py = a.get(0).cloned().unwrap_or(Py::None);
            let v_visit: Py = py::rec_get(&rec4952)?;
            let mut v_p: Py = Py::None;
            if ({ let t4953 = v_x.clone(); Py::Bool({ let t4954 = v_done.clone(); if py::contains(&t4954, &t4953)? { true } else { false } }) }).truthy() {
                return Ok(Py::None);
            }
            let _ = f__need({ let t4955 = v_x.clone(); Py::Bool({ let t4956 = v_seen.clone(); if !py::contains(&t4956, &t4955)? { true } else { false } }) }, Py::str("the pedigree has a cycle"))?;
            let _ = py::method(&v_seen.clone(), "add", &[v_x.clone()])?;
            for t4957 in py::iter(&(py::method(&v_parents.clone(), "get", &[v_x.clone(), Py::list(vec![])])?))? {
                v_p = t4957;
                let _ = py::call(&v_visit.clone(), &[v_p.clone()])?;
            }
            let _ = py::method(&v_done.clone(), "add", &[v_x.clone()])?;
            #[allow(unreachable_code)]
            Ok(Py::None)
        }
        }) };
    py::rec_set(&rec4952, &v_visit);
    for t4958 in py::iter(&(Py::list((v_parents.clone()).items()?)))? {
        v_c = t4958;
        let _ = py::call(&v_visit.clone(), &[v_c.clone()])?;
    }
    return Ok(v_parents.clone());
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn f__pos(mut v_name: Py, mut v_vals: Py) -> R {
    let mut v_v: Py = Py::None;
    for t6664 in py::iter(&(v_vals.clone()))? {
        v_v = t6664;
        let _ = f__need({ let t6669 = Py::Bool(matches!(v_v.clone(), Py::Int(_) | Py::Big(_) | Py::Bool(_) | Py::Float(_))); if !t6669.truthy() { t6669 } else { { let t6668 = py::not(&Py::Bool(matches!(v_v.clone(), Py::Bool(_)))); if !t6668.truthy() { t6668 } else { { let t6667 = py::math1("isfinite", &v_v.clone())?; if !t6667.truthy() { t6667 } else { { let t6665 = v_v.clone(); Py::Bool({ let t6666 = Py::Int(0i64); if py::gt(&t6665, &t6666)? { true } else { false } }) } } } } } } }, Py::str(&[py::format(&v_name.clone(), "")?, " must be a positive number".to_string()].concat()))?;
    }
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn f__r(mut v_x: Py) -> R {
    if ({ let t7125 = Py::Bool(matches!(v_x.clone(), Py::Bool(_))); if t7125.truthy() { t7125 } else { { let t7124 = Py::Bool(matches!(v_x.clone(), Py::Str(_))); if t7124.truthy() { t7124 } else { { let t7122 = v_x.clone(); Py::Bool({ let t7123 = Py::None; if py::is(&t7122, &t7123) { true } else { false } }) } } } } }).truthy() {
        return Ok(v_x.clone());
    }
    if (Py::Bool(matches!(v_x.clone(), Py::Int(_) | Py::Big(_) | Py::Bool(_)))).truthy() {
        return Ok(v_x.clone());
    }
    if (Py::Bool(matches!(v_x.clone(), Py::Dict(_)))).truthy() {
        return Ok({ let mut out = vec![]; for t7126 in py::iter(&(py::method(&v_x.clone(), "items", &[])?))? { let mut v_k: Py; let mut v_v: Py; 
let t7127 = (t7126).items()?;
if t7127.len() != 2 { return Err(PyErr::Value("unpack".into())); }
v_k = t7127[0].clone();
v_v = t7127[1].clone();
out.push((v_k.clone(), f__r(v_v.clone())?)); } Py::dict(out) });
    }
    if (Py::Bool(matches!(v_x.clone(), Py::List(_) | Py::Tuple(_)))).truthy() {
        return Ok({ let mut out = vec![]; for t7128 in py::iter(&(v_x.clone()))? { let mut v_v: Py; 
v_v = t7128;
out.push(f__r(v_v.clone())?); } Py::list(out) });
    }
    if ({ let t7131 = { let t7129 = v_x.clone(); Py::Bool({ let t7130 = Py::Int(0i64); if py::eq(&t7129, &t7130) { true } else { false } }) }; if t7131.truthy() { t7131 } else { py::not(&py::math1("isfinite", &v_x.clone())?) } }).truthy() {
        return Ok((if ({ let t7132 = v_x.clone(); Py::Bool({ let t7133 = Py::Int(0i64); if py::eq(&t7132, &t7133) { true } else { false } }) }).truthy() { Py::Float(0.0) } else { v_x.clone() }));
    }
    return Ok(py::add(&py::float(&Py::str(&[py::format(&v_x.clone(), ".12g")?].concat()))?, &Py::Float(0.0))?);
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn f__ratio_text(mut v_probs: Py) -> R {
    let mut v_den: Py = Py::None;
    let mut v_g: Py = Py::None;
    let mut v_ints: Py = Py::None;
    let mut v_p: Py = Py::None;
    let mut v_v: Py = Py::None;
    v_probs = { let mut out = vec![]; for t931 in py::iter(&(v_probs.clone()))? { let mut v_p: Py; 
v_p = t931;
if ({ let t932 = v_p.clone(); Py::Bool({ let t933 = Py::Int(0i64); if py::gt(&t932, &t933)? { true } else { false } }) }).truthy() { out.push(v_p.clone()); }} Py::list(out) };
    v_den = Py::Int(1i64);
    for t934 in py::iter(&(v_probs.clone()))? {
        v_p = t934;
        v_den = py::floordiv(&py::mul(&v_den.clone(), &py::denominator(&v_p.clone())?)?, &py::gcd(&v_den.clone(), &py::denominator(&v_p.clone())?)?)?;
    }
    v_ints = { let mut out = vec![]; for t935 in py::iter(&(v_probs.clone()))? { let mut v_p: Py; 
v_p = t935;
out.push(py::int(&py::mul(&v_p.clone(), &v_den.clone())?)?); } Py::list(out) };
    v_g = Py::Int(0i64);
    for t936 in py::iter(&(v_ints.clone()))? {
        v_v = t936;
        v_g = py::gcd(&v_g.clone(), &v_v.clone())?;
    }
    return Ok(py::method(&Py::str(":"), "join", &[{ let mut out = vec![]; for t937 in py::iter(&(py::sorted(&v_ints.clone(), (Py::Bool(true)).truthy())?))? { let mut v_v: Py; 
v_v = t937;
out.push(Py::str(&py::to_str(&py::floordiv(&v_v.clone(), &v_g.clone())?))); } Py::list(out) }])?);
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn f__rel(mut v_a: Py, mut v_b: Py) -> R {
    let mut v_d: Py = Py::None;
    if ({ let t7140 = Py::Bool(matches!(v_a.clone(), Py::Dict(_))); if t7140.truthy() { t7140 } else { Py::Bool(matches!(v_b.clone(), Py::Dict(_))) } }).truthy() {
        if ({ let t7144 = py::not(&{ let t7141 = Py::Bool(matches!(v_a.clone(), Py::Dict(_))); if !t7141.truthy() { t7141 } else { Py::Bool(matches!(v_b.clone(), Py::Dict(_))) } }); if t7144.truthy() { t7144 } else { { let t7142 = py::method(&v_a.clone(), "keys", &[])?; Py::Bool({ let t7143 = py::method(&v_b.clone(), "keys", &[])?; if py::ne(&t7142, &t7143) { true } else { false } }) } } }).truthy() {
            return Ok(Py::Float(f64::INFINITY));
        }
        return Ok(py::minmax(&[py::add(&{ let mut out = vec![]; for t7145 in py::iter(&(v_a.clone()))? { let mut v_k: Py; 
v_k = t7145;
out.push(f__rel(py::getitem(&v_a.clone(), &v_k.clone())?, py::getitem(&v_b.clone(), &v_k.clone())?)?); } Py::list(out) }, &Py::list(vec![Py::Float(0.0)]))?], true)?);
    }
    if ({ let t7146 = Py::Bool(matches!(v_a.clone(), Py::List(_) | Py::Tuple(_))); if t7146.truthy() { t7146 } else { Py::Bool(matches!(v_b.clone(), Py::List(_) | Py::Tuple(_))) } }).truthy() {
        if ({ let t7150 = py::not(&{ let t7147 = Py::Bool(matches!(v_a.clone(), Py::List(_) | Py::Tuple(_))); if !t7147.truthy() { t7147 } else { Py::Bool(matches!(v_b.clone(), Py::List(_) | Py::Tuple(_))) } }); if t7150.truthy() { t7150 } else { { let t7148 = py::len(&v_a.clone())?; Py::Bool({ let t7149 = py::len(&v_b.clone())?; if py::ne(&t7148, &t7149) { true } else { false } }) } } }).truthy() {
            return Ok(Py::Float(f64::INFINITY));
        }
        return Ok(py::minmax(&[py::add(&{ let mut out = vec![]; for t7151 in py::iter(&(py::zip(&[v_a.clone(), v_b.clone()])?))? { let mut v_x: Py; let mut v_y: Py; 
let t7152 = (t7151).items()?;
if t7152.len() != 2 { return Err(PyErr::Value("unpack".into())); }
v_x = t7152[0].clone();
v_y = t7152[1].clone();
out.push(f__rel(v_x.clone(), v_y.clone())?); } Py::list(out) }, &Py::list(vec![Py::Float(0.0)]))?], true)?);
    }
    if ({ let t7159 = Py::Bool(matches!(v_a.clone(), Py::Str(_) | Py::Bool(_))); if t7159.truthy() { t7159 } else { { let t7158 = Py::Bool(matches!(v_b.clone(), Py::Str(_) | Py::Bool(_))); if t7158.truthy() { t7158 } else { { let t7157 = { let t7153 = v_a.clone(); Py::Bool({ let t7154 = Py::None; if py::is(&t7153, &t7154) { true } else { false } }) }; if t7157.truthy() { t7157 } else { { let t7155 = v_b.clone(); Py::Bool({ let t7156 = Py::None; if py::is(&t7155, &t7156) { true } else { false } }) } } } } } } }).truthy() {
        return Ok((if ({ let t7160 = v_a.clone(); Py::Bool({ let t7161 = v_b.clone(); if py::eq(&t7160, &t7161) { true } else { false } }) }).truthy() { Py::Float(0.0) } else { Py::Float(f64::INFINITY) }));
    }
    let t7163: Result<Option<Py>, PyErr> = (|| -> Result<Option<Py>, PyErr> {
            let t7162 = (Py::tuple(vec![py::float(&v_a.clone())?, py::float(&v_b.clone())?])).items()?;
            if t7162.len() != 2 { return Err(PyErr::Value("unpack".into())); }
            v_a = t7162[0].clone();
            v_b = t7162[1].clone();
            #[allow(unreachable_code)]
            Ok(None)
        })();
    match t7163 {
        Ok(Some(r)) => { return Ok(r); }
        Ok(None) => {}
        Err(e) if ((matches!(e, PyErr::Type(_)) || matches!(&e, PyErr::Other(c) if c.as_str() == "TypeError"))) || (matches!(e, PyErr::Value(_))) => {
            return Ok(Py::Float(f64::INFINITY));
        }
        Err(e) => return Err(e),
    }
    if (py::not(&{ let t7164 = py::math1("isfinite", &v_a.clone())?; if !t7164.truthy() { t7164 } else { py::math1("isfinite", &v_b.clone())? } })).truthy() {
        return Ok((if ({ let t7165 = v_a.clone(); Py::Bool({ let t7166 = v_b.clone(); if py::eq(&t7165, &t7166) { true } else { false } }) }).truthy() { Py::Float(0.0) } else { Py::Float(f64::INFINITY) }));
    }
    v_d = py::abs(&py::sub(&v_a.clone(), &v_b.clone())?)?;
    if ({ let t7167 = v_d.clone(); Py::Bool({ let t7168 = Py::Int(0i64); if py::eq(&t7167, &t7168) { true } else { false } }) }).truthy() {
        return Ok(Py::Float(0.0));
    }
    if ({ let t7173 = { let t7169 = py::abs(&v_a.clone())?; Py::Bool({ let t7170 = Py::Float(1e-15); if py::lt(&t7169, &t7170)? { true } else { false } }) }; if !t7173.truthy() { t7173 } else { { let t7171 = py::abs(&v_b.clone())?; Py::Bool({ let t7172 = Py::Float(1e-15); if py::lt(&t7171, &t7172)? { true } else { false } }) } } }).truthy() {
        return Ok(Py::Float(0.0));
    }
    return Ok(py::truediv(&v_d.clone(), &py::minmax(&[py::abs(&v_a.clone())?, py::abs(&v_b.clone())?], true)?)?);
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn f__revcomp(mut v_s: Py) -> R {
    return Ok(py::slice(&py::method(&v_s.clone(), "translate", &[c__COMP()?])?, None, None, Some((py::neg(&Py::Int(1i64))?).i()?))?);
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn f__same_loci(mut v_gs: Py) -> R {
    let mut v_g: Py = Py::None;
    let mut v_ref: Py = Py::None;
    v_ref = { let mut out = vec![]; for t905 in py::iter(&(py::getitem(&v_gs.clone(), &Py::Int(0i64))?))? { let mut v_p: Py; 
v_p = t905;
out.push(py::method(&py::getitem(&v_p.clone(), &Py::Int(0i64))?, "upper", &[])?); } Py::list(out) };
    for t906 in py::iter(&(py::slice(&v_gs.clone(), Some((Py::Int(1i64)).i()?), None, None)?))? {
        v_g = t906;
        let _ = f__need({ let t907 = { let mut out = vec![]; for t909 in py::iter(&(v_g.clone()))? { let mut v_p: Py; 
v_p = t909;
out.push(py::method(&py::getitem(&v_p.clone(), &Py::Int(0i64))?, "upper", &[])?); } Py::list(out) }; Py::Bool({ let t908 = v_ref.clone(); if py::eq(&t907, &t908) { true } else { false } }) }, Py::str("every genotype must list the same loci in the same order"))?;
    }
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn f__seq(mut v_s: Py, mut v_name: Py, mut v_kind: Py) -> R {
    let mut v_u: Py = Py::None;
    let _ = f__need({ let t6880 = Py::Bool(matches!(v_s.clone(), Py::Str(_))); if !t6880.truthy() { t6880 } else { { let t6878 = v_s.clone(); Py::Bool({ let t6879 = Py::str(""); if py::ne(&t6878, &t6879) { true } else { false } }) } } }, Py::str(&[py::format(&v_name.clone(), "")?, " must be a sequence of letters".to_string()].concat()))?;
    v_u = py::method(&v_s.clone(), "upper", &[])?;
    let _ = f__need({ let t6881 = py::len(&v_u.clone())?; Py::Bool({ let t6882 = Py::Int(100000i64); if py::le(&t6881, &t6882)? { true } else { false } }) }, Py::str(&[py::format(&v_name.clone(), "")?, ": at most 100000 letters".to_string()].concat()))?;
    if ({ let t6883 = v_kind.clone(); Py::Bool({ let t6884 = Py::str("DNA"); if py::eq(&t6883, &t6884) { true } else { false } }) }).truthy() {
        let _ = f__need({ let t6885 = ax_engine::pyre::module_call("fullmatch", &Py::str("[ACGT]+"), &[v_u.clone()], 0)?; Py::Bool({ let t6886 = Py::None; if !py::is(&t6885, &t6886) { true } else { false } }) }, Py::str(&[py::format(&v_name.clone(), "")?, " must be a DNA sequence of A, C, G and T only".to_string()].concat()))?;
    } else {
        if ({ let t6887 = v_kind.clone(); Py::Bool({ let t6888 = Py::str("RNA"); if py::eq(&t6887, &t6888) { true } else { false } }) }).truthy() {
            let _ = f__need({ let t6889 = ax_engine::pyre::module_call("fullmatch", &Py::str("[ACGU]+"), &[v_u.clone()], 0)?; Py::Bool({ let t6890 = Py::None; if !py::is(&t6889, &t6890) { true } else { false } }) }, Py::str(&[py::format(&v_name.clone(), "")?, " must be an RNA sequence of A, C, G and U only".to_string()].concat()))?;
        } else {
            if ({ let t6891 = v_kind.clone(); Py::Bool({ let t6892 = Py::str("NUC"); if py::eq(&t6891, &t6892) { true } else { false } }) }).truthy() {
                let _ = f__need({ let t6893 = ax_engine::pyre::module_call("fullmatch", &Py::str("[ACGT]+|[ACGU]+"), &[v_u.clone()], 0)?; Py::Bool({ let t6894 = Py::None; if !py::is(&t6893, &t6894) { true } else { false } }) }, Py::str(&[py::format(&v_name.clone(), "")?, " must be a DNA (A, C, G, T) or RNA (A, C, G, U) sequence, not a mixture".to_string()].concat()))?;
            } else {
                if ({ let t6895 = v_kind.clone(); Py::Bool({ let t6896 = Py::str("PROTEIN"); if py::eq(&t6895, &t6896) { true } else { false } }) }).truthy() {
                    let _ = f__need({ let t6897 = ax_engine::pyre::module_call("fullmatch", &py::add(&py::add(&Py::str("["), &c_AMINO()?)?, &Py::str("]+"))?, &[v_u.clone()], 0)?; Py::Bool({ let t6898 = Py::None; if !py::is(&t6897, &t6898) { true } else { false } }) }, Py::str(&[py::format(&v_name.clone(), "")?, " must be a protein sequence of the 20 standard one-letter amino-acid codes".to_string()].concat()))?;
                } else {
                    let _ = f__need({ let t6899 = ax_engine::pyre::module_call("fullmatch", &Py::str("[A-Z]+"), &[v_u.clone()], 0)?; Py::Bool({ let t6900 = Py::None; if !py::is(&t6899, &t6900) { true } else { false } }) }, Py::str(&[py::format(&v_name.clone(), "")?, " must be letters only".to_string()].concat()))?;
                }
            }
        }
    }
    return Ok(v_u.clone());
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn f__simplex_max(mut v_A: Py, mut v_b: Py, mut v_c: Py) -> R {
    let mut v_N: Py = Py::None;
    let mut v_T: Py = Py::None;
    let mut v_basis: Py = Py::None;
    let mut v_col: Py = Py::None;
    let mut v_cost: Py = Py::None;
    let mut v_i: Py = Py::None;
    let mut v_m: Py = Py::None;
    let mut v_n: Py = Py::None;
    let mut v_pivot: Py = Py::None;
    let mut v_run: Py = Py::None;
    let mut v_x: Py = Py::None;
    let t1546 = (Py::tuple(vec![py::len(&v_A.clone())?, py::len(&py::getitem(&v_A.clone(), &Py::Int(0i64))?)?])).items()?;
    if t1546.len() != 2 { return Err(PyErr::Value("unpack".into())); }
    v_m = t1546[0].clone();
    v_n = t1546[1].clone();
    v_A = { let mut out = vec![]; for t1547 in py::iter(&(v_A.clone()))? { let mut v_row: Py; 
v_row = t1547;
out.push({ let mut out = vec![]; for t1548 in py::iter(&(v_row.clone()))? { let mut v_x: Py; 
v_x = t1548;
out.push(py::fraction(&[v_x.clone()])?); } Py::list(out) }); } Py::list(out) };
    v_b = { let mut out = vec![]; for t1549 in py::iter(&(v_b.clone()))? { let mut v_x: Py; 
v_x = t1549;
out.push(py::fraction(&[v_x.clone()])?); } Py::list(out) };
    for t1550 in py::range_iter(&[(v_m.clone()).i()?]) {
        v_i = t1550;
        if ({ let t1551 = py::getitem(&v_b.clone(), &v_i.clone())?; Py::Bool({ let t1552 = Py::Int(0i64); if py::lt(&t1551, &t1552)? { true } else { false } }) }).truthy() {
            py::setitem(&v_A.clone(), &v_i.clone(), { let mut out = vec![]; for t1553 in py::iter(&(py::getitem(&v_A.clone(), &v_i.clone())?))? { let mut v_x: Py; 
v_x = t1553;
out.push(py::neg(&v_x.clone())?); } Py::list(out) })?;
            py::setitem(&v_b.clone(), &v_i.clone(), py::neg(&py::getitem(&v_b.clone(), &v_i.clone())?)?)?;
        }
    }
    v_T = { let mut out = vec![]; for t1554 in py::range_iter(&[(v_m.clone()).i()?]) { let mut v_i: Py; 
v_i = t1554;
out.push(py::add(&py::add(&py::getitem(&v_A.clone(), &v_i.clone())?, &{ let mut out = vec![]; for t1555 in py::range_iter(&[(v_m.clone()).i()?]) { let mut v_k: Py; 
v_k = t1555;
out.push(py::fraction(&[(if ({ let t1556 = v_k.clone(); Py::Bool({ let t1557 = v_i.clone(); if py::eq(&t1556, &t1557) { true } else { false } }) }).truthy() { Py::Int(1i64) } else { Py::Int(0i64) })])?); } Py::list(out) })?, &Py::list(vec![py::getitem(&v_b.clone(), &v_i.clone())?]))?); } Py::list(out) };
    v_basis = { let mut out = vec![]; for t1558 in py::range_iter(&[(v_m.clone()).i()?]) { let mut v_i: Py; 
v_i = t1558;
out.push(py::add(&v_n.clone(), &v_i.clone())?); } Py::list(out) };
    v_N = py::add(&v_n.clone(), &v_m.clone())?;
    v_pivot = { let v_A = v_A.clone(); let v_N = v_N.clone(); let v_T = v_T.clone(); let v_b = v_b.clone(); let v_basis = v_basis.clone(); let v_c = v_c.clone(); let v_i = v_i.clone(); let v_m = v_m.clone(); let v_n = v_n.clone(); let v_pivot = v_pivot.clone(); py::func(move |a: &[Py]| -> R {
            #[allow(unused_mut, unused_variables, unused_assignments)]
{
            let mut v_r: Py = a.get(0).cloned().unwrap_or(Py::None);
            let mut v_col: Py = a.get(1).cloned().unwrap_or(Py::None);
            let mut v_f: Py = Py::None;
            let mut v_i: Py = Py::None;
            let mut v_pv: Py = Py::None;
            v_pv = py::getitem(&py::getitem(&v_T.clone(), &v_r.clone())?, &v_col.clone())?;
            py::setitem(&v_T.clone(), &v_r.clone(), { let mut out = vec![]; for t1559 in py::iter(&(py::getitem(&v_T.clone(), &v_r.clone())?))? { let mut v_x: Py; 
v_x = t1559;
out.push(py::truediv(&v_x.clone(), &v_pv.clone())?); } Py::list(out) })?;
            for t1560 in py::range_iter(&[(v_m.clone()).i()?]) {
                v_i = t1560;
                if ({ let t1565 = { let t1561 = v_i.clone(); Py::Bool({ let t1562 = v_r.clone(); if py::ne(&t1561, &t1562) { true } else { false } }) }; if !t1565.truthy() { t1565 } else { { let t1563 = py::getitem(&py::getitem(&v_T.clone(), &v_i.clone())?, &v_col.clone())?; Py::Bool({ let t1564 = Py::Int(0i64); if py::ne(&t1563, &t1564) { true } else { false } }) } } }).truthy() {
                    v_f = py::getitem(&py::getitem(&v_T.clone(), &v_i.clone())?, &v_col.clone())?;
                    py::setitem(&v_T.clone(), &v_i.clone(), { let mut out = vec![]; for t1566 in py::iter(&(py::zip(&[py::getitem(&v_T.clone(), &v_i.clone())?, py::getitem(&v_T.clone(), &v_r.clone())?])?))? { let mut v_x: Py; let mut v_y: Py; 
let t1567 = (t1566).items()?;
if t1567.len() != 2 { return Err(PyErr::Value("unpack".into())); }
v_x = t1567[0].clone();
v_y = t1567[1].clone();
out.push(py::sub(&v_x.clone(), &py::mul(&v_f.clone(), &v_y.clone())?)?); } Py::list(out) })?;
                }
            }
            py::setitem(&v_basis.clone(), &v_r.clone(), v_col.clone())?;
            #[allow(unreachable_code)]
            Ok(Py::None)
        }
        }) };
    v_run = { let v_A = v_A.clone(); let v_N = v_N.clone(); let v_T = v_T.clone(); let v_b = v_b.clone(); let v_basis = v_basis.clone(); let v_c = v_c.clone(); let v_i = v_i.clone(); let v_m = v_m.clone(); let v_n = v_n.clone(); let v_pivot = v_pivot.clone(); let v_run = v_run.clone(); py::func(move |a: &[Py]| -> R {
            #[allow(unused_mut, unused_variables, unused_assignments)]
{
            let mut v_cost: Py = a.get(0).cloned().unwrap_or(Py::None);
            let mut v_allowed: Py = a.get(1).cloned().unwrap_or(Py::None);
            let mut v_enter: Py = Py::None;
            let mut v_r: Py = Py::None;
            let mut v_red: Py = Py::None;
            let mut v_rows: Py = Py::None;
            while (Py::Bool(true)).truthy() {
                v_red = { let mut out = vec![]; for t1568 in py::range_iter(&[(v_N.clone()).i()?]) { let mut v_j: Py; 
v_j = t1568;
out.push(py::sub(&py::getitem(&v_cost.clone(), &v_j.clone())?, &py::sum(&{ let mut out = vec![]; for t1569 in py::range_iter(&[(v_m.clone()).i()?]) { let mut v_i: Py; 
v_i = t1569;
out.push(py::mul(&py::getitem(&v_cost.clone(), &py::getitem(&v_basis.clone(), &v_i.clone())?)?, &py::getitem(&py::getitem(&v_T.clone(), &v_i.clone())?, &v_j.clone())?)?); } Py::list(out) }, None)?)?); } Py::list(out) };
                v_enter = (match ('g1570: { let mut out = vec![]; for t1571 in py::range_iter(&[(v_N.clone()).i()?]) { let mut v_j: Py; 
v_j = t1571;
if ({ let t1574 = py::call(&v_allowed.clone(), &[v_j.clone()])?; if !t1574.truthy() { t1574 } else { { let t1572 = py::getitem(&v_red.clone(), &v_j.clone())?; Py::Bool({ let t1573 = Py::Int(0i64); if py::gt(&t1572, &t1573)? { true } else { false } }) } } }).truthy() { out.push(v_j.clone()); break 'g1570 Py::list(out); }} Py::list(out) }).items()?.into_iter().next() { Some(x) => x, None => Py::None });
                if ({ let t1575 = v_enter.clone(); Py::Bool({ let t1576 = Py::None; if py::is(&t1575, &t1576) { true } else { false } }) }).truthy() {
                    return Ok(Py::None);
                }
                v_rows = { let mut out = vec![]; for t1577 in py::range_iter(&[(v_m.clone()).i()?]) { let mut v_i: Py; 
v_i = t1577;
if ({ let t1578 = py::getitem(&py::getitem(&v_T.clone(), &v_i.clone())?, &v_enter.clone())?; Py::Bool({ let t1579 = Py::Int(0i64); if py::gt(&t1578, &t1579)? { true } else { false } }) }).truthy() { out.push(Py::tuple(vec![py::truediv(&py::getitem(&py::getitem(&v_T.clone(), &v_i.clone())?, &py::neg(&Py::Int(1i64))?)?, &py::getitem(&py::getitem(&v_T.clone(), &v_i.clone())?, &v_enter.clone())?)?, py::getitem(&v_basis.clone(), &v_i.clone())?, v_i.clone()])); }} Py::list(out) };
                if (py::not(&v_rows.clone())).truthy() {
                    return Err(py::value_error(py::to_str(&Py::str("the objective is unbounded"))));
                }
                v_r = py::getitem(&py::minmax(&[v_rows.clone()], false)?, &Py::Int(2i64))?;
                let _ = py::call(&v_pivot.clone(), &[v_r.clone(), v_enter.clone()])?;
            }
            #[allow(unreachable_code)]
            Ok(Py::None)
        }
        }) };
    let _ = py::call(&v_run.clone(), &[py::add(&py::mul(&Py::list(vec![py::fraction(&[Py::Int(0i64)])?]), &v_n.clone())?, &py::mul(&Py::list(vec![py::fraction(&[py::neg(&Py::Int(1i64))?])?]), &v_m.clone())?)?, { let v_A = v_A.clone(); let v_N = v_N.clone(); let v_T = v_T.clone(); let v_b = v_b.clone(); let v_basis = v_basis.clone(); let v_c = v_c.clone(); let v_i = v_i.clone(); let v_m = v_m.clone(); let v_n = v_n.clone(); let v_pivot = v_pivot.clone(); let v_run = v_run.clone(); py::func(move |a: &[Py]| -> R {
            #[allow(unused_mut, unused_variables, unused_assignments)]
{
            let mut v_j: Py = a.get(0).cloned().unwrap_or(Py::None);
            return Ok(Py::Bool(true));
            #[allow(unreachable_code)]
            Ok(Py::None)
        }
        }) }])?;
    if ({ let t1580 = py::sum(&{ let mut out = vec![]; for t1582 in py::range_iter(&[(v_m.clone()).i()?]) { let mut v_i: Py; 
v_i = t1582;
if ({ let t1583 = py::getitem(&v_basis.clone(), &v_i.clone())?; Py::Bool({ let t1584 = v_n.clone(); if py::ge(&t1583, &t1584)? { true } else { false } }) }).truthy() { out.push(py::getitem(&py::getitem(&v_T.clone(), &v_i.clone())?, &py::neg(&Py::Int(1i64))?)?); }} Py::list(out) }, None)?; Py::Bool({ let t1581 = Py::Int(0i64); if py::ne(&t1580, &t1581) { true } else { false } }) }).truthy() {
        return Err(py::value_error(py::to_str(&Py::str("no flux distribution satisfies the steady state and the bounds (infeasible)"))));
    }
    for t1585 in py::range_iter(&[(v_m.clone()).i()?]) {
        v_i = t1585;
        if ({ let t1586 = py::getitem(&v_basis.clone(), &v_i.clone())?; Py::Bool({ let t1587 = v_n.clone(); if py::ge(&t1586, &t1587)? { true } else { false } }) }).truthy() {
            v_col = (match ('g1588: { let mut out = vec![]; for t1589 in py::range_iter(&[(v_n.clone()).i()?]) { let mut v_j: Py; 
v_j = t1589;
if ({ let t1590 = py::getitem(&py::getitem(&v_T.clone(), &v_i.clone())?, &v_j.clone())?; Py::Bool({ let t1591 = Py::Int(0i64); if py::ne(&t1590, &t1591) { true } else { false } }) }).truthy() { out.push(v_j.clone()); break 'g1588 Py::list(out); }} Py::list(out) }).items()?.into_iter().next() { Some(x) => x, None => Py::None });
            if ({ let t1592 = v_col.clone(); Py::Bool({ let t1593 = Py::None; if !py::is(&t1592, &t1593) { true } else { false } }) }).truthy() {
                let _ = py::call(&v_pivot.clone(), &[v_i.clone(), v_col.clone()])?;
            }
        }
    }
    v_cost = py::add(&{ let mut out = vec![]; for t1594 in py::iter(&(v_c.clone()))? { let mut v_x: Py; 
v_x = t1594;
out.push(py::fraction(&[v_x.clone()])?); } Py::list(out) }, &py::mul(&Py::list(vec![py::fraction(&[Py::Int(0i64)])?]), &v_m.clone())?)?;
    let _ = py::call(&v_run.clone(), &[v_cost.clone(), { let v_A = v_A.clone(); let v_N = v_N.clone(); let v_T = v_T.clone(); let v_b = v_b.clone(); let v_basis = v_basis.clone(); let v_c = v_c.clone(); let v_col = v_col.clone(); let v_cost = v_cost.clone(); let v_i = v_i.clone(); let v_m = v_m.clone(); let v_n = v_n.clone(); let v_pivot = v_pivot.clone(); let v_run = v_run.clone(); py::func(move |a: &[Py]| -> R {
            #[allow(unused_mut, unused_variables, unused_assignments)]
{
            let mut v_j: Py = a.get(0).cloned().unwrap_or(Py::None);
            return Ok({ let t1595 = v_j.clone(); Py::Bool({ let t1596 = v_n.clone(); if py::lt(&t1595, &t1596)? { true } else { false } }) });
            #[allow(unreachable_code)]
            Ok(Py::None)
        }
        }) }])?;
    v_x = py::mul(&Py::list(vec![py::fraction(&[Py::Int(0i64)])?]), &v_N.clone())?;
    for t1597 in py::range_iter(&[(v_m.clone()).i()?]) {
        v_i = t1597;
        py::setitem(&v_x.clone(), &py::getitem(&v_basis.clone(), &v_i.clone())?, py::getitem(&py::getitem(&v_T.clone(), &v_i.clone())?, &py::neg(&Py::Int(1i64))?)?)?;
    }
    return Ok(py::sum(&{ let mut out = vec![]; for t1598 in py::iter(&(py::zip(&[v_c.clone(), py::slice(&v_x.clone(), None, Some((v_n.clone()).i()?), None)?])?))? { let mut v_cf: Py; let mut v_xv: Py; 
let t1599 = (t1598).items()?;
if t1599.len() != 2 { return Err(PyErr::Value("unpack".into())); }
v_cf = t1599[0].clone();
v_xv = t1599[1].clone();
out.push(py::mul(&v_cf.clone(), &v_xv.clone())?); } Py::list(out) }, None)?);
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn f__solve3(mut v_A: Py, mut v_b: Py) -> R {
    let mut v_M: Py = Py::None;
    let mut v_c: Py = Py::None;
    let mut v_f: Py = Py::None;
    let mut v_k: Py = Py::None;
    let mut v_p: Py = Py::None;
    let mut v_r: Py = Py::None;
    let mut v_x: Py = Py::None;
    v_M = { let mut out = vec![]; for t2701 in py::range_iter(&[(Py::Int(3i64)).i()?]) { let mut v_i: Py; 
v_i = t2701;
out.push(py::add(&Py::list((py::getitem(&v_A.clone(), &v_i.clone())?).items()?), &Py::list(vec![py::getitem(&v_b.clone(), &v_i.clone())?]))?); } Py::list(out) };
    for t2702 in py::range_iter(&[(Py::Int(3i64)).i()?]) {
        v_c = t2702;
        v_p = py::minmax_key(&[Py::list((py::range_iter(&[(v_c.clone()).i()?, (Py::Int(3i64)).i()?])).into_iter().collect::<Vec<Py>>())], &{ let v_A = v_A.clone(); let v_M = v_M.clone(); let v_b = v_b.clone(); let v_c = v_c.clone(); py::func(move |a: &[Py]| -> R {
            #[allow(unused_mut, unused_variables, unused_assignments)]
{
            let mut v_r: Py = a.get(0).cloned().unwrap_or(Py::None);
            return Ok(py::abs(&py::getitem(&py::getitem(&v_M.clone(), &v_r.clone())?, &v_c.clone())?)?);
            #[allow(unreachable_code)]
            Ok(Py::None)
        }
        }) }, true)?;
        if ({ let t2703 = py::abs(&py::getitem(&py::getitem(&v_M.clone(), &v_p.clone())?, &v_c.clone())?)?; Py::Bool({ let t2704 = Py::Float(1e-300); if py::lt(&t2703, &t2704)? { true } else { false } }) }).truthy() {
            return Err(py::value_error(py::to_str(&Py::str("singular normal equations"))));
        }
        let t2705 = (Py::tuple(vec![py::getitem(&v_M.clone(), &v_p.clone())?, py::getitem(&v_M.clone(), &v_c.clone())?])).items()?;
        if t2705.len() != 2 { return Err(PyErr::Value("unpack".into())); }
        py::setitem(&v_M.clone(), &v_c.clone(), t2705[0].clone())?;
        py::setitem(&v_M.clone(), &v_p.clone(), t2705[1].clone())?;
        for t2706 in py::range_iter(&[(py::add(&v_c.clone(), &Py::Int(1i64))?).i()?, (Py::Int(3i64)).i()?]) {
            v_r = t2706;
            v_f = py::truediv(&py::getitem(&py::getitem(&v_M.clone(), &v_r.clone())?, &v_c.clone())?, &py::getitem(&py::getitem(&v_M.clone(), &v_c.clone())?, &v_c.clone())?)?;
            for t2707 in py::range_iter(&[(v_c.clone()).i()?, (Py::Int(4i64)).i()?]) {
                v_k = t2707;
                let t2708 = py::getitem(&v_M.clone(), &v_r.clone())?; let t2709 = v_k.clone();
                py::setitem(&t2708, &t2709, py::sub(&py::getitem(&t2708, &t2709)?, &py::mul(&v_f.clone(), &py::getitem(&py::getitem(&v_M.clone(), &v_c.clone())?, &v_k.clone())?)?)?)?;
            }
        }
    }
    v_x = py::mul(&Py::list(vec![Py::Float(0.0)]), &Py::Int(3i64))?;
    for t2710 in py::iter(&(Py::tuple(vec![Py::Int(2i64), Py::Int(1i64), Py::Int(0i64)])))? {
        v_r = t2710;
        py::setitem(&v_x.clone(), &v_r.clone(), py::truediv(&py::sub(&py::getitem(&py::getitem(&v_M.clone(), &v_r.clone())?, &Py::Int(3i64))?, &py::sum(&{ let mut out = vec![]; for t2711 in py::range_iter(&[(py::add(&v_r.clone(), &Py::Int(1i64))?).i()?, (Py::Int(3i64)).i()?]) { let mut v_k: Py; 
v_k = t2711;
out.push(py::mul(&py::getitem(&py::getitem(&v_M.clone(), &v_r.clone())?, &v_k.clone())?, &py::getitem(&v_x.clone(), &v_k.clone())?)?); } Py::list(out) }, None)?)?, &py::getitem(&py::getitem(&v_M.clone(), &v_r.clone())?, &v_r.clone())?)?)?;
    }
    return Ok(v_x.clone());
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn f__splits(mut v_root: Py, mut v_leaves: Py) -> R {
    let mut v_allx: Py = Py::None;
    let mut v_first: Py = Py::None;
    let mut v_out: Py = Py::None;
    let mut v_walk: Py = Py::None;
    v_allx = py::set(&v_leaves.clone())?;
    v_first = py::minmax(&[v_leaves.clone()], false)?;
    v_out = py::set(&Py::list(vec![]))?;
    let rec6024 = py::rec_cell();
    v_walk = { let v_allx = v_allx.clone(); let v_first = v_first.clone(); let v_leaves = v_leaves.clone(); let v_out = v_out.clone(); let v_root = v_root.clone(); let rec6024 = rec6024.clone(); py::func(move |a: &[Py]| -> R {
            #[allow(unused_mut, unused_variables, unused_assignments)]
{
            let mut v_t: Py = a.get(0).cloned().unwrap_or(Py::None);
            let v_walk: Py = py::rec_get(&rec6024)?;
            let mut v_s: Py = Py::None;
            let mut v_side: Py = Py::None;
            if (py::not(&py::getitem(&v_t.clone(), &Py::Int(2i64))?)).truthy() {
                return Ok(py::set(&Py::list(vec![py::getitem(&v_t.clone(), &Py::Int(0i64))?]))?);
            }
            v_s = py::method(&py::set(&Py::list(vec![]))?, "union", &({ let mut v: Vec<Py> = vec![]; v.extend(({ let mut out = vec![]; for t6025 in py::iter(&(py::getitem(&v_t.clone(), &Py::Int(2i64))?))? { let mut v_k: Py; 
v_k = t6025;
out.push(py::call(&v_walk.clone(), &[v_k.clone()])?); } Py::list(out) }).items()?); v })[..])?;
            v_side = (if ({ let t6026 = v_first.clone(); Py::Bool({ let t6027 = v_s.clone(); if py::contains(&t6027, &t6026)? { true } else { false } }) }).truthy() { py::sub(&v_allx.clone(), &v_s.clone())? } else { v_s.clone() });
            if ({ let t6028 = Py::Int(2i64); Py::Bool({ let t6029 = py::len(&v_side.clone())?; if py::le(&t6028, &t6029)? { { let t6030 = py::sub(&py::len(&v_allx.clone())?, &Py::Int(2i64))?; if py::le(&t6029, &t6030)? { true } else { false } } } else { false } }) }).truthy() {
                let _ = py::method(&v_out.clone(), "add", &[v_side.clone()])?;
            }
            return Ok(v_s.clone());
            #[allow(unreachable_code)]
            Ok(Py::None)
        }
        }) };
    py::rec_set(&rec6024, &v_walk);
    let _ = py::call(&v_walk.clone(), &[v_root.clone()])?;
    return Ok(v_out.clone());
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn f__states(mut v_text: Py, mut v_need: Py) -> R {
    let mut v_lens: Py = Py::None;
    let mut v_m: Py = Py::None;
    let mut v_out: Py = Py::None;
    let mut v_tok: Py = Py::None;
    let _ = f__need(Py::Bool(matches!(v_text.clone(), Py::Str(_))), Py::str("states must be written as name:SEQUENCE pairs"))?;
    v_out = Py::dict(vec![]);
    for t1439 in py::iter(&(py::method(&v_text.clone(), "split", &[])?))? {
        v_tok = t1439;
        v_m = ax_engine::pyre::module_call("fullmatch", &Py::str("([A-Za-z][A-Za-z0-9_]*):([A-Za-z]+)"), &[v_tok.clone()], 0)?;
        let _ = f__need({ let t1440 = v_m.clone(); Py::Bool({ let t1441 = Py::None; if !py::is(&t1440, &t1441) { true } else { false } }) }, Py::str(&["cannot read '".to_string(), py::format(&v_tok.clone(), "")?, "': write name:SEQUENCE".to_string()].concat()))?;
        let _ = f__need({ let t1442 = py::method(&v_m.clone(), "group", &[Py::Int(1i64)])?; Py::Bool({ let t1443 = v_out.clone(); if !py::contains(&t1443, &t1442)? { true } else { false } }) }, Py::str(&[py::format(&py::method(&v_m.clone(), "group", &[Py::Int(1i64)])?, "")?, " is given twice".to_string()].concat()))?;
        py::setitem(&v_out.clone(), &py::method(&v_m.clone(), "group", &[Py::Int(1i64)])?, py::method(&py::method(&v_m.clone(), "group", &[Py::Int(2i64)])?, "upper", &[])?)?;
    }
    let _ = f__need({ let t1444 = py::set(&v_out.clone())?; Py::Bool({ let t1445 = py::set(&v_need.clone())?; if py::eq(&t1444, &t1445) { true } else { false } }) }, Py::str("give one sequence for every leaf of the tree, and no others"))?;
    v_lens = py::set(&{ let mut out = vec![]; for t1446 in py::iter(&(py::method(&v_out.clone(), "values", &[])?))? { let mut v_v: Py; 
v_v = t1446;
out.push(py::len(&v_v.clone())?); } Py::list(out) })?;
    let _ = f__need({ let t1447 = py::len(&v_lens.clone())?; Py::Bool({ let t1448 = Py::Int(1i64); if py::eq(&t1447, &t1448) { true } else { false } }) }, Py::str("all sequences must have the same length (aligned)"))?;
    return Ok(v_out.clone());
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn f__tie(mut v_vals: Py, mut v_best: Py, mut v_scale: Py) -> R {
    return Ok({ let t7026 = py::sum(&{ let mut out = vec![]; for t7028 in py::iter(&(v_vals.clone()))? { let mut v_v: Py; 
v_v = t7028;
if ({ let t7029 = py::abs(&py::sub(&v_v.clone(), &v_best.clone())?)?; Py::Bool({ let t7030 = py::mul(&Py::Float(1e-12), &v_scale.clone())?; if py::le(&t7029, &t7030)? { true } else { false } }) }).truthy() { out.push(Py::Int(1i64)); }} Py::list(out) }, None)?; Py::Bool({ let t7027 = Py::Int(1i64); if py::gt(&t7026, &t7027)? { true } else { false } }) });
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn f__tm_nn(mut v_s: Py, mut v_Na: Py, mut v_dnac1: Py, mut v_dnac2: Py, mut v_saltcorr: Py) -> R {
    let mut v_at: Py = Py::None;
    let mut v_c: Py = Py::None;
    let mut v_dh: Py = Py::None;
    let mut v_ds: Py = Py::None;
    let mut v_ends: Py = Py::None;
    let mut v_gc: Py = Py::None;
    let mut v_i: Py = Py::None;
    let mut v_k: Py = Py::None;
    let mut v_key: Py = Py::None;
    let mut v_nb: Py = Py::None;
    let mut v_selfcomp: Py = Py::None;
    let mut v_t: Py = Py::None;
    let mut v_v: Py = Py::None;
    v_c = py::method(&v_s.clone(), "translate", &[c__COMP()?])?;
    let t6694 = (Py::tuple(vec![Py::Float(0.0), Py::Float(0.0)])).items()?;
    if t6694.len() != 2 { return Err(PyErr::Value("unpack".into())); }
    v_dh = t6694[0].clone();
    v_ds = t6694[1].clone();
    v_t = c_DNA_NN3()?;
    v_dh = py::iadd(&v_dh, &py::getitem(&py::getitem(&v_t.clone(), &Py::str("init"))?, &Py::Int(0i64))?)?;
    v_ds = py::iadd(&v_ds, &py::getitem(&py::getitem(&v_t.clone(), &Py::str("init"))?, &Py::Int(1i64))?)?;
    v_key = (if ({ let t6699 = { let t6695 = Py::str("G"); Py::Bool({ let t6696 = v_s.clone(); if py::contains(&t6696, &t6695)? { true } else { false } }) }; if t6699.truthy() { t6699 } else { { let t6697 = Py::str("C"); Py::Bool({ let t6698 = v_s.clone(); if py::contains(&t6698, &t6697)? { true } else { false } }) } } }).truthy() { Py::str("init_oneG/C") } else { Py::str("init_allA/T") });
    v_dh = py::iadd(&v_dh, &py::getitem(&py::getitem(&v_t.clone(), &v_key.clone())?, &Py::Int(0i64))?)?;
    v_ds = py::iadd(&v_ds, &py::getitem(&py::getitem(&v_t.clone(), &v_key.clone())?, &Py::Int(1i64))?)?;
    if (py::method(&v_s.clone(), "startswith", &[Py::str("T")])?).truthy() {
        v_dh = py::iadd(&v_dh, &py::getitem(&py::getitem(&v_t.clone(), &Py::str("init_5T/A"))?, &Py::Int(0i64))?)?;
        v_ds = py::iadd(&v_ds, &py::getitem(&py::getitem(&v_t.clone(), &Py::str("init_5T/A"))?, &Py::Int(1i64))?)?;
    }
    if (py::method(&v_s.clone(), "endswith", &[Py::str("A")])?).truthy() {
        v_dh = py::iadd(&v_dh, &py::getitem(&py::getitem(&v_t.clone(), &Py::str("init_5T/A"))?, &Py::Int(0i64))?)?;
        v_ds = py::iadd(&v_ds, &py::getitem(&py::getitem(&v_t.clone(), &Py::str("init_5T/A"))?, &Py::Int(1i64))?)?;
    }
    v_ends = py::add(&py::getitem(&v_s.clone(), &Py::Int(0i64))?, &py::getitem(&v_s.clone(), &py::neg(&Py::Int(1i64))?)?)?;
    v_at = py::add(&py::method(&v_ends.clone(), "count", &[Py::str("A")])?, &py::method(&v_ends.clone(), "count", &[Py::str("T")])?)?;
    v_gc = py::add(&py::method(&v_ends.clone(), "count", &[Py::str("G")])?, &py::method(&v_ends.clone(), "count", &[Py::str("C")])?)?;
    v_dh = py::iadd(&v_dh, &py::mul(&py::getitem(&py::getitem(&v_t.clone(), &Py::str("init_A/T"))?, &Py::Int(0i64))?, &v_at.clone())?)?;
    v_ds = py::iadd(&v_ds, &py::mul(&py::getitem(&py::getitem(&v_t.clone(), &Py::str("init_A/T"))?, &Py::Int(1i64))?, &v_at.clone())?)?;
    v_dh = py::iadd(&v_dh, &py::mul(&py::getitem(&py::getitem(&v_t.clone(), &Py::str("init_G/C"))?, &Py::Int(0i64))?, &v_gc.clone())?)?;
    v_ds = py::iadd(&v_ds, &py::mul(&py::getitem(&py::getitem(&v_t.clone(), &Py::str("init_G/C"))?, &Py::Int(1i64))?, &v_gc.clone())?)?;
    for t6700 in py::range_iter(&[(py::sub(&py::len(&v_s.clone())?, &Py::Int(1i64))?).i()?]) {
        v_i = t6700;
        v_nb = py::add(&py::add(&py::slice(&v_s.clone(), Some((v_i.clone()).i()?), Some((py::add(&v_i.clone(), &Py::Int(2i64))?).i()?), None)?, &Py::str("/"))?, &py::slice(&v_c.clone(), Some((v_i.clone()).i()?), Some((py::add(&v_i.clone(), &Py::Int(2i64))?).i()?), None)?)?;
        v_v = (if ({ let t6701 = v_nb.clone(); Py::Bool({ let t6702 = v_t.clone(); if py::contains(&t6702, &t6701)? { true } else { false } }) }).truthy() { py::getitem(&v_t.clone(), &v_nb.clone())? } else { py::getitem(&v_t.clone(), &py::slice(&v_nb.clone(), None, None, Some((py::neg(&Py::Int(1i64))?).i()?))?)? });
        v_dh = py::iadd(&v_dh, &py::getitem(&v_v.clone(), &Py::Int(0i64))?)?;
        v_ds = py::iadd(&v_ds, &py::getitem(&v_v.clone(), &Py::Int(1i64))?)?;
    }
    v_selfcomp = { let t6703 = v_s.clone(); Py::Bool({ let t6704 = f__revcomp(v_s.clone())?; if py::eq(&t6703, &t6704) { true } else { false } }) };
    v_k = py::mul(&py::sub(&v_dnac1.clone(), &py::truediv(&v_dnac2.clone(), &Py::Float(2.0))?)?, &Py::Float(1e-09))?;
    if (v_selfcomp.clone()).truthy() {
        v_k = py::mul(&v_dnac1.clone(), &Py::Float(1e-09))?;
        v_dh = py::iadd(&v_dh, &py::getitem(&py::getitem(&v_t.clone(), &Py::str("sym"))?, &Py::Int(0i64))?)?;
        v_ds = py::iadd(&v_ds, &py::getitem(&py::getitem(&v_t.clone(), &Py::str("sym"))?, &Py::Int(1i64))?)?;
    }
    if (v_saltcorr.clone()).truthy() {
        v_ds = py::iadd(&v_ds, &py::mul(&py::mul(&Py::Float(0.368), &py::sub(&py::len(&v_s.clone())?, &Py::Int(1i64))?)?, &py::math1("log", &py::mul(&v_Na.clone(), &Py::Float(0.001))?)?)?)?;
    }
    return Ok(py::sub(&py::truediv(&py::mul(&Py::Int(1000i64), &v_dh.clone())?, &py::add(&v_ds.clone(), &py::mul(&Py::Float(1.987), &py::math1("log", &v_k.clone())?)?)?)?, &Py::Float(273.15))?);
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn f__translate(mut v_s: Py, mut v_table: Py, mut v_truncate: Py) -> R {
    let mut v_p: Py = Py::None;
    v_p = py::method(&Py::str(""), "join", &[{ let mut out = vec![]; for t6905 in py::range_iter(&[(Py::Int(0i64)).i()?, (py::len(&v_s.clone())?).i()?, (Py::Int(3i64)).i()?]) { let mut v_i: Py; 
v_i = t6905;
out.push(py::getitem(&v_table.clone(), &py::slice(&v_s.clone(), Some((v_i.clone()).i()?), Some((py::add(&v_i.clone(), &Py::Int(3i64))?).i()?), None)?)?); } Py::list(out) }])?;
    if ({ let t6908 = v_truncate.clone(); if !t6908.truthy() { t6908 } else { { let t6906 = Py::str("*"); Py::Bool({ let t6907 = v_p.clone(); if py::contains(&t6907, &t6906)? { true } else { false } }) } } }).truthy() {
        v_p = py::slice(&v_p.clone(), None, Some((py::method(&v_p.clone(), "index", &[Py::str("*")])?).i()?), None)?;
    }
    return Ok(v_p.clone());
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn f__whole(mut v_name: Py, mut v_v: Py, mut v_lo: Py) -> R {
    let _ = f__need({ let t6424 = Py::Bool(matches!(v_v.clone(), Py::Int(_) | Py::Big(_) | Py::Bool(_) | Py::Float(_))); if !t6424.truthy() { t6424 } else { { let t6423 = py::not(&Py::Bool(matches!(v_v.clone(), Py::Bool(_)))); if !t6423.truthy() { t6423 } else { { let t6422 = py::math1("isfinite", &v_v.clone())?; if !t6422.truthy() { t6422 } else { { let t6421 = { let t6417 = py::float(&v_v.clone())?; Py::Bool({ let t6418 = py::int(&v_v.clone())?; if py::eq(&t6417, &t6418) { true } else { false } }) }; if !t6421.truthy() { t6421 } else { { let t6419 = v_v.clone(); Py::Bool({ let t6420 = v_lo.clone(); if py::ge(&t6419, &t6420)? { true } else { false } }) } } } } } } } } }, Py::str(&[py::format(&v_name.clone(), "")?, " must be a whole number of at least ".to_string(), py::format(&v_lo.clone(), "")?].concat()))?;
    return Ok(py::int(&v_v.clone())?);
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn f_aa_composition(mut v_seq: Py) -> R {
    let mut v_pct: Py = Py::None;
    let mut v_s: Py = Py::None;
    let mut v_top: Py = Py::None;
    v_s = f__seq(v_seq.clone(), Py::str("seq"), Py::str("PROTEIN"))?;
    v_pct = { let mut out = vec![]; for t100 in py::iter(&(c_AMINO()?))? { let mut v_a: Py; 
v_a = t100;
out.push((v_a.clone(), py::truediv(&py::mul(&Py::Float(100.0), &py::method(&v_s.clone(), "count", &[v_a.clone()])?)?, &py::len(&v_s.clone())?)?)); } Py::dict(out) };
    v_top = py::minmax(&[py::method(&v_pct.clone(), "values", &[])?], true)?;
    return Ok(f__out(Py::dict(vec![(Py::str("percent"), v_pct.clone()), (Py::str("most_common"), py::minmax(&[{ let mut out = vec![]; for t115 in py::iter(&(c_AMINO()?))? { let mut v_a: Py; 
v_a = t115;
if ({ let t116 = py::getitem(&v_pct.clone(), &v_a.clone())?; Py::Bool({ let t117 = v_top.clone(); if py::eq(&t116, &t117) { true } else { false } }) }).truthy() { out.push(v_a.clone()); }} Py::list(out) }], false)?), (Py::str("length"), py::len(&v_s.clone())?)]))?);
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn f_ace(mut v_counts: Py, mut v_rare_threshold: Py, mut v_no_gamma: Py) -> R {
    let mut v_c: Py = Py::None;
    let mut v_c_ace: Py = Py::None;
    let mut v_gamma: Py = Py::None;
    let mut v_n_rare: Py = Py::None;
    let mut v_rare: Py = Py::None;
    let mut v_s_abun: Py = Py::None;
    let mut v_s_rare: Py = Py::None;
    let mut v_singles: Py = Py::None;
    let mut v_top: Py = Py::None;
    v_c = f__abund(v_counts.clone(), Py::str("counts"))?;
    let _ = f__need({ let t21 = Py::Bool(matches!(v_rare_threshold.clone(), Py::Int(_) | Py::Big(_) | Py::Bool(_))); if !t21.truthy() { t21 } else { { let t19 = v_rare_threshold.clone(); Py::Bool({ let t20 = Py::Int(1i64); if py::ge(&t19, &t20)? { true } else { false } }) } } }, Py::str("rare_threshold must be a whole number >= 1"))?;
    v_rare = { let mut out = vec![]; for t22 in py::iter(&(v_c.clone()))? { let mut v_x: Py; 
v_x = t22;
if ({ let t23 = v_x.clone(); Py::Bool({ let t24 = v_rare_threshold.clone(); if py::le(&t23, &t24)? { true } else { false } }) }).truthy() { out.push(v_x.clone()); }} Py::list(out) };
    v_s_rare = py::len(&v_rare.clone())?;
    v_singles = py::sum(&{ let mut out = vec![]; for t25 in py::iter(&(v_c.clone()))? { let mut v_x: Py; 
v_x = t25;
if ({ let t26 = v_x.clone(); Py::Bool({ let t27 = Py::Int(1i64); if py::eq(&t26, &t27) { true } else { false } }) }).truthy() { out.push(Py::Int(1i64)); }} Py::list(out) }, None)?;
    let _ = f__need(py::not(&{ let t32 = { let t28 = v_singles.clone(); Py::Bool({ let t29 = Py::Int(0i64); if py::gt(&t28, &t29)? { true } else { false } }) }; if !t32.truthy() { t32 } else { { let t30 = v_singles.clone(); Py::Bool({ let t31 = v_s_rare.clone(); if py::eq(&t30, &t31) { true } else { false } }) } } }), Py::str("the only rare species are singletons, so ACE is undefined (use the bias-corrected Chao1)"))?;
    v_s_abun = py::sub(&py::len(&v_c.clone())?, &v_s_rare.clone())?;
    if ({ let t33 = v_s_rare.clone(); Py::Bool({ let t34 = Py::Int(0i64); if py::eq(&t33, &t34) { true } else { false } }) }).truthy() {
        return Ok(Py::dict(vec![(Py::str("ace"), py::float(&v_s_abun.clone())?)]));
    }
    v_n_rare = py::sum(&v_rare.clone(), None)?;
    v_c_ace = py::sub(&Py::Int(1i64), &py::fraction(&[v_singles.clone(), v_n_rare.clone()])?)?;
    v_top = py::mul(&v_s_rare.clone(), &py::sum(&{ let mut out = vec![]; for t35 in py::iter(&(v_rare.clone()))? { let mut v_x: Py; 
v_x = t35;
out.push(py::mul(&v_x.clone(), &py::sub(&v_x.clone(), &Py::Int(1i64))?)?); } Py::list(out) }, None)?)?;
    v_gamma = py::sub(&py::truediv(&py::fraction(&[v_top.clone()])?, &py::mul(&py::mul(&v_c_ace.clone(), &v_n_rare.clone())?, &py::sub(&v_n_rare.clone(), &Py::Int(1i64))?)?)?, &Py::Int(1i64))?;
    v_gamma = py::minmax(&[v_gamma.clone(), py::fraction(&[Py::Int(0i64)])?], true)?;
    if (v_no_gamma.clone()).truthy() {
        v_gamma = py::fraction(&[Py::Int(0i64)])?;
    }
    return Ok(Py::dict(vec![(Py::str("ace"), py::float(&py::add(&py::add(&v_s_abun.clone(), &py::truediv(&v_s_rare.clone(), &v_c_ace.clone())?)?, &py::mul(&py::truediv(&v_singles.clone(), &v_c_ace.clone())?, &v_gamma.clone())?)?)?)]));
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn f_aromaticity(mut v_seq: Py) -> R {
    let mut v_n: Py = Py::None;
    let mut v_s: Py = Py::None;
    v_s = f__seq(v_seq.clone(), Py::str("seq"), Py::str("PROTEIN"))?;
    v_n = py::add(&py::add(&py::method(&v_s.clone(), "count", &[Py::str("F")])?, &py::method(&v_s.clone(), "count", &[Py::str("W")])?)?, &py::method(&v_s.clone(), "count", &[Py::str("Y")])?)?;
    return Ok(Py::dict(vec![(Py::str("aromaticity"), py::truediv(&v_n.clone(), &py::len(&v_s.clone())?)?), (Py::str("aromatic_count"), v_n.clone())]));
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn f_beta_diversity(mut v_x: Py, mut v_y: Py, mut v_similarity: Py) -> R {
    let mut v_A: Py = Py::None;
    let mut v_B: Py = Py::None;
    let mut v_bc: Py = Py::None;
    let mut v_j: Py = Py::None;
    let mut v_s: Py = Py::None;
    let mut v_v: Py = Py::None;
    let _ = f__need({ let t229 = Py::Bool(matches!(v_x.clone(), Py::List(_))); if !t229.truthy() { t229 } else { { let t228 = Py::Bool(matches!(v_y.clone(), Py::List(_))); if !t228.truthy() { t228 } else { { let t227 = { let t223 = py::len(&v_x.clone())?; Py::Bool({ let t224 = py::len(&v_y.clone())?; if py::eq(&t223, &t224) { true } else { false } }) }; if !t227.truthy() { t227 } else { { let t225 = py::len(&v_x.clone())?; Py::Bool({ let t226 = Py::Int(1i64); if py::ge(&t225, &t226)? { true } else { false } }) } } } } } } }, Py::str("x and y must list the same species in the same order"))?;
    for t230 in py::iter(&(py::add(&v_x.clone(), &v_y.clone())?))? {
        v_v = t230;
        let _ = f__fin(Py::str("abundance"), Py::tuple(vec![v_v.clone()]))?;
        let _ = f__need({ let t234 = v_v.clone(); Py::Bool({ let t235 = Py::Int(0i64); if py::ge(&t234, &t235)? { true } else { false } }) }, Py::str("abundances must be zero or positive"))?;
    }
    let _ = f__need({ let t242 = Py::Bool(({ let mut out = vec![]; for t236 in py::iter(&(v_x.clone()))? { let mut v_v: Py; 
v_v = t236;
out.push({ let t237 = v_v.clone(); Py::Bool({ let t238 = Py::Int(0i64); if py::gt(&t237, &t238)? { true } else { false } }) }); } Py::list(out) }).items()?.iter().any(|x| x.truthy())); if !t242.truthy() { t242 } else { Py::Bool(({ let mut out = vec![]; for t239 in py::iter(&(v_y.clone()))? { let mut v_v: Py; 
v_v = t239;
out.push({ let t240 = v_v.clone(); Py::Bool({ let t241 = Py::Int(0i64); if py::gt(&t240, &t241)? { true } else { false } }) }); } Py::list(out) }).items()?.iter().any(|x| x.truthy())) } }, Py::str("each sample needs at least one species"))?;
    v_bc = py::truediv(&py::fsum(&{ let mut out = vec![]; for t243 in py::iter(&(py::zip(&[v_x.clone(), v_y.clone()])?))? { let mut v_a: Py; let mut v_b: Py; 
let t244 = (t243).items()?;
if t244.len() != 2 { return Err(PyErr::Value("unpack".into())); }
v_a = t244[0].clone();
v_b = t244[1].clone();
out.push(py::abs(&py::sub(&v_a.clone(), &v_b.clone())?)?); } Py::list(out) })?, &py::fsum(&{ let mut out = vec![]; for t245 in py::iter(&(py::zip(&[v_x.clone(), v_y.clone()])?))? { let mut v_a: Py; let mut v_b: Py; 
let t246 = (t245).items()?;
if t246.len() != 2 { return Err(PyErr::Value("unpack".into())); }
v_a = t246[0].clone();
v_b = t246[1].clone();
out.push(py::add(&v_a.clone(), &v_b.clone())?); } Py::list(out) })?)?;
    v_A = py::set(&{ let mut out = vec![]; for t247 in py::iter(&(py::enumerate(&v_x.clone(), 0)?))? { let mut v_i: Py; let mut v_v: Py; 
let t248 = (t247).items()?;
if t248.len() != 2 { return Err(PyErr::Value("unpack".into())); }
v_i = t248[0].clone();
v_v = t248[1].clone();
if ({ let t249 = v_v.clone(); Py::Bool({ let t250 = Py::Int(0i64); if py::gt(&t249, &t250)? { true } else { false } }) }).truthy() { out.push(v_i.clone()); }} Py::list(out) })?;
    v_B = py::set(&{ let mut out = vec![]; for t251 in py::iter(&(py::enumerate(&v_y.clone(), 0)?))? { let mut v_i: Py; let mut v_v: Py; 
let t252 = (t251).items()?;
if t252.len() != 2 { return Err(PyErr::Value("unpack".into())); }
v_i = t252[0].clone();
v_v = t252[1].clone();
if ({ let t253 = v_v.clone(); Py::Bool({ let t254 = Py::Int(0i64); if py::gt(&t253, &t254)? { true } else { false } }) }).truthy() { out.push(v_i.clone()); }} Py::list(out) })?;
    v_j = py::fraction(&[py::len(&py::bitxor(&v_A.clone(), &v_B.clone())?)?, py::len(&py::bitor(&v_A.clone(), &v_B.clone())?)?])?;
    v_s = py::fraction(&[py::len(&py::bitxor(&v_A.clone(), &v_B.clone())?)?, py::add(&py::len(&v_A.clone())?, &py::len(&v_B.clone())?)?])?;
    if (v_similarity.clone()).truthy() {
        let t255 = (Py::tuple(vec![py::sub(&Py::Int(1i64), &v_j.clone())?, py::sub(&Py::Int(1i64), &v_s.clone())?])).items()?;
        if t255.len() != 2 { return Err(PyErr::Value("unpack".into())); }
        v_j = t255[0].clone();
        v_s = t255[1].clone();
    }
    return Ok(Py::dict(vec![(Py::str("bray_curtis"), v_bc.clone()), (Py::str("jaccard"), py::float(&v_j.clone())?), (Py::str("sorensen"), py::float(&v_s.clone())?)]));
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn f_brillouin(mut v_counts: Py, mut v_log10: Py) -> R {
    let mut v_N: Py = Py::None;
    let mut v_c: Py = Py::None;
    let mut v_m: Py = Py::None;
    let mut v_x: Py = Py::None;
    v_c = f__abund(v_counts.clone(), Py::str("counts"))?;
    v_N = py::sum(&v_c.clone(), None)?;
    v_m = py::factorial(&v_N.clone())?;
    for t315 in py::iter(&(v_c.clone()))? {
        v_x = t315;
        v_m = py::floordiv(&v_m, &py::factorial(&v_x.clone())?)?;
    }
    return Ok(f__out(Py::dict(vec![(Py::str("H_B"), py::truediv(&(if (v_log10.clone()).truthy() { py::math1("log10", &v_m.clone())? } else { py::math1("log", &v_m.clone())? }), &v_N.clone())?)]))?);
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn f_chao1(mut v_counts: Py, mut v_form: Py, mut v_swap: Py) -> R {
    let mut v_c: Py = Py::None;
    let mut v_d: Py = Py::None;
    let mut v_est: Py = Py::None;
    let mut v_o: Py = Py::None;
    let mut v_s: Py = Py::None;
    v_c = f__abund(v_counts.clone(), Py::str("counts"))?;
    let _ = f__need({ let t389 = v_form.clone(); Py::Bool({ let t390 = Py::tuple(vec![Py::str("bias_corrected"), Py::str("classic")]); if py::contains(&t390, &t389)? { true } else { false } }) }, Py::str("form must be bias_corrected or classic"))?;
    if (v_swap.clone()).truthy() {
        v_form = (if ({ let t391 = v_form.clone(); Py::Bool({ let t392 = Py::str("bias_corrected"); if py::eq(&t391, &t392) { true } else { false } }) }).truthy() { Py::str("classic") } else { Py::str("bias_corrected") });
    }
    let t399 = (f__osd(v_c.clone())?).items()?;
    if t399.len() != 3 { return Err(PyErr::Value("unpack".into())); }
    v_o = t399[0].clone();
    v_s = t399[1].clone();
    v_d = t399[2].clone();
    if ({ let t403 = { let t400 = v_form.clone(); Py::Bool({ let t401 = Py::str("classic"); if py::eq(&t400, &t401) { true } else { false } }) }; if !t403.truthy() { t403 } else { { let t402 = v_s.clone(); if !t402.truthy() { t402 } else { v_d.clone() } } } }).truthy() {
        v_est = py::add(&v_o.clone(), &py::fraction(&[py::mul(&v_s.clone(), &v_s.clone())?, py::mul(&Py::Int(2i64), &v_d.clone())?])?)?;
    } else {
        v_est = py::add(&v_o.clone(), &py::fraction(&[py::mul(&v_s.clone(), &py::sub(&v_s.clone(), &Py::Int(1i64))?)?, py::mul(&Py::Int(2i64), &py::add(&v_d.clone(), &Py::Int(1i64))?)?])?)?;
    }
    return Ok(Py::dict(vec![(Py::str("chao1"), py::float(&v_est.clone())?), (Py::str("s_obs"), v_o.clone()), (Py::str("singletons"), v_s.clone()), (Py::str("doubletons"), v_d.clone())]));
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn f_chemostat(mut v_mu_max: Py, mut v_Ks: Py, mut v_Y: Py, mut v_D: Py, mut v_S0: Py, mut v_no_residual: Py) -> R {
    let mut v_Ss: Py = Py::None;
    let _ = f__pos(Py::str("mu_max, Ks, Y, D and S0"), Py::tuple(vec![v_mu_max.clone(), v_Ks.clone(), v_Y.clone(), v_D.clone(), v_S0.clone()]))?;
    let _ = f__need({ let t451 = v_D.clone(); Py::Bool({ let t452 = v_mu_max.clone(); if py::lt(&t451, &t452)? { true } else { false } }) }, Py::str("the dilution rate must be below mu_max (otherwise the culture washes out)"))?;
    v_Ss = py::truediv(&py::mul(&v_Ks.clone(), &v_D.clone())?, &py::sub(&v_mu_max.clone(), &v_D.clone())?)?;
    let _ = f__need({ let t453 = v_Ss.clone(); Py::Bool({ let t454 = v_S0.clone(); if py::lt(&t453, &t454)? { true } else { false } }) }, Py::str("washout: the steady-state substrate would exceed the feed concentration"))?;
    return Ok(f__out(Py::dict(vec![(Py::str("S"), v_Ss.clone()), (Py::str("X"), (if (v_no_residual.clone()).truthy() { py::mul(&v_Y.clone(), &v_S0.clone())? } else { py::mul(&v_Y.clone(), &py::sub(&v_S0.clone(), &v_Ss.clone())?)? }))]))?);
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn f_coalescent_times(mut v_n: Py, mut v_Ne: Py) -> R {
    let mut v_a: Py = Py::None;
    let _ = f__need({ let t513 = Py::Bool(matches!(v_n.clone(), Py::Int(_) | Py::Big(_) | Py::Bool(_))); if !t513.truthy() { t513 } else { { let t510 = Py::Int(2i64); Py::Bool({ let t511 = v_n.clone(); if py::le(&t510, &t511)? { { let t512 = Py::Int(100000i64); if py::le(&t511, &t512)? { true } else { false } } } else { false } }) } } }, Py::str("n must be a whole number of sampled gene copies (2 to 100000)"))?;
    let _ = f__pos(Py::str("Ne"), Py::tuple(vec![v_Ne.clone()]))?;
    v_a = py::fsum(&{ let mut out = vec![]; for t520 in py::range_iter(&[(Py::Int(1i64)).i()?, (v_n.clone()).i()?]) { let mut v_i: Py; 
v_i = t520;
out.push(py::truediv(&Py::Float(1.0), &v_i.clone())?); } Py::list(out) })?;
    return Ok(f__out(Py::dict(vec![(Py::str("t_mrca"), py::mul(&py::mul(&Py::Int(4i64), &v_Ne.clone())?, &py::sub(&Py::Int(1i64), &py::truediv(&Py::Int(1i64), &v_n.clone())?)?)?), (Py::str("total_length"), py::mul(&py::mul(&Py::Int(4i64), &v_Ne.clone())?, &v_a.clone())?)]))?);
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn f_competitive_fit(mut v_S: Py, mut v_I: Py, mut v_v: Py, mut v_model: Py) -> R {
    let mut v_J: Py = Py::None;
    let mut v_base: Py = Py::None;
    let mut v_f: Py = Py::None;
    let mut v_g: Py = Py::None;
    let mut v_mdl: Py = Py::None;
    let mut v_p: Py = Py::None;
    let mut v_r: Py = Py::None;
    let mut v_r0: Py = Py::None;
    let mut v_x: Py = Py::None;
    let _ = f__need({ let t582 = Py::Bool(matches!(v_S.clone(), Py::List(_))); if !t582.truthy() { t582 } else { { let t581 = Py::Bool(matches!(v_I.clone(), Py::List(_))); if !t581.truthy() { t581 } else { { let t580 = Py::Bool(matches!(v_v.clone(), Py::List(_))); if !t580.truthy() { t580 } else { { let t576 = py::len(&v_S.clone())?; Py::Bool({ let t577 = py::len(&v_I.clone())?; if py::eq(&t576, &t577) { { let t578 = py::len(&v_v.clone())?; if py::eq(&t577, &t578) { { let t579 = Py::Int(4i64); if py::ge(&t578, &t579)? { true } else { false } } } else { false } } } else { false } }) } } } } } } }, Py::str("S, I and v must be lists of equal length (at least 4 points)"))?;
    for t583 in py::iter(&(py::add(&v_S.clone(), &v_v.clone())?))? {
        v_x = t583;
        let _ = f__pos(Py::str("every substrate concentration and rate"), Py::tuple(vec![v_x.clone()]))?;
    }
    for t590 in py::iter(&(v_I.clone()))? {
        v_x = t590;
        let _ = f__fin(Py::str("I"), Py::tuple(vec![v_x.clone()]))?;
        let _ = f__need({ let t594 = v_x.clone(); Py::Bool({ let t595 = Py::Int(0i64); if py::ge(&t594, &t595)? { true } else { false } }) }, Py::str("inhibitor concentrations must be zero or positive"))?;
    }
    let _ = f__need({ let t602 = Py::Bool(({ let mut out = vec![]; for t596 in py::iter(&(v_I.clone()))? { let mut v_x: Py; 
v_x = t596;
out.push({ let t597 = v_x.clone(); Py::Bool({ let t598 = Py::Int(0i64); if py::eq(&t597, &t598) { true } else { false } }) }); } Py::list(out) }).items()?.iter().any(|x| x.truthy())); if !t602.truthy() { t602 } else { Py::Bool(({ let mut out = vec![]; for t599 in py::iter(&(v_I.clone()))? { let mut v_x: Py; 
v_x = t599;
out.push({ let t600 = v_x.clone(); Py::Bool({ let t601 = Py::Int(0i64); if py::gt(&t600, &t601)? { true } else { false } }) }); } Py::list(out) }).items()?.iter().any(|x| x.truthy())) } }, Py::str("include rates without inhibitor and with inhibitor"))?;
    let t606 = (Py::tuple(vec![{ let mut out = vec![]; for t603 in py::iter(&(v_S.clone()))? { let mut v_x: Py; 
v_x = t603;
out.push(py::float(&v_x.clone())?); } Py::list(out) }, { let mut out = vec![]; for t604 in py::iter(&(v_I.clone()))? { let mut v_x: Py; 
v_x = t604;
out.push(py::float(&v_x.clone())?); } Py::list(out) }, { let mut out = vec![]; for t605 in py::iter(&(v_v.clone()))? { let mut v_x: Py; 
v_x = t605;
out.push(py::float(&v_x.clone())?); } Py::list(out) }])).items()?;
    if t606.len() != 3 { return Err(PyErr::Value("unpack".into())); }
    v_S = t606[0].clone();
    v_I = t606[1].clone();
    v_v = t606[2].clone();
    v_base = { let mut out = vec![]; for t607 in py::iter(&(py::zip(&[v_S.clone(), v_I.clone(), v_v.clone()])?))? { let mut v_i: Py; let mut v_s: Py; let mut v_y: Py; 
let t608 = (t607).items()?;
if t608.len() != 3 { return Err(PyErr::Value("unpack".into())); }
v_s = t608[0].clone();
v_i = t608[1].clone();
v_y = t608[2].clone();
if ({ let t609 = v_i.clone(); Py::Bool({ let t610 = Py::Int(0i64); if py::eq(&t609, &t610) { true } else { false } }) }).truthy() { out.push(Py::tuple(vec![v_s.clone(), v_y.clone()])); }} Py::list(out) };
    let _ = f__need({ let t611 = py::len(&py::set(&{ let mut out = vec![]; for t613 in py::iter(&(v_base.clone()))? { let mut v__: Py; let mut v_s: Py; 
let t614 = (t613).items()?;
if t614.len() != 2 { return Err(PyErr::Value("unpack".into())); }
v_s = t614[0].clone();
v__ = t614[1].clone();
out.push(v_s.clone()); } Py::list(out) })?)?; Py::Bool({ let t612 = Py::Int(3i64); if py::ge(&t611, &t612)? { true } else { false } }) }, Py::str("at least 3 different S without inhibitor are needed"))?;
    v_r0 = f_mm_fit({ let mut out = vec![]; for t685 in py::iter(&(v_base.clone()))? { let mut v__: Py; let mut v_s: Py; 
let t686 = (t685).items()?;
if t686.len() != 2 { return Err(PyErr::Value("unpack".into())); }
v_s = t686[0].clone();
v__ = t686[1].clone();
out.push(v_s.clone()); } Py::list(out) }, { let mut out = vec![]; for t687 in py::iter(&(v_base.clone()))? { let mut v__: Py; let mut v_y: Py; 
let t688 = (t687).items()?;
if t688.len() != 2 { return Err(PyErr::Value("unpack".into())); }
v__ = t688[0].clone();
v_y = t688[1].clone();
out.push(v_y.clone()); } Py::list(out) })?;
    v_mdl = (if ({ let t689 = v_model.clone(); Py::Bool({ let t690 = Py::str("competitive"); if py::eq(&t689, &t690) { true } else { false } }) }).truthy() { f__ci_model(v_S.clone(), v_I.clone())? } else { f__nc_model(v_S.clone(), v_I.clone())? });
    v_p = f__lm(v_mdl.clone(), Py::list(vec![py::getitem(&v_r0.clone(), &Py::str("Vmax"))?, py::getitem(&v_r0.clone(), &Py::str("Km"))?, py::minmax(&[v_I.clone()], true)?]), v_v.clone(), { let v_I = v_I.clone(); let v_S = v_S.clone(); let v_base = v_base.clone(); let v_mdl = v_mdl.clone(); let v_model = v_model.clone(); let v_r0 = v_r0.clone(); let v_v = v_v.clone(); let v_x = v_x.clone(); py::func(move |a: &[Py]| -> R {
            #[allow(unused_mut, unused_variables, unused_assignments)]
{
            let mut v_q: Py = a.get(0).cloned().unwrap_or(Py::None);
            return Ok({ let t764 = { let t757 = py::getitem(&v_q.clone(), &Py::Int(0i64))?; Py::Bool({ let t758 = Py::Int(0i64); if py::gt(&t757, &t758)? { true } else { false } }) }; if !t764.truthy() { t764 } else { { let t763 = { let t759 = py::getitem(&v_q.clone(), &Py::Int(1i64))?; Py::Bool({ let t760 = Py::Int(0i64); if py::gt(&t759, &t760)? { true } else { false } }) }; if !t763.truthy() { t763 } else { { let t761 = py::getitem(&v_q.clone(), &Py::Int(2i64))?; Py::Bool({ let t762 = Py::Int(0i64); if py::gt(&t761, &t762)? { true } else { false } }) } } } } });
            #[allow(unreachable_code)]
            Ok(Py::None)
        }
        }) })?;
    let t765 = (py::call(&v_mdl.clone(), &[v_p.clone()])?).items()?;
    if t765.len() != 2 { return Err(PyErr::Value("unpack".into())); }
    v_f = t765[0].clone();
    v_J = t765[1].clone();
    v_r = { let mut out = vec![]; for t766 in py::iter(&(py::zip(&[v_v.clone(), v_f.clone()])?))? { let mut v_a: Py; let mut v_y: Py; 
let t767 = (t766).items()?;
if t767.len() != 2 { return Err(PyErr::Value("unpack".into())); }
v_y = t767[0].clone();
v_a = t767[1].clone();
out.push(py::sub(&v_y.clone(), &v_a.clone())?); } Py::list(out) };
    v_g = { let mut out = vec![]; for t768 in py::range_iter(&[(Py::Int(3i64)).i()?]) { let mut v_i: Py; 
v_i = t768;
out.push(py::abs(&py::fsum(&{ let mut out = vec![]; for t769 in py::iter(&(py::zip(&[v_J.clone(), v_r.clone()])?))? { let mut v_ri: Py; let mut v_row: Py; 
let t770 = (t769).items()?;
if t770.len() != 2 { return Err(PyErr::Value("unpack".into())); }
v_row = t770[0].clone();
v_ri = t770[1].clone();
out.push(py::mul(&py::getitem(&v_row.clone(), &v_i.clone())?, &v_ri.clone())?); } Py::list(out) })?)?); } Py::list(out) };
    let _ = f__need({ let t771 = py::minmax(&[{ let mut out = vec![]; for t773 in py::iter(&(py::zip(&[v_g.clone(), v_p.clone()])?))? { let mut v_gi: Py; let mut v_pi: Py; 
let t774 = (t773).items()?;
if t774.len() != 2 { return Err(PyErr::Value("unpack".into())); }
v_gi = t774[0].clone();
v_pi = t774[1].clone();
out.push(py::mul(&v_gi.clone(), &py::abs(&v_pi.clone())?)?); } Py::list(out) }], true)?; Py::Bool({ let t772 = py::mul(&Py::Float(1e-09), &py::fsum(&{ let mut out = vec![]; for t775 in py::iter(&(v_v.clone()))? { let mut v_y: Py; 
v_y = t775;
out.push(py::mul(&v_y.clone(), &v_y.clone())?); } Py::list(out) })?)?; if py::le(&t771, &t772)? { true } else { false } }) }, Py::str("the inhibition fit did not converge"))?;
    return Ok(f__out(Py::dict(vec![(Py::str("Vmax"), py::getitem(&v_p.clone(), &Py::Int(0i64))?), (Py::str("Km"), py::getitem(&v_p.clone(), &Py::Int(1i64))?), (Py::str("Ki"), py::getitem(&v_p.clone(), &Py::Int(2i64))?)]))?);
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn f_cross_genotype(mut v_parent1: Py, mut v_parent2: Py, mut v_offspring: Py) -> R {
    let mut v_a: Py = Py::None;
    let mut v_b: Py = Py::None;
    let mut v_hits: Py = Py::None;
    let mut v_o: Py = Py::None;
    let mut v_p: Py = Py::None;
    let mut v_want: Py = Py::None;
    let mut v_x1: Py = Py::None;
    let mut v_x2: Py = Py::None;
    let mut v_y1: Py = Py::None;
    let mut v_y2: Py = Py::None;
    let t833 = (Py::tuple(vec![f__genotype(v_parent1.clone(), Py::str("parent1"))?, f__genotype(v_parent2.clone(), Py::str("parent2"))?, f__genotype(v_offspring.clone(), Py::str("offspring"))?])).items()?;
    if t833.len() != 3 { return Err(PyErr::Value("unpack".into())); }
    v_a = t833[0].clone();
    v_b = t833[1].clone();
    v_o = t833[2].clone();
    let _ = f__same_loci(Py::tuple(vec![v_a.clone(), v_b.clone(), v_o.clone()]))?;
    v_p = py::fraction(&[Py::Int(1i64)])?;
    for t839 in py::iter(&(py::zip(&[v_a.clone(), v_b.clone(), v_o.clone()])?))? {
        let t844 = (t839).items()?;
        if t844.len() != 3 { return Err(PyErr::Value("unpack".into())); }
        let t845 = (t844[0].clone()).items()?;
        if t845.len() != 2 { return Err(PyErr::Value("unpack".into())); }
        v_x1 = t845[0].clone();
        v_y1 = t845[1].clone();
        let t846 = (t844[1].clone()).items()?;
        if t846.len() != 2 { return Err(PyErr::Value("unpack".into())); }
        v_x2 = t846[0].clone();
        v_y2 = t846[1].clone();
        v_want = t844[2].clone();
        v_hits = py::sum(&{ let mut out = vec![]; for t840 in py::iter(&(Py::tuple(vec![v_x1.clone(), v_y1.clone()])))? { let mut v_u: Py; 
v_u = t840;
for t841 in py::iter(&(Py::tuple(vec![v_x2.clone(), v_y2.clone()])))? { let mut v_v: Py; 
v_v = t841;
if ({ let t842 = Py::tuple((py::sorted_key(&Py::tuple(vec![v_u.clone(), v_v.clone()]), &{ let v_a = v_a.clone(); let v_b = v_b.clone(); let v_o = v_o.clone(); let v_offspring = v_offspring.clone(); let v_p = v_p.clone(); let v_parent1 = v_parent1.clone(); let v_parent2 = v_parent2.clone(); let v_u = v_u.clone(); let v_v = v_v.clone(); let v_want = v_want.clone(); let v_x1 = v_x1.clone(); let v_x2 = v_x2.clone(); let v_y1 = v_y1.clone(); let v_y2 = v_y2.clone(); py::func(move |a: &[Py]| -> R {
            #[allow(unused_mut, unused_variables, unused_assignments)]
{
            let mut v_c: Py = a.get(0).cloned().unwrap_or(Py::None);
            return Ok(Py::tuple(vec![py::method(&v_c.clone(), "islower", &[])?, v_c.clone()]));
            #[allow(unreachable_code)]
            Ok(Py::None)
        }
        }) }, (Py::Bool(false)).truthy())?).items()?); Py::Bool({ let t843 = v_want.clone(); if py::eq(&t842, &t843) { true } else { false } }) }).truthy() { out.push(Py::Int(1i64)); }}} Py::list(out) }, None)?;
        v_p = py::mul(&v_p, &py::fraction(&[v_hits.clone(), Py::Int(4i64)])?)?;
    }
    return Ok(Py::dict(vec![(Py::str("probability"), py::float(&v_p.clone())?), (Py::str("fraction"), f__fr_text(v_p.clone())?)]));
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn f_cross_phenotypes(mut v_parent1: Py, mut v_parent2: Py, mut v_complete: Py) -> R {
    let mut v__c: Py = Py::None;
    let mut v__dom: Py = Py::None;
    let mut v_a: Py = Py::None;
    let mut v_b: Py = Py::None;
    let mut v_classes: Py = Py::None;
    let mut v_d: Py = Py::None;
    let mut v_dom: Py = Py::None;
    let mut v_n_geno: Py = Py::None;
    let mut v_p_dom: Py = Py::None;
    let mut v_p_rec: Py = Py::None;
    let mut v_per: Py = Py::None;
    let mut v_probs: Py = Py::None;
    let mut v_x1: Py = Py::None;
    let mut v_x2: Py = Py::None;
    let mut v_y1: Py = Py::None;
    let mut v_y2: Py = Py::None;
    let t904 = (Py::tuple(vec![f__genotype(v_parent1.clone(), Py::str("parent1"))?, f__genotype(v_parent2.clone(), Py::str("parent2"))?])).items()?;
    if t904.len() != 2 { return Err(PyErr::Value("unpack".into())); }
    v_a = t904[0].clone();
    v_b = t904[1].clone();
    let _ = f__same_loci(Py::tuple(vec![v_a.clone(), v_b.clone()]))?;
    v_per = Py::list(vec![]);
    v_n_geno = Py::Int(1i64);
    for t910 in py::iter(&(py::zip(&[v_a.clone(), v_b.clone()])?))? {
        let t918 = (t910).items()?;
        if t918.len() != 2 { return Err(PyErr::Value("unpack".into())); }
        let t919 = (t918[0].clone()).items()?;
        if t919.len() != 2 { return Err(PyErr::Value("unpack".into())); }
        v_x1 = t919[0].clone();
        v_y1 = t919[1].clone();
        let t920 = (t918[1].clone()).items()?;
        if t920.len() != 2 { return Err(PyErr::Value("unpack".into())); }
        v_x2 = t920[0].clone();
        v_y2 = t920[1].clone();
        v_d = f__locus_dist(Py::tuple(vec![v_x1.clone(), v_y1.clone()]), Py::tuple(vec![v_x2.clone(), v_y2.clone()]))?;
        v_n_geno = py::mul(&v_n_geno, &py::len(&v_d.clone())?)?;
        if (v_complete.clone()).truthy() {
            v_dom = py::sum(&{ let mut out = vec![]; for t913 in py::iter(&(py::method(&v_d.clone(), "items", &[])?))? { let mut v_k: Py; let mut v_p: Py; 
let t914 = (t913).items()?;
if t914.len() != 2 { return Err(PyErr::Value("unpack".into())); }
v_k = t914[0].clone();
v_p = t914[1].clone();
if (py::method(&py::getitem(&v_k.clone(), &Py::Int(0i64))?, "isupper", &[])?).truthy() { out.push(v_p.clone()); }} Py::list(out) }, None)?;
            v_classes = { let mut out = vec![]; for t915 in py::iter(&(Py::tuple(vec![v_dom.clone(), py::sub(&Py::Int(1i64), &v_dom.clone())?])))? { let mut v_p: Py; 
v_p = t915;
if ({ let t916 = v_p.clone(); Py::Bool({ let t917 = Py::Int(0i64); if py::gt(&t916, &t917)? { true } else { false } }) }).truthy() { out.push(v_p.clone()); }} Py::list(out) };
            let _ = py::method(&v_per.clone(), "append", &[Py::tuple(vec![v_dom.clone(), v_classes.clone()])])?;
        } else {
            let _ = py::method(&v_per.clone(), "append", &[Py::tuple(vec![py::method(&v_d.clone(), "get", &[Py::tuple(vec![py::method(&v_x1.clone(), "upper", &[])?, py::method(&v_x1.clone(), "upper", &[])?]), py::fraction(&[Py::Int(0i64)])?])?, Py::list((py::method(&v_d.clone(), "values", &[])?).items()?)])])?;
        }
    }
    v_p_dom = py::fraction(&[Py::Int(1i64)])?;
    for t921 in py::iter(&(v_per.clone()))? {
        let t922 = (t921).items()?;
        if t922.len() != 2 { return Err(PyErr::Value("unpack".into())); }
        v_dom = t922[0].clone();
        v__c = t922[1].clone();
        v_p_dom = py::mul(&v_p_dom, &v_dom.clone())?;
    }
    v_probs = Py::list(vec![py::fraction(&[Py::Int(1i64)])?]);
    for t923 in py::iter(&(v_per.clone()))? {
        let t926 = (t923).items()?;
        if t926.len() != 2 { return Err(PyErr::Value("unpack".into())); }
        v__dom = t926[0].clone();
        v_classes = t926[1].clone();
        v_probs = { let mut out = vec![]; for t924 in py::iter(&(v_probs.clone()))? { let mut v_p: Py; 
v_p = t924;
for t925 in py::iter(&(v_classes.clone()))? { let mut v_c: Py; 
v_c = t925;
out.push(py::mul(&v_p.clone(), &v_c.clone())?); }} Py::list(out) };
    }
    v_p_rec = py::fraction(&[Py::Int(1i64)])?;
    for t927 in py::iter(&(py::zip(&[v_a.clone(), v_b.clone()])?))? {
        let t928 = (t927).items()?;
        if t928.len() != 2 { return Err(PyErr::Value("unpack".into())); }
        let t929 = (t928[0].clone()).items()?;
        if t929.len() != 2 { return Err(PyErr::Value("unpack".into())); }
        v_x1 = t929[0].clone();
        v_y1 = t929[1].clone();
        let t930 = (t928[1].clone()).items()?;
        if t930.len() != 2 { return Err(PyErr::Value("unpack".into())); }
        v_x2 = t930[0].clone();
        v_y2 = t930[1].clone();
        v_p_rec = py::mul(&v_p_rec, &py::method(&f__locus_dist(Py::tuple(vec![v_x1.clone(), v_y1.clone()]), Py::tuple(vec![v_x2.clone(), v_y2.clone()]))?, "get", &[Py::tuple(vec![py::method(&v_x1.clone(), "lower", &[])?, py::method(&v_x1.clone(), "lower", &[])?]), py::fraction(&[Py::Int(0i64)])?])?)?;
    }
    return Ok(Py::dict(vec![(Py::str("p_all_dominant"), py::float(&v_p_dom.clone())?), (Py::str("p_all_recessive"), py::float(&v_p_rec.clone())?), (Py::str("n_genotypes"), v_n_geno.clone()), (Py::str("n_phenotypes"), py::len(&v_probs.clone())?), (Py::str("ratio"), f__ratio_text(v_probs.clone())?)]));
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn f_dominance_indices(mut v_counts: Py, mut v_wrong: Py) -> R {
    let mut v_N: Py = Py::None;
    let mut v_U: Py = Py::None;
    let mut v_c: Py = Py::None;
    v_c = f__abund(v_counts.clone(), Py::str("counts"))?;
    v_N = py::sum(&v_c.clone(), None)?;
    let _ = f__need({ let t997 = v_N.clone(); Py::Bool({ let t998 = Py::Int(2i64); if py::ge(&t997, &t998)? { true } else { false } }) }, Py::str("at least 2 individuals are needed"))?;
    v_U = py::math1("sqrt", &py::sum(&{ let mut out = vec![]; for t999 in py::iter(&(v_c.clone()))? { let mut v_x: Py; 
v_x = t999;
out.push(py::mul(&v_x.clone(), &v_x.clone())?); } Py::list(out) }, None)?)?;
    return Ok(f__out(Py::dict(vec![(Py::str("berger_parker"), py::truediv(&py::minmax(&[v_c.clone()], true)?, &v_N.clone())?), (Py::str("mcintosh"), py::truediv(&py::sub(&v_N.clone(), &v_U.clone())?, &(if (v_wrong.clone()).truthy() { v_N.clone() } else { py::sub(&v_N.clone(), &py::math1("sqrt", &v_N.clone())?)? }))?)]))?);
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn f_drift_heterozygosity(mut v_N: Py, mut v_i: Py, mut v_t: Py) -> R {
    let mut v_h0: Py = Py::None;
    let mut v_p: Py = Py::None;
    let mut v_ret: Py = Py::None;
    let _ = f__need({ let t1058 = Py::Bool(matches!(v_N.clone(), Py::Int(_) | Py::Big(_) | Py::Bool(_))); if !t1058.truthy() { t1058 } else { { let t1055 = Py::Int(1i64); Py::Bool({ let t1056 = v_N.clone(); if py::le(&t1055, &t1056)? { { let t1057 = Py::Int(100000i64); if py::le(&t1056, &t1057)? { true } else { false } } } else { false } }) } } }, Py::str("N must be a whole number of diploid individuals (1 to 100000)"))?;
    let _ = f__need({ let t1062 = Py::Bool(matches!(v_i.clone(), Py::Int(_) | Py::Big(_) | Py::Bool(_))); if !t1062.truthy() { t1062 } else { { let t1059 = Py::Int(0i64); Py::Bool({ let t1060 = v_i.clone(); if py::lt(&t1059, &t1060)? { { let t1061 = py::mul(&Py::Int(2i64), &v_N.clone())?; if py::lt(&t1060, &t1061)? { true } else { false } } } else { false } }) } } }, Py::str("i must be a whole number of copies between 1 and 2N - 1"))?;
    let _ = f__need({ let t1066 = Py::Bool(matches!(v_t.clone(), Py::Int(_) | Py::Big(_) | Py::Bool(_))); if !t1066.truthy() { t1066 } else { { let t1063 = Py::Int(0i64); Py::Bool({ let t1064 = v_t.clone(); if py::le(&t1063, &t1064)? { { let t1065 = Py::Int(100000i64); if py::le(&t1064, &t1065)? { true } else { false } } } else { false } }) } } }, Py::str("t must be a whole number of generations"))?;
    v_p = py::fraction(&[v_i.clone(), py::mul(&Py::Int(2i64), &v_N.clone())?])?;
    v_h0 = py::float(&py::mul(&py::mul(&Py::Int(2i64), &v_p.clone())?, &py::sub(&Py::Int(1i64), &v_p.clone())?)?)?;
    v_ret = py::pow(&py::sub(&Py::Int(1i64), &py::truediv(&Py::Int(1i64), &py::mul(&Py::Int(2i64), &v_N.clone())?)?)?, &v_t.clone())?;
    return Ok(f__out(Py::dict(vec![(Py::str("H0"), v_h0.clone()), (Py::str("H_t"), py::mul(&v_h0.clone(), &v_ret.clone())?), (Py::str("retained"), v_ret.clone())]))?);
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn f_eadie_hofstee(mut v_S: Py, mut v_v: Py) -> R {
    let mut v_a: Py = Py::None;
    let mut v_b: Py = Py::None;
    let mut v_r2: Py = Py::None;
    let t1138 = (f__kin_data(v_S.clone(), v_v.clone(), Py::Int(2i64))?).items()?;
    if t1138.len() != 2 { return Err(PyErr::Value("unpack".into())); }
    v_S = t1138[0].clone();
    v_v = t1138[1].clone();
    let t1150 = (f__ols({ let mut out = vec![]; for t1148 in py::iter(&(py::zip(&[v_S.clone(), v_v.clone()])?))? { let mut v_x: Py; let mut v_y: Py; 
let t1149 = (t1148).items()?;
if t1149.len() != 2 { return Err(PyErr::Value("unpack".into())); }
v_x = t1149[0].clone();
v_y = t1149[1].clone();
out.push(py::truediv(&v_y.clone(), &v_x.clone())?); } Py::list(out) }, v_v.clone())?).items()?;
    if t1150.len() != 3 { return Err(PyErr::Value("unpack".into())); }
    v_b = t1150[0].clone();
    v_a = t1150[1].clone();
    v_r2 = t1150[2].clone();
    let _ = f__need({ let t1155 = { let t1151 = v_b.clone(); Py::Bool({ let t1152 = Py::Int(0i64); if py::lt(&t1151, &t1152)? { true } else { false } }) }; if !t1155.truthy() { t1155 } else { { let t1153 = v_a.clone(); Py::Bool({ let t1154 = Py::Int(0i64); if py::gt(&t1153, &t1154)? { true } else { false } }) } } }, Py::str("the Eadie-Hofstee line must have a negative slope and a positive intercept"))?;
    return Ok(f__out(Py::dict(vec![(Py::str("Vmax"), v_a.clone()), (Py::str("Km"), py::neg(&v_b.clone())?), (Py::str("r2"), v_r2.clone())]))?);
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn f_edit_distance(mut v_seq1: Py, mut v_seq2: Py) -> R {
    let mut v_a: Py = Py::None;
    let mut v_b: Py = Py::None;
    let t1234 = (Py::tuple(vec![f__seq(v_seq1.clone(), Py::str("seq1"), Py::str("ANY"))?, f__seq(v_seq2.clone(), Py::str("seq2"), Py::str("ANY"))?])).items()?;
    if t1234.len() != 2 { return Err(PyErr::Value("unpack".into())); }
    v_a = t1234[0].clone();
    v_b = t1234[1].clone();
    let _ = f__need({ let t1235 = py::mul(&py::len(&v_a.clone())?, &py::len(&v_b.clone())?)?; Py::Bool({ let t1236 = Py::Int(4000000i64); if py::le(&t1235, &t1236)? { true } else { false } }) }, Py::str("sequences too long (length product at most 4 000 000)"))?;
    return Ok(Py::dict(vec![(Py::str("distance"), f__levenshtein(v_a.clone(), v_b.clone())?)]));
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn f_extinction(mut v_seq: Py) -> R {
    let mut v_m: Py = Py::None;
    let mut v_ox: Py = Py::None;
    let mut v_red: Py = Py::None;
    let mut v_s: Py = Py::None;
    v_s = f__seq(v_seq.clone(), Py::str("seq"), Py::str("PROTEIN"))?;
    v_red = py::add(&py::mul(&Py::Int(5500i64), &py::method(&v_s.clone(), "count", &[Py::str("W")])?)?, &py::mul(&Py::Int(1490i64), &py::method(&v_s.clone(), "count", &[Py::str("Y")])?)?)?;
    v_ox = py::add(&v_red.clone(), &py::mul(&py::floordiv(&py::method(&v_s.clone(), "count", &[Py::str("C")])?, &Py::Int(2i64))?, &Py::Int(125i64))?)?;
    v_m = f__chain_mass(v_s.clone(), c_PROTEIN_W()?)?;
    return Ok(f__out(Py::dict(vec![(Py::str("eps_reduced"), v_red.clone()), (Py::str("eps_cystines"), v_ox.clone()), (Py::str("abs_reduced"), py::truediv(&v_red.clone(), &v_m.clone())?), (Py::str("abs_cystines"), py::truediv(&v_ox.clone(), &v_m.clone())?)]))?);
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn f_faith_pd(mut v_tree: Py, mut v_taxa: Py, mut v_unrooted: Py) -> R {
    let mut v_keep: Py = Py::None;
    let mut v_kids: Py = Py::None;
    let mut v_leaves: Py = Py::None;
    let mut v_names: Py = Py::None;
    let mut v_node: Py = Py::None;
    let mut v_root: Py = Py::None;
    let mut v_walk: Py = Py::None;
    let mut v_want: Py = Py::None;
    let mut v_x: Py = Py::None;
    let c_total = py::cell_new();
    let t1326 = (f__parse_newick(v_tree.clone())?).items()?;
    if t1326.len() != 2 { return Err(PyErr::Value("unpack".into())); }
    v_root = t1326[0].clone();
    v_leaves = t1326[1].clone();
    v_names = (if (Py::Bool(matches!(v_taxa.clone(), Py::Str(_)))).truthy() { py::method(&v_taxa.clone(), "split", &[])? } else { Py::list(vec![]) });
    let _ = f__need({ let t1329 = v_names.clone(); if !t1329.truthy() { t1329 } else { { let t1327 = py::len(&py::set(&v_names.clone())?)?; Py::Bool({ let t1328 = py::len(&v_names.clone())?; if py::eq(&t1327, &t1328) { true } else { false } }) } } }, Py::str("give the observed taxa as distinct names separated by spaces"))?;
    for t1330 in py::iter(&(v_names.clone()))? {
        v_x = t1330;
        let _ = f__need({ let t1331 = v_x.clone(); Py::Bool({ let t1332 = v_leaves.clone(); if py::contains(&t1332, &t1331)? { true } else { false } }) }, Py::str(&[py::format(&v_x.clone(), "")?, " is not a leaf of the tree".to_string()].concat()))?;
    }
    v_want = py::set(&v_names.clone())?;
    py::cell_set(&c_total, py::fraction(&[Py::Int(0i64)])?);
    v_keep = Py::list(vec![]);
    let rec1333 = py::rec_cell();
    v_walk = { let v_keep = v_keep.clone(); let v_leaves = v_leaves.clone(); let v_names = v_names.clone(); let v_root = v_root.clone(); let v_taxa = v_taxa.clone(); let v_tree = v_tree.clone(); let v_unrooted = v_unrooted.clone(); let v_want = v_want.clone(); let v_x = v_x.clone(); let c_total = c_total.clone(); let rec1333 = rec1333.clone(); py::func(move |a: &[Py]| -> R {
            #[allow(unused_mut, unused_variables, unused_assignments)]
{
            let mut v_t: Py = a.get(0).cloned().unwrap_or(Py::None);
            let mut v_is_root: Py = a.get(1).cloned().unwrap_or(Py::None);
            let v_walk: Py = py::rec_get(&rec1333)?;
            let mut v_hit: Py = Py::None;
            let mut v_hits: Py = Py::None;
            if (py::not(&py::getitem(&v_t.clone(), &Py::Int(2i64))?)).truthy() {
                v_hit = { let t1334 = py::getitem(&v_t.clone(), &Py::Int(0i64))?; Py::Bool({ let t1335 = v_want.clone(); if py::contains(&t1335, &t1334)? { true } else { false } }) };
            } else {
                v_hits = { let mut out = vec![]; for t1336 in py::iter(&(py::getitem(&v_t.clone(), &Py::Int(2i64))?))? { let mut v_k: Py; 
v_k = t1336;
out.push(py::call(&v_walk.clone(), &[v_k.clone(), Py::Bool(false)])?); } Py::list(out) };
                v_hit = Py::Bool((v_hits.clone()).items()?.iter().any(|x| x.truthy()));
                let _ = py::method(&v_keep.clone(), "append", &[Py::tuple(vec![v_t.clone(), py::sum(&v_hits.clone(), None)?])])?;
            }
            if ({ let t1337 = v_hit.clone(); if !t1337.truthy() { t1337 } else { py::not(&v_is_root.clone()) } }).truthy() {
                let _ = f__need({ let t1338 = py::getitem(&v_t.clone(), &Py::Int(1i64))?; Py::Bool({ let t1339 = Py::None; if !py::is(&t1338, &t1339) { true } else { false } }) }, Py::str("every branch needs a length"))?;
                py::cell_set(&c_total, py::iadd(&py::cell_get(&c_total)?, &py::fraction(&[Py::str(&py::repr(&py::getitem(&v_t.clone(), &Py::Int(1i64))?))])?)?);
            }
            return Ok(v_hit.clone());
            #[allow(unreachable_code)]
            Ok(Py::None)
        }
        }) };
    py::rec_set(&rec1333, &v_walk);
    let _ = py::call(&v_walk.clone(), &[v_root.clone(), Py::Bool(true)])?;
    if (v_unrooted.clone()).truthy() {
        v_node = v_root.clone();
        while (py::getitem(&v_node.clone(), &Py::Int(2i64))?).truthy() {
            v_kids = { let mut out = vec![]; for t1340 in py::iter(&(py::getitem(&v_node.clone(), &Py::Int(2i64))?))? { let mut v_k: Py; 
v_k = t1340;
if (f__has_any(v_k.clone(), v_want.clone())?).truthy() { out.push(v_k.clone()); }} Py::list(out) };
            if ({ let t1344 = py::len(&v_kids.clone())?; Py::Bool({ let t1345 = Py::Int(1i64); if py::ne(&t1344, &t1345) { true } else { false } }) }).truthy() {
                break;
            }
            v_node = py::getitem(&v_kids.clone(), &Py::Int(0i64))?;
            if ({ let t1348 = py::getitem(&v_node.clone(), &Py::Int(2i64))?; if t1348.truthy() { t1348 } else { { let t1346 = py::getitem(&v_node.clone(), &Py::Int(0i64))?; Py::Bool({ let t1347 = v_want.clone(); if !py::contains(&t1347, &t1346)? { true } else { false } }) } } }).truthy() {
                py::cell_set(&c_total, py::sub(&py::cell_get(&c_total)?, &py::fraction(&[Py::str(&py::repr(&py::getitem(&v_node.clone(), &Py::Int(1i64))?))])?)?);
            } else {
                break;
            }
        }
    }
    return Ok(Py::dict(vec![(Py::str("pd"), py::float(&py::cell_get(&c_total)?)?)]));
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn f_fba(mut v_S: Py, mut v_lb: Py, mut v_ub: Py, mut v_c: Py, mut v_box_only: Py) -> R {
    let mut v_A: Py = Py::None;
    let mut v_F: Py = Py::None;
    let mut v_L: Py = Py::None;
    let mut v_Sf: Py = Py::None;
    let mut v_U: Py = Py::None;
    let mut v_a: Py = Py::None;
    let mut v_bnd: Py = Py::None;
    let mut v_bvec: Py = Py::None;
    let mut v_cvec: Py = Py::None;
    let mut v_j: Py = Py::None;
    let mut v_m: Py = Py::None;
    let mut v_n: Py = Py::None;
    let mut v_name: Py = Py::None;
    let mut v_row: Py = Py::None;
    let mut v_v: Py = Py::None;
    let mut v_w: Py = Py::None;
    let mut v_x: Py = Py::None;
    let mut v_z: Py = Py::None;
    let _ = f__need({ let t1507 = Py::Bool(matches!(v_S.clone(), Py::List(_))); if !t1507.truthy() { t1507 } else { { let t1506 = v_S.clone(); if !t1506.truthy() { t1506 } else { Py::Bool(({ let mut out = vec![]; for t1503 in py::iter(&(v_S.clone()))? { let mut v_r: Py; 
v_r = t1503;
out.push({ let t1504 = py::len(&v_r.clone())?; Py::Bool({ let t1505 = py::len(&py::getitem(&v_S.clone(), &Py::Int(0i64))?)?; if py::eq(&t1504, &t1505) { true } else { false } }) }); } Py::list(out) }).items()?.iter().all(|x| x.truthy())) } } } }, Py::str("S must be a stoichiometric matrix"))?;
    let t1508 = (Py::tuple(vec![py::len(&v_S.clone())?, py::len(&py::getitem(&v_S.clone(), &Py::Int(0i64))?)?])).items()?;
    if t1508.len() != 2 { return Err(PyErr::Value("unpack".into())); }
    v_m = t1508[0].clone();
    v_n = t1508[1].clone();
    let _ = f__need({ let t1513 = { let t1509 = v_m.clone(); Py::Bool({ let t1510 = Py::Int(15i64); if py::le(&t1509, &t1510)? { true } else { false } }) }; if !t1513.truthy() { t1513 } else { { let t1511 = v_n.clone(); Py::Bool({ let t1512 = Py::Int(20i64); if py::le(&t1511, &t1512)? { true } else { false } }) } } }, Py::str("at most 15 metabolites and 20 reactions"))?;
    for t1514 in py::iter(&(Py::tuple(vec![Py::tuple(vec![Py::str("lb"), v_lb.clone()]), Py::tuple(vec![Py::str("ub"), v_ub.clone()]), Py::tuple(vec![Py::str("c"), v_c.clone()])])))? {
        let t1518 = (t1514).items()?;
        if t1518.len() != 2 { return Err(PyErr::Value("unpack".into())); }
        v_name = t1518[0].clone();
        v_v = t1518[1].clone();
        let _ = f__need({ let t1517 = Py::Bool(matches!(v_v.clone(), Py::List(_))); if !t1517.truthy() { t1517 } else { { let t1515 = py::len(&v_v.clone())?; Py::Bool({ let t1516 = v_n.clone(); if py::eq(&t1515, &t1516) { true } else { false } }) } } }, Py::str(&[py::format(&v_name.clone(), "")?, " needs one value per reaction (".to_string(), py::format(&v_n.clone(), "")?, ")".to_string()].concat()))?;
    }
    for t1519 in py::iter(&(v_S.clone()))? {
        v_row = t1519;
        for t1520 in py::iter(&(v_row.clone()))? {
            v_x = t1520;
            let _ = f__fin(Py::str("S"), Py::tuple(vec![v_x.clone()]))?;
        }
    }
    for t1524 in py::iter(&(py::zip(&[v_lb.clone(), v_ub.clone(), v_c.clone()])?))? {
        let t1527 = (t1524).items()?;
        if t1527.len() != 3 { return Err(PyErr::Value("unpack".into())); }
        v_a = t1527[0].clone();
        v_bnd = t1527[1].clone();
        v_w = t1527[2].clone();
        let _ = f__fin(Py::str("bounds and objective"), Py::tuple(vec![v_a.clone(), v_bnd.clone(), v_w.clone()]))?;
        let _ = f__need({ let t1525 = v_a.clone(); Py::Bool({ let t1526 = v_bnd.clone(); if py::le(&t1525, &t1526)? { true } else { false } }) }, Py::str("every lower bound must not exceed its upper bound"))?;
    }
    v_F = { let v_S = v_S.clone(); let v_a = v_a.clone(); let v_bnd = v_bnd.clone(); let v_box_only = v_box_only.clone(); let v_c = v_c.clone(); let v_lb = v_lb.clone(); let v_m = v_m.clone(); let v_n = v_n.clone(); let v_name = v_name.clone(); let v_row = v_row.clone(); let v_ub = v_ub.clone(); let v_v = v_v.clone(); let v_w = v_w.clone(); py::func(move |a: &[Py]| -> R {
            #[allow(unused_mut, unused_variables, unused_assignments)]
{
            let mut v_x: Py = a.get(0).cloned().unwrap_or(Py::None);
            return Ok(py::fraction(&[Py::str(&py::repr(&py::float(&v_x.clone())?))])?);
            #[allow(unreachable_code)]
            Ok(Py::None)
        }
        }) };
    if (v_box_only.clone()).truthy() {
        return Ok(Py::dict(vec![(Py::str("objective"), py::float(&py::sum(&{ let mut out = vec![]; for t1528 in py::iter(&(py::zip(&[v_lb.clone(), v_ub.clone(), v_c.clone()])?))? { let mut v_a: Py; let mut v_bnd: Py; let mut v_w: Py; 
let t1529 = (t1528).items()?;
if t1529.len() != 3 { return Err(PyErr::Value("unpack".into())); }
v_a = t1529[0].clone();
v_bnd = t1529[1].clone();
v_w = t1529[2].clone();
out.push(py::minmax(&[py::mul(&py::call(&v_F.clone(), &[v_w.clone()])?, &py::call(&v_F.clone(), &[v_a.clone()])?)?, py::mul(&py::call(&v_F.clone(), &[v_w.clone()])?, &py::call(&v_F.clone(), &[v_bnd.clone()])?)?], true)?); } Py::list(out) }, None)?)?)]));
    }
    let t1532 = (Py::tuple(vec![{ let mut out = vec![]; for t1530 in py::iter(&(v_lb.clone()))? { let mut v_a: Py; 
v_a = t1530;
out.push(py::call(&v_F.clone(), &[v_a.clone()])?); } Py::list(out) }, { let mut out = vec![]; for t1531 in py::iter(&(v_ub.clone()))? { let mut v_bnd: Py; 
v_bnd = t1531;
out.push(py::call(&v_F.clone(), &[v_bnd.clone()])?); } Py::list(out) }])).items()?;
    if t1532.len() != 2 { return Err(PyErr::Value("unpack".into())); }
    v_L = t1532[0].clone();
    v_U = t1532[1].clone();
    v_Sf = { let mut out = vec![]; for t1533 in py::iter(&(v_S.clone()))? { let mut v_row: Py; 
v_row = t1533;
out.push({ let mut out = vec![]; for t1534 in py::iter(&(v_row.clone()))? { let mut v_x: Py; 
v_x = t1534;
out.push(py::call(&v_F.clone(), &[v_x.clone()])?); } Py::list(out) }); } Py::list(out) };
    v_A = { let mut out = vec![]; for t1535 in py::iter(&(v_Sf.clone()))? { let mut v_row: Py; 
v_row = t1535;
out.push(py::add(&v_row.clone(), &py::mul(&Py::list(vec![py::fraction(&[Py::Int(0i64)])?]), &v_n.clone())?)?); } Py::list(out) };
    v_bvec = { let mut out = vec![]; for t1536 in py::iter(&(v_Sf.clone()))? { let mut v_row: Py; 
v_row = t1536;
out.push(py::neg(&py::sum(&{ let mut out = vec![]; for t1537 in py::range_iter(&[(v_n.clone()).i()?]) { let mut v_j: Py; 
v_j = t1537;
out.push(py::mul(&py::getitem(&v_row.clone(), &v_j.clone())?, &py::getitem(&v_L.clone(), &v_j.clone())?)?); } Py::list(out) }, None)?)?); } Py::list(out) };
    for t1538 in py::range_iter(&[(v_n.clone()).i()?]) {
        v_j = t1538;
        let _ = py::method(&v_A.clone(), "append", &[py::add(&{ let mut out = vec![]; for t1539 in py::range_iter(&[(v_n.clone()).i()?]) { let mut v_k: Py; 
v_k = t1539;
out.push(py::fraction(&[(if ({ let t1540 = v_k.clone(); Py::Bool({ let t1541 = v_j.clone(); if py::eq(&t1540, &t1541) { true } else { false } }) }).truthy() { Py::Int(1i64) } else { Py::Int(0i64) })])?); } Py::list(out) }, &{ let mut out = vec![]; for t1542 in py::range_iter(&[(v_n.clone()).i()?]) { let mut v_k: Py; 
v_k = t1542;
out.push(py::fraction(&[(if ({ let t1543 = v_k.clone(); Py::Bool({ let t1544 = v_j.clone(); if py::eq(&t1543, &t1544) { true } else { false } }) }).truthy() { Py::Int(1i64) } else { Py::Int(0i64) })])?); } Py::list(out) })?])?;
        let _ = py::method(&v_bvec.clone(), "append", &[py::sub(&py::getitem(&v_U.clone(), &v_j.clone())?, &py::getitem(&v_L.clone(), &v_j.clone())?)?])?;
    }
    v_cvec = py::add(&{ let mut out = vec![]; for t1545 in py::iter(&(v_c.clone()))? { let mut v_w: Py; 
v_w = t1545;
out.push(py::call(&v_F.clone(), &[v_w.clone()])?); } Py::list(out) }, &py::mul(&Py::list(vec![py::fraction(&[Py::Int(0i64)])?]), &v_n.clone())?)?;
    v_z = py::add(&f__simplex_max(v_A.clone(), v_bvec.clone(), v_cvec.clone())?, &py::sum(&{ let mut out = vec![]; for t1600 in py::iter(&(py::zip(&[v_c.clone(), v_L.clone()])?))? { let mut v_l: Py; let mut v_w: Py; 
let t1601 = (t1600).items()?;
if t1601.len() != 2 { return Err(PyErr::Value("unpack".into())); }
v_w = t1601[0].clone();
v_l = t1601[1].clone();
out.push(py::mul(&py::call(&v_F.clone(), &[v_w.clone()])?, &v_l.clone())?); } Py::list(out) }, None)?)?;
    return Ok(Py::dict(vec![(Py::str("objective"), py::float(&v_z.clone())?)]));
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn f_fitch_parsimony(mut v_tree: Py, mut v_states: Py, mut v_lower_bound: Py) -> R {
    let mut v_L: Py = Py::None;
    let mut v_alphabet: Py = Py::None;
    let mut v_binary: Py = Py::None;
    let mut v_col: Py = Py::None;
    let mut v_cost: Py = Py::None;
    let mut v_leaves: Py = Py::None;
    let mut v_root: Py = Py::None;
    let mut v_seqs: Py = Py::None;
    let mut v_site: Py = Py::None;
    let mut v_total: Py = Py::None;
    let t1434 = (f__parse_newick(v_tree.clone())?).items()?;
    if t1434.len() != 2 { return Err(PyErr::Value("unpack".into())); }
    v_root = t1434[0].clone();
    v_leaves = t1434[1].clone();
    let rec1435 = py::rec_cell();
    v_binary = { let v_leaves = v_leaves.clone(); let v_lower_bound = v_lower_bound.clone(); let v_root = v_root.clone(); let v_states = v_states.clone(); let v_tree = v_tree.clone(); let rec1435 = rec1435.clone(); py::func(move |a: &[Py]| -> R {
            #[allow(unused_mut, unused_variables, unused_assignments)]
{
            let mut v_t: Py = a.get(0).cloned().unwrap_or(Py::None);
            let mut v_top: Py = a.get(1).cloned().unwrap_or(Py::None);
            let v_binary: Py = py::rec_get(&rec1435)?;
            let mut v_k: Py = Py::None;
            if (py::getitem(&v_t.clone(), &Py::Int(2i64))?).truthy() {
                let _ = f__need({ let t1436 = py::len(&py::getitem(&v_t.clone(), &Py::Int(2i64))?)?; Py::Bool({ let t1437 = Py::Int(2i64); if py::eq(&t1436, &t1437) { true } else { false } }) }, Py::str("the tree must be rooted and strictly bifurcating (two children per internal node)"))?;
                for t1438 in py::iter(&(py::getitem(&v_t.clone(), &Py::Int(2i64))?))? {
                    v_k = t1438;
                    let _ = py::call(&v_binary.clone(), &[v_k.clone(), Py::Bool(false)])?;
                }
            }
            #[allow(unreachable_code)]
            Ok(Py::None)
        }
        }) };
    py::rec_set(&rec1435, &v_binary);
    let _ = py::call(&v_binary.clone(), &[v_root.clone(), Py::Bool(true)])?;
    v_seqs = f__states(v_states.clone(), v_leaves.clone())?;
    v_L = py::len(&(match (py::method(&v_seqs.clone(), "values", &[])?).items()?.into_iter().next() { Some(x) => x, None => return Err(PyErr::Other("StopIteration".into())) }))?;
    v_total = Py::Int(0i64);
    for t1449 in py::range_iter(&[(v_L.clone()).i()?]) {
        v_site = t1449;
        v_col = { let mut out = vec![]; for t1450 in py::iter(&(py::method(&v_seqs.clone(), "items", &[])?))? { let mut v_n: Py; let mut v_s: Py; 
let t1451 = (t1450).items()?;
if t1451.len() != 2 { return Err(PyErr::Value("unpack".into())); }
v_n = t1451[0].clone();
v_s = t1451[1].clone();
out.push((v_n.clone(), py::getitem(&v_s.clone(), &v_site.clone())?)); } Py::dict(out) };
        v_alphabet = py::sorted(&py::set(&py::method(&v_col.clone(), "values", &[])?)?, (Py::Bool(false)).truthy())?;
        if (v_lower_bound.clone()).truthy() {
            v_total = py::iadd(&v_total, &py::sub(&py::len(&v_alphabet.clone())?, &Py::Int(1i64))?)?;
            continue;
        }
        let rec1452 = py::rec_cell();
        v_cost = { let v_L = v_L.clone(); let v_alphabet = v_alphabet.clone(); let v_binary = v_binary.clone(); let v_col = v_col.clone(); let v_leaves = v_leaves.clone(); let v_lower_bound = v_lower_bound.clone(); let v_root = v_root.clone(); let v_seqs = v_seqs.clone(); let v_site = v_site.clone(); let v_states = v_states.clone(); let v_total = v_total.clone(); let v_tree = v_tree.clone(); let rec1452 = rec1452.clone(); py::func(move |a: &[Py]| -> R {
            #[allow(unused_mut, unused_variables, unused_assignments)]
{
            let mut v_t: Py = a.get(0).cloned().unwrap_or(Py::None);
            let v_cost: Py = py::rec_get(&rec1452)?;
            let mut v_cs: Py = Py::None;
            if (py::not(&py::getitem(&v_t.clone(), &Py::Int(2i64))?)).truthy() {
                return Ok({ let mut out = vec![]; for t1453 in py::iter(&(v_alphabet.clone()))? { let mut v_a: Py; 
v_a = t1453;
out.push((v_a.clone(), (if ({ let t1454 = v_a.clone(); Py::Bool({ let t1455 = py::getitem(&v_col.clone(), &py::getitem(&v_t.clone(), &Py::Int(0i64))?)?; if py::eq(&t1454, &t1455) { true } else { false } }) }).truthy() { Py::Int(0i64) } else { Py::Float(f64::INFINITY) }))); } Py::dict(out) });
            }
            v_cs = { let mut out = vec![]; for t1456 in py::iter(&(py::getitem(&v_t.clone(), &Py::Int(2i64))?))? { let mut v_k: Py; 
v_k = t1456;
out.push(py::call(&v_cost.clone(), &[v_k.clone()])?); } Py::list(out) };
            return Ok({ let mut out = vec![]; for t1457 in py::iter(&(v_alphabet.clone()))? { let mut v_a: Py; 
v_a = t1457;
out.push((v_a.clone(), py::sum(&{ let mut out = vec![]; for t1458 in py::iter(&(v_cs.clone()))? { let mut v_c: Py; 
v_c = t1458;
out.push(py::minmax(&[{ let mut out = vec![]; for t1459 in py::iter(&(v_alphabet.clone()))? { let mut v_b: Py; 
v_b = t1459;
out.push(py::add(&py::getitem(&v_c.clone(), &v_b.clone())?, &{ let t1460 = v_a.clone(); Py::Bool({ let t1461 = v_b.clone(); if py::ne(&t1460, &t1461) { true } else { false } }) })?); } Py::list(out) }], false)?); } Py::list(out) }, None)?)); } Py::dict(out) });
            #[allow(unreachable_code)]
            Ok(Py::None)
        }
        }) };
        py::rec_set(&rec1452, &v_cost);
        v_total = py::iadd(&v_total, &py::int(&py::minmax(&[py::method(&py::call(&v_cost.clone(), &[v_root.clone()])?, "values", &[])?], false)?)?)?;
    }
    return Ok(Py::dict(vec![(Py::str("score"), v_total.clone())]));
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn f_fst(mut v_p: Py) -> R {
    let mut v_hs: Py = Py::None;
    let mut v_ht: Py = Py::None;
    let mut v_pb: Py = Py::None;
    let mut v_v: Py = Py::None;
    let _ = f__need({ let t1645 = Py::Bool(matches!(v_p.clone(), Py::List(_))); if !t1645.truthy() { t1645 } else { { let t1643 = py::len(&v_p.clone())?; Py::Bool({ let t1644 = Py::Int(2i64); if py::ge(&t1643, &t1644)? { true } else { false } }) } } }, Py::str("give the allele frequency in at least 2 subpopulations"))?;
    for t1646 in py::iter(&(v_p.clone()))? {
        v_v = t1646;
        let _ = f__frac(Py::str("p"), v_v.clone(), Py::Bool(false), Py::Bool(false))?;
    }
    v_pb = py::truediv(&py::sum(&v_p.clone(), None)?, &py::len(&v_p.clone())?)?;
    let _ = f__need({ let t1659 = Py::Int(0i64); Py::Bool({ let t1660 = v_pb.clone(); if py::lt(&t1659, &t1660)? { { let t1661 = Py::Int(1i64); if py::lt(&t1660, &t1661)? { true } else { false } } } else { false } }) }, Py::str("the allele must be polymorphic overall (0 < mean p < 1)"))?;
    v_hs = py::truediv(&py::sum(&{ let mut out = vec![]; for t1662 in py::iter(&(v_p.clone()))? { let mut v_v: Py; 
v_v = t1662;
out.push(py::mul(&py::mul(&Py::Int(2i64), &v_v.clone())?, &py::sub(&Py::Int(1i64), &v_v.clone())?)?); } Py::list(out) }, None)?, &py::len(&v_p.clone())?)?;
    v_ht = py::mul(&py::mul(&Py::Int(2i64), &v_pb.clone())?, &py::sub(&Py::Int(1i64), &v_pb.clone())?)?;
    return Ok(f__out(Py::dict(vec![(Py::str("HS"), v_hs.clone()), (Py::str("HT"), v_ht.clone()), (Py::str("FST"), py::truediv(&py::sub(&v_ht.clone(), &v_hs.clone())?, &v_ht.clone())?)]))?);
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn f_gc_content(mut v_seq: Py) -> R {
    let mut v_gc: Py = Py::None;
    let mut v_s: Py = Py::None;
    v_s = f__seq(v_seq.clone(), Py::str("seq"), Py::str("NUC"))?;
    v_gc = py::add(&py::method(&v_s.clone(), "count", &[Py::str("G")])?, &py::method(&v_s.clone(), "count", &[Py::str("C")])?)?;
    return Ok(Py::dict(vec![(Py::str("gc_fraction"), py::truediv(&v_gc.clone(), &py::len(&v_s.clone())?)?), (Py::str("gc_count"), v_gc.clone())]));
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn f_global_affine(mut v_seq1: Py, mut v_seq2: Py, mut v_match: Py, mut v_mismatch: Py, mut v_open: Py, mut v_extend: Py) -> R {
    let mut v_a: Py = Py::None;
    let mut v_b: Py = Py::None;
    let t1805 = (Py::tuple(vec![f__seq(v_seq1.clone(), Py::str("seq1"), Py::str("ANY"))?, f__seq(v_seq2.clone(), Py::str("seq2"), Py::str("ANY"))?])).items()?;
    if t1805.len() != 2 { return Err(PyErr::Value("unpack".into())); }
    v_a = t1805[0].clone();
    v_b = t1805[1].clone();
    let _ = f__check_align(v_a.clone(), v_b.clone(), Py::dict(vec![(Py::str("match"), v_match.clone()), (Py::str("mismatch"), v_mismatch.clone()), (Py::str("open"), v_open.clone()), (Py::str("extend"), v_extend.clone())]))?;
    let _ = f__need({ let t1817 = { let t1813 = v_open.clone(); Py::Bool({ let t1814 = Py::Int(0i64); if py::le(&t1813, &t1814)? { true } else { false } }) }; if !t1817.truthy() { t1817 } else { { let t1815 = v_extend.clone(); Py::Bool({ let t1816 = Py::Int(0i64); if py::le(&t1815, &t1816)? { true } else { false } }) } } }, Py::str("gap scores must be zero or negative"))?;
    return Ok(f__out(Py::dict(vec![(Py::str("score"), f__align(v_a.clone(), v_b.clone(), f__mm(v_match.clone(), v_mismatch.clone())?, v_open.clone(), v_extend.clone(), Py::str("global"), Py::Bool(false))?)]))?);
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn f_global_linear(mut v_seq1: Py, mut v_seq2: Py, mut v_match: Py, mut v_mismatch: Py, mut v_gap: Py) -> R {
    let mut v_a: Py = Py::None;
    let mut v_b: Py = Py::None;
    let t1927 = (Py::tuple(vec![f__seq(v_seq1.clone(), Py::str("seq1"), Py::str("ANY"))?, f__seq(v_seq2.clone(), Py::str("seq2"), Py::str("ANY"))?])).items()?;
    if t1927.len() != 2 { return Err(PyErr::Value("unpack".into())); }
    v_a = t1927[0].clone();
    v_b = t1927[1].clone();
    let _ = f__check_align(v_a.clone(), v_b.clone(), Py::dict(vec![(Py::str("match"), v_match.clone()), (Py::str("mismatch"), v_mismatch.clone()), (Py::str("gap"), v_gap.clone())]))?;
    let _ = f__need({ let t1935 = v_gap.clone(); Py::Bool({ let t1936 = Py::Int(0i64); if py::le(&t1935, &t1936)? { true } else { false } }) }, Py::str("the gap score must be zero or negative"))?;
    return Ok(f__out(Py::dict(vec![(Py::str("score"), f__align(v_a.clone(), v_b.clone(), f__mm(v_match.clone(), v_mismatch.clone())?, v_gap.clone(), v_gap.clone(), Py::str("global"), Py::Bool(false))?)]))?);
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn f_goods_coverage(mut v_counts: Py) -> R {
    let mut v_c: Py = Py::None;
    v_c = f__abund(v_counts.clone(), Py::str("counts"))?;
    return Ok(Py::dict(vec![(Py::str("coverage"), py::float(&py::sub(&Py::Int(1i64), &py::fraction(&[py::sum(&{ let mut out = vec![]; for t2041 in py::iter(&(v_c.clone()))? { let mut v_x: Py; 
v_x = t2041;
if ({ let t2042 = v_x.clone(); Py::Bool({ let t2043 = Py::Int(1i64); if py::eq(&t2042, &t2043) { true } else { false } }) }).truthy() { out.push(Py::Int(1i64)); }} Py::list(out) }, None)?, py::sum(&v_c.clone(), None)?])?)?)?)]));
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn f_gravy(mut v_seq: Py) -> R {
    let mut v_s: Py = Py::None;
    v_s = f__seq(v_seq.clone(), Py::str("seq"), Py::str("PROTEIN"))?;
    return Ok(f__out(Py::dict(vec![(Py::str("gravy"), py::truediv(&py::sum(&{ let mut out = vec![]; for t2122 in py::iter(&(v_s.clone()))? { let mut v_a: Py; 
v_a = t2122;
out.push(py::getitem(&c_KD()?, &v_a.clone())?); } Py::list(out) }, None)?, &py::len(&v_s.clone())?)?)]))?);
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn f_hamming(mut v_seq1: Py, mut v_seq2: Py) -> R {
    let mut v_a: Py = Py::None;
    let mut v_b: Py = Py::None;
    let mut v_d: Py = Py::None;
    let t2187 = (Py::tuple(vec![f__seq(v_seq1.clone(), Py::str("seq1"), Py::str("ANY"))?, f__seq(v_seq2.clone(), Py::str("seq2"), Py::str("ANY"))?])).items()?;
    if t2187.len() != 2 { return Err(PyErr::Value("unpack".into())); }
    v_a = t2187[0].clone();
    v_b = t2187[1].clone();
    let _ = f__need({ let t2188 = py::len(&v_a.clone())?; Py::Bool({ let t2189 = py::len(&v_b.clone())?; if py::eq(&t2188, &t2189) { true } else { false } }) }, Py::str("the two aligned sequences must have the same length"))?;
    v_d = py::sum(&{ let mut out = vec![]; for t2190 in py::iter(&(py::zip(&[v_a.clone(), v_b.clone()])?))? { let mut v_x: Py; let mut v_y: Py; 
let t2191 = (t2190).items()?;
if t2191.len() != 2 { return Err(PyErr::Value("unpack".into())); }
v_x = t2191[0].clone();
v_y = t2191[1].clone();
if ({ let t2192 = v_x.clone(); Py::Bool({ let t2193 = v_y.clone(); if py::ne(&t2192, &t2193) { true } else { false } }) }).truthy() { out.push(Py::Int(1i64)); }} Py::list(out) }, None)?;
    return Ok(Py::dict(vec![(Py::str("mismatches"), v_d.clone()), (Py::str("p_distance"), py::truediv(&v_d.clone(), &py::len(&v_a.clone())?)?), (Py::str("identity"), py::truediv(&py::sub(&py::len(&v_a.clone())?, &v_d.clone())?, &py::len(&v_a.clone())?)?)]));
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn f_hanes_woolf(mut v_S: Py, mut v_v: Py) -> R {
    let mut v_a: Py = Py::None;
    let mut v_b: Py = Py::None;
    let mut v_r2: Py = Py::None;
    let t2251 = (f__kin_data(v_S.clone(), v_v.clone(), Py::Int(2i64))?).items()?;
    if t2251.len() != 2 { return Err(PyErr::Value("unpack".into())); }
    v_S = t2251[0].clone();
    v_v = t2251[1].clone();
    let t2263 = (f__ols(v_S.clone(), { let mut out = vec![]; for t2261 in py::iter(&(py::zip(&[v_S.clone(), v_v.clone()])?))? { let mut v_x: Py; let mut v_y: Py; 
let t2262 = (t2261).items()?;
if t2262.len() != 2 { return Err(PyErr::Value("unpack".into())); }
v_x = t2262[0].clone();
v_y = t2262[1].clone();
out.push(py::truediv(&v_x.clone(), &v_y.clone())?); } Py::list(out) })?).items()?;
    if t2263.len() != 3 { return Err(PyErr::Value("unpack".into())); }
    v_b = t2263[0].clone();
    v_a = t2263[1].clone();
    v_r2 = t2263[2].clone();
    let _ = f__need({ let t2268 = { let t2264 = v_a.clone(); Py::Bool({ let t2265 = Py::Int(0i64); if py::gt(&t2264, &t2265)? { true } else { false } }) }; if !t2268.truthy() { t2268 } else { { let t2266 = v_b.clone(); Py::Bool({ let t2267 = Py::Int(0i64); if py::gt(&t2266, &t2267)? { true } else { false } }) } } }, Py::str("the Hanes-Woolf line must have a positive slope and intercept"))?;
    return Ok(f__out(Py::dict(vec![(Py::str("Vmax"), py::truediv(&Py::Int(1i64), &v_b.clone())?), (Py::str("Km"), py::truediv(&v_a.clone(), &v_b.clone())?), (Py::str("r2"), v_r2.clone())]))?);
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn f_haploid_selection(mut v_p0: Py, mut v_s: Py, mut v_t: Py) -> R {
    let mut v__: Py = Py::None;
    let mut v_p: Py = Py::None;
    let _ = f__frac(Py::str("p0"), v_p0.clone(), Py::Bool(true), Py::Bool(true))?;
    let _ = f__fin(Py::str("s"), Py::tuple(vec![v_s.clone()]))?;
    let _ = f__need({ let t2336 = v_s.clone(); Py::Bool({ let t2337 = py::neg(&Py::Int(1i64))?; if py::gt(&t2336, &t2337)? { true } else { false } }) }, Py::str("s must be greater than -1 (fitness 1 + s must be positive)"))?;
    let _ = f__need({ let t2341 = Py::Bool(matches!(v_t.clone(), Py::Int(_) | Py::Big(_) | Py::Bool(_))); if !t2341.truthy() { t2341 } else { { let t2338 = Py::Int(0i64); Py::Bool({ let t2339 = v_t.clone(); if py::le(&t2338, &t2339)? { { let t2340 = Py::Int(100000i64); if py::le(&t2339, &t2340)? { true } else { false } } } else { false } }) } } }, Py::str("t must be a whole number of generations, 0 to 100000"))?;
    v_p = v_p0.clone();
    for t2342 in py::range_iter(&[(v_t.clone()).i()?]) {
        v__ = t2342;
        v_p = py::truediv(&py::mul(&v_p.clone(), &py::add(&Py::Int(1i64), &v_s.clone())?)?, &py::add(&Py::Int(1i64), &py::mul(&v_s.clone(), &v_p.clone())?)?)?;
    }
    return Ok(f__out(Py::dict(vec![(Py::str("p_t"), v_p.clone())]))?);
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn f_hardy_weinberg_test(mut v_n_AA: Py, mut v_n_Aa: Py, mut v_n_aa: Py) -> R {
    let mut v_N: Py = Py::None;
    let mut v_chi2: Py = Py::None;
    let mut v_exp: Py = Py::None;
    let mut v_he: Py = Py::None;
    let mut v_ho: Py = Py::None;
    let mut v_nm: Py = Py::None;
    let mut v_obs: Py = Py::None;
    let mut v_p: Py = Py::None;
    let mut v_q: Py = Py::None;
    let mut v_v: Py = Py::None;
    for t2455 in py::iter(&(Py::tuple(vec![Py::tuple(vec![Py::str("n_AA"), v_n_AA.clone()]), Py::tuple(vec![Py::str("n_Aa"), v_n_Aa.clone()]), Py::tuple(vec![Py::str("n_aa"), v_n_aa.clone()])])))? {
        let t2459 = (t2455).items()?;
        if t2459.len() != 2 { return Err(PyErr::Value("unpack".into())); }
        v_nm = t2459[0].clone();
        v_v = t2459[1].clone();
        let _ = f__need({ let t2458 = Py::Bool(matches!(v_v.clone(), Py::Int(_) | Py::Big(_) | Py::Bool(_))); if !t2458.truthy() { t2458 } else { { let t2456 = v_v.clone(); Py::Bool({ let t2457 = Py::Int(0i64); if py::ge(&t2456, &t2457)? { true } else { false } }) } } }, Py::str(&[py::format(&v_nm.clone(), "")?, " must be a whole number >= 0".to_string()].concat()))?;
    }
    v_N = py::add(&py::add(&v_n_AA.clone(), &v_n_Aa.clone())?, &v_n_aa.clone())?;
    let _ = f__need({ let t2460 = v_N.clone(); Py::Bool({ let t2461 = Py::Int(1i64); if py::ge(&t2460, &t2461)? { true } else { false } }) }, Py::str("at least one individual is needed"))?;
    v_p = py::fraction(&[py::add(&py::mul(&Py::Int(2i64), &v_n_AA.clone())?, &v_n_Aa.clone())?, py::mul(&Py::Int(2i64), &v_N.clone())?])?;
    let _ = f__need({ let t2462 = Py::Int(0i64); Py::Bool({ let t2463 = v_p.clone(); if py::lt(&t2462, &t2463)? { { let t2464 = Py::Int(1i64); if py::lt(&t2463, &t2464)? { true } else { false } } } else { false } }) }, Py::str("both alleles must be present (otherwise the expected counts are zero)"))?;
    v_q = py::sub(&Py::Int(1i64), &v_p.clone())?;
    v_exp = Py::list(vec![py::mul(&py::mul(&v_N.clone(), &v_p.clone())?, &v_p.clone())?, py::mul(&py::mul(&py::mul(&Py::Int(2i64), &v_N.clone())?, &v_p.clone())?, &v_q.clone())?, py::mul(&py::mul(&v_N.clone(), &v_q.clone())?, &v_q.clone())?]);
    v_obs = Py::list(vec![v_n_AA.clone(), v_n_Aa.clone(), v_n_aa.clone()]);
    v_chi2 = py::sum(&{ let mut out = vec![]; for t2465 in py::iter(&(py::zip(&[v_obs.clone(), v_exp.clone()])?))? { let mut v_e: Py; let mut v_o: Py; 
let t2466 = (t2465).items()?;
if t2466.len() != 2 { return Err(PyErr::Value("unpack".into())); }
v_o = t2466[0].clone();
v_e = t2466[1].clone();
out.push(py::truediv(&py::pow(&py::sub(&v_o.clone(), &v_e.clone())?, &Py::Int(2i64))?, &v_e.clone())?); } Py::list(out) }, None)?;
    let t2467 = (Py::tuple(vec![py::fraction(&[v_n_Aa.clone(), v_N.clone()])?, py::mul(&py::mul(&Py::Int(2i64), &v_p.clone())?, &v_q.clone())?])).items()?;
    if t2467.len() != 2 { return Err(PyErr::Value("unpack".into())); }
    v_ho = t2467[0].clone();
    v_he = t2467[1].clone();
    return Ok(f__out(Py::dict(vec![(Py::str("p"), py::float(&v_p.clone())?), (Py::str("q"), py::float(&v_q.clone())?), (Py::str("exp_AA"), py::float(&py::getitem(&v_exp.clone(), &Py::Int(0i64))?)?), (Py::str("exp_Aa"), py::float(&py::getitem(&v_exp.clone(), &Py::Int(1i64))?)?), (Py::str("exp_aa"), py::float(&py::getitem(&v_exp.clone(), &Py::Int(2i64))?)?), (Py::str("chi2"), py::float(&v_chi2.clone())?), (Py::str("p_value"), py::math1("erfc", &py::math1("sqrt", &py::truediv(&py::float(&v_chi2.clone())?, &Py::Int(2i64))?)?)?), (Py::str("Ho"), py::float(&v_ho.clone())?), (Py::str("He"), py::float(&v_he.clone())?), (Py::str("F"), py::float(&py::sub(&Py::Int(1i64), &py::truediv(&v_ho.clone(), &v_he.clone())?)?)?)]))?);
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn f_het_advantage(mut v_w_AA: Py, mut v_w_Aa: Py, mut v_w_aa: Py) -> R {
    let _ = f__pos(Py::str("fitness"), Py::tuple(vec![v_w_AA.clone(), v_w_Aa.clone(), v_w_aa.clone()]))?;
    let _ = f__need({ let t2533 = { let t2529 = v_w_Aa.clone(); Py::Bool({ let t2530 = v_w_AA.clone(); if py::gt(&t2529, &t2530)? { true } else { false } }) }; if !t2533.truthy() { t2533 } else { { let t2531 = v_w_Aa.clone(); Py::Bool({ let t2532 = v_w_aa.clone(); if py::gt(&t2531, &t2532)? { true } else { false } }) } } }, Py::str("heterozygote advantage needs w_Aa greater than both homozygote fitnesses"))?;
    return Ok(f__out(Py::dict(vec![(Py::str("p_eq"), py::truediv(&py::sub(&v_w_Aa.clone(), &v_w_aa.clone())?, &py::sub(&py::sub(&py::mul(&Py::Int(2i64), &v_w_Aa.clone())?, &v_w_AA.clone())?, &v_w_aa.clone())?)?)]))?);
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn f_hill_fit(mut v_S: Py, mut v_v: Py, mut v_fixed_n: Py) -> R {
    let mut v_J: Py = Py::None;
    let mut v_Y: Py = Py::None;
    let mut v__: Py = Py::None;
    let mut v_a: Py = Py::None;
    let mut v_b: Py = Py::None;
    let mut v_f: Py = Py::None;
    let mut v_g: Py = Py::None;
    let mut v_p: Py = Py::None;
    let mut v_p0: Py = Py::None;
    let mut v_r: Py = Py::None;
    let mut v_scale: Py = Py::None;
    let mut v_vmax0: Py = Py::None;
    let t2605 = (f__kin_data(v_S.clone(), v_v.clone(), Py::Int(4i64))?).items()?;
    if t2605.len() != 2 { return Err(PyErr::Value("unpack".into())); }
    v_S = t2605[0].clone();
    v_v = t2605[1].clone();
    let _ = f__need({ let t2606 = py::len(&py::set(&v_S.clone())?)?; Py::Bool({ let t2607 = Py::Int(4i64); if py::ge(&t2606, &t2607)? { true } else { false } }) }, Py::str("at least 4 different concentrations are needed"))?;
    if (v_fixed_n.clone()).truthy() {
        v_r = f_mm_fit(v_S.clone(), v_v.clone())?;
        return Ok(Py::dict(vec![(Py::str("Vmax"), py::getitem(&v_r.clone(), &Py::str("Vmax"))?), (Py::str("K"), py::getitem(&v_r.clone(), &Py::str("Km"))?), (Py::str("n"), Py::Float(1.0))]));
    }
    v_vmax0 = py::mul(&py::minmax(&[v_v.clone()], true)?, &Py::Float(1.1))?;
    v_Y = { let mut out = vec![]; for t2668 in py::iter(&(v_v.clone()))? { let mut v_y: Py; 
v_y = t2668;
out.push(py::minmax(&[py::minmax(&[py::truediv(&v_y.clone(), &v_vmax0.clone())?, Py::Float(0.001)], true)?, py::sub(&Py::Int(1i64), &Py::Float(0.001))?], false)?); } Py::list(out) };
    let t2680 = (f__ols({ let mut out = vec![]; for t2678 in py::iter(&(v_S.clone()))? { let mut v_s: Py; 
v_s = t2678;
out.push(py::math1("log", &v_s.clone())?); } Py::list(out) }, { let mut out = vec![]; for t2679 in py::iter(&(v_Y.clone()))? { let mut v_y: Py; 
v_y = t2679;
out.push(py::math1("log", &py::truediv(&v_y.clone(), &py::sub(&Py::Int(1i64), &v_y.clone())?)?)?); } Py::list(out) })?).items()?;
    if t2680.len() != 3 { return Err(PyErr::Value("unpack".into())); }
    v_b = t2680[0].clone();
    v_a = t2680[1].clone();
    v__ = t2680[2].clone();
    v_p0 = Py::list(vec![v_vmax0.clone(), (if ({ let t2681 = v_b.clone(); Py::Bool({ let t2682 = Py::Int(0i64); if py::gt(&t2681, &t2682)? { true } else { false } }) }).truthy() { py::math1("exp", &py::truediv(&py::neg(&v_a.clone())?, &v_b.clone())?)? } else { py::getitem(&py::sorted(&v_S.clone(), (Py::Bool(false)).truthy())?, &py::floordiv(&py::len(&v_S.clone())?, &Py::Int(2i64))?)? }), py::minmax(&[v_b.clone(), Py::Float(0.5)], true)?]);
    v_p = f__lm(f__hill_model(v_S.clone())?, v_p0.clone(), v_v.clone(), { let v_S = v_S.clone(); let v_Y = v_Y.clone(); let v__ = v__.clone(); let v_a = v_a.clone(); let v_b = v_b.clone(); let v_fixed_n = v_fixed_n.clone(); let v_p0 = v_p0.clone(); let v_r = v_r.clone(); let v_v = v_v.clone(); let v_vmax0 = v_vmax0.clone(); py::func(move |a: &[Py]| -> R {
            #[allow(unused_mut, unused_variables, unused_assignments)]
{
            let mut v_q: Py = a.get(0).cloned().unwrap_or(Py::None);
            return Ok({ let t2752 = { let t2744 = py::getitem(&v_q.clone(), &Py::Int(0i64))?; Py::Bool({ let t2745 = Py::Int(0i64); if py::gt(&t2744, &t2745)? { true } else { false } }) }; if !t2752.truthy() { t2752 } else { { let t2751 = { let t2746 = py::getitem(&v_q.clone(), &Py::Int(1i64))?; Py::Bool({ let t2747 = Py::Int(0i64); if py::gt(&t2746, &t2747)? { true } else { false } }) }; if !t2751.truthy() { t2751 } else { { let t2748 = Py::Float(0.05); Py::Bool({ let t2749 = py::getitem(&v_q.clone(), &Py::Int(2i64))?; if py::lt(&t2748, &t2749)? { { let t2750 = Py::Int(20i64); if py::lt(&t2749, &t2750)? { true } else { false } } } else { false } }) } } } } });
            #[allow(unreachable_code)]
            Ok(Py::None)
        }
        }) })?;
    let t2753 = (py::call(&f__hill_model(v_S.clone())?, &[v_p.clone()])?).items()?;
    if t2753.len() != 2 { return Err(PyErr::Value("unpack".into())); }
    v_f = t2753[0].clone();
    v_J = t2753[1].clone();
    v_r = { let mut out = vec![]; for t2754 in py::iter(&(py::zip(&[v_v.clone(), v_f.clone()])?))? { let mut v_a: Py; let mut v_y: Py; 
let t2755 = (t2754).items()?;
if t2755.len() != 2 { return Err(PyErr::Value("unpack".into())); }
v_y = t2755[0].clone();
v_a = t2755[1].clone();
out.push(py::sub(&v_y.clone(), &v_a.clone())?); } Py::list(out) };
    v_g = { let mut out = vec![]; for t2756 in py::range_iter(&[(Py::Int(3i64)).i()?]) { let mut v_i: Py; 
v_i = t2756;
out.push(py::abs(&py::fsum(&{ let mut out = vec![]; for t2757 in py::iter(&(py::zip(&[v_J.clone(), v_r.clone()])?))? { let mut v_ri: Py; let mut v_row: Py; 
let t2758 = (t2757).items()?;
if t2758.len() != 2 { return Err(PyErr::Value("unpack".into())); }
v_row = t2758[0].clone();
v_ri = t2758[1].clone();
out.push(py::mul(&py::getitem(&v_row.clone(), &v_i.clone())?, &v_ri.clone())?); } Py::list(out) })?)?); } Py::list(out) };
    v_scale = py::fsum(&{ let mut out = vec![]; for t2759 in py::iter(&(v_v.clone()))? { let mut v_y: Py; 
v_y = t2759;
out.push(py::mul(&v_y.clone(), &v_y.clone())?); } Py::list(out) })?;
    let _ = f__need({ let t2760 = py::minmax(&[{ let mut out = vec![]; for t2762 in py::iter(&(py::zip(&[v_g.clone(), v_p.clone()])?))? { let mut v_gi: Py; let mut v_pi: Py; 
let t2763 = (t2762).items()?;
if t2763.len() != 2 { return Err(PyErr::Value("unpack".into())); }
v_gi = t2763[0].clone();
v_pi = t2763[1].clone();
out.push(py::mul(&v_gi.clone(), &py::abs(&v_pi.clone())?)?); } Py::list(out) }], true)?; Py::Bool({ let t2761 = py::mul(&Py::Float(1e-09), &v_scale.clone())?; if py::le(&t2760, &t2761)? { true } else { false } }) }, Py::str("the Hill fit did not converge"))?;
    return Ok(f__out(Py::dict(vec![(Py::str("Vmax"), py::getitem(&v_p.clone(), &Py::Int(0i64))?), (Py::str("K"), py::getitem(&v_p.clone(), &Py::Int(1i64))?), (Py::str("n"), py::getitem(&v_p.clone(), &Py::Int(2i64))?)]))?);
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn f_hill_plot(mut v_S: Py, mut v_Y: Py) -> R {
    let mut v_a: Py = Py::None;
    let mut v_b: Py = Py::None;
    let mut v_r2: Py = Py::None;
    let mut v_s: Py = Py::None;
    let mut v_y: Py = Py::None;
    let _ = f__need({ let t2811 = Py::Bool(matches!(v_S.clone(), Py::List(_))); if !t2811.truthy() { t2811 } else { { let t2810 = Py::Bool(matches!(v_Y.clone(), Py::List(_))); if !t2810.truthy() { t2810 } else { { let t2809 = { let t2805 = py::len(&v_S.clone())?; Py::Bool({ let t2806 = py::len(&v_Y.clone())?; if py::eq(&t2805, &t2806) { true } else { false } }) }; if !t2809.truthy() { t2809 } else { { let t2807 = py::len(&v_S.clone())?; Py::Bool({ let t2808 = Py::Int(2i64); if py::ge(&t2807, &t2808)? { true } else { false } }) } } } } } } }, Py::str("S and Y must be lists of equal length (at least 2 points)"))?;
    for t2812 in py::iter(&(v_S.clone()))? {
        v_s = t2812;
        let _ = f__pos(Py::str("S"), Py::tuple(vec![v_s.clone()]))?;
    }
    for t2819 in py::iter(&(v_Y.clone()))? {
        v_y = t2819;
        let _ = f__frac(Py::str("Y"), v_y.clone(), Py::Bool(true), Py::Bool(true))?;
    }
    let t2843 = (f__ols({ let mut out = vec![]; for t2841 in py::iter(&(v_S.clone()))? { let mut v_s: Py; 
v_s = t2841;
out.push(py::math1("log10", &v_s.clone())?); } Py::list(out) }, { let mut out = vec![]; for t2842 in py::iter(&(v_Y.clone()))? { let mut v_y: Py; 
v_y = t2842;
out.push(py::math1("log10", &py::truediv(&v_y.clone(), &py::sub(&Py::Int(1i64), &v_y.clone())?)?)?); } Py::list(out) })?).items()?;
    if t2843.len() != 3 { return Err(PyErr::Value("unpack".into())); }
    v_b = t2843[0].clone();
    v_a = t2843[1].clone();
    v_r2 = t2843[2].clone();
    let _ = f__need({ let t2844 = v_b.clone(); Py::Bool({ let t2845 = Py::Int(0i64); if py::gt(&t2844, &t2845)? { true } else { false } }) }, Py::str("the Hill plot must rise (positive Hill coefficient)"))?;
    return Ok(f__out(Py::dict(vec![(Py::str("n"), v_b.clone()), (Py::str("K"), py::pow(&Py::Int(10i64), &py::truediv(&py::neg(&v_a.clone())?, &v_b.clone())?)?), (Py::str("r2"), v_r2.clone())]))?);
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn f_hwe_exact(mut v_n_AA: Py, mut v_n_Aa: Py, mut v_n_aa: Py) -> R {
    let mut v_N: Py = Py::None;
    let mut v_nA: Py = Py::None;
    let mut v_nm: Py = Py::None;
    let mut v_po: Py = Py::None;
    let mut v_pr: Py = Py::None;
    let mut v_v: Py = Py::None;
    for t2398 in py::iter(&(Py::tuple(vec![Py::tuple(vec![Py::str("n_AA"), v_n_AA.clone()]), Py::tuple(vec![Py::str("n_Aa"), v_n_Aa.clone()]), Py::tuple(vec![Py::str("n_aa"), v_n_aa.clone()])])))? {
        let t2402 = (t2398).items()?;
        if t2402.len() != 2 { return Err(PyErr::Value("unpack".into())); }
        v_nm = t2402[0].clone();
        v_v = t2402[1].clone();
        let _ = f__need({ let t2401 = Py::Bool(matches!(v_v.clone(), Py::Int(_) | Py::Big(_) | Py::Bool(_))); if !t2401.truthy() { t2401 } else { { let t2399 = v_v.clone(); Py::Bool({ let t2400 = Py::Int(0i64); if py::ge(&t2399, &t2400)? { true } else { false } }) } } }, Py::str(&[py::format(&v_nm.clone(), "")?, " must be a whole number >= 0".to_string()].concat()))?;
    }
    v_N = py::add(&py::add(&v_n_AA.clone(), &v_n_Aa.clone())?, &v_n_aa.clone())?;
    let _ = f__need({ let t2403 = Py::Int(1i64); Py::Bool({ let t2404 = v_N.clone(); if py::le(&t2403, &t2404)? { { let t2405 = Py::Int(1000i64); if py::le(&t2404, &t2405)? { true } else { false } } } else { false } }) }, Py::str("between 1 and 1000 individuals"))?;
    v_nA = py::add(&py::mul(&Py::Int(2i64), &v_n_AA.clone())?, &v_n_Aa.clone())?;
    let _ = f__need({ let t2406 = Py::Int(0i64); Py::Bool({ let t2407 = v_nA.clone(); if py::lt(&t2406, &t2407)? { { let t2408 = py::mul(&Py::Int(2i64), &v_N.clone())?; if py::lt(&t2407, &t2408)? { true } else { false } } } else { false } }) }, Py::str("both alleles must be present"))?;
    v_pr = f__hwe_probs(v_nA.clone(), v_N.clone())?;
    v_po = py::getitem(&v_pr.clone(), &v_n_Aa.clone())?;
    return Ok(Py::dict(vec![(Py::str("p_value"), py::float(&py::sum(&{ let mut out = vec![]; for t2411 in py::iter(&(py::method(&v_pr.clone(), "values", &[])?))? { let mut v_v: Py; 
v_v = t2411;
if ({ let t2412 = v_v.clone(); Py::Bool({ let t2413 = v_po.clone(); if py::le(&t2412, &t2413)? { true } else { false } }) }).truthy() { out.push(v_v.clone()); }} Py::list(out) }, None)?)?), (Py::str("p_observed"), py::float(&v_po.clone())?)]));
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn f_instability(mut v_seq: Py) -> R {
    let mut v_i: Py = Py::None;
    let mut v_s: Py = Py::None;
    let mut v_score: Py = Py::None;
    v_s = f__seq(v_seq.clone(), Py::str("seq"), Py::str("PROTEIN"))?;
    let _ = f__need({ let t2924 = py::len(&v_s.clone())?; Py::Bool({ let t2925 = Py::Int(2i64); if py::ge(&t2924, &t2925)? { true } else { false } }) }, Py::str("the protein needs at least 2 residues"))?;
    v_score = Py::Float(0.0);
    for t2926 in py::range_iter(&[(py::sub(&py::len(&v_s.clone())?, &Py::Int(1i64))?).i()?]) {
        v_i = t2926;
        v_score = py::iadd(&v_score, &py::getitem(&py::getitem(&c_DIWV()?, &py::getitem(&v_s.clone(), &v_i.clone())?)?, &py::getitem(&v_s.clone(), &py::add(&v_i.clone(), &Py::Int(1i64))?)?)?)?;
    }
    return Ok(f__out(Py::dict(vec![(Py::str("index"), py::mul(&py::truediv(&Py::Float(10.0), &py::len(&v_s.clone())?)?, &v_score.clone())?)]))?);
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn f_isoelectric_point(mut v_seq: Py) -> R {
    let mut v_pi: Py = Py::None;
    let mut v_s: Py = Py::None;
    v_s = f__seq(v_seq.clone(), Py::str("seq"), Py::str("PROTEIN"))?;
    v_pi = f__bisect({ let v_s = v_s.clone(); let v_seq = v_seq.clone(); py::func(move |a: &[Py]| -> R {
            #[allow(unused_mut, unused_variables, unused_assignments)]
{
            let mut v_x: Py = a.get(0).cloned().unwrap_or(Py::None);
            return Ok(f__charge(v_s.clone(), v_x.clone(), Py::Bool(true))?);
            #[allow(unreachable_code)]
            Ok(Py::None)
        }
        }) }, Py::Float(0.0), Py::Float(14.0), Py::Int(400i64))?;
    return Ok(f__out(Py::dict(vec![(Py::str("pI"), v_pi.clone()), (Py::str("charge_at_7"), f__charge(v_s.clone(), Py::Float(7.0), Py::Bool(true))?)]))?);
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn f_jukes_cantor(mut v_seq1: Py, mut v_seq2: Py) -> R {
    let mut v_a: Py = Py::None;
    let mut v_b: Py = Py::None;
    let mut v_p: Py = Py::None;
    let t3123 = (f__aligned_dna(v_seq1.clone(), v_seq2.clone())?).items()?;
    if t3123.len() != 2 { return Err(PyErr::Value("unpack".into())); }
    v_a = t3123[0].clone();
    v_b = t3123[1].clone();
    v_p = py::truediv(&py::sum(&{ let mut out = vec![]; for t3124 in py::iter(&(py::zip(&[v_a.clone(), v_b.clone()])?))? { let mut v_x: Py; let mut v_y: Py; 
let t3125 = (t3124).items()?;
if t3125.len() != 2 { return Err(PyErr::Value("unpack".into())); }
v_x = t3125[0].clone();
v_y = t3125[1].clone();
if ({ let t3126 = v_x.clone(); Py::Bool({ let t3127 = v_y.clone(); if py::ne(&t3126, &t3127) { true } else { false } }) }).truthy() { out.push(Py::Int(1i64)); }} Py::list(out) }, None)?, &py::len(&v_a.clone())?)?;
    let _ = f__need({ let t3128 = v_p.clone(); Py::Bool({ let t3129 = Py::Float(0.75); if py::lt(&t3128, &t3129)? { true } else { false } }) }, Py::str("p >= 0.75: the sequences are saturated and the Jukes-Cantor distance is undefined"))?;
    return Ok(f__out(Py::dict(vec![(Py::str("p"), v_p.clone()), (Py::str("d"), py::mul(&py::neg(&Py::Float(0.75))?, &py::math1("log1p", &py::truediv(&py::mul(&py::neg(&Py::Float(4.0))?, &v_p.clone())?, &Py::Float(3.0))?)?)?)]))?);
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn f_kimura(mut v_seq1: Py, mut v_seq2: Py) -> R {
    let mut v_P: Py = Py::None;
    let mut v_Q: Py = Py::None;
    let mut v_a: Py = Py::None;
    let mut v_b: Py = Py::None;
    let mut v_n: Py = Py::None;
    let mut v_ts: Py = Py::None;
    let mut v_tv: Py = Py::None;
    let t3211 = (f__aligned_dna(v_seq1.clone(), v_seq2.clone())?).items()?;
    if t3211.len() != 2 { return Err(PyErr::Value("unpack".into())); }
    v_a = t3211[0].clone();
    v_b = t3211[1].clone();
    v_n = py::len(&v_a.clone())?;
    v_ts = py::sum(&{ let mut out = vec![]; for t3212 in py::iter(&(py::zip(&[v_a.clone(), v_b.clone()])?))? { let mut v_x: Py; let mut v_y: Py; 
let t3213 = (t3212).items()?;
if t3213.len() != 2 { return Err(PyErr::Value("unpack".into())); }
v_x = t3213[0].clone();
v_y = t3213[1].clone();
if ({ let t3222 = { let t3214 = v_x.clone(); Py::Bool({ let t3215 = v_y.clone(); if py::ne(&t3214, &t3215) { true } else { false } }) }; if !t3222.truthy() { t3222 } else { { let t3216 = { let t3218 = v_x.clone(); Py::Bool({ let t3219 = c__PUR()?; if py::contains(&t3219, &t3218)? { true } else { false } }) }; Py::Bool({ let t3217 = { let t3220 = v_y.clone(); Py::Bool({ let t3221 = c__PUR()?; if py::contains(&t3221, &t3220)? { true } else { false } }) }; if py::eq(&t3216, &t3217) { true } else { false } }) } } }).truthy() { out.push(Py::Int(1i64)); }} Py::list(out) }, None)?;
    v_tv = py::sum(&{ let mut out = vec![]; for t3223 in py::iter(&(py::zip(&[v_a.clone(), v_b.clone()])?))? { let mut v_x: Py; let mut v_y: Py; 
let t3224 = (t3223).items()?;
if t3224.len() != 2 { return Err(PyErr::Value("unpack".into())); }
v_x = t3224[0].clone();
v_y = t3224[1].clone();
if ({ let t3225 = { let t3227 = v_x.clone(); Py::Bool({ let t3228 = c__PUR()?; if py::contains(&t3228, &t3227)? { true } else { false } }) }; Py::Bool({ let t3226 = { let t3229 = v_y.clone(); Py::Bool({ let t3230 = c__PUR()?; if py::contains(&t3230, &t3229)? { true } else { false } }) }; if py::ne(&t3225, &t3226) { true } else { false } }) }).truthy() { out.push(Py::Int(1i64)); }} Py::list(out) }, None)?;
    let t3231 = (Py::tuple(vec![py::truediv(&v_ts.clone(), &v_n.clone())?, py::truediv(&v_tv.clone(), &v_n.clone())?])).items()?;
    if t3231.len() != 2 { return Err(PyErr::Value("unpack".into())); }
    v_P = t3231[0].clone();
    v_Q = t3231[1].clone();
    let _ = f__need({ let t3236 = { let t3232 = py::sub(&py::sub(&Py::Int(1i64), &py::mul(&Py::Int(2i64), &v_P.clone())?)?, &v_Q.clone())?; Py::Bool({ let t3233 = Py::Int(0i64); if py::gt(&t3232, &t3233)? { true } else { false } }) }; if !t3236.truthy() { t3236 } else { { let t3234 = py::sub(&Py::Int(1i64), &py::mul(&Py::Int(2i64), &v_Q.clone())?)?; Py::Bool({ let t3235 = Py::Int(0i64); if py::gt(&t3234, &t3235)? { true } else { false } }) } } }, Py::str("too many differences: the Kimura distance is undefined"))?;
    return Ok(f__out(Py::dict(vec![(Py::str("d"), py::sub(&py::mul(&py::neg(&Py::Float(0.5))?, &py::math1("log", &py::sub(&py::sub(&Py::Int(1i64), &py::mul(&Py::Int(2i64), &v_P.clone())?)?, &v_Q.clone())?)?)?, &py::mul(&Py::Float(0.25), &py::math1("log", &py::sub(&Py::Int(1i64), &py::mul(&Py::Int(2i64), &v_Q.clone())?)?)?)?)?)]))?);
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn f_kmer_counts(mut v_seq: Py, mut v_k: Py) -> R {
    let mut v_counts: Py = Py::None;
    let mut v_i: Py = Py::None;
    let mut v_s: Py = Py::None;
    let mut v_top: Py = Py::None;
    let mut v_w: Py = Py::None;
    v_s = f__seq(v_seq.clone(), Py::str("seq"), Py::str("NUC"))?;
    let _ = f__need({ let t3318 = Py::Bool(matches!(v_k.clone(), Py::Int(_) | Py::Big(_) | Py::Bool(_))); if !t3318.truthy() { t3318 } else { { let t3315 = Py::Int(1i64); Py::Bool({ let t3316 = v_k.clone(); if py::le(&t3315, &t3316)? { { let t3317 = py::len(&v_s.clone())?; if py::le(&t3316, &t3317)? { true } else { false } } } else { false } }) } } }, Py::str("k must be a whole number between 1 and the sequence length"))?;
    v_counts = Py::dict(vec![]);
    for t3319 in py::range_iter(&[(py::add(&py::sub(&py::len(&v_s.clone())?, &v_k.clone())?, &Py::Int(1i64))?).i()?]) {
        v_i = t3319;
        v_w = py::slice(&v_s.clone(), Some((v_i.clone()).i()?), Some((py::add(&v_i.clone(), &v_k.clone())?).i()?), None)?;
        py::setitem(&v_counts.clone(), &v_w.clone(), py::add(&py::method(&v_counts.clone(), "get", &[v_w.clone(), Py::Int(0i64)])?, &Py::Int(1i64))?)?;
    }
    v_top = py::minmax(&[py::method(&v_counts.clone(), "values", &[])?], true)?;
    return Ok(Py::dict(vec![(Py::str("total"), py::add(&py::sub(&py::len(&v_s.clone())?, &v_k.clone())?, &Py::Int(1i64))?), (Py::str("distinct"), py::len(&v_counts.clone())?), (Py::str("max_count"), v_top.clone()), (Py::str("most_frequent"), py::minmax(&[{ let mut out = vec![]; for t3320 in py::iter(&(py::method(&v_counts.clone(), "items", &[])?))? { let mut v_c: Py; let mut v_w: Py; 
let t3321 = (t3320).items()?;
if t3321.len() != 2 { return Err(PyErr::Value("unpack".into())); }
v_w = t3321[0].clone();
v_c = t3321[1].clone();
if ({ let t3322 = v_c.clone(); Py::Bool({ let t3323 = v_top.clone(); if py::eq(&t3322, &t3323) { true } else { false } }) }).truthy() { out.push(v_w.clone()); }} Py::list(out) }], false)?)]));
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn f_leslie_growth(mut v_F: Py, mut v_S: Py, mut v_use_r0: Py) -> R {
    let mut v_R0: Py = Py::None;
    let mut v_g: Py = Py::None;
    let mut v_hi: Py = Py::None;
    let mut v_l: Py = Py::None;
    let mut v_lam: Py = Py::None;
    let mut v_lo: Py = Py::None;
    let mut v_tot: Py = Py::None;
    let mut v_w: Py = Py::None;
    let _ = f__leslie_in(v_F.clone(), v_S.clone())?;
    v_l = f__euler_lotka(v_F.clone(), v_S.clone())?;
    v_R0 = py::fsum(&{ let mut out = vec![]; for t3386 in py::iter(&(py::zip(&[v_F.clone(), v_l.clone()])?))? { let mut v_f: Py; let mut v_li: Py; 
let t3387 = (t3386).items()?;
if t3387.len() != 2 { return Err(PyErr::Value("unpack".into())); }
v_f = t3387[0].clone();
v_li = t3387[1].clone();
out.push(py::mul(&v_f.clone(), &v_li.clone())?); } Py::list(out) })?;
    if (v_use_r0.clone()).truthy() {
        v_lam = v_R0.clone();
    } else {
        v_g = { let v_F = v_F.clone(); let v_R0 = v_R0.clone(); let v_S = v_S.clone(); let v_l = v_l.clone(); let v_lam = v_lam.clone(); let v_use_r0 = v_use_r0.clone(); py::func(move |a: &[Py]| -> R {
            #[allow(unused_mut, unused_variables, unused_assignments)]
{
            let mut v_x: Py = a.get(0).cloned().unwrap_or(Py::None);
            return Ok(py::sub(&py::fsum(&{ let mut out = vec![]; for t3388 in py::iter(&(py::enumerate(&py::zip(&[v_F.clone(), v_l.clone()])?, 0)?))? { let mut v_f: Py; let mut v_i: Py; let mut v_li: Py; 
let t3389 = (t3388).items()?;
if t3389.len() != 2 { return Err(PyErr::Value("unpack".into())); }
v_i = t3389[0].clone();
let t3390 = (t3389[1].clone()).items()?;
if t3390.len() != 2 { return Err(PyErr::Value("unpack".into())); }
v_f = t3390[0].clone();
v_li = t3390[1].clone();
out.push(py::mul(&py::mul(&v_f.clone(), &v_li.clone())?, &py::math1("exp", &py::mul(&py::neg(&py::add(&v_i.clone(), &Py::Int(1i64))?)?, &v_x.clone())?)?)?); } Py::list(out) })?, &Py::Float(1.0))?);
            #[allow(unreachable_code)]
            Ok(Py::None)
        }
        }) };
        v_hi = py::add(&py::math1("log", &py::minmax(&[v_R0.clone(), Py::Float(1.0)], true)?)?, &Py::Float(1.0))?;
        v_lo = py::neg(&Py::Float(1.0))?;
        while ({ let t3391 = py::call(&v_g.clone(), &[v_lo.clone()])?; Py::Bool({ let t3392 = Py::Int(0i64); if py::le(&t3391, &t3392)? { true } else { false } }) }).truthy() {
            v_lo = py::sub(&v_lo, &Py::Float(1.0))?;
        }
        while ({ let t3393 = py::call(&v_g.clone(), &[v_hi.clone()])?; Py::Bool({ let t3394 = Py::Int(0i64); if py::ge(&t3393, &t3394)? { true } else { false } }) }).truthy() {
            v_hi = py::iadd(&v_hi, &Py::Float(1.0))?;
        }
        v_lam = py::math1("exp", &f__bisect(v_g.clone(), v_lo.clone(), v_hi.clone(), Py::Int(400i64))?)?;
    }
    v_w = { let mut out = vec![]; for t3421 in py::iter(&(py::enumerate(&v_l.clone(), 0)?))? { let mut v_i: Py; let mut v_li: Py; 
let t3422 = (t3421).items()?;
if t3422.len() != 2 { return Err(PyErr::Value("unpack".into())); }
v_i = t3422[0].clone();
v_li = t3422[1].clone();
out.push(py::mul(&v_li.clone(), &py::pow(&v_lam.clone(), &py::neg(&v_i.clone())?)?)?); } Py::list(out) };
    v_tot = py::fsum(&v_w.clone())?;
    return Ok(f__out(Py::dict(vec![(Py::str("lambda1"), v_lam.clone()), (Py::str("r"), py::math1("log", &v_lam.clone())?), (Py::str("R0"), v_R0.clone()), (Py::str("stable"), { let mut out = vec![]; for t3437 in py::iter(&(v_w.clone()))? { let mut v_x: Py; 
v_x = t3437;
out.push(py::truediv(&v_x.clone(), &v_tot.clone())?); } Py::list(out) })]))?);
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn f_leslie_projection(mut v_F: Py, mut v_S: Py, mut v_n0: Py, mut v_t: Py, mut v_transpose: Py) -> R {
    let mut v_L: Py = Py::None;
    let mut v__: Py = Py::None;
    let mut v_i: Py = Py::None;
    let mut v_k: Py = Py::None;
    let mut v_n: Py = Py::None;
    let mut v_s: Py = Py::None;
    let mut v_x: Py = Py::None;
    let _ = f__leslie_in(v_F.clone(), v_S.clone())?;
    let _ = f__need({ let t3501 = Py::Bool(matches!(v_n0.clone(), Py::List(_))); if !t3501.truthy() { t3501 } else { { let t3499 = py::len(&v_n0.clone())?; Py::Bool({ let t3500 = py::len(&v_F.clone())?; if py::eq(&t3499, &t3500) { true } else { false } }) } } }, Py::str("n0 needs one number per age class"))?;
    for t3502 in py::iter(&(v_n0.clone()))? {
        v_x = t3502;
        let _ = f__fin(Py::str("n0"), Py::tuple(vec![v_x.clone()]))?;
        let _ = f__need({ let t3503 = v_x.clone(); Py::Bool({ let t3504 = Py::Int(0i64); if py::ge(&t3503, &t3504)? { true } else { false } }) }, Py::str("population numbers must be zero or positive"))?;
    }
    let _ = f__need({ let t3508 = Py::Bool(matches!(v_t.clone(), Py::Int(_) | Py::Big(_) | Py::Bool(_))); if !t3508.truthy() { t3508 } else { { let t3505 = Py::Int(0i64); Py::Bool({ let t3506 = v_t.clone(); if py::le(&t3505, &t3506)? { { let t3507 = Py::Int(500i64); if py::le(&t3506, &t3507)? { true } else { false } } } else { false } }) } } }, Py::str("t must be a whole number of time steps (0 to 500)"))?;
    v_k = py::len(&v_F.clone())?;
    v_L = { let mut out = vec![]; for t3509 in py::range_iter(&[(v_k.clone()).i()?]) { let mut v__: Py; 
v__ = t3509;
out.push(py::mul(&Py::list(vec![Py::Float(0.0)]), &v_k.clone())?); } Py::list(out) };
    py::setitem(&v_L.clone(), &Py::Int(0i64), { let mut out = vec![]; for t3510 in py::iter(&(v_F.clone()))? { let mut v_f: Py; 
v_f = t3510;
out.push(py::float(&v_f.clone())?); } Py::list(out) })?;
    for t3511 in py::iter(&(py::enumerate(&v_S.clone(), 0)?))? {
        let t3512 = (t3511).items()?;
        if t3512.len() != 2 { return Err(PyErr::Value("unpack".into())); }
        v_i = t3512[0].clone();
        v_s = t3512[1].clone();
        py::setitem(&py::getitem(&v_L.clone(), &py::add(&v_i.clone(), &Py::Int(1i64))?)?, &v_i.clone(), py::float(&v_s.clone())?)?;
    }
    if (v_transpose.clone()).truthy() {
        v_L = { let mut out = vec![]; for t3513 in py::iter(&(py::zip(&{ let mut v: Vec<Py> = vec![]; v.extend((v_L.clone()).items()?); v })?))? { let mut v_r: Py; 
v_r = t3513;
out.push(Py::list((v_r.clone()).items()?)); } Py::list(out) };
    }
    v_n = { let mut out = vec![]; for t3514 in py::iter(&(v_n0.clone()))? { let mut v_x: Py; 
v_x = t3514;
out.push(py::float(&v_x.clone())?); } Py::list(out) };
    for t3515 in py::range_iter(&[(v_t.clone()).i()?]) {
        v__ = t3515;
        v_n = { let mut out = vec![]; for t3516 in py::range_iter(&[(v_k.clone()).i()?]) { let mut v_i: Py; 
v_i = t3516;
out.push(py::fsum(&{ let mut out = vec![]; for t3517 in py::range_iter(&[(v_k.clone()).i()?]) { let mut v_j: Py; 
v_j = t3517;
out.push(py::mul(&py::getitem(&py::getitem(&v_L.clone(), &v_i.clone())?, &v_j.clone())?, &py::getitem(&v_n.clone(), &v_j.clone())?)?); } Py::list(out) })?); } Py::list(out) };
    }
    return Ok(f__out(Py::dict(vec![(Py::str("n_t"), v_n.clone()), (Py::str("total"), py::fsum(&v_n.clone())?)]))?);
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn f_lineweaver_burk(mut v_S: Py, mut v_v: Py) -> R {
    let mut v_a: Py = Py::None;
    let mut v_b: Py = Py::None;
    let mut v_r2: Py = Py::None;
    let t3589 = (f__kin_data(v_S.clone(), v_v.clone(), Py::Int(2i64))?).items()?;
    if t3589.len() != 2 { return Err(PyErr::Value("unpack".into())); }
    v_S = t3589[0].clone();
    v_v = t3589[1].clone();
    let t3601 = (f__ols({ let mut out = vec![]; for t3599 in py::iter(&(v_S.clone()))? { let mut v_x: Py; 
v_x = t3599;
out.push(py::truediv(&Py::Int(1i64), &v_x.clone())?); } Py::list(out) }, { let mut out = vec![]; for t3600 in py::iter(&(v_v.clone()))? { let mut v_y: Py; 
v_y = t3600;
out.push(py::truediv(&Py::Int(1i64), &v_y.clone())?); } Py::list(out) })?).items()?;
    if t3601.len() != 3 { return Err(PyErr::Value("unpack".into())); }
    v_b = t3601[0].clone();
    v_a = t3601[1].clone();
    v_r2 = t3601[2].clone();
    let _ = f__need({ let t3606 = { let t3602 = v_a.clone(); Py::Bool({ let t3603 = Py::Int(0i64); if py::gt(&t3602, &t3603)? { true } else { false } }) }; if !t3606.truthy() { t3606 } else { { let t3604 = v_b.clone(); Py::Bool({ let t3605 = Py::Int(0i64); if py::gt(&t3604, &t3605)? { true } else { false } }) } } }, Py::str("the double-reciprocal line must have a positive slope and intercept"))?;
    return Ok(f__out(Py::dict(vec![(Py::str("Vmax"), py::truediv(&Py::Int(1i64), &v_a.clone())?), (Py::str("Km"), py::truediv(&v_b.clone(), &v_a.clone())?), (Py::str("r2"), v_r2.clone())]))?);
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn f_linkage_disequilibrium(mut v_n_AB: Py, mut v_n_Ab: Py, mut v_n_aB: Py, mut v_n_ab: Py) -> R {
    let mut v_D: Py = Py::None;
    let mut v_c: Py = Py::None;
    let mut v_dmax: Py = Py::None;
    let mut v_n: Py = Py::None;
    let mut v_pA: Py = Py::None;
    let mut v_pB: Py = Py::None;
    let mut v_v: Py = Py::None;
    v_c = Py::list(vec![v_n_AB.clone(), v_n_Ab.clone(), v_n_aB.clone(), v_n_ab.clone()]);
    for t3662 in py::iter(&(v_c.clone()))? {
        v_v = t3662;
        let _ = f__need({ let t3665 = Py::Bool(matches!(v_v.clone(), Py::Int(_) | Py::Big(_) | Py::Bool(_))); if !t3665.truthy() { t3665 } else { { let t3663 = v_v.clone(); Py::Bool({ let t3664 = Py::Int(0i64); if py::ge(&t3663, &t3664)? { true } else { false } }) } } }, Py::str("haplotype counts must be whole numbers >= 0"))?;
    }
    v_n = py::sum(&v_c.clone(), None)?;
    let _ = f__need({ let t3666 = v_n.clone(); Py::Bool({ let t3667 = Py::Int(1i64); if py::ge(&t3666, &t3667)? { true } else { false } }) }, Py::str("at least one haplotype is needed"))?;
    let t3668 = (Py::tuple(vec![py::fraction(&[py::add(&v_n_AB.clone(), &v_n_Ab.clone())?, v_n.clone()])?, py::fraction(&[py::add(&v_n_AB.clone(), &v_n_aB.clone())?, v_n.clone()])?])).items()?;
    if t3668.len() != 2 { return Err(PyErr::Value("unpack".into())); }
    v_pA = t3668[0].clone();
    v_pB = t3668[1].clone();
    let _ = f__need({ let t3675 = { let t3669 = Py::Int(0i64); Py::Bool({ let t3670 = v_pA.clone(); if py::lt(&t3669, &t3670)? { { let t3671 = Py::Int(1i64); if py::lt(&t3670, &t3671)? { true } else { false } } } else { false } }) }; if !t3675.truthy() { t3675 } else { { let t3672 = Py::Int(0i64); Py::Bool({ let t3673 = v_pB.clone(); if py::lt(&t3672, &t3673)? { { let t3674 = Py::Int(1i64); if py::lt(&t3673, &t3674)? { true } else { false } } } else { false } }) } } }, Py::str("both loci must be polymorphic"))?;
    v_D = py::sub(&py::fraction(&[v_n_AB.clone(), v_n.clone()])?, &py::mul(&v_pA.clone(), &v_pB.clone())?)?;
    if ({ let t3676 = v_D.clone(); Py::Bool({ let t3677 = Py::Int(0i64); if py::gt(&t3676, &t3677)? { true } else { false } }) }).truthy() {
        v_dmax = py::minmax(&[py::mul(&v_pA.clone(), &py::sub(&Py::Int(1i64), &v_pB.clone())?)?, py::mul(&py::sub(&Py::Int(1i64), &v_pA.clone())?, &v_pB.clone())?], false)?;
    } else {
        v_dmax = py::minmax(&[py::mul(&v_pA.clone(), &v_pB.clone())?, py::mul(&py::sub(&Py::Int(1i64), &v_pA.clone())?, &py::sub(&Py::Int(1i64), &v_pB.clone())?)?], false)?;
    }
    return Ok(Py::dict(vec![(Py::str("D"), py::float(&v_D.clone())?), (Py::str("D_prime"), (if ({ let t3678 = v_D.clone(); Py::Bool({ let t3679 = Py::Int(0i64); if py::ne(&t3678, &t3679) { true } else { false } }) }).truthy() { py::float(&py::truediv(&v_D.clone(), &v_dmax.clone())?)? } else { Py::Float(0.0) })), (Py::str("r2"), py::float(&py::truediv(&py::mul(&v_D.clone(), &v_D.clone())?, &py::mul(&py::mul(&py::mul(&v_pA.clone(), &py::sub(&Py::Int(1i64), &v_pA.clone())?)?, &v_pB.clone())?, &py::sub(&Py::Int(1i64), &v_pB.clone())?)?)?)?)]));
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn f_local_affine(mut v_seq1: Py, mut v_seq2: Py, mut v_match: Py, mut v_mismatch: Py, mut v_open: Py, mut v_extend: Py) -> R {
    let mut v_a: Py = Py::None;
    let mut v_b: Py = Py::None;
    let t3744 = (Py::tuple(vec![f__seq(v_seq1.clone(), Py::str("seq1"), Py::str("ANY"))?, f__seq(v_seq2.clone(), Py::str("seq2"), Py::str("ANY"))?])).items()?;
    if t3744.len() != 2 { return Err(PyErr::Value("unpack".into())); }
    v_a = t3744[0].clone();
    v_b = t3744[1].clone();
    let _ = f__check_align(v_a.clone(), v_b.clone(), Py::dict(vec![(Py::str("match"), v_match.clone()), (Py::str("mismatch"), v_mismatch.clone()), (Py::str("open"), v_open.clone()), (Py::str("extend"), v_extend.clone())]))?;
    let _ = f__need({ let t3756 = { let t3752 = v_open.clone(); Py::Bool({ let t3753 = Py::Int(0i64); if py::le(&t3752, &t3753)? { true } else { false } }) }; if !t3756.truthy() { t3756 } else { { let t3754 = v_extend.clone(); Py::Bool({ let t3755 = Py::Int(0i64); if py::le(&t3754, &t3755)? { true } else { false } }) } } }, Py::str("gap scores must be zero or negative"))?;
    let _ = f__need({ let t3757 = v_match.clone(); Py::Bool({ let t3758 = Py::Int(0i64); if py::gt(&t3757, &t3758)? { true } else { false } }) }, Py::str("the match score must be positive for a local alignment"))?;
    return Ok(f__out(Py::dict(vec![(Py::str("score"), f__align(v_a.clone(), v_b.clone(), f__mm(v_match.clone(), v_mismatch.clone())?, v_open.clone(), v_extend.clone(), Py::str("local"), Py::Bool(false))?)]))?);
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn f_local_linear(mut v_seq1: Py, mut v_seq2: Py, mut v_match: Py, mut v_mismatch: Py, mut v_gap: Py) -> R {
    let mut v_a: Py = Py::None;
    let mut v_b: Py = Py::None;
    let t3868 = (Py::tuple(vec![f__seq(v_seq1.clone(), Py::str("seq1"), Py::str("ANY"))?, f__seq(v_seq2.clone(), Py::str("seq2"), Py::str("ANY"))?])).items()?;
    if t3868.len() != 2 { return Err(PyErr::Value("unpack".into())); }
    v_a = t3868[0].clone();
    v_b = t3868[1].clone();
    let _ = f__check_align(v_a.clone(), v_b.clone(), Py::dict(vec![(Py::str("match"), v_match.clone()), (Py::str("mismatch"), v_mismatch.clone()), (Py::str("gap"), v_gap.clone())]))?;
    let _ = f__need({ let t3876 = v_gap.clone(); Py::Bool({ let t3877 = Py::Int(0i64); if py::le(&t3876, &t3877)? { true } else { false } }) }, Py::str("the gap score must be zero or negative"))?;
    let _ = f__need({ let t3878 = v_match.clone(); Py::Bool({ let t3879 = Py::Int(0i64); if py::gt(&t3878, &t3879)? { true } else { false } }) }, Py::str("the match score must be positive for a local alignment"))?;
    return Ok(f__out(Py::dict(vec![(Py::str("score"), f__align(v_a.clone(), v_b.clone(), f__mm(v_match.clone(), v_mismatch.clone())?, v_gap.clone(), v_gap.clone(), Py::str("local"), Py::Bool(false))?)]))?);
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn f_longest_orf(mut v_seq: Py) -> R {
    let mut v_L: Py = Py::None;
    let mut v_res: Py = Py::None;
    let mut v_s: Py = Py::None;
    v_s = f__seq(v_seq.clone(), Py::str("seq"), Py::str("DNA"))?;
    v_res = py::add(&f__orf_lengths(v_s.clone())?, &f__orf_lengths(f__revcomp(v_s.clone())?)?)?;
    v_L = (if (v_res.clone()).truthy() { py::minmax(&[v_res.clone()], true)? } else { Py::Int(0i64) });
    return Ok(Py::dict(vec![(Py::str("length_nt"), v_L.clone()), (Py::str("protein_length"), (if (v_L.clone()).truthy() { py::sub(&py::floordiv(&v_L.clone(), &Py::Int(3i64))?, &Py::Int(1i64))? } else { Py::Int(0i64) })), (Py::str("orf_count"), py::len(&v_res.clone())?)]));
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn f_map_distance(mut v_r: Py, mut v_map: Py) -> R {
    let mut v_m: Py = Py::None;
    let _ = f__fin(Py::str("r"), Py::tuple(vec![v_r.clone()]))?;
    let _ = f__need({ let t4041 = Py::Int(0i64); Py::Bool({ let t4042 = v_r.clone(); if py::le(&t4041, &t4042)? { { let t4043 = Py::Float(0.5); if py::lt(&t4042, &t4043)? { true } else { false } } } else { false } }) }, Py::str("the recombination fraction must lie in [0, 0.5)"))?;
    let _ = f__need({ let t4044 = v_map.clone(); Py::Bool({ let t4045 = Py::tuple(vec![Py::str("haldane"), Py::str("kosambi")]); if py::contains(&t4045, &t4044)? { true } else { false } }) }, Py::str("map must be haldane or kosambi"))?;
    v_m = (if ({ let t4046 = v_map.clone(); Py::Bool({ let t4047 = Py::str("haldane"); if py::eq(&t4046, &t4047) { true } else { false } }) }).truthy() { py::mul(&py::neg(&Py::Float(0.5))?, &py::math1("log1p", &py::mul(&py::neg(&Py::Int(2i64))?, &v_r.clone())?)?)? } else { py::mul(&Py::Float(0.25), &py::math1("log", &py::truediv(&py::add(&Py::Int(1i64), &py::mul(&Py::Int(2i64), &v_r.clone())?)?, &py::sub(&Py::Int(1i64), &py::mul(&Py::Int(2i64), &v_r.clone())?)?)?)?)? });
    return Ok(f__out(Py::dict(vec![(Py::str("d_cM"), py::mul(&Py::Int(100i64), &v_m.clone())?)]))?);
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn f_mm_fit(mut v_S: Py, mut v_v: Py) -> R {
    let mut v_F: Py = Py::None;
    let mut v_K: Py = Py::None;
    let mut v__: Py = Py::None;
    let mut v_grid: Py = Py::None;
    let mut v_hi: Py = Py::None;
    let mut v_i: Py = Py::None;
    let mut v_lo: Py = Py::None;
    let mut v_svx: Py = Py::None;
    let mut v_sxx: Py = Py::None;
    let t4119 = (f__kin_data(v_S.clone(), v_v.clone(), Py::Int(3i64))?).items()?;
    if t4119.len() != 2 { return Err(PyErr::Value("unpack".into())); }
    v_S = t4119[0].clone();
    v_v = t4119[1].clone();
    let _ = f__need({ let t4120 = py::len(&py::set(&v_S.clone())?)?; Py::Bool({ let t4121 = Py::Int(3i64); if py::ge(&t4120, &t4121)? { true } else { false } }) }, Py::str("at least 3 different substrate concentrations are needed"))?;
    let t4122 = (Py::tuple(vec![py::mul(&py::minmax(&[v_S.clone()], false)?, &Py::Float(0.001))?, py::mul(&py::minmax(&[v_S.clone()], true)?, &Py::Float(1000.0))?])).items()?;
    if t4122.len() != 2 { return Err(PyErr::Value("unpack".into())); }
    v_lo = t4122[0].clone();
    v_hi = t4122[1].clone();
    v_grid = { let mut out = vec![]; for t4123 in py::range_iter(&[(Py::Int(301i64)).i()?]) { let mut v_i: Py; 
v_i = t4123;
out.push(py::mul(&v_lo.clone(), &py::pow(&py::truediv(&v_hi.clone(), &v_lo.clone())?, &py::truediv(&v_i.clone(), &Py::Int(300i64))?)?)?); } Py::list(out) };
    v_F = Py::list(vec![]);
    for t4124 in py::iter(&(v_grid.clone()))? {
        v_K = t4124;
        let t4134 = (f__mm_varpro(v_S.clone(), v_v.clone(), v_K.clone())?).items()?;
        if t4134.len() != 3 { return Err(PyErr::Value("unpack".into())); }
        v_svx = t4134[0].clone();
        v_sxx = t4134[1].clone();
        v__ = t4134[2].clone();
        let _ = py::method(&v_F.clone(), "append", &[py::truediv(&py::mul(&v_svx.clone(), &v_svx.clone())?, &v_sxx.clone())?])?;
    }
    v_i = py::minmax_key(&[Py::list((py::range_iter(&[(py::len(&v_F.clone())?).i()?])).into_iter().collect::<Vec<Py>>())], &{ let v_F = v_F.clone(); let v_K = v_K.clone(); let v_S = v_S.clone(); let v__ = v__.clone(); let v_grid = v_grid.clone(); let v_hi = v_hi.clone(); let v_lo = v_lo.clone(); let v_svx = v_svx.clone(); let v_sxx = v_sxx.clone(); let v_v = v_v.clone(); py::func(move |a: &[Py]| -> R {
            #[allow(unused_mut, unused_variables, unused_assignments)]
{
            let mut v_k: Py = a.get(0).cloned().unwrap_or(Py::None);
            return Ok(py::getitem(&v_F.clone(), &v_k.clone())?);
            #[allow(unreachable_code)]
            Ok(Py::None)
        }
        }) }, true)?;
    let _ = f__need({ let t4135 = Py::Int(0i64); Py::Bool({ let t4136 = v_i.clone(); if py::lt(&t4135, &t4136)? { { let t4137 = py::sub(&py::len(&v_F.clone())?, &Py::Int(1i64))?; if py::lt(&t4136, &t4137)? { true } else { false } } } else { false } }) }, Py::str("the data do not determine Km (the best fit lies at the edge of the search range)"))?;
    v_K = f__bisect({ let v_F = v_F.clone(); let v_K = v_K.clone(); let v_S = v_S.clone(); let v__ = v__.clone(); let v_grid = v_grid.clone(); let v_hi = v_hi.clone(); let v_i = v_i.clone(); let v_lo = v_lo.clone(); let v_svx = v_svx.clone(); let v_sxx = v_sxx.clone(); let v_v = v_v.clone(); py::func(move |a: &[Py]| -> R {
            #[allow(unused_mut, unused_variables, unused_assignments)]
{
            let mut v_k: Py = a.get(0).cloned().unwrap_or(Py::None);
            return Ok(py::getitem(&f__mm_varpro(v_S.clone(), v_v.clone(), v_k.clone())?, &Py::Int(2i64))?);
            #[allow(unreachable_code)]
            Ok(Py::None)
        }
        }) }, py::getitem(&v_grid.clone(), &py::sub(&v_i.clone(), &Py::Int(1i64))?)?, py::getitem(&v_grid.clone(), &py::add(&v_i.clone(), &Py::Int(1i64))?)?, Py::Int(400i64))?;
    let t4164 = (f__mm_varpro(v_S.clone(), v_v.clone(), v_K.clone())?).items()?;
    if t4164.len() != 3 { return Err(PyErr::Value("unpack".into())); }
    v_svx = t4164[0].clone();
    v_sxx = t4164[1].clone();
    v__ = t4164[2].clone();
    return Ok(f__out(Py::dict(vec![(Py::str("Vmax"), py::truediv(&v_svx.clone(), &v_sxx.clone())?), (Py::str("Km"), v_K.clone())]))?);
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn f_mm_progress(mut v_S0: Py, mut v_Vmax: Py, mut v_Km: Py, mut v_t: Py, mut v_first_order: Py) -> R {
    let mut v_P: Py = Py::None;
    let mut v_Ph: Py = Py::None;
    let mut v_S: Py = Py::None;
    let mut v_fP: Py = Py::None;
    let mut v_fx: Py = Py::None;
    let mut v_vt: Py = Py::None;
    let mut v_x: Py = Py::None;
    let _ = f__mm_check(v_S0.clone(), v_Vmax.clone(), v_Km.clone())?;
    let _ = f__fin(Py::str("t"), Py::tuple(vec![v_t.clone()]))?;
    let _ = f__need({ let t4229 = v_t.clone(); Py::Bool({ let t4230 = Py::Int(0i64); if py::ge(&t4229, &t4230)? { true } else { false } }) }, Py::str("t must be zero or positive"))?;
    v_vt = py::mul(&v_Vmax.clone(), &v_t.clone())?;
    let _ = f__need({ let t4231 = py::truediv(&v_vt.clone(), &v_Km.clone())?; Py::Bool({ let t4232 = Py::Int(600i64); if py::le(&t4231, &t4232)? { true } else { false } }) }, Py::str("the substrate is exhausted beyond representable precision (Vmax t / Km > 600)"))?;
    if (v_first_order.clone()).truthy() {
        v_S = py::mul(&v_S0.clone(), &py::math1("exp", &py::truediv(&py::neg(&v_vt.clone())?, &v_Km.clone())?)?)?;
        return Ok(Py::dict(vec![(Py::str("S"), v_S.clone()), (Py::str("P"), py::sub(&v_S0.clone(), &v_S.clone())?)]));
    }
    if ({ let t4233 = v_vt.clone(); Py::Bool({ let t4234 = Py::Int(0i64); if py::eq(&t4233, &t4234) { true } else { false } }) }).truthy() {
        return Ok(Py::dict(vec![(Py::str("S"), py::float(&v_S0.clone())?), (Py::str("P"), Py::Float(0.0))]));
    }
    v_fP = { let v_Km = v_Km.clone(); let v_S = v_S.clone(); let v_S0 = v_S0.clone(); let v_Vmax = v_Vmax.clone(); let v_first_order = v_first_order.clone(); let v_t = v_t.clone(); let v_vt = v_vt.clone(); py::func(move |a: &[Py]| -> R {
            #[allow(unused_mut, unused_variables, unused_assignments)]
{
            let mut v_P: Py = a.get(0).cloned().unwrap_or(Py::None);
            return Ok(py::sub(&py::add(&py::mul(&py::neg(&v_Km.clone())?, &py::math1("log1p", &py::truediv(&py::neg(&v_P.clone())?, &v_S0.clone())?)?)?, &v_P.clone())?, &v_vt.clone())?);
            #[allow(unreachable_code)]
            Ok(Py::None)
        }
        }) };
    v_Ph = py::truediv(&v_S0.clone(), &Py::Int(2i64))?;
    if ({ let t4235 = py::call(&v_fP.clone(), &[v_Ph.clone()])?; Py::Bool({ let t4236 = Py::Int(0i64); if py::ge(&t4235, &t4236)? { true } else { false } }) }).truthy() {
        v_P = f__bisect(v_fP.clone(), Py::Float(0.0), v_Ph.clone(), Py::Int(400i64))?;
        return Ok(f__out(Py::dict(vec![(Py::str("S"), py::sub(&v_S0.clone(), &v_P.clone())?), (Py::str("P"), v_P.clone())]))?);
    }
    v_fx = { let v_Km = v_Km.clone(); let v_P = v_P.clone(); let v_Ph = v_Ph.clone(); let v_S = v_S.clone(); let v_S0 = v_S0.clone(); let v_Vmax = v_Vmax.clone(); let v_fP = v_fP.clone(); let v_first_order = v_first_order.clone(); let v_t = v_t.clone(); let v_vt = v_vt.clone(); py::func(move |a: &[Py]| -> R {
            #[allow(unused_mut, unused_variables, unused_assignments)]
{
            let mut v_x: Py = a.get(0).cloned().unwrap_or(Py::None);
            return Ok(py::sub(&py::sub(&py::add(&py::mul(&v_Km.clone(), &py::sub(&py::math1("log", &v_S0.clone())?, &v_x.clone())?)?, &v_S0.clone())?, &py::math1("exp", &v_x.clone())?)?, &v_vt.clone())?);
            #[allow(unreachable_code)]
            Ok(Py::None)
        }
        }) };
    v_x = f__bisect(v_fx.clone(), py::sub(&py::sub(&py::math1("log", &v_S0.clone())?, &py::truediv(&v_vt.clone(), &v_Km.clone())?)?, &Py::Float(1e-09))?, py::math1("log", &py::truediv(&v_S0.clone(), &Py::Int(2i64))?)?, Py::Int(400i64))?;
    v_S = py::math1("exp", &v_x.clone())?;
    return Ok(f__out(Py::dict(vec![(Py::str("S"), v_S.clone()), (Py::str("P"), py::sub(&v_S0.clone(), &v_S.clone())?)]))?);
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn f_mm_time(mut v_S0: Py, mut v_Vmax: Py, mut v_Km: Py, mut v_X: Py, mut v_zero_order: Py) -> R {
    let _ = f__mm_check(v_S0.clone(), v_Vmax.clone(), v_Km.clone())?;
    let _ = f__fin(Py::str("X"), Py::tuple(vec![v_X.clone()]))?;
    let _ = f__need({ let t4327 = Py::Int(0i64); Py::Bool({ let t4328 = v_X.clone(); if py::lt(&t4327, &t4328)? { { let t4329 = Py::Int(1i64); if py::lt(&t4328, &t4329)? { true } else { false } } } else { false } }) }, Py::str("the fraction converted must lie strictly between 0 and 1"))?;
    if (v_zero_order.clone()).truthy() {
        return Ok(Py::dict(vec![(Py::str("t"), py::truediv(&py::mul(&v_S0.clone(), &v_X.clone())?, &v_Vmax.clone())?)]));
    }
    return Ok(f__out(Py::dict(vec![(Py::str("t"), py::truediv(&py::add(&py::mul(&py::neg(&v_Km.clone())?, &py::math1("log1p", &py::neg(&v_X.clone())?)?)?, &py::mul(&v_S0.clone(), &v_X.clone())?)?, &v_Vmax.clone())?)]))?);
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn f_moran_fixation(mut v_N: Py, mut v_i: Py, mut v_r: Py) -> R {
    let mut v_x: Py = Py::None;
    let _ = f__need({ let t4388 = Py::Bool(matches!(v_N.clone(), Py::Int(_) | Py::Big(_) | Py::Bool(_))); if !t4388.truthy() { t4388 } else { { let t4385 = Py::Int(2i64); Py::Bool({ let t4386 = v_N.clone(); if py::le(&t4385, &t4386)? { { let t4387 = Py::Int(100000i64); if py::le(&t4386, &t4387)? { true } else { false } } } else { false } }) } } }, Py::str("N must be a whole number of individuals (2 to 100000)"))?;
    let _ = f__need({ let t4392 = Py::Bool(matches!(v_i.clone(), Py::Int(_) | Py::Big(_) | Py::Bool(_))); if !t4392.truthy() { t4392 } else { { let t4389 = Py::Int(0i64); Py::Bool({ let t4390 = v_i.clone(); if py::lt(&t4389, &t4390)? { { let t4391 = v_N.clone(); if py::lt(&t4390, &t4391)? { true } else { false } } } else { false } }) } } }, Py::str("i must be a whole number of mutants between 1 and N - 1"))?;
    let _ = f__pos(Py::str("r"), Py::tuple(vec![v_r.clone()]))?;
    if ({ let t4399 = v_r.clone(); Py::Bool({ let t4400 = Py::Int(1i64); if py::eq(&t4399, &t4400) { true } else { false } }) }).truthy() {
        return Ok(f__out(Py::dict(vec![(Py::str("rho"), py::truediv(&v_i.clone(), &v_N.clone())?)]))?);
    }
    v_x = py::math1("log", &v_r.clone())?;
    return Ok(f__out(Py::dict(vec![(Py::str("rho"), py::math1("exp", &py::sub(&f__log_abs_expm1(py::mul(&py::neg(&v_i.clone())?, &v_x.clone())?)?, &f__log_abs_expm1(py::mul(&py::neg(&v_N.clone())?, &v_x.clone())?)?)?)?)]))?);
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn f_mutation_selection(mut v_mu: Py, mut v_s: Py, mut v_h: Py) -> R {
    let mut v_f: Py = Py::None;
    let _ = f__fin(Py::str("mu"), Py::tuple(vec![v_mu.clone()]))?;
    let _ = f__need({ let t4463 = Py::Int(0i64); Py::Bool({ let t4464 = v_mu.clone(); if py::lt(&t4463, &t4464)? { { let t4465 = Py::Float(0.01); if py::le(&t4464, &t4465)? { true } else { false } } } else { false } }) }, Py::str("mu must lie in (0, 0.01]"))?;
    let _ = f__need({ let t4469 = Py::Bool(matches!(v_s.clone(), Py::Int(_) | Py::Big(_) | Py::Bool(_) | Py::Float(_))); if !t4469.truthy() { t4469 } else { { let t4466 = Py::Float(0.01); Py::Bool({ let t4467 = v_s.clone(); if py::le(&t4466, &t4467)? { { let t4468 = Py::Int(1i64); if py::le(&t4467, &t4468)? { true } else { false } } } else { false } }) } } }, Py::str("s must lie in [0.01, 1]"))?;
    let _ = f__need({ let t4473 = Py::Bool(matches!(v_h.clone(), Py::Int(_) | Py::Big(_) | Py::Bool(_) | Py::Float(_))); if !t4473.truthy() { t4473 } else { { let t4470 = Py::Int(0i64); Py::Bool({ let t4471 = v_h.clone(); if py::le(&t4470, &t4471)? { { let t4472 = Py::Int(1i64); if py::le(&t4471, &t4472)? { true } else { false } } } else { false } }) } } }, Py::str("h must lie in [0, 1]"))?;
    v_f = { let v_h = v_h.clone(); let v_mu = v_mu.clone(); let v_s = v_s.clone(); py::func(move |a: &[Py]| -> R {
            #[allow(unused_mut, unused_variables, unused_assignments)]
{
            let mut v_q: Py = a.get(0).cloned().unwrap_or(Py::None);
            return Ok(py::sub(&f__mutsel_map(v_q.clone(), v_mu.clone(), v_s.clone(), v_h.clone())?, &v_q.clone())?);
            #[allow(unreachable_code)]
            Ok(Py::None)
        }
        }) };
    let _ = f__need({ let t4474 = py::call(&v_f.clone(), &[Py::Float(0.5)])?; Py::Bool({ let t4475 = Py::Int(0i64); if py::lt(&t4474, &t4475)? { true } else { false } }) }, Py::str("mutation is too strong relative to selection: no low-frequency balance"))?;
    return Ok(f__out(Py::dict(vec![(Py::str("q_eq"), f__bisect(v_f.clone(), Py::Float(0.0), Py::Float(0.5), Py::Int(400i64))?)]))?);
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn f_ne_fluctuating(mut v_sizes: Py) -> R {
    let mut v_n: Py = Py::None;
    let mut v_ne: Py = Py::None;
    let mut v_ns: Py = Py::None;
    let mut v_ret: Py = Py::None;
    let _ = f__need({ let t4559 = Py::Bool(matches!(v_sizes.clone(), Py::List(_))); if !t4559.truthy() { t4559 } else { { let t4557 = py::len(&v_sizes.clone())?; Py::Bool({ let t4558 = Py::Int(2i64); if py::ge(&t4557, &t4558)? { true } else { false } }) } } }, Py::str("give the sizes of at least 2 generations"))?;
    v_ns = { let mut out = vec![]; for t4560 in py::iter(&(v_sizes.clone()))? { let mut v_v: Py; 
v_v = t4560;
out.push(f__whole(Py::str("sizes"), v_v.clone(), Py::Int(1i64))?); } Py::list(out) };
    v_ne = py::truediv(&py::fraction(&[py::len(&v_ns.clone())?])?, &py::sum(&{ let mut out = vec![]; for t4569 in py::iter(&(v_ns.clone()))? { let mut v_n: Py; 
v_n = t4569;
out.push(py::fraction(&[Py::Int(1i64), v_n.clone()])?); } Py::list(out) }, None)?)?;
    v_ret = Py::Float(1.0);
    for t4570 in py::iter(&(v_ns.clone()))? {
        v_n = t4570;
        v_ret = py::mul(&v_ret, &py::sub(&Py::Int(1i64), &py::truediv(&Py::Int(1i64), &py::mul(&Py::Int(2i64), &v_n.clone())?)?)?)?;
    }
    return Ok(f__out(Py::dict(vec![(Py::str("Ne"), py::float(&v_ne.clone())?), (Py::str("retained"), v_ret.clone())]))?);
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn f_ne_sex_ratio(mut v_N_m: Py, mut v_N_f: Py) -> R {
    let mut v_ne: Py = Py::None;
    let _ = f__need({ let t4632 = Py::Bool(matches!(v_N_m.clone(), Py::Int(_) | Py::Big(_) | Py::Bool(_))); if !t4632.truthy() { t4632 } else { { let t4631 = Py::Bool(matches!(v_N_f.clone(), Py::Int(_) | Py::Big(_) | Py::Bool(_))); if !t4631.truthy() { t4631 } else { { let t4630 = { let t4626 = v_N_m.clone(); Py::Bool({ let t4627 = Py::Int(1i64); if py::ge(&t4626, &t4627)? { true } else { false } }) }; if !t4630.truthy() { t4630 } else { { let t4628 = v_N_f.clone(); Py::Bool({ let t4629 = Py::Int(1i64); if py::ge(&t4628, &t4629)? { true } else { false } }) } } } } } } }, Py::str("N_m and N_f must be whole numbers of breeding males and females (at least 1)"))?;
    v_ne = py::fraction(&[py::mul(&py::mul(&Py::Int(4i64), &v_N_m.clone())?, &v_N_f.clone())?, py::add(&v_N_m.clone(), &v_N_f.clone())?])?;
    return Ok(Py::dict(vec![(Py::str("Ne"), py::float(&v_ne.clone())?), (Py::str("Ne_over_N"), py::float(&py::truediv(&v_ne.clone(), &py::add(&v_N_m.clone(), &v_N_f.clone())?)?)?)]));
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn f_neighbor_joining(mut v_D: Py, mut v_corrected: Py) -> R {
    let mut v_M: Py = Py::None;
    let mut v_Qv: Py = Py::None;
    let mut v_R: Py = Py::None;
    let mut v_a: Py = Py::None;
    let mut v_b: Py = Py::None;
    let mut v_best: Py = Py::None;
    let mut v_branch: Py = Py::None;
    let mut v_d: Py = Py::None;
    let mut v_edges: Py = Py::None;
    let mut v_harmless: Py = Py::None;
    let mut v_i: Py = Py::None;
    let mut v_j: Py = Py::None;
    let mut v_k: Py = Py::None;
    let mut v_li: Py = Py::None;
    let mut v_lj: Py = Py::None;
    let mut v_lx: Py = Py::None;
    let mut v_n: Py = Py::None;
    let mut v_nodes: Py = Py::None;
    let mut v_nxt: Py = Py::None;
    let mut v_r: Py = Py::None;
    let mut v_scale: Py = Py::None;
    let mut v_tied: Py = Py::None;
    let mut v_u: Py = Py::None;
    let mut v_x: Py = Py::None;
    v_M = f__distance_matrix(v_D.clone(), Py::str("D"), Py::Int(3i64), Py::Int(30i64))?;
    v_n = py::len(&v_M.clone())?;
    v_d = { let mut out = vec![]; for t4691 in py::range_iter(&[(v_n.clone()).i()?]) { let mut v_i: Py; 
v_i = t4691;
for t4692 in py::range_iter(&[(v_n.clone()).i()?]) { let mut v_j: Py; 
v_j = t4692;
if ({ let t4693 = v_i.clone(); Py::Bool({ let t4694 = v_j.clone(); if py::ne(&t4693, &t4694) { true } else { false } }) }).truthy() { out.push((Py::tuple(vec![v_i.clone(), v_j.clone()]), py::getitem(&py::getitem(&v_M.clone(), &v_i.clone())?, &v_j.clone())?)); }}} Py::dict(out) };
    v_nodes = Py::list((Py::list((py::range_iter(&[(v_n.clone()).i()?])).into_iter().collect::<Vec<Py>>())).items()?);
    v_branch = Py::dict(vec![]);
    v_edges = Py::list(vec![]);
    v_nxt = v_n.clone();
    v_scale = py::minmax(&[py::method(&v_d.clone(), "values", &[])?], true)?;
    while ({ let t4732 = py::len(&v_nodes.clone())?; Py::Bool({ let t4733 = Py::Int(2i64); if py::gt(&t4732, &t4733)? { true } else { false } }) }).truthy() {
        v_r = py::len(&v_nodes.clone())?;
        v_R = { let mut out = vec![]; for t4695 in py::iter(&(v_nodes.clone()))? { let mut v_i: Py; 
v_i = t4695;
out.push((v_i.clone(), py::sum(&{ let mut out = vec![]; for t4696 in py::iter(&(v_nodes.clone()))? { let mut v_k: Py; 
v_k = t4696;
if ({ let t4697 = v_k.clone(); Py::Bool({ let t4698 = v_i.clone(); if py::ne(&t4697, &t4698) { true } else { false } }) }).truthy() { out.push(py::getitem(&v_d.clone(), &Py::tuple(vec![v_i.clone(), v_k.clone()]))?); }} Py::list(out) }, None)?)); } Py::dict(out) };
        v_Qv = { let mut out = vec![]; for t4699 in py::iter(&(py::enumerate(&v_nodes.clone(), 0)?))? { let mut v_a: Py; let mut v_i: Py; 
let t4700 = (t4699).items()?;
if t4700.len() != 2 { return Err(PyErr::Value("unpack".into())); }
v_a = t4700[0].clone();
v_i = t4700[1].clone();
for t4701 in py::iter(&(py::slice(&v_nodes.clone(), Some((py::add(&v_a.clone(), &Py::Int(1i64))?).i()?), None, None)?))? { let mut v_j: Py; 
v_j = t4701;
out.push((Py::tuple(vec![v_i.clone(), v_j.clone()]), py::sub(&py::sub(&py::mul(&py::sub(&v_r.clone(), &Py::Int(2i64))?, &py::getitem(&v_d.clone(), &Py::tuple(vec![v_i.clone(), v_j.clone()]))?)?, &py::getitem(&v_R.clone(), &v_i.clone())?)?, &py::getitem(&v_R.clone(), &v_j.clone())?)?)); }} Py::dict(out) };
        v_best = py::minmax(&[py::method(&v_Qv.clone(), "values", &[])?], false)?;
        let t4706 = (py::minmax(&[{ let mut out = vec![]; for t4702 in py::iter(&(py::method(&v_Qv.clone(), "items", &[])?))? { let mut v_k: Py; let mut v_v: Py; 
let t4703 = (t4702).items()?;
if t4703.len() != 2 { return Err(PyErr::Value("unpack".into())); }
v_k = t4703[0].clone();
v_v = t4703[1].clone();
if ({ let t4704 = v_v.clone(); Py::Bool({ let t4705 = v_best.clone(); if py::eq(&t4704, &t4705) { true } else { false } }) }).truthy() { out.push(v_k.clone()); }} Py::list(out) }], false)?).items()?;
        if t4706.len() != 2 { return Err(PyErr::Value("unpack".into())); }
        v_i = t4706[0].clone();
        v_j = t4706[1].clone();
        v_tied = py::set(&{ let mut out = vec![]; for t4707 in py::iter(&(py::method(&v_Qv.clone(), "items", &[])?))? { let mut v_k: Py; let mut v_v: Py; 
let t4708 = (t4707).items()?;
if t4708.len() != 2 { return Err(PyErr::Value("unpack".into())); }
v_k = t4708[0].clone();
v_v = t4708[1].clone();
if ({ let t4709 = py::abs(&py::sub(&v_v.clone(), &v_best.clone())?)?; Py::Bool({ let t4710 = py::mul(&py::mul(&Py::Float(1e-12), &v_scale.clone())?, &v_r.clone())?; if py::le(&t4709, &t4710)? { true } else { false } }) }).truthy() { out.push(v_k.clone()); }} Py::list(out) })?;
        v_harmless = (if ({ let t4711 = v_r.clone(); Py::Bool({ let t4712 = Py::Int(4i64); if py::eq(&t4711, &t4712) { true } else { false } }) }).truthy() { py::set(&Py::list(vec![Py::tuple(vec![v_i.clone(), v_j.clone()]), Py::tuple((py::sorted(&{ let mut out = vec![]; for t4713 in py::iter(&(v_nodes.clone()))? { let mut v_k: Py; 
v_k = t4713;
if ({ let t4714 = v_k.clone(); Py::Bool({ let t4715 = Py::tuple(vec![v_i.clone(), v_j.clone()]); if !py::contains(&t4715, &t4714)? { true } else { false } }) }).truthy() { out.push(v_k.clone()); }} Py::list(out) }, (Py::Bool(false)).truthy())?).items()?)]))? } else { py::set(&Py::list(vec![Py::tuple(vec![v_i.clone(), v_j.clone()])]))? });
        let _ = f__need({ let t4720 = { let t4716 = v_r.clone(); Py::Bool({ let t4717 = Py::Int(3i64); if py::eq(&t4716, &t4717) { true } else { false } }) }; if t4720.truthy() { t4720 } else { { let t4718 = v_tied.clone(); Py::Bool({ let t4719 = v_harmless.clone(); if py::le(&t4718, &t4719)? { true } else { false } }) } } }, Py::str("two pairs tie in the neighbour-joining criterion: the tree is not unique"))?;
        if (v_corrected.clone()).truthy() {
            v_li = py::add(&py::truediv(&py::getitem(&v_d.clone(), &Py::tuple(vec![v_i.clone(), v_j.clone()]))?, &Py::Int(2i64))?, &py::truediv(&py::sub(&py::getitem(&v_R.clone(), &v_i.clone())?, &py::getitem(&v_R.clone(), &v_j.clone())?)?, &py::mul(&Py::Int(2i64), &py::sub(&v_r.clone(), &Py::Int(2i64))?)?)?)?;
        } else {
            v_li = py::truediv(&py::getitem(&v_d.clone(), &Py::tuple(vec![v_i.clone(), v_j.clone()]))?, &Py::Int(2i64))?;
        }
        v_lj = py::sub(&py::getitem(&v_d.clone(), &Py::tuple(vec![v_i.clone(), v_j.clone()]))?, &v_li.clone())?;
        for t4721 in py::iter(&(Py::tuple(vec![Py::tuple(vec![v_i.clone(), v_li.clone()]), Py::tuple(vec![v_j.clone(), v_lj.clone()])])))? {
            let t4724 = (t4721).items()?;
            if t4724.len() != 2 { return Err(PyErr::Value("unpack".into())); }
            v_x = t4724[0].clone();
            v_lx = t4724[1].clone();
            let _ = py::method(&v_edges.clone(), "append", &[v_lx.clone()])?;
            if ({ let t4722 = v_x.clone(); Py::Bool({ let t4723 = v_n.clone(); if py::lt(&t4722, &t4723)? { true } else { false } }) }).truthy() {
                py::setitem(&v_branch.clone(), &v_x.clone(), v_lx.clone())?;
            }
        }
        v_u = v_nxt.clone();
        v_nxt = py::iadd(&v_nxt, &Py::Int(1i64))?;
        for t4725 in py::iter(&(v_nodes.clone()))? {
            v_k = t4725;
            if ({ let t4726 = v_k.clone(); Py::Bool({ let t4727 = Py::tuple(vec![v_i.clone(), v_j.clone()]); if !py::contains(&t4727, &t4726)? { true } else { false } }) }).truthy() {
                let t4728 = py::truediv(&py::sub(&py::add(&py::getitem(&v_d.clone(), &Py::tuple(vec![v_i.clone(), v_k.clone()]))?, &py::getitem(&v_d.clone(), &Py::tuple(vec![v_j.clone(), v_k.clone()]))?)?, &py::getitem(&v_d.clone(), &Py::tuple(vec![v_i.clone(), v_j.clone()]))?)?, &Py::Int(2i64))?;
                py::setitem(&v_d.clone(), &Py::tuple(vec![v_u.clone(), v_k.clone()]), t4728.clone())?;
                py::setitem(&v_d.clone(), &Py::tuple(vec![v_k.clone(), v_u.clone()]), t4728.clone())?;
            }
        }
        v_nodes = py::add(&{ let mut out = vec![]; for t4729 in py::iter(&(v_nodes.clone()))? { let mut v_k: Py; 
v_k = t4729;
if ({ let t4730 = v_k.clone(); Py::Bool({ let t4731 = Py::tuple(vec![v_i.clone(), v_j.clone()]); if !py::contains(&t4731, &t4730)? { true } else { false } }) }).truthy() { out.push(v_k.clone()); }} Py::list(out) }, &Py::list(vec![v_u.clone()]))?;
    }
    let t4734 = (v_nodes.clone()).items()?;
    if t4734.len() != 2 { return Err(PyErr::Value("unpack".into())); }
    v_a = t4734[0].clone();
    v_b = t4734[1].clone();
    let _ = py::method(&v_edges.clone(), "append", &[py::getitem(&v_d.clone(), &Py::tuple(vec![v_a.clone(), v_b.clone()]))?])?;
    for t4735 in py::iter(&(Py::tuple(vec![v_a.clone(), v_b.clone()])))? {
        v_x = t4735;
        if ({ let t4736 = v_x.clone(); Py::Bool({ let t4737 = v_n.clone(); if py::lt(&t4736, &t4737)? { true } else { false } }) }).truthy() {
            py::setitem(&v_branch.clone(), &v_x.clone(), py::getitem(&v_d.clone(), &Py::tuple(vec![v_a.clone(), v_b.clone()]))?)?;
        }
    }
    return Ok(f__out(Py::dict(vec![(Py::str("total_length"), py::sum(&v_edges.clone(), None)?), (Py::str("branch"), { let mut out = vec![]; for t4752 in py::range_iter(&[(v_n.clone()).i()?]) { let mut v_i: Py; 
v_i = t4752;
out.push(py::getitem(&v_branch.clone(), &v_i.clone())?); } Py::list(out) })]))?);
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn f_net_charge(mut v_seq: Py, mut v_pH: Py) -> R {
    let mut v_s: Py = Py::None;
    v_s = f__seq(v_seq.clone(), Py::str("seq"), Py::str("PROTEIN"))?;
    let _ = f__fin(Py::str("pH"), Py::tuple(vec![v_pH.clone()]))?;
    let _ = f__need({ let t5493 = Py::Int(0i64); Py::Bool({ let t5494 = v_pH.clone(); if py::le(&t5493, &t5494)? { { let t5495 = Py::Int(14i64); if py::le(&t5494, &t5495)? { true } else { false } } } else { false } }) }, Py::str("pH must lie between 0 and 14"))?;
    return Ok(f__out(Py::dict(vec![(Py::str("charge"), f__charge(v_s.clone(), v_pH.clone(), Py::Bool(true))?)]))?);
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn f_nucleic_mass(mut v_seq: Py, mut v_molecule: Py, mut v_strands: Py) -> R {
    let mut v_s: Py = Py::None;
    let mut v_table: Py = Py::None;
    let mut v_w: Py = Py::None;
    let _ = f__need({ let t4794 = v_molecule.clone(); Py::Bool({ let t4795 = Py::tuple(vec![Py::str("DNA"), Py::str("RNA")]); if py::contains(&t4795, &t4794)? { true } else { false } }) }, Py::str("molecule must be DNA or RNA"))?;
    let _ = f__need({ let t4796 = v_strands.clone(); Py::Bool({ let t4797 = Py::tuple(vec![Py::str("single"), Py::str("double")]); if py::contains(&t4797, &t4796)? { true } else { false } }) }, Py::str("strands must be single or double"))?;
    v_s = f__seq(v_seq.clone(), Py::str("seq"), v_molecule.clone())?;
    v_table = (if ({ let t4821 = v_molecule.clone(); Py::Bool({ let t4822 = Py::str("DNA"); if py::eq(&t4821, &t4822) { true } else { false } }) }).truthy() { c_DNA_W()? } else { c_RNA_W()? });
    v_w = f__chain_mass(v_s.clone(), v_table.clone())?;
    if ({ let t4824 = v_strands.clone(); Py::Bool({ let t4825 = Py::str("double"); if py::eq(&t4824, &t4825) { true } else { false } }) }).truthy() {
        v_w = py::iadd(&v_w, &f__chain_mass(py::method(&v_s.clone(), "translate", &[(if ({ let t4826 = v_molecule.clone(); Py::Bool({ let t4827 = Py::str("DNA"); if py::eq(&t4826, &t4827) { true } else { false } }) }).truthy() { c__COMP()? } else { c__COMP_RNA()? })])?, v_table.clone())?)?;
    }
    return Ok(f__out(Py::dict(vec![(Py::str("mass"), v_w.clone())]))?);
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn f_nucleotide_diversity(mut v_seqs: Py) -> R {
    let mut v_k: Py = Py::None;
    let mut v_rows: Py = Py::None;
    v_rows = f__alignment(v_seqs.clone(), Py::Int(2i64))?;
    v_k = f__pairwise_k(v_rows.clone())?;
    return Ok(Py::dict(vec![(Py::str("k"), py::float(&v_k.clone())?), (Py::str("pi"), py::float(&py::truediv(&v_k.clone(), &py::len(&py::getitem(&v_rows.clone(), &Py::Int(0i64))?)?)?)?)]));
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn f_pedigree_inbreeding(mut v_pedigree: Py, mut v_individual: Py) -> R {
    let mut v_F: Py = Py::None;
    let mut v_half: Py = Py::None;
    let mut v_memo: Py = Py::None;
    let mut v_names: Py = Py::None;
    let mut v_parents: Py = Py::None;
    let mut v_up_paths: Py = Py::None;
    v_half = py::fraction(&[Py::Int(1i64), Py::Int(2i64)])?;
    v_parents = f__pedigree(v_pedigree.clone())?;
    v_names = py::bitor(&py::set(&v_parents.clone())?, &py::set(&{ let mut out = vec![]; for t4959 in py::iter(&(py::method(&v_parents.clone(), "values", &[])?))? { let mut v_ps: Py; 
v_ps = t4959;
for t4960 in py::iter(&(v_ps.clone()))? { let mut v_p: Py; 
v_p = t4960;
out.push(v_p.clone()); }} Py::list(out) })?)?;
    let _ = f__need({ let t4961 = v_individual.clone(); Py::Bool({ let t4962 = v_names.clone(); if py::contains(&t4962, &t4961)? { true } else { false } }) }, Py::str(&[py::format(&v_individual.clone(), "")?, " is not in the pedigree".to_string()].concat()))?;
    v_memo = Py::dict(vec![]);
    let rec4963 = py::rec_cell();
    v_up_paths = { let v_half = v_half.clone(); let v_individual = v_individual.clone(); let v_memo = v_memo.clone(); let v_names = v_names.clone(); let v_parents = v_parents.clone(); let v_pedigree = v_pedigree.clone(); let rec4963 = rec4963.clone(); py::func(move |a: &[Py]| -> R {
            #[allow(unused_mut, unused_variables, unused_assignments)]
{
            let mut v_x: Py = a.get(0).cloned().unwrap_or(Py::None);
            let v_up_paths: Py = py::rec_get(&rec4963)?;
            let mut v_out: Py = Py::None;
            let mut v_p: Py = Py::None;
            v_out = Py::list(vec![Py::list(vec![v_x.clone()])]);
            for t4964 in py::iter(&(py::method(&v_parents.clone(), "get", &[v_x.clone(), Py::list(vec![])])?))? {
                v_p = t4964;
                v_out = py::iadd(&v_out, &{ let mut out = vec![]; for t4965 in py::iter(&(py::call(&v_up_paths.clone(), &[v_p.clone()])?))? { let mut v_q: Py; 
v_q = t4965;
out.push(py::add(&Py::list(vec![v_x.clone()]), &v_q.clone())?); } Py::list(out) })?;
            }
            return Ok(v_out.clone());
            #[allow(unreachable_code)]
            Ok(Py::None)
        }
        }) };
    py::rec_set(&rec4963, &v_up_paths);
    let rec4966 = py::rec_cell();
    v_F = { let v_half = v_half.clone(); let v_individual = v_individual.clone(); let v_memo = v_memo.clone(); let v_names = v_names.clone(); let v_parents = v_parents.clone(); let v_pedigree = v_pedigree.clone(); let v_up_paths = v_up_paths.clone(); let rec4966 = rec4966.clone(); py::func(move |a: &[Py]| -> R {
            #[allow(unused_mut, unused_variables, unused_assignments)]
{
            let mut v_x: Py = a.get(0).cloned().unwrap_or(Py::None);
            let v_F: Py = py::rec_get(&rec4966)?;
            let mut v_a: Py = Py::None;
            let mut v_b: Py = Py::None;
            let mut v_d: Py = Py::None;
            let mut v_pd: Py = Py::None;
            let mut v_ps: Py = Py::None;
            let mut v_s: Py = Py::None;
            let mut v_tot: Py = Py::None;
            if ({ let t4967 = v_x.clone(); Py::Bool({ let t4968 = v_parents.clone(); if !py::contains(&t4968, &t4967)? { true } else { false } }) }).truthy() {
                return Ok(py::fraction(&[Py::Int(0i64)])?);
            }
            if ({ let t4969 = v_x.clone(); Py::Bool({ let t4970 = v_memo.clone(); if py::contains(&t4970, &t4969)? { true } else { false } }) }).truthy() {
                return Ok(py::getitem(&v_memo.clone(), &v_x.clone())?);
            }
            let t4971 = (py::getitem(&v_parents.clone(), &v_x.clone())?).items()?;
            if t4971.len() != 2 { return Err(PyErr::Value("unpack".into())); }
            v_s = t4971[0].clone();
            v_d = t4971[1].clone();
            v_tot = py::fraction(&[Py::Int(0i64)])?;
            let t4972 = (Py::tuple(vec![py::call(&v_up_paths.clone(), &[v_s.clone()])?, py::call(&v_up_paths.clone(), &[v_d.clone()])?])).items()?;
            if t4972.len() != 2 { return Err(PyErr::Value("unpack".into())); }
            v_ps = t4972[0].clone();
            v_pd = t4972[1].clone();
            for t4973 in py::iter(&(v_ps.clone()))? {
                v_a = t4973;
                for t4974 in py::iter(&(v_pd.clone()))? {
                    v_b = t4974;
                    if ({ let t4979 = { let t4975 = py::getitem(&v_a.clone(), &py::neg(&Py::Int(1i64))?)?; Py::Bool({ let t4976 = py::getitem(&v_b.clone(), &py::neg(&Py::Int(1i64))?)?; if py::eq(&t4975, &t4976) { true } else { false } }) }; if !t4979.truthy() { t4979 } else { { let t4977 = py::bitand(&py::set(&v_a.clone())?, &py::set(&v_b.clone())?)?; Py::Bool({ let t4978 = py::set(&Py::list(vec![py::getitem(&v_a.clone(), &py::neg(&Py::Int(1i64))?)?]))?; if py::eq(&t4977, &t4978) { true } else { false } }) } } }).truthy() {
                        v_tot = py::iadd(&v_tot, &py::mul(&py::pow(&v_half.clone(), &py::sub(&py::add(&py::len(&v_a.clone())?, &py::len(&v_b.clone())?)?, &Py::Int(1i64))?)?, &py::add(&Py::Int(1i64), &py::call(&v_F.clone(), &[py::getitem(&v_a.clone(), &py::neg(&Py::Int(1i64))?)?])?)?)?)?;
                    }
                }
            }
            py::setitem(&v_memo.clone(), &v_x.clone(), v_tot.clone())?;
            return Ok(v_tot.clone());
            #[allow(unreachable_code)]
            Ok(Py::None)
        }
        }) };
    py::rec_set(&rec4966, &v_F);
    return Ok(Py::dict(vec![(Py::str("F"), py::float(&py::call(&v_F.clone(), &[v_individual.clone()])?)?)]));
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn f_petersen_mle(mut v_M: Py, mut v_C: Py, mut v_R: Py, mut v_ceil: Py) -> R {
    let mut v_nm: Py = Py::None;
    let mut v_v: Py = Py::None;
    for t5021 in py::iter(&(Py::tuple(vec![Py::tuple(vec![Py::str("M"), v_M.clone()]), Py::tuple(vec![Py::str("C"), v_C.clone()]), Py::tuple(vec![Py::str("R"), v_R.clone()])])))? {
        let t5025 = (t5021).items()?;
        if t5025.len() != 2 { return Err(PyErr::Value("unpack".into())); }
        v_nm = t5025[0].clone();
        v_v = t5025[1].clone();
        let _ = f__need({ let t5024 = Py::Bool(matches!(v_v.clone(), Py::Int(_) | Py::Big(_) | Py::Bool(_))); if !t5024.truthy() { t5024 } else { { let t5022 = v_v.clone(); Py::Bool({ let t5023 = Py::Int(1i64); if py::ge(&t5022, &t5023)? { true } else { false } }) } } }, Py::str(&[py::format(&v_nm.clone(), "")?, " must be a whole number >= 1".to_string()].concat()))?;
    }
    let _ = f__need({ let t5026 = v_R.clone(); Py::Bool({ let t5027 = py::minmax(&[v_M.clone(), v_C.clone()], false)?; if py::le(&t5026, &t5027)? { true } else { false } }) }, Py::str("R cannot exceed M or C"))?;
    return Ok(Py::dict(vec![(Py::str("N_mle"), (if (v_ceil.clone()).truthy() { py::neg(&py::floordiv(&py::mul(&py::neg(&v_M.clone())?, &v_C.clone())?, &v_R.clone())?)? } else { py::floordiv(&py::mul(&v_M.clone(), &v_C.clone())?, &v_R.clone())? }))]));
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn f_protein_global(mut v_seq1: Py, mut v_seq2: Py, mut v_open: Py, mut v_extend: Py) -> R {
    let mut v_a: Py = Py::None;
    let mut v_b: Py = Py::None;
    let t5171 = (Py::tuple(vec![f__seq(v_seq1.clone(), Py::str("seq1"), Py::str("PROTEIN"))?, f__seq(v_seq2.clone(), Py::str("seq2"), Py::str("PROTEIN"))?])).items()?;
    if t5171.len() != 2 { return Err(PyErr::Value("unpack".into())); }
    v_a = t5171[0].clone();
    v_b = t5171[1].clone();
    let _ = f__check_align(v_a.clone(), v_b.clone(), Py::dict(vec![(Py::str("open"), v_open.clone()), (Py::str("extend"), v_extend.clone())]))?;
    let _ = f__need({ let t5183 = { let t5179 = v_open.clone(); Py::Bool({ let t5180 = Py::Int(0i64); if py::le(&t5179, &t5180)? { true } else { false } }) }; if !t5183.truthy() { t5183 } else { { let t5181 = v_extend.clone(); Py::Bool({ let t5182 = Py::Int(0i64); if py::le(&t5181, &t5182)? { true } else { false } }) } } }, Py::str("gap scores must be zero or negative"))?;
    return Ok(f__out(Py::dict(vec![(Py::str("score"), f__align(v_a.clone(), v_b.clone(), py::func(|a: &[Py]| -> R { f__blo(a.get(0).cloned().unwrap_or(Py::None), a.get(1).cloned().unwrap_or(Py::None)) }), v_open.clone(), v_extend.clone(), Py::str("global"), Py::Bool(false))?)]))?);
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn f_protein_local(mut v_seq1: Py, mut v_seq2: Py, mut v_open: Py, mut v_extend: Py) -> R {
    let mut v_a: Py = Py::None;
    let mut v_b: Py = Py::None;
    let t5291 = (Py::tuple(vec![f__seq(v_seq1.clone(), Py::str("seq1"), Py::str("PROTEIN"))?, f__seq(v_seq2.clone(), Py::str("seq2"), Py::str("PROTEIN"))?])).items()?;
    if t5291.len() != 2 { return Err(PyErr::Value("unpack".into())); }
    v_a = t5291[0].clone();
    v_b = t5291[1].clone();
    let _ = f__check_align(v_a.clone(), v_b.clone(), Py::dict(vec![(Py::str("open"), v_open.clone()), (Py::str("extend"), v_extend.clone())]))?;
    let _ = f__need({ let t5303 = { let t5299 = v_open.clone(); Py::Bool({ let t5300 = Py::Int(0i64); if py::le(&t5299, &t5300)? { true } else { false } }) }; if !t5303.truthy() { t5303 } else { { let t5301 = v_extend.clone(); Py::Bool({ let t5302 = Py::Int(0i64); if py::le(&t5301, &t5302)? { true } else { false } }) } } }, Py::str("gap scores must be zero or negative"))?;
    return Ok(f__out(Py::dict(vec![(Py::str("score"), f__align(v_a.clone(), v_b.clone(), py::func(|a: &[Py]| -> R { f__blo(a.get(0).cloned().unwrap_or(Py::None), a.get(1).cloned().unwrap_or(Py::None)) }), v_open.clone(), v_extend.clone(), Py::str("local"), Py::Bool(false))?)]))?);
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn f_protein_mass(mut v_seq: Py) -> R {
    let mut v_s: Py = Py::None;
    v_s = f__seq(v_seq.clone(), Py::str("seq"), Py::str("PROTEIN"))?;
    return Ok(f__out(Py::dict(vec![(Py::str("mass"), f__chain_mass(v_s.clone(), c_PROTEIN_W()?)?)]))?);
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn f_rarefaction(mut v_counts: Py, mut v_n: Py, mut v_with_replacement: Py) -> R {
    let mut v_N: Py = Py::None;
    let mut v_c: Py = Py::None;
    let mut v_e: Py = Py::None;
    let mut v_tot: Py = Py::None;
    v_c = f__abund(v_counts.clone(), Py::str("counts"))?;
    v_N = py::sum(&v_c.clone(), None)?;
    let _ = f__need({ let t5583 = Py::Bool(matches!(v_n.clone(), Py::Int(_) | Py::Big(_) | Py::Bool(_))); if !t5583.truthy() { t5583 } else { { let t5580 = Py::Int(1i64); Py::Bool({ let t5581 = v_n.clone(); if py::le(&t5580, &t5581)? { { let t5582 = v_N.clone(); if py::le(&t5581, &t5582)? { true } else { false } } } else { false } }) } } }, Py::str("n must be a whole number from 1 to the number of individuals"))?;
    let _ = f__need({ let t5584 = v_N.clone(); Py::Bool({ let t5585 = Py::Int(100000i64); if py::le(&t5584, &t5585)? { true } else { false } }) }, Py::str("at most 100000 individuals"))?;
    if (v_with_replacement.clone()).truthy() {
        return Ok(Py::dict(vec![(Py::str("expected_species"), py::fsum(&{ let mut out = vec![]; for t5586 in py::iter(&(v_c.clone()))? { let mut v_x: Py; 
v_x = t5586;
out.push(py::sub(&Py::Int(1i64), &py::pow(&py::sub(&Py::Int(1i64), &py::truediv(&v_x.clone(), &v_N.clone())?)?, &v_n.clone())?)?); } Py::list(out) })?)]));
    }
    v_tot = py::comb(&v_N.clone(), &v_n.clone())?;
    v_e = py::sum(&{ let mut out = vec![]; for t5587 in py::iter(&(v_c.clone()))? { let mut v_x: Py; 
v_x = t5587;
out.push(py::sub(&Py::Int(1i64), &py::fraction(&[py::comb(&py::sub(&v_N.clone(), &v_x.clone())?, &v_n.clone())?, v_tot.clone()])?)?); } Py::list(out) }, None)?;
    return Ok(Py::dict(vec![(Py::str("expected_species"), py::float(&v_e.clone())?)]));
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn f_recessive_lethal(mut v_q0: Py, mut v_t: Py) -> R {
    let mut v__: Py = Py::None;
    let mut v_q: Py = Py::None;
    let _ = f__frac(Py::str("q0"), v_q0.clone(), Py::Bool(true), Py::Bool(false))?;
    let _ = f__need({ let t5644 = Py::Bool(matches!(v_t.clone(), Py::Int(_) | Py::Big(_) | Py::Bool(_))); if !t5644.truthy() { t5644 } else { { let t5641 = Py::Int(0i64); Py::Bool({ let t5642 = v_t.clone(); if py::le(&t5641, &t5642)? { { let t5643 = Py::Int(100000i64); if py::le(&t5642, &t5643)? { true } else { false } } } else { false } }) } } }, Py::str("t must be a whole number of generations, 0 to 100000"))?;
    v_q = v_q0.clone();
    for t5645 in py::range_iter(&[(v_t.clone()).i()?]) {
        v__ = t5645;
        v_q = py::truediv(&v_q.clone(), &py::add(&Py::Int(1i64), &v_q.clone())?)?;
    }
    return Ok(f__out(Py::dict(vec![(Py::str("q_t"), v_q.clone())]))?);
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn f_recombination_fraction(mut v_d_cM: Py, mut v_map: Py) -> R {
    let mut v_m: Py = Py::None;
    let _ = f__fin(Py::str("d_cM"), Py::tuple(vec![v_d_cM.clone()]))?;
    let _ = f__need({ let t5704 = v_d_cM.clone(); Py::Bool({ let t5705 = Py::Int(0i64); if py::ge(&t5704, &t5705)? { true } else { false } }) }, Py::str("the map distance must be zero or positive"))?;
    let _ = f__need({ let t5706 = v_map.clone(); Py::Bool({ let t5707 = Py::tuple(vec![Py::str("haldane"), Py::str("kosambi")]); if py::contains(&t5707, &t5706)? { true } else { false } }) }, Py::str("map must be haldane or kosambi"))?;
    v_m = py::truediv(&v_d_cM.clone(), &Py::Int(100i64))?;
    return Ok(f__out(Py::dict(vec![(Py::str("r"), (if ({ let t5722 = v_map.clone(); Py::Bool({ let t5723 = Py::str("haldane"); if py::eq(&t5722, &t5723) { true } else { false } }) }).truthy() { py::mul(&py::neg(&Py::Float(0.5))?, &py::math1("expm1", &py::mul(&py::neg(&Py::Int(2i64))?, &v_m.clone())?)?)? } else { py::mul(&Py::Float(0.5), &py::math1("tanh", &py::mul(&Py::Int(2i64), &v_m.clone())?)?)? }))]))?);
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn f_restriction_sites(mut v_seq: Py, mut v_enzyme: Py) -> R {
    let mut v_bounds: Py = Py::None;
    let mut v_cut: Py = Py::None;
    let mut v_pos: Py = Py::None;
    let mut v_s: Py = Py::None;
    let mut v_site: Py = Py::None;
    v_s = f__seq(v_seq.clone(), Py::str("seq"), Py::str("DNA"))?;
    let _ = f__need({ let t5788 = v_enzyme.clone(); Py::Bool({ let t5789 = c_ENZYMES()?; if py::contains(&t5789, &t5788)? { true } else { false } }) }, py::add(&Py::str("enzyme must be one of "), &py::method(&Py::str(", "), "join", &[c_ENZYMES()?])?)?)?;
    let t5790 = (py::getitem(&c_ENZYMES()?, &v_enzyme.clone())?).items()?;
    if t5790.len() != 2 { return Err(PyErr::Value("unpack".into())); }
    v_site = t5790[0].clone();
    v_cut = t5790[1].clone();
    v_pos = { let mut out = vec![]; for t5791 in py::range_iter(&[(py::add(&py::sub(&py::len(&v_s.clone())?, &py::len(&v_site.clone())?)?, &Py::Int(1i64))?).i()?]) { let mut v_i: Py; 
v_i = t5791;
if (py::method(&v_s.clone(), "startswith", &[v_site.clone(), v_i.clone()])?).truthy() { out.push(py::add(&py::add(&v_i.clone(), &v_cut.clone())?, &Py::Int(1i64))?); }} Py::list(out) };
    v_bounds = py::add(&py::add(&Py::list(vec![Py::Int(1i64)]), &v_pos.clone())?, &Py::list(vec![py::add(&py::len(&v_s.clone())?, &Py::Int(1i64))?]))?;
    return Ok(Py::dict(vec![(Py::str("n_sites"), py::len(&v_pos.clone())?), (Py::str("positions"), v_pos.clone()), (Py::str("fragments"), { let mut out = vec![]; for t5792 in py::range_iter(&[(py::sub(&py::len(&v_bounds.clone())?, &Py::Int(1i64))?).i()?]) { let mut v_i: Py; 
v_i = t5792;
out.push(py::sub(&py::getitem(&v_bounds.clone(), &py::add(&v_i.clone(), &Py::Int(1i64))?)?, &py::getitem(&v_bounds.clone(), &v_i.clone())?)?); } Py::list(out) })]));
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn f_reverse_complement(mut v_seq: Py) -> R {
    let mut v_s: Py = Py::None;
    v_s = f__seq(v_seq.clone(), Py::str("seq"), Py::str("DNA"))?;
    return Ok(Py::dict(vec![(Py::str("revcomp"), f__revcomp(v_s.clone())?)]));
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn f_richness_indices(mut v_counts: Py, mut v_log10: Py) -> R {
    let mut v_N: Py = Py::None;
    let mut v_S: Py = Py::None;
    let mut v_c: Py = Py::None;
    v_c = f__abund(v_counts.clone(), Py::str("counts"))?;
    let t5916 = (Py::tuple(vec![py::sum(&v_c.clone(), None)?, py::len(&v_c.clone())?])).items()?;
    if t5916.len() != 2 { return Err(PyErr::Value("unpack".into())); }
    v_N = t5916[0].clone();
    v_S = t5916[1].clone();
    let _ = f__need({ let t5917 = v_N.clone(); Py::Bool({ let t5918 = Py::Int(2i64); if py::ge(&t5917, &t5918)? { true } else { false } }) }, Py::str("at least 2 individuals are needed (ln N must be positive)"))?;
    return Ok(f__out(Py::dict(vec![(Py::str("margalef"), py::truediv(&py::sub(&v_S.clone(), &Py::Int(1i64))?, &(if (v_log10.clone()).truthy() { py::math1("log10", &v_N.clone())? } else { py::math1("log", &v_N.clone())? }))?), (Py::str("menhinick"), py::truediv(&v_S.clone(), &py::math1("sqrt", &v_N.clone())?)?)]))?);
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn f_robinson_foulds(mut v_tree1: Py, mut v_tree2: Py, mut v_one_way: Py) -> R {
    let mut v_l1: Py = Py::None;
    let mut v_l2: Py = Py::None;
    let mut v_r1: Py = Py::None;
    let mut v_r2: Py = Py::None;
    let mut v_rf: Py = Py::None;
    let mut v_s1: Py = Py::None;
    let mut v_s2: Py = Py::None;
    let t6018 = (f__parse_newick(v_tree1.clone())?).items()?;
    if t6018.len() != 2 { return Err(PyErr::Value("unpack".into())); }
    v_r1 = t6018[0].clone();
    v_l1 = t6018[1].clone();
    let t6019 = (f__parse_newick(v_tree2.clone())?).items()?;
    if t6019.len() != 2 { return Err(PyErr::Value("unpack".into())); }
    v_r2 = t6019[0].clone();
    v_l2 = t6019[1].clone();
    let _ = f__need({ let t6020 = py::set(&v_l1.clone())?; Py::Bool({ let t6021 = py::set(&v_l2.clone())?; if py::eq(&t6020, &t6021) { true } else { false } }) }, Py::str("both trees must have the same taxa"))?;
    let _ = f__need({ let t6022 = py::len(&v_l1.clone())?; Py::Bool({ let t6023 = Py::Int(4i64); if py::ge(&t6022, &t6023)? { true } else { false } }) }, Py::str("at least 4 taxa are needed"))?;
    let t6031 = (Py::tuple(vec![f__splits(v_r1.clone(), v_l1.clone())?, f__splits(v_r2.clone(), v_l2.clone())?])).items()?;
    if t6031.len() != 2 { return Err(PyErr::Value("unpack".into())); }
    v_s1 = t6031[0].clone();
    v_s2 = t6031[1].clone();
    let _ = f__need({ let t6032 = v_s1.clone(); if t6032.truthy() { t6032 } else { v_s2.clone() } }, Py::str("both trees are unresolved stars"))?;
    v_rf = py::add(&py::len(&py::sub(&v_s1.clone(), &v_s2.clone())?)?, &(if (v_one_way.clone()).truthy() { Py::Int(0i64) } else { py::len(&py::sub(&v_s2.clone(), &v_s1.clone())?)? }))?;
    return Ok(Py::dict(vec![(Py::str("rf"), v_rf.clone()), (Py::str("rf_normalized"), py::truediv(&v_rf.clone(), &py::add(&py::len(&v_s1.clone())?, &py::len(&v_s2.clone())?)?)?)]));
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn f_schumacher(mut v_C: Py, mut v_R: Py, mut v_M: Py, mut v_schnabel: Py) -> R {
    let mut v_c: Py = Py::None;
    let mut v_m: Py = Py::None;
    let mut v_r: Py = Py::None;
    let mut v_srm: Py = Py::None;
    let _ = f__need({ let t6080 = Py::Bool(matches!(v_C.clone(), Py::List(_))); if !t6080.truthy() { t6080 } else { { let t6079 = Py::Bool(matches!(v_R.clone(), Py::List(_))); if !t6079.truthy() { t6079 } else { { let t6078 = Py::Bool(matches!(v_M.clone(), Py::List(_))); if !t6078.truthy() { t6078 } else { { let t6074 = py::len(&v_C.clone())?; Py::Bool({ let t6075 = py::len(&v_R.clone())?; if py::eq(&t6074, &t6075) { { let t6076 = py::len(&v_M.clone())?; if py::eq(&t6075, &t6076) { { let t6077 = Py::Int(2i64); if py::ge(&t6076, &t6077)? { true } else { false } } } else { false } } } else { false } }) } } } } } } }, Py::str("C, R and M need one value per sampling occasion (at least 2)"))?;
    for t6081 in py::iter(&(py::zip(&[v_C.clone(), v_R.clone(), v_M.clone()])?))? {
        let t6098 = (t6081).items()?;
        if t6098.len() != 3 { return Err(PyErr::Value("unpack".into())); }
        v_c = t6098[0].clone();
        v_r = t6098[1].clone();
        v_m = t6098[2].clone();
        let _ = f__need({ let t6086 = { let t6082 = py::float(&v_c.clone())?; Py::Bool({ let t6083 = py::int(&v_c.clone())?; if py::eq(&t6082, &t6083) { true } else { false } }) }; if !t6086.truthy() { t6086 } else { { let t6084 = v_c.clone(); Py::Bool({ let t6085 = Py::Int(1i64); if py::ge(&t6084, &t6085)? { true } else { false } }) } } }, Py::str("catches must be whole numbers >= 1"))?;
        let _ = f__need({ let t6092 = { let t6087 = py::float(&v_r.clone())?; Py::Bool({ let t6088 = py::int(&v_r.clone())?; if py::eq(&t6087, &t6088) { true } else { false } }) }; if !t6092.truthy() { t6092 } else { { let t6089 = Py::Int(0i64); Py::Bool({ let t6090 = v_r.clone(); if py::le(&t6089, &t6090)? { { let t6091 = v_c.clone(); if py::le(&t6090, &t6091)? { true } else { false } } } else { false } }) } } }, Py::str("recaptures must be whole numbers from 0 to the catch"))?;
        let _ = f__need({ let t6097 = { let t6093 = py::float(&v_m.clone())?; Py::Bool({ let t6094 = py::int(&v_m.clone())?; if py::eq(&t6093, &t6094) { true } else { false } }) }; if !t6097.truthy() { t6097 } else { { let t6095 = v_m.clone(); Py::Bool({ let t6096 = Py::Int(0i64); if py::ge(&t6095, &t6096)? { true } else { false } }) } } }, Py::str("marked animals at large must be whole numbers >= 0"))?;
    }
    v_srm = py::sum(&{ let mut out = vec![]; for t6099 in py::iter(&(py::zip(&[v_R.clone(), v_M.clone()])?))? { let mut v_m: Py; let mut v_r: Py; 
let t6100 = (t6099).items()?;
if t6100.len() != 2 { return Err(PyErr::Value("unpack".into())); }
v_r = t6100[0].clone();
v_m = t6100[1].clone();
out.push(py::fraction(&[py::mul(&py::int(&v_r.clone())?, &py::int(&v_m.clone())?)?])?); } Py::list(out) }, None)?;
    let _ = f__need({ let t6101 = v_srm.clone(); Py::Bool({ let t6102 = Py::Int(0i64); if py::gt(&t6101, &t6102)? { true } else { false } }) }, Py::str("no marked animal was ever recaptured: the population cannot be estimated"))?;
    if (v_schnabel.clone()).truthy() {
        return Ok(Py::dict(vec![(Py::str("N"), py::float(&py::truediv(&py::sum(&{ let mut out = vec![]; for t6103 in py::iter(&(py::zip(&[v_C.clone(), v_M.clone()])?))? { let mut v_c: Py; let mut v_m: Py; 
let t6104 = (t6103).items()?;
if t6104.len() != 2 { return Err(PyErr::Value("unpack".into())); }
v_c = t6104[0].clone();
v_m = t6104[1].clone();
out.push(py::fraction(&[py::mul(&py::int(&v_c.clone())?, &py::int(&v_m.clone())?)?])?); } Py::list(out) }, None)?, &py::sum(&{ let mut out = vec![]; for t6105 in py::iter(&(v_R.clone()))? { let mut v_r: Py; 
v_r = t6105;
out.push(py::fraction(&[py::int(&v_r.clone())?])?); } Py::list(out) }, None)?)?)?)]));
    }
    return Ok(Py::dict(vec![(Py::str("N"), py::float(&py::truediv(&py::sum(&{ let mut out = vec![]; for t6106 in py::iter(&(py::zip(&[v_C.clone(), v_M.clone()])?))? { let mut v_c: Py; let mut v_m: Py; 
let t6107 = (t6106).items()?;
if t6107.len() != 2 { return Err(PyErr::Value("unpack".into())); }
v_c = t6107[0].clone();
v_m = t6107[1].clone();
out.push(py::fraction(&[py::mul(&py::int(&v_c.clone())?, &py::pow(&py::int(&v_m.clone())?, &Py::Int(2i64))?)?])?); } Py::list(out) }, None)?, &v_srm.clone())?)?)]));
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn f_semiglobal(mut v_seq1: Py, mut v_seq2: Py, mut v_match: Py, mut v_mismatch: Py, mut v_gap: Py) -> R {
    let mut v_a: Py = Py::None;
    let mut v_b: Py = Py::None;
    let t6239 = (Py::tuple(vec![f__seq(v_seq1.clone(), Py::str("seq1"), Py::str("ANY"))?, f__seq(v_seq2.clone(), Py::str("seq2"), Py::str("ANY"))?])).items()?;
    if t6239.len() != 2 { return Err(PyErr::Value("unpack".into())); }
    v_a = t6239[0].clone();
    v_b = t6239[1].clone();
    let _ = f__check_align(v_a.clone(), v_b.clone(), Py::dict(vec![(Py::str("match"), v_match.clone()), (Py::str("mismatch"), v_mismatch.clone()), (Py::str("gap"), v_gap.clone())]))?;
    let _ = f__need({ let t6247 = v_gap.clone(); Py::Bool({ let t6248 = Py::Int(0i64); if py::le(&t6247, &t6248)? { true } else { false } }) }, Py::str("the gap score must be zero or negative"))?;
    return Ok(f__out(Py::dict(vec![(Py::str("score"), f__align(v_a.clone(), v_b.clone(), f__mm(v_match.clone(), v_mismatch.clone())?, v_gap.clone(), v_gap.clone(), Py::str("semiglobal"), Py::Bool(false))?)]))?);
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn f_shannon(mut v_counts: Py, mut v_base: Py, mut v_wrong: Py) -> R {
    let mut v_H: Py = Py::None;
    let mut v_Hb: Py = Py::None;
    let mut v_J: Py = Py::None;
    let mut v_N: Py = Py::None;
    let mut v_S: Py = Py::None;
    let mut v_c: Py = Py::None;
    v_c = f__abund(v_counts.clone(), Py::str("counts"))?;
    let _ = f__need({ let t6353 = v_base.clone(); Py::Bool({ let t6354 = c__LOGS()?; if py::contains(&t6354, &t6353)? { true } else { false } }) }, Py::str("base must be e, 2 or 10"))?;
    v_N = py::sum(&v_c.clone(), None)?;
    v_H = py::neg(&py::fsum(&{ let mut out = vec![]; for t6355 in py::iter(&(v_c.clone()))? { let mut v_x: Py; 
v_x = t6355;
out.push(py::mul(&py::truediv(&v_x.clone(), &v_N.clone())?, &py::math1("log", &py::truediv(&v_x.clone(), &v_N.clone())?)?)?); } Py::list(out) })?)?;
    v_S = py::len(&v_c.clone())?;
    v_Hb = py::truediv(&v_H.clone(), &py::math1("log", &py::getitem(&c__LOGS()?, &v_base.clone())?)?)?;
    v_J = (if ({ let t6356 = v_S.clone(); Py::Bool({ let t6357 = Py::Int(1i64); if py::eq(&t6356, &t6357) { true } else { false } }) }).truthy() { Py::Float(1.0) } else { (if (v_wrong.clone()).truthy() { py::truediv(&v_H.clone(), &v_S.clone())? } else { py::truediv(&v_H.clone(), &py::math1("log", &v_S.clone())?)? }) });
    return Ok(f__out(Py::dict(vec![(Py::str("H"), v_Hb.clone()), (Py::str("evenness"), v_J.clone()), (Py::str("effective_number"), py::math1("exp", &v_H.clone())?)]))?);
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn f_simpson(mut v_counts: Py, mut v_estimator: Py, mut v_swap: Py) -> R {
    let mut v_D: Py = Py::None;
    let mut v_N: Py = Py::None;
    let mut v_c: Py = Py::None;
    let mut v_est: Py = Py::None;
    v_c = f__abund(v_counts.clone(), Py::str("counts"))?;
    let _ = f__need({ let t6431 = v_estimator.clone(); Py::Bool({ let t6432 = Py::tuple(vec![Py::str("infinite"), Py::str("finite")]); if py::contains(&t6432, &t6431)? { true } else { false } }) }, Py::str("estimator must be infinite or finite"))?;
    v_est = (if (v_swap.clone()).truthy() { py::getitem(&Py::dict(vec![(Py::str("infinite"), Py::str("finite")), (Py::str("finite"), Py::str("infinite"))]), &v_estimator.clone())? } else { v_estimator.clone() });
    v_N = py::sum(&v_c.clone(), None)?;
    if ({ let t6433 = v_est.clone(); Py::Bool({ let t6434 = Py::str("finite"); if py::eq(&t6433, &t6434) { true } else { false } }) }).truthy() {
        let _ = f__need({ let t6435 = v_N.clone(); Py::Bool({ let t6436 = Py::Int(2i64); if py::ge(&t6435, &t6436)? { true } else { false } }) }, Py::str("the finite-sample estimator needs at least 2 individuals"))?;
        v_D = py::fraction(&[py::sum(&{ let mut out = vec![]; for t6437 in py::iter(&(v_c.clone()))? { let mut v_x: Py; 
v_x = t6437;
out.push(py::mul(&v_x.clone(), &py::sub(&v_x.clone(), &Py::Int(1i64))?)?); } Py::list(out) }, None)?, py::mul(&v_N.clone(), &py::sub(&v_N.clone(), &Py::Int(1i64))?)?])?;
        let _ = f__need({ let t6438 = v_D.clone(); Py::Bool({ let t6439 = Py::Int(0i64); if py::gt(&t6438, &t6439)? { true } else { false } }) }, Py::str("every species is a singleton: the finite-sample dominance is 0 and its inverse undefined"))?;
    } else {
        v_D = py::fraction(&[py::sum(&{ let mut out = vec![]; for t6440 in py::iter(&(v_c.clone()))? { let mut v_x: Py; 
v_x = t6440;
out.push(py::mul(&v_x.clone(), &v_x.clone())?); } Py::list(out) }, None)?, py::mul(&v_N.clone(), &v_N.clone())?])?;
    }
    return Ok(Py::dict(vec![(Py::str("D"), py::float(&v_D.clone())?), (Py::str("simpson"), py::float(&py::sub(&Py::Int(1i64), &v_D.clone())?)?), (Py::str("inverse"), py::float(&py::truediv(&Py::Int(1i64), &v_D.clone())?)?)]));
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn f_species_area(mut v_A: Py, mut v_S: Py, mut v_no_exp: Py) -> R {
    let mut v_a: Py = Py::None;
    let mut v_b: Py = Py::None;
    let mut v_r2: Py = Py::None;
    let mut v_v: Py = Py::None;
    let _ = f__need({ let t6486 = Py::Bool(matches!(v_A.clone(), Py::List(_))); if !t6486.truthy() { t6486 } else { { let t6485 = Py::Bool(matches!(v_S.clone(), Py::List(_))); if !t6485.truthy() { t6485 } else { { let t6482 = py::len(&v_A.clone())?; Py::Bool({ let t6483 = py::len(&v_S.clone())?; if py::eq(&t6482, &t6483) { { let t6484 = Py::Int(2i64); if py::ge(&t6483, &t6484)? { true } else { false } } } else { false } }) } } } } }, Py::str("A and S must be lists of equal length"))?;
    for t6487 in py::iter(&(py::add(&v_A.clone(), &v_S.clone())?))? {
        v_v = t6487;
        let _ = f__pos(Py::str("areas and species counts"), Py::tuple(vec![v_v.clone()]))?;
    }
    let t6505 = (f__ols({ let mut out = vec![]; for t6503 in py::iter(&(v_A.clone()))? { let mut v_x: Py; 
v_x = t6503;
out.push(py::math1("log", &v_x.clone())?); } Py::list(out) }, { let mut out = vec![]; for t6504 in py::iter(&(v_S.clone()))? { let mut v_y: Py; 
v_y = t6504;
out.push(py::math1("log", &v_y.clone())?); } Py::list(out) })?).items()?;
    if t6505.len() != 3 { return Err(PyErr::Value("unpack".into())); }
    v_b = t6505[0].clone();
    v_a = t6505[1].clone();
    v_r2 = t6505[2].clone();
    return Ok(f__out(Py::dict(vec![(Py::str("c"), (if (v_no_exp.clone()).truthy() { v_a.clone() } else { py::math1("exp", &v_a.clone())? })), (Py::str("z"), v_b.clone()), (Py::str("r2"), v_r2.clone())]))?);
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn f_ss_fraction(mut v_seq: Py, mut v_order: Py) -> R {
    let mut v_f: Py = Py::None;
    let mut v_s: Py = Py::None;
    v_s = f__seq(v_seq.clone(), Py::str("seq"), Py::str("PROTEIN"))?;
    v_f = { let mut out = vec![]; for t6172 in py::iter(&(py::zip(&[v_order.clone(), Py::tuple(vec![c__HELIX()?, c__TURN()?, c__SHEET()?])])?))? { let mut v_g: Py; let mut v_k: Py; 
let t6173 = (t6172).items()?;
if t6173.len() != 2 { return Err(PyErr::Value("unpack".into())); }
v_k = t6173[0].clone();
v_g = t6173[1].clone();
out.push((v_k.clone(), py::truediv(&py::sum(&{ let mut out = vec![]; for t6174 in py::iter(&(v_g.clone()))? { let mut v_a: Py; 
v_a = t6174;
out.push(py::method(&v_s.clone(), "count", &[v_a.clone()])?); } Py::list(out) }, None)?, &py::len(&v_s.clone())?)?)); } Py::dict(out) };
    return Ok(v_f.clone());
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn f_tajima_d(mut v_seqs: Py, mut v_sqrt_denominator: Py) -> R {
    let mut v_S: Py = Py::None;
    let mut v_a1: Py = Py::None;
    let mut v_a2: Py = Py::None;
    let mut v_b1: Py = Py::None;
    let mut v_b2: Py = Py::None;
    let mut v_c1: Py = Py::None;
    let mut v_c2: Py = Py::None;
    let mut v_e1: Py = Py::None;
    let mut v_e2: Py = Py::None;
    let mut v_k: Py = Py::None;
    let mut v_n: Py = Py::None;
    let mut v_rows: Py = Py::None;
    let mut v_var: Py = Py::None;
    v_rows = f__alignment(v_seqs.clone(), Py::Int(4i64))?;
    v_n = py::len(&v_rows.clone())?;
    v_S = py::sum(&{ let mut out = vec![]; for t6570 in py::iter(&(py::zip(&{ let mut v: Vec<Py> = vec![]; v.extend((v_rows.clone()).items()?); v })?))? { let mut v_col: Py; 
v_col = t6570;
if ({ let t6571 = py::len(&py::set(&v_col.clone())?)?; Py::Bool({ let t6572 = Py::Int(1i64); if py::gt(&t6571, &t6572)? { true } else { false } }) }).truthy() { out.push(Py::Int(1i64)); }} Py::list(out) }, None)?;
    let _ = f__need({ let t6573 = v_S.clone(); Py::Bool({ let t6574 = Py::Int(1i64); if py::ge(&t6573, &t6574)? { true } else { false } }) }, Py::str("Tajima's D needs at least one segregating site"))?;
    v_k = py::float(&f__pairwise_k(v_rows.clone())?)?;
    v_a1 = py::fsum(&{ let mut out = vec![]; for t6581 in py::range_iter(&[(Py::Int(1i64)).i()?, (v_n.clone()).i()?]) { let mut v_i: Py; 
v_i = t6581;
out.push(py::truediv(&Py::Float(1.0), &v_i.clone())?); } Py::list(out) })?;
    v_a2 = py::fsum(&{ let mut out = vec![]; for t6582 in py::range_iter(&[(Py::Int(1i64)).i()?, (v_n.clone()).i()?]) { let mut v_i: Py; 
v_i = t6582;
out.push(py::truediv(&Py::Float(1.0), &py::mul(&v_i.clone(), &v_i.clone())?)?); } Py::list(out) })?;
    v_b1 = py::truediv(&py::add(&v_n.clone(), &Py::Int(1i64))?, &py::mul(&Py::Int(3i64), &py::sub(&v_n.clone(), &Py::Int(1i64))?)?)?;
    v_b2 = py::truediv(&py::mul(&Py::Int(2i64), &py::add(&py::add(&py::mul(&v_n.clone(), &v_n.clone())?, &v_n.clone())?, &Py::Int(3i64))?)?, &py::mul(&py::mul(&Py::Int(9i64), &v_n.clone())?, &py::sub(&v_n.clone(), &Py::Int(1i64))?)?)?;
    v_c1 = py::sub(&v_b1.clone(), &py::truediv(&Py::Int(1i64), &v_a1.clone())?)?;
    v_c2 = py::add(&py::sub(&v_b2.clone(), &py::truediv(&py::add(&v_n.clone(), &Py::Int(2i64))?, &py::mul(&v_a1.clone(), &v_n.clone())?)?)?, &py::truediv(&v_a2.clone(), &py::mul(&v_a1.clone(), &v_a1.clone())?)?)?;
    let t6583 = (Py::tuple(vec![py::truediv(&v_c1.clone(), &v_a1.clone())?, py::truediv(&v_c2.clone(), &py::add(&py::mul(&v_a1.clone(), &v_a1.clone())?, &v_a2.clone())?)?])).items()?;
    if t6583.len() != 2 { return Err(PyErr::Value("unpack".into())); }
    v_e1 = t6583[0].clone();
    v_e2 = t6583[1].clone();
    v_var = py::add(&py::mul(&v_e1.clone(), &v_S.clone())?, &py::mul(&py::mul(&v_e2.clone(), &v_S.clone())?, &py::sub(&v_S.clone(), &Py::Int(1i64))?)?)?;
    return Ok(f__out(Py::dict(vec![(Py::str("D"), py::truediv(&py::sub(&v_k.clone(), &py::truediv(&v_S.clone(), &v_a1.clone())?)?, &(if (v_sqrt_denominator.clone()).truthy() { py::math1("sqrt", &v_var.clone())? } else { v_var.clone() }))?), (Py::str("k"), v_k.clone())]))?);
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn f_tm_nn(mut v_seq: Py, mut v_Na: Py, mut v_dnac1: Py, mut v_dnac2: Py) -> R {
    let mut v_s: Py = Py::None;
    v_s = f__seq(v_seq.clone(), Py::str("seq"), Py::str("DNA"))?;
    let _ = f__need({ let t6662 = py::len(&v_s.clone())?; Py::Bool({ let t6663 = Py::Int(2i64); if py::ge(&t6662, &t6663)? { true } else { false } }) }, Py::str("the sequence needs at least 2 bases"))?;
    let _ = f__pos(Py::str("Na"), Py::tuple(vec![v_Na.clone()]))?;
    let _ = f__pos(Py::str("dnac1"), Py::tuple(vec![v_dnac1.clone()]))?;
    let _ = f__fin(Py::str("dnac2"), Py::tuple(vec![v_dnac2.clone()]))?;
    let _ = f__need({ let t6673 = v_dnac2.clone(); Py::Bool({ let t6674 = Py::Int(0i64); if py::ge(&t6673, &t6674)? { true } else { false } }) }, Py::str("dnac2 must be zero or positive"))?;
    let _ = f__need({ let t6679 = { let t6675 = v_s.clone(); Py::Bool({ let t6676 = f__revcomp(v_s.clone())?; if py::eq(&t6675, &t6676) { true } else { false } }) }; if t6679.truthy() { t6679 } else { { let t6677 = v_dnac1.clone(); Py::Bool({ let t6678 = py::truediv(&v_dnac2.clone(), &Py::Int(2i64))?; if py::gt(&t6677, &t6678)? { true } else { false } }) } } }, Py::str("dnac1 must exceed dnac2/2 (k = c1 - c2/2 must be positive)"))?;
    return Ok(f__out(Py::dict(vec![(Py::str("tm"), f__tm_nn(v_s.clone(), v_Na.clone(), v_dnac1.clone(), v_dnac2.clone(), Py::Bool(true))?)]))?);
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn f_tm_wallace(mut v_seq: Py) -> R {
    let mut v_s: Py = Py::None;
    v_s = f__seq(v_seq.clone(), Py::str("seq"), Py::str("DNA"))?;
    return Ok(Py::dict(vec![(Py::str("tm"), py::add(&py::mul(&Py::Int(2i64), &py::add(&py::method(&v_s.clone(), "count", &[Py::str("A")])?, &py::method(&v_s.clone(), "count", &[Py::str("T")])?)?)?, &py::mul(&Py::Int(4i64), &py::add(&py::method(&v_s.clone(), "count", &[Py::str("G")])?, &py::method(&v_s.clone(), "count", &[Py::str("C")])?)?)?)?)]));
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn f_transcription(mut v_seq: Py, mut v_strand: Py) -> R {
    let mut v_m: Py = Py::None;
    let mut v_s: Py = Py::None;
    v_s = f__seq(v_seq.clone(), Py::str("seq"), Py::str("DNA"))?;
    let _ = f__need({ let t6833 = v_strand.clone(); Py::Bool({ let t6834 = Py::tuple(vec![Py::str("coding"), Py::str("template")]); if py::contains(&t6834, &t6833)? { true } else { false } }) }, Py::str("strand must be coding or template"))?;
    v_m = (if ({ let t6835 = v_strand.clone(); Py::Bool({ let t6836 = Py::str("coding"); if py::eq(&t6835, &t6836) { true } else { false } }) }).truthy() { v_s.clone() } else { f__revcomp(v_s.clone())? });
    return Ok(Py::dict(vec![(Py::str("mrna"), py::method(&v_m.clone(), "replace", &[Py::str("T"), Py::str("U")])?)]));
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn f_translation(mut v_seq: Py, mut v_stops: Py) -> R {
    let mut v_p: Py = Py::None;
    let mut v_s: Py = Py::None;
    v_s = py::method(&f__seq(v_seq.clone(), Py::str("seq"), Py::str("NUC"))?, "replace", &[Py::str("U"), Py::str("T")])?;
    let _ = f__need({ let t6901 = py::modulo(&py::len(&v_s.clone())?, &Py::Int(3i64))?; Py::Bool({ let t6902 = Py::Int(0i64); if py::eq(&t6901, &t6902) { true } else { false } }) }, Py::str("the sequence length must be a multiple of 3 (whole codons read in frame 1)"))?;
    let _ = f__need({ let t6903 = v_stops.clone(); Py::Bool({ let t6904 = Py::tuple(vec![Py::str("show"), Py::str("truncate")]); if py::contains(&t6904, &t6903)? { true } else { false } }) }, Py::str("stops must be show or truncate"))?;
    v_p = f__translate(v_s.clone(), c_CODON()?, { let t6909 = v_stops.clone(); Py::Bool({ let t6910 = Py::str("truncate"); if py::eq(&t6909, &t6910) { true } else { false } }) })?;
    return Ok(Py::dict(vec![(Py::str("protein"), v_p.clone()), (Py::str("length"), py::len(&v_p.clone())?)]));
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn f_tree_count(mut v_n: Py, mut v_rooted: Py) -> R {
    let mut v_lo: Py = Py::None;
    let _ = f__need({ let t6952 = v_rooted.clone(); Py::Bool({ let t6953 = Py::tuple(vec![Py::str("rooted"), Py::str("unrooted")]); if py::contains(&t6953, &t6952)? { true } else { false } }) }, Py::str("rooted must be rooted or unrooted"))?;
    v_lo = (if ({ let t6954 = v_rooted.clone(); Py::Bool({ let t6955 = Py::str("rooted"); if py::eq(&t6954, &t6955) { true } else { false } }) }).truthy() { Py::Int(2i64) } else { Py::Int(3i64) });
    let _ = f__need({ let t6959 = Py::Bool(matches!(v_n.clone(), Py::Int(_) | Py::Big(_) | Py::Bool(_))); if !t6959.truthy() { t6959 } else { { let t6956 = v_lo.clone(); Py::Bool({ let t6957 = v_n.clone(); if py::le(&t6956, &t6957)? { { let t6958 = Py::Int(100i64); if py::le(&t6957, &t6958)? { true } else { false } } } else { false } }) } } }, Py::str(&["n must be a whole number of taxa from ".to_string(), py::format(&v_lo.clone(), "")?, " to 100".to_string()].concat()))?;
    return Ok(Py::dict(vec![(Py::str("count"), (if ({ let t6960 = v_rooted.clone(); Py::Bool({ let t6961 = Py::str("rooted"); if py::eq(&t6960, &t6961) { true } else { false } }) }).truthy() { f__dfact(py::sub(&py::mul(&Py::Int(2i64), &v_n.clone())?, &Py::Int(3i64))?)? } else { f__dfact(py::sub(&py::mul(&Py::Int(2i64), &v_n.clone())?, &Py::Int(5i64))?)? }))]));
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn f_upgma(mut v_D: Py, mut v_linkage: Py) -> R {
    let mut v_M: Py = Py::None;
    let mut v_best: Py = Py::None;
    let mut v_d: Py = Py::None;
    let mut v_dik: Py = Py::None;
    let mut v_djk: Py = Py::None;
    let mut v_heights: Py = Py::None;
    let mut v_i: Py = Py::None;
    let mut v_j: Py = Py::None;
    let mut v_k: Py = Py::None;
    let mut v_n: Py = Py::None;
    let mut v_new: Py = Py::None;
    let mut v_nxt: Py = Py::None;
    let mut v_scale: Py = Py::None;
    let mut v_size: Py = Py::None;
    v_M = f__distance_matrix(v_D.clone(), Py::str("D"), Py::Int(2i64), Py::Int(30i64))?;
    v_n = py::len(&v_M.clone())?;
    v_size = { let mut out = vec![]; for t7022 in py::range_iter(&[(v_n.clone()).i()?]) { let mut v_i: Py; 
v_i = t7022;
out.push((v_i.clone(), Py::Int(1i64))); } Py::dict(out) };
    v_d = { let mut out = vec![]; for t7023 in py::range_iter(&[(v_n.clone()).i()?]) { let mut v_i: Py; 
v_i = t7023;
for t7024 in py::range_iter(&[(py::add(&v_i.clone(), &Py::Int(1i64))?).i()?, (v_n.clone()).i()?]) { let mut v_j: Py; 
v_j = t7024;
out.push((Py::tuple(vec![v_i.clone(), v_j.clone()]), py::getitem(&py::getitem(&v_M.clone(), &v_i.clone())?, &v_j.clone())?)); }} Py::dict(out) };
    v_scale = py::minmax(&[py::method(&v_d.clone(), "values", &[])?], true)?;
    let t7025 = (Py::tuple(vec![Py::list(vec![]), v_n.clone()])).items()?;
    if t7025.len() != 2 { return Err(PyErr::Value("unpack".into())); }
    v_heights = t7025[0].clone();
    v_nxt = t7025[1].clone();
    while ({ let t7049 = py::len(&v_size.clone())?; Py::Bool({ let t7050 = Py::Int(1i64); if py::gt(&t7049, &t7050)? { true } else { false } }) }).truthy() {
        v_best = py::minmax(&[py::method(&v_d.clone(), "values", &[])?], false)?;
        let _ = f__need(py::not(&f__tie(Py::list((py::method(&v_d.clone(), "values", &[])?).items()?), v_best.clone(), v_scale.clone())?), Py::str("two pairs tie at the closest distance: the tree is not unique"))?;
        let t7035 = (py::minmax(&[{ let mut out = vec![]; for t7031 in py::iter(&(py::method(&v_d.clone(), "items", &[])?))? { let mut v_k: Py; let mut v_v: Py; 
let t7032 = (t7031).items()?;
if t7032.len() != 2 { return Err(PyErr::Value("unpack".into())); }
v_k = t7032[0].clone();
v_v = t7032[1].clone();
if ({ let t7033 = v_v.clone(); Py::Bool({ let t7034 = v_best.clone(); if py::eq(&t7033, &t7034) { true } else { false } }) }).truthy() { out.push(v_k.clone()); }} Py::list(out) }], false)?).items()?;
        if t7035.len() != 2 { return Err(PyErr::Value("unpack".into())); }
        v_i = t7035[0].clone();
        v_j = t7035[1].clone();
        let _ = py::method(&v_heights.clone(), "append", &[py::truediv(&v_best.clone(), &Py::Int(2i64))?])?;
        v_new = v_nxt.clone();
        v_nxt = py::iadd(&v_nxt, &Py::Int(1i64))?;
        for t7036 in py::iter(&(Py::list((v_size.clone()).items()?)))? {
            v_k = t7036;
            if ({ let t7037 = v_k.clone(); Py::Bool({ let t7038 = Py::tuple(vec![v_i.clone(), v_j.clone()]); if py::contains(&t7038, &t7037)? { true } else { false } }) }).truthy() {
                continue;
            }
            let t7039 = (Py::tuple(vec![py::getitem(&v_d.clone(), &Py::tuple((py::sorted(&Py::tuple(vec![v_i.clone(), v_k.clone()]), (Py::Bool(false)).truthy())?).items()?))?, py::getitem(&v_d.clone(), &Py::tuple((py::sorted(&Py::tuple(vec![v_j.clone(), v_k.clone()]), (Py::Bool(false)).truthy())?).items()?))?])).items()?;
            if t7039.len() != 2 { return Err(PyErr::Value("unpack".into())); }
            v_dik = t7039[0].clone();
            v_djk = t7039[1].clone();
            if ({ let t7040 = v_linkage.clone(); Py::Bool({ let t7041 = Py::str("average"); if py::eq(&t7040, &t7041) { true } else { false } }) }).truthy() {
                py::setitem(&v_d.clone(), &Py::tuple(vec![v_k.clone(), v_new.clone()]), py::truediv(&py::add(&py::mul(&py::getitem(&v_size.clone(), &v_i.clone())?, &v_dik.clone())?, &py::mul(&py::getitem(&v_size.clone(), &v_j.clone())?, &v_djk.clone())?)?, &py::add(&py::getitem(&v_size.clone(), &v_i.clone())?, &py::getitem(&v_size.clone(), &v_j.clone())?)?)?)?;
            } else {
                py::setitem(&v_d.clone(), &Py::tuple(vec![v_k.clone(), v_new.clone()]), py::minmax(&[v_dik.clone(), v_djk.clone()], false)?)?;
            }
        }
        py::setitem(&v_size.clone(), &v_new.clone(), py::add(&py::method(&v_size.clone(), "pop", &[v_i.clone()])?, &py::method(&v_size.clone(), "pop", &[v_j.clone()])?)?)?;
        v_d = { let mut out = vec![]; for t7042 in py::iter(&(py::method(&v_d.clone(), "items", &[])?))? { let mut v_k: Py; let mut v_v: Py; 
let t7043 = (t7042).items()?;
if t7043.len() != 2 { return Err(PyErr::Value("unpack".into())); }
v_k = t7043[0].clone();
v_v = t7043[1].clone();
if ({ let t7048 = { let t7044 = v_i.clone(); Py::Bool({ let t7045 = v_k.clone(); if !py::contains(&t7045, &t7044)? { true } else { false } }) }; if !t7048.truthy() { t7048 } else { { let t7046 = v_j.clone(); Py::Bool({ let t7047 = v_k.clone(); if !py::contains(&t7047, &t7046)? { true } else { false } }) } } }).truthy() { out.push((v_k.clone(), v_v.clone())); }} Py::dict(out) };
    }
    return Ok(f__out(Py::dict(vec![(Py::str("root_height"), py::getitem(&v_heights.clone(), &py::neg(&Py::Int(1i64))?)?), (Py::str("heights"), py::sorted(&v_heights.clone(), (Py::Bool(false)).truthy())?)]))?);
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn f_watterson(mut v_seqs: Py, mut v_harmonic_to_n: Py) -> R {
    let mut v_L: Py = Py::None;
    let mut v_S: Py = Py::None;
    let mut v_a: Py = Py::None;
    let mut v_n: Py = Py::None;
    let mut v_rows: Py = Py::None;
    v_rows = f__alignment(v_seqs.clone(), Py::Int(2i64))?;
    let t7115 = (Py::tuple(vec![py::len(&v_rows.clone())?, py::len(&py::getitem(&v_rows.clone(), &Py::Int(0i64))?)?])).items()?;
    if t7115.len() != 2 { return Err(PyErr::Value("unpack".into())); }
    v_n = t7115[0].clone();
    v_L = t7115[1].clone();
    v_S = py::sum(&{ let mut out = vec![]; for t7116 in py::iter(&(py::zip(&{ let mut v: Vec<Py> = vec![]; v.extend((v_rows.clone()).items()?); v })?))? { let mut v_col: Py; 
v_col = t7116;
if ({ let t7117 = py::len(&py::set(&v_col.clone())?)?; Py::Bool({ let t7118 = Py::Int(1i64); if py::gt(&t7117, &t7118)? { true } else { false } }) }).truthy() { out.push(Py::Int(1i64)); }} Py::list(out) }, None)?;
    v_a = py::fsum(&{ let mut out = vec![]; for t7119 in py::range_iter(&[(Py::Int(1i64)).i()?, (py::add(&v_n.clone(), &(if (v_harmonic_to_n.clone()).truthy() { Py::Int(1i64) } else { Py::Int(0i64) }))?).i()?]) { let mut v_i: Py; 
v_i = t7119;
out.push(py::truediv(&Py::Float(1.0), &v_i.clone())?); } Py::list(out) })?;
    return Ok(f__out(Py::dict(vec![(Py::str("S"), v_S.clone()), (Py::str("theta"), py::truediv(&v_S.clone(), &v_a.clone())?), (Py::str("theta_site"), py::truediv(&py::truediv(&v_S.clone(), &v_a.clone())?, &v_L.clone())?)]))?);
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn s_1014(mut v_inputs: Py, mut v_outputs: Py) -> R {
    let mut v_k: Py = Py::None;
    let mut v_v: Py = Py::None;
    let mut v_want: Py = Py::None;
    let mut v_worst: Py = Py::None;
    let t1016: Result<Option<Py>, PyErr> = (|| -> Result<Option<Py>, PyErr> {
            v_want = { let t1015 = v_inputs.clone(); f_dominance_indices(py::getitem(&t1015, &Py::str("counts"))?, py::method(&t1015, "get", &[Py::str("wrong"), Py::Bool(false)])?)? };
            #[allow(unreachable_code)]
            Ok(None)
        })();
    match t1016 {
        Ok(Some(r)) => { return Ok(r); }
        Ok(None) => {}
        Err(e) if (!matches!(e, PyErr::Unsupported(_))) => {
            return Ok(Py::Float(f64::INFINITY));
        }
        Err(e) => return Err(e),
    }
    v_worst = Py::Float(0.0);
    for t1017 in py::iter(&(py::method(&v_want.clone(), "items", &[])?))? {
        let t1054 = (t1017).items()?;
        if t1054.len() != 2 { return Err(PyErr::Value("unpack".into())); }
        v_k = t1054[0].clone();
        v_v = t1054[1].clone();
        if ({ let t1018 = v_k.clone(); Py::Bool({ let t1019 = v_outputs.clone(); if !py::contains(&t1019, &t1018)? { true } else { false } }) }).truthy() {
            return Ok(Py::Float(f64::INFINITY));
        }
        v_worst = py::minmax(&[v_worst.clone(), f__rel(py::getitem(&v_outputs.clone(), &v_k.clone())?, v_v.clone())?], true)?;
    }
    return Ok(v_worst.clone());
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn s_1081(mut v_inputs: Py, mut v_outputs: Py) -> R {
    let mut v_k: Py = Py::None;
    let mut v_v: Py = Py::None;
    let mut v_want: Py = Py::None;
    let mut v_worst: Py = Py::None;
    let t1083: Result<Option<Py>, PyErr> = (|| -> Result<Option<Py>, PyErr> {
            v_want = { let t1082 = v_inputs.clone(); f_drift_heterozygosity(py::getitem(&t1082, &Py::str("N"))?, py::getitem(&t1082, &Py::str("i"))?, py::getitem(&t1082, &Py::str("t"))?)? };
            #[allow(unreachable_code)]
            Ok(None)
        })();
    match t1083 {
        Ok(Some(r)) => { return Ok(r); }
        Ok(None) => {}
        Err(e) if (!matches!(e, PyErr::Unsupported(_))) => {
            return Ok(Py::Float(f64::INFINITY));
        }
        Err(e) => return Err(e),
    }
    v_worst = Py::Float(0.0);
    for t1084 in py::iter(&(py::method(&v_want.clone(), "items", &[])?))? {
        let t1121 = (t1084).items()?;
        if t1121.len() != 2 { return Err(PyErr::Value("unpack".into())); }
        v_k = t1121[0].clone();
        v_v = t1121[1].clone();
        if ({ let t1085 = v_k.clone(); Py::Bool({ let t1086 = v_outputs.clone(); if !py::contains(&t1086, &t1085)? { true } else { false } }) }).truthy() {
            return Ok(Py::Float(f64::INFINITY));
        }
        v_worst = py::minmax(&[v_worst.clone(), f__rel(py::getitem(&v_outputs.clone(), &v_k.clone())?, v_v.clone())?], true)?;
    }
    return Ok(v_worst.clone());
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn s_1170(mut v_inputs: Py, mut v_outputs: Py) -> R {
    let mut v_k: Py = Py::None;
    let mut v_v: Py = Py::None;
    let mut v_want: Py = Py::None;
    let mut v_worst: Py = Py::None;
    let t1172: Result<Option<Py>, PyErr> = (|| -> Result<Option<Py>, PyErr> {
            v_want = { let t1171 = v_inputs.clone(); f_eadie_hofstee(py::getitem(&t1171, &Py::str("S"))?, py::getitem(&t1171, &Py::str("v"))?)? };
            #[allow(unreachable_code)]
            Ok(None)
        })();
    match t1172 {
        Ok(Some(r)) => { return Ok(r); }
        Ok(None) => {}
        Err(e) if (!matches!(e, PyErr::Unsupported(_))) => {
            return Ok(Py::Float(f64::INFINITY));
        }
        Err(e) => return Err(e),
    }
    v_worst = Py::Float(0.0);
    for t1173 in py::iter(&(py::method(&v_want.clone(), "items", &[])?))? {
        let t1210 = (t1173).items()?;
        if t1210.len() != 2 { return Err(PyErr::Value("unpack".into())); }
        v_k = t1210[0].clone();
        v_v = t1210[1].clone();
        if ({ let t1174 = v_k.clone(); Py::Bool({ let t1175 = v_outputs.clone(); if !py::contains(&t1175, &t1174)? { true } else { false } }) }).truthy() {
            return Ok(Py::Float(f64::INFINITY));
        }
        v_worst = py::minmax(&[v_worst.clone(), f__rel(py::getitem(&v_outputs.clone(), &v_k.clone())?, v_v.clone())?], true)?;
    }
    return Ok(v_worst.clone());
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn s_118(mut v_inputs: Py, mut v_outputs: Py) -> R {
    let mut v_k: Py = Py::None;
    let mut v_v: Py = Py::None;
    let mut v_want: Py = Py::None;
    let mut v_worst: Py = Py::None;
    let t120: Result<Option<Py>, PyErr> = (|| -> Result<Option<Py>, PyErr> {
            v_want = { let t119 = v_inputs.clone(); f_aa_composition(py::getitem(&t119, &Py::str("seq"))?)? };
            #[allow(unreachable_code)]
            Ok(None)
        })();
    match t120 {
        Ok(Some(r)) => { return Ok(r); }
        Ok(None) => {}
        Err(e) if (!matches!(e, PyErr::Unsupported(_))) => {
            return Ok(Py::Float(f64::INFINITY));
        }
        Err(e) => return Err(e),
    }
    v_worst = Py::Float(0.0);
    for t121 in py::iter(&(py::method(&v_want.clone(), "items", &[])?))? {
        let t158 = (t121).items()?;
        if t158.len() != 2 { return Err(PyErr::Value("unpack".into())); }
        v_k = t158[0].clone();
        v_v = t158[1].clone();
        if ({ let t122 = v_k.clone(); Py::Bool({ let t123 = v_outputs.clone(); if !py::contains(&t123, &t122)? { true } else { false } }) }).truthy() {
            return Ok(Py::Float(f64::INFINITY));
        }
        v_worst = py::minmax(&[v_worst.clone(), f__rel(py::getitem(&v_outputs.clone(), &v_k.clone())?, v_v.clone())?], true)?;
    }
    return Ok(v_worst.clone());
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn s_1241(mut v_inputs: Py, mut v_outputs: Py) -> R {
    let mut v_k: Py = Py::None;
    let mut v_v: Py = Py::None;
    let mut v_want: Py = Py::None;
    let mut v_worst: Py = Py::None;
    let t1243: Result<Option<Py>, PyErr> = (|| -> Result<Option<Py>, PyErr> {
            v_want = { let t1242 = v_inputs.clone(); f_edit_distance(py::getitem(&t1242, &Py::str("seq1"))?, py::getitem(&t1242, &Py::str("seq2"))?)? };
            #[allow(unreachable_code)]
            Ok(None)
        })();
    match t1243 {
        Ok(Some(r)) => { return Ok(r); }
        Ok(None) => {}
        Err(e) if (!matches!(e, PyErr::Unsupported(_))) => {
            return Ok(Py::Float(f64::INFINITY));
        }
        Err(e) => return Err(e),
    }
    v_worst = Py::Float(0.0);
    for t1244 in py::iter(&(py::method(&v_want.clone(), "items", &[])?))? {
        let t1281 = (t1244).items()?;
        if t1281.len() != 2 { return Err(PyErr::Value("unpack".into())); }
        v_k = t1281[0].clone();
        v_v = t1281[1].clone();
        if ({ let t1245 = v_k.clone(); Py::Bool({ let t1246 = v_outputs.clone(); if !py::contains(&t1246, &t1245)? { true } else { false } }) }).truthy() {
            return Ok(Py::Float(f64::INFINITY));
        }
        v_worst = py::minmax(&[v_worst.clone(), f__rel(py::getitem(&v_outputs.clone(), &v_k.clone())?, v_v.clone())?], true)?;
    }
    return Ok(v_worst.clone());
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn s_1349(mut v_inputs: Py, mut v_outputs: Py) -> R {
    let mut v_k: Py = Py::None;
    let mut v_v: Py = Py::None;
    let mut v_want: Py = Py::None;
    let mut v_worst: Py = Py::None;
    let t1351: Result<Option<Py>, PyErr> = (|| -> Result<Option<Py>, PyErr> {
            v_want = { let t1350 = v_inputs.clone(); f_faith_pd(py::getitem(&t1350, &Py::str("tree"))?, py::getitem(&t1350, &Py::str("taxa"))?, py::method(&t1350, "get", &[Py::str("unrooted"), Py::Bool(false)])?)? };
            #[allow(unreachable_code)]
            Ok(None)
        })();
    match t1351 {
        Ok(Some(r)) => { return Ok(r); }
        Ok(None) => {}
        Err(e) if (!matches!(e, PyErr::Unsupported(_))) => {
            return Ok(Py::Float(f64::INFINITY));
        }
        Err(e) => return Err(e),
    }
    v_worst = Py::Float(0.0);
    for t1352 in py::iter(&(py::method(&v_want.clone(), "items", &[])?))? {
        let t1389 = (t1352).items()?;
        if t1389.len() != 2 { return Err(PyErr::Value("unpack".into())); }
        v_k = t1389[0].clone();
        v_v = t1389[1].clone();
        if ({ let t1353 = v_k.clone(); Py::Bool({ let t1354 = v_outputs.clone(); if !py::contains(&t1354, &t1353)? { true } else { false } }) }).truthy() {
            return Ok(Py::Float(f64::INFINITY));
        }
        v_worst = py::minmax(&[v_worst.clone(), f__rel(py::getitem(&v_outputs.clone(), &v_k.clone())?, v_v.clone())?], true)?;
    }
    return Ok(v_worst.clone());
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn s_1462(mut v_inputs: Py, mut v_outputs: Py) -> R {
    let mut v_k: Py = Py::None;
    let mut v_v: Py = Py::None;
    let mut v_want: Py = Py::None;
    let mut v_worst: Py = Py::None;
    let t1464: Result<Option<Py>, PyErr> = (|| -> Result<Option<Py>, PyErr> {
            v_want = { let t1463 = v_inputs.clone(); f_fitch_parsimony(py::getitem(&t1463, &Py::str("tree"))?, py::getitem(&t1463, &Py::str("states"))?, py::method(&t1463, "get", &[Py::str("lower_bound"), Py::Bool(false)])?)? };
            #[allow(unreachable_code)]
            Ok(None)
        })();
    match t1464 {
        Ok(Some(r)) => { return Ok(r); }
        Ok(None) => {}
        Err(e) if (!matches!(e, PyErr::Unsupported(_))) => {
            return Ok(Py::Float(f64::INFINITY));
        }
        Err(e) => return Err(e),
    }
    v_worst = Py::Float(0.0);
    for t1465 in py::iter(&(py::method(&v_want.clone(), "items", &[])?))? {
        let t1502 = (t1465).items()?;
        if t1502.len() != 2 { return Err(PyErr::Value("unpack".into())); }
        v_k = t1502[0].clone();
        v_v = t1502[1].clone();
        if ({ let t1466 = v_k.clone(); Py::Bool({ let t1467 = v_outputs.clone(); if !py::contains(&t1467, &t1466)? { true } else { false } }) }).truthy() {
            return Ok(Py::Float(f64::INFINITY));
        }
        v_worst = py::minmax(&[v_worst.clone(), f__rel(py::getitem(&v_outputs.clone(), &v_k.clone())?, v_v.clone())?], true)?;
    }
    return Ok(v_worst.clone());
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn s_1602(mut v_inputs: Py, mut v_outputs: Py) -> R {
    let mut v_k: Py = Py::None;
    let mut v_v: Py = Py::None;
    let mut v_want: Py = Py::None;
    let mut v_worst: Py = Py::None;
    let t1604: Result<Option<Py>, PyErr> = (|| -> Result<Option<Py>, PyErr> {
            v_want = { let t1603 = v_inputs.clone(); f_fba(py::getitem(&t1603, &Py::str("S"))?, py::getitem(&t1603, &Py::str("lb"))?, py::getitem(&t1603, &Py::str("ub"))?, py::getitem(&t1603, &Py::str("c"))?, py::method(&t1603, "get", &[Py::str("box_only"), Py::Bool(false)])?)? };
            #[allow(unreachable_code)]
            Ok(None)
        })();
    match t1604 {
        Ok(Some(r)) => { return Ok(r); }
        Ok(None) => {}
        Err(e) if (!matches!(e, PyErr::Unsupported(_))) => {
            return Ok(Py::Float(f64::INFINITY));
        }
        Err(e) => return Err(e),
    }
    v_worst = Py::Float(0.0);
    for t1605 in py::iter(&(py::method(&v_want.clone(), "items", &[])?))? {
        let t1642 = (t1605).items()?;
        if t1642.len() != 2 { return Err(PyErr::Value("unpack".into())); }
        v_k = t1642[0].clone();
        v_v = t1642[1].clone();
        if ({ let t1606 = v_k.clone(); Py::Bool({ let t1607 = v_outputs.clone(); if !py::contains(&t1607, &t1606)? { true } else { false } }) }).truthy() {
            return Ok(Py::Float(f64::INFINITY));
        }
        v_worst = py::minmax(&[v_worst.clone(), f__rel(py::getitem(&v_outputs.clone(), &v_k.clone())?, v_v.clone())?], true)?;
    }
    return Ok(v_worst.clone());
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn s_1677(mut v_inputs: Py, mut v_outputs: Py) -> R {
    let mut v_k: Py = Py::None;
    let mut v_v: Py = Py::None;
    let mut v_want: Py = Py::None;
    let mut v_worst: Py = Py::None;
    let t1679: Result<Option<Py>, PyErr> = (|| -> Result<Option<Py>, PyErr> {
            v_want = { let t1678 = v_inputs.clone(); f_fst(py::getitem(&t1678, &Py::str("p"))?)? };
            #[allow(unreachable_code)]
            Ok(None)
        })();
    match t1679 {
        Ok(Some(r)) => { return Ok(r); }
        Ok(None) => {}
        Err(e) if (!matches!(e, PyErr::Unsupported(_))) => {
            return Ok(Py::Float(f64::INFINITY));
        }
        Err(e) => return Err(e),
    }
    v_worst = Py::Float(0.0);
    for t1680 in py::iter(&(py::method(&v_want.clone(), "items", &[])?))? {
        let t1717 = (t1680).items()?;
        if t1717.len() != 2 { return Err(PyErr::Value("unpack".into())); }
        v_k = t1717[0].clone();
        v_v = t1717[1].clone();
        if ({ let t1681 = v_k.clone(); Py::Bool({ let t1682 = v_outputs.clone(); if !py::contains(&t1682, &t1681)? { true } else { false } }) }).truthy() {
            return Ok(Py::Float(f64::INFINITY));
        }
        v_worst = py::minmax(&[v_worst.clone(), f__rel(py::getitem(&v_outputs.clone(), &v_k.clone())?, v_v.clone())?], true)?;
    }
    return Ok(v_worst.clone());
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn s_1741(mut v_inputs: Py, mut v_outputs: Py) -> R {
    let mut v_k: Py = Py::None;
    let mut v_v: Py = Py::None;
    let mut v_want: Py = Py::None;
    let mut v_worst: Py = Py::None;
    let t1743: Result<Option<Py>, PyErr> = (|| -> Result<Option<Py>, PyErr> {
            v_want = { let t1742 = v_inputs.clone(); f_gc_content(py::getitem(&t1742, &Py::str("seq"))?)? };
            #[allow(unreachable_code)]
            Ok(None)
        })();
    match t1743 {
        Ok(Some(r)) => { return Ok(r); }
        Ok(None) => {}
        Err(e) if (!matches!(e, PyErr::Unsupported(_))) => {
            return Ok(Py::Float(f64::INFINITY));
        }
        Err(e) => return Err(e),
    }
    v_worst = Py::Float(0.0);
    for t1744 in py::iter(&(py::method(&v_want.clone(), "items", &[])?))? {
        let t1781 = (t1744).items()?;
        if t1781.len() != 2 { return Err(PyErr::Value("unpack".into())); }
        v_k = t1781[0].clone();
        v_v = t1781[1].clone();
        if ({ let t1745 = v_k.clone(); Py::Bool({ let t1746 = v_outputs.clone(); if !py::contains(&t1746, &t1745)? { true } else { false } }) }).truthy() {
            return Ok(Py::Float(f64::INFINITY));
        }
        v_worst = py::minmax(&[v_worst.clone(), f__rel(py::getitem(&v_outputs.clone(), &v_k.clone())?, v_v.clone())?], true)?;
    }
    return Ok(v_worst.clone());
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn s_182(mut v_inputs: Py, mut v_outputs: Py) -> R {
    let mut v_k: Py = Py::None;
    let mut v_v: Py = Py::None;
    let mut v_want: Py = Py::None;
    let mut v_worst: Py = Py::None;
    let t184: Result<Option<Py>, PyErr> = (|| -> Result<Option<Py>, PyErr> {
            v_want = { let t183 = v_inputs.clone(); f_aromaticity(py::getitem(&t183, &Py::str("seq"))?)? };
            #[allow(unreachable_code)]
            Ok(None)
        })();
    match t184 {
        Ok(Some(r)) => { return Ok(r); }
        Ok(None) => {}
        Err(e) if (!matches!(e, PyErr::Unsupported(_))) => {
            return Ok(Py::Float(f64::INFINITY));
        }
        Err(e) => return Err(e),
    }
    v_worst = Py::Float(0.0);
    for t185 in py::iter(&(py::method(&v_want.clone(), "items", &[])?))? {
        let t222 = (t185).items()?;
        if t222.len() != 2 { return Err(PyErr::Value("unpack".into())); }
        v_k = t222[0].clone();
        v_v = t222[1].clone();
        if ({ let t186 = v_k.clone(); Py::Bool({ let t187 = v_outputs.clone(); if !py::contains(&t187, &t186)? { true } else { false } }) }).truthy() {
            return Ok(Py::Float(f64::INFINITY));
        }
        v_worst = py::minmax(&[v_worst.clone(), f__rel(py::getitem(&v_outputs.clone(), &v_k.clone())?, v_v.clone())?], true)?;
    }
    return Ok(v_worst.clone());
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn s_1863(mut v_inputs: Py, mut v_outputs: Py) -> R {
    let mut v_k: Py = Py::None;
    let mut v_v: Py = Py::None;
    let mut v_want: Py = Py::None;
    let mut v_worst: Py = Py::None;
    let t1865: Result<Option<Py>, PyErr> = (|| -> Result<Option<Py>, PyErr> {
            v_want = { let t1864 = v_inputs.clone(); f_global_affine(py::getitem(&t1864, &Py::str("seq1"))?, py::getitem(&t1864, &Py::str("seq2"))?, py::getitem(&t1864, &Py::str("match"))?, py::getitem(&t1864, &Py::str("mismatch"))?, py::getitem(&t1864, &Py::str("open"))?, py::getitem(&t1864, &Py::str("extend"))?)? };
            #[allow(unreachable_code)]
            Ok(None)
        })();
    match t1865 {
        Ok(Some(r)) => { return Ok(r); }
        Ok(None) => {}
        Err(e) if (!matches!(e, PyErr::Unsupported(_))) => {
            return Ok(Py::Float(f64::INFINITY));
        }
        Err(e) => return Err(e),
    }
    v_worst = Py::Float(0.0);
    for t1866 in py::iter(&(py::method(&v_want.clone(), "items", &[])?))? {
        let t1903 = (t1866).items()?;
        if t1903.len() != 2 { return Err(PyErr::Value("unpack".into())); }
        v_k = t1903[0].clone();
        v_v = t1903[1].clone();
        if ({ let t1867 = v_k.clone(); Py::Bool({ let t1868 = v_outputs.clone(); if !py::contains(&t1868, &t1867)? { true } else { false } }) }).truthy() {
            return Ok(Py::Float(f64::INFINITY));
        }
        v_worst = py::minmax(&[v_worst.clone(), f__rel(py::getitem(&v_outputs.clone(), &v_k.clone())?, v_v.clone())?], true)?;
    }
    return Ok(v_worst.clone());
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn s_1982(mut v_inputs: Py, mut v_outputs: Py) -> R {
    let mut v_k: Py = Py::None;
    let mut v_v: Py = Py::None;
    let mut v_want: Py = Py::None;
    let mut v_worst: Py = Py::None;
    let t1984: Result<Option<Py>, PyErr> = (|| -> Result<Option<Py>, PyErr> {
            v_want = { let t1983 = v_inputs.clone(); f_global_linear(py::getitem(&t1983, &Py::str("seq1"))?, py::getitem(&t1983, &Py::str("seq2"))?, py::getitem(&t1983, &Py::str("match"))?, py::getitem(&t1983, &Py::str("mismatch"))?, py::getitem(&t1983, &Py::str("gap"))?)? };
            #[allow(unreachable_code)]
            Ok(None)
        })();
    match t1984 {
        Ok(Some(r)) => { return Ok(r); }
        Ok(None) => {}
        Err(e) if (!matches!(e, PyErr::Unsupported(_))) => {
            return Ok(Py::Float(f64::INFINITY));
        }
        Err(e) => return Err(e),
    }
    v_worst = Py::Float(0.0);
    for t1985 in py::iter(&(py::method(&v_want.clone(), "items", &[])?))? {
        let t2022 = (t1985).items()?;
        if t2022.len() != 2 { return Err(PyErr::Value("unpack".into())); }
        v_k = t2022[0].clone();
        v_v = t2022[1].clone();
        if ({ let t1986 = v_k.clone(); Py::Bool({ let t1987 = v_outputs.clone(); if !py::contains(&t1987, &t1986)? { true } else { false } }) }).truthy() {
            return Ok(Py::Float(f64::INFINITY));
        }
        v_worst = py::minmax(&[v_worst.clone(), f__rel(py::getitem(&v_outputs.clone(), &v_k.clone())?, v_v.clone())?], true)?;
    }
    return Ok(v_worst.clone());
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn s_2044(mut v_inputs: Py, mut v_outputs: Py) -> R {
    let mut v_k: Py = Py::None;
    let mut v_v: Py = Py::None;
    let mut v_want: Py = Py::None;
    let mut v_worst: Py = Py::None;
    let t2046: Result<Option<Py>, PyErr> = (|| -> Result<Option<Py>, PyErr> {
            v_want = { let t2045 = v_inputs.clone(); f_goods_coverage(py::getitem(&t2045, &Py::str("counts"))?)? };
            #[allow(unreachable_code)]
            Ok(None)
        })();
    match t2046 {
        Ok(Some(r)) => { return Ok(r); }
        Ok(None) => {}
        Err(e) if (!matches!(e, PyErr::Unsupported(_))) => {
            return Ok(Py::Float(f64::INFINITY));
        }
        Err(e) => return Err(e),
    }
    v_worst = Py::Float(0.0);
    for t2047 in py::iter(&(py::method(&v_want.clone(), "items", &[])?))? {
        let t2084 = (t2047).items()?;
        if t2084.len() != 2 { return Err(PyErr::Value("unpack".into())); }
        v_k = t2084[0].clone();
        v_v = t2084[1].clone();
        if ({ let t2048 = v_k.clone(); Py::Bool({ let t2049 = v_outputs.clone(); if !py::contains(&t2049, &t2048)? { true } else { false } }) }).truthy() {
            return Ok(Py::Float(f64::INFINITY));
        }
        v_worst = py::minmax(&[v_worst.clone(), f__rel(py::getitem(&v_outputs.clone(), &v_k.clone())?, v_v.clone())?], true)?;
    }
    return Ok(v_worst.clone());
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn s_2123(mut v_inputs: Py, mut v_outputs: Py) -> R {
    let mut v_k: Py = Py::None;
    let mut v_v: Py = Py::None;
    let mut v_want: Py = Py::None;
    let mut v_worst: Py = Py::None;
    let t2125: Result<Option<Py>, PyErr> = (|| -> Result<Option<Py>, PyErr> {
            v_want = { let t2124 = v_inputs.clone(); f_gravy(py::getitem(&t2124, &Py::str("seq"))?)? };
            #[allow(unreachable_code)]
            Ok(None)
        })();
    match t2125 {
        Ok(Some(r)) => { return Ok(r); }
        Ok(None) => {}
        Err(e) if (!matches!(e, PyErr::Unsupported(_))) => {
            return Ok(Py::Float(f64::INFINITY));
        }
        Err(e) => return Err(e),
    }
    v_worst = Py::Float(0.0);
    for t2126 in py::iter(&(py::method(&v_want.clone(), "items", &[])?))? {
        let t2163 = (t2126).items()?;
        if t2163.len() != 2 { return Err(PyErr::Value("unpack".into())); }
        v_k = t2163[0].clone();
        v_v = t2163[1].clone();
        if ({ let t2127 = v_k.clone(); Py::Bool({ let t2128 = v_outputs.clone(); if !py::contains(&t2128, &t2127)? { true } else { false } }) }).truthy() {
            return Ok(Py::Float(f64::INFINITY));
        }
        v_worst = py::minmax(&[v_worst.clone(), f__rel(py::getitem(&v_outputs.clone(), &v_k.clone())?, v_v.clone())?], true)?;
    }
    return Ok(v_worst.clone());
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn s_2194(mut v_inputs: Py, mut v_outputs: Py) -> R {
    let mut v_k: Py = Py::None;
    let mut v_v: Py = Py::None;
    let mut v_want: Py = Py::None;
    let mut v_worst: Py = Py::None;
    let t2196: Result<Option<Py>, PyErr> = (|| -> Result<Option<Py>, PyErr> {
            v_want = { let t2195 = v_inputs.clone(); f_hamming(py::getitem(&t2195, &Py::str("seq1"))?, py::getitem(&t2195, &Py::str("seq2"))?)? };
            #[allow(unreachable_code)]
            Ok(None)
        })();
    match t2196 {
        Ok(Some(r)) => { return Ok(r); }
        Ok(None) => {}
        Err(e) if (!matches!(e, PyErr::Unsupported(_))) => {
            return Ok(Py::Float(f64::INFINITY));
        }
        Err(e) => return Err(e),
    }
    v_worst = Py::Float(0.0);
    for t2197 in py::iter(&(py::method(&v_want.clone(), "items", &[])?))? {
        let t2234 = (t2197).items()?;
        if t2234.len() != 2 { return Err(PyErr::Value("unpack".into())); }
        v_k = t2234[0].clone();
        v_v = t2234[1].clone();
        if ({ let t2198 = v_k.clone(); Py::Bool({ let t2199 = v_outputs.clone(); if !py::contains(&t2199, &t2198)? { true } else { false } }) }).truthy() {
            return Ok(Py::Float(f64::INFINITY));
        }
        v_worst = py::minmax(&[v_worst.clone(), f__rel(py::getitem(&v_outputs.clone(), &v_k.clone())?, v_v.clone())?], true)?;
    }
    return Ok(v_worst.clone());
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn s_2283(mut v_inputs: Py, mut v_outputs: Py) -> R {
    let mut v_k: Py = Py::None;
    let mut v_v: Py = Py::None;
    let mut v_want: Py = Py::None;
    let mut v_worst: Py = Py::None;
    let t2285: Result<Option<Py>, PyErr> = (|| -> Result<Option<Py>, PyErr> {
            v_want = { let t2284 = v_inputs.clone(); f_hanes_woolf(py::getitem(&t2284, &Py::str("S"))?, py::getitem(&t2284, &Py::str("v"))?)? };
            #[allow(unreachable_code)]
            Ok(None)
        })();
    match t2285 {
        Ok(Some(r)) => { return Ok(r); }
        Ok(None) => {}
        Err(e) if (!matches!(e, PyErr::Unsupported(_))) => {
            return Ok(Py::Float(f64::INFINITY));
        }
        Err(e) => return Err(e),
    }
    v_worst = Py::Float(0.0);
    for t2286 in py::iter(&(py::method(&v_want.clone(), "items", &[])?))? {
        let t2323 = (t2286).items()?;
        if t2323.len() != 2 { return Err(PyErr::Value("unpack".into())); }
        v_k = t2323[0].clone();
        v_v = t2323[1].clone();
        if ({ let t2287 = v_k.clone(); Py::Bool({ let t2288 = v_outputs.clone(); if !py::contains(&t2288, &t2287)? { true } else { false } }) }).truthy() {
            return Ok(Py::Float(f64::INFINITY));
        }
        v_worst = py::minmax(&[v_worst.clone(), f__rel(py::getitem(&v_outputs.clone(), &v_k.clone())?, v_v.clone())?], true)?;
    }
    return Ok(v_worst.clone());
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn s_2357(mut v_inputs: Py, mut v_outputs: Py) -> R {
    let mut v_k: Py = Py::None;
    let mut v_v: Py = Py::None;
    let mut v_want: Py = Py::None;
    let mut v_worst: Py = Py::None;
    let t2359: Result<Option<Py>, PyErr> = (|| -> Result<Option<Py>, PyErr> {
            v_want = { let t2358 = v_inputs.clone(); f_haploid_selection(py::getitem(&t2358, &Py::str("p0"))?, py::getitem(&t2358, &Py::str("s"))?, py::getitem(&t2358, &Py::str("t"))?)? };
            #[allow(unreachable_code)]
            Ok(None)
        })();
    match t2359 {
        Ok(Some(r)) => { return Ok(r); }
        Ok(None) => {}
        Err(e) if (!matches!(e, PyErr::Unsupported(_))) => {
            return Ok(Py::Float(f64::INFINITY));
        }
        Err(e) => return Err(e),
    }
    v_worst = Py::Float(0.0);
    for t2360 in py::iter(&(py::method(&v_want.clone(), "items", &[])?))? {
        let t2397 = (t2360).items()?;
        if t2397.len() != 2 { return Err(PyErr::Value("unpack".into())); }
        v_k = t2397[0].clone();
        v_v = t2397[1].clone();
        if ({ let t2361 = v_k.clone(); Py::Bool({ let t2362 = v_outputs.clone(); if !py::contains(&t2362, &t2361)? { true } else { false } }) }).truthy() {
            return Ok(Py::Float(f64::INFINITY));
        }
        v_worst = py::minmax(&[v_worst.clone(), f__rel(py::getitem(&v_outputs.clone(), &v_k.clone())?, v_v.clone())?], true)?;
    }
    return Ok(v_worst.clone());
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn s_2414(mut v_inputs: Py, mut v_outputs: Py) -> R {
    let mut v_k: Py = Py::None;
    let mut v_v: Py = Py::None;
    let mut v_want: Py = Py::None;
    let mut v_worst: Py = Py::None;
    let t2416: Result<Option<Py>, PyErr> = (|| -> Result<Option<Py>, PyErr> {
            v_want = { let t2415 = v_inputs.clone(); f_hwe_exact(py::getitem(&t2415, &Py::str("n_AA"))?, py::getitem(&t2415, &Py::str("n_Aa"))?, py::getitem(&t2415, &Py::str("n_aa"))?)? };
            #[allow(unreachable_code)]
            Ok(None)
        })();
    match t2416 {
        Ok(Some(r)) => { return Ok(r); }
        Ok(None) => {}
        Err(e) if (!matches!(e, PyErr::Unsupported(_))) => {
            return Ok(Py::Float(f64::INFINITY));
        }
        Err(e) => return Err(e),
    }
    v_worst = Py::Float(0.0);
    for t2417 in py::iter(&(py::method(&v_want.clone(), "items", &[])?))? {
        let t2454 = (t2417).items()?;
        if t2454.len() != 2 { return Err(PyErr::Value("unpack".into())); }
        v_k = t2454[0].clone();
        v_v = t2454[1].clone();
        if ({ let t2418 = v_k.clone(); Py::Bool({ let t2419 = v_outputs.clone(); if !py::contains(&t2419, &t2418)? { true } else { false } }) }).truthy() {
            return Ok(Py::Float(f64::INFINITY));
        }
        v_worst = py::minmax(&[v_worst.clone(), f__rel(py::getitem(&v_outputs.clone(), &v_k.clone())?, v_v.clone())?], true)?;
    }
    return Ok(v_worst.clone());
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn s_2482(mut v_inputs: Py, mut v_outputs: Py) -> R {
    let mut v_k: Py = Py::None;
    let mut v_v: Py = Py::None;
    let mut v_want: Py = Py::None;
    let mut v_worst: Py = Py::None;
    let t2484: Result<Option<Py>, PyErr> = (|| -> Result<Option<Py>, PyErr> {
            v_want = { let t2483 = v_inputs.clone(); f_hardy_weinberg_test(py::getitem(&t2483, &Py::str("n_AA"))?, py::getitem(&t2483, &Py::str("n_Aa"))?, py::getitem(&t2483, &Py::str("n_aa"))?)? };
            #[allow(unreachable_code)]
            Ok(None)
        })();
    match t2484 {
        Ok(Some(r)) => { return Ok(r); }
        Ok(None) => {}
        Err(e) if (!matches!(e, PyErr::Unsupported(_))) => {
            return Ok(Py::Float(f64::INFINITY));
        }
        Err(e) => return Err(e),
    }
    v_worst = Py::Float(0.0);
    for t2485 in py::iter(&(py::method(&v_want.clone(), "items", &[])?))? {
        let t2522 = (t2485).items()?;
        if t2522.len() != 2 { return Err(PyErr::Value("unpack".into())); }
        v_k = t2522[0].clone();
        v_v = t2522[1].clone();
        if ({ let t2486 = v_k.clone(); Py::Bool({ let t2487 = v_outputs.clone(); if !py::contains(&t2487, &t2486)? { true } else { false } }) }).truthy() {
            return Ok(Py::Float(f64::INFINITY));
        }
        v_worst = py::minmax(&[v_worst.clone(), f__rel(py::getitem(&v_outputs.clone(), &v_k.clone())?, v_v.clone())?], true)?;
    }
    return Ok(v_worst.clone());
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn s_2548(mut v_inputs: Py, mut v_outputs: Py) -> R {
    let mut v_k: Py = Py::None;
    let mut v_v: Py = Py::None;
    let mut v_want: Py = Py::None;
    let mut v_worst: Py = Py::None;
    let t2550: Result<Option<Py>, PyErr> = (|| -> Result<Option<Py>, PyErr> {
            v_want = { let t2549 = v_inputs.clone(); f_het_advantage(py::getitem(&t2549, &Py::str("w_AA"))?, py::getitem(&t2549, &Py::str("w_Aa"))?, py::getitem(&t2549, &Py::str("w_aa"))?)? };
            #[allow(unreachable_code)]
            Ok(None)
        })();
    match t2550 {
        Ok(Some(r)) => { return Ok(r); }
        Ok(None) => {}
        Err(e) if (!matches!(e, PyErr::Unsupported(_))) => {
            return Ok(Py::Float(f64::INFINITY));
        }
        Err(e) => return Err(e),
    }
    v_worst = Py::Float(0.0);
    for t2551 in py::iter(&(py::method(&v_want.clone(), "items", &[])?))? {
        let t2588 = (t2551).items()?;
        if t2588.len() != 2 { return Err(PyErr::Value("unpack".into())); }
        v_k = t2588[0].clone();
        v_v = t2588[1].clone();
        if ({ let t2552 = v_k.clone(); Py::Bool({ let t2553 = v_outputs.clone(); if !py::contains(&t2553, &t2552)? { true } else { false } }) }).truthy() {
            return Ok(Py::Float(f64::INFINITY));
        }
        v_worst = py::minmax(&[v_worst.clone(), f__rel(py::getitem(&v_outputs.clone(), &v_k.clone())?, v_v.clone())?], true)?;
    }
    return Ok(v_worst.clone());
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn s_256(mut v_inputs: Py, mut v_outputs: Py) -> R {
    let mut v_k: Py = Py::None;
    let mut v_v: Py = Py::None;
    let mut v_want: Py = Py::None;
    let mut v_worst: Py = Py::None;
    let t258: Result<Option<Py>, PyErr> = (|| -> Result<Option<Py>, PyErr> {
            v_want = { let t257 = v_inputs.clone(); f_beta_diversity(py::getitem(&t257, &Py::str("x"))?, py::getitem(&t257, &Py::str("y"))?, py::method(&t257, "get", &[Py::str("similarity"), Py::Bool(false)])?)? };
            #[allow(unreachable_code)]
            Ok(None)
        })();
    match t258 {
        Ok(Some(r)) => { return Ok(r); }
        Ok(None) => {}
        Err(e) if (!matches!(e, PyErr::Unsupported(_))) => {
            return Ok(Py::Float(f64::INFINITY));
        }
        Err(e) => return Err(e),
    }
    v_worst = Py::Float(0.0);
    for t259 in py::iter(&(py::method(&v_want.clone(), "items", &[])?))? {
        let t296 = (t259).items()?;
        if t296.len() != 2 { return Err(PyErr::Value("unpack".into())); }
        v_k = t296[0].clone();
        v_v = t296[1].clone();
        if ({ let t260 = v_k.clone(); Py::Bool({ let t261 = v_outputs.clone(); if !py::contains(&t261, &t260)? { true } else { false } }) }).truthy() {
            return Ok(Py::Float(f64::INFINITY));
        }
        v_worst = py::minmax(&[v_worst.clone(), f__rel(py::getitem(&v_outputs.clone(), &v_k.clone())?, v_v.clone())?], true)?;
    }
    return Ok(v_worst.clone());
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn s_2764(mut v_inputs: Py, mut v_outputs: Py) -> R {
    let mut v_k: Py = Py::None;
    let mut v_v: Py = Py::None;
    let mut v_want: Py = Py::None;
    let mut v_worst: Py = Py::None;
    let t2766: Result<Option<Py>, PyErr> = (|| -> Result<Option<Py>, PyErr> {
            v_want = { let t2765 = v_inputs.clone(); f_hill_fit(py::getitem(&t2765, &Py::str("S"))?, py::getitem(&t2765, &Py::str("v"))?, py::method(&t2765, "get", &[Py::str("fixed_n"), Py::Bool(false)])?)? };
            #[allow(unreachable_code)]
            Ok(None)
        })();
    match t2766 {
        Ok(Some(r)) => { return Ok(r); }
        Ok(None) => {}
        Err(e) if (!matches!(e, PyErr::Unsupported(_))) => {
            return Ok(Py::Float(f64::INFINITY));
        }
        Err(e) => return Err(e),
    }
    v_worst = Py::Float(0.0);
    for t2767 in py::iter(&(py::method(&v_want.clone(), "items", &[])?))? {
        let t2804 = (t2767).items()?;
        if t2804.len() != 2 { return Err(PyErr::Value("unpack".into())); }
        v_k = t2804[0].clone();
        v_v = t2804[1].clone();
        if ({ let t2768 = v_k.clone(); Py::Bool({ let t2769 = v_outputs.clone(); if !py::contains(&t2769, &t2768)? { true } else { false } }) }).truthy() {
            return Ok(Py::Float(f64::INFINITY));
        }
        v_worst = py::minmax(&[v_worst.clone(), f__rel(py::getitem(&v_outputs.clone(), &v_k.clone())?, v_v.clone())?], true)?;
    }
    return Ok(v_worst.clone());
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn s_2860(mut v_inputs: Py, mut v_outputs: Py) -> R {
    let mut v_k: Py = Py::None;
    let mut v_v: Py = Py::None;
    let mut v_want: Py = Py::None;
    let mut v_worst: Py = Py::None;
    let t2862: Result<Option<Py>, PyErr> = (|| -> Result<Option<Py>, PyErr> {
            v_want = { let t2861 = v_inputs.clone(); f_hill_plot(py::getitem(&t2861, &Py::str("S"))?, py::getitem(&t2861, &Py::str("Y"))?)? };
            #[allow(unreachable_code)]
            Ok(None)
        })();
    match t2862 {
        Ok(Some(r)) => { return Ok(r); }
        Ok(None) => {}
        Err(e) if (!matches!(e, PyErr::Unsupported(_))) => {
            return Ok(Py::Float(f64::INFINITY));
        }
        Err(e) => return Err(e),
    }
    v_worst = Py::Float(0.0);
    for t2863 in py::iter(&(py::method(&v_want.clone(), "items", &[])?))? {
        let t2900 = (t2863).items()?;
        if t2900.len() != 2 { return Err(PyErr::Value("unpack".into())); }
        v_k = t2900[0].clone();
        v_v = t2900[1].clone();
        if ({ let t2864 = v_k.clone(); Py::Bool({ let t2865 = v_outputs.clone(); if !py::contains(&t2865, &t2864)? { true } else { false } }) }).truthy() {
            return Ok(Py::Float(f64::INFINITY));
        }
        v_worst = py::minmax(&[v_worst.clone(), f__rel(py::getitem(&v_outputs.clone(), &v_k.clone())?, v_v.clone())?], true)?;
    }
    return Ok(v_worst.clone());
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn s_2941(mut v_inputs: Py, mut v_outputs: Py) -> R {
    let mut v_k: Py = Py::None;
    let mut v_v: Py = Py::None;
    let mut v_want: Py = Py::None;
    let mut v_worst: Py = Py::None;
    let t2943: Result<Option<Py>, PyErr> = (|| -> Result<Option<Py>, PyErr> {
            v_want = { let t2942 = v_inputs.clone(); f_instability(py::getitem(&t2942, &Py::str("seq"))?)? };
            #[allow(unreachable_code)]
            Ok(None)
        })();
    match t2943 {
        Ok(Some(r)) => { return Ok(r); }
        Ok(None) => {}
        Err(e) if (!matches!(e, PyErr::Unsupported(_))) => {
            return Ok(Py::Float(f64::INFINITY));
        }
        Err(e) => return Err(e),
    }
    v_worst = Py::Float(0.0);
    for t2944 in py::iter(&(py::method(&v_want.clone(), "items", &[])?))? {
        let t2981 = (t2944).items()?;
        if t2981.len() != 2 { return Err(PyErr::Value("unpack".into())); }
        v_k = t2981[0].clone();
        v_v = t2981[1].clone();
        if ({ let t2945 = v_k.clone(); Py::Bool({ let t2946 = v_outputs.clone(); if !py::contains(&t2946, &t2945)? { true } else { false } }) }).truthy() {
            return Ok(Py::Float(f64::INFINITY));
        }
        v_worst = py::minmax(&[v_worst.clone(), f__rel(py::getitem(&v_outputs.clone(), &v_k.clone())?, v_v.clone())?], true)?;
    }
    return Ok(v_worst.clone());
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn s_3056(mut v_inputs: Py, mut v_outputs: Py) -> R {
    let mut v_k: Py = Py::None;
    let mut v_v: Py = Py::None;
    let mut v_want: Py = Py::None;
    let mut v_worst: Py = Py::None;
    let t3058: Result<Option<Py>, PyErr> = (|| -> Result<Option<Py>, PyErr> {
            v_want = { let t3057 = v_inputs.clone(); f_isoelectric_point(py::getitem(&t3057, &Py::str("seq"))?)? };
            #[allow(unreachable_code)]
            Ok(None)
        })();
    match t3058 {
        Ok(Some(r)) => { return Ok(r); }
        Ok(None) => {}
        Err(e) if (!matches!(e, PyErr::Unsupported(_))) => {
            return Ok(Py::Float(f64::INFINITY));
        }
        Err(e) => return Err(e),
    }
    v_worst = Py::Float(0.0);
    for t3059 in py::iter(&(py::method(&v_want.clone(), "items", &[])?))? {
        let t3096 = (t3059).items()?;
        if t3096.len() != 2 { return Err(PyErr::Value("unpack".into())); }
        v_k = t3096[0].clone();
        v_v = t3096[1].clone();
        if ({ let t3060 = v_k.clone(); Py::Bool({ let t3061 = v_outputs.clone(); if !py::contains(&t3061, &t3060)? { true } else { false } }) }).truthy() {
            return Ok(Py::Float(f64::INFINITY));
        }
        v_worst = py::minmax(&[v_worst.clone(), f__rel(py::getitem(&v_outputs.clone(), &v_k.clone())?, v_v.clone())?], true)?;
    }
    return Ok(v_worst.clone());
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn s_3144(mut v_inputs: Py, mut v_outputs: Py) -> R {
    let mut v_k: Py = Py::None;
    let mut v_v: Py = Py::None;
    let mut v_want: Py = Py::None;
    let mut v_worst: Py = Py::None;
    let t3146: Result<Option<Py>, PyErr> = (|| -> Result<Option<Py>, PyErr> {
            v_want = { let t3145 = v_inputs.clone(); f_jukes_cantor(py::getitem(&t3145, &Py::str("seq1"))?, py::getitem(&t3145, &Py::str("seq2"))?)? };
            #[allow(unreachable_code)]
            Ok(None)
        })();
    match t3146 {
        Ok(Some(r)) => { return Ok(r); }
        Ok(None) => {}
        Err(e) if (!matches!(e, PyErr::Unsupported(_))) => {
            return Ok(Py::Float(f64::INFINITY));
        }
        Err(e) => return Err(e),
    }
    v_worst = Py::Float(0.0);
    for t3147 in py::iter(&(py::method(&v_want.clone(), "items", &[])?))? {
        let t3184 = (t3147).items()?;
        if t3184.len() != 2 { return Err(PyErr::Value("unpack".into())); }
        v_k = t3184[0].clone();
        v_v = t3184[1].clone();
        if ({ let t3148 = v_k.clone(); Py::Bool({ let t3149 = v_outputs.clone(); if !py::contains(&t3149, &t3148)? { true } else { false } }) }).truthy() {
            return Ok(Py::Float(f64::INFINITY));
        }
        v_worst = py::minmax(&[v_worst.clone(), f__rel(py::getitem(&v_outputs.clone(), &v_k.clone())?, v_v.clone())?], true)?;
    }
    return Ok(v_worst.clone());
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn s_3251(mut v_inputs: Py, mut v_outputs: Py) -> R {
    let mut v_k: Py = Py::None;
    let mut v_v: Py = Py::None;
    let mut v_want: Py = Py::None;
    let mut v_worst: Py = Py::None;
    let t3253: Result<Option<Py>, PyErr> = (|| -> Result<Option<Py>, PyErr> {
            v_want = { let t3252 = v_inputs.clone(); f_kimura(py::getitem(&t3252, &Py::str("seq1"))?, py::getitem(&t3252, &Py::str("seq2"))?)? };
            #[allow(unreachable_code)]
            Ok(None)
        })();
    match t3253 {
        Ok(Some(r)) => { return Ok(r); }
        Ok(None) => {}
        Err(e) if (!matches!(e, PyErr::Unsupported(_))) => {
            return Ok(Py::Float(f64::INFINITY));
        }
        Err(e) => return Err(e),
    }
    v_worst = Py::Float(0.0);
    for t3254 in py::iter(&(py::method(&v_want.clone(), "items", &[])?))? {
        let t3291 = (t3254).items()?;
        if t3291.len() != 2 { return Err(PyErr::Value("unpack".into())); }
        v_k = t3291[0].clone();
        v_v = t3291[1].clone();
        if ({ let t3255 = v_k.clone(); Py::Bool({ let t3256 = v_outputs.clone(); if !py::contains(&t3256, &t3255)? { true } else { false } }) }).truthy() {
            return Ok(Py::Float(f64::INFINITY));
        }
        v_worst = py::minmax(&[v_worst.clone(), f__rel(py::getitem(&v_outputs.clone(), &v_k.clone())?, v_v.clone())?], true)?;
    }
    return Ok(v_worst.clone());
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn s_330(mut v_inputs: Py, mut v_outputs: Py) -> R {
    let mut v_k: Py = Py::None;
    let mut v_v: Py = Py::None;
    let mut v_want: Py = Py::None;
    let mut v_worst: Py = Py::None;
    let t332: Result<Option<Py>, PyErr> = (|| -> Result<Option<Py>, PyErr> {
            v_want = { let t331 = v_inputs.clone(); f_brillouin(py::getitem(&t331, &Py::str("counts"))?, py::method(&t331, "get", &[Py::str("log10"), Py::Bool(false)])?)? };
            #[allow(unreachable_code)]
            Ok(None)
        })();
    match t332 {
        Ok(Some(r)) => { return Ok(r); }
        Ok(None) => {}
        Err(e) if (!matches!(e, PyErr::Unsupported(_))) => {
            return Ok(Py::Float(f64::INFINITY));
        }
        Err(e) => return Err(e),
    }
    v_worst = Py::Float(0.0);
    for t333 in py::iter(&(py::method(&v_want.clone(), "items", &[])?))? {
        let t370 = (t333).items()?;
        if t370.len() != 2 { return Err(PyErr::Value("unpack".into())); }
        v_k = t370[0].clone();
        v_v = t370[1].clone();
        if ({ let t334 = v_k.clone(); Py::Bool({ let t335 = v_outputs.clone(); if !py::contains(&t335, &t334)? { true } else { false } }) }).truthy() {
            return Ok(Py::Float(f64::INFINITY));
        }
        v_worst = py::minmax(&[v_worst.clone(), f__rel(py::getitem(&v_outputs.clone(), &v_k.clone())?, v_v.clone())?], true)?;
    }
    return Ok(v_worst.clone());
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn s_3324(mut v_inputs: Py, mut v_outputs: Py) -> R {
    let mut v_k: Py = Py::None;
    let mut v_v: Py = Py::None;
    let mut v_want: Py = Py::None;
    let mut v_worst: Py = Py::None;
    let t3326: Result<Option<Py>, PyErr> = (|| -> Result<Option<Py>, PyErr> {
            v_want = { let t3325 = v_inputs.clone(); f_kmer_counts(py::getitem(&t3325, &Py::str("seq"))?, py::getitem(&t3325, &Py::str("k"))?)? };
            #[allow(unreachable_code)]
            Ok(None)
        })();
    match t3326 {
        Ok(Some(r)) => { return Ok(r); }
        Ok(None) => {}
        Err(e) if (!matches!(e, PyErr::Unsupported(_))) => {
            return Ok(Py::Float(f64::INFINITY));
        }
        Err(e) => return Err(e),
    }
    v_worst = Py::Float(0.0);
    for t3327 in py::iter(&(py::method(&v_want.clone(), "items", &[])?))? {
        let t3364 = (t3327).items()?;
        if t3364.len() != 2 { return Err(PyErr::Value("unpack".into())); }
        v_k = t3364[0].clone();
        v_v = t3364[1].clone();
        if ({ let t3328 = v_k.clone(); Py::Bool({ let t3329 = v_outputs.clone(); if !py::contains(&t3329, &t3328)? { true } else { false } }) }).truthy() {
            return Ok(Py::Float(f64::INFINITY));
        }
        v_worst = py::minmax(&[v_worst.clone(), f__rel(py::getitem(&v_outputs.clone(), &v_k.clone())?, v_v.clone())?], true)?;
    }
    return Ok(v_worst.clone());
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn s_3438(mut v_inputs: Py, mut v_outputs: Py) -> R {
    let mut v_k: Py = Py::None;
    let mut v_v: Py = Py::None;
    let mut v_want: Py = Py::None;
    let mut v_worst: Py = Py::None;
    let t3440: Result<Option<Py>, PyErr> = (|| -> Result<Option<Py>, PyErr> {
            v_want = { let t3439 = v_inputs.clone(); f_leslie_growth(py::getitem(&t3439, &Py::str("F"))?, py::getitem(&t3439, &Py::str("S"))?, py::method(&t3439, "get", &[Py::str("use_r0"), Py::Bool(false)])?)? };
            #[allow(unreachable_code)]
            Ok(None)
        })();
    match t3440 {
        Ok(Some(r)) => { return Ok(r); }
        Ok(None) => {}
        Err(e) if (!matches!(e, PyErr::Unsupported(_))) => {
            return Ok(Py::Float(f64::INFINITY));
        }
        Err(e) => return Err(e),
    }
    v_worst = Py::Float(0.0);
    for t3441 in py::iter(&(py::method(&v_want.clone(), "items", &[])?))? {
        let t3478 = (t3441).items()?;
        if t3478.len() != 2 { return Err(PyErr::Value("unpack".into())); }
        v_k = t3478[0].clone();
        v_v = t3478[1].clone();
        if ({ let t3442 = v_k.clone(); Py::Bool({ let t3443 = v_outputs.clone(); if !py::contains(&t3443, &t3442)? { true } else { false } }) }).truthy() {
            return Ok(Py::Float(f64::INFINITY));
        }
        v_worst = py::minmax(&[v_worst.clone(), f__rel(py::getitem(&v_outputs.clone(), &v_k.clone())?, v_v.clone())?], true)?;
    }
    return Ok(v_worst.clone());
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn s_3532(mut v_inputs: Py, mut v_outputs: Py) -> R {
    let mut v_k: Py = Py::None;
    let mut v_v: Py = Py::None;
    let mut v_want: Py = Py::None;
    let mut v_worst: Py = Py::None;
    let t3534: Result<Option<Py>, PyErr> = (|| -> Result<Option<Py>, PyErr> {
            v_want = { let t3533 = v_inputs.clone(); f_leslie_projection(py::getitem(&t3533, &Py::str("F"))?, py::getitem(&t3533, &Py::str("S"))?, py::getitem(&t3533, &Py::str("n0"))?, py::getitem(&t3533, &Py::str("t"))?, py::method(&t3533, "get", &[Py::str("transpose"), Py::Bool(false)])?)? };
            #[allow(unreachable_code)]
            Ok(None)
        })();
    match t3534 {
        Ok(Some(r)) => { return Ok(r); }
        Ok(None) => {}
        Err(e) if (!matches!(e, PyErr::Unsupported(_))) => {
            return Ok(Py::Float(f64::INFINITY));
        }
        Err(e) => return Err(e),
    }
    v_worst = Py::Float(0.0);
    for t3535 in py::iter(&(py::method(&v_want.clone(), "items", &[])?))? {
        let t3572 = (t3535).items()?;
        if t3572.len() != 2 { return Err(PyErr::Value("unpack".into())); }
        v_k = t3572[0].clone();
        v_v = t3572[1].clone();
        if ({ let t3536 = v_k.clone(); Py::Bool({ let t3537 = v_outputs.clone(); if !py::contains(&t3537, &t3536)? { true } else { false } }) }).truthy() {
            return Ok(Py::Float(f64::INFINITY));
        }
        v_worst = py::minmax(&[v_worst.clone(), f__rel(py::getitem(&v_outputs.clone(), &v_k.clone())?, v_v.clone())?], true)?;
    }
    return Ok(v_worst.clone());
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn s_36(mut v_inputs: Py, mut v_outputs: Py) -> R {
    let mut v_k: Py = Py::None;
    let mut v_v: Py = Py::None;
    let mut v_want: Py = Py::None;
    let mut v_worst: Py = Py::None;
    let t38: Result<Option<Py>, PyErr> = (|| -> Result<Option<Py>, PyErr> {
            v_want = { let t37 = v_inputs.clone(); f_ace(py::getitem(&t37, &Py::str("counts"))?, py::getitem(&t37, &Py::str("rare_threshold"))?, py::method(&t37, "get", &[Py::str("no_gamma"), Py::Bool(false)])?)? };
            #[allow(unreachable_code)]
            Ok(None)
        })();
    match t38 {
        Ok(Some(r)) => { return Ok(r); }
        Ok(None) => {}
        Err(e) if (!matches!(e, PyErr::Unsupported(_))) => {
            return Ok(Py::Float(f64::INFINITY));
        }
        Err(e) => return Err(e),
    }
    v_worst = Py::Float(0.0);
    for t39 in py::iter(&(py::method(&v_want.clone(), "items", &[])?))? {
        let t76 = (t39).items()?;
        if t76.len() != 2 { return Err(PyErr::Value("unpack".into())); }
        v_k = t76[0].clone();
        v_v = t76[1].clone();
        if ({ let t40 = v_k.clone(); Py::Bool({ let t41 = v_outputs.clone(); if !py::contains(&t41, &t40)? { true } else { false } }) }).truthy() {
            return Ok(Py::Float(f64::INFINITY));
        }
        v_worst = py::minmax(&[v_worst.clone(), f__rel(py::getitem(&v_outputs.clone(), &v_k.clone())?, v_v.clone())?], true)?;
    }
    return Ok(v_worst.clone());
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn s_3621(mut v_inputs: Py, mut v_outputs: Py) -> R {
    let mut v_k: Py = Py::None;
    let mut v_v: Py = Py::None;
    let mut v_want: Py = Py::None;
    let mut v_worst: Py = Py::None;
    let t3623: Result<Option<Py>, PyErr> = (|| -> Result<Option<Py>, PyErr> {
            v_want = { let t3622 = v_inputs.clone(); f_lineweaver_burk(py::getitem(&t3622, &Py::str("S"))?, py::getitem(&t3622, &Py::str("v"))?)? };
            #[allow(unreachable_code)]
            Ok(None)
        })();
    match t3623 {
        Ok(Some(r)) => { return Ok(r); }
        Ok(None) => {}
        Err(e) if (!matches!(e, PyErr::Unsupported(_))) => {
            return Ok(Py::Float(f64::INFINITY));
        }
        Err(e) => return Err(e),
    }
    v_worst = Py::Float(0.0);
    for t3624 in py::iter(&(py::method(&v_want.clone(), "items", &[])?))? {
        let t3661 = (t3624).items()?;
        if t3661.len() != 2 { return Err(PyErr::Value("unpack".into())); }
        v_k = t3661[0].clone();
        v_v = t3661[1].clone();
        if ({ let t3625 = v_k.clone(); Py::Bool({ let t3626 = v_outputs.clone(); if !py::contains(&t3626, &t3625)? { true } else { false } }) }).truthy() {
            return Ok(Py::Float(f64::INFINITY));
        }
        v_worst = py::minmax(&[v_worst.clone(), f__rel(py::getitem(&v_outputs.clone(), &v_k.clone())?, v_v.clone())?], true)?;
    }
    return Ok(v_worst.clone());
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn s_3680(mut v_inputs: Py, mut v_outputs: Py) -> R {
    let mut v_k: Py = Py::None;
    let mut v_v: Py = Py::None;
    let mut v_want: Py = Py::None;
    let mut v_worst: Py = Py::None;
    let t3682: Result<Option<Py>, PyErr> = (|| -> Result<Option<Py>, PyErr> {
            v_want = { let t3681 = v_inputs.clone(); f_linkage_disequilibrium(py::getitem(&t3681, &Py::str("n_AB"))?, py::getitem(&t3681, &Py::str("n_Ab"))?, py::getitem(&t3681, &Py::str("n_aB"))?, py::getitem(&t3681, &Py::str("n_ab"))?)? };
            #[allow(unreachable_code)]
            Ok(None)
        })();
    match t3682 {
        Ok(Some(r)) => { return Ok(r); }
        Ok(None) => {}
        Err(e) if (!matches!(e, PyErr::Unsupported(_))) => {
            return Ok(Py::Float(f64::INFINITY));
        }
        Err(e) => return Err(e),
    }
    v_worst = Py::Float(0.0);
    for t3683 in py::iter(&(py::method(&v_want.clone(), "items", &[])?))? {
        let t3720 = (t3683).items()?;
        if t3720.len() != 2 { return Err(PyErr::Value("unpack".into())); }
        v_k = t3720[0].clone();
        v_v = t3720[1].clone();
        if ({ let t3684 = v_k.clone(); Py::Bool({ let t3685 = v_outputs.clone(); if !py::contains(&t3685, &t3684)? { true } else { false } }) }).truthy() {
            return Ok(Py::Float(f64::INFINITY));
        }
        v_worst = py::minmax(&[v_worst.clone(), f__rel(py::getitem(&v_outputs.clone(), &v_k.clone())?, v_v.clone())?], true)?;
    }
    return Ok(v_worst.clone());
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn s_3804(mut v_inputs: Py, mut v_outputs: Py) -> R {
    let mut v_k: Py = Py::None;
    let mut v_v: Py = Py::None;
    let mut v_want: Py = Py::None;
    let mut v_worst: Py = Py::None;
    let t3806: Result<Option<Py>, PyErr> = (|| -> Result<Option<Py>, PyErr> {
            v_want = { let t3805 = v_inputs.clone(); f_local_affine(py::getitem(&t3805, &Py::str("seq1"))?, py::getitem(&t3805, &Py::str("seq2"))?, py::getitem(&t3805, &Py::str("match"))?, py::getitem(&t3805, &Py::str("mismatch"))?, py::getitem(&t3805, &Py::str("open"))?, py::getitem(&t3805, &Py::str("extend"))?)? };
            #[allow(unreachable_code)]
            Ok(None)
        })();
    match t3806 {
        Ok(Some(r)) => { return Ok(r); }
        Ok(None) => {}
        Err(e) if (!matches!(e, PyErr::Unsupported(_))) => {
            return Ok(Py::Float(f64::INFINITY));
        }
        Err(e) => return Err(e),
    }
    v_worst = Py::Float(0.0);
    for t3807 in py::iter(&(py::method(&v_want.clone(), "items", &[])?))? {
        let t3844 = (t3807).items()?;
        if t3844.len() != 2 { return Err(PyErr::Value("unpack".into())); }
        v_k = t3844[0].clone();
        v_v = t3844[1].clone();
        if ({ let t3808 = v_k.clone(); Py::Bool({ let t3809 = v_outputs.clone(); if !py::contains(&t3809, &t3808)? { true } else { false } }) }).truthy() {
            return Ok(Py::Float(f64::INFINITY));
        }
        v_worst = py::minmax(&[v_worst.clone(), f__rel(py::getitem(&v_outputs.clone(), &v_k.clone())?, v_v.clone())?], true)?;
    }
    return Ok(v_worst.clone());
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn s_3925(mut v_inputs: Py, mut v_outputs: Py) -> R {
    let mut v_k: Py = Py::None;
    let mut v_v: Py = Py::None;
    let mut v_want: Py = Py::None;
    let mut v_worst: Py = Py::None;
    let t3927: Result<Option<Py>, PyErr> = (|| -> Result<Option<Py>, PyErr> {
            v_want = { let t3926 = v_inputs.clone(); f_local_linear(py::getitem(&t3926, &Py::str("seq1"))?, py::getitem(&t3926, &Py::str("seq2"))?, py::getitem(&t3926, &Py::str("match"))?, py::getitem(&t3926, &Py::str("mismatch"))?, py::getitem(&t3926, &Py::str("gap"))?)? };
            #[allow(unreachable_code)]
            Ok(None)
        })();
    match t3927 {
        Ok(Some(r)) => { return Ok(r); }
        Ok(None) => {}
        Err(e) if (!matches!(e, PyErr::Unsupported(_))) => {
            return Ok(Py::Float(f64::INFINITY));
        }
        Err(e) => return Err(e),
    }
    v_worst = Py::Float(0.0);
    for t3928 in py::iter(&(py::method(&v_want.clone(), "items", &[])?))? {
        let t3965 = (t3928).items()?;
        if t3965.len() != 2 { return Err(PyErr::Value("unpack".into())); }
        v_k = t3965[0].clone();
        v_v = t3965[1].clone();
        if ({ let t3929 = v_k.clone(); Py::Bool({ let t3930 = v_outputs.clone(); if !py::contains(&t3930, &t3929)? { true } else { false } }) }).truthy() {
            return Ok(Py::Float(f64::INFINITY));
        }
        v_worst = py::minmax(&[v_worst.clone(), f__rel(py::getitem(&v_outputs.clone(), &v_k.clone())?, v_v.clone())?], true)?;
    }
    return Ok(v_worst.clone());
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn s_3997(mut v_inputs: Py, mut v_outputs: Py) -> R {
    let mut v_k: Py = Py::None;
    let mut v_v: Py = Py::None;
    let mut v_want: Py = Py::None;
    let mut v_worst: Py = Py::None;
    let t3999: Result<Option<Py>, PyErr> = (|| -> Result<Option<Py>, PyErr> {
            v_want = { let t3998 = v_inputs.clone(); f_longest_orf(py::getitem(&t3998, &Py::str("seq"))?)? };
            #[allow(unreachable_code)]
            Ok(None)
        })();
    match t3999 {
        Ok(Some(r)) => { return Ok(r); }
        Ok(None) => {}
        Err(e) if (!matches!(e, PyErr::Unsupported(_))) => {
            return Ok(Py::Float(f64::INFINITY));
        }
        Err(e) => return Err(e),
    }
    v_worst = Py::Float(0.0);
    for t4000 in py::iter(&(py::method(&v_want.clone(), "items", &[])?))? {
        let t4037 = (t4000).items()?;
        if t4037.len() != 2 { return Err(PyErr::Value("unpack".into())); }
        v_k = t4037[0].clone();
        v_v = t4037[1].clone();
        if ({ let t4001 = v_k.clone(); Py::Bool({ let t4002 = v_outputs.clone(); if !py::contains(&t4002, &t4001)? { true } else { false } }) }).truthy() {
            return Ok(Py::Float(f64::INFINITY));
        }
        v_worst = py::minmax(&[v_worst.clone(), f__rel(py::getitem(&v_outputs.clone(), &v_k.clone())?, v_v.clone())?], true)?;
    }
    return Ok(v_worst.clone());
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn s_404(mut v_inputs: Py, mut v_outputs: Py) -> R {
    let mut v_k: Py = Py::None;
    let mut v_v: Py = Py::None;
    let mut v_want: Py = Py::None;
    let mut v_worst: Py = Py::None;
    let t406: Result<Option<Py>, PyErr> = (|| -> Result<Option<Py>, PyErr> {
            v_want = { let t405 = v_inputs.clone(); f_chao1(py::getitem(&t405, &Py::str("counts"))?, py::getitem(&t405, &Py::str("form"))?, py::method(&t405, "get", &[Py::str("swap"), Py::Bool(false)])?)? };
            #[allow(unreachable_code)]
            Ok(None)
        })();
    match t406 {
        Ok(Some(r)) => { return Ok(r); }
        Ok(None) => {}
        Err(e) if (!matches!(e, PyErr::Unsupported(_))) => {
            return Ok(Py::Float(f64::INFINITY));
        }
        Err(e) => return Err(e),
    }
    v_worst = Py::Float(0.0);
    for t407 in py::iter(&(py::method(&v_want.clone(), "items", &[])?))? {
        let t444 = (t407).items()?;
        if t444.len() != 2 { return Err(PyErr::Value("unpack".into())); }
        v_k = t444[0].clone();
        v_v = t444[1].clone();
        if ({ let t408 = v_k.clone(); Py::Bool({ let t409 = v_outputs.clone(); if !py::contains(&t409, &t408)? { true } else { false } }) }).truthy() {
            return Ok(Py::Float(f64::INFINITY));
        }
        v_worst = py::minmax(&[v_worst.clone(), f__rel(py::getitem(&v_outputs.clone(), &v_k.clone())?, v_v.clone())?], true)?;
    }
    return Ok(v_worst.clone());
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn s_4062(mut v_inputs: Py, mut v_outputs: Py) -> R {
    let mut v_k: Py = Py::None;
    let mut v_v: Py = Py::None;
    let mut v_want: Py = Py::None;
    let mut v_worst: Py = Py::None;
    let t4064: Result<Option<Py>, PyErr> = (|| -> Result<Option<Py>, PyErr> {
            v_want = { let t4063 = v_inputs.clone(); f_map_distance(py::getitem(&t4063, &Py::str("r"))?, py::getitem(&t4063, &Py::str("map"))?)? };
            #[allow(unreachable_code)]
            Ok(None)
        })();
    match t4064 {
        Ok(Some(r)) => { return Ok(r); }
        Ok(None) => {}
        Err(e) if (!matches!(e, PyErr::Unsupported(_))) => {
            return Ok(Py::Float(f64::INFINITY));
        }
        Err(e) => return Err(e),
    }
    v_worst = Py::Float(0.0);
    for t4065 in py::iter(&(py::method(&v_want.clone(), "items", &[])?))? {
        let t4102 = (t4065).items()?;
        if t4102.len() != 2 { return Err(PyErr::Value("unpack".into())); }
        v_k = t4102[0].clone();
        v_v = t4102[1].clone();
        if ({ let t4066 = v_k.clone(); Py::Bool({ let t4067 = v_outputs.clone(); if !py::contains(&t4067, &t4066)? { true } else { false } }) }).truthy() {
            return Ok(Py::Float(f64::INFINITY));
        }
        v_worst = py::minmax(&[v_worst.clone(), f__rel(py::getitem(&v_outputs.clone(), &v_k.clone())?, v_v.clone())?], true)?;
    }
    return Ok(v_worst.clone());
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn s_4179(mut v_inputs: Py, mut v_outputs: Py) -> R {
    let mut v_k: Py = Py::None;
    let mut v_v: Py = Py::None;
    let mut v_want: Py = Py::None;
    let mut v_worst: Py = Py::None;
    let t4181: Result<Option<Py>, PyErr> = (|| -> Result<Option<Py>, PyErr> {
            v_want = { let t4180 = v_inputs.clone(); f_mm_fit(py::getitem(&t4180, &Py::str("S"))?, py::getitem(&t4180, &Py::str("v"))?)? };
            #[allow(unreachable_code)]
            Ok(None)
        })();
    match t4181 {
        Ok(Some(r)) => { return Ok(r); }
        Ok(None) => {}
        Err(e) if (!matches!(e, PyErr::Unsupported(_))) => {
            return Ok(Py::Float(f64::INFINITY));
        }
        Err(e) => return Err(e),
    }
    v_worst = Py::Float(0.0);
    for t4182 in py::iter(&(py::method(&v_want.clone(), "items", &[])?))? {
        let t4219 = (t4182).items()?;
        if t4219.len() != 2 { return Err(PyErr::Value("unpack".into())); }
        v_k = t4219[0].clone();
        v_v = t4219[1].clone();
        if ({ let t4183 = v_k.clone(); Py::Bool({ let t4184 = v_outputs.clone(); if !py::contains(&t4184, &t4183)? { true } else { false } }) }).truthy() {
            return Ok(Py::Float(f64::INFINITY));
        }
        v_worst = py::minmax(&[v_worst.clone(), f__rel(py::getitem(&v_outputs.clone(), &v_k.clone())?, v_v.clone())?], true)?;
    }
    return Ok(v_worst.clone());
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn s_4277(mut v_inputs: Py, mut v_outputs: Py) -> R {
    let mut v_k: Py = Py::None;
    let mut v_v: Py = Py::None;
    let mut v_want: Py = Py::None;
    let mut v_worst: Py = Py::None;
    let t4279: Result<Option<Py>, PyErr> = (|| -> Result<Option<Py>, PyErr> {
            v_want = { let t4278 = v_inputs.clone(); f_mm_progress(py::getitem(&t4278, &Py::str("S0"))?, py::getitem(&t4278, &Py::str("Vmax"))?, py::getitem(&t4278, &Py::str("Km"))?, py::getitem(&t4278, &Py::str("t"))?, py::method(&t4278, "get", &[Py::str("first_order"), Py::Bool(false)])?)? };
            #[allow(unreachable_code)]
            Ok(None)
        })();
    match t4279 {
        Ok(Some(r)) => { return Ok(r); }
        Ok(None) => {}
        Err(e) if (!matches!(e, PyErr::Unsupported(_))) => {
            return Ok(Py::Float(f64::INFINITY));
        }
        Err(e) => return Err(e),
    }
    v_worst = Py::Float(0.0);
    for t4280 in py::iter(&(py::method(&v_want.clone(), "items", &[])?))? {
        let t4317 = (t4280).items()?;
        if t4317.len() != 2 { return Err(PyErr::Value("unpack".into())); }
        v_k = t4317[0].clone();
        v_v = t4317[1].clone();
        if ({ let t4281 = v_k.clone(); Py::Bool({ let t4282 = v_outputs.clone(); if !py::contains(&t4282, &t4281)? { true } else { false } }) }).truthy() {
            return Ok(Py::Float(f64::INFINITY));
        }
        v_worst = py::minmax(&[v_worst.clone(), f__rel(py::getitem(&v_outputs.clone(), &v_k.clone())?, v_v.clone())?], true)?;
    }
    return Ok(v_worst.clone());
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn s_4344(mut v_inputs: Py, mut v_outputs: Py) -> R {
    let mut v_k: Py = Py::None;
    let mut v_v: Py = Py::None;
    let mut v_want: Py = Py::None;
    let mut v_worst: Py = Py::None;
    let t4346: Result<Option<Py>, PyErr> = (|| -> Result<Option<Py>, PyErr> {
            v_want = { let t4345 = v_inputs.clone(); f_mm_time(py::getitem(&t4345, &Py::str("S0"))?, py::getitem(&t4345, &Py::str("Vmax"))?, py::getitem(&t4345, &Py::str("Km"))?, py::getitem(&t4345, &Py::str("X"))?, py::method(&t4345, "get", &[Py::str("zero_order"), Py::Bool(false)])?)? };
            #[allow(unreachable_code)]
            Ok(None)
        })();
    match t4346 {
        Ok(Some(r)) => { return Ok(r); }
        Ok(None) => {}
        Err(e) if (!matches!(e, PyErr::Unsupported(_))) => {
            return Ok(Py::Float(f64::INFINITY));
        }
        Err(e) => return Err(e),
    }
    v_worst = Py::Float(0.0);
    for t4347 in py::iter(&(py::method(&v_want.clone(), "items", &[])?))? {
        let t4384 = (t4347).items()?;
        if t4384.len() != 2 { return Err(PyErr::Value("unpack".into())); }
        v_k = t4384[0].clone();
        v_v = t4384[1].clone();
        if ({ let t4348 = v_k.clone(); Py::Bool({ let t4349 = v_outputs.clone(); if !py::contains(&t4349, &t4348)? { true } else { false } }) }).truthy() {
            return Ok(Py::Float(f64::INFINITY));
        }
        v_worst = py::minmax(&[v_worst.clone(), f__rel(py::getitem(&v_outputs.clone(), &v_k.clone())?, v_v.clone())?], true)?;
    }
    return Ok(v_worst.clone());
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn s_4419(mut v_inputs: Py, mut v_outputs: Py) -> R {
    let mut v_k: Py = Py::None;
    let mut v_v: Py = Py::None;
    let mut v_want: Py = Py::None;
    let mut v_worst: Py = Py::None;
    let t4421: Result<Option<Py>, PyErr> = (|| -> Result<Option<Py>, PyErr> {
            v_want = { let t4420 = v_inputs.clone(); f_moran_fixation(py::getitem(&t4420, &Py::str("N"))?, py::getitem(&t4420, &Py::str("i"))?, py::getitem(&t4420, &Py::str("r"))?)? };
            #[allow(unreachable_code)]
            Ok(None)
        })();
    match t4421 {
        Ok(Some(r)) => { return Ok(r); }
        Ok(None) => {}
        Err(e) if (!matches!(e, PyErr::Unsupported(_))) => {
            return Ok(Py::Float(f64::INFINITY));
        }
        Err(e) => return Err(e),
    }
    v_worst = Py::Float(0.0);
    for t4422 in py::iter(&(py::method(&v_want.clone(), "items", &[])?))? {
        let t4459 = (t4422).items()?;
        if t4459.len() != 2 { return Err(PyErr::Value("unpack".into())); }
        v_k = t4459[0].clone();
        v_v = t4459[1].clone();
        if ({ let t4423 = v_k.clone(); Py::Bool({ let t4424 = v_outputs.clone(); if !py::contains(&t4424, &t4423)? { true } else { false } }) }).truthy() {
            return Ok(Py::Float(f64::INFINITY));
        }
        v_worst = py::minmax(&[v_worst.clone(), f__rel(py::getitem(&v_outputs.clone(), &v_k.clone())?, v_v.clone())?], true)?;
    }
    return Ok(v_worst.clone());
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn s_4516(mut v_inputs: Py, mut v_outputs: Py) -> R {
    let mut v_k: Py = Py::None;
    let mut v_v: Py = Py::None;
    let mut v_want: Py = Py::None;
    let mut v_worst: Py = Py::None;
    let t4518: Result<Option<Py>, PyErr> = (|| -> Result<Option<Py>, PyErr> {
            v_want = { let t4517 = v_inputs.clone(); f_mutation_selection(py::getitem(&t4517, &Py::str("mu"))?, py::getitem(&t4517, &Py::str("s"))?, py::getitem(&t4517, &Py::str("h"))?)? };
            #[allow(unreachable_code)]
            Ok(None)
        })();
    match t4518 {
        Ok(Some(r)) => { return Ok(r); }
        Ok(None) => {}
        Err(e) if (!matches!(e, PyErr::Unsupported(_))) => {
            return Ok(Py::Float(f64::INFINITY));
        }
        Err(e) => return Err(e),
    }
    v_worst = Py::Float(0.0);
    for t4519 in py::iter(&(py::method(&v_want.clone(), "items", &[])?))? {
        let t4556 = (t4519).items()?;
        if t4556.len() != 2 { return Err(PyErr::Value("unpack".into())); }
        v_k = t4556[0].clone();
        v_v = t4556[1].clone();
        if ({ let t4520 = v_k.clone(); Py::Bool({ let t4521 = v_outputs.clone(); if !py::contains(&t4521, &t4520)? { true } else { false } }) }).truthy() {
            return Ok(Py::Float(f64::INFINITY));
        }
        v_worst = py::minmax(&[v_worst.clone(), f__rel(py::getitem(&v_outputs.clone(), &v_k.clone())?, v_v.clone())?], true)?;
    }
    return Ok(v_worst.clone());
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn s_4585(mut v_inputs: Py, mut v_outputs: Py) -> R {
    let mut v_k: Py = Py::None;
    let mut v_v: Py = Py::None;
    let mut v_want: Py = Py::None;
    let mut v_worst: Py = Py::None;
    let t4587: Result<Option<Py>, PyErr> = (|| -> Result<Option<Py>, PyErr> {
            v_want = { let t4586 = v_inputs.clone(); f_ne_fluctuating(py::getitem(&t4586, &Py::str("sizes"))?)? };
            #[allow(unreachable_code)]
            Ok(None)
        })();
    match t4587 {
        Ok(Some(r)) => { return Ok(r); }
        Ok(None) => {}
        Err(e) if (!matches!(e, PyErr::Unsupported(_))) => {
            return Ok(Py::Float(f64::INFINITY));
        }
        Err(e) => return Err(e),
    }
    v_worst = Py::Float(0.0);
    for t4588 in py::iter(&(py::method(&v_want.clone(), "items", &[])?))? {
        let t4625 = (t4588).items()?;
        if t4625.len() != 2 { return Err(PyErr::Value("unpack".into())); }
        v_k = t4625[0].clone();
        v_v = t4625[1].clone();
        if ({ let t4589 = v_k.clone(); Py::Bool({ let t4590 = v_outputs.clone(); if !py::contains(&t4590, &t4589)? { true } else { false } }) }).truthy() {
            return Ok(Py::Float(f64::INFINITY));
        }
        v_worst = py::minmax(&[v_worst.clone(), f__rel(py::getitem(&v_outputs.clone(), &v_k.clone())?, v_v.clone())?], true)?;
    }
    return Ok(v_worst.clone());
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn s_4633(mut v_inputs: Py, mut v_outputs: Py) -> R {
    let mut v_k: Py = Py::None;
    let mut v_v: Py = Py::None;
    let mut v_want: Py = Py::None;
    let mut v_worst: Py = Py::None;
    let t4635: Result<Option<Py>, PyErr> = (|| -> Result<Option<Py>, PyErr> {
            v_want = { let t4634 = v_inputs.clone(); f_ne_sex_ratio(py::getitem(&t4634, &Py::str("N_m"))?, py::getitem(&t4634, &Py::str("N_f"))?)? };
            #[allow(unreachable_code)]
            Ok(None)
        })();
    match t4635 {
        Ok(Some(r)) => { return Ok(r); }
        Ok(None) => {}
        Err(e) if (!matches!(e, PyErr::Unsupported(_))) => {
            return Ok(Py::Float(f64::INFINITY));
        }
        Err(e) => return Err(e),
    }
    v_worst = Py::Float(0.0);
    for t4636 in py::iter(&(py::method(&v_want.clone(), "items", &[])?))? {
        let t4673 = (t4636).items()?;
        if t4673.len() != 2 { return Err(PyErr::Value("unpack".into())); }
        v_k = t4673[0].clone();
        v_v = t4673[1].clone();
        if ({ let t4637 = v_k.clone(); Py::Bool({ let t4638 = v_outputs.clone(); if !py::contains(&t4638, &t4637)? { true } else { false } }) }).truthy() {
            return Ok(Py::Float(f64::INFINITY));
        }
        v_worst = py::minmax(&[v_worst.clone(), f__rel(py::getitem(&v_outputs.clone(), &v_k.clone())?, v_v.clone())?], true)?;
    }
    return Ok(v_worst.clone());
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn s_469(mut v_inputs: Py, mut v_outputs: Py) -> R {
    let mut v_k: Py = Py::None;
    let mut v_v: Py = Py::None;
    let mut v_want: Py = Py::None;
    let mut v_worst: Py = Py::None;
    let t471: Result<Option<Py>, PyErr> = (|| -> Result<Option<Py>, PyErr> {
            v_want = { let t470 = v_inputs.clone(); f_chemostat(py::getitem(&t470, &Py::str("mu_max"))?, py::getitem(&t470, &Py::str("Ks"))?, py::getitem(&t470, &Py::str("Y"))?, py::getitem(&t470, &Py::str("D"))?, py::getitem(&t470, &Py::str("S0"))?, py::method(&t470, "get", &[Py::str("no_residual"), Py::Bool(false)])?)? };
            #[allow(unreachable_code)]
            Ok(None)
        })();
    match t471 {
        Ok(Some(r)) => { return Ok(r); }
        Ok(None) => {}
        Err(e) if (!matches!(e, PyErr::Unsupported(_))) => {
            return Ok(Py::Float(f64::INFINITY));
        }
        Err(e) => return Err(e),
    }
    v_worst = Py::Float(0.0);
    for t472 in py::iter(&(py::method(&v_want.clone(), "items", &[])?))? {
        let t509 = (t472).items()?;
        if t509.len() != 2 { return Err(PyErr::Value("unpack".into())); }
        v_k = t509[0].clone();
        v_v = t509[1].clone();
        if ({ let t473 = v_k.clone(); Py::Bool({ let t474 = v_outputs.clone(); if !py::contains(&t474, &t473)? { true } else { false } }) }).truthy() {
            return Ok(Py::Float(f64::INFINITY));
        }
        v_worst = py::minmax(&[v_worst.clone(), f__rel(py::getitem(&v_outputs.clone(), &v_k.clone())?, v_v.clone())?], true)?;
    }
    return Ok(v_worst.clone());
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn s_4753(mut v_inputs: Py, mut v_outputs: Py) -> R {
    let mut v_k: Py = Py::None;
    let mut v_v: Py = Py::None;
    let mut v_want: Py = Py::None;
    let mut v_worst: Py = Py::None;
    let t4755: Result<Option<Py>, PyErr> = (|| -> Result<Option<Py>, PyErr> {
            v_want = { let t4754 = v_inputs.clone(); f_neighbor_joining(py::getitem(&t4754, &Py::str("D"))?, py::method(&t4754, "get", &[Py::str("corrected"), Py::Bool(true)])?)? };
            #[allow(unreachable_code)]
            Ok(None)
        })();
    match t4755 {
        Ok(Some(r)) => { return Ok(r); }
        Ok(None) => {}
        Err(e) if (!matches!(e, PyErr::Unsupported(_))) => {
            return Ok(Py::Float(f64::INFINITY));
        }
        Err(e) => return Err(e),
    }
    v_worst = Py::Float(0.0);
    for t4756 in py::iter(&(py::method(&v_want.clone(), "items", &[])?))? {
        let t4793 = (t4756).items()?;
        if t4793.len() != 2 { return Err(PyErr::Value("unpack".into())); }
        v_k = t4793[0].clone();
        v_v = t4793[1].clone();
        if ({ let t4757 = v_k.clone(); Py::Bool({ let t4758 = v_outputs.clone(); if !py::contains(&t4758, &t4757)? { true } else { false } }) }).truthy() {
            return Ok(Py::Float(f64::INFINITY));
        }
        v_worst = py::minmax(&[v_worst.clone(), f__rel(py::getitem(&v_outputs.clone(), &v_k.clone())?, v_v.clone())?], true)?;
    }
    return Ok(v_worst.clone());
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn s_4842(mut v_inputs: Py, mut v_outputs: Py) -> R {
    let mut v_k: Py = Py::None;
    let mut v_v: Py = Py::None;
    let mut v_want: Py = Py::None;
    let mut v_worst: Py = Py::None;
    let t4844: Result<Option<Py>, PyErr> = (|| -> Result<Option<Py>, PyErr> {
            v_want = { let t4843 = v_inputs.clone(); f_nucleic_mass(py::getitem(&t4843, &Py::str("seq"))?, py::getitem(&t4843, &Py::str("molecule"))?, py::getitem(&t4843, &Py::str("strands"))?)? };
            #[allow(unreachable_code)]
            Ok(None)
        })();
    match t4844 {
        Ok(Some(r)) => { return Ok(r); }
        Ok(None) => {}
        Err(e) if (!matches!(e, PyErr::Unsupported(_))) => {
            return Ok(Py::Float(f64::INFINITY));
        }
        Err(e) => return Err(e),
    }
    v_worst = Py::Float(0.0);
    for t4845 in py::iter(&(py::method(&v_want.clone(), "items", &[])?))? {
        let t4882 = (t4845).items()?;
        if t4882.len() != 2 { return Err(PyErr::Value("unpack".into())); }
        v_k = t4882[0].clone();
        v_v = t4882[1].clone();
        if ({ let t4846 = v_k.clone(); Py::Bool({ let t4847 = v_outputs.clone(); if !py::contains(&t4847, &t4846)? { true } else { false } }) }).truthy() {
            return Ok(Py::Float(f64::INFINITY));
        }
        v_worst = py::minmax(&[v_worst.clone(), f__rel(py::getitem(&v_outputs.clone(), &v_k.clone())?, v_v.clone())?], true)?;
    }
    return Ok(v_worst.clone());
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn s_4898(mut v_inputs: Py, mut v_outputs: Py) -> R {
    let mut v_k: Py = Py::None;
    let mut v_v: Py = Py::None;
    let mut v_want: Py = Py::None;
    let mut v_worst: Py = Py::None;
    let t4900: Result<Option<Py>, PyErr> = (|| -> Result<Option<Py>, PyErr> {
            v_want = { let t4899 = v_inputs.clone(); f_nucleotide_diversity(py::getitem(&t4899, &Py::str("seqs"))?)? };
            #[allow(unreachable_code)]
            Ok(None)
        })();
    match t4900 {
        Ok(Some(r)) => { return Ok(r); }
        Ok(None) => {}
        Err(e) if (!matches!(e, PyErr::Unsupported(_))) => {
            return Ok(Py::Float(f64::INFINITY));
        }
        Err(e) => return Err(e),
    }
    v_worst = Py::Float(0.0);
    for t4901 in py::iter(&(py::method(&v_want.clone(), "items", &[])?))? {
        let t4938 = (t4901).items()?;
        if t4938.len() != 2 { return Err(PyErr::Value("unpack".into())); }
        v_k = t4938[0].clone();
        v_v = t4938[1].clone();
        if ({ let t4902 = v_k.clone(); Py::Bool({ let t4903 = v_outputs.clone(); if !py::contains(&t4903, &t4902)? { true } else { false } }) }).truthy() {
            return Ok(Py::Float(f64::INFINITY));
        }
        v_worst = py::minmax(&[v_worst.clone(), f__rel(py::getitem(&v_outputs.clone(), &v_k.clone())?, v_v.clone())?], true)?;
    }
    return Ok(v_worst.clone());
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn s_4980(mut v_inputs: Py, mut v_outputs: Py) -> R {
    let mut v_k: Py = Py::None;
    let mut v_v: Py = Py::None;
    let mut v_want: Py = Py::None;
    let mut v_worst: Py = Py::None;
    let t4982: Result<Option<Py>, PyErr> = (|| -> Result<Option<Py>, PyErr> {
            v_want = { let t4981 = v_inputs.clone(); f_pedigree_inbreeding(py::getitem(&t4981, &Py::str("pedigree"))?, py::getitem(&t4981, &Py::str("individual"))?)? };
            #[allow(unreachable_code)]
            Ok(None)
        })();
    match t4982 {
        Ok(Some(r)) => { return Ok(r); }
        Ok(None) => {}
        Err(e) if (!matches!(e, PyErr::Unsupported(_))) => {
            return Ok(Py::Float(f64::INFINITY));
        }
        Err(e) => return Err(e),
    }
    v_worst = Py::Float(0.0);
    for t4983 in py::iter(&(py::method(&v_want.clone(), "items", &[])?))? {
        let t5020 = (t4983).items()?;
        if t5020.len() != 2 { return Err(PyErr::Value("unpack".into())); }
        v_k = t5020[0].clone();
        v_v = t5020[1].clone();
        if ({ let t4984 = v_k.clone(); Py::Bool({ let t4985 = v_outputs.clone(); if !py::contains(&t4985, &t4984)? { true } else { false } }) }).truthy() {
            return Ok(Py::Float(f64::INFINITY));
        }
        v_worst = py::minmax(&[v_worst.clone(), f__rel(py::getitem(&v_outputs.clone(), &v_k.clone())?, v_v.clone())?], true)?;
    }
    return Ok(v_worst.clone());
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn s_5028(mut v_inputs: Py, mut v_outputs: Py) -> R {
    let mut v_k: Py = Py::None;
    let mut v_v: Py = Py::None;
    let mut v_want: Py = Py::None;
    let mut v_worst: Py = Py::None;
    let t5030: Result<Option<Py>, PyErr> = (|| -> Result<Option<Py>, PyErr> {
            v_want = { let t5029 = v_inputs.clone(); f_petersen_mle(py::getitem(&t5029, &Py::str("M"))?, py::getitem(&t5029, &Py::str("C"))?, py::getitem(&t5029, &Py::str("R"))?, py::method(&t5029, "get", &[Py::str("ceil"), Py::Bool(false)])?)? };
            #[allow(unreachable_code)]
            Ok(None)
        })();
    match t5030 {
        Ok(Some(r)) => { return Ok(r); }
        Ok(None) => {}
        Err(e) if (!matches!(e, PyErr::Unsupported(_))) => {
            return Ok(Py::Float(f64::INFINITY));
        }
        Err(e) => return Err(e),
    }
    v_worst = Py::Float(0.0);
    for t5031 in py::iter(&(py::method(&v_want.clone(), "items", &[])?))? {
        let t5068 = (t5031).items()?;
        if t5068.len() != 2 { return Err(PyErr::Value("unpack".into())); }
        v_k = t5068[0].clone();
        v_v = t5068[1].clone();
        if ({ let t5032 = v_k.clone(); Py::Bool({ let t5033 = v_outputs.clone(); if !py::contains(&t5033, &t5032)? { true } else { false } }) }).truthy() {
            return Ok(Py::Float(f64::INFINITY));
        }
        v_worst = py::minmax(&[v_worst.clone(), f__rel(py::getitem(&v_outputs.clone(), &v_k.clone())?, v_v.clone())?], true)?;
    }
    return Ok(v_worst.clone());
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn s_5107(mut v_inputs: Py, mut v_outputs: Py) -> R {
    let mut v_k: Py = Py::None;
    let mut v_v: Py = Py::None;
    let mut v_want: Py = Py::None;
    let mut v_worst: Py = Py::None;
    let t5109: Result<Option<Py>, PyErr> = (|| -> Result<Option<Py>, PyErr> {
            v_want = { let t5108 = v_inputs.clone(); f_extinction(py::getitem(&t5108, &Py::str("seq"))?)? };
            #[allow(unreachable_code)]
            Ok(None)
        })();
    match t5109 {
        Ok(Some(r)) => { return Ok(r); }
        Ok(None) => {}
        Err(e) if (!matches!(e, PyErr::Unsupported(_))) => {
            return Ok(Py::Float(f64::INFINITY));
        }
        Err(e) => return Err(e),
    }
    v_worst = Py::Float(0.0);
    for t5110 in py::iter(&(py::method(&v_want.clone(), "items", &[])?))? {
        let t5147 = (t5110).items()?;
        if t5147.len() != 2 { return Err(PyErr::Value("unpack".into())); }
        v_k = t5147[0].clone();
        v_v = t5147[1].clone();
        if ({ let t5111 = v_k.clone(); Py::Bool({ let t5112 = v_outputs.clone(); if !py::contains(&t5112, &t5111)? { true } else { false } }) }).truthy() {
            return Ok(Py::Float(f64::INFINITY));
        }
        v_worst = py::minmax(&[v_worst.clone(), f__rel(py::getitem(&v_outputs.clone(), &v_k.clone())?, v_v.clone())?], true)?;
    }
    return Ok(v_worst.clone());
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn s_5227(mut v_inputs: Py, mut v_outputs: Py) -> R {
    let mut v_k: Py = Py::None;
    let mut v_v: Py = Py::None;
    let mut v_want: Py = Py::None;
    let mut v_worst: Py = Py::None;
    let t5229: Result<Option<Py>, PyErr> = (|| -> Result<Option<Py>, PyErr> {
            v_want = { let t5228 = v_inputs.clone(); f_protein_global(py::getitem(&t5228, &Py::str("seq1"))?, py::getitem(&t5228, &Py::str("seq2"))?, py::getitem(&t5228, &Py::str("open"))?, py::getitem(&t5228, &Py::str("extend"))?)? };
            #[allow(unreachable_code)]
            Ok(None)
        })();
    match t5229 {
        Ok(Some(r)) => { return Ok(r); }
        Ok(None) => {}
        Err(e) if (!matches!(e, PyErr::Unsupported(_))) => {
            return Ok(Py::Float(f64::INFINITY));
        }
        Err(e) => return Err(e),
    }
    v_worst = Py::Float(0.0);
    for t5230 in py::iter(&(py::method(&v_want.clone(), "items", &[])?))? {
        let t5267 = (t5230).items()?;
        if t5267.len() != 2 { return Err(PyErr::Value("unpack".into())); }
        v_k = t5267[0].clone();
        v_v = t5267[1].clone();
        if ({ let t5231 = v_k.clone(); Py::Bool({ let t5232 = v_outputs.clone(); if !py::contains(&t5232, &t5231)? { true } else { false } }) }).truthy() {
            return Ok(Py::Float(f64::INFINITY));
        }
        v_worst = py::minmax(&[v_worst.clone(), f__rel(py::getitem(&v_outputs.clone(), &v_k.clone())?, v_v.clone())?], true)?;
    }
    return Ok(v_worst.clone());
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn s_5347(mut v_inputs: Py, mut v_outputs: Py) -> R {
    let mut v_k: Py = Py::None;
    let mut v_v: Py = Py::None;
    let mut v_want: Py = Py::None;
    let mut v_worst: Py = Py::None;
    let t5349: Result<Option<Py>, PyErr> = (|| -> Result<Option<Py>, PyErr> {
            v_want = { let t5348 = v_inputs.clone(); f_protein_local(py::getitem(&t5348, &Py::str("seq1"))?, py::getitem(&t5348, &Py::str("seq2"))?, py::getitem(&t5348, &Py::str("open"))?, py::getitem(&t5348, &Py::str("extend"))?)? };
            #[allow(unreachable_code)]
            Ok(None)
        })();
    match t5349 {
        Ok(Some(r)) => { return Ok(r); }
        Ok(None) => {}
        Err(e) if (!matches!(e, PyErr::Unsupported(_))) => {
            return Ok(Py::Float(f64::INFINITY));
        }
        Err(e) => return Err(e),
    }
    v_worst = Py::Float(0.0);
    for t5350 in py::iter(&(py::method(&v_want.clone(), "items", &[])?))? {
        let t5387 = (t5350).items()?;
        if t5387.len() != 2 { return Err(PyErr::Value("unpack".into())); }
        v_k = t5387[0].clone();
        v_v = t5387[1].clone();
        if ({ let t5351 = v_k.clone(); Py::Bool({ let t5352 = v_outputs.clone(); if !py::contains(&t5352, &t5351)? { true } else { false } }) }).truthy() {
            return Ok(Py::Float(f64::INFINITY));
        }
        v_worst = py::minmax(&[v_worst.clone(), f__rel(py::getitem(&v_outputs.clone(), &v_k.clone())?, v_v.clone())?], true)?;
    }
    return Ok(v_worst.clone());
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn s_535(mut v_inputs: Py, mut v_outputs: Py) -> R {
    let mut v_k: Py = Py::None;
    let mut v_v: Py = Py::None;
    let mut v_want: Py = Py::None;
    let mut v_worst: Py = Py::None;
    let t537: Result<Option<Py>, PyErr> = (|| -> Result<Option<Py>, PyErr> {
            v_want = { let t536 = v_inputs.clone(); f_coalescent_times(py::getitem(&t536, &Py::str("n"))?, py::getitem(&t536, &Py::str("Ne"))?)? };
            #[allow(unreachable_code)]
            Ok(None)
        })();
    match t537 {
        Ok(Some(r)) => { return Ok(r); }
        Ok(None) => {}
        Err(e) if (!matches!(e, PyErr::Unsupported(_))) => {
            return Ok(Py::Float(f64::INFINITY));
        }
        Err(e) => return Err(e),
    }
    v_worst = Py::Float(0.0);
    for t538 in py::iter(&(py::method(&v_want.clone(), "items", &[])?))? {
        let t575 = (t538).items()?;
        if t575.len() != 2 { return Err(PyErr::Value("unpack".into())); }
        v_k = t575[0].clone();
        v_v = t575[1].clone();
        if ({ let t539 = v_k.clone(); Py::Bool({ let t540 = v_outputs.clone(); if !py::contains(&t540, &t539)? { true } else { false } }) }).truthy() {
            return Ok(Py::Float(f64::INFINITY));
        }
        v_worst = py::minmax(&[v_worst.clone(), f__rel(py::getitem(&v_outputs.clone(), &v_k.clone())?, v_v.clone())?], true)?;
    }
    return Ok(v_worst.clone());
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn s_5426(mut v_inputs: Py, mut v_outputs: Py) -> R {
    let mut v_k: Py = Py::None;
    let mut v_v: Py = Py::None;
    let mut v_want: Py = Py::None;
    let mut v_worst: Py = Py::None;
    let t5428: Result<Option<Py>, PyErr> = (|| -> Result<Option<Py>, PyErr> {
            v_want = { let t5427 = v_inputs.clone(); f_protein_mass(py::getitem(&t5427, &Py::str("seq"))?)? };
            #[allow(unreachable_code)]
            Ok(None)
        })();
    match t5428 {
        Ok(Some(r)) => { return Ok(r); }
        Ok(None) => {}
        Err(e) if (!matches!(e, PyErr::Unsupported(_))) => {
            return Ok(Py::Float(f64::INFINITY));
        }
        Err(e) => return Err(e),
    }
    v_worst = Py::Float(0.0);
    for t5429 in py::iter(&(py::method(&v_want.clone(), "items", &[])?))? {
        let t5466 = (t5429).items()?;
        if t5466.len() != 2 { return Err(PyErr::Value("unpack".into())); }
        v_k = t5466[0].clone();
        v_v = t5466[1].clone();
        if ({ let t5430 = v_k.clone(); Py::Bool({ let t5431 = v_outputs.clone(); if !py::contains(&t5431, &t5430)? { true } else { false } }) }).truthy() {
            return Ok(Py::Float(f64::INFINITY));
        }
        v_worst = py::minmax(&[v_worst.clone(), f__rel(py::getitem(&v_outputs.clone(), &v_k.clone())?, v_v.clone())?], true)?;
    }
    return Ok(v_worst.clone());
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn s_5521(mut v_inputs: Py, mut v_outputs: Py) -> R {
    let mut v_k: Py = Py::None;
    let mut v_v: Py = Py::None;
    let mut v_want: Py = Py::None;
    let mut v_worst: Py = Py::None;
    let t5523: Result<Option<Py>, PyErr> = (|| -> Result<Option<Py>, PyErr> {
            v_want = { let t5522 = v_inputs.clone(); f_net_charge(py::getitem(&t5522, &Py::str("seq"))?, py::getitem(&t5522, &Py::str("pH"))?)? };
            #[allow(unreachable_code)]
            Ok(None)
        })();
    match t5523 {
        Ok(Some(r)) => { return Ok(r); }
        Ok(None) => {}
        Err(e) if (!matches!(e, PyErr::Unsupported(_))) => {
            return Ok(Py::Float(f64::INFINITY));
        }
        Err(e) => return Err(e),
    }
    v_worst = Py::Float(0.0);
    for t5524 in py::iter(&(py::method(&v_want.clone(), "items", &[])?))? {
        let t5561 = (t5524).items()?;
        if t5561.len() != 2 { return Err(PyErr::Value("unpack".into())); }
        v_k = t5561[0].clone();
        v_v = t5561[1].clone();
        if ({ let t5525 = v_k.clone(); Py::Bool({ let t5526 = v_outputs.clone(); if !py::contains(&t5526, &t5525)? { true } else { false } }) }).truthy() {
            return Ok(Py::Float(f64::INFINITY));
        }
        v_worst = py::minmax(&[v_worst.clone(), f__rel(py::getitem(&v_outputs.clone(), &v_k.clone())?, v_v.clone())?], true)?;
    }
    return Ok(v_worst.clone());
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn s_5588(mut v_inputs: Py, mut v_outputs: Py) -> R {
    let mut v_k: Py = Py::None;
    let mut v_v: Py = Py::None;
    let mut v_want: Py = Py::None;
    let mut v_worst: Py = Py::None;
    let t5590: Result<Option<Py>, PyErr> = (|| -> Result<Option<Py>, PyErr> {
            v_want = { let t5589 = v_inputs.clone(); f_rarefaction(py::getitem(&t5589, &Py::str("counts"))?, py::getitem(&t5589, &Py::str("n"))?, py::method(&t5589, "get", &[Py::str("with_replacement"), Py::Bool(false)])?)? };
            #[allow(unreachable_code)]
            Ok(None)
        })();
    match t5590 {
        Ok(Some(r)) => { return Ok(r); }
        Ok(None) => {}
        Err(e) if (!matches!(e, PyErr::Unsupported(_))) => {
            return Ok(Py::Float(f64::INFINITY));
        }
        Err(e) => return Err(e),
    }
    v_worst = Py::Float(0.0);
    for t5591 in py::iter(&(py::method(&v_want.clone(), "items", &[])?))? {
        let t5628 = (t5591).items()?;
        if t5628.len() != 2 { return Err(PyErr::Value("unpack".into())); }
        v_k = t5628[0].clone();
        v_v = t5628[1].clone();
        if ({ let t5592 = v_k.clone(); Py::Bool({ let t5593 = v_outputs.clone(); if !py::contains(&t5593, &t5592)? { true } else { false } }) }).truthy() {
            return Ok(Py::Float(f64::INFINITY));
        }
        v_worst = py::minmax(&[v_worst.clone(), f__rel(py::getitem(&v_outputs.clone(), &v_k.clone())?, v_v.clone())?], true)?;
    }
    return Ok(v_worst.clone());
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn s_5660(mut v_inputs: Py, mut v_outputs: Py) -> R {
    let mut v_k: Py = Py::None;
    let mut v_v: Py = Py::None;
    let mut v_want: Py = Py::None;
    let mut v_worst: Py = Py::None;
    let t5662: Result<Option<Py>, PyErr> = (|| -> Result<Option<Py>, PyErr> {
            v_want = { let t5661 = v_inputs.clone(); f_recessive_lethal(py::getitem(&t5661, &Py::str("q0"))?, py::getitem(&t5661, &Py::str("t"))?)? };
            #[allow(unreachable_code)]
            Ok(None)
        })();
    match t5662 {
        Ok(Some(r)) => { return Ok(r); }
        Ok(None) => {}
        Err(e) if (!matches!(e, PyErr::Unsupported(_))) => {
            return Ok(Py::Float(f64::INFINITY));
        }
        Err(e) => return Err(e),
    }
    v_worst = Py::Float(0.0);
    for t5663 in py::iter(&(py::method(&v_want.clone(), "items", &[])?))? {
        let t5700 = (t5663).items()?;
        if t5700.len() != 2 { return Err(PyErr::Value("unpack".into())); }
        v_k = t5700[0].clone();
        v_v = t5700[1].clone();
        if ({ let t5664 = v_k.clone(); Py::Bool({ let t5665 = v_outputs.clone(); if !py::contains(&t5665, &t5664)? { true } else { false } }) }).truthy() {
            return Ok(Py::Float(f64::INFINITY));
        }
        v_worst = py::minmax(&[v_worst.clone(), f__rel(py::getitem(&v_outputs.clone(), &v_k.clone())?, v_v.clone())?], true)?;
    }
    return Ok(v_worst.clone());
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn s_5724(mut v_inputs: Py, mut v_outputs: Py) -> R {
    let mut v_k: Py = Py::None;
    let mut v_v: Py = Py::None;
    let mut v_want: Py = Py::None;
    let mut v_worst: Py = Py::None;
    let t5726: Result<Option<Py>, PyErr> = (|| -> Result<Option<Py>, PyErr> {
            v_want = { let t5725 = v_inputs.clone(); f_recombination_fraction(py::getitem(&t5725, &Py::str("d_cM"))?, py::getitem(&t5725, &Py::str("map"))?)? };
            #[allow(unreachable_code)]
            Ok(None)
        })();
    match t5726 {
        Ok(Some(r)) => { return Ok(r); }
        Ok(None) => {}
        Err(e) if (!matches!(e, PyErr::Unsupported(_))) => {
            return Ok(Py::Float(f64::INFINITY));
        }
        Err(e) => return Err(e),
    }
    v_worst = Py::Float(0.0);
    for t5727 in py::iter(&(py::method(&v_want.clone(), "items", &[])?))? {
        let t5764 = (t5727).items()?;
        if t5764.len() != 2 { return Err(PyErr::Value("unpack".into())); }
        v_k = t5764[0].clone();
        v_v = t5764[1].clone();
        if ({ let t5728 = v_k.clone(); Py::Bool({ let t5729 = v_outputs.clone(); if !py::contains(&t5729, &t5728)? { true } else { false } }) }).truthy() {
            return Ok(Py::Float(f64::INFINITY));
        }
        v_worst = py::minmax(&[v_worst.clone(), f__rel(py::getitem(&v_outputs.clone(), &v_k.clone())?, v_v.clone())?], true)?;
    }
    return Ok(v_worst.clone());
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn s_5793(mut v_inputs: Py, mut v_outputs: Py) -> R {
    let mut v_k: Py = Py::None;
    let mut v_v: Py = Py::None;
    let mut v_want: Py = Py::None;
    let mut v_worst: Py = Py::None;
    let t5795: Result<Option<Py>, PyErr> = (|| -> Result<Option<Py>, PyErr> {
            v_want = { let t5794 = v_inputs.clone(); f_restriction_sites(py::getitem(&t5794, &Py::str("seq"))?, py::getitem(&t5794, &Py::str("enzyme"))?)? };
            #[allow(unreachable_code)]
            Ok(None)
        })();
    match t5795 {
        Ok(Some(r)) => { return Ok(r); }
        Ok(None) => {}
        Err(e) if (!matches!(e, PyErr::Unsupported(_))) => {
            return Ok(Py::Float(f64::INFINITY));
        }
        Err(e) => return Err(e),
    }
    v_worst = Py::Float(0.0);
    for t5796 in py::iter(&(py::method(&v_want.clone(), "items", &[])?))? {
        let t5833 = (t5796).items()?;
        if t5833.len() != 2 { return Err(PyErr::Value("unpack".into())); }
        v_k = t5833[0].clone();
        v_v = t5833[1].clone();
        if ({ let t5797 = v_k.clone(); Py::Bool({ let t5798 = v_outputs.clone(); if !py::contains(&t5798, &t5797)? { true } else { false } }) }).truthy() {
            return Ok(Py::Float(f64::INFINITY));
        }
        v_worst = py::minmax(&[v_worst.clone(), f__rel(py::getitem(&v_outputs.clone(), &v_k.clone())?, v_v.clone())?], true)?;
    }
    return Ok(v_worst.clone());
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn s_5857(mut v_inputs: Py, mut v_outputs: Py) -> R {
    let mut v_k: Py = Py::None;
    let mut v_v: Py = Py::None;
    let mut v_want: Py = Py::None;
    let mut v_worst: Py = Py::None;
    let t5859: Result<Option<Py>, PyErr> = (|| -> Result<Option<Py>, PyErr> {
            v_want = { let t5858 = v_inputs.clone(); f_reverse_complement(py::getitem(&t5858, &Py::str("seq"))?)? };
            #[allow(unreachable_code)]
            Ok(None)
        })();
    match t5859 {
        Ok(Some(r)) => { return Ok(r); }
        Ok(None) => {}
        Err(e) if (!matches!(e, PyErr::Unsupported(_))) => {
            return Ok(Py::Float(f64::INFINITY));
        }
        Err(e) => return Err(e),
    }
    v_worst = Py::Float(0.0);
    for t5860 in py::iter(&(py::method(&v_want.clone(), "items", &[])?))? {
        let t5897 = (t5860).items()?;
        if t5897.len() != 2 { return Err(PyErr::Value("unpack".into())); }
        v_k = t5897[0].clone();
        v_v = t5897[1].clone();
        if ({ let t5861 = v_k.clone(); Py::Bool({ let t5862 = v_outputs.clone(); if !py::contains(&t5862, &t5861)? { true } else { false } }) }).truthy() {
            return Ok(Py::Float(f64::INFINITY));
        }
        v_worst = py::minmax(&[v_worst.clone(), f__rel(py::getitem(&v_outputs.clone(), &v_k.clone())?, v_v.clone())?], true)?;
    }
    return Ok(v_worst.clone());
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn s_5933(mut v_inputs: Py, mut v_outputs: Py) -> R {
    let mut v_k: Py = Py::None;
    let mut v_v: Py = Py::None;
    let mut v_want: Py = Py::None;
    let mut v_worst: Py = Py::None;
    let t5935: Result<Option<Py>, PyErr> = (|| -> Result<Option<Py>, PyErr> {
            v_want = { let t5934 = v_inputs.clone(); f_richness_indices(py::getitem(&t5934, &Py::str("counts"))?, py::method(&t5934, "get", &[Py::str("log10"), Py::Bool(false)])?)? };
            #[allow(unreachable_code)]
            Ok(None)
        })();
    match t5935 {
        Ok(Some(r)) => { return Ok(r); }
        Ok(None) => {}
        Err(e) if (!matches!(e, PyErr::Unsupported(_))) => {
            return Ok(Py::Float(f64::INFINITY));
        }
        Err(e) => return Err(e),
    }
    v_worst = Py::Float(0.0);
    for t5936 in py::iter(&(py::method(&v_want.clone(), "items", &[])?))? {
        let t5973 = (t5936).items()?;
        if t5973.len() != 2 { return Err(PyErr::Value("unpack".into())); }
        v_k = t5973[0].clone();
        v_v = t5973[1].clone();
        if ({ let t5937 = v_k.clone(); Py::Bool({ let t5938 = v_outputs.clone(); if !py::contains(&t5938, &t5937)? { true } else { false } }) }).truthy() {
            return Ok(Py::Float(f64::INFINITY));
        }
        v_worst = py::minmax(&[v_worst.clone(), f__rel(py::getitem(&v_outputs.clone(), &v_k.clone())?, v_v.clone())?], true)?;
    }
    return Ok(v_worst.clone());
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn s_6033(mut v_inputs: Py, mut v_outputs: Py) -> R {
    let mut v_k: Py = Py::None;
    let mut v_v: Py = Py::None;
    let mut v_want: Py = Py::None;
    let mut v_worst: Py = Py::None;
    let t6035: Result<Option<Py>, PyErr> = (|| -> Result<Option<Py>, PyErr> {
            v_want = { let t6034 = v_inputs.clone(); f_robinson_foulds(py::getitem(&t6034, &Py::str("tree1"))?, py::getitem(&t6034, &Py::str("tree2"))?, py::method(&t6034, "get", &[Py::str("one_way"), Py::Bool(false)])?)? };
            #[allow(unreachable_code)]
            Ok(None)
        })();
    match t6035 {
        Ok(Some(r)) => { return Ok(r); }
        Ok(None) => {}
        Err(e) if (!matches!(e, PyErr::Unsupported(_))) => {
            return Ok(Py::Float(f64::INFINITY));
        }
        Err(e) => return Err(e),
    }
    v_worst = Py::Float(0.0);
    for t6036 in py::iter(&(py::method(&v_want.clone(), "items", &[])?))? {
        let t6073 = (t6036).items()?;
        if t6073.len() != 2 { return Err(PyErr::Value("unpack".into())); }
        v_k = t6073[0].clone();
        v_v = t6073[1].clone();
        if ({ let t6037 = v_k.clone(); Py::Bool({ let t6038 = v_outputs.clone(); if !py::contains(&t6038, &t6037)? { true } else { false } }) }).truthy() {
            return Ok(Py::Float(f64::INFINITY));
        }
        v_worst = py::minmax(&[v_worst.clone(), f__rel(py::getitem(&v_outputs.clone(), &v_k.clone())?, v_v.clone())?], true)?;
    }
    return Ok(v_worst.clone());
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn s_6108(mut v_inputs: Py, mut v_outputs: Py) -> R {
    let mut v_k: Py = Py::None;
    let mut v_v: Py = Py::None;
    let mut v_want: Py = Py::None;
    let mut v_worst: Py = Py::None;
    let t6110: Result<Option<Py>, PyErr> = (|| -> Result<Option<Py>, PyErr> {
            v_want = { let t6109 = v_inputs.clone(); f_schumacher(py::getitem(&t6109, &Py::str("C"))?, py::getitem(&t6109, &Py::str("R"))?, py::getitem(&t6109, &Py::str("M"))?, py::method(&t6109, "get", &[Py::str("schnabel"), Py::Bool(false)])?)? };
            #[allow(unreachable_code)]
            Ok(None)
        })();
    match t6110 {
        Ok(Some(r)) => { return Ok(r); }
        Ok(None) => {}
        Err(e) if (!matches!(e, PyErr::Unsupported(_))) => {
            return Ok(Py::Float(f64::INFINITY));
        }
        Err(e) => return Err(e),
    }
    v_worst = Py::Float(0.0);
    for t6111 in py::iter(&(py::method(&v_want.clone(), "items", &[])?))? {
        let t6148 = (t6111).items()?;
        if t6148.len() != 2 { return Err(PyErr::Value("unpack".into())); }
        v_k = t6148[0].clone();
        v_v = t6148[1].clone();
        if ({ let t6112 = v_k.clone(); Py::Bool({ let t6113 = v_outputs.clone(); if !py::contains(&t6113, &t6112)? { true } else { false } }) }).truthy() {
            return Ok(Py::Float(f64::INFINITY));
        }
        v_worst = py::minmax(&[v_worst.clone(), f__rel(py::getitem(&v_outputs.clone(), &v_k.clone())?, v_v.clone())?], true)?;
    }
    return Ok(v_worst.clone());
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn s_6175(mut v_inputs: Py, mut v_outputs: Py) -> R {
    let mut v_k: Py = Py::None;
    let mut v_v: Py = Py::None;
    let mut v_want: Py = Py::None;
    let mut v_worst: Py = Py::None;
    let t6177: Result<Option<Py>, PyErr> = (|| -> Result<Option<Py>, PyErr> {
            v_want = { let t6176 = v_inputs.clone(); f_ss_fraction(py::getitem(&t6176, &Py::str("seq"))?, py::method(&t6176, "get", &[Py::str("order"), Py::tuple(vec![Py::str("helix"), Py::str("turn"), Py::str("sheet")])])?)? };
            #[allow(unreachable_code)]
            Ok(None)
        })();
    match t6177 {
        Ok(Some(r)) => { return Ok(r); }
        Ok(None) => {}
        Err(e) if (!matches!(e, PyErr::Unsupported(_))) => {
            return Ok(Py::Float(f64::INFINITY));
        }
        Err(e) => return Err(e),
    }
    v_worst = Py::Float(0.0);
    for t6178 in py::iter(&(py::method(&v_want.clone(), "items", &[])?))? {
        let t6215 = (t6178).items()?;
        if t6215.len() != 2 { return Err(PyErr::Value("unpack".into())); }
        v_k = t6215[0].clone();
        v_v = t6215[1].clone();
        if ({ let t6179 = v_k.clone(); Py::Bool({ let t6180 = v_outputs.clone(); if !py::contains(&t6180, &t6179)? { true } else { false } }) }).truthy() {
            return Ok(Py::Float(f64::INFINITY));
        }
        v_worst = py::minmax(&[v_worst.clone(), f__rel(py::getitem(&v_outputs.clone(), &v_k.clone())?, v_v.clone())?], true)?;
    }
    return Ok(v_worst.clone());
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn s_6294(mut v_inputs: Py, mut v_outputs: Py) -> R {
    let mut v_k: Py = Py::None;
    let mut v_v: Py = Py::None;
    let mut v_want: Py = Py::None;
    let mut v_worst: Py = Py::None;
    let t6296: Result<Option<Py>, PyErr> = (|| -> Result<Option<Py>, PyErr> {
            v_want = { let t6295 = v_inputs.clone(); f_semiglobal(py::getitem(&t6295, &Py::str("seq1"))?, py::getitem(&t6295, &Py::str("seq2"))?, py::getitem(&t6295, &Py::str("match"))?, py::getitem(&t6295, &Py::str("mismatch"))?, py::getitem(&t6295, &Py::str("gap"))?)? };
            #[allow(unreachable_code)]
            Ok(None)
        })();
    match t6296 {
        Ok(Some(r)) => { return Ok(r); }
        Ok(None) => {}
        Err(e) if (!matches!(e, PyErr::Unsupported(_))) => {
            return Ok(Py::Float(f64::INFINITY));
        }
        Err(e) => return Err(e),
    }
    v_worst = Py::Float(0.0);
    for t6297 in py::iter(&(py::method(&v_want.clone(), "items", &[])?))? {
        let t6334 = (t6297).items()?;
        if t6334.len() != 2 { return Err(PyErr::Value("unpack".into())); }
        v_k = t6334[0].clone();
        v_v = t6334[1].clone();
        if ({ let t6298 = v_k.clone(); Py::Bool({ let t6299 = v_outputs.clone(); if !py::contains(&t6299, &t6298)? { true } else { false } }) }).truthy() {
            return Ok(Py::Float(f64::INFINITY));
        }
        v_worst = py::minmax(&[v_worst.clone(), f__rel(py::getitem(&v_outputs.clone(), &v_k.clone())?, v_v.clone())?], true)?;
    }
    return Ok(v_worst.clone());
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn s_6372(mut v_inputs: Py, mut v_outputs: Py) -> R {
    let mut v_k: Py = Py::None;
    let mut v_v: Py = Py::None;
    let mut v_want: Py = Py::None;
    let mut v_worst: Py = Py::None;
    let t6374: Result<Option<Py>, PyErr> = (|| -> Result<Option<Py>, PyErr> {
            v_want = { let t6373 = v_inputs.clone(); f_shannon(py::getitem(&t6373, &Py::str("counts"))?, py::getitem(&t6373, &Py::str("base"))?, py::method(&t6373, "get", &[Py::str("wrong"), Py::Bool(false)])?)? };
            #[allow(unreachable_code)]
            Ok(None)
        })();
    match t6374 {
        Ok(Some(r)) => { return Ok(r); }
        Ok(None) => {}
        Err(e) if (!matches!(e, PyErr::Unsupported(_))) => {
            return Ok(Py::Float(f64::INFINITY));
        }
        Err(e) => return Err(e),
    }
    v_worst = Py::Float(0.0);
    for t6375 in py::iter(&(py::method(&v_want.clone(), "items", &[])?))? {
        let t6412 = (t6375).items()?;
        if t6412.len() != 2 { return Err(PyErr::Value("unpack".into())); }
        v_k = t6412[0].clone();
        v_v = t6412[1].clone();
        if ({ let t6376 = v_k.clone(); Py::Bool({ let t6377 = v_outputs.clone(); if !py::contains(&t6377, &t6376)? { true } else { false } }) }).truthy() {
            return Ok(Py::Float(f64::INFINITY));
        }
        v_worst = py::minmax(&[v_worst.clone(), f__rel(py::getitem(&v_outputs.clone(), &v_k.clone())?, v_v.clone())?], true)?;
    }
    return Ok(v_worst.clone());
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn s_6441(mut v_inputs: Py, mut v_outputs: Py) -> R {
    let mut v_k: Py = Py::None;
    let mut v_v: Py = Py::None;
    let mut v_want: Py = Py::None;
    let mut v_worst: Py = Py::None;
    let t6443: Result<Option<Py>, PyErr> = (|| -> Result<Option<Py>, PyErr> {
            v_want = { let t6442 = v_inputs.clone(); f_simpson(py::getitem(&t6442, &Py::str("counts"))?, py::getitem(&t6442, &Py::str("estimator"))?, py::method(&t6442, "get", &[Py::str("swap"), Py::Bool(false)])?)? };
            #[allow(unreachable_code)]
            Ok(None)
        })();
    match t6443 {
        Ok(Some(r)) => { return Ok(r); }
        Ok(None) => {}
        Err(e) if (!matches!(e, PyErr::Unsupported(_))) => {
            return Ok(Py::Float(f64::INFINITY));
        }
        Err(e) => return Err(e),
    }
    v_worst = Py::Float(0.0);
    for t6444 in py::iter(&(py::method(&v_want.clone(), "items", &[])?))? {
        let t6481 = (t6444).items()?;
        if t6481.len() != 2 { return Err(PyErr::Value("unpack".into())); }
        v_k = t6481[0].clone();
        v_v = t6481[1].clone();
        if ({ let t6445 = v_k.clone(); Py::Bool({ let t6446 = v_outputs.clone(); if !py::contains(&t6446, &t6445)? { true } else { false } }) }).truthy() {
            return Ok(Py::Float(f64::INFINITY));
        }
        v_worst = py::minmax(&[v_worst.clone(), f__rel(py::getitem(&v_outputs.clone(), &v_k.clone())?, v_v.clone())?], true)?;
    }
    return Ok(v_worst.clone());
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn s_6520(mut v_inputs: Py, mut v_outputs: Py) -> R {
    let mut v_k: Py = Py::None;
    let mut v_v: Py = Py::None;
    let mut v_want: Py = Py::None;
    let mut v_worst: Py = Py::None;
    let t6522: Result<Option<Py>, PyErr> = (|| -> Result<Option<Py>, PyErr> {
            v_want = { let t6521 = v_inputs.clone(); f_species_area(py::getitem(&t6521, &Py::str("A"))?, py::getitem(&t6521, &Py::str("S"))?, py::method(&t6521, "get", &[Py::str("no_exp"), Py::Bool(false)])?)? };
            #[allow(unreachable_code)]
            Ok(None)
        })();
    match t6522 {
        Ok(Some(r)) => { return Ok(r); }
        Ok(None) => {}
        Err(e) if (!matches!(e, PyErr::Unsupported(_))) => {
            return Ok(Py::Float(f64::INFINITY));
        }
        Err(e) => return Err(e),
    }
    v_worst = Py::Float(0.0);
    for t6523 in py::iter(&(py::method(&v_want.clone(), "items", &[])?))? {
        let t6560 = (t6523).items()?;
        if t6560.len() != 2 { return Err(PyErr::Value("unpack".into())); }
        v_k = t6560[0].clone();
        v_v = t6560[1].clone();
        if ({ let t6524 = v_k.clone(); Py::Bool({ let t6525 = v_outputs.clone(); if !py::contains(&t6525, &t6524)? { true } else { false } }) }).truthy() {
            return Ok(Py::Float(f64::INFINITY));
        }
        v_worst = py::minmax(&[v_worst.clone(), f__rel(py::getitem(&v_outputs.clone(), &v_k.clone())?, v_v.clone())?], true)?;
    }
    return Ok(v_worst.clone());
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn s_6598(mut v_inputs: Py, mut v_outputs: Py) -> R {
    let mut v_k: Py = Py::None;
    let mut v_v: Py = Py::None;
    let mut v_want: Py = Py::None;
    let mut v_worst: Py = Py::None;
    let t6600: Result<Option<Py>, PyErr> = (|| -> Result<Option<Py>, PyErr> {
            v_want = { let t6599 = v_inputs.clone(); f_tajima_d(py::getitem(&t6599, &Py::str("seqs"))?, py::method(&t6599, "get", &[Py::str("sqrt_denominator"), Py::Bool(true)])?)? };
            #[allow(unreachable_code)]
            Ok(None)
        })();
    match t6600 {
        Ok(Some(r)) => { return Ok(r); }
        Ok(None) => {}
        Err(e) if (!matches!(e, PyErr::Unsupported(_))) => {
            return Ok(Py::Float(f64::INFINITY));
        }
        Err(e) => return Err(e),
    }
    v_worst = Py::Float(0.0);
    for t6601 in py::iter(&(py::method(&v_want.clone(), "items", &[])?))? {
        let t6638 = (t6601).items()?;
        if t6638.len() != 2 { return Err(PyErr::Value("unpack".into())); }
        v_k = t6638[0].clone();
        v_v = t6638[1].clone();
        if ({ let t6602 = v_k.clone(); Py::Bool({ let t6603 = v_outputs.clone(); if !py::contains(&t6603, &t6602)? { true } else { false } }) }).truthy() {
            return Ok(Py::Float(f64::INFINITY));
        }
        v_worst = py::minmax(&[v_worst.clone(), f__rel(py::getitem(&v_outputs.clone(), &v_k.clone())?, v_v.clone())?], true)?;
    }
    return Ok(v_worst.clone());
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn s_6705(mut v_inputs: Py, mut v_outputs: Py) -> R {
    let mut v_k: Py = Py::None;
    let mut v_v: Py = Py::None;
    let mut v_want: Py = Py::None;
    let mut v_worst: Py = Py::None;
    let t6707: Result<Option<Py>, PyErr> = (|| -> Result<Option<Py>, PyErr> {
            v_want = { let t6706 = v_inputs.clone(); f_tm_nn(py::getitem(&t6706, &Py::str("seq"))?, py::getitem(&t6706, &Py::str("Na"))?, py::getitem(&t6706, &Py::str("dnac1"))?, py::getitem(&t6706, &Py::str("dnac2"))?)? };
            #[allow(unreachable_code)]
            Ok(None)
        })();
    match t6707 {
        Ok(Some(r)) => { return Ok(r); }
        Ok(None) => {}
        Err(e) if (!matches!(e, PyErr::Unsupported(_))) => {
            return Ok(Py::Float(f64::INFINITY));
        }
        Err(e) => return Err(e),
    }
    v_worst = Py::Float(0.0);
    for t6708 in py::iter(&(py::method(&v_want.clone(), "items", &[])?))? {
        let t6745 = (t6708).items()?;
        if t6745.len() != 2 { return Err(PyErr::Value("unpack".into())); }
        v_k = t6745[0].clone();
        v_v = t6745[1].clone();
        if ({ let t6709 = v_k.clone(); Py::Bool({ let t6710 = v_outputs.clone(); if !py::contains(&t6710, &t6709)? { true } else { false } }) }).truthy() {
            return Ok(Py::Float(f64::INFINITY));
        }
        v_worst = py::minmax(&[v_worst.clone(), f__rel(py::getitem(&v_outputs.clone(), &v_k.clone())?, v_v.clone())?], true)?;
    }
    return Ok(v_worst.clone());
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn s_6769(mut v_inputs: Py, mut v_outputs: Py) -> R {
    let mut v_k: Py = Py::None;
    let mut v_v: Py = Py::None;
    let mut v_want: Py = Py::None;
    let mut v_worst: Py = Py::None;
    let t6771: Result<Option<Py>, PyErr> = (|| -> Result<Option<Py>, PyErr> {
            v_want = { let t6770 = v_inputs.clone(); f_tm_wallace(py::getitem(&t6770, &Py::str("seq"))?)? };
            #[allow(unreachable_code)]
            Ok(None)
        })();
    match t6771 {
        Ok(Some(r)) => { return Ok(r); }
        Ok(None) => {}
        Err(e) if (!matches!(e, PyErr::Unsupported(_))) => {
            return Ok(Py::Float(f64::INFINITY));
        }
        Err(e) => return Err(e),
    }
    v_worst = Py::Float(0.0);
    for t6772 in py::iter(&(py::method(&v_want.clone(), "items", &[])?))? {
        let t6809 = (t6772).items()?;
        if t6809.len() != 2 { return Err(PyErr::Value("unpack".into())); }
        v_k = t6809[0].clone();
        v_v = t6809[1].clone();
        if ({ let t6773 = v_k.clone(); Py::Bool({ let t6774 = v_outputs.clone(); if !py::contains(&t6774, &t6773)? { true } else { false } }) }).truthy() {
            return Ok(Py::Float(f64::INFINITY));
        }
        v_worst = py::minmax(&[v_worst.clone(), f__rel(py::getitem(&v_outputs.clone(), &v_k.clone())?, v_v.clone())?], true)?;
    }
    return Ok(v_worst.clone());
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn s_6837(mut v_inputs: Py, mut v_outputs: Py) -> R {
    let mut v_k: Py = Py::None;
    let mut v_v: Py = Py::None;
    let mut v_want: Py = Py::None;
    let mut v_worst: Py = Py::None;
    let t6839: Result<Option<Py>, PyErr> = (|| -> Result<Option<Py>, PyErr> {
            v_want = { let t6838 = v_inputs.clone(); f_transcription(py::getitem(&t6838, &Py::str("seq"))?, py::getitem(&t6838, &Py::str("strand"))?)? };
            #[allow(unreachable_code)]
            Ok(None)
        })();
    match t6839 {
        Ok(Some(r)) => { return Ok(r); }
        Ok(None) => {}
        Err(e) if (!matches!(e, PyErr::Unsupported(_))) => {
            return Ok(Py::Float(f64::INFINITY));
        }
        Err(e) => return Err(e),
    }
    v_worst = Py::Float(0.0);
    for t6840 in py::iter(&(py::method(&v_want.clone(), "items", &[])?))? {
        let t6877 = (t6840).items()?;
        if t6877.len() != 2 { return Err(PyErr::Value("unpack".into())); }
        v_k = t6877[0].clone();
        v_v = t6877[1].clone();
        if ({ let t6841 = v_k.clone(); Py::Bool({ let t6842 = v_outputs.clone(); if !py::contains(&t6842, &t6841)? { true } else { false } }) }).truthy() {
            return Ok(Py::Float(f64::INFINITY));
        }
        v_worst = py::minmax(&[v_worst.clone(), f__rel(py::getitem(&v_outputs.clone(), &v_k.clone())?, v_v.clone())?], true)?;
    }
    return Ok(v_worst.clone());
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn s_6911(mut v_inputs: Py, mut v_outputs: Py) -> R {
    let mut v_k: Py = Py::None;
    let mut v_v: Py = Py::None;
    let mut v_want: Py = Py::None;
    let mut v_worst: Py = Py::None;
    let t6913: Result<Option<Py>, PyErr> = (|| -> Result<Option<Py>, PyErr> {
            v_want = { let t6912 = v_inputs.clone(); f_translation(py::getitem(&t6912, &Py::str("seq"))?, py::getitem(&t6912, &Py::str("stops"))?)? };
            #[allow(unreachable_code)]
            Ok(None)
        })();
    match t6913 {
        Ok(Some(r)) => { return Ok(r); }
        Ok(None) => {}
        Err(e) if (!matches!(e, PyErr::Unsupported(_))) => {
            return Ok(Py::Float(f64::INFINITY));
        }
        Err(e) => return Err(e),
    }
    v_worst = Py::Float(0.0);
    for t6914 in py::iter(&(py::method(&v_want.clone(), "items", &[])?))? {
        let t6951 = (t6914).items()?;
        if t6951.len() != 2 { return Err(PyErr::Value("unpack".into())); }
        v_k = t6951[0].clone();
        v_v = t6951[1].clone();
        if ({ let t6915 = v_k.clone(); Py::Bool({ let t6916 = v_outputs.clone(); if !py::contains(&t6916, &t6915)? { true } else { false } }) }).truthy() {
            return Ok(Py::Float(f64::INFINITY));
        }
        v_worst = py::minmax(&[v_worst.clone(), f__rel(py::getitem(&v_outputs.clone(), &v_k.clone())?, v_v.clone())?], true)?;
    }
    return Ok(v_worst.clone());
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn s_6964(mut v_inputs: Py, mut v_outputs: Py) -> R {
    let mut v_k: Py = Py::None;
    let mut v_v: Py = Py::None;
    let mut v_want: Py = Py::None;
    let mut v_worst: Py = Py::None;
    let t6966: Result<Option<Py>, PyErr> = (|| -> Result<Option<Py>, PyErr> {
            v_want = { let t6965 = v_inputs.clone(); f_tree_count(py::getitem(&t6965, &Py::str("n"))?, py::getitem(&t6965, &Py::str("rooted"))?)? };
            #[allow(unreachable_code)]
            Ok(None)
        })();
    match t6966 {
        Ok(Some(r)) => { return Ok(r); }
        Ok(None) => {}
        Err(e) if (!matches!(e, PyErr::Unsupported(_))) => {
            return Ok(Py::Float(f64::INFINITY));
        }
        Err(e) => return Err(e),
    }
    v_worst = Py::Float(0.0);
    for t6967 in py::iter(&(py::method(&v_want.clone(), "items", &[])?))? {
        let t7004 = (t6967).items()?;
        if t7004.len() != 2 { return Err(PyErr::Value("unpack".into())); }
        v_k = t7004[0].clone();
        v_v = t7004[1].clone();
        if ({ let t6968 = v_k.clone(); Py::Bool({ let t6969 = v_outputs.clone(); if !py::contains(&t6969, &t6968)? { true } else { false } }) }).truthy() {
            return Ok(Py::Float(f64::INFINITY));
        }
        v_worst = py::minmax(&[v_worst.clone(), f__rel(py::getitem(&v_outputs.clone(), &v_k.clone())?, v_v.clone())?], true)?;
    }
    return Ok(v_worst.clone());
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn s_7065(mut v_inputs: Py, mut v_outputs: Py) -> R {
    let mut v_k: Py = Py::None;
    let mut v_v: Py = Py::None;
    let mut v_want: Py = Py::None;
    let mut v_worst: Py = Py::None;
    let t7067: Result<Option<Py>, PyErr> = (|| -> Result<Option<Py>, PyErr> {
            v_want = { let t7066 = v_inputs.clone(); f_upgma(py::getitem(&t7066, &Py::str("D"))?, py::method(&t7066, "get", &[Py::str("linkage"), Py::str("average")])?)? };
            #[allow(unreachable_code)]
            Ok(None)
        })();
    match t7067 {
        Ok(Some(r)) => { return Ok(r); }
        Ok(None) => {}
        Err(e) if (!matches!(e, PyErr::Unsupported(_))) => {
            return Ok(Py::Float(f64::INFINITY));
        }
        Err(e) => return Err(e),
    }
    v_worst = Py::Float(0.0);
    for t7068 in py::iter(&(py::method(&v_want.clone(), "items", &[])?))? {
        let t7105 = (t7068).items()?;
        if t7105.len() != 2 { return Err(PyErr::Value("unpack".into())); }
        v_k = t7105[0].clone();
        v_v = t7105[1].clone();
        if ({ let t7069 = v_k.clone(); Py::Bool({ let t7070 = v_outputs.clone(); if !py::contains(&t7070, &t7069)? { true } else { false } }) }).truthy() {
            return Ok(Py::Float(f64::INFINITY));
        }
        v_worst = py::minmax(&[v_worst.clone(), f__rel(py::getitem(&v_outputs.clone(), &v_k.clone())?, v_v.clone())?], true)?;
    }
    return Ok(v_worst.clone());
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn s_7134(mut v_inputs: Py, mut v_outputs: Py) -> R {
    let mut v_k: Py = Py::None;
    let mut v_v: Py = Py::None;
    let mut v_want: Py = Py::None;
    let mut v_worst: Py = Py::None;
    let t7136: Result<Option<Py>, PyErr> = (|| -> Result<Option<Py>, PyErr> {
            v_want = { let t7135 = v_inputs.clone(); f_watterson(py::getitem(&t7135, &Py::str("seqs"))?, py::method(&t7135, "get", &[Py::str("harmonic_to_n"), Py::Bool(false)])?)? };
            #[allow(unreachable_code)]
            Ok(None)
        })();
    match t7136 {
        Ok(Some(r)) => { return Ok(r); }
        Ok(None) => {}
        Err(e) if (!matches!(e, PyErr::Unsupported(_))) => {
            return Ok(Py::Float(f64::INFINITY));
        }
        Err(e) => return Err(e),
    }
    v_worst = Py::Float(0.0);
    for t7137 in py::iter(&(py::method(&v_want.clone(), "items", &[])?))? {
        let t7174 = (t7137).items()?;
        if t7174.len() != 2 { return Err(PyErr::Value("unpack".into())); }
        v_k = t7174[0].clone();
        v_v = t7174[1].clone();
        if ({ let t7138 = v_k.clone(); Py::Bool({ let t7139 = v_outputs.clone(); if !py::contains(&t7139, &t7138)? { true } else { false } }) }).truthy() {
            return Ok(Py::Float(f64::INFINITY));
        }
        v_worst = py::minmax(&[v_worst.clone(), f__rel(py::getitem(&v_outputs.clone(), &v_k.clone())?, v_v.clone())?], true)?;
    }
    return Ok(v_worst.clone());
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn s_776(mut v_inputs: Py, mut v_outputs: Py) -> R {
    let mut v_k: Py = Py::None;
    let mut v_v: Py = Py::None;
    let mut v_want: Py = Py::None;
    let mut v_worst: Py = Py::None;
    let t778: Result<Option<Py>, PyErr> = (|| -> Result<Option<Py>, PyErr> {
            v_want = { let t777 = v_inputs.clone(); f_competitive_fit(py::getitem(&t777, &Py::str("S"))?, py::getitem(&t777, &Py::str("I"))?, py::getitem(&t777, &Py::str("v"))?, py::method(&t777, "get", &[Py::str("model"), Py::str("competitive")])?)? };
            #[allow(unreachable_code)]
            Ok(None)
        })();
    match t778 {
        Ok(Some(r)) => { return Ok(r); }
        Ok(None) => {}
        Err(e) if (!matches!(e, PyErr::Unsupported(_))) => {
            return Ok(Py::Float(f64::INFINITY));
        }
        Err(e) => return Err(e),
    }
    v_worst = Py::Float(0.0);
    for t779 in py::iter(&(py::method(&v_want.clone(), "items", &[])?))? {
        let t816 = (t779).items()?;
        if t816.len() != 2 { return Err(PyErr::Value("unpack".into())); }
        v_k = t816[0].clone();
        v_v = t816[1].clone();
        if ({ let t780 = v_k.clone(); Py::Bool({ let t781 = v_outputs.clone(); if !py::contains(&t781, &t780)? { true } else { false } }) }).truthy() {
            return Ok(Py::Float(f64::INFINITY));
        }
        v_worst = py::minmax(&[v_worst.clone(), f__rel(py::getitem(&v_outputs.clone(), &v_k.clone())?, v_v.clone())?], true)?;
    }
    return Ok(v_worst.clone());
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn s_847(mut v_inputs: Py, mut v_outputs: Py) -> R {
    let mut v_k: Py = Py::None;
    let mut v_v: Py = Py::None;
    let mut v_want: Py = Py::None;
    let mut v_worst: Py = Py::None;
    let t849: Result<Option<Py>, PyErr> = (|| -> Result<Option<Py>, PyErr> {
            v_want = { let t848 = v_inputs.clone(); f_cross_genotype(py::getitem(&t848, &Py::str("parent1"))?, py::getitem(&t848, &Py::str("parent2"))?, py::getitem(&t848, &Py::str("offspring"))?)? };
            #[allow(unreachable_code)]
            Ok(None)
        })();
    match t849 {
        Ok(Some(r)) => { return Ok(r); }
        Ok(None) => {}
        Err(e) if (!matches!(e, PyErr::Unsupported(_))) => {
            return Ok(Py::Float(f64::INFINITY));
        }
        Err(e) => return Err(e),
    }
    v_worst = Py::Float(0.0);
    for t850 in py::iter(&(py::method(&v_want.clone(), "items", &[])?))? {
        let t887 = (t850).items()?;
        if t887.len() != 2 { return Err(PyErr::Value("unpack".into())); }
        v_k = t887[0].clone();
        v_v = t887[1].clone();
        if ({ let t851 = v_k.clone(); Py::Bool({ let t852 = v_outputs.clone(); if !py::contains(&t852, &t851)? { true } else { false } }) }).truthy() {
            return Ok(Py::Float(f64::INFINITY));
        }
        v_worst = py::minmax(&[v_worst.clone(), f__rel(py::getitem(&v_outputs.clone(), &v_k.clone())?, v_v.clone())?], true)?;
    }
    return Ok(v_worst.clone());
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn s_938(mut v_inputs: Py, mut v_outputs: Py) -> R {
    let mut v_k: Py = Py::None;
    let mut v_v: Py = Py::None;
    let mut v_want: Py = Py::None;
    let mut v_worst: Py = Py::None;
    let t940: Result<Option<Py>, PyErr> = (|| -> Result<Option<Py>, PyErr> {
            v_want = { let t939 = v_inputs.clone(); f_cross_phenotypes(py::getitem(&t939, &Py::str("parent1"))?, py::getitem(&t939, &Py::str("parent2"))?, py::method(&t939, "get", &[Py::str("complete"), Py::Bool(true)])?)? };
            #[allow(unreachable_code)]
            Ok(None)
        })();
    match t940 {
        Ok(Some(r)) => { return Ok(r); }
        Ok(None) => {}
        Err(e) if (!matches!(e, PyErr::Unsupported(_))) => {
            return Ok(Py::Float(f64::INFINITY));
        }
        Err(e) => return Err(e),
    }
    v_worst = Py::Float(0.0);
    for t941 in py::iter(&(py::method(&v_want.clone(), "items", &[])?))? {
        let t978 = (t941).items()?;
        if t978.len() != 2 { return Err(PyErr::Value("unpack".into())); }
        v_k = t978[0].clone();
        v_v = t978[1].clone();
        if ({ let t942 = v_k.clone(); Py::Bool({ let t943 = v_outputs.clone(); if !py::contains(&t943, &t942)? { true } else { false } }) }).truthy() {
            return Ok(Py::Float(f64::INFINITY));
        }
        v_worst = py::minmax(&[v_worst.clone(), f__rel(py::getitem(&v_outputs.clone(), &v_k.clone())?, v_v.clone())?], true)?;
    }
    return Ok(v_worst.clone());
    #[allow(unreachable_code)]
    Ok(Py::None)
}
