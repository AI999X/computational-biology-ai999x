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
    "bio_hardy_weinberg_frequencies",
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
    "bio_testcross_recombination",
    "bio_three_point_testcross",
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
        "bio_hardy_weinberg_frequencies" => {
            let r = f_hardy_weinberg_frequencies(match v.get("p") { Some(x) => py::from_val(x), None => { (|| -> R { Ok(Py::None) })().ok()? } }, match v.get("q") { Some(x) => py::from_val(x), None => { (|| -> R { Ok(Py::None) })().ok()? } }, match v.get("q2") { Some(x) => py::from_val(x), None => { (|| -> R { Ok(Py::None) })().ok()? } }, match v.get("N") { Some(x) => py::from_val(x), None => { (|| -> R { Ok(Py::None) })().ok()? } });
            finish(r, |i: &Py, o: &Py| -> f64 { match s_2496(i.clone(), o.clone()).and_then(|r| r.f()) { Ok(x) => x, Err(e) => { if std::env::var_os("AX_DEBUG").is_some() { eprintln!("check error: {e:?}"); } f64::NAN } } }, v)
        }
        "bio_hardy_weinberg_test" => {
            let r = f_hardy_weinberg_test(match v.get("n_AA") { Some(x) => py::from_val(x), None => { return None } }, match v.get("n_Aa") { Some(x) => py::from_val(x), None => { return None } }, match v.get("n_aa") { Some(x) => py::from_val(x), None => { return None } });
            finish(r, |i: &Py, o: &Py| -> f64 { match s_2564(i.clone(), o.clone()).and_then(|r| r.f()) { Ok(x) => x, Err(e) => { if std::env::var_os("AX_DEBUG").is_some() { eprintln!("check error: {e:?}"); } f64::NAN } } }, v)
        }
        "bio_heterozygote_advantage" => {
            let r = f_het_advantage(match v.get("w_AA") { Some(x) => py::from_val(x), None => { return None } }, match v.get("w_Aa") { Some(x) => py::from_val(x), None => { return None } }, match v.get("w_aa") { Some(x) => py::from_val(x), None => { return None } });
            finish(r, |i: &Py, o: &Py| -> f64 { match s_2630(i.clone(), o.clone()).and_then(|r| r.f()) { Ok(x) => x, Err(e) => { if std::env::var_os("AX_DEBUG").is_some() { eprintln!("check error: {e:?}"); } f64::NAN } } }, v)
        }
        "bio_hill_fit" => {
            let r = f_hill_fit(match v.get("S") { Some(x) => py::from_val(x), None => { return None } }, match v.get("v") { Some(x) => py::from_val(x), None => { return None } }, (|| -> R { Ok(Py::Bool(false)) })().ok()?);
            finish(r, |i: &Py, o: &Py| -> f64 { match s_2846(i.clone(), o.clone()).and_then(|r| r.f()) { Ok(x) => x, Err(e) => { if std::env::var_os("AX_DEBUG").is_some() { eprintln!("check error: {e:?}"); } f64::NAN } } }, v)
        }
        "bio_hill_plot" => {
            let r = f_hill_plot(match v.get("S") { Some(x) => py::from_val(x), None => { return None } }, match v.get("Y") { Some(x) => py::from_val(x), None => { return None } });
            finish(r, |i: &Py, o: &Py| -> f64 { match s_2942(i.clone(), o.clone()).and_then(|r| r.f()) { Ok(x) => x, Err(e) => { if std::env::var_os("AX_DEBUG").is_some() { eprintln!("check error: {e:?}"); } f64::NAN } } }, v)
        }
        "bio_instability_index" => {
            let r = f_instability(match v.get("seq") { Some(x) => py::from_val(x), None => { return None } });
            finish(r, |i: &Py, o: &Py| -> f64 { match s_3023(i.clone(), o.clone()).and_then(|r| r.f()) { Ok(x) => x, Err(e) => { if std::env::var_os("AX_DEBUG").is_some() { eprintln!("check error: {e:?}"); } f64::NAN } } }, v)
        }
        "bio_isoelectric_point" => {
            let r = f_isoelectric_point(match v.get("seq") { Some(x) => py::from_val(x), None => { return None } });
            finish(r, |i: &Py, o: &Py| -> f64 { match s_3138(i.clone(), o.clone()).and_then(|r| r.f()) { Ok(x) => x, Err(e) => { if std::env::var_os("AX_DEBUG").is_some() { eprintln!("check error: {e:?}"); } f64::NAN } } }, v)
        }
        "bio_jukes_cantor_distance" => {
            let r = f_jukes_cantor(match v.get("seq1") { Some(x) => py::from_val(x), None => { return None } }, match v.get("seq2") { Some(x) => py::from_val(x), None => { return None } });
            finish(r, |i: &Py, o: &Py| -> f64 { match s_3226(i.clone(), o.clone()).and_then(|r| r.f()) { Ok(x) => x, Err(e) => { if std::env::var_os("AX_DEBUG").is_some() { eprintln!("check error: {e:?}"); } f64::NAN } } }, v)
        }
        "bio_kimura_distance" => {
            let r = f_kimura(match v.get("seq1") { Some(x) => py::from_val(x), None => { return None } }, match v.get("seq2") { Some(x) => py::from_val(x), None => { return None } });
            finish(r, |i: &Py, o: &Py| -> f64 { match s_3333(i.clone(), o.clone()).and_then(|r| r.f()) { Ok(x) => x, Err(e) => { if std::env::var_os("AX_DEBUG").is_some() { eprintln!("check error: {e:?}"); } f64::NAN } } }, v)
        }
        "bio_kmer_counts" => {
            let r = f_kmer_counts(match v.get("seq") { Some(x) => py::from_val(x), None => { return None } }, match v.get("k") { Some(x) => py::from_val(x), None => { return None } });
            finish(r, |i: &Py, o: &Py| -> f64 { match s_3406(i.clone(), o.clone()).and_then(|r| r.f()) { Ok(x) => x, Err(e) => { if std::env::var_os("AX_DEBUG").is_some() { eprintln!("check error: {e:?}"); } f64::NAN } } }, v)
        }
        "bio_leslie_growth_rate" => {
            let r = f_leslie_growth(match v.get("F") { Some(x) => py::from_val(x), None => { return None } }, match v.get("S") { Some(x) => py::from_val(x), None => { return None } }, (|| -> R { Ok(Py::Bool(false)) })().ok()?);
            finish(r, |i: &Py, o: &Py| -> f64 { match s_3520(i.clone(), o.clone()).and_then(|r| r.f()) { Ok(x) => x, Err(e) => { if std::env::var_os("AX_DEBUG").is_some() { eprintln!("check error: {e:?}"); } f64::NAN } } }, v)
        }
        "bio_leslie_projection" => {
            let r = f_leslie_projection(match v.get("F") { Some(x) => py::from_val(x), None => { return None } }, match v.get("S") { Some(x) => py::from_val(x), None => { return None } }, match v.get("n0") { Some(x) => py::from_val(x), None => { return None } }, match v.get("t") { Some(x) => py::from_val(x), None => { return None } }, (|| -> R { Ok(Py::Bool(false)) })().ok()?);
            finish(r, |i: &Py, o: &Py| -> f64 { match s_3614(i.clone(), o.clone()).and_then(|r| r.f()) { Ok(x) => x, Err(e) => { if std::env::var_os("AX_DEBUG").is_some() { eprintln!("check error: {e:?}"); } f64::NAN } } }, v)
        }
        "bio_lineweaver_burk" => {
            let r = f_lineweaver_burk(match v.get("S") { Some(x) => py::from_val(x), None => { return None } }, match v.get("v") { Some(x) => py::from_val(x), None => { return None } });
            finish(r, |i: &Py, o: &Py| -> f64 { match s_3703(i.clone(), o.clone()).and_then(|r| r.f()) { Ok(x) => x, Err(e) => { if std::env::var_os("AX_DEBUG").is_some() { eprintln!("check error: {e:?}"); } f64::NAN } } }, v)
        }
        "bio_linkage_disequilibrium" => {
            let r = f_linkage_disequilibrium(match v.get("n_AB") { Some(x) => py::from_val(x), None => { return None } }, match v.get("n_Ab") { Some(x) => py::from_val(x), None => { return None } }, match v.get("n_aB") { Some(x) => py::from_val(x), None => { return None } }, match v.get("n_ab") { Some(x) => py::from_val(x), None => { return None } });
            finish(r, |i: &Py, o: &Py| -> f64 { match s_3762(i.clone(), o.clone()).and_then(|r| r.f()) { Ok(x) => x, Err(e) => { if std::env::var_os("AX_DEBUG").is_some() { eprintln!("check error: {e:?}"); } f64::NAN } } }, v)
        }
        "bio_local_alignment_affine" => {
            let r = f_local_affine(match v.get("seq1") { Some(x) => py::from_val(x), None => { return None } }, match v.get("seq2") { Some(x) => py::from_val(x), None => { return None } }, match v.get("match") { Some(x) => py::from_val(x), None => { return None } }, match v.get("mismatch") { Some(x) => py::from_val(x), None => { return None } }, match v.get("open") { Some(x) => py::from_val(x), None => { return None } }, match v.get("extend") { Some(x) => py::from_val(x), None => { return None } });
            finish(r, |i: &Py, o: &Py| -> f64 { match s_3886(i.clone(), o.clone()).and_then(|r| r.f()) { Ok(x) => x, Err(e) => { if std::env::var_os("AX_DEBUG").is_some() { eprintln!("check error: {e:?}"); } f64::NAN } } }, v)
        }
        "bio_local_alignment_score" => {
            let r = f_local_linear(match v.get("seq1") { Some(x) => py::from_val(x), None => { return None } }, match v.get("seq2") { Some(x) => py::from_val(x), None => { return None } }, match v.get("match") { Some(x) => py::from_val(x), None => { return None } }, match v.get("mismatch") { Some(x) => py::from_val(x), None => { return None } }, match v.get("gap") { Some(x) => py::from_val(x), None => { return None } });
            finish(r, |i: &Py, o: &Py| -> f64 { match s_4007(i.clone(), o.clone()).and_then(|r| r.f()) { Ok(x) => x, Err(e) => { if std::env::var_os("AX_DEBUG").is_some() { eprintln!("check error: {e:?}"); } f64::NAN } } }, v)
        }
        "bio_longest_orf" => {
            let r = f_longest_orf(match v.get("seq") { Some(x) => py::from_val(x), None => { return None } });
            finish(r, |i: &Py, o: &Py| -> f64 { match s_4079(i.clone(), o.clone()).and_then(|r| r.f()) { Ok(x) => x, Err(e) => { if std::env::var_os("AX_DEBUG").is_some() { eprintln!("check error: {e:?}"); } f64::NAN } } }, v)
        }
        "bio_map_distance" => {
            let r = f_map_distance(match v.get("r") { Some(x) => py::from_val(x), None => { return None } }, match v.get("map") { Some(x) => py::from_val(x), None => { return None } });
            finish(r, |i: &Py, o: &Py| -> f64 { match s_4144(i.clone(), o.clone()).and_then(|r| r.f()) { Ok(x) => x, Err(e) => { if std::env::var_os("AX_DEBUG").is_some() { eprintln!("check error: {e:?}"); } f64::NAN } } }, v)
        }
        "bio_michaelis_menten_fit" => {
            let r = f_mm_fit(match v.get("S") { Some(x) => py::from_val(x), None => { return None } }, match v.get("v") { Some(x) => py::from_val(x), None => { return None } });
            finish(r, |i: &Py, o: &Py| -> f64 { match s_4261(i.clone(), o.clone()).and_then(|r| r.f()) { Ok(x) => x, Err(e) => { if std::env::var_os("AX_DEBUG").is_some() { eprintln!("check error: {e:?}"); } f64::NAN } } }, v)
        }
        "bio_michaelis_menten_progress" => {
            let r = f_mm_progress(match v.get("S0") { Some(x) => py::from_val(x), None => { return None } }, match v.get("Vmax") { Some(x) => py::from_val(x), None => { return None } }, match v.get("Km") { Some(x) => py::from_val(x), None => { return None } }, match v.get("t") { Some(x) => py::from_val(x), None => { return None } }, (|| -> R { Ok(Py::Bool(false)) })().ok()?);
            finish(r, |i: &Py, o: &Py| -> f64 { match s_4359(i.clone(), o.clone()).and_then(|r| r.f()) { Ok(x) => x, Err(e) => { if std::env::var_os("AX_DEBUG").is_some() { eprintln!("check error: {e:?}"); } f64::NAN } } }, v)
        }
        "bio_michaelis_menten_time" => {
            let r = f_mm_time(match v.get("S0") { Some(x) => py::from_val(x), None => { return None } }, match v.get("Vmax") { Some(x) => py::from_val(x), None => { return None } }, match v.get("Km") { Some(x) => py::from_val(x), None => { return None } }, match v.get("X") { Some(x) => py::from_val(x), None => { return None } }, (|| -> R { Ok(Py::Bool(false)) })().ok()?);
            finish(r, |i: &Py, o: &Py| -> f64 { match s_4426(i.clone(), o.clone()).and_then(|r| r.f()) { Ok(x) => x, Err(e) => { if std::env::var_os("AX_DEBUG").is_some() { eprintln!("check error: {e:?}"); } f64::NAN } } }, v)
        }
        "bio_moran_fixation" => {
            let r = f_moran_fixation(match v.get("N") { Some(x) => py::from_val(x), None => { return None } }, match v.get("i") { Some(x) => py::from_val(x), None => { return None } }, match v.get("r") { Some(x) => py::from_val(x), None => { return None } });
            finish(r, |i: &Py, o: &Py| -> f64 { match s_4501(i.clone(), o.clone()).and_then(|r| r.f()) { Ok(x) => x, Err(e) => { if std::env::var_os("AX_DEBUG").is_some() { eprintln!("check error: {e:?}"); } f64::NAN } } }, v)
        }
        "bio_mutation_selection_balance" => {
            let r = f_mutation_selection(match v.get("mu") { Some(x) => py::from_val(x), None => { return None } }, match v.get("s") { Some(x) => py::from_val(x), None => { return None } }, match v.get("h") { Some(x) => py::from_val(x), None => { return None } });
            finish(r, |i: &Py, o: &Py| -> f64 { match s_4598(i.clone(), o.clone()).and_then(|r| r.f()) { Ok(x) => x, Err(e) => { if std::env::var_os("AX_DEBUG").is_some() { eprintln!("check error: {e:?}"); } f64::NAN } } }, v)
        }
        "bio_ne_fluctuating_size" => {
            let r = f_ne_fluctuating(match v.get("sizes") { Some(x) => py::from_val(x), None => { return None } });
            finish(r, |i: &Py, o: &Py| -> f64 { match s_4667(i.clone(), o.clone()).and_then(|r| r.f()) { Ok(x) => x, Err(e) => { if std::env::var_os("AX_DEBUG").is_some() { eprintln!("check error: {e:?}"); } f64::NAN } } }, v)
        }
        "bio_ne_sex_ratio" => {
            let r = f_ne_sex_ratio(match v.get("N_m") { Some(x) => py::from_val(x), None => { return None } }, match v.get("N_f") { Some(x) => py::from_val(x), None => { return None } });
            finish(r, |i: &Py, o: &Py| -> f64 { match s_4715(i.clone(), o.clone()).and_then(|r| r.f()) { Ok(x) => x, Err(e) => { if std::env::var_os("AX_DEBUG").is_some() { eprintln!("check error: {e:?}"); } f64::NAN } } }, v)
        }
        "bio_neighbor_joining" => {
            let r = f_neighbor_joining(match v.get("D") { Some(x) => py::from_val(x), None => { return None } }, (|| -> R { Ok(Py::Bool(true)) })().ok()?);
            finish(r, |i: &Py, o: &Py| -> f64 { match s_4835(i.clone(), o.clone()).and_then(|r| r.f()) { Ok(x) => x, Err(e) => { if std::env::var_os("AX_DEBUG").is_some() { eprintln!("check error: {e:?}"); } f64::NAN } } }, v)
        }
        "bio_nucleic_acid_mass" => {
            let r = f_nucleic_mass(match v.get("seq") { Some(x) => py::from_val(x), None => { return None } }, match v.get("molecule") { Some(x) => py::from_val(x), None => { return None } }, match v.get("strands") { Some(x) => py::from_val(x), None => { return None } });
            finish(r, |i: &Py, o: &Py| -> f64 { match s_4924(i.clone(), o.clone()).and_then(|r| r.f()) { Ok(x) => x, Err(e) => { if std::env::var_os("AX_DEBUG").is_some() { eprintln!("check error: {e:?}"); } f64::NAN } } }, v)
        }
        "bio_nucleotide_diversity" => {
            let r = f_nucleotide_diversity(match v.get("seqs") { Some(x) => py::from_val(x), None => { return None } });
            finish(r, |i: &Py, o: &Py| -> f64 { match s_4980(i.clone(), o.clone()).and_then(|r| r.f()) { Ok(x) => x, Err(e) => { if std::env::var_os("AX_DEBUG").is_some() { eprintln!("check error: {e:?}"); } f64::NAN } } }, v)
        }
        "bio_pedigree_inbreeding" => {
            let r = f_pedigree_inbreeding(match v.get("pedigree") { Some(x) => py::from_val(x), None => { return None } }, match v.get("individual") { Some(x) => py::from_val(x), None => { return None } });
            finish(r, |i: &Py, o: &Py| -> f64 { match s_5062(i.clone(), o.clone()).and_then(|r| r.f()) { Ok(x) => x, Err(e) => { if std::env::var_os("AX_DEBUG").is_some() { eprintln!("check error: {e:?}"); } f64::NAN } } }, v)
        }
        "bio_petersen_estimate" => {
            let r = f_petersen_mle(match v.get("M") { Some(x) => py::from_val(x), None => { return None } }, match v.get("C") { Some(x) => py::from_val(x), None => { return None } }, match v.get("R") { Some(x) => py::from_val(x), None => { return None } }, (|| -> R { Ok(Py::Bool(false)) })().ok()?);
            finish(r, |i: &Py, o: &Py| -> f64 { match s_5110(i.clone(), o.clone()).and_then(|r| r.f()) { Ok(x) => x, Err(e) => { if std::env::var_os("AX_DEBUG").is_some() { eprintln!("check error: {e:?}"); } f64::NAN } } }, v)
        }
        "bio_protein_extinction" => {
            let r = f_extinction(match v.get("seq") { Some(x) => py::from_val(x), None => { return None } });
            finish(r, |i: &Py, o: &Py| -> f64 { match s_5189(i.clone(), o.clone()).and_then(|r| r.f()) { Ok(x) => x, Err(e) => { if std::env::var_os("AX_DEBUG").is_some() { eprintln!("check error: {e:?}"); } f64::NAN } } }, v)
        }
        "bio_protein_global_blosum62" => {
            let r = f_protein_global(match v.get("seq1") { Some(x) => py::from_val(x), None => { return None } }, match v.get("seq2") { Some(x) => py::from_val(x), None => { return None } }, match v.get("open") { Some(x) => py::from_val(x), None => { return None } }, match v.get("extend") { Some(x) => py::from_val(x), None => { return None } });
            finish(r, |i: &Py, o: &Py| -> f64 { match s_5309(i.clone(), o.clone()).and_then(|r| r.f()) { Ok(x) => x, Err(e) => { if std::env::var_os("AX_DEBUG").is_some() { eprintln!("check error: {e:?}"); } f64::NAN } } }, v)
        }
        "bio_protein_local_blosum62" => {
            let r = f_protein_local(match v.get("seq1") { Some(x) => py::from_val(x), None => { return None } }, match v.get("seq2") { Some(x) => py::from_val(x), None => { return None } }, match v.get("open") { Some(x) => py::from_val(x), None => { return None } }, match v.get("extend") { Some(x) => py::from_val(x), None => { return None } });
            finish(r, |i: &Py, o: &Py| -> f64 { match s_5429(i.clone(), o.clone()).and_then(|r| r.f()) { Ok(x) => x, Err(e) => { if std::env::var_os("AX_DEBUG").is_some() { eprintln!("check error: {e:?}"); } f64::NAN } } }, v)
        }
        "bio_protein_mass" => {
            let r = f_protein_mass(match v.get("seq") { Some(x) => py::from_val(x), None => { return None } });
            finish(r, |i: &Py, o: &Py| -> f64 { match s_5508(i.clone(), o.clone()).and_then(|r| r.f()) { Ok(x) => x, Err(e) => { if std::env::var_os("AX_DEBUG").is_some() { eprintln!("check error: {e:?}"); } f64::NAN } } }, v)
        }
        "bio_protein_net_charge" => {
            let r = f_net_charge(match v.get("seq") { Some(x) => py::from_val(x), None => { return None } }, match v.get("pH") { Some(x) => py::from_val(x), None => { return None } });
            finish(r, |i: &Py, o: &Py| -> f64 { match s_5603(i.clone(), o.clone()).and_then(|r| r.f()) { Ok(x) => x, Err(e) => { if std::env::var_os("AX_DEBUG").is_some() { eprintln!("check error: {e:?}"); } f64::NAN } } }, v)
        }
        "bio_rarefaction" => {
            let r = f_rarefaction(match v.get("counts") { Some(x) => py::from_val(x), None => { return None } }, match v.get("n") { Some(x) => py::from_val(x), None => { return None } }, (|| -> R { Ok(Py::Bool(false)) })().ok()?);
            finish(r, |i: &Py, o: &Py| -> f64 { match s_5670(i.clone(), o.clone()).and_then(|r| r.f()) { Ok(x) => x, Err(e) => { if std::env::var_os("AX_DEBUG").is_some() { eprintln!("check error: {e:?}"); } f64::NAN } } }, v)
        }
        "bio_recessive_lethal" => {
            let r = f_recessive_lethal(match v.get("q0") { Some(x) => py::from_val(x), None => { return None } }, match v.get("t") { Some(x) => py::from_val(x), None => { return None } });
            finish(r, |i: &Py, o: &Py| -> f64 { match s_5742(i.clone(), o.clone()).and_then(|r| r.f()) { Ok(x) => x, Err(e) => { if std::env::var_os("AX_DEBUG").is_some() { eprintln!("check error: {e:?}"); } f64::NAN } } }, v)
        }
        "bio_recombination_fraction" => {
            let r = f_recombination_fraction(match v.get("d_cM") { Some(x) => py::from_val(x), None => { return None } }, match v.get("map") { Some(x) => py::from_val(x), None => { return None } });
            finish(r, |i: &Py, o: &Py| -> f64 { match s_5806(i.clone(), o.clone()).and_then(|r| r.f()) { Ok(x) => x, Err(e) => { if std::env::var_os("AX_DEBUG").is_some() { eprintln!("check error: {e:?}"); } f64::NAN } } }, v)
        }
        "bio_restriction_sites" => {
            let r = f_restriction_sites(match v.get("seq") { Some(x) => py::from_val(x), None => { return None } }, match v.get("enzyme") { Some(x) => py::from_val(x), None => { return None } });
            finish(r, |i: &Py, o: &Py| -> f64 { match s_5875(i.clone(), o.clone()).and_then(|r| r.f()) { Ok(x) => x, Err(e) => { if std::env::var_os("AX_DEBUG").is_some() { eprintln!("check error: {e:?}"); } f64::NAN } } }, v)
        }
        "bio_reverse_complement" => {
            let r = f_reverse_complement(match v.get("seq") { Some(x) => py::from_val(x), None => { return None } });
            finish(r, |i: &Py, o: &Py| -> f64 { match s_5939(i.clone(), o.clone()).and_then(|r| r.f()) { Ok(x) => x, Err(e) => { if std::env::var_os("AX_DEBUG").is_some() { eprintln!("check error: {e:?}"); } f64::NAN } } }, v)
        }
        "bio_richness_indices" => {
            let r = f_richness_indices(match v.get("counts") { Some(x) => py::from_val(x), None => { return None } }, (|| -> R { Ok(Py::Bool(false)) })().ok()?);
            finish(r, |i: &Py, o: &Py| -> f64 { match s_6015(i.clone(), o.clone()).and_then(|r| r.f()) { Ok(x) => x, Err(e) => { if std::env::var_os("AX_DEBUG").is_some() { eprintln!("check error: {e:?}"); } f64::NAN } } }, v)
        }
        "bio_robinson_foulds" => {
            let r = f_robinson_foulds(match v.get("tree1") { Some(x) => py::from_val(x), None => { return None } }, match v.get("tree2") { Some(x) => py::from_val(x), None => { return None } }, (|| -> R { Ok(Py::Bool(false)) })().ok()?);
            finish(r, |i: &Py, o: &Py| -> f64 { match s_6115(i.clone(), o.clone()).and_then(|r| r.f()) { Ok(x) => x, Err(e) => { if std::env::var_os("AX_DEBUG").is_some() { eprintln!("check error: {e:?}"); } f64::NAN } } }, v)
        }
        "bio_schumacher_eschmeyer" => {
            let r = f_schumacher(match v.get("C") { Some(x) => py::from_val(x), None => { return None } }, match v.get("R") { Some(x) => py::from_val(x), None => { return None } }, match v.get("M") { Some(x) => py::from_val(x), None => { return None } }, (|| -> R { Ok(Py::Bool(false)) })().ok()?);
            finish(r, |i: &Py, o: &Py| -> f64 { match s_6190(i.clone(), o.clone()).and_then(|r| r.f()) { Ok(x) => x, Err(e) => { if std::env::var_os("AX_DEBUG").is_some() { eprintln!("check error: {e:?}"); } f64::NAN } } }, v)
        }
        "bio_secondary_structure_fraction" => {
            let r = f_ss_fraction(match v.get("seq") { Some(x) => py::from_val(x), None => { return None } }, (|| -> R { Ok(Py::tuple(vec![Py::str("helix"), Py::str("turn"), Py::str("sheet")])) })().ok()?);
            finish(r, |i: &Py, o: &Py| -> f64 { match s_6257(i.clone(), o.clone()).and_then(|r| r.f()) { Ok(x) => x, Err(e) => { if std::env::var_os("AX_DEBUG").is_some() { eprintln!("check error: {e:?}"); } f64::NAN } } }, v)
        }
        "bio_semiglobal_alignment_score" => {
            let r = f_semiglobal(match v.get("seq1") { Some(x) => py::from_val(x), None => { return None } }, match v.get("seq2") { Some(x) => py::from_val(x), None => { return None } }, match v.get("match") { Some(x) => py::from_val(x), None => { return None } }, match v.get("mismatch") { Some(x) => py::from_val(x), None => { return None } }, match v.get("gap") { Some(x) => py::from_val(x), None => { return None } });
            finish(r, |i: &Py, o: &Py| -> f64 { match s_6376(i.clone(), o.clone()).and_then(|r| r.f()) { Ok(x) => x, Err(e) => { if std::env::var_os("AX_DEBUG").is_some() { eprintln!("check error: {e:?}"); } f64::NAN } } }, v)
        }
        "bio_shannon_diversity" => {
            let r = f_shannon(match v.get("counts") { Some(x) => py::from_val(x), None => { return None } }, match v.get("base") { Some(x) => py::from_val(x), None => { return None } }, (|| -> R { Ok(Py::Bool(false)) })().ok()?);
            finish(r, |i: &Py, o: &Py| -> f64 { match s_6454(i.clone(), o.clone()).and_then(|r| r.f()) { Ok(x) => x, Err(e) => { if std::env::var_os("AX_DEBUG").is_some() { eprintln!("check error: {e:?}"); } f64::NAN } } }, v)
        }
        "bio_simpson_diversity" => {
            let r = f_simpson(match v.get("counts") { Some(x) => py::from_val(x), None => { return None } }, match v.get("estimator") { Some(x) => py::from_val(x), None => { return None } }, (|| -> R { Ok(Py::Bool(false)) })().ok()?);
            finish(r, |i: &Py, o: &Py| -> f64 { match s_6523(i.clone(), o.clone()).and_then(|r| r.f()) { Ok(x) => x, Err(e) => { if std::env::var_os("AX_DEBUG").is_some() { eprintln!("check error: {e:?}"); } f64::NAN } } }, v)
        }
        "bio_species_area" => {
            let r = f_species_area(match v.get("A") { Some(x) => py::from_val(x), None => { return None } }, match v.get("S") { Some(x) => py::from_val(x), None => { return None } }, (|| -> R { Ok(Py::Bool(false)) })().ok()?);
            finish(r, |i: &Py, o: &Py| -> f64 { match s_6602(i.clone(), o.clone()).and_then(|r| r.f()) { Ok(x) => x, Err(e) => { if std::env::var_os("AX_DEBUG").is_some() { eprintln!("check error: {e:?}"); } f64::NAN } } }, v)
        }
        "bio_tajima_d" => {
            let r = f_tajima_d(match v.get("seqs") { Some(x) => py::from_val(x), None => { return None } }, (|| -> R { Ok(Py::Bool(true)) })().ok()?);
            finish(r, |i: &Py, o: &Py| -> f64 { match s_6680(i.clone(), o.clone()).and_then(|r| r.f()) { Ok(x) => x, Err(e) => { if std::env::var_os("AX_DEBUG").is_some() { eprintln!("check error: {e:?}"); } f64::NAN } } }, v)
        }
        "bio_testcross_recombination" => {
            let r = f_testcross_recombination(match v.get("parental") { Some(x) => py::from_val(x), None => { return None } }, match v.get("recombinant") { Some(x) => py::from_val(x), None => { return None } });
            finish(r, |i: &Py, o: &Py| -> f64 { match s_6749(i.clone(), o.clone()).and_then(|r| r.f()) { Ok(x) => x, Err(e) => { if std::env::var_os("AX_DEBUG").is_some() { eprintln!("check error: {e:?}"); } f64::NAN } } }, v)
        }
        "bio_three_point_testcross" => {
            let r = f_three_point_testcross(match v.get("classes") { Some(x) => py::from_val(x), None => { return None } }, match v.get("counts") { Some(x) => py::from_val(x), None => { return None } });
            finish(r, |i: &Py, o: &Py| -> f64 { match s_6893(i.clone(), o.clone()).and_then(|r| r.f()) { Ok(x) => x, Err(e) => { if std::env::var_os("AX_DEBUG").is_some() { eprintln!("check error: {e:?}"); } f64::NAN } } }, v)
        }
        "bio_tm_nearest_neighbor" => {
            let r = f_tm_nn(match v.get("seq") { Some(x) => py::from_val(x), None => { return None } }, match v.get("Na") { Some(x) => py::from_val(x), None => { return None } }, match v.get("dnac1") { Some(x) => py::from_val(x), None => { return None } }, match v.get("dnac2") { Some(x) => py::from_val(x), None => { return None } });
            finish(r, |i: &Py, o: &Py| -> f64 { match s_7000(i.clone(), o.clone()).and_then(|r| r.f()) { Ok(x) => x, Err(e) => { if std::env::var_os("AX_DEBUG").is_some() { eprintln!("check error: {e:?}"); } f64::NAN } } }, v)
        }
        "bio_tm_wallace" => {
            let r = f_tm_wallace(match v.get("seq") { Some(x) => py::from_val(x), None => { return None } });
            finish(r, |i: &Py, o: &Py| -> f64 { match s_7064(i.clone(), o.clone()).and_then(|r| r.f()) { Ok(x) => x, Err(e) => { if std::env::var_os("AX_DEBUG").is_some() { eprintln!("check error: {e:?}"); } f64::NAN } } }, v)
        }
        "bio_transcription" => {
            let r = f_transcription(match v.get("seq") { Some(x) => py::from_val(x), None => { return None } }, match v.get("strand") { Some(x) => py::from_val(x), None => { return None } });
            finish(r, |i: &Py, o: &Py| -> f64 { match s_7132(i.clone(), o.clone()).and_then(|r| r.f()) { Ok(x) => x, Err(e) => { if std::env::var_os("AX_DEBUG").is_some() { eprintln!("check error: {e:?}"); } f64::NAN } } }, v)
        }
        "bio_translation" => {
            let r = f_translation(match v.get("seq") { Some(x) => py::from_val(x), None => { return None } }, match v.get("stops") { Some(x) => py::from_val(x), None => { return None } });
            finish(r, |i: &Py, o: &Py| -> f64 { match s_7206(i.clone(), o.clone()).and_then(|r| r.f()) { Ok(x) => x, Err(e) => { if std::env::var_os("AX_DEBUG").is_some() { eprintln!("check error: {e:?}"); } f64::NAN } } }, v)
        }
        "bio_tree_count" => {
            let r = f_tree_count(match v.get("n") { Some(x) => py::from_val(x), None => { return None } }, match v.get("rooted") { Some(x) => py::from_val(x), None => { return None } });
            finish(r, |i: &Py, o: &Py| -> f64 { match s_7259(i.clone(), o.clone()).and_then(|r| r.f()) { Ok(x) => x, Err(e) => { if std::env::var_os("AX_DEBUG").is_some() { eprintln!("check error: {e:?}"); } f64::NAN } } }, v)
        }
        "bio_upgma" => {
            let r = f_upgma(match v.get("D") { Some(x) => py::from_val(x), None => { return None } }, (|| -> R { Ok(Py::str("average")) })().ok()?);
            finish(r, |i: &Py, o: &Py| -> f64 { match s_7360(i.clone(), o.clone()).and_then(|r| r.f()) { Ok(x) => x, Err(e) => { if std::env::var_os("AX_DEBUG").is_some() { eprintln!("check error: {e:?}"); } f64::NAN } } }, v)
        }
        "bio_watterson_theta" => {
            let r = f_watterson(match v.get("seqs") { Some(x) => py::from_val(x), None => { return None } }, (|| -> R { Ok(Py::Bool(false)) })().ok()?);
            finish(r, |i: &Py, o: &Py| -> f64 { match s_7429(i.clone(), o.clone()).and_then(|r| r.f()) { Ok(x) => x, Err(e) => { if std::env::var_os("AX_DEBUG").is_some() { eprintln!("check error: {e:?}"); } f64::NAN } } }, v)
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
    v_pos = { let mut out = vec![]; for t6510 in py::iter(&(v_c.clone()))? { let mut v_x: Py; 
v_x = t6510;
if ({ let t6511 = v_x.clone(); Py::Bool({ let t6512 = Py::Int(0i64); if py::gt(&t6511, &t6512)? { true } else { false } }) }).truthy() { out.push(v_x.clone()); }} Py::list(out) };
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
    let t6345 = (Py::tuple(vec![py::len(&v_a.clone())?, py::len(&v_b.clone())?])).items()?;
    if t6345.len() != 2 { return Err(PyErr::Value("unpack".into())); }
    v_n = t6345[0].clone();
    v_m = t6345[1].clone();
    v_o2 = (if (v_wrong_gaps.clone()).truthy() { py::add(&v_op.clone(), &v_ex.clone())? } else { v_op.clone() });
    v_H = { let mut out = vec![]; for t6346 in py::range_iter(&[(py::add(&v_n.clone(), &Py::Int(1i64))?).i()?]) { let mut v__: Py; 
v__ = t6346;
out.push(py::mul(&Py::list(vec![Py::Float(0.0)]), &py::add(&v_m.clone(), &Py::Int(1i64))?)?); } Py::list(out) };
    v_E = { let mut out = vec![]; for t6347 in py::range_iter(&[(py::add(&v_n.clone(), &Py::Int(1i64))?).i()?]) { let mut v__: Py; 
v__ = t6347;
out.push(py::mul(&Py::list(vec![c__NEG()?]), &py::add(&v_m.clone(), &Py::Int(1i64))?)?); } Py::list(out) };
    v_F = { let mut out = vec![]; for t6348 in py::range_iter(&[(py::add(&v_n.clone(), &Py::Int(1i64))?).i()?]) { let mut v__: Py; 
v__ = t6348;
out.push(py::mul(&Py::list(vec![c__NEG()?]), &py::add(&v_m.clone(), &Py::Int(1i64))?)?); } Py::list(out) };
    for t6349 in py::range_iter(&[(Py::Int(1i64)).i()?, (py::add(&v_n.clone(), &Py::Int(1i64))?).i()?]) {
        v_i = t6349;
        if ({ let t6350 = v_mode.clone(); Py::Bool({ let t6351 = Py::str("global"); if py::eq(&t6350, &t6351) { true } else { false } }) }).truthy() {
            let t6352 = py::add(&v_o2.clone(), &py::mul(&py::sub(&v_i.clone(), &Py::Int(1i64))?, &v_ex.clone())?)?;
            py::setitem(&py::getitem(&v_H.clone(), &v_i.clone())?, &Py::Int(0i64), t6352.clone())?;
            py::setitem(&py::getitem(&v_F.clone(), &v_i.clone())?, &Py::Int(0i64), t6352.clone())?;
        }
    }
    for t6353 in py::range_iter(&[(Py::Int(1i64)).i()?, (py::add(&v_m.clone(), &Py::Int(1i64))?).i()?]) {
        v_j = t6353;
        if ({ let t6354 = v_mode.clone(); Py::Bool({ let t6355 = Py::str("global"); if py::eq(&t6354, &t6355) { true } else { false } }) }).truthy() {
            let t6356 = py::add(&v_o2.clone(), &py::mul(&py::sub(&v_j.clone(), &Py::Int(1i64))?, &v_ex.clone())?)?;
            py::setitem(&py::getitem(&v_H.clone(), &Py::Int(0i64))?, &v_j.clone(), t6356.clone())?;
            py::setitem(&py::getitem(&v_E.clone(), &Py::Int(0i64))?, &v_j.clone(), t6356.clone())?;
        }
    }
    v_best = Py::Float(0.0);
    for t6357 in py::range_iter(&[(Py::Int(1i64)).i()?, (py::add(&v_n.clone(), &Py::Int(1i64))?).i()?]) {
        v_i = t6357;
        let t6358 = (Py::tuple(vec![py::getitem(&v_H.clone(), &v_i.clone())?, py::getitem(&v_H.clone(), &py::sub(&v_i.clone(), &Py::Int(1i64))?)?, py::getitem(&v_E.clone(), &v_i.clone())?, py::getitem(&v_F.clone(), &v_i.clone())?, py::getitem(&v_F.clone(), &py::sub(&v_i.clone(), &Py::Int(1i64))?)?])).items()?;
        if t6358.len() != 5 { return Err(PyErr::Value("unpack".into())); }
        v_Hi = t6358[0].clone();
        v_Hp = t6358[1].clone();
        v_Ei = t6358[2].clone();
        v_Fi = t6358[3].clone();
        v_Fp = t6358[4].clone();
        v_ai = py::getitem(&v_a.clone(), &py::sub(&v_i.clone(), &Py::Int(1i64))?)?;
        for t6359 in py::range_iter(&[(Py::Int(1i64)).i()?, (py::add(&v_m.clone(), &Py::Int(1i64))?).i()?]) {
            v_j = t6359;
            v_e = py::minmax(&[py::add(&py::getitem(&v_Ei.clone(), &py::sub(&v_j.clone(), &Py::Int(1i64))?)?, &v_ex.clone())?, py::add(&py::getitem(&v_Hi.clone(), &py::sub(&v_j.clone(), &Py::Int(1i64))?)?, &v_o2.clone())?], true)?;
            v_f = py::minmax(&[py::add(&py::getitem(&v_Fp.clone(), &v_j.clone())?, &v_ex.clone())?, py::add(&py::getitem(&v_Hp.clone(), &v_j.clone())?, &v_o2.clone())?], true)?;
            v_h = py::minmax(&[py::add(&py::getitem(&v_Hp.clone(), &py::sub(&v_j.clone(), &Py::Int(1i64))?)?, &py::call(&v_sc.clone(), &[v_ai.clone(), py::getitem(&v_b.clone(), &py::sub(&v_j.clone(), &Py::Int(1i64))?)?])?)?, v_e.clone(), v_f.clone()], true)?;
            if ({ let t6364 = { let t6360 = v_mode.clone(); Py::Bool({ let t6361 = Py::str("local"); if py::eq(&t6360, &t6361) { true } else { false } }) }; if !t6364.truthy() { t6364 } else { { let t6362 = v_h.clone(); Py::Bool({ let t6363 = Py::Int(0i64); if py::lt(&t6362, &t6363)? { true } else { false } }) } } }).truthy() {
                v_h = Py::Float(0.0);
            }
            let t6365 = (Py::tuple(vec![v_e.clone(), v_f.clone(), v_h.clone()])).items()?;
            if t6365.len() != 3 { return Err(PyErr::Value("unpack".into())); }
            py::setitem(&v_Ei.clone(), &v_j.clone(), t6365[0].clone())?;
            py::setitem(&v_Fi.clone(), &v_j.clone(), t6365[1].clone())?;
            py::setitem(&v_Hi.clone(), &v_j.clone(), t6365[2].clone())?;
            if ({ let t6366 = v_h.clone(); Py::Bool({ let t6367 = v_best.clone(); if py::gt(&t6366, &t6367)? { true } else { false } }) }).truthy() {
                v_best = v_h.clone();
            }
        }
    }
    if ({ let t6368 = v_mode.clone(); Py::Bool({ let t6369 = Py::str("global"); if py::eq(&t6368, &t6369) { true } else { false } }) }).truthy() {
        return Ok(py::getitem(&py::getitem(&v_H.clone(), &v_n.clone())?, &v_m.clone())?);
    }
    if ({ let t6370 = v_mode.clone(); Py::Bool({ let t6371 = Py::str("local"); if py::eq(&t6370, &t6371) { true } else { false } }) }).truthy() {
        return Ok(v_best.clone());
    }
    return Ok(py::minmax(&[py::add(&{ let mut out = vec![]; for t6372 in py::range_iter(&[(py::add(&v_m.clone(), &Py::Int(1i64))?).i()?]) { let mut v_j: Py; 
v_j = t6372;
out.push(py::getitem(&py::getitem(&v_H.clone(), &v_n.clone())?, &v_j.clone())?); } Py::list(out) }, &{ let mut out = vec![]; for t6373 in py::range_iter(&[(py::add(&v_n.clone(), &Py::Int(1i64))?).i()?]) { let mut v_i: Py; 
v_i = t6373;
out.push(py::getitem(&py::getitem(&v_H.clone(), &v_i.clone())?, &v_m.clone())?); } Py::list(out) })?], true)?);
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn f__aligned_dna(mut v_seq1: Py, mut v_seq2: Py) -> R {
    let mut v_a: Py = Py::None;
    let mut v_b: Py = Py::None;
    let t3290 = (Py::tuple(vec![f__seq(v_seq1.clone(), Py::str("seq1"), Py::str("DNA"))?, f__seq(v_seq2.clone(), Py::str("seq2"), Py::str("DNA"))?])).items()?;
    if t3290.len() != 2 { return Err(PyErr::Value("unpack".into())); }
    v_a = t3290[0].clone();
    v_b = t3290[1].clone();
    let _ = f__need({ let t3291 = py::len(&v_a.clone())?; Py::Bool({ let t3292 = py::len(&v_b.clone())?; if py::eq(&t3291, &t3292) { true } else { false } }) }, Py::str("the two aligned sequences must have the same length (no gaps)"))?;
    return Ok(Py::tuple(vec![v_a.clone(), v_b.clone()]));
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn f__alignment(mut v_seqs: Py, mut v_lo: Py) -> R {
    let mut v_r: Py = Py::None;
    let mut v_rows: Py = Py::None;
    let _ = f__need(Py::Bool(matches!(v_seqs.clone(), Py::Str(_))), Py::str("give the aligned sequences separated by spaces"))?;
    v_rows = { let mut out = vec![]; for t7401 in py::iter(&(py::method(&v_seqs.clone(), "split", &[])?))? { let mut v_s: Py; 
v_s = t7401;
out.push(py::method(&v_s.clone(), "upper", &[])?); } Py::list(out) };
    let _ = f__need({ let t7402 = py::len(&v_rows.clone())?; Py::Bool({ let t7403 = v_lo.clone(); if py::ge(&t7402, &t7403)? { true } else { false } }) }, Py::str(&["at least ".to_string(), py::format(&v_lo.clone(), "")?, " sequences are needed".to_string()].concat()))?;
    for t7404 in py::iter(&(v_rows.clone()))? {
        v_r = t7404;
        let _ = f__need({ let t7405 = ax_engine::pyre::module_call("fullmatch", &Py::str("[ACGT]+"), &[v_r.clone()], 0)?; Py::Bool({ let t7406 = Py::None; if !py::is(&t7405, &t7406) { true } else { false } }) }, Py::str("sequences must be DNA (A, C, G, T) without gaps"))?;
    }
    let _ = f__need({ let t7407 = py::len(&py::set(&{ let mut out = vec![]; for t7409 in py::iter(&(v_rows.clone()))? { let mut v_r: Py; 
v_r = t7409;
out.push(py::len(&v_r.clone())?); } Py::list(out) })?)?; Py::Bool({ let t7408 = Py::Int(1i64); if py::eq(&t7407, &t7408) { true } else { false } }) }, Py::str("the sequences must be aligned (equal lengths)"))?;
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
    let t4572 = (Py::tuple(vec![py::call(&v_f.clone(), &[v_a.clone()])?, py::call(&v_f.clone(), &[v_b.clone()])?])).items()?;
    if t4572.len() != 2 { return Err(PyErr::Value("unpack".into())); }
    v_fa = t4572[0].clone();
    v_fb = t4572[1].clone();
    if ({ let t4573 = v_fa.clone(); Py::Bool({ let t4574 = Py::Int(0i64); if py::eq(&t4573, &t4574) { true } else { false } }) }).truthy() {
        return Ok(v_a.clone());
    }
    if ({ let t4575 = v_fb.clone(); Py::Bool({ let t4576 = Py::Int(0i64); if py::eq(&t4575, &t4576) { true } else { false } }) }).truthy() {
        return Ok(v_b.clone());
    }
    let _ = f__need({ let t4577 = { let t4579 = v_fa.clone(); Py::Bool({ let t4580 = Py::Int(0i64); if py::gt(&t4579, &t4580)? { true } else { false } }) }; Py::Bool({ let t4578 = { let t4581 = v_fb.clone(); Py::Bool({ let t4582 = Py::Int(0i64); if py::gt(&t4581, &t4582)? { true } else { false } }) }; if py::ne(&t4577, &t4578) { true } else { false } }) }, Py::str("no root in the interval"))?;
    for t4583 in py::range_iter(&[(v_n.clone()).i()?]) {
        v__ = t4583;
        v_m = py::mul(&Py::Float(0.5), &py::add(&v_a.clone(), &v_b.clone())?)?;
        if ({ let t4588 = { let t4584 = v_m.clone(); Py::Bool({ let t4585 = v_a.clone(); if py::eq(&t4584, &t4585) { true } else { false } }) }; if t4588.truthy() { t4588 } else { { let t4586 = v_m.clone(); Py::Bool({ let t4587 = v_b.clone(); if py::eq(&t4586, &t4587) { true } else { false } }) } } }).truthy() {
            break;
        }
        v_fm = py::call(&v_f.clone(), &[v_m.clone()])?;
        if ({ let t4589 = v_fm.clone(); Py::Bool({ let t4590 = Py::Int(0i64); if py::eq(&t4589, &t4590) { true } else { false } }) }).truthy() {
            return Ok(v_m.clone());
        }
        if ({ let t4591 = { let t4593 = v_fm.clone(); Py::Bool({ let t4594 = Py::Int(0i64); if py::gt(&t4593, &t4594)? { true } else { false } }) }; Py::Bool({ let t4592 = { let t4595 = v_fa.clone(); Py::Bool({ let t4596 = Py::Int(0i64); if py::gt(&t4595, &t4596)? { true } else { false } }) }; if py::eq(&t4591, &t4592) { true } else { false } }) }).truthy() {
            let t4597 = (Py::tuple(vec![v_m.clone(), v_fm.clone()])).items()?;
            if t4597.len() != 2 { return Err(PyErr::Value("unpack".into())); }
            v_a = t4597[0].clone();
            v_fa = t4597[1].clone();
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
    return Ok(py::sub(&py::sum(&{ let mut out = vec![]; for t5507 in py::iter(&(v_s.clone()))? { let mut v_x: Py; 
v_x = t5507;
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
    let t5592 = (Py::tuple(vec![py::dict_from(&c_PK_POS()?, vec![])?, py::dict_from(&c_PK_NEG()?, vec![])?])).items()?;
    if t5592.len() != 2 { return Err(PyErr::Value("unpack".into())); }
    v_pos = t5592[0].clone();
    v_neg = t5592[1].clone();
    if ({ let t5593 = py::getitem(&v_s.clone(), &Py::Int(0i64))?; Py::Bool({ let t5594 = c_PK_NTERM()?; if py::contains(&t5594, &t5593)? { true } else { false } }) }).truthy() {
        py::setitem(&v_pos.clone(), &Py::str("Nterm"), py::getitem(&c_PK_NTERM()?, &py::getitem(&v_s.clone(), &Py::Int(0i64))?)?)?;
    }
    if ({ let t5595 = py::getitem(&v_s.clone(), &py::neg(&Py::Int(1i64))?)?; Py::Bool({ let t5596 = c_PK_CTERM()?; if py::contains(&t5596, &t5595)? { true } else { false } }) }).truthy() {
        py::setitem(&v_neg.clone(), &Py::str("Cterm"), py::getitem(&c_PK_CTERM()?, &py::getitem(&v_s.clone(), &py::neg(&Py::Int(1i64))?)?)?)?;
    }
    v_n = { let mut out = vec![]; for t5597 in py::iter(&(Py::str("KRHDECY")))? { let mut v_aa: Py; 
v_aa = t5597;
out.push((v_aa.clone(), py::float(&py::method(&v_s.clone(), "count", &[v_aa.clone()])?)?)); } Py::dict(out) };
    let t5598 = (if (v_termini.clone()).truthy() { Py::Float(1.0) } else { Py::Float(0.0) });
    py::setitem(&v_n.clone(), &Py::str("Nterm"), t5598.clone())?;
    py::setitem(&v_n.clone(), &Py::str("Cterm"), t5598.clone())?;
    v_plus = Py::Float(0.0);
    for t5599 in py::iter(&(py::method(&v_pos.clone(), "items", &[])?))? {
        let t5600 = (t5599).items()?;
        if t5600.len() != 2 { return Err(PyErr::Value("unpack".into())); }
        v_aa = t5600[0].clone();
        v_pk = t5600[1].clone();
        v_plus = py::iadd(&v_plus, &py::mul(&py::getitem(&v_n.clone(), &v_aa.clone())?, &py::truediv(&Py::Float(1.0), &py::add(&py::pow(&Py::Int(10i64), &py::sub(&v_pH.clone(), &v_pk.clone())?)?, &Py::Float(1.0))?)?)?)?;
    }
    v_minus = Py::Float(0.0);
    for t5601 in py::iter(&(py::method(&v_neg.clone(), "items", &[])?))? {
        let t5602 = (t5601).items()?;
        if t5602.len() != 2 { return Err(PyErr::Value("unpack".into())); }
        v_aa = t5602[0].clone();
        v_pk = t5602[1].clone();
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
    let _ = f__need({ let t6322 = py::mul(&py::len(&v_a.clone())?, &py::len(&v_b.clone())?)?; Py::Bool({ let t6323 = Py::Int(1000000i64); if py::le(&t6322, &t6323)? { true } else { false } }) }, Py::str("sequences too long (length product at most 1 000 000)"))?;
    for t6324 in py::iter(&(py::method(&v_scores.clone(), "items", &[])?))? {
        let t6328 = (t6324).items()?;
        if t6328.len() != 2 { return Err(PyErr::Value("unpack".into())); }
        v_name = t6328[0].clone();
        v_v = t6328[1].clone();
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
    let _ = f__need({ let t6794 = Py::Bool(matches!(v_v.clone(), Py::List(_) | Py::Tuple(_))); if !t6794.truthy() { t6794 } else { { let t6792 = py::len(&v_v.clone())?; Py::Bool({ let t6793 = Py::Int(1i64); if py::ge(&t6792, &t6793)? { true } else { false } }) } } }, Py::str(&[py::format(&v_name.clone(), "")?, " must be a list of counts".to_string()].concat()))?;
    v_out = { let mut out = vec![]; for t6795 in py::iter(&(v_v.clone()))? { let mut v_x: Py; 
v_x = t6795;
out.push(f__whole(v_name.clone(), v_x.clone(), v_lo.clone())?); } Py::list(out) };
    if (v_positive.clone()).truthy() {
        let _ = f__need(Py::Bool(({ let mut out = vec![]; for t6804 in py::iter(&(v_out.clone()))? { let mut v_x: Py; 
v_x = t6804;
out.push({ let t6805 = v_x.clone(); Py::Bool({ let t6806 = Py::Int(0i64); if py::gt(&t6805, &t6806)? { true } else { false } }) }); } Py::list(out) }).items()?.iter().all(|x| x.truthy())), Py::str(&[py::format(&v_name.clone(), "")?, ": every count must be positive".to_string()].concat()))?;
    }
    return Ok(v_out.clone());
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn f__dfact(mut v_k: Py) -> R {
    let mut v_out: Py = Py::None;
    v_out = Py::Int(1i64);
    while ({ let t7257 = v_k.clone(); Py::Bool({ let t7258 = Py::Int(1i64); if py::gt(&t7257, &t7258)? { true } else { false } }) }).truthy() {
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
    let _ = f__need({ let t7306 = { let t7300 = v_lo.clone(); Py::Bool({ let t7301 = v_n.clone(); if py::le(&t7300, &t7301)? { { let t7302 = v_hi.clone(); if py::le(&t7301, &t7302)? { true } else { false } } } else { false } }) }; if !t7306.truthy() { t7306 } else { Py::Bool(({ let mut out = vec![]; for t7303 in py::iter(&(v_D.clone()))? { let mut v_r: Py; 
v_r = t7303;
out.push({ let t7304 = py::len(&v_r.clone())?; Py::Bool({ let t7305 = v_n.clone(); if py::eq(&t7304, &t7305) { true } else { false } }) }); } Py::list(out) }).items()?.iter().all(|x| x.truthy())) } }, Py::str(&[py::format(&v_name.clone(), "")?, " must be a square matrix of ".to_string(), py::format(&v_lo.clone(), "")?, " to ".to_string(), py::format(&v_hi.clone(), "")?, " taxa".to_string()].concat()))?;
    for t7307 in py::range_iter(&[(v_n.clone()).i()?]) {
        v_i = t7307;
        let _ = f__need({ let t7308 = py::getitem(&py::getitem(&v_D.clone(), &v_i.clone())?, &v_i.clone())?; Py::Bool({ let t7309 = Py::Int(0i64); if py::eq(&t7308, &t7309) { true } else { false } }) }, Py::str(&[py::format(&v_name.clone(), "")?, ": the diagonal must be 0".to_string()].concat()))?;
        for t7310 in py::range_iter(&[(py::add(&v_i.clone(), &Py::Int(1i64))?).i()?, (v_n.clone()).i()?]) {
            v_j = t7310;
            let _ = f__need({ let t7311 = py::getitem(&py::getitem(&v_D.clone(), &v_i.clone())?, &v_j.clone())?; Py::Bool({ let t7312 = py::getitem(&py::getitem(&v_D.clone(), &v_j.clone())?, &v_i.clone())?; if py::eq(&t7311, &t7312) { true } else { false } }) }, Py::str(&[py::format(&v_name.clone(), "")?, " must be symmetric".to_string()].concat()))?;
            let _ = f__need({ let t7313 = py::getitem(&py::getitem(&v_D.clone(), &v_i.clone())?, &v_j.clone())?; Py::Bool({ let t7314 = Py::Int(0i64); if py::gt(&t7313, &t7314)? { true } else { false } }) }, Py::str(&[py::format(&v_name.clone(), "")?, ": distances between different taxa must be positive".to_string()].concat()))?;
        }
    }
    return Ok({ let mut out = vec![]; for t7315 in py::iter(&(v_D.clone()))? { let mut v_r: Py; 
v_r = t7315;
out.push({ let mut out = vec![]; for t7316 in py::iter(&(v_r.clone()))? { let mut v_x: Py; 
v_x = t7316;
out.push(py::float(&v_x.clone())?); } Py::list(out) }); } Py::list(out) });
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn f__euler_lotka(mut v_F: Py, mut v_S: Py) -> R {
    let mut v_l: Py = Py::None;
    let mut v_s: Py = Py::None;
    v_l = Py::list(vec![Py::Float(1.0)]);
    for t3467 in py::iter(&(v_S.clone()))? {
        v_s = t3467;
        let _ = py::method(&v_l.clone(), "append", &[py::mul(&py::getitem(&v_l.clone(), &py::neg(&Py::Int(1i64))?)?, &v_s.clone())?])?;
    }
    return Ok(v_l.clone());
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn f__fin(mut v_name: Py, mut v_vals: Py) -> R {
    let mut v_v: Py = Py::None;
    for t6965 in py::iter(&(v_vals.clone()))? {
        v_v = t6965;
        let _ = f__need({ let t6967 = Py::Bool(matches!(v_v.clone(), Py::Int(_) | Py::Big(_) | Py::Bool(_) | Py::Float(_))); if !t6967.truthy() { t6967 } else { { let t6966 = py::not(&Py::Bool(matches!(v_v.clone(), Py::Bool(_)))); if !t6966.truthy() { t6966 } else { py::math1("isfinite", &v_v.clone())? } } } }, Py::str(&[py::format(&v_name.clone(), "")?, " must be a finite number".to_string()].concat()))?;
    }
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn f__flip_case(mut v_s: Py) -> R {
    let mut v_ch: Py = Py::None;
    let mut v_out: Py = Py::None;
    v_out = Py::str("");
    for t6845 in py::iter(&(v_s.clone()))? {
        v_ch = t6845;
        if (py::method(&v_ch.clone(), "isupper", &[])?).truthy() {
            v_out = py::iadd(&v_out, &py::method(&v_ch.clone(), "lower", &[])?)?;
        } else {
            v_out = py::iadd(&v_out, &py::method(&v_ch.clone(), "upper", &[])?)?;
        }
    }
    return Ok(v_out.clone());
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
    let _ = f__need({ let t5722 = (if (v_lo_open.clone()).truthy() { { let t5714 = v_v.clone(); Py::Bool({ let t5715 = Py::Int(0i64); if py::gt(&t5714, &t5715)? { true } else { false } }) } } else { { let t5716 = v_v.clone(); Py::Bool({ let t5717 = Py::Int(0i64); if py::ge(&t5716, &t5717)? { true } else { false } }) } }); if !t5722.truthy() { t5722 } else { (if (v_hi_open.clone()).truthy() { { let t5718 = v_v.clone(); Py::Bool({ let t5719 = Py::Int(1i64); if py::lt(&t5718, &t5719)? { true } else { false } }) } } else { { let t5720 = v_v.clone(); Py::Bool({ let t5721 = Py::Int(1i64); if py::le(&t5720, &t5721)? { true } else { false } }) } }) } }, Py::str(&[py::format(&v_name.clone(), "")?, " must be a fraction between 0 and 1".to_string()].concat()))?;
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
    for t2803 in py::range_iter(&[(Py::Int(60i64)).i()?]) {
        v__ = t2803;
        let t2804 = (py::call(&v_model.clone(), &[v_p.clone()])?).items()?;
        if t2804.len() != 2 { return Err(PyErr::Value("unpack".into())); }
        v_f = t2804[0].clone();
        v_J = t2804[1].clone();
        v_r = { let mut out = vec![]; for t2805 in py::iter(&(py::zip(&[v_ys.clone(), v_f.clone()])?))? { let mut v_a: Py; let mut v_y: Py; 
let t2806 = (t2805).items()?;
if t2806.len() != 2 { return Err(PyErr::Value("unpack".into())); }
v_y = t2806[0].clone();
v_a = t2806[1].clone();
out.push(py::sub(&v_y.clone(), &v_a.clone())?); } Py::list(out) };
        v_A = { let mut out = vec![]; for t2807 in py::range_iter(&[(Py::Int(3i64)).i()?]) { let mut v_i: Py; 
v_i = t2807;
out.push({ let mut out = vec![]; for t2808 in py::range_iter(&[(Py::Int(3i64)).i()?]) { let mut v_j: Py; 
v_j = t2808;
out.push(py::fsum(&{ let mut out = vec![]; for t2809 in py::iter(&(v_J.clone()))? { let mut v_row: Py; 
v_row = t2809;
out.push(py::mul(&py::getitem(&v_row.clone(), &v_i.clone())?, &py::getitem(&v_row.clone(), &v_j.clone())?)?); } Py::list(out) })?); } Py::list(out) }); } Py::list(out) };
        v_g = { let mut out = vec![]; for t2810 in py::range_iter(&[(Py::Int(3i64)).i()?]) { let mut v_i: Py; 
v_i = t2810;
out.push(py::fsum(&{ let mut out = vec![]; for t2811 in py::iter(&(py::zip(&[v_J.clone(), v_r.clone()])?))? { let mut v_ri: Py; let mut v_row: Py; 
let t2812 = (t2811).items()?;
if t2812.len() != 2 { return Err(PyErr::Value("unpack".into())); }
v_row = t2812[0].clone();
v_ri = t2812[1].clone();
out.push(py::mul(&py::getitem(&v_row.clone(), &v_i.clone())?, &v_ri.clone())?); } Py::list(out) })?); } Py::list(out) };
        let t2813: Result<Option<Py>, PyErr> = (|| -> Result<Option<Py>, PyErr> {
                v_d = f__solve3(v_A.clone(), v_g.clone())?;
                #[allow(unreachable_code)]
                Ok(None)
            })();
        match t2813 {
            Ok(Some(r)) => { return Ok(r); }
            Ok(None) => {}
            Err(e) if (matches!(e, PyErr::Value(_))) => {
                return Ok(v_p.clone());
            }
            Err(e) => return Err(e),
        }
        v_q = { let mut out = vec![]; for t2814 in py::range_iter(&[(Py::Int(3i64)).i()?]) { let mut v_i: Py; 
v_i = t2814;
out.push(py::add(&py::getitem(&v_p.clone(), &v_i.clone())?, &py::getitem(&v_d.clone(), &v_i.clone())?)?); } Py::list(out) };
        if (py::not(&py::call(&v_ok.clone(), &[v_q.clone()])?)).truthy() {
            return Ok(v_p.clone());
        }
        v_size = py::minmax(&[{ let mut out = vec![]; for t2815 in py::range_iter(&[(Py::Int(3i64)).i()?]) { let mut v_i: Py; 
v_i = t2815;
out.push(py::truediv(&py::abs(&py::getitem(&v_d.clone(), &v_i.clone())?)?, &py::minmax(&[py::abs(&py::getitem(&v_q.clone(), &v_i.clone())?)?, Py::Float(1e-300)], true)?)?); } Py::list(out) }], true)?;
        v_p = v_q.clone();
        if ({ let t2820 = { let t2816 = v_size.clone(); Py::Bool({ let t2817 = Py::Float(1e-16); if py::lt(&t2816, &t2817)? { true } else { false } }) }; if t2820.truthy() { t2820 } else { { let t2818 = v_size.clone(); Py::Bool({ let t2819 = v_last.clone(); if py::ge(&t2818, &t2819)? { true } else { false } }) } } }).truthy() {
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
            let t2823 = (v_p.clone()).items()?;
            if t2823.len() != 3 { return Err(PyErr::Value("unpack".into())); }
            v_Vm = t2823[0].clone();
            v_K = t2823[1].clone();
            v_n = t2823[2].clone();
            let t2824 = (Py::tuple(vec![Py::list(vec![]), Py::list(vec![])])).items()?;
            if t2824.len() != 2 { return Err(PyErr::Value("unpack".into())); }
            v_f = t2824[0].clone();
            v_J = t2824[1].clone();
            for t2825 in py::iter(&(v_S.clone()))? {
                v_s = t2825;
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
    let _ = f__need({ let t4191 = Py::Bool(matches!(v_S.clone(), Py::List(_))); if !t4191.truthy() { t4191 } else { { let t4190 = Py::Bool(matches!(v_v.clone(), Py::List(_))); if !t4190.truthy() { t4190 } else { { let t4189 = { let t4185 = py::len(&v_S.clone())?; Py::Bool({ let t4186 = py::len(&v_v.clone())?; if py::eq(&t4185, &t4186) { true } else { false } }) }; if !t4189.truthy() { t4189 } else { { let t4187 = py::len(&v_S.clone())?; Py::Bool({ let t4188 = v_lo.clone(); if py::ge(&t4187, &t4188)? { true } else { false } }) } } } } } } }, Py::str(&["S and v must be lists of equal length (at least ".to_string(), py::format(&v_lo.clone(), "")?, " points)".to_string()].concat()))?;
    for t4192 in py::iter(&(py::add(&v_S.clone(), &v_v.clone())?))? {
        v_x = t4192;
        let _ = f__pos(Py::str("every concentration and rate"), Py::tuple(vec![v_x.clone()]))?;
    }
    return Ok(Py::tuple(vec![{ let mut out = vec![]; for t4199 in py::iter(&(v_S.clone()))? { let mut v_x: Py; 
v_x = t4199;
out.push(py::float(&v_x.clone())?); } Py::list(out) }, { let mut out = vec![]; for t4200 in py::iter(&(v_v.clone()))? { let mut v_x: Py; 
v_x = t4200;
out.push(py::float(&v_x.clone())?); } Py::list(out) }]));
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn f__leslie_in(mut v_F: Py, mut v_S: Py) -> R {
    let mut v_f: Py = Py::None;
    let mut v_s: Py = Py::None;
    let _ = f__need({ let t3567 = Py::Bool(matches!(v_F.clone(), Py::List(_))); if !t3567.truthy() { t3567 } else { { let t3566 = Py::Bool(matches!(v_S.clone(), Py::List(_))); if !t3566.truthy() { t3566 } else { { let t3565 = { let t3561 = py::len(&v_F.clone())?; Py::Bool({ let t3562 = Py::Int(2i64); if py::ge(&t3561, &t3562)? { true } else { false } }) }; if !t3565.truthy() { t3565 } else { { let t3563 = py::len(&v_S.clone())?; Py::Bool({ let t3564 = py::sub(&py::len(&v_F.clone())?, &Py::Int(1i64))?; if py::eq(&t3563, &t3564) { true } else { false } }) } } } } } } }, Py::str("F needs one fecundity per age class (at least 2) and S one survival fewer"))?;
    for t3568 in py::iter(&(v_F.clone()))? {
        v_f = t3568;
        let _ = f__fin(Py::str("F"), Py::tuple(vec![v_f.clone()]))?;
        let _ = f__need({ let t3572 = v_f.clone(); Py::Bool({ let t3573 = Py::Int(0i64); if py::ge(&t3572, &t3573)? { true } else { false } }) }, Py::str("fecundities must be zero or positive"))?;
    }
    for t3574 in py::iter(&(v_S.clone()))? {
        v_s = t3574;
        let _ = f__fin(Py::str("S"), Py::tuple(vec![v_s.clone()]))?;
        let _ = f__need({ let t3575 = Py::Int(0i64); Py::Bool({ let t3576 = v_s.clone(); if py::lt(&t3575, &t3576)? { { let t3577 = Py::Int(1i64); if py::le(&t3576, &t3577)? { true } else { false } } } else { false } }) }, Py::str("survival probabilities must lie in (0, 1]"))?;
    }
    let _ = f__need(Py::Bool(({ let mut out = vec![]; for t3578 in py::iter(&(v_F.clone()))? { let mut v_f: Py; 
v_f = t3578;
out.push({ let t3579 = v_f.clone(); Py::Bool({ let t3580 = Py::Int(0i64); if py::gt(&t3579, &t3580)? { true } else { false } }) }); } Py::list(out) }).items()?.iter().any(|x| x.truthy())), Py::str("at least one age class must reproduce"))?;
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
            let t2765 = (py::call(&v_model.clone(), &[v_pp.clone()])?).items()?;
            if t2765.len() != 2 { return Err(PyErr::Value("unpack".into())); }
            v_f = t2765[0].clone();
            v__ = t2765[1].clone();
            return Ok(py::fsum(&{ let mut out = vec![]; for t2766 in py::iter(&(py::zip(&[v_ys.clone(), v_f.clone()])?))? { let mut v_a: Py; let mut v_y: Py; 
let t2767 = (t2766).items()?;
if t2767.len() != 2 { return Err(PyErr::Value("unpack".into())); }
v_y = t2767[0].clone();
v_a = t2767[1].clone();
out.push(py::pow(&py::sub(&v_y.clone(), &v_a.clone())?, &Py::Int(2i64))?); } Py::list(out) })?);
            #[allow(unreachable_code)]
            Ok(Py::None)
        }
        }) };
    let t2768 = (Py::tuple(vec![Py::Float(0.001), py::call(&v_cost.clone(), &[v_p.clone()])?])).items()?;
    if t2768.len() != 2 { return Err(PyErr::Value("unpack".into())); }
    v_lam = t2768[0].clone();
    v_c = t2768[1].clone();
    for t2769 in py::range_iter(&[(Py::Int(2000i64)).i()?]) {
        v__ = t2769;
        let t2770 = (py::call(&v_model.clone(), &[v_p.clone()])?).items()?;
        if t2770.len() != 2 { return Err(PyErr::Value("unpack".into())); }
        v_f = t2770[0].clone();
        v_J = t2770[1].clone();
        v_r = { let mut out = vec![]; for t2771 in py::iter(&(py::zip(&[v_ys.clone(), v_f.clone()])?))? { let mut v_a: Py; let mut v_y: Py; 
let t2772 = (t2771).items()?;
if t2772.len() != 2 { return Err(PyErr::Value("unpack".into())); }
v_y = t2772[0].clone();
v_a = t2772[1].clone();
out.push(py::sub(&v_y.clone(), &v_a.clone())?); } Py::list(out) };
        v_A = { let mut out = vec![]; for t2773 in py::range_iter(&[(Py::Int(3i64)).i()?]) { let mut v_i: Py; 
v_i = t2773;
out.push({ let mut out = vec![]; for t2774 in py::range_iter(&[(Py::Int(3i64)).i()?]) { let mut v_j: Py; 
v_j = t2774;
out.push(py::fsum(&{ let mut out = vec![]; for t2775 in py::iter(&(v_J.clone()))? { let mut v_row: Py; 
v_row = t2775;
out.push(py::mul(&py::getitem(&v_row.clone(), &v_i.clone())?, &py::getitem(&v_row.clone(), &v_j.clone())?)?); } Py::list(out) })?); } Py::list(out) }); } Py::list(out) };
        v_g = { let mut out = vec![]; for t2776 in py::range_iter(&[(Py::Int(3i64)).i()?]) { let mut v_i: Py; 
v_i = t2776;
out.push(py::fsum(&{ let mut out = vec![]; for t2777 in py::iter(&(py::zip(&[v_J.clone(), v_r.clone()])?))? { let mut v_ri: Py; let mut v_row: Py; 
let t2778 = (t2777).items()?;
if t2778.len() != 2 { return Err(PyErr::Value("unpack".into())); }
v_row = t2778[0].clone();
v_ri = t2778[1].clone();
out.push(py::mul(&py::getitem(&v_row.clone(), &v_i.clone())?, &v_ri.clone())?); } Py::list(out) })?); } Py::list(out) };
        while (Py::Bool(true)).truthy() {
            v_Ad = { let mut out = vec![]; for t2779 in py::range_iter(&[(Py::Int(3i64)).i()?]) { let mut v_i: Py; 
v_i = t2779;
out.push({ let mut out = vec![]; for t2780 in py::range_iter(&[(Py::Int(3i64)).i()?]) { let mut v_j: Py; 
v_j = t2780;
out.push(py::add(&py::getitem(&py::getitem(&v_A.clone(), &v_i.clone())?, &v_j.clone())?, &(if ({ let t2781 = v_i.clone(); Py::Bool({ let t2782 = v_j.clone(); if py::eq(&t2781, &t2782) { true } else { false } }) }).truthy() { py::mul(&v_lam.clone(), &py::getitem(&py::getitem(&v_A.clone(), &v_i.clone())?, &v_i.clone())?)? } else { Py::Float(0.0) }))?); } Py::list(out) }); } Py::list(out) };
            let t2794: Result<Option<Py>, PyErr> = (|| -> Result<Option<Py>, PyErr> {
                    v_d = f__solve3(v_Ad.clone(), v_g.clone())?;
                    #[allow(unreachable_code)]
                    Ok(None)
                })();
            match t2794 {
                Ok(Some(r)) => { return Ok(r); }
                Ok(None) => {}
                Err(e) if (matches!(e, PyErr::Value(_))) => {
                    v_d = Py::None;
                }
                Err(e) => return Err(e),
            }
            v_q = (if (v_d.clone()).truthy() { { let mut out = vec![]; for t2795 in py::range_iter(&[(Py::Int(3i64)).i()?]) { let mut v_i: Py; 
v_i = t2795;
out.push(py::add(&py::getitem(&v_p.clone(), &v_i.clone())?, &py::getitem(&v_d.clone(), &v_i.clone())?)?); } Py::list(out) } } else { Py::None });
            if ({ let t2796 = v_q.clone(); if !t2796.truthy() { t2796 } else { py::call(&v_ok.clone(), &[v_q.clone()])? } }).truthy() {
                v_cq = py::call(&v_cost.clone(), &[v_q.clone()])?;
                if ({ let t2797 = v_cq.clone(); Py::Bool({ let t2798 = v_c.clone(); if py::lt(&t2797, &t2798)? { true } else { false } }) }).truthy() {
                    v_small = { let t2799 = py::minmax(&[{ let mut out = vec![]; for t2801 in py::range_iter(&[(Py::Int(3i64)).i()?]) { let mut v_i: Py; 
v_i = t2801;
out.push(py::truediv(&py::abs(&py::getitem(&v_d.clone(), &v_i.clone())?)?, &py::minmax(&[py::abs(&py::getitem(&v_p.clone(), &v_i.clone())?)?, Py::Float(1e-300)], true)?)?); } Py::list(out) }], true)?; Py::Bool({ let t2800 = Py::Float(1e-15); if py::lt(&t2799, &t2800)? { true } else { false } }) };
                    let t2802 = (Py::tuple(vec![v_q.clone(), v_cq.clone()])).items()?;
                    if t2802.len() != 2 { return Err(PyErr::Value("unpack".into())); }
                    v_p = t2802[0].clone();
                    v_c = t2802[1].clone();
                    v_lam = py::minmax(&[py::truediv(&v_lam.clone(), &Py::Int(10i64))?, Py::Float(1e-15)], true)?;
                    if (v_small.clone()).truthy() {
                        return Ok(f__gn(v_model.clone(), v_p.clone(), v_ys.clone(), v_ok.clone())?);
                    }
                    break;
                }
            }
            v_lam = py::mul(&v_lam, &Py::Int(10i64))?;
            if ({ let t2821 = v_lam.clone(); Py::Bool({ let t2822 = Py::Float(1e+20); if py::gt(&t2821, &t2822)? { true } else { false } }) }).truthy() {
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
    if ({ let t4497 = v_z.clone(); Py::Bool({ let t4498 = Py::Int(0i64); if py::gt(&t4497, &t4498)? { true } else { false } }) }).truthy() {
        return Ok((if ({ let t4499 = v_z.clone(); Py::Bool({ let t4500 = Py::Int(30i64); if py::gt(&t4499, &t4500)? { true } else { false } }) }).truthy() { py::add(&v_z.clone(), &py::math1("log1p", &py::neg(&py::math1("exp", &py::neg(&v_z.clone())?)?)?)?)? } else { py::math1("log", &py::math1("expm1", &v_z.clone())?)? }));
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
            return Ok((if ({ let t6374 = v_x.clone(); Py::Bool({ let t6375 = v_y.clone(); if py::eq(&t6374, &t6375) { true } else { false } }) }).truthy() { v_match.clone() } else { v_mismatch.clone() }));
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
    v_x = { let mut out = vec![]; for t4207 in py::iter(&(v_S.clone()))? { let mut v_s: Py; 
v_s = t4207;
out.push(py::truediv(&v_s.clone(), &py::add(&v_K.clone(), &v_s.clone())?)?); } Py::list(out) };
    v_xd = { let mut out = vec![]; for t4208 in py::iter(&(v_S.clone()))? { let mut v_s: Py; 
v_s = t4208;
out.push(py::truediv(&py::neg(&v_s.clone())?, &py::pow(&py::add(&v_K.clone(), &v_s.clone())?, &Py::Int(2i64))?)?); } Py::list(out) };
    v_svx = py::fsum(&{ let mut out = vec![]; for t4209 in py::iter(&(py::zip(&[v_v.clone(), v_x.clone()])?))? { let mut v_a: Py; let mut v_b: Py; 
let t4210 = (t4209).items()?;
if t4210.len() != 2 { return Err(PyErr::Value("unpack".into())); }
v_a = t4210[0].clone();
v_b = t4210[1].clone();
out.push(py::mul(&v_a.clone(), &v_b.clone())?); } Py::list(out) })?;
    v_sxx = py::fsum(&{ let mut out = vec![]; for t4211 in py::iter(&(v_x.clone()))? { let mut v_a: Py; 
v_a = t4211;
out.push(py::mul(&v_a.clone(), &v_a.clone())?); } Py::list(out) })?;
    v_svxd = py::fsum(&{ let mut out = vec![]; for t4212 in py::iter(&(py::zip(&[v_v.clone(), v_xd.clone()])?))? { let mut v_a: Py; let mut v_b: Py; 
let t4213 = (t4212).items()?;
if t4213.len() != 2 { return Err(PyErr::Value("unpack".into())); }
v_a = t4213[0].clone();
v_b = t4213[1].clone();
out.push(py::mul(&v_a.clone(), &v_b.clone())?); } Py::list(out) })?;
    v_sxxd = py::fsum(&{ let mut out = vec![]; for t4214 in py::iter(&(py::zip(&[v_x.clone(), v_xd.clone()])?))? { let mut v_a: Py; let mut v_b: Py; 
let t4215 = (t4214).items()?;
if t4215.len() != 2 { return Err(PyErr::Value("unpack".into())); }
v_a = t4215[0].clone();
v_b = t4215[1].clone();
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
    let t6576 = (Py::tuple(vec![py::truediv(&py::fsum(&v_x.clone())?, &v_n.clone())?, py::truediv(&py::fsum(&v_y.clone())?, &v_n.clone())?])).items()?;
    if t6576.len() != 2 { return Err(PyErr::Value("unpack".into())); }
    v_mx = t6576[0].clone();
    v_my = t6576[1].clone();
    v_sxx = py::fsum(&{ let mut out = vec![]; for t6577 in py::iter(&(v_x.clone()))? { let mut v_a: Py; 
v_a = t6577;
out.push(py::pow(&py::sub(&v_a.clone(), &v_mx.clone())?, &Py::Int(2i64))?); } Py::list(out) })?;
    v_syy = py::fsum(&{ let mut out = vec![]; for t6578 in py::iter(&(v_y.clone()))? { let mut v_b: Py; 
v_b = t6578;
out.push(py::pow(&py::sub(&v_b.clone(), &v_my.clone())?, &Py::Int(2i64))?); } Py::list(out) })?;
    v_sxy = py::fsum(&{ let mut out = vec![]; for t6579 in py::iter(&(py::zip(&[v_x.clone(), v_y.clone()])?))? { let mut v_a: Py; let mut v_b: Py; 
let t6580 = (t6579).items()?;
if t6580.len() != 2 { return Err(PyErr::Value("unpack".into())); }
v_a = t6580[0].clone();
v_b = t6580[1].clone();
out.push(py::mul(&py::sub(&v_a.clone(), &v_mx.clone())?, &py::sub(&v_b.clone(), &v_my.clone())?)?); } Py::list(out) })?;
    let _ = f__need({ let t6581 = v_sxx.clone(); Py::Bool({ let t6582 = Py::Int(0i64); if py::gt(&t6581, &t6582)? { true } else { false } }) }, Py::str("the x values must not all be equal"))?;
    v_b = py::truediv(&v_sxy.clone(), &v_sxx.clone())?;
    v_r2 = (if ({ let t6583 = v_syy.clone(); Py::Bool({ let t6584 = Py::Int(0i64); if py::gt(&t6583, &t6584)? { true } else { false } }) }).truthy() { py::truediv(&py::mul(&v_sxy.clone(), &v_sxy.clone())?, &py::mul(&v_sxx.clone(), &v_syy.clone())?)? } else { Py::Float(1.0) });
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
    for t4071 in py::range_iter(&[(Py::Int(3i64)).i()?]) {
        v_f = t4071;
        v_start = Py::None;
        for t4072 in py::range_iter(&[(v_f.clone()).i()?, (py::sub(&py::len(&v_s.clone())?, &Py::Int(2i64))?).i()?, (Py::Int(3i64)).i()?]) {
            v_i = t4072;
            v_c = py::slice(&v_s.clone(), Some((v_i.clone()).i()?), Some((py::add(&v_i.clone(), &Py::Int(3i64))?).i()?), None)?;
            if ({ let t4073 = v_start.clone(); Py::Bool({ let t4074 = Py::None; if py::is(&t4073, &t4074) { true } else { false } }) }).truthy() {
                if ({ let t4075 = v_c.clone(); Py::Bool({ let t4076 = Py::str("ATG"); if py::eq(&t4075, &t4076) { true } else { false } }) }).truthy() {
                    v_start = v_i.clone();
                }
            } else {
                if ({ let t4077 = v_c.clone(); Py::Bool({ let t4078 = c_STOPS()?; if py::contains(&t4078, &t4077)? { true } else { false } }) }).truthy() {
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
    return Ok({ let mut out = vec![]; for t7415 in py::iter(&(py::method(&v_kw.clone(), "items", &[])?))? { let mut v_k: Py; let mut v_v: Py; 
let t7416 = (t7415).items()?;
if t7416.len() != 2 { return Err(PyErr::Value("unpack".into())); }
v_k = t7416[0].clone();
v_v = t7416[1].clone();
out.push((v_k.clone(), f__r(v_v.clone())?)); } Py::dict(out) });
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn f__pairwise_k(mut v_rows: Py) -> R {
    let mut v_n: Py = Py::None;
    let mut v_tot: Py = Py::None;
    v_n = py::len(&v_rows.clone())?;
    v_tot = py::sum(&{ let mut out = vec![]; for t6657 in py::range_iter(&[(v_n.clone()).i()?]) { let mut v_i: Py; 
v_i = t6657;
for t6658 in py::range_iter(&[(py::add(&v_i.clone(), &Py::Int(1i64))?).i()?, (v_n.clone()).i()?]) { let mut v_j: Py; 
v_j = t6658;
out.push(py::sum(&{ let mut out = vec![]; for t6659 in py::iter(&(py::zip(&[py::getitem(&v_rows.clone(), &v_i.clone())?, py::getitem(&v_rows.clone(), &v_j.clone())?])?))? { let mut v_x: Py; let mut v_y: Py; 
let t6660 = (t6659).items()?;
if t6660.len() != 2 { return Err(PyErr::Value("unpack".into())); }
v_x = t6660[0].clone();
v_y = t6660[1].clone();
if ({ let t6661 = v_x.clone(); Py::Bool({ let t6662 = v_y.clone(); if py::ne(&t6661, &t6662) { true } else { false } }) }).truthy() { out.push(Py::Int(1i64)); }} Py::list(out) }, None)?); }} Py::list(out) }, None)?;
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
            let _ = f__need({ let t6056 = v_m.clone(); Py::Bool({ let t6057 = Py::None; if !py::is(&t6056, &t6057) { true } else { false } }) }, Py::str("a branch length must be a number"))?;
            let t6058 = v_pos.clone(); let t6059 = Py::Int(0i64);
            py::setitem(&t6058, &t6059, py::iadd(&py::getitem(&t6058, &t6059)?, &py::method(&v_m.clone(), "end", &[])?)?)?;
            return Ok(py::float(&py::method(&v_m.clone(), "group", &[Py::Int(0i64)])?)?);
            #[allow(unreachable_code)]
            Ok(Py::None)
        }
        }) };
    let rec6060 = py::rec_cell();
    v_node = { let v_num = v_num.clone(); let v_pos = v_pos.clone(); let v_s = v_s.clone(); let v_text = v_text.clone(); let rec6060 = rec6060.clone(); py::func(move |a: &[Py]| -> R {
            #[allow(unused_mut, unused_variables, unused_assignments)]
{
            let v_node: Py = py::rec_get(&rec6060)?;
            let mut v_kids: Py = Py::None;
            let mut v_length: Py = Py::None;
            let mut v_m: Py = Py::None;
            let mut v_name: Py = Py::None;
            v_kids = Py::list(vec![]);
            if ({ let t6065 = { let t6061 = py::getitem(&v_pos.clone(), &Py::Int(0i64))?; Py::Bool({ let t6062 = py::len(&v_s.clone())?; if py::lt(&t6061, &t6062)? { true } else { false } }) }; if !t6065.truthy() { t6065 } else { { let t6063 = py::getitem(&v_s.clone(), &py::getitem(&v_pos.clone(), &Py::Int(0i64))?)?; Py::Bool({ let t6064 = Py::str("("); if py::eq(&t6063, &t6064) { true } else { false } }) } } }).truthy() {
                let t6066 = v_pos.clone(); let t6067 = Py::Int(0i64);
                py::setitem(&t6066, &t6067, py::iadd(&py::getitem(&t6066, &t6067)?, &Py::Int(1i64))?)?;
                let _ = py::method(&v_kids.clone(), "append", &[py::call(&v_node.clone(), &[])?])?;
                while ({ let t6074 = { let t6070 = py::getitem(&v_pos.clone(), &Py::Int(0i64))?; Py::Bool({ let t6071 = py::len(&v_s.clone())?; if py::lt(&t6070, &t6071)? { true } else { false } }) }; if !t6074.truthy() { t6074 } else { { let t6072 = py::getitem(&v_s.clone(), &py::getitem(&v_pos.clone(), &Py::Int(0i64))?)?; Py::Bool({ let t6073 = Py::str(","); if py::eq(&t6072, &t6073) { true } else { false } }) } } }).truthy() {
                    let t6068 = v_pos.clone(); let t6069 = Py::Int(0i64);
                    py::setitem(&t6068, &t6069, py::iadd(&py::getitem(&t6068, &t6069)?, &Py::Int(1i64))?)?;
                    let _ = py::method(&v_kids.clone(), "append", &[py::call(&v_node.clone(), &[])?])?;
                }
                let _ = f__need({ let t6079 = { let t6075 = py::getitem(&v_pos.clone(), &Py::Int(0i64))?; Py::Bool({ let t6076 = py::len(&v_s.clone())?; if py::lt(&t6075, &t6076)? { true } else { false } }) }; if !t6079.truthy() { t6079 } else { { let t6077 = py::getitem(&v_s.clone(), &py::getitem(&v_pos.clone(), &Py::Int(0i64))?)?; Py::Bool({ let t6078 = Py::str(")"); if py::eq(&t6077, &t6078) { true } else { false } }) } } }, Py::str("unbalanced parentheses in the tree"))?;
                let t6080 = v_pos.clone(); let t6081 = Py::Int(0i64);
                py::setitem(&t6080, &t6081, py::iadd(&py::getitem(&t6080, &t6081)?, &Py::Int(1i64))?)?;
            }
            v_m = ax_engine::pyre::module_call("match", &Py::str("[A-Za-z0-9_.]*"), &[py::slice(&v_s.clone(), Some((py::getitem(&v_pos.clone(), &Py::Int(0i64))?).i()?), None, None)?], 0)?;
            v_name = py::method(&v_m.clone(), "group", &[Py::Int(0i64)])?;
            let t6082 = v_pos.clone(); let t6083 = Py::Int(0i64);
            py::setitem(&t6082, &t6083, py::iadd(&py::getitem(&t6082, &t6083)?, &py::method(&v_m.clone(), "end", &[])?)?)?;
            v_length = Py::None;
            if ({ let t6088 = { let t6084 = py::getitem(&v_pos.clone(), &Py::Int(0i64))?; Py::Bool({ let t6085 = py::len(&v_s.clone())?; if py::lt(&t6084, &t6085)? { true } else { false } }) }; if !t6088.truthy() { t6088 } else { { let t6086 = py::getitem(&v_s.clone(), &py::getitem(&v_pos.clone(), &Py::Int(0i64))?)?; Py::Bool({ let t6087 = Py::str(":"); if py::eq(&t6086, &t6087) { true } else { false } }) } } }).truthy() {
                let t6089 = v_pos.clone(); let t6090 = Py::Int(0i64);
                py::setitem(&t6089, &t6090, py::iadd(&py::getitem(&t6089, &t6090)?, &Py::Int(1i64))?)?;
                v_length = py::call(&v_num.clone(), &[])?;
            }
            if (py::not(&v_kids.clone())).truthy() {
                let _ = f__need({ let t6091 = ax_engine::pyre::module_call("fullmatch", &Py::str("[A-Za-z][A-Za-z0-9_]*"), &[{ let t6093 = v_name.clone(); if t6093.truthy() { t6093 } else { Py::str("") } }], 0)?; Py::Bool({ let t6092 = Py::None; if !py::is(&t6091, &t6092) { true } else { false } }) }, Py::str("every leaf needs a name"))?;
            }
            return Ok(Py::tuple(vec![v_name.clone(), v_length.clone(), v_kids.clone()]));
            #[allow(unreachable_code)]
            Ok(Py::None)
        }
        }) };
    py::rec_set(&rec6060, &v_node);
    v_root = py::call(&v_node.clone(), &[])?;
    let _ = f__need({ let t6094 = py::getitem(&v_pos.clone(), &Py::Int(0i64))?; Py::Bool({ let t6095 = py::len(&v_s.clone())?; if py::eq(&t6094, &t6095) { true } else { false } }) }, Py::str("unreadable text after the tree"))?;
    v_leaves = Py::list(vec![]);
    let rec6096 = py::rec_cell();
    v_walk = { let v_leaves = v_leaves.clone(); let v_node = v_node.clone(); let v_num = v_num.clone(); let v_pos = v_pos.clone(); let v_root = v_root.clone(); let v_s = v_s.clone(); let v_text = v_text.clone(); let rec6096 = rec6096.clone(); py::func(move |a: &[Py]| -> R {
            #[allow(unused_mut, unused_variables, unused_assignments)]
{
            let mut v_t: Py = a.get(0).cloned().unwrap_or(Py::None);
            let v_walk: Py = py::rec_get(&rec6096)?;
            let mut v_k: Py = Py::None;
            if (py::not(&py::getitem(&v_t.clone(), &Py::Int(2i64))?)).truthy() {
                let _ = py::method(&v_leaves.clone(), "append", &[py::getitem(&v_t.clone(), &Py::Int(0i64))?])?;
            }
            for t6097 in py::iter(&(py::getitem(&v_t.clone(), &Py::Int(2i64))?))? {
                v_k = t6097;
                let _ = py::call(&v_walk.clone(), &[v_k.clone()])?;
            }
            #[allow(unreachable_code)]
            Ok(Py::None)
        }
        }) };
    py::rec_set(&rec6096, &v_walk);
    let _ = py::call(&v_walk.clone(), &[v_root.clone()])?;
    let _ = f__need({ let t6098 = py::len(&py::set(&v_leaves.clone())?)?; Py::Bool({ let t6099 = py::len(&v_leaves.clone())?; if py::eq(&t6098, &t6099) { true } else { false } }) }, Py::str("leaf names must be unique"))?;
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
    for t5021 in py::iter(&(v_pedigree.clone()))? {
        v_e = t5021;
        let _ = f__need({ let t5022 = py::len(&v_e.clone())?; Py::Bool({ let t5023 = Py::Int(2i64); if py::eq(&t5022, &t5023) { true } else { false } }) }, Py::str("pedigree edges are parent-child pairs without weights"))?;
        let t5024 = (v_e.clone()).items()?;
        if t5024.len() != 2 { return Err(PyErr::Value("unpack".into())); }
        v_p = t5024[0].clone();
        v_c = t5024[1].clone();
        let _ = f__need({ let t5025 = v_p.clone(); Py::Bool({ let t5026 = v_c.clone(); if py::ne(&t5025, &t5026) { true } else { false } }) }, Py::str("an individual cannot be its own parent"))?;
        let _ = py::method(&v_parents.clone(), "setdefault", &[v_c.clone(), Py::list(vec![])])?;
        let _ = f__need({ let t5027 = v_p.clone(); Py::Bool({ let t5028 = py::getitem(&v_parents.clone(), &v_c.clone())?; if !py::contains(&t5028, &t5027)? { true } else { false } }) }, Py::str(&[py::format(&v_p.clone(), "")?, " is listed twice as a parent of ".to_string(), py::format(&v_c.clone(), "")?].concat()))?;
        let _ = py::method(&py::getitem(&v_parents.clone(), &v_c.clone())?, "append", &[v_p.clone()])?;
    }
    for t5029 in py::iter(&(py::method(&v_parents.clone(), "items", &[])?))? {
        let t5032 = (t5029).items()?;
        if t5032.len() != 2 { return Err(PyErr::Value("unpack".into())); }
        v_c = t5032[0].clone();
        v_ps = t5032[1].clone();
        let _ = f__need({ let t5030 = py::len(&v_ps.clone())?; Py::Bool({ let t5031 = Py::Int(2i64); if py::eq(&t5030, &t5031) { true } else { false } }) }, Py::str(&[py::format(&v_c.clone(), "")?, " must have exactly two parents listed (founders have none)".to_string()].concat()))?;
    }
    let t5033 = (Py::tuple(vec![py::set(&Py::list(vec![]))?, py::set(&Py::list(vec![]))?])).items()?;
    if t5033.len() != 2 { return Err(PyErr::Value("unpack".into())); }
    v_seen = t5033[0].clone();
    v_done = t5033[1].clone();
    let rec5034 = py::rec_cell();
    v_visit = { let v_c = v_c.clone(); let v_done = v_done.clone(); let v_e = v_e.clone(); let v_p = v_p.clone(); let v_parents = v_parents.clone(); let v_pedigree = v_pedigree.clone(); let v_ps = v_ps.clone(); let v_seen = v_seen.clone(); let rec5034 = rec5034.clone(); py::func(move |a: &[Py]| -> R {
            #[allow(unused_mut, unused_variables, unused_assignments)]
{
            let mut v_x: Py = a.get(0).cloned().unwrap_or(Py::None);
            let v_visit: Py = py::rec_get(&rec5034)?;
            let mut v_p: Py = Py::None;
            if ({ let t5035 = v_x.clone(); Py::Bool({ let t5036 = v_done.clone(); if py::contains(&t5036, &t5035)? { true } else { false } }) }).truthy() {
                return Ok(Py::None);
            }
            let _ = f__need({ let t5037 = v_x.clone(); Py::Bool({ let t5038 = v_seen.clone(); if !py::contains(&t5038, &t5037)? { true } else { false } }) }, Py::str("the pedigree has a cycle"))?;
            let _ = py::method(&v_seen.clone(), "add", &[v_x.clone()])?;
            for t5039 in py::iter(&(py::method(&v_parents.clone(), "get", &[v_x.clone(), Py::list(vec![])])?))? {
                v_p = t5039;
                let _ = py::call(&v_visit.clone(), &[v_p.clone()])?;
            }
            let _ = py::method(&v_done.clone(), "add", &[v_x.clone()])?;
            #[allow(unreachable_code)]
            Ok(Py::None)
        }
        }) };
    py::rec_set(&rec5034, &v_visit);
    for t5040 in py::iter(&(Py::list((v_parents.clone()).items()?)))? {
        v_c = t5040;
        let _ = py::call(&v_visit.clone(), &[v_c.clone()])?;
    }
    return Ok(v_parents.clone());
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn f__pos(mut v_name: Py, mut v_vals: Py) -> R {
    let mut v_v: Py = Py::None;
    for t6959 in py::iter(&(v_vals.clone()))? {
        v_v = t6959;
        let _ = f__need({ let t6964 = Py::Bool(matches!(v_v.clone(), Py::Int(_) | Py::Big(_) | Py::Bool(_) | Py::Float(_))); if !t6964.truthy() { t6964 } else { { let t6963 = py::not(&Py::Bool(matches!(v_v.clone(), Py::Bool(_)))); if !t6963.truthy() { t6963 } else { { let t6962 = py::math1("isfinite", &v_v.clone())?; if !t6962.truthy() { t6962 } else { { let t6960 = v_v.clone(); Py::Bool({ let t6961 = Py::Int(0i64); if py::gt(&t6960, &t6961)? { true } else { false } }) } } } } } } }, Py::str(&[py::format(&v_name.clone(), "")?, " must be a positive number".to_string()].concat()))?;
    }
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn f__r(mut v_x: Py) -> R {
    if ({ let t7420 = Py::Bool(matches!(v_x.clone(), Py::Bool(_))); if t7420.truthy() { t7420 } else { { let t7419 = Py::Bool(matches!(v_x.clone(), Py::Str(_))); if t7419.truthy() { t7419 } else { { let t7417 = v_x.clone(); Py::Bool({ let t7418 = Py::None; if py::is(&t7417, &t7418) { true } else { false } }) } } } } }).truthy() {
        return Ok(v_x.clone());
    }
    if (Py::Bool(matches!(v_x.clone(), Py::Int(_) | Py::Big(_) | Py::Bool(_)))).truthy() {
        return Ok(v_x.clone());
    }
    if (Py::Bool(matches!(v_x.clone(), Py::Dict(_)))).truthy() {
        return Ok({ let mut out = vec![]; for t7421 in py::iter(&(py::method(&v_x.clone(), "items", &[])?))? { let mut v_k: Py; let mut v_v: Py; 
let t7422 = (t7421).items()?;
if t7422.len() != 2 { return Err(PyErr::Value("unpack".into())); }
v_k = t7422[0].clone();
v_v = t7422[1].clone();
out.push((v_k.clone(), f__r(v_v.clone())?)); } Py::dict(out) });
    }
    if (Py::Bool(matches!(v_x.clone(), Py::List(_) | Py::Tuple(_)))).truthy() {
        return Ok({ let mut out = vec![]; for t7423 in py::iter(&(v_x.clone()))? { let mut v_v: Py; 
v_v = t7423;
out.push(f__r(v_v.clone())?); } Py::list(out) });
    }
    if ({ let t7426 = { let t7424 = v_x.clone(); Py::Bool({ let t7425 = Py::Int(0i64); if py::eq(&t7424, &t7425) { true } else { false } }) }; if t7426.truthy() { t7426 } else { py::not(&py::math1("isfinite", &v_x.clone())?) } }).truthy() {
        return Ok((if ({ let t7427 = v_x.clone(); Py::Bool({ let t7428 = Py::Int(0i64); if py::eq(&t7427, &t7428) { true } else { false } }) }).truthy() { Py::Float(0.0) } else { v_x.clone() }));
    }
    return Ok(py::add(&py::float(&Py::str(&[py::format(&v_x.clone(), ".12g")?].concat()))?, &Py::Float(0.0))?);
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn f__rank_desc(mut v_cnt: Py) -> R {
    let mut v_i: Py = Py::None;
    let mut v_idx: Py = Py::None;
    let mut v_k: Py = Py::None;
    v_idx = Py::list(vec![]);
    for t6833 in py::range_iter(&[(py::len(&v_cnt.clone())?).i()?]) {
        v_i = t6833;
        v_k = py::len(&v_idx.clone())?;
        while ({ let t6838 = { let t6834 = v_k.clone(); Py::Bool({ let t6835 = Py::Int(0i64); if py::gt(&t6834, &t6835)? { true } else { false } }) }; if !t6838.truthy() { t6838 } else { { let t6836 = py::getitem(&v_cnt.clone(), &py::getitem(&v_idx.clone(), &py::sub(&v_k.clone(), &Py::Int(1i64))?)?)?; Py::Bool({ let t6837 = py::getitem(&v_cnt.clone(), &v_i.clone())?; if py::lt(&t6836, &t6837)? { true } else { false } }) } } }).truthy() {
            v_k = py::sub(&v_k, &Py::Int(1i64))?;
        }
        let _ = py::method(&v_idx.clone(), "insert", &[v_k.clone(), v_i.clone()])?;
    }
    return Ok(v_idx.clone());
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
    if ({ let t7435 = Py::Bool(matches!(v_a.clone(), Py::Dict(_))); if t7435.truthy() { t7435 } else { Py::Bool(matches!(v_b.clone(), Py::Dict(_))) } }).truthy() {
        if ({ let t7439 = py::not(&{ let t7436 = Py::Bool(matches!(v_a.clone(), Py::Dict(_))); if !t7436.truthy() { t7436 } else { Py::Bool(matches!(v_b.clone(), Py::Dict(_))) } }); if t7439.truthy() { t7439 } else { { let t7437 = py::method(&v_a.clone(), "keys", &[])?; Py::Bool({ let t7438 = py::method(&v_b.clone(), "keys", &[])?; if py::ne(&t7437, &t7438) { true } else { false } }) } } }).truthy() {
            return Ok(Py::Float(f64::INFINITY));
        }
        return Ok(py::minmax(&[py::add(&{ let mut out = vec![]; for t7440 in py::iter(&(v_a.clone()))? { let mut v_k: Py; 
v_k = t7440;
out.push(f__rel(py::getitem(&v_a.clone(), &v_k.clone())?, py::getitem(&v_b.clone(), &v_k.clone())?)?); } Py::list(out) }, &Py::list(vec![Py::Float(0.0)]))?], true)?);
    }
    if ({ let t7441 = Py::Bool(matches!(v_a.clone(), Py::List(_) | Py::Tuple(_))); if t7441.truthy() { t7441 } else { Py::Bool(matches!(v_b.clone(), Py::List(_) | Py::Tuple(_))) } }).truthy() {
        if ({ let t7445 = py::not(&{ let t7442 = Py::Bool(matches!(v_a.clone(), Py::List(_) | Py::Tuple(_))); if !t7442.truthy() { t7442 } else { Py::Bool(matches!(v_b.clone(), Py::List(_) | Py::Tuple(_))) } }); if t7445.truthy() { t7445 } else { { let t7443 = py::len(&v_a.clone())?; Py::Bool({ let t7444 = py::len(&v_b.clone())?; if py::ne(&t7443, &t7444) { true } else { false } }) } } }).truthy() {
            return Ok(Py::Float(f64::INFINITY));
        }
        return Ok(py::minmax(&[py::add(&{ let mut out = vec![]; for t7446 in py::iter(&(py::zip(&[v_a.clone(), v_b.clone()])?))? { let mut v_x: Py; let mut v_y: Py; 
let t7447 = (t7446).items()?;
if t7447.len() != 2 { return Err(PyErr::Value("unpack".into())); }
v_x = t7447[0].clone();
v_y = t7447[1].clone();
out.push(f__rel(v_x.clone(), v_y.clone())?); } Py::list(out) }, &Py::list(vec![Py::Float(0.0)]))?], true)?);
    }
    if ({ let t7454 = Py::Bool(matches!(v_a.clone(), Py::Str(_) | Py::Bool(_))); if t7454.truthy() { t7454 } else { { let t7453 = Py::Bool(matches!(v_b.clone(), Py::Str(_) | Py::Bool(_))); if t7453.truthy() { t7453 } else { { let t7452 = { let t7448 = v_a.clone(); Py::Bool({ let t7449 = Py::None; if py::is(&t7448, &t7449) { true } else { false } }) }; if t7452.truthy() { t7452 } else { { let t7450 = v_b.clone(); Py::Bool({ let t7451 = Py::None; if py::is(&t7450, &t7451) { true } else { false } }) } } } } } } }).truthy() {
        return Ok((if ({ let t7455 = v_a.clone(); Py::Bool({ let t7456 = v_b.clone(); if py::eq(&t7455, &t7456) { true } else { false } }) }).truthy() { Py::Float(0.0) } else { Py::Float(f64::INFINITY) }));
    }
    let t7458: Result<Option<Py>, PyErr> = (|| -> Result<Option<Py>, PyErr> {
            let t7457 = (Py::tuple(vec![py::float(&v_a.clone())?, py::float(&v_b.clone())?])).items()?;
            if t7457.len() != 2 { return Err(PyErr::Value("unpack".into())); }
            v_a = t7457[0].clone();
            v_b = t7457[1].clone();
            #[allow(unreachable_code)]
            Ok(None)
        })();
    match t7458 {
        Ok(Some(r)) => { return Ok(r); }
        Ok(None) => {}
        Err(e) if ((matches!(e, PyErr::Type(_)) || matches!(&e, PyErr::Other(c) if c.as_str() == "TypeError"))) || (matches!(e, PyErr::Value(_))) => {
            return Ok(Py::Float(f64::INFINITY));
        }
        Err(e) => return Err(e),
    }
    if (py::not(&{ let t7459 = py::math1("isfinite", &v_a.clone())?; if !t7459.truthy() { t7459 } else { py::math1("isfinite", &v_b.clone())? } })).truthy() {
        return Ok((if ({ let t7460 = v_a.clone(); Py::Bool({ let t7461 = v_b.clone(); if py::eq(&t7460, &t7461) { true } else { false } }) }).truthy() { Py::Float(0.0) } else { Py::Float(f64::INFINITY) }));
    }
    v_d = py::abs(&py::sub(&v_a.clone(), &v_b.clone())?)?;
    if ({ let t7462 = v_d.clone(); Py::Bool({ let t7463 = Py::Int(0i64); if py::eq(&t7462, &t7463) { true } else { false } }) }).truthy() {
        return Ok(Py::Float(0.0));
    }
    if ({ let t7468 = { let t7464 = py::abs(&v_a.clone())?; Py::Bool({ let t7465 = Py::Float(1e-15); if py::lt(&t7464, &t7465)? { true } else { false } }) }; if !t7468.truthy() { t7468 } else { { let t7466 = py::abs(&v_b.clone())?; Py::Bool({ let t7467 = Py::Float(1e-15); if py::lt(&t7466, &t7467)? { true } else { false } }) } } }).truthy() {
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
    let _ = f__need({ let t7175 = Py::Bool(matches!(v_s.clone(), Py::Str(_))); if !t7175.truthy() { t7175 } else { { let t7173 = v_s.clone(); Py::Bool({ let t7174 = Py::str(""); if py::ne(&t7173, &t7174) { true } else { false } }) } } }, Py::str(&[py::format(&v_name.clone(), "")?, " must be a sequence of letters".to_string()].concat()))?;
    v_u = py::method(&v_s.clone(), "upper", &[])?;
    let _ = f__need({ let t7176 = py::len(&v_u.clone())?; Py::Bool({ let t7177 = Py::Int(100000i64); if py::le(&t7176, &t7177)? { true } else { false } }) }, Py::str(&[py::format(&v_name.clone(), "")?, ": at most 100000 letters".to_string()].concat()))?;
    if ({ let t7178 = v_kind.clone(); Py::Bool({ let t7179 = Py::str("DNA"); if py::eq(&t7178, &t7179) { true } else { false } }) }).truthy() {
        let _ = f__need({ let t7180 = ax_engine::pyre::module_call("fullmatch", &Py::str("[ACGT]+"), &[v_u.clone()], 0)?; Py::Bool({ let t7181 = Py::None; if !py::is(&t7180, &t7181) { true } else { false } }) }, Py::str(&[py::format(&v_name.clone(), "")?, " must be a DNA sequence of A, C, G and T only".to_string()].concat()))?;
    } else {
        if ({ let t7182 = v_kind.clone(); Py::Bool({ let t7183 = Py::str("RNA"); if py::eq(&t7182, &t7183) { true } else { false } }) }).truthy() {
            let _ = f__need({ let t7184 = ax_engine::pyre::module_call("fullmatch", &Py::str("[ACGU]+"), &[v_u.clone()], 0)?; Py::Bool({ let t7185 = Py::None; if !py::is(&t7184, &t7185) { true } else { false } }) }, Py::str(&[py::format(&v_name.clone(), "")?, " must be an RNA sequence of A, C, G and U only".to_string()].concat()))?;
        } else {
            if ({ let t7186 = v_kind.clone(); Py::Bool({ let t7187 = Py::str("NUC"); if py::eq(&t7186, &t7187) { true } else { false } }) }).truthy() {
                let _ = f__need({ let t7188 = ax_engine::pyre::module_call("fullmatch", &Py::str("[ACGT]+|[ACGU]+"), &[v_u.clone()], 0)?; Py::Bool({ let t7189 = Py::None; if !py::is(&t7188, &t7189) { true } else { false } }) }, Py::str(&[py::format(&v_name.clone(), "")?, " must be a DNA (A, C, G, T) or RNA (A, C, G, U) sequence, not a mixture".to_string()].concat()))?;
            } else {
                if ({ let t7190 = v_kind.clone(); Py::Bool({ let t7191 = Py::str("PROTEIN"); if py::eq(&t7190, &t7191) { true } else { false } }) }).truthy() {
                    let _ = f__need({ let t7192 = ax_engine::pyre::module_call("fullmatch", &py::add(&py::add(&Py::str("["), &c_AMINO()?)?, &Py::str("]+"))?, &[v_u.clone()], 0)?; Py::Bool({ let t7193 = Py::None; if !py::is(&t7192, &t7193) { true } else { false } }) }, Py::str(&[py::format(&v_name.clone(), "")?, " must be a protein sequence of the 20 standard one-letter amino-acid codes".to_string()].concat()))?;
                } else {
                    let _ = f__need({ let t7194 = ax_engine::pyre::module_call("fullmatch", &Py::str("[A-Z]+"), &[v_u.clone()], 0)?; Py::Bool({ let t7195 = Py::None; if !py::is(&t7194, &t7195) { true } else { false } }) }, Py::str(&[py::format(&v_name.clone(), "")?, " must be letters only".to_string()].concat()))?;
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
    v_M = { let mut out = vec![]; for t2783 in py::range_iter(&[(Py::Int(3i64)).i()?]) { let mut v_i: Py; 
v_i = t2783;
out.push(py::add(&Py::list((py::getitem(&v_A.clone(), &v_i.clone())?).items()?), &Py::list(vec![py::getitem(&v_b.clone(), &v_i.clone())?]))?); } Py::list(out) };
    for t2784 in py::range_iter(&[(Py::Int(3i64)).i()?]) {
        v_c = t2784;
        v_p = py::minmax_key(&[Py::list((py::range_iter(&[(v_c.clone()).i()?, (Py::Int(3i64)).i()?])).into_iter().collect::<Vec<Py>>())], &{ let v_A = v_A.clone(); let v_M = v_M.clone(); let v_b = v_b.clone(); let v_c = v_c.clone(); py::func(move |a: &[Py]| -> R {
            #[allow(unused_mut, unused_variables, unused_assignments)]
{
            let mut v_r: Py = a.get(0).cloned().unwrap_or(Py::None);
            return Ok(py::abs(&py::getitem(&py::getitem(&v_M.clone(), &v_r.clone())?, &v_c.clone())?)?);
            #[allow(unreachable_code)]
            Ok(Py::None)
        }
        }) }, true)?;
        if ({ let t2785 = py::abs(&py::getitem(&py::getitem(&v_M.clone(), &v_p.clone())?, &v_c.clone())?)?; Py::Bool({ let t2786 = Py::Float(1e-300); if py::lt(&t2785, &t2786)? { true } else { false } }) }).truthy() {
            return Err(py::value_error(py::to_str(&Py::str("singular normal equations"))));
        }
        let t2787 = (Py::tuple(vec![py::getitem(&v_M.clone(), &v_p.clone())?, py::getitem(&v_M.clone(), &v_c.clone())?])).items()?;
        if t2787.len() != 2 { return Err(PyErr::Value("unpack".into())); }
        py::setitem(&v_M.clone(), &v_c.clone(), t2787[0].clone())?;
        py::setitem(&v_M.clone(), &v_p.clone(), t2787[1].clone())?;
        for t2788 in py::range_iter(&[(py::add(&v_c.clone(), &Py::Int(1i64))?).i()?, (Py::Int(3i64)).i()?]) {
            v_r = t2788;
            v_f = py::truediv(&py::getitem(&py::getitem(&v_M.clone(), &v_r.clone())?, &v_c.clone())?, &py::getitem(&py::getitem(&v_M.clone(), &v_c.clone())?, &v_c.clone())?)?;
            for t2789 in py::range_iter(&[(v_c.clone()).i()?, (Py::Int(4i64)).i()?]) {
                v_k = t2789;
                let t2790 = py::getitem(&v_M.clone(), &v_r.clone())?; let t2791 = v_k.clone();
                py::setitem(&t2790, &t2791, py::sub(&py::getitem(&t2790, &t2791)?, &py::mul(&v_f.clone(), &py::getitem(&py::getitem(&v_M.clone(), &v_c.clone())?, &v_k.clone())?)?)?)?;
            }
        }
    }
    v_x = py::mul(&Py::list(vec![Py::Float(0.0)]), &Py::Int(3i64))?;
    for t2792 in py::iter(&(Py::tuple(vec![Py::Int(2i64), Py::Int(1i64), Py::Int(0i64)])))? {
        v_r = t2792;
        py::setitem(&v_x.clone(), &v_r.clone(), py::truediv(&py::sub(&py::getitem(&py::getitem(&v_M.clone(), &v_r.clone())?, &Py::Int(3i64))?, &py::sum(&{ let mut out = vec![]; for t2793 in py::range_iter(&[(py::add(&v_r.clone(), &Py::Int(1i64))?).i()?, (Py::Int(3i64)).i()?]) { let mut v_k: Py; 
v_k = t2793;
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
    let rec6106 = py::rec_cell();
    v_walk = { let v_allx = v_allx.clone(); let v_first = v_first.clone(); let v_leaves = v_leaves.clone(); let v_out = v_out.clone(); let v_root = v_root.clone(); let rec6106 = rec6106.clone(); py::func(move |a: &[Py]| -> R {
            #[allow(unused_mut, unused_variables, unused_assignments)]
{
            let mut v_t: Py = a.get(0).cloned().unwrap_or(Py::None);
            let v_walk: Py = py::rec_get(&rec6106)?;
            let mut v_s: Py = Py::None;
            let mut v_side: Py = Py::None;
            if (py::not(&py::getitem(&v_t.clone(), &Py::Int(2i64))?)).truthy() {
                return Ok(py::set(&Py::list(vec![py::getitem(&v_t.clone(), &Py::Int(0i64))?]))?);
            }
            v_s = py::method(&py::set(&Py::list(vec![]))?, "union", &({ let mut v: Vec<Py> = vec![]; v.extend(({ let mut out = vec![]; for t6107 in py::iter(&(py::getitem(&v_t.clone(), &Py::Int(2i64))?))? { let mut v_k: Py; 
v_k = t6107;
out.push(py::call(&v_walk.clone(), &[v_k.clone()])?); } Py::list(out) }).items()?); v })[..])?;
            v_side = (if ({ let t6108 = v_first.clone(); Py::Bool({ let t6109 = v_s.clone(); if py::contains(&t6109, &t6108)? { true } else { false } }) }).truthy() { py::sub(&v_allx.clone(), &v_s.clone())? } else { v_s.clone() });
            if ({ let t6110 = Py::Int(2i64); Py::Bool({ let t6111 = py::len(&v_side.clone())?; if py::le(&t6110, &t6111)? { { let t6112 = py::sub(&py::len(&v_allx.clone())?, &Py::Int(2i64))?; if py::le(&t6111, &t6112)? { true } else { false } } } else { false } }) }).truthy() {
                let _ = py::method(&v_out.clone(), "add", &[v_side.clone()])?;
            }
            return Ok(v_s.clone());
            #[allow(unreachable_code)]
            Ok(Py::None)
        }
        }) };
    py::rec_set(&rec6106, &v_walk);
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
    return Ok({ let t7321 = py::sum(&{ let mut out = vec![]; for t7323 in py::iter(&(v_vals.clone()))? { let mut v_v: Py; 
v_v = t7323;
if ({ let t7324 = py::abs(&py::sub(&v_v.clone(), &v_best.clone())?)?; Py::Bool({ let t7325 = py::mul(&Py::Float(1e-12), &v_scale.clone())?; if py::le(&t7324, &t7325)? { true } else { false } }) }).truthy() { out.push(Py::Int(1i64)); }} Py::list(out) }, None)?; Py::Bool({ let t7322 = Py::Int(1i64); if py::gt(&t7321, &t7322)? { true } else { false } }) });
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
    let t6989 = (Py::tuple(vec![Py::Float(0.0), Py::Float(0.0)])).items()?;
    if t6989.len() != 2 { return Err(PyErr::Value("unpack".into())); }
    v_dh = t6989[0].clone();
    v_ds = t6989[1].clone();
    v_t = c_DNA_NN3()?;
    v_dh = py::iadd(&v_dh, &py::getitem(&py::getitem(&v_t.clone(), &Py::str("init"))?, &Py::Int(0i64))?)?;
    v_ds = py::iadd(&v_ds, &py::getitem(&py::getitem(&v_t.clone(), &Py::str("init"))?, &Py::Int(1i64))?)?;
    v_key = (if ({ let t6994 = { let t6990 = Py::str("G"); Py::Bool({ let t6991 = v_s.clone(); if py::contains(&t6991, &t6990)? { true } else { false } }) }; if t6994.truthy() { t6994 } else { { let t6992 = Py::str("C"); Py::Bool({ let t6993 = v_s.clone(); if py::contains(&t6993, &t6992)? { true } else { false } }) } } }).truthy() { Py::str("init_oneG/C") } else { Py::str("init_allA/T") });
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
    for t6995 in py::range_iter(&[(py::sub(&py::len(&v_s.clone())?, &Py::Int(1i64))?).i()?]) {
        v_i = t6995;
        v_nb = py::add(&py::add(&py::slice(&v_s.clone(), Some((v_i.clone()).i()?), Some((py::add(&v_i.clone(), &Py::Int(2i64))?).i()?), None)?, &Py::str("/"))?, &py::slice(&v_c.clone(), Some((v_i.clone()).i()?), Some((py::add(&v_i.clone(), &Py::Int(2i64))?).i()?), None)?)?;
        v_v = (if ({ let t6996 = v_nb.clone(); Py::Bool({ let t6997 = v_t.clone(); if py::contains(&t6997, &t6996)? { true } else { false } }) }).truthy() { py::getitem(&v_t.clone(), &v_nb.clone())? } else { py::getitem(&v_t.clone(), &py::slice(&v_nb.clone(), None, None, Some((py::neg(&Py::Int(1i64))?).i()?))?)? });
        v_dh = py::iadd(&v_dh, &py::getitem(&v_v.clone(), &Py::Int(0i64))?)?;
        v_ds = py::iadd(&v_ds, &py::getitem(&v_v.clone(), &Py::Int(1i64))?)?;
    }
    v_selfcomp = { let t6998 = v_s.clone(); Py::Bool({ let t6999 = f__revcomp(v_s.clone())?; if py::eq(&t6998, &t6999) { true } else { false } }) };
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
fn f__tp_classes(mut v_classes: Py, mut v_counts: Py) -> R {
    let mut v_c: Py = Py::None;
    let mut v_cnt: Py = Py::None;
    let mut v_i: Py = Py::None;
    let mut v_j: Py = Py::None;
    let mut v_letters: Py = Py::None;
    let mut v_n: Py = Py::None;
    let mut v_names: Py = Py::None;
    let mut v_nm: Py = Py::None;
    v_names = (if (Py::Bool(matches!(v_classes.clone(), Py::Str(_)))).truthy() { py::method(&v_classes.clone(), "split", &[])? } else { Py::list(vec![]) });
    let _ = f__need({ let t6790 = py::len(&v_names.clone())?; Py::Bool({ let t6791 = Py::Int(8i64); if py::eq(&t6790, &t6791) { true } else { false } }) }, Py::str("give the 8 offspring classes, three letters each, separated by spaces"))?;
    v_cnt = f__counts(Py::str("counts"), v_counts.clone(), Py::Int(0i64), Py::Bool(false))?;
    let _ = f__need({ let t6807 = py::len(&v_cnt.clone())?; Py::Bool({ let t6808 = Py::Int(8i64); if py::eq(&t6807, &t6808) { true } else { false } }) }, Py::str("give one count per class (8 counts)"))?;
    v_letters = Py::str("");
    for t6809 in py::iter(&(v_names.clone()))? {
        v_nm = t6809;
        let _ = f__need({ let t6812 = { let t6810 = py::len(&v_nm.clone())?; Py::Bool({ let t6811 = Py::Int(3i64); if py::eq(&t6810, &t6811) { true } else { false } }) }; if !t6812.truthy() { t6812 } else { py::method(&v_nm.clone(), "isalpha", &[])? } }, Py::str("each class is written as three letters, one per gene"))?;
        if ({ let t6813 = v_letters.clone(); Py::Bool({ let t6814 = Py::str(""); if py::eq(&t6813, &t6814) { true } else { false } }) }).truthy() {
            v_letters = py::method(&v_nm.clone(), "upper", &[])?;
        }
        let _ = f__need({ let t6815 = py::method(&v_nm.clone(), "upper", &[])?; Py::Bool({ let t6816 = v_letters.clone(); if py::eq(&t6815, &t6816) { true } else { false } }) }, Py::str("every class must name the same three genes in the same order"))?;
    }
    let _ = f__need({ let t6824 = { let t6817 = py::getitem(&v_letters.clone(), &Py::Int(0i64))?; Py::Bool({ let t6818 = py::getitem(&v_letters.clone(), &Py::Int(1i64))?; if py::ne(&t6817, &t6818) { true } else { false } }) }; if !t6824.truthy() { t6824 } else { { let t6823 = { let t6819 = py::getitem(&v_letters.clone(), &Py::Int(0i64))?; Py::Bool({ let t6820 = py::getitem(&v_letters.clone(), &Py::Int(2i64))?; if py::ne(&t6819, &t6820) { true } else { false } }) }; if !t6823.truthy() { t6823 } else { { let t6821 = py::getitem(&v_letters.clone(), &Py::Int(1i64))?; Py::Bool({ let t6822 = py::getitem(&v_letters.clone(), &Py::Int(2i64))?; if py::ne(&t6821, &t6822) { true } else { false } }) } } } } }, Py::str("the three genes need different letters"))?;
    for t6825 in py::range_iter(&[(Py::Int(8i64)).i()?]) {
        v_i = t6825;
        for t6826 in py::range_iter(&[(py::add(&v_i.clone(), &Py::Int(1i64))?).i()?, (Py::Int(8i64)).i()?]) {
            v_j = t6826;
            let _ = f__need({ let t6827 = py::getitem(&v_names.clone(), &v_i.clone())?; Py::Bool({ let t6828 = py::getitem(&v_names.clone(), &v_j.clone())?; if py::ne(&t6827, &t6828) { true } else { false } }) }, Py::str("each of the 8 classes must be given once"))?;
        }
    }
    v_n = Py::Int(0i64);
    for t6829 in py::iter(&(v_cnt.clone()))? {
        v_c = t6829;
        v_n = py::iadd(&v_n, &v_c.clone())?;
    }
    let _ = f__need({ let t6830 = v_n.clone(); Py::Bool({ let t6831 = Py::Int(1i64); if py::ge(&t6830, &t6831)? { true } else { false } }) }, Py::str("at least one offspring is needed"))?;
    return Ok(Py::tuple(vec![v_names.clone(), v_cnt.clone(), v_letters.clone(), v_n.clone()]));
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn f__tp_rec(mut v_names: Py, mut v_cnt: Py, mut v_P: Py, mut v_i: Py, mut v_j: Py) -> R {
    let mut v_k: Py = Py::None;
    let mut v_s: Py = Py::None;
    v_s = Py::Int(0i64);
    for t6867 in py::range_iter(&[(Py::Int(8i64)).i()?]) {
        v_k = t6867;
        if ({ let t6868 = { let t6870 = py::getitem(&py::getitem(&v_names.clone(), &v_k.clone())?, &v_i.clone())?; Py::Bool({ let t6871 = py::getitem(&v_P.clone(), &v_i.clone())?; if py::eq(&t6870, &t6871) { true } else { false } }) }; Py::Bool({ let t6869 = { let t6872 = py::getitem(&py::getitem(&v_names.clone(), &v_k.clone())?, &v_j.clone())?; Py::Bool({ let t6873 = py::getitem(&v_P.clone(), &v_j.clone())?; if py::eq(&t6872, &t6873) { true } else { false } }) }; if py::ne(&t6868, &t6869) { true } else { false } }) }).truthy() {
            v_s = py::iadd(&v_s, &py::getitem(&v_cnt.clone(), &v_k.clone())?)?;
        }
    }
    return Ok(v_s.clone());
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn f__translate(mut v_s: Py, mut v_table: Py, mut v_truncate: Py) -> R {
    let mut v_p: Py = Py::None;
    v_p = py::method(&Py::str(""), "join", &[{ let mut out = vec![]; for t7200 in py::range_iter(&[(Py::Int(0i64)).i()?, (py::len(&v_s.clone())?).i()?, (Py::Int(3i64)).i()?]) { let mut v_i: Py; 
v_i = t7200;
out.push(py::getitem(&v_table.clone(), &py::slice(&v_s.clone(), Some((v_i.clone()).i()?), Some((py::add(&v_i.clone(), &Py::Int(3i64))?).i()?), None)?)?); } Py::list(out) }])?;
    if ({ let t7203 = v_truncate.clone(); if !t7203.truthy() { t7203 } else { { let t7201 = Py::str("*"); Py::Bool({ let t7202 = v_p.clone(); if py::contains(&t7202, &t7201)? { true } else { false } }) } } }).truthy() {
        v_p = py::slice(&v_p.clone(), None, Some((py::method(&v_p.clone(), "index", &[Py::str("*")])?).i()?), None)?;
    }
    return Ok(v_p.clone());
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn f__whole(mut v_name: Py, mut v_v: Py, mut v_lo: Py) -> R {
    let _ = f__need({ let t6803 = Py::Bool(matches!(v_v.clone(), Py::Int(_) | Py::Big(_) | Py::Bool(_) | Py::Float(_))); if !t6803.truthy() { t6803 } else { { let t6802 = py::not(&Py::Bool(matches!(v_v.clone(), Py::Bool(_)))); if !t6802.truthy() { t6802 } else { { let t6801 = py::math1("isfinite", &v_v.clone())?; if !t6801.truthy() { t6801 } else { { let t6800 = { let t6796 = py::float(&v_v.clone())?; Py::Bool({ let t6797 = py::int(&v_v.clone())?; if py::eq(&t6796, &t6797) { true } else { false } }) }; if !t6800.truthy() { t6800 } else { { let t6798 = v_v.clone(); Py::Bool({ let t6799 = v_lo.clone(); if py::ge(&t6798, &t6799)? { true } else { false } }) } } } } } } } } }, Py::str(&[py::format(&v_name.clone(), "")?, " must be a whole number of at least ".to_string(), py::format(&v_lo.clone(), "")?].concat()))?;
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
fn f_hardy_weinberg_frequencies(mut v_p: Py, mut v_q: Py, mut v_q2: Py, mut v_N: Py) -> R {
    let mut v_given: Py = Py::None;
    let mut v_het: Py = Py::None;
    let mut v_out: Py = Py::None;
    let mut v_pp: Py = Py::None;
    let mut v_qq: Py = Py::None;
    let mut v_v: Py = Py::None;
    v_given = Py::Int(0i64);
    for t2455 in py::iter(&(Py::tuple(vec![v_p.clone(), v_q.clone(), v_q2.clone()])))? {
        v_v = t2455;
        if ({ let t2456 = v_v.clone(); Py::Bool({ let t2457 = Py::None; if !py::is(&t2456, &t2457) { true } else { false } }) }).truthy() {
            v_given = py::iadd(&v_given, &Py::Int(1i64))?;
        }
    }
    let _ = f__need({ let t2458 = v_given.clone(); Py::Bool({ let t2459 = Py::Int(1i64); if py::eq(&t2458, &t2459) { true } else { false } }) }, Py::str("give exactly one of p (frequency of A), q (frequency of a) or q2 (frequency of the recessive phenotype aa)"))?;
    if ({ let t2460 = v_p.clone(); Py::Bool({ let t2461 = Py::None; if !py::is(&t2460, &t2461) { true } else { false } }) }).truthy() {
        let _ = f__frac(Py::str("p"), v_p.clone(), Py::Bool(true), Py::Bool(true))?;
        v_qq = py::sub(&Py::Float(1.0), &v_p.clone())?;
    } else {
        if ({ let t2474 = v_q.clone(); Py::Bool({ let t2475 = Py::None; if !py::is(&t2474, &t2475) { true } else { false } }) }).truthy() {
            let _ = f__frac(Py::str("q"), v_q.clone(), Py::Bool(true), Py::Bool(true))?;
            v_qq = v_q.clone();
        } else {
            let _ = f__frac(Py::str("q2"), v_q2.clone(), Py::Bool(true), Py::Bool(true))?;
            v_qq = py::math1("sqrt", &v_q2.clone())?;
        }
    }
    v_pp = py::sub(&Py::Float(1.0), &v_qq.clone())?;
    v_het = py::mul(&py::mul(&Py::Float(2.0), &v_pp.clone())?, &v_qq.clone())?;
    v_out = Py::dict(vec![(Py::str("p"), v_pp.clone()), (Py::str("q"), v_qq.clone()), (Py::str("p2"), py::mul(&v_pp.clone(), &v_pp.clone())?), (Py::str("het_2pq"), v_het.clone()), (Py::str("q2"), py::mul(&v_qq.clone(), &v_qq.clone())?), (Py::str("carrier_frequency"), v_het.clone()), (Py::str("carrier_among_unaffected"), py::truediv(&v_het.clone(), &py::sub(&Py::Float(1.0), &py::mul(&v_qq.clone(), &v_qq.clone())?)?)?), (Py::str("dominant_phenotype"), py::sub(&Py::Float(1.0), &py::mul(&v_qq.clone(), &v_qq.clone())?)?)]);
    if ({ let t2476 = v_N.clone(); Py::Bool({ let t2477 = Py::None; if !py::is(&t2476, &t2477) { true } else { false } }) }).truthy() {
        let _ = f__need({ let t2480 = Py::Bool(matches!(v_N.clone(), Py::Int(_) | Py::Big(_) | Py::Bool(_))); if !t2480.truthy() { t2480 } else { { let t2478 = v_N.clone(); Py::Bool({ let t2479 = Py::Int(1i64); if py::ge(&t2478, &t2479)? { true } else { false } }) } } }, Py::str("N must be a whole number of at least 1"))?;
        py::setitem(&v_out.clone(), &Py::str("exp_AA"), py::mul(&py::mul(&v_N.clone(), &v_pp.clone())?, &v_pp.clone())?)?;
        py::setitem(&v_out.clone(), &Py::str("exp_Aa"), py::mul(&v_N.clone(), &v_het.clone())?)?;
        py::setitem(&v_out.clone(), &Py::str("exp_aa"), py::mul(&py::mul(&v_N.clone(), &v_qq.clone())?, &v_qq.clone())?)?;
    }
    return Ok(f__out(v_out.clone())?);
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
    for t2537 in py::iter(&(Py::tuple(vec![Py::tuple(vec![Py::str("n_AA"), v_n_AA.clone()]), Py::tuple(vec![Py::str("n_Aa"), v_n_Aa.clone()]), Py::tuple(vec![Py::str("n_aa"), v_n_aa.clone()])])))? {
        let t2541 = (t2537).items()?;
        if t2541.len() != 2 { return Err(PyErr::Value("unpack".into())); }
        v_nm = t2541[0].clone();
        v_v = t2541[1].clone();
        let _ = f__need({ let t2540 = Py::Bool(matches!(v_v.clone(), Py::Int(_) | Py::Big(_) | Py::Bool(_))); if !t2540.truthy() { t2540 } else { { let t2538 = v_v.clone(); Py::Bool({ let t2539 = Py::Int(0i64); if py::ge(&t2538, &t2539)? { true } else { false } }) } } }, Py::str(&[py::format(&v_nm.clone(), "")?, " must be a whole number >= 0".to_string()].concat()))?;
    }
    v_N = py::add(&py::add(&v_n_AA.clone(), &v_n_Aa.clone())?, &v_n_aa.clone())?;
    let _ = f__need({ let t2542 = v_N.clone(); Py::Bool({ let t2543 = Py::Int(1i64); if py::ge(&t2542, &t2543)? { true } else { false } }) }, Py::str("at least one individual is needed"))?;
    v_p = py::fraction(&[py::add(&py::mul(&Py::Int(2i64), &v_n_AA.clone())?, &v_n_Aa.clone())?, py::mul(&Py::Int(2i64), &v_N.clone())?])?;
    let _ = f__need({ let t2544 = Py::Int(0i64); Py::Bool({ let t2545 = v_p.clone(); if py::lt(&t2544, &t2545)? { { let t2546 = Py::Int(1i64); if py::lt(&t2545, &t2546)? { true } else { false } } } else { false } }) }, Py::str("both alleles must be present (otherwise the expected counts are zero)"))?;
    v_q = py::sub(&Py::Int(1i64), &v_p.clone())?;
    v_exp = Py::list(vec![py::mul(&py::mul(&v_N.clone(), &v_p.clone())?, &v_p.clone())?, py::mul(&py::mul(&py::mul(&Py::Int(2i64), &v_N.clone())?, &v_p.clone())?, &v_q.clone())?, py::mul(&py::mul(&v_N.clone(), &v_q.clone())?, &v_q.clone())?]);
    v_obs = Py::list(vec![v_n_AA.clone(), v_n_Aa.clone(), v_n_aa.clone()]);
    v_chi2 = py::sum(&{ let mut out = vec![]; for t2547 in py::iter(&(py::zip(&[v_obs.clone(), v_exp.clone()])?))? { let mut v_e: Py; let mut v_o: Py; 
let t2548 = (t2547).items()?;
if t2548.len() != 2 { return Err(PyErr::Value("unpack".into())); }
v_o = t2548[0].clone();
v_e = t2548[1].clone();
out.push(py::truediv(&py::pow(&py::sub(&v_o.clone(), &v_e.clone())?, &Py::Int(2i64))?, &v_e.clone())?); } Py::list(out) }, None)?;
    let t2549 = (Py::tuple(vec![py::fraction(&[v_n_Aa.clone(), v_N.clone()])?, py::mul(&py::mul(&Py::Int(2i64), &v_p.clone())?, &v_q.clone())?])).items()?;
    if t2549.len() != 2 { return Err(PyErr::Value("unpack".into())); }
    v_ho = t2549[0].clone();
    v_he = t2549[1].clone();
    return Ok(f__out(Py::dict(vec![(Py::str("p"), py::float(&v_p.clone())?), (Py::str("q"), py::float(&v_q.clone())?), (Py::str("exp_AA"), py::float(&py::getitem(&v_exp.clone(), &Py::Int(0i64))?)?), (Py::str("exp_Aa"), py::float(&py::getitem(&v_exp.clone(), &Py::Int(1i64))?)?), (Py::str("exp_aa"), py::float(&py::getitem(&v_exp.clone(), &Py::Int(2i64))?)?), (Py::str("chi2"), py::float(&v_chi2.clone())?), (Py::str("p_value"), py::math1("erfc", &py::math1("sqrt", &py::truediv(&py::float(&v_chi2.clone())?, &Py::Int(2i64))?)?)?), (Py::str("Ho"), py::float(&v_ho.clone())?), (Py::str("He"), py::float(&v_he.clone())?), (Py::str("F"), py::float(&py::sub(&Py::Int(1i64), &py::truediv(&v_ho.clone(), &v_he.clone())?)?)?)]))?);
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn f_het_advantage(mut v_w_AA: Py, mut v_w_Aa: Py, mut v_w_aa: Py) -> R {
    let _ = f__pos(Py::str("fitness"), Py::tuple(vec![v_w_AA.clone(), v_w_Aa.clone(), v_w_aa.clone()]))?;
    let _ = f__need({ let t2615 = { let t2611 = v_w_Aa.clone(); Py::Bool({ let t2612 = v_w_AA.clone(); if py::gt(&t2611, &t2612)? { true } else { false } }) }; if !t2615.truthy() { t2615 } else { { let t2613 = v_w_Aa.clone(); Py::Bool({ let t2614 = v_w_aa.clone(); if py::gt(&t2613, &t2614)? { true } else { false } }) } } }, Py::str("heterozygote advantage needs w_Aa greater than both homozygote fitnesses"))?;
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
    let t2687 = (f__kin_data(v_S.clone(), v_v.clone(), Py::Int(4i64))?).items()?;
    if t2687.len() != 2 { return Err(PyErr::Value("unpack".into())); }
    v_S = t2687[0].clone();
    v_v = t2687[1].clone();
    let _ = f__need({ let t2688 = py::len(&py::set(&v_S.clone())?)?; Py::Bool({ let t2689 = Py::Int(4i64); if py::ge(&t2688, &t2689)? { true } else { false } }) }, Py::str("at least 4 different concentrations are needed"))?;
    if (v_fixed_n.clone()).truthy() {
        v_r = f_mm_fit(v_S.clone(), v_v.clone())?;
        return Ok(Py::dict(vec![(Py::str("Vmax"), py::getitem(&v_r.clone(), &Py::str("Vmax"))?), (Py::str("K"), py::getitem(&v_r.clone(), &Py::str("Km"))?), (Py::str("n"), Py::Float(1.0))]));
    }
    v_vmax0 = py::mul(&py::minmax(&[v_v.clone()], true)?, &Py::Float(1.1))?;
    v_Y = { let mut out = vec![]; for t2750 in py::iter(&(v_v.clone()))? { let mut v_y: Py; 
v_y = t2750;
out.push(py::minmax(&[py::minmax(&[py::truediv(&v_y.clone(), &v_vmax0.clone())?, Py::Float(0.001)], true)?, py::sub(&Py::Int(1i64), &Py::Float(0.001))?], false)?); } Py::list(out) };
    let t2762 = (f__ols({ let mut out = vec![]; for t2760 in py::iter(&(v_S.clone()))? { let mut v_s: Py; 
v_s = t2760;
out.push(py::math1("log", &v_s.clone())?); } Py::list(out) }, { let mut out = vec![]; for t2761 in py::iter(&(v_Y.clone()))? { let mut v_y: Py; 
v_y = t2761;
out.push(py::math1("log", &py::truediv(&v_y.clone(), &py::sub(&Py::Int(1i64), &v_y.clone())?)?)?); } Py::list(out) })?).items()?;
    if t2762.len() != 3 { return Err(PyErr::Value("unpack".into())); }
    v_b = t2762[0].clone();
    v_a = t2762[1].clone();
    v__ = t2762[2].clone();
    v_p0 = Py::list(vec![v_vmax0.clone(), (if ({ let t2763 = v_b.clone(); Py::Bool({ let t2764 = Py::Int(0i64); if py::gt(&t2763, &t2764)? { true } else { false } }) }).truthy() { py::math1("exp", &py::truediv(&py::neg(&v_a.clone())?, &v_b.clone())?)? } else { py::getitem(&py::sorted(&v_S.clone(), (Py::Bool(false)).truthy())?, &py::floordiv(&py::len(&v_S.clone())?, &Py::Int(2i64))?)? }), py::minmax(&[v_b.clone(), Py::Float(0.5)], true)?]);
    v_p = f__lm(f__hill_model(v_S.clone())?, v_p0.clone(), v_v.clone(), { let v_S = v_S.clone(); let v_Y = v_Y.clone(); let v__ = v__.clone(); let v_a = v_a.clone(); let v_b = v_b.clone(); let v_fixed_n = v_fixed_n.clone(); let v_p0 = v_p0.clone(); let v_r = v_r.clone(); let v_v = v_v.clone(); let v_vmax0 = v_vmax0.clone(); py::func(move |a: &[Py]| -> R {
            #[allow(unused_mut, unused_variables, unused_assignments)]
{
            let mut v_q: Py = a.get(0).cloned().unwrap_or(Py::None);
            return Ok({ let t2834 = { let t2826 = py::getitem(&v_q.clone(), &Py::Int(0i64))?; Py::Bool({ let t2827 = Py::Int(0i64); if py::gt(&t2826, &t2827)? { true } else { false } }) }; if !t2834.truthy() { t2834 } else { { let t2833 = { let t2828 = py::getitem(&v_q.clone(), &Py::Int(1i64))?; Py::Bool({ let t2829 = Py::Int(0i64); if py::gt(&t2828, &t2829)? { true } else { false } }) }; if !t2833.truthy() { t2833 } else { { let t2830 = Py::Float(0.05); Py::Bool({ let t2831 = py::getitem(&v_q.clone(), &Py::Int(2i64))?; if py::lt(&t2830, &t2831)? { { let t2832 = Py::Int(20i64); if py::lt(&t2831, &t2832)? { true } else { false } } } else { false } }) } } } } });
            #[allow(unreachable_code)]
            Ok(Py::None)
        }
        }) })?;
    let t2835 = (py::call(&f__hill_model(v_S.clone())?, &[v_p.clone()])?).items()?;
    if t2835.len() != 2 { return Err(PyErr::Value("unpack".into())); }
    v_f = t2835[0].clone();
    v_J = t2835[1].clone();
    v_r = { let mut out = vec![]; for t2836 in py::iter(&(py::zip(&[v_v.clone(), v_f.clone()])?))? { let mut v_a: Py; let mut v_y: Py; 
let t2837 = (t2836).items()?;
if t2837.len() != 2 { return Err(PyErr::Value("unpack".into())); }
v_y = t2837[0].clone();
v_a = t2837[1].clone();
out.push(py::sub(&v_y.clone(), &v_a.clone())?); } Py::list(out) };
    v_g = { let mut out = vec![]; for t2838 in py::range_iter(&[(Py::Int(3i64)).i()?]) { let mut v_i: Py; 
v_i = t2838;
out.push(py::abs(&py::fsum(&{ let mut out = vec![]; for t2839 in py::iter(&(py::zip(&[v_J.clone(), v_r.clone()])?))? { let mut v_ri: Py; let mut v_row: Py; 
let t2840 = (t2839).items()?;
if t2840.len() != 2 { return Err(PyErr::Value("unpack".into())); }
v_row = t2840[0].clone();
v_ri = t2840[1].clone();
out.push(py::mul(&py::getitem(&v_row.clone(), &v_i.clone())?, &v_ri.clone())?); } Py::list(out) })?)?); } Py::list(out) };
    v_scale = py::fsum(&{ let mut out = vec![]; for t2841 in py::iter(&(v_v.clone()))? { let mut v_y: Py; 
v_y = t2841;
out.push(py::mul(&v_y.clone(), &v_y.clone())?); } Py::list(out) })?;
    let _ = f__need({ let t2842 = py::minmax(&[{ let mut out = vec![]; for t2844 in py::iter(&(py::zip(&[v_g.clone(), v_p.clone()])?))? { let mut v_gi: Py; let mut v_pi: Py; 
let t2845 = (t2844).items()?;
if t2845.len() != 2 { return Err(PyErr::Value("unpack".into())); }
v_gi = t2845[0].clone();
v_pi = t2845[1].clone();
out.push(py::mul(&v_gi.clone(), &py::abs(&v_pi.clone())?)?); } Py::list(out) }], true)?; Py::Bool({ let t2843 = py::mul(&Py::Float(1e-09), &v_scale.clone())?; if py::le(&t2842, &t2843)? { true } else { false } }) }, Py::str("the Hill fit did not converge"))?;
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
    let _ = f__need({ let t2893 = Py::Bool(matches!(v_S.clone(), Py::List(_))); if !t2893.truthy() { t2893 } else { { let t2892 = Py::Bool(matches!(v_Y.clone(), Py::List(_))); if !t2892.truthy() { t2892 } else { { let t2891 = { let t2887 = py::len(&v_S.clone())?; Py::Bool({ let t2888 = py::len(&v_Y.clone())?; if py::eq(&t2887, &t2888) { true } else { false } }) }; if !t2891.truthy() { t2891 } else { { let t2889 = py::len(&v_S.clone())?; Py::Bool({ let t2890 = Py::Int(2i64); if py::ge(&t2889, &t2890)? { true } else { false } }) } } } } } } }, Py::str("S and Y must be lists of equal length (at least 2 points)"))?;
    for t2894 in py::iter(&(v_S.clone()))? {
        v_s = t2894;
        let _ = f__pos(Py::str("S"), Py::tuple(vec![v_s.clone()]))?;
    }
    for t2901 in py::iter(&(v_Y.clone()))? {
        v_y = t2901;
        let _ = f__frac(Py::str("Y"), v_y.clone(), Py::Bool(true), Py::Bool(true))?;
    }
    let t2925 = (f__ols({ let mut out = vec![]; for t2923 in py::iter(&(v_S.clone()))? { let mut v_s: Py; 
v_s = t2923;
out.push(py::math1("log10", &v_s.clone())?); } Py::list(out) }, { let mut out = vec![]; for t2924 in py::iter(&(v_Y.clone()))? { let mut v_y: Py; 
v_y = t2924;
out.push(py::math1("log10", &py::truediv(&v_y.clone(), &py::sub(&Py::Int(1i64), &v_y.clone())?)?)?); } Py::list(out) })?).items()?;
    if t2925.len() != 3 { return Err(PyErr::Value("unpack".into())); }
    v_b = t2925[0].clone();
    v_a = t2925[1].clone();
    v_r2 = t2925[2].clone();
    let _ = f__need({ let t2926 = v_b.clone(); Py::Bool({ let t2927 = Py::Int(0i64); if py::gt(&t2926, &t2927)? { true } else { false } }) }, Py::str("the Hill plot must rise (positive Hill coefficient)"))?;
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
    let _ = f__need({ let t3006 = py::len(&v_s.clone())?; Py::Bool({ let t3007 = Py::Int(2i64); if py::ge(&t3006, &t3007)? { true } else { false } }) }, Py::str("the protein needs at least 2 residues"))?;
    v_score = Py::Float(0.0);
    for t3008 in py::range_iter(&[(py::sub(&py::len(&v_s.clone())?, &Py::Int(1i64))?).i()?]) {
        v_i = t3008;
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
    let t3205 = (f__aligned_dna(v_seq1.clone(), v_seq2.clone())?).items()?;
    if t3205.len() != 2 { return Err(PyErr::Value("unpack".into())); }
    v_a = t3205[0].clone();
    v_b = t3205[1].clone();
    v_p = py::truediv(&py::sum(&{ let mut out = vec![]; for t3206 in py::iter(&(py::zip(&[v_a.clone(), v_b.clone()])?))? { let mut v_x: Py; let mut v_y: Py; 
let t3207 = (t3206).items()?;
if t3207.len() != 2 { return Err(PyErr::Value("unpack".into())); }
v_x = t3207[0].clone();
v_y = t3207[1].clone();
if ({ let t3208 = v_x.clone(); Py::Bool({ let t3209 = v_y.clone(); if py::ne(&t3208, &t3209) { true } else { false } }) }).truthy() { out.push(Py::Int(1i64)); }} Py::list(out) }, None)?, &py::len(&v_a.clone())?)?;
    let _ = f__need({ let t3210 = v_p.clone(); Py::Bool({ let t3211 = Py::Float(0.75); if py::lt(&t3210, &t3211)? { true } else { false } }) }, Py::str("p >= 0.75: the sequences are saturated and the Jukes-Cantor distance is undefined"))?;
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
    let t3293 = (f__aligned_dna(v_seq1.clone(), v_seq2.clone())?).items()?;
    if t3293.len() != 2 { return Err(PyErr::Value("unpack".into())); }
    v_a = t3293[0].clone();
    v_b = t3293[1].clone();
    v_n = py::len(&v_a.clone())?;
    v_ts = py::sum(&{ let mut out = vec![]; for t3294 in py::iter(&(py::zip(&[v_a.clone(), v_b.clone()])?))? { let mut v_x: Py; let mut v_y: Py; 
let t3295 = (t3294).items()?;
if t3295.len() != 2 { return Err(PyErr::Value("unpack".into())); }
v_x = t3295[0].clone();
v_y = t3295[1].clone();
if ({ let t3304 = { let t3296 = v_x.clone(); Py::Bool({ let t3297 = v_y.clone(); if py::ne(&t3296, &t3297) { true } else { false } }) }; if !t3304.truthy() { t3304 } else { { let t3298 = { let t3300 = v_x.clone(); Py::Bool({ let t3301 = c__PUR()?; if py::contains(&t3301, &t3300)? { true } else { false } }) }; Py::Bool({ let t3299 = { let t3302 = v_y.clone(); Py::Bool({ let t3303 = c__PUR()?; if py::contains(&t3303, &t3302)? { true } else { false } }) }; if py::eq(&t3298, &t3299) { true } else { false } }) } } }).truthy() { out.push(Py::Int(1i64)); }} Py::list(out) }, None)?;
    v_tv = py::sum(&{ let mut out = vec![]; for t3305 in py::iter(&(py::zip(&[v_a.clone(), v_b.clone()])?))? { let mut v_x: Py; let mut v_y: Py; 
let t3306 = (t3305).items()?;
if t3306.len() != 2 { return Err(PyErr::Value("unpack".into())); }
v_x = t3306[0].clone();
v_y = t3306[1].clone();
if ({ let t3307 = { let t3309 = v_x.clone(); Py::Bool({ let t3310 = c__PUR()?; if py::contains(&t3310, &t3309)? { true } else { false } }) }; Py::Bool({ let t3308 = { let t3311 = v_y.clone(); Py::Bool({ let t3312 = c__PUR()?; if py::contains(&t3312, &t3311)? { true } else { false } }) }; if py::ne(&t3307, &t3308) { true } else { false } }) }).truthy() { out.push(Py::Int(1i64)); }} Py::list(out) }, None)?;
    let t3313 = (Py::tuple(vec![py::truediv(&v_ts.clone(), &v_n.clone())?, py::truediv(&v_tv.clone(), &v_n.clone())?])).items()?;
    if t3313.len() != 2 { return Err(PyErr::Value("unpack".into())); }
    v_P = t3313[0].clone();
    v_Q = t3313[1].clone();
    let _ = f__need({ let t3318 = { let t3314 = py::sub(&py::sub(&Py::Int(1i64), &py::mul(&Py::Int(2i64), &v_P.clone())?)?, &v_Q.clone())?; Py::Bool({ let t3315 = Py::Int(0i64); if py::gt(&t3314, &t3315)? { true } else { false } }) }; if !t3318.truthy() { t3318 } else { { let t3316 = py::sub(&Py::Int(1i64), &py::mul(&Py::Int(2i64), &v_Q.clone())?)?; Py::Bool({ let t3317 = Py::Int(0i64); if py::gt(&t3316, &t3317)? { true } else { false } }) } } }, Py::str("too many differences: the Kimura distance is undefined"))?;
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
    let _ = f__need({ let t3400 = Py::Bool(matches!(v_k.clone(), Py::Int(_) | Py::Big(_) | Py::Bool(_))); if !t3400.truthy() { t3400 } else { { let t3397 = Py::Int(1i64); Py::Bool({ let t3398 = v_k.clone(); if py::le(&t3397, &t3398)? { { let t3399 = py::len(&v_s.clone())?; if py::le(&t3398, &t3399)? { true } else { false } } } else { false } }) } } }, Py::str("k must be a whole number between 1 and the sequence length"))?;
    v_counts = Py::dict(vec![]);
    for t3401 in py::range_iter(&[(py::add(&py::sub(&py::len(&v_s.clone())?, &v_k.clone())?, &Py::Int(1i64))?).i()?]) {
        v_i = t3401;
        v_w = py::slice(&v_s.clone(), Some((v_i.clone()).i()?), Some((py::add(&v_i.clone(), &v_k.clone())?).i()?), None)?;
        py::setitem(&v_counts.clone(), &v_w.clone(), py::add(&py::method(&v_counts.clone(), "get", &[v_w.clone(), Py::Int(0i64)])?, &Py::Int(1i64))?)?;
    }
    v_top = py::minmax(&[py::method(&v_counts.clone(), "values", &[])?], true)?;
    return Ok(Py::dict(vec![(Py::str("total"), py::add(&py::sub(&py::len(&v_s.clone())?, &v_k.clone())?, &Py::Int(1i64))?), (Py::str("distinct"), py::len(&v_counts.clone())?), (Py::str("max_count"), v_top.clone()), (Py::str("most_frequent"), py::minmax(&[{ let mut out = vec![]; for t3402 in py::iter(&(py::method(&v_counts.clone(), "items", &[])?))? { let mut v_c: Py; let mut v_w: Py; 
let t3403 = (t3402).items()?;
if t3403.len() != 2 { return Err(PyErr::Value("unpack".into())); }
v_w = t3403[0].clone();
v_c = t3403[1].clone();
if ({ let t3404 = v_c.clone(); Py::Bool({ let t3405 = v_top.clone(); if py::eq(&t3404, &t3405) { true } else { false } }) }).truthy() { out.push(v_w.clone()); }} Py::list(out) }], false)?)]));
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
    v_R0 = py::fsum(&{ let mut out = vec![]; for t3468 in py::iter(&(py::zip(&[v_F.clone(), v_l.clone()])?))? { let mut v_f: Py; let mut v_li: Py; 
let t3469 = (t3468).items()?;
if t3469.len() != 2 { return Err(PyErr::Value("unpack".into())); }
v_f = t3469[0].clone();
v_li = t3469[1].clone();
out.push(py::mul(&v_f.clone(), &v_li.clone())?); } Py::list(out) })?;
    if (v_use_r0.clone()).truthy() {
        v_lam = v_R0.clone();
    } else {
        v_g = { let v_F = v_F.clone(); let v_R0 = v_R0.clone(); let v_S = v_S.clone(); let v_l = v_l.clone(); let v_lam = v_lam.clone(); let v_use_r0 = v_use_r0.clone(); py::func(move |a: &[Py]| -> R {
            #[allow(unused_mut, unused_variables, unused_assignments)]
{
            let mut v_x: Py = a.get(0).cloned().unwrap_or(Py::None);
            return Ok(py::sub(&py::fsum(&{ let mut out = vec![]; for t3470 in py::iter(&(py::enumerate(&py::zip(&[v_F.clone(), v_l.clone()])?, 0)?))? { let mut v_f: Py; let mut v_i: Py; let mut v_li: Py; 
let t3471 = (t3470).items()?;
if t3471.len() != 2 { return Err(PyErr::Value("unpack".into())); }
v_i = t3471[0].clone();
let t3472 = (t3471[1].clone()).items()?;
if t3472.len() != 2 { return Err(PyErr::Value("unpack".into())); }
v_f = t3472[0].clone();
v_li = t3472[1].clone();
out.push(py::mul(&py::mul(&v_f.clone(), &v_li.clone())?, &py::math1("exp", &py::mul(&py::neg(&py::add(&v_i.clone(), &Py::Int(1i64))?)?, &v_x.clone())?)?)?); } Py::list(out) })?, &Py::Float(1.0))?);
            #[allow(unreachable_code)]
            Ok(Py::None)
        }
        }) };
        v_hi = py::add(&py::math1("log", &py::minmax(&[v_R0.clone(), Py::Float(1.0)], true)?)?, &Py::Float(1.0))?;
        v_lo = py::neg(&Py::Float(1.0))?;
        while ({ let t3473 = py::call(&v_g.clone(), &[v_lo.clone()])?; Py::Bool({ let t3474 = Py::Int(0i64); if py::le(&t3473, &t3474)? { true } else { false } }) }).truthy() {
            v_lo = py::sub(&v_lo, &Py::Float(1.0))?;
        }
        while ({ let t3475 = py::call(&v_g.clone(), &[v_hi.clone()])?; Py::Bool({ let t3476 = Py::Int(0i64); if py::ge(&t3475, &t3476)? { true } else { false } }) }).truthy() {
            v_hi = py::iadd(&v_hi, &Py::Float(1.0))?;
        }
        v_lam = py::math1("exp", &f__bisect(v_g.clone(), v_lo.clone(), v_hi.clone(), Py::Int(400i64))?)?;
    }
    v_w = { let mut out = vec![]; for t3503 in py::iter(&(py::enumerate(&v_l.clone(), 0)?))? { let mut v_i: Py; let mut v_li: Py; 
let t3504 = (t3503).items()?;
if t3504.len() != 2 { return Err(PyErr::Value("unpack".into())); }
v_i = t3504[0].clone();
v_li = t3504[1].clone();
out.push(py::mul(&v_li.clone(), &py::pow(&v_lam.clone(), &py::neg(&v_i.clone())?)?)?); } Py::list(out) };
    v_tot = py::fsum(&v_w.clone())?;
    return Ok(f__out(Py::dict(vec![(Py::str("lambda1"), v_lam.clone()), (Py::str("r"), py::math1("log", &v_lam.clone())?), (Py::str("R0"), v_R0.clone()), (Py::str("stable"), { let mut out = vec![]; for t3519 in py::iter(&(v_w.clone()))? { let mut v_x: Py; 
v_x = t3519;
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
    let _ = f__need({ let t3583 = Py::Bool(matches!(v_n0.clone(), Py::List(_))); if !t3583.truthy() { t3583 } else { { let t3581 = py::len(&v_n0.clone())?; Py::Bool({ let t3582 = py::len(&v_F.clone())?; if py::eq(&t3581, &t3582) { true } else { false } }) } } }, Py::str("n0 needs one number per age class"))?;
    for t3584 in py::iter(&(v_n0.clone()))? {
        v_x = t3584;
        let _ = f__fin(Py::str("n0"), Py::tuple(vec![v_x.clone()]))?;
        let _ = f__need({ let t3585 = v_x.clone(); Py::Bool({ let t3586 = Py::Int(0i64); if py::ge(&t3585, &t3586)? { true } else { false } }) }, Py::str("population numbers must be zero or positive"))?;
    }
    let _ = f__need({ let t3590 = Py::Bool(matches!(v_t.clone(), Py::Int(_) | Py::Big(_) | Py::Bool(_))); if !t3590.truthy() { t3590 } else { { let t3587 = Py::Int(0i64); Py::Bool({ let t3588 = v_t.clone(); if py::le(&t3587, &t3588)? { { let t3589 = Py::Int(500i64); if py::le(&t3588, &t3589)? { true } else { false } } } else { false } }) } } }, Py::str("t must be a whole number of time steps (0 to 500)"))?;
    v_k = py::len(&v_F.clone())?;
    v_L = { let mut out = vec![]; for t3591 in py::range_iter(&[(v_k.clone()).i()?]) { let mut v__: Py; 
v__ = t3591;
out.push(py::mul(&Py::list(vec![Py::Float(0.0)]), &v_k.clone())?); } Py::list(out) };
    py::setitem(&v_L.clone(), &Py::Int(0i64), { let mut out = vec![]; for t3592 in py::iter(&(v_F.clone()))? { let mut v_f: Py; 
v_f = t3592;
out.push(py::float(&v_f.clone())?); } Py::list(out) })?;
    for t3593 in py::iter(&(py::enumerate(&v_S.clone(), 0)?))? {
        let t3594 = (t3593).items()?;
        if t3594.len() != 2 { return Err(PyErr::Value("unpack".into())); }
        v_i = t3594[0].clone();
        v_s = t3594[1].clone();
        py::setitem(&py::getitem(&v_L.clone(), &py::add(&v_i.clone(), &Py::Int(1i64))?)?, &v_i.clone(), py::float(&v_s.clone())?)?;
    }
    if (v_transpose.clone()).truthy() {
        v_L = { let mut out = vec![]; for t3595 in py::iter(&(py::zip(&{ let mut v: Vec<Py> = vec![]; v.extend((v_L.clone()).items()?); v })?))? { let mut v_r: Py; 
v_r = t3595;
out.push(Py::list((v_r.clone()).items()?)); } Py::list(out) };
    }
    v_n = { let mut out = vec![]; for t3596 in py::iter(&(v_n0.clone()))? { let mut v_x: Py; 
v_x = t3596;
out.push(py::float(&v_x.clone())?); } Py::list(out) };
    for t3597 in py::range_iter(&[(v_t.clone()).i()?]) {
        v__ = t3597;
        v_n = { let mut out = vec![]; for t3598 in py::range_iter(&[(v_k.clone()).i()?]) { let mut v_i: Py; 
v_i = t3598;
out.push(py::fsum(&{ let mut out = vec![]; for t3599 in py::range_iter(&[(v_k.clone()).i()?]) { let mut v_j: Py; 
v_j = t3599;
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
    let t3671 = (f__kin_data(v_S.clone(), v_v.clone(), Py::Int(2i64))?).items()?;
    if t3671.len() != 2 { return Err(PyErr::Value("unpack".into())); }
    v_S = t3671[0].clone();
    v_v = t3671[1].clone();
    let t3683 = (f__ols({ let mut out = vec![]; for t3681 in py::iter(&(v_S.clone()))? { let mut v_x: Py; 
v_x = t3681;
out.push(py::truediv(&Py::Int(1i64), &v_x.clone())?); } Py::list(out) }, { let mut out = vec![]; for t3682 in py::iter(&(v_v.clone()))? { let mut v_y: Py; 
v_y = t3682;
out.push(py::truediv(&Py::Int(1i64), &v_y.clone())?); } Py::list(out) })?).items()?;
    if t3683.len() != 3 { return Err(PyErr::Value("unpack".into())); }
    v_b = t3683[0].clone();
    v_a = t3683[1].clone();
    v_r2 = t3683[2].clone();
    let _ = f__need({ let t3688 = { let t3684 = v_a.clone(); Py::Bool({ let t3685 = Py::Int(0i64); if py::gt(&t3684, &t3685)? { true } else { false } }) }; if !t3688.truthy() { t3688 } else { { let t3686 = v_b.clone(); Py::Bool({ let t3687 = Py::Int(0i64); if py::gt(&t3686, &t3687)? { true } else { false } }) } } }, Py::str("the double-reciprocal line must have a positive slope and intercept"))?;
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
    for t3744 in py::iter(&(v_c.clone()))? {
        v_v = t3744;
        let _ = f__need({ let t3747 = Py::Bool(matches!(v_v.clone(), Py::Int(_) | Py::Big(_) | Py::Bool(_))); if !t3747.truthy() { t3747 } else { { let t3745 = v_v.clone(); Py::Bool({ let t3746 = Py::Int(0i64); if py::ge(&t3745, &t3746)? { true } else { false } }) } } }, Py::str("haplotype counts must be whole numbers >= 0"))?;
    }
    v_n = py::sum(&v_c.clone(), None)?;
    let _ = f__need({ let t3748 = v_n.clone(); Py::Bool({ let t3749 = Py::Int(1i64); if py::ge(&t3748, &t3749)? { true } else { false } }) }, Py::str("at least one haplotype is needed"))?;
    let t3750 = (Py::tuple(vec![py::fraction(&[py::add(&v_n_AB.clone(), &v_n_Ab.clone())?, v_n.clone()])?, py::fraction(&[py::add(&v_n_AB.clone(), &v_n_aB.clone())?, v_n.clone()])?])).items()?;
    if t3750.len() != 2 { return Err(PyErr::Value("unpack".into())); }
    v_pA = t3750[0].clone();
    v_pB = t3750[1].clone();
    let _ = f__need({ let t3757 = { let t3751 = Py::Int(0i64); Py::Bool({ let t3752 = v_pA.clone(); if py::lt(&t3751, &t3752)? { { let t3753 = Py::Int(1i64); if py::lt(&t3752, &t3753)? { true } else { false } } } else { false } }) }; if !t3757.truthy() { t3757 } else { { let t3754 = Py::Int(0i64); Py::Bool({ let t3755 = v_pB.clone(); if py::lt(&t3754, &t3755)? { { let t3756 = Py::Int(1i64); if py::lt(&t3755, &t3756)? { true } else { false } } } else { false } }) } } }, Py::str("both loci must be polymorphic"))?;
    v_D = py::sub(&py::fraction(&[v_n_AB.clone(), v_n.clone()])?, &py::mul(&v_pA.clone(), &v_pB.clone())?)?;
    if ({ let t3758 = v_D.clone(); Py::Bool({ let t3759 = Py::Int(0i64); if py::gt(&t3758, &t3759)? { true } else { false } }) }).truthy() {
        v_dmax = py::minmax(&[py::mul(&v_pA.clone(), &py::sub(&Py::Int(1i64), &v_pB.clone())?)?, py::mul(&py::sub(&Py::Int(1i64), &v_pA.clone())?, &v_pB.clone())?], false)?;
    } else {
        v_dmax = py::minmax(&[py::mul(&v_pA.clone(), &v_pB.clone())?, py::mul(&py::sub(&Py::Int(1i64), &v_pA.clone())?, &py::sub(&Py::Int(1i64), &v_pB.clone())?)?], false)?;
    }
    return Ok(Py::dict(vec![(Py::str("D"), py::float(&v_D.clone())?), (Py::str("D_prime"), (if ({ let t3760 = v_D.clone(); Py::Bool({ let t3761 = Py::Int(0i64); if py::ne(&t3760, &t3761) { true } else { false } }) }).truthy() { py::float(&py::truediv(&v_D.clone(), &v_dmax.clone())?)? } else { Py::Float(0.0) })), (Py::str("r2"), py::float(&py::truediv(&py::mul(&v_D.clone(), &v_D.clone())?, &py::mul(&py::mul(&py::mul(&v_pA.clone(), &py::sub(&Py::Int(1i64), &v_pA.clone())?)?, &v_pB.clone())?, &py::sub(&Py::Int(1i64), &v_pB.clone())?)?)?)?)]));
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn f_local_affine(mut v_seq1: Py, mut v_seq2: Py, mut v_match: Py, mut v_mismatch: Py, mut v_open: Py, mut v_extend: Py) -> R {
    let mut v_a: Py = Py::None;
    let mut v_b: Py = Py::None;
    let t3826 = (Py::tuple(vec![f__seq(v_seq1.clone(), Py::str("seq1"), Py::str("ANY"))?, f__seq(v_seq2.clone(), Py::str("seq2"), Py::str("ANY"))?])).items()?;
    if t3826.len() != 2 { return Err(PyErr::Value("unpack".into())); }
    v_a = t3826[0].clone();
    v_b = t3826[1].clone();
    let _ = f__check_align(v_a.clone(), v_b.clone(), Py::dict(vec![(Py::str("match"), v_match.clone()), (Py::str("mismatch"), v_mismatch.clone()), (Py::str("open"), v_open.clone()), (Py::str("extend"), v_extend.clone())]))?;
    let _ = f__need({ let t3838 = { let t3834 = v_open.clone(); Py::Bool({ let t3835 = Py::Int(0i64); if py::le(&t3834, &t3835)? { true } else { false } }) }; if !t3838.truthy() { t3838 } else { { let t3836 = v_extend.clone(); Py::Bool({ let t3837 = Py::Int(0i64); if py::le(&t3836, &t3837)? { true } else { false } }) } } }, Py::str("gap scores must be zero or negative"))?;
    let _ = f__need({ let t3839 = v_match.clone(); Py::Bool({ let t3840 = Py::Int(0i64); if py::gt(&t3839, &t3840)? { true } else { false } }) }, Py::str("the match score must be positive for a local alignment"))?;
    return Ok(f__out(Py::dict(vec![(Py::str("score"), f__align(v_a.clone(), v_b.clone(), f__mm(v_match.clone(), v_mismatch.clone())?, v_open.clone(), v_extend.clone(), Py::str("local"), Py::Bool(false))?)]))?);
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn f_local_linear(mut v_seq1: Py, mut v_seq2: Py, mut v_match: Py, mut v_mismatch: Py, mut v_gap: Py) -> R {
    let mut v_a: Py = Py::None;
    let mut v_b: Py = Py::None;
    let t3950 = (Py::tuple(vec![f__seq(v_seq1.clone(), Py::str("seq1"), Py::str("ANY"))?, f__seq(v_seq2.clone(), Py::str("seq2"), Py::str("ANY"))?])).items()?;
    if t3950.len() != 2 { return Err(PyErr::Value("unpack".into())); }
    v_a = t3950[0].clone();
    v_b = t3950[1].clone();
    let _ = f__check_align(v_a.clone(), v_b.clone(), Py::dict(vec![(Py::str("match"), v_match.clone()), (Py::str("mismatch"), v_mismatch.clone()), (Py::str("gap"), v_gap.clone())]))?;
    let _ = f__need({ let t3958 = v_gap.clone(); Py::Bool({ let t3959 = Py::Int(0i64); if py::le(&t3958, &t3959)? { true } else { false } }) }, Py::str("the gap score must be zero or negative"))?;
    let _ = f__need({ let t3960 = v_match.clone(); Py::Bool({ let t3961 = Py::Int(0i64); if py::gt(&t3960, &t3961)? { true } else { false } }) }, Py::str("the match score must be positive for a local alignment"))?;
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
    let _ = f__need({ let t4123 = Py::Int(0i64); Py::Bool({ let t4124 = v_r.clone(); if py::le(&t4123, &t4124)? { { let t4125 = Py::Float(0.5); if py::lt(&t4124, &t4125)? { true } else { false } } } else { false } }) }, Py::str("the recombination fraction must lie in [0, 0.5)"))?;
    let _ = f__need({ let t4126 = v_map.clone(); Py::Bool({ let t4127 = Py::tuple(vec![Py::str("haldane"), Py::str("kosambi")]); if py::contains(&t4127, &t4126)? { true } else { false } }) }, Py::str("map must be haldane or kosambi"))?;
    v_m = (if ({ let t4128 = v_map.clone(); Py::Bool({ let t4129 = Py::str("haldane"); if py::eq(&t4128, &t4129) { true } else { false } }) }).truthy() { py::mul(&py::neg(&Py::Float(0.5))?, &py::math1("log1p", &py::mul(&py::neg(&Py::Int(2i64))?, &v_r.clone())?)?)? } else { py::mul(&Py::Float(0.25), &py::math1("log", &py::truediv(&py::add(&Py::Int(1i64), &py::mul(&Py::Int(2i64), &v_r.clone())?)?, &py::sub(&Py::Int(1i64), &py::mul(&Py::Int(2i64), &v_r.clone())?)?)?)?)? });
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
    let t4201 = (f__kin_data(v_S.clone(), v_v.clone(), Py::Int(3i64))?).items()?;
    if t4201.len() != 2 { return Err(PyErr::Value("unpack".into())); }
    v_S = t4201[0].clone();
    v_v = t4201[1].clone();
    let _ = f__need({ let t4202 = py::len(&py::set(&v_S.clone())?)?; Py::Bool({ let t4203 = Py::Int(3i64); if py::ge(&t4202, &t4203)? { true } else { false } }) }, Py::str("at least 3 different substrate concentrations are needed"))?;
    let t4204 = (Py::tuple(vec![py::mul(&py::minmax(&[v_S.clone()], false)?, &Py::Float(0.001))?, py::mul(&py::minmax(&[v_S.clone()], true)?, &Py::Float(1000.0))?])).items()?;
    if t4204.len() != 2 { return Err(PyErr::Value("unpack".into())); }
    v_lo = t4204[0].clone();
    v_hi = t4204[1].clone();
    v_grid = { let mut out = vec![]; for t4205 in py::range_iter(&[(Py::Int(301i64)).i()?]) { let mut v_i: Py; 
v_i = t4205;
out.push(py::mul(&v_lo.clone(), &py::pow(&py::truediv(&v_hi.clone(), &v_lo.clone())?, &py::truediv(&v_i.clone(), &Py::Int(300i64))?)?)?); } Py::list(out) };
    v_F = Py::list(vec![]);
    for t4206 in py::iter(&(v_grid.clone()))? {
        v_K = t4206;
        let t4216 = (f__mm_varpro(v_S.clone(), v_v.clone(), v_K.clone())?).items()?;
        if t4216.len() != 3 { return Err(PyErr::Value("unpack".into())); }
        v_svx = t4216[0].clone();
        v_sxx = t4216[1].clone();
        v__ = t4216[2].clone();
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
    let _ = f__need({ let t4217 = Py::Int(0i64); Py::Bool({ let t4218 = v_i.clone(); if py::lt(&t4217, &t4218)? { { let t4219 = py::sub(&py::len(&v_F.clone())?, &Py::Int(1i64))?; if py::lt(&t4218, &t4219)? { true } else { false } } } else { false } }) }, Py::str("the data do not determine Km (the best fit lies at the edge of the search range)"))?;
    v_K = f__bisect({ let v_F = v_F.clone(); let v_K = v_K.clone(); let v_S = v_S.clone(); let v__ = v__.clone(); let v_grid = v_grid.clone(); let v_hi = v_hi.clone(); let v_i = v_i.clone(); let v_lo = v_lo.clone(); let v_svx = v_svx.clone(); let v_sxx = v_sxx.clone(); let v_v = v_v.clone(); py::func(move |a: &[Py]| -> R {
            #[allow(unused_mut, unused_variables, unused_assignments)]
{
            let mut v_k: Py = a.get(0).cloned().unwrap_or(Py::None);
            return Ok(py::getitem(&f__mm_varpro(v_S.clone(), v_v.clone(), v_k.clone())?, &Py::Int(2i64))?);
            #[allow(unreachable_code)]
            Ok(Py::None)
        }
        }) }, py::getitem(&v_grid.clone(), &py::sub(&v_i.clone(), &Py::Int(1i64))?)?, py::getitem(&v_grid.clone(), &py::add(&v_i.clone(), &Py::Int(1i64))?)?, Py::Int(400i64))?;
    let t4246 = (f__mm_varpro(v_S.clone(), v_v.clone(), v_K.clone())?).items()?;
    if t4246.len() != 3 { return Err(PyErr::Value("unpack".into())); }
    v_svx = t4246[0].clone();
    v_sxx = t4246[1].clone();
    v__ = t4246[2].clone();
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
    let _ = f__need({ let t4311 = v_t.clone(); Py::Bool({ let t4312 = Py::Int(0i64); if py::ge(&t4311, &t4312)? { true } else { false } }) }, Py::str("t must be zero or positive"))?;
    v_vt = py::mul(&v_Vmax.clone(), &v_t.clone())?;
    let _ = f__need({ let t4313 = py::truediv(&v_vt.clone(), &v_Km.clone())?; Py::Bool({ let t4314 = Py::Int(600i64); if py::le(&t4313, &t4314)? { true } else { false } }) }, Py::str("the substrate is exhausted beyond representable precision (Vmax t / Km > 600)"))?;
    if (v_first_order.clone()).truthy() {
        v_S = py::mul(&v_S0.clone(), &py::math1("exp", &py::truediv(&py::neg(&v_vt.clone())?, &v_Km.clone())?)?)?;
        return Ok(Py::dict(vec![(Py::str("S"), v_S.clone()), (Py::str("P"), py::sub(&v_S0.clone(), &v_S.clone())?)]));
    }
    if ({ let t4315 = v_vt.clone(); Py::Bool({ let t4316 = Py::Int(0i64); if py::eq(&t4315, &t4316) { true } else { false } }) }).truthy() {
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
    if ({ let t4317 = py::call(&v_fP.clone(), &[v_Ph.clone()])?; Py::Bool({ let t4318 = Py::Int(0i64); if py::ge(&t4317, &t4318)? { true } else { false } }) }).truthy() {
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
    let _ = f__need({ let t4409 = Py::Int(0i64); Py::Bool({ let t4410 = v_X.clone(); if py::lt(&t4409, &t4410)? { { let t4411 = Py::Int(1i64); if py::lt(&t4410, &t4411)? { true } else { false } } } else { false } }) }, Py::str("the fraction converted must lie strictly between 0 and 1"))?;
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
    let _ = f__need({ let t4470 = Py::Bool(matches!(v_N.clone(), Py::Int(_) | Py::Big(_) | Py::Bool(_))); if !t4470.truthy() { t4470 } else { { let t4467 = Py::Int(2i64); Py::Bool({ let t4468 = v_N.clone(); if py::le(&t4467, &t4468)? { { let t4469 = Py::Int(100000i64); if py::le(&t4468, &t4469)? { true } else { false } } } else { false } }) } } }, Py::str("N must be a whole number of individuals (2 to 100000)"))?;
    let _ = f__need({ let t4474 = Py::Bool(matches!(v_i.clone(), Py::Int(_) | Py::Big(_) | Py::Bool(_))); if !t4474.truthy() { t4474 } else { { let t4471 = Py::Int(0i64); Py::Bool({ let t4472 = v_i.clone(); if py::lt(&t4471, &t4472)? { { let t4473 = v_N.clone(); if py::lt(&t4472, &t4473)? { true } else { false } } } else { false } }) } } }, Py::str("i must be a whole number of mutants between 1 and N - 1"))?;
    let _ = f__pos(Py::str("r"), Py::tuple(vec![v_r.clone()]))?;
    if ({ let t4481 = v_r.clone(); Py::Bool({ let t4482 = Py::Int(1i64); if py::eq(&t4481, &t4482) { true } else { false } }) }).truthy() {
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
    let _ = f__need({ let t4545 = Py::Int(0i64); Py::Bool({ let t4546 = v_mu.clone(); if py::lt(&t4545, &t4546)? { { let t4547 = Py::Float(0.01); if py::le(&t4546, &t4547)? { true } else { false } } } else { false } }) }, Py::str("mu must lie in (0, 0.01]"))?;
    let _ = f__need({ let t4551 = Py::Bool(matches!(v_s.clone(), Py::Int(_) | Py::Big(_) | Py::Bool(_) | Py::Float(_))); if !t4551.truthy() { t4551 } else { { let t4548 = Py::Float(0.01); Py::Bool({ let t4549 = v_s.clone(); if py::le(&t4548, &t4549)? { { let t4550 = Py::Int(1i64); if py::le(&t4549, &t4550)? { true } else { false } } } else { false } }) } } }, Py::str("s must lie in [0.01, 1]"))?;
    let _ = f__need({ let t4555 = Py::Bool(matches!(v_h.clone(), Py::Int(_) | Py::Big(_) | Py::Bool(_) | Py::Float(_))); if !t4555.truthy() { t4555 } else { { let t4552 = Py::Int(0i64); Py::Bool({ let t4553 = v_h.clone(); if py::le(&t4552, &t4553)? { { let t4554 = Py::Int(1i64); if py::le(&t4553, &t4554)? { true } else { false } } } else { false } }) } } }, Py::str("h must lie in [0, 1]"))?;
    v_f = { let v_h = v_h.clone(); let v_mu = v_mu.clone(); let v_s = v_s.clone(); py::func(move |a: &[Py]| -> R {
            #[allow(unused_mut, unused_variables, unused_assignments)]
{
            let mut v_q: Py = a.get(0).cloned().unwrap_or(Py::None);
            return Ok(py::sub(&f__mutsel_map(v_q.clone(), v_mu.clone(), v_s.clone(), v_h.clone())?, &v_q.clone())?);
            #[allow(unreachable_code)]
            Ok(Py::None)
        }
        }) };
    let _ = f__need({ let t4556 = py::call(&v_f.clone(), &[Py::Float(0.5)])?; Py::Bool({ let t4557 = Py::Int(0i64); if py::lt(&t4556, &t4557)? { true } else { false } }) }, Py::str("mutation is too strong relative to selection: no low-frequency balance"))?;
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
    let _ = f__need({ let t4641 = Py::Bool(matches!(v_sizes.clone(), Py::List(_))); if !t4641.truthy() { t4641 } else { { let t4639 = py::len(&v_sizes.clone())?; Py::Bool({ let t4640 = Py::Int(2i64); if py::ge(&t4639, &t4640)? { true } else { false } }) } } }, Py::str("give the sizes of at least 2 generations"))?;
    v_ns = { let mut out = vec![]; for t4642 in py::iter(&(v_sizes.clone()))? { let mut v_v: Py; 
v_v = t4642;
out.push(f__whole(Py::str("sizes"), v_v.clone(), Py::Int(1i64))?); } Py::list(out) };
    v_ne = py::truediv(&py::fraction(&[py::len(&v_ns.clone())?])?, &py::sum(&{ let mut out = vec![]; for t4651 in py::iter(&(v_ns.clone()))? { let mut v_n: Py; 
v_n = t4651;
out.push(py::fraction(&[Py::Int(1i64), v_n.clone()])?); } Py::list(out) }, None)?)?;
    v_ret = Py::Float(1.0);
    for t4652 in py::iter(&(v_ns.clone()))? {
        v_n = t4652;
        v_ret = py::mul(&v_ret, &py::sub(&Py::Int(1i64), &py::truediv(&Py::Int(1i64), &py::mul(&Py::Int(2i64), &v_n.clone())?)?)?)?;
    }
    return Ok(f__out(Py::dict(vec![(Py::str("Ne"), py::float(&v_ne.clone())?), (Py::str("retained"), v_ret.clone())]))?);
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn f_ne_sex_ratio(mut v_N_m: Py, mut v_N_f: Py) -> R {
    let mut v_ne: Py = Py::None;
    let _ = f__need({ let t4714 = Py::Bool(matches!(v_N_m.clone(), Py::Int(_) | Py::Big(_) | Py::Bool(_))); if !t4714.truthy() { t4714 } else { { let t4713 = Py::Bool(matches!(v_N_f.clone(), Py::Int(_) | Py::Big(_) | Py::Bool(_))); if !t4713.truthy() { t4713 } else { { let t4712 = { let t4708 = v_N_m.clone(); Py::Bool({ let t4709 = Py::Int(1i64); if py::ge(&t4708, &t4709)? { true } else { false } }) }; if !t4712.truthy() { t4712 } else { { let t4710 = v_N_f.clone(); Py::Bool({ let t4711 = Py::Int(1i64); if py::ge(&t4710, &t4711)? { true } else { false } }) } } } } } } }, Py::str("N_m and N_f must be whole numbers of breeding males and females (at least 1)"))?;
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
    v_d = { let mut out = vec![]; for t4773 in py::range_iter(&[(v_n.clone()).i()?]) { let mut v_i: Py; 
v_i = t4773;
for t4774 in py::range_iter(&[(v_n.clone()).i()?]) { let mut v_j: Py; 
v_j = t4774;
if ({ let t4775 = v_i.clone(); Py::Bool({ let t4776 = v_j.clone(); if py::ne(&t4775, &t4776) { true } else { false } }) }).truthy() { out.push((Py::tuple(vec![v_i.clone(), v_j.clone()]), py::getitem(&py::getitem(&v_M.clone(), &v_i.clone())?, &v_j.clone())?)); }}} Py::dict(out) };
    v_nodes = Py::list((Py::list((py::range_iter(&[(v_n.clone()).i()?])).into_iter().collect::<Vec<Py>>())).items()?);
    v_branch = Py::dict(vec![]);
    v_edges = Py::list(vec![]);
    v_nxt = v_n.clone();
    v_scale = py::minmax(&[py::method(&v_d.clone(), "values", &[])?], true)?;
    while ({ let t4814 = py::len(&v_nodes.clone())?; Py::Bool({ let t4815 = Py::Int(2i64); if py::gt(&t4814, &t4815)? { true } else { false } }) }).truthy() {
        v_r = py::len(&v_nodes.clone())?;
        v_R = { let mut out = vec![]; for t4777 in py::iter(&(v_nodes.clone()))? { let mut v_i: Py; 
v_i = t4777;
out.push((v_i.clone(), py::sum(&{ let mut out = vec![]; for t4778 in py::iter(&(v_nodes.clone()))? { let mut v_k: Py; 
v_k = t4778;
if ({ let t4779 = v_k.clone(); Py::Bool({ let t4780 = v_i.clone(); if py::ne(&t4779, &t4780) { true } else { false } }) }).truthy() { out.push(py::getitem(&v_d.clone(), &Py::tuple(vec![v_i.clone(), v_k.clone()]))?); }} Py::list(out) }, None)?)); } Py::dict(out) };
        v_Qv = { let mut out = vec![]; for t4781 in py::iter(&(py::enumerate(&v_nodes.clone(), 0)?))? { let mut v_a: Py; let mut v_i: Py; 
let t4782 = (t4781).items()?;
if t4782.len() != 2 { return Err(PyErr::Value("unpack".into())); }
v_a = t4782[0].clone();
v_i = t4782[1].clone();
for t4783 in py::iter(&(py::slice(&v_nodes.clone(), Some((py::add(&v_a.clone(), &Py::Int(1i64))?).i()?), None, None)?))? { let mut v_j: Py; 
v_j = t4783;
out.push((Py::tuple(vec![v_i.clone(), v_j.clone()]), py::sub(&py::sub(&py::mul(&py::sub(&v_r.clone(), &Py::Int(2i64))?, &py::getitem(&v_d.clone(), &Py::tuple(vec![v_i.clone(), v_j.clone()]))?)?, &py::getitem(&v_R.clone(), &v_i.clone())?)?, &py::getitem(&v_R.clone(), &v_j.clone())?)?)); }} Py::dict(out) };
        v_best = py::minmax(&[py::method(&v_Qv.clone(), "values", &[])?], false)?;
        let t4788 = (py::minmax(&[{ let mut out = vec![]; for t4784 in py::iter(&(py::method(&v_Qv.clone(), "items", &[])?))? { let mut v_k: Py; let mut v_v: Py; 
let t4785 = (t4784).items()?;
if t4785.len() != 2 { return Err(PyErr::Value("unpack".into())); }
v_k = t4785[0].clone();
v_v = t4785[1].clone();
if ({ let t4786 = v_v.clone(); Py::Bool({ let t4787 = v_best.clone(); if py::eq(&t4786, &t4787) { true } else { false } }) }).truthy() { out.push(v_k.clone()); }} Py::list(out) }], false)?).items()?;
        if t4788.len() != 2 { return Err(PyErr::Value("unpack".into())); }
        v_i = t4788[0].clone();
        v_j = t4788[1].clone();
        v_tied = py::set(&{ let mut out = vec![]; for t4789 in py::iter(&(py::method(&v_Qv.clone(), "items", &[])?))? { let mut v_k: Py; let mut v_v: Py; 
let t4790 = (t4789).items()?;
if t4790.len() != 2 { return Err(PyErr::Value("unpack".into())); }
v_k = t4790[0].clone();
v_v = t4790[1].clone();
if ({ let t4791 = py::abs(&py::sub(&v_v.clone(), &v_best.clone())?)?; Py::Bool({ let t4792 = py::mul(&py::mul(&Py::Float(1e-12), &v_scale.clone())?, &v_r.clone())?; if py::le(&t4791, &t4792)? { true } else { false } }) }).truthy() { out.push(v_k.clone()); }} Py::list(out) })?;
        v_harmless = (if ({ let t4793 = v_r.clone(); Py::Bool({ let t4794 = Py::Int(4i64); if py::eq(&t4793, &t4794) { true } else { false } }) }).truthy() { py::set(&Py::list(vec![Py::tuple(vec![v_i.clone(), v_j.clone()]), Py::tuple((py::sorted(&{ let mut out = vec![]; for t4795 in py::iter(&(v_nodes.clone()))? { let mut v_k: Py; 
v_k = t4795;
if ({ let t4796 = v_k.clone(); Py::Bool({ let t4797 = Py::tuple(vec![v_i.clone(), v_j.clone()]); if !py::contains(&t4797, &t4796)? { true } else { false } }) }).truthy() { out.push(v_k.clone()); }} Py::list(out) }, (Py::Bool(false)).truthy())?).items()?)]))? } else { py::set(&Py::list(vec![Py::tuple(vec![v_i.clone(), v_j.clone()])]))? });
        let _ = f__need({ let t4802 = { let t4798 = v_r.clone(); Py::Bool({ let t4799 = Py::Int(3i64); if py::eq(&t4798, &t4799) { true } else { false } }) }; if t4802.truthy() { t4802 } else { { let t4800 = v_tied.clone(); Py::Bool({ let t4801 = v_harmless.clone(); if py::le(&t4800, &t4801)? { true } else { false } }) } } }, Py::str("two pairs tie in the neighbour-joining criterion: the tree is not unique"))?;
        if (v_corrected.clone()).truthy() {
            v_li = py::add(&py::truediv(&py::getitem(&v_d.clone(), &Py::tuple(vec![v_i.clone(), v_j.clone()]))?, &Py::Int(2i64))?, &py::truediv(&py::sub(&py::getitem(&v_R.clone(), &v_i.clone())?, &py::getitem(&v_R.clone(), &v_j.clone())?)?, &py::mul(&Py::Int(2i64), &py::sub(&v_r.clone(), &Py::Int(2i64))?)?)?)?;
        } else {
            v_li = py::truediv(&py::getitem(&v_d.clone(), &Py::tuple(vec![v_i.clone(), v_j.clone()]))?, &Py::Int(2i64))?;
        }
        v_lj = py::sub(&py::getitem(&v_d.clone(), &Py::tuple(vec![v_i.clone(), v_j.clone()]))?, &v_li.clone())?;
        for t4803 in py::iter(&(Py::tuple(vec![Py::tuple(vec![v_i.clone(), v_li.clone()]), Py::tuple(vec![v_j.clone(), v_lj.clone()])])))? {
            let t4806 = (t4803).items()?;
            if t4806.len() != 2 { return Err(PyErr::Value("unpack".into())); }
            v_x = t4806[0].clone();
            v_lx = t4806[1].clone();
            let _ = py::method(&v_edges.clone(), "append", &[v_lx.clone()])?;
            if ({ let t4804 = v_x.clone(); Py::Bool({ let t4805 = v_n.clone(); if py::lt(&t4804, &t4805)? { true } else { false } }) }).truthy() {
                py::setitem(&v_branch.clone(), &v_x.clone(), v_lx.clone())?;
            }
        }
        v_u = v_nxt.clone();
        v_nxt = py::iadd(&v_nxt, &Py::Int(1i64))?;
        for t4807 in py::iter(&(v_nodes.clone()))? {
            v_k = t4807;
            if ({ let t4808 = v_k.clone(); Py::Bool({ let t4809 = Py::tuple(vec![v_i.clone(), v_j.clone()]); if !py::contains(&t4809, &t4808)? { true } else { false } }) }).truthy() {
                let t4810 = py::truediv(&py::sub(&py::add(&py::getitem(&v_d.clone(), &Py::tuple(vec![v_i.clone(), v_k.clone()]))?, &py::getitem(&v_d.clone(), &Py::tuple(vec![v_j.clone(), v_k.clone()]))?)?, &py::getitem(&v_d.clone(), &Py::tuple(vec![v_i.clone(), v_j.clone()]))?)?, &Py::Int(2i64))?;
                py::setitem(&v_d.clone(), &Py::tuple(vec![v_u.clone(), v_k.clone()]), t4810.clone())?;
                py::setitem(&v_d.clone(), &Py::tuple(vec![v_k.clone(), v_u.clone()]), t4810.clone())?;
            }
        }
        v_nodes = py::add(&{ let mut out = vec![]; for t4811 in py::iter(&(v_nodes.clone()))? { let mut v_k: Py; 
v_k = t4811;
if ({ let t4812 = v_k.clone(); Py::Bool({ let t4813 = Py::tuple(vec![v_i.clone(), v_j.clone()]); if !py::contains(&t4813, &t4812)? { true } else { false } }) }).truthy() { out.push(v_k.clone()); }} Py::list(out) }, &Py::list(vec![v_u.clone()]))?;
    }
    let t4816 = (v_nodes.clone()).items()?;
    if t4816.len() != 2 { return Err(PyErr::Value("unpack".into())); }
    v_a = t4816[0].clone();
    v_b = t4816[1].clone();
    let _ = py::method(&v_edges.clone(), "append", &[py::getitem(&v_d.clone(), &Py::tuple(vec![v_a.clone(), v_b.clone()]))?])?;
    for t4817 in py::iter(&(Py::tuple(vec![v_a.clone(), v_b.clone()])))? {
        v_x = t4817;
        if ({ let t4818 = v_x.clone(); Py::Bool({ let t4819 = v_n.clone(); if py::lt(&t4818, &t4819)? { true } else { false } }) }).truthy() {
            py::setitem(&v_branch.clone(), &v_x.clone(), py::getitem(&v_d.clone(), &Py::tuple(vec![v_a.clone(), v_b.clone()]))?)?;
        }
    }
    return Ok(f__out(Py::dict(vec![(Py::str("total_length"), py::sum(&v_edges.clone(), None)?), (Py::str("branch"), { let mut out = vec![]; for t4834 in py::range_iter(&[(v_n.clone()).i()?]) { let mut v_i: Py; 
v_i = t4834;
out.push(py::getitem(&v_branch.clone(), &v_i.clone())?); } Py::list(out) })]))?);
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn f_net_charge(mut v_seq: Py, mut v_pH: Py) -> R {
    let mut v_s: Py = Py::None;
    v_s = f__seq(v_seq.clone(), Py::str("seq"), Py::str("PROTEIN"))?;
    let _ = f__fin(Py::str("pH"), Py::tuple(vec![v_pH.clone()]))?;
    let _ = f__need({ let t5575 = Py::Int(0i64); Py::Bool({ let t5576 = v_pH.clone(); if py::le(&t5575, &t5576)? { { let t5577 = Py::Int(14i64); if py::le(&t5576, &t5577)? { true } else { false } } } else { false } }) }, Py::str("pH must lie between 0 and 14"))?;
    return Ok(f__out(Py::dict(vec![(Py::str("charge"), f__charge(v_s.clone(), v_pH.clone(), Py::Bool(true))?)]))?);
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn f_nucleic_mass(mut v_seq: Py, mut v_molecule: Py, mut v_strands: Py) -> R {
    let mut v_s: Py = Py::None;
    let mut v_table: Py = Py::None;
    let mut v_w: Py = Py::None;
    let _ = f__need({ let t4876 = v_molecule.clone(); Py::Bool({ let t4877 = Py::tuple(vec![Py::str("DNA"), Py::str("RNA")]); if py::contains(&t4877, &t4876)? { true } else { false } }) }, Py::str("molecule must be DNA or RNA"))?;
    let _ = f__need({ let t4878 = v_strands.clone(); Py::Bool({ let t4879 = Py::tuple(vec![Py::str("single"), Py::str("double")]); if py::contains(&t4879, &t4878)? { true } else { false } }) }, Py::str("strands must be single or double"))?;
    v_s = f__seq(v_seq.clone(), Py::str("seq"), v_molecule.clone())?;
    v_table = (if ({ let t4903 = v_molecule.clone(); Py::Bool({ let t4904 = Py::str("DNA"); if py::eq(&t4903, &t4904) { true } else { false } }) }).truthy() { c_DNA_W()? } else { c_RNA_W()? });
    v_w = f__chain_mass(v_s.clone(), v_table.clone())?;
    if ({ let t4906 = v_strands.clone(); Py::Bool({ let t4907 = Py::str("double"); if py::eq(&t4906, &t4907) { true } else { false } }) }).truthy() {
        v_w = py::iadd(&v_w, &f__chain_mass(py::method(&v_s.clone(), "translate", &[(if ({ let t4908 = v_molecule.clone(); Py::Bool({ let t4909 = Py::str("DNA"); if py::eq(&t4908, &t4909) { true } else { false } }) }).truthy() { c__COMP()? } else { c__COMP_RNA()? })])?, v_table.clone())?)?;
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
    v_names = py::bitor(&py::set(&v_parents.clone())?, &py::set(&{ let mut out = vec![]; for t5041 in py::iter(&(py::method(&v_parents.clone(), "values", &[])?))? { let mut v_ps: Py; 
v_ps = t5041;
for t5042 in py::iter(&(v_ps.clone()))? { let mut v_p: Py; 
v_p = t5042;
out.push(v_p.clone()); }} Py::list(out) })?)?;
    let _ = f__need({ let t5043 = v_individual.clone(); Py::Bool({ let t5044 = v_names.clone(); if py::contains(&t5044, &t5043)? { true } else { false } }) }, Py::str(&[py::format(&v_individual.clone(), "")?, " is not in the pedigree".to_string()].concat()))?;
    v_memo = Py::dict(vec![]);
    let rec5045 = py::rec_cell();
    v_up_paths = { let v_half = v_half.clone(); let v_individual = v_individual.clone(); let v_memo = v_memo.clone(); let v_names = v_names.clone(); let v_parents = v_parents.clone(); let v_pedigree = v_pedigree.clone(); let rec5045 = rec5045.clone(); py::func(move |a: &[Py]| -> R {
            #[allow(unused_mut, unused_variables, unused_assignments)]
{
            let mut v_x: Py = a.get(0).cloned().unwrap_or(Py::None);
            let v_up_paths: Py = py::rec_get(&rec5045)?;
            let mut v_out: Py = Py::None;
            let mut v_p: Py = Py::None;
            v_out = Py::list(vec![Py::list(vec![v_x.clone()])]);
            for t5046 in py::iter(&(py::method(&v_parents.clone(), "get", &[v_x.clone(), Py::list(vec![])])?))? {
                v_p = t5046;
                v_out = py::iadd(&v_out, &{ let mut out = vec![]; for t5047 in py::iter(&(py::call(&v_up_paths.clone(), &[v_p.clone()])?))? { let mut v_q: Py; 
v_q = t5047;
out.push(py::add(&Py::list(vec![v_x.clone()]), &v_q.clone())?); } Py::list(out) })?;
            }
            return Ok(v_out.clone());
            #[allow(unreachable_code)]
            Ok(Py::None)
        }
        }) };
    py::rec_set(&rec5045, &v_up_paths);
    let rec5048 = py::rec_cell();
    v_F = { let v_half = v_half.clone(); let v_individual = v_individual.clone(); let v_memo = v_memo.clone(); let v_names = v_names.clone(); let v_parents = v_parents.clone(); let v_pedigree = v_pedigree.clone(); let v_up_paths = v_up_paths.clone(); let rec5048 = rec5048.clone(); py::func(move |a: &[Py]| -> R {
            #[allow(unused_mut, unused_variables, unused_assignments)]
{
            let mut v_x: Py = a.get(0).cloned().unwrap_or(Py::None);
            let v_F: Py = py::rec_get(&rec5048)?;
            let mut v_a: Py = Py::None;
            let mut v_b: Py = Py::None;
            let mut v_d: Py = Py::None;
            let mut v_pd: Py = Py::None;
            let mut v_ps: Py = Py::None;
            let mut v_s: Py = Py::None;
            let mut v_tot: Py = Py::None;
            if ({ let t5049 = v_x.clone(); Py::Bool({ let t5050 = v_parents.clone(); if !py::contains(&t5050, &t5049)? { true } else { false } }) }).truthy() {
                return Ok(py::fraction(&[Py::Int(0i64)])?);
            }
            if ({ let t5051 = v_x.clone(); Py::Bool({ let t5052 = v_memo.clone(); if py::contains(&t5052, &t5051)? { true } else { false } }) }).truthy() {
                return Ok(py::getitem(&v_memo.clone(), &v_x.clone())?);
            }
            let t5053 = (py::getitem(&v_parents.clone(), &v_x.clone())?).items()?;
            if t5053.len() != 2 { return Err(PyErr::Value("unpack".into())); }
            v_s = t5053[0].clone();
            v_d = t5053[1].clone();
            v_tot = py::fraction(&[Py::Int(0i64)])?;
            let t5054 = (Py::tuple(vec![py::call(&v_up_paths.clone(), &[v_s.clone()])?, py::call(&v_up_paths.clone(), &[v_d.clone()])?])).items()?;
            if t5054.len() != 2 { return Err(PyErr::Value("unpack".into())); }
            v_ps = t5054[0].clone();
            v_pd = t5054[1].clone();
            for t5055 in py::iter(&(v_ps.clone()))? {
                v_a = t5055;
                for t5056 in py::iter(&(v_pd.clone()))? {
                    v_b = t5056;
                    if ({ let t5061 = { let t5057 = py::getitem(&v_a.clone(), &py::neg(&Py::Int(1i64))?)?; Py::Bool({ let t5058 = py::getitem(&v_b.clone(), &py::neg(&Py::Int(1i64))?)?; if py::eq(&t5057, &t5058) { true } else { false } }) }; if !t5061.truthy() { t5061 } else { { let t5059 = py::bitand(&py::set(&v_a.clone())?, &py::set(&v_b.clone())?)?; Py::Bool({ let t5060 = py::set(&Py::list(vec![py::getitem(&v_a.clone(), &py::neg(&Py::Int(1i64))?)?]))?; if py::eq(&t5059, &t5060) { true } else { false } }) } } }).truthy() {
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
    py::rec_set(&rec5048, &v_F);
    return Ok(Py::dict(vec![(Py::str("F"), py::float(&py::call(&v_F.clone(), &[v_individual.clone()])?)?)]));
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn f_petersen_mle(mut v_M: Py, mut v_C: Py, mut v_R: Py, mut v_ceil: Py) -> R {
    let mut v_nm: Py = Py::None;
    let mut v_v: Py = Py::None;
    for t5103 in py::iter(&(Py::tuple(vec![Py::tuple(vec![Py::str("M"), v_M.clone()]), Py::tuple(vec![Py::str("C"), v_C.clone()]), Py::tuple(vec![Py::str("R"), v_R.clone()])])))? {
        let t5107 = (t5103).items()?;
        if t5107.len() != 2 { return Err(PyErr::Value("unpack".into())); }
        v_nm = t5107[0].clone();
        v_v = t5107[1].clone();
        let _ = f__need({ let t5106 = Py::Bool(matches!(v_v.clone(), Py::Int(_) | Py::Big(_) | Py::Bool(_))); if !t5106.truthy() { t5106 } else { { let t5104 = v_v.clone(); Py::Bool({ let t5105 = Py::Int(1i64); if py::ge(&t5104, &t5105)? { true } else { false } }) } } }, Py::str(&[py::format(&v_nm.clone(), "")?, " must be a whole number >= 1".to_string()].concat()))?;
    }
    let _ = f__need({ let t5108 = v_R.clone(); Py::Bool({ let t5109 = py::minmax(&[v_M.clone(), v_C.clone()], false)?; if py::le(&t5108, &t5109)? { true } else { false } }) }, Py::str("R cannot exceed M or C"))?;
    return Ok(Py::dict(vec![(Py::str("N_mle"), (if (v_ceil.clone()).truthy() { py::neg(&py::floordiv(&py::mul(&py::neg(&v_M.clone())?, &v_C.clone())?, &v_R.clone())?)? } else { py::floordiv(&py::mul(&v_M.clone(), &v_C.clone())?, &v_R.clone())? }))]));
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn f_protein_global(mut v_seq1: Py, mut v_seq2: Py, mut v_open: Py, mut v_extend: Py) -> R {
    let mut v_a: Py = Py::None;
    let mut v_b: Py = Py::None;
    let t5253 = (Py::tuple(vec![f__seq(v_seq1.clone(), Py::str("seq1"), Py::str("PROTEIN"))?, f__seq(v_seq2.clone(), Py::str("seq2"), Py::str("PROTEIN"))?])).items()?;
    if t5253.len() != 2 { return Err(PyErr::Value("unpack".into())); }
    v_a = t5253[0].clone();
    v_b = t5253[1].clone();
    let _ = f__check_align(v_a.clone(), v_b.clone(), Py::dict(vec![(Py::str("open"), v_open.clone()), (Py::str("extend"), v_extend.clone())]))?;
    let _ = f__need({ let t5265 = { let t5261 = v_open.clone(); Py::Bool({ let t5262 = Py::Int(0i64); if py::le(&t5261, &t5262)? { true } else { false } }) }; if !t5265.truthy() { t5265 } else { { let t5263 = v_extend.clone(); Py::Bool({ let t5264 = Py::Int(0i64); if py::le(&t5263, &t5264)? { true } else { false } }) } } }, Py::str("gap scores must be zero or negative"))?;
    return Ok(f__out(Py::dict(vec![(Py::str("score"), f__align(v_a.clone(), v_b.clone(), py::func(|a: &[Py]| -> R { f__blo(a.get(0).cloned().unwrap_or(Py::None), a.get(1).cloned().unwrap_or(Py::None)) }), v_open.clone(), v_extend.clone(), Py::str("global"), Py::Bool(false))?)]))?);
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn f_protein_local(mut v_seq1: Py, mut v_seq2: Py, mut v_open: Py, mut v_extend: Py) -> R {
    let mut v_a: Py = Py::None;
    let mut v_b: Py = Py::None;
    let t5373 = (Py::tuple(vec![f__seq(v_seq1.clone(), Py::str("seq1"), Py::str("PROTEIN"))?, f__seq(v_seq2.clone(), Py::str("seq2"), Py::str("PROTEIN"))?])).items()?;
    if t5373.len() != 2 { return Err(PyErr::Value("unpack".into())); }
    v_a = t5373[0].clone();
    v_b = t5373[1].clone();
    let _ = f__check_align(v_a.clone(), v_b.clone(), Py::dict(vec![(Py::str("open"), v_open.clone()), (Py::str("extend"), v_extend.clone())]))?;
    let _ = f__need({ let t5385 = { let t5381 = v_open.clone(); Py::Bool({ let t5382 = Py::Int(0i64); if py::le(&t5381, &t5382)? { true } else { false } }) }; if !t5385.truthy() { t5385 } else { { let t5383 = v_extend.clone(); Py::Bool({ let t5384 = Py::Int(0i64); if py::le(&t5383, &t5384)? { true } else { false } }) } } }, Py::str("gap scores must be zero or negative"))?;
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
    let _ = f__need({ let t5665 = Py::Bool(matches!(v_n.clone(), Py::Int(_) | Py::Big(_) | Py::Bool(_))); if !t5665.truthy() { t5665 } else { { let t5662 = Py::Int(1i64); Py::Bool({ let t5663 = v_n.clone(); if py::le(&t5662, &t5663)? { { let t5664 = v_N.clone(); if py::le(&t5663, &t5664)? { true } else { false } } } else { false } }) } } }, Py::str("n must be a whole number from 1 to the number of individuals"))?;
    let _ = f__need({ let t5666 = v_N.clone(); Py::Bool({ let t5667 = Py::Int(100000i64); if py::le(&t5666, &t5667)? { true } else { false } }) }, Py::str("at most 100000 individuals"))?;
    if (v_with_replacement.clone()).truthy() {
        return Ok(Py::dict(vec![(Py::str("expected_species"), py::fsum(&{ let mut out = vec![]; for t5668 in py::iter(&(v_c.clone()))? { let mut v_x: Py; 
v_x = t5668;
out.push(py::sub(&Py::Int(1i64), &py::pow(&py::sub(&Py::Int(1i64), &py::truediv(&v_x.clone(), &v_N.clone())?)?, &v_n.clone())?)?); } Py::list(out) })?)]));
    }
    v_tot = py::comb(&v_N.clone(), &v_n.clone())?;
    v_e = py::sum(&{ let mut out = vec![]; for t5669 in py::iter(&(v_c.clone()))? { let mut v_x: Py; 
v_x = t5669;
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
    let _ = f__need({ let t5726 = Py::Bool(matches!(v_t.clone(), Py::Int(_) | Py::Big(_) | Py::Bool(_))); if !t5726.truthy() { t5726 } else { { let t5723 = Py::Int(0i64); Py::Bool({ let t5724 = v_t.clone(); if py::le(&t5723, &t5724)? { { let t5725 = Py::Int(100000i64); if py::le(&t5724, &t5725)? { true } else { false } } } else { false } }) } } }, Py::str("t must be a whole number of generations, 0 to 100000"))?;
    v_q = v_q0.clone();
    for t5727 in py::range_iter(&[(v_t.clone()).i()?]) {
        v__ = t5727;
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
    let _ = f__need({ let t5786 = v_d_cM.clone(); Py::Bool({ let t5787 = Py::Int(0i64); if py::ge(&t5786, &t5787)? { true } else { false } }) }, Py::str("the map distance must be zero or positive"))?;
    let _ = f__need({ let t5788 = v_map.clone(); Py::Bool({ let t5789 = Py::tuple(vec![Py::str("haldane"), Py::str("kosambi")]); if py::contains(&t5789, &t5788)? { true } else { false } }) }, Py::str("map must be haldane or kosambi"))?;
    v_m = py::truediv(&v_d_cM.clone(), &Py::Int(100i64))?;
    return Ok(f__out(Py::dict(vec![(Py::str("r"), (if ({ let t5804 = v_map.clone(); Py::Bool({ let t5805 = Py::str("haldane"); if py::eq(&t5804, &t5805) { true } else { false } }) }).truthy() { py::mul(&py::neg(&Py::Float(0.5))?, &py::math1("expm1", &py::mul(&py::neg(&Py::Int(2i64))?, &v_m.clone())?)?)? } else { py::mul(&Py::Float(0.5), &py::math1("tanh", &py::mul(&Py::Int(2i64), &v_m.clone())?)?)? }))]))?);
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
    let _ = f__need({ let t5870 = v_enzyme.clone(); Py::Bool({ let t5871 = c_ENZYMES()?; if py::contains(&t5871, &t5870)? { true } else { false } }) }, py::add(&Py::str("enzyme must be one of "), &py::method(&Py::str(", "), "join", &[c_ENZYMES()?])?)?)?;
    let t5872 = (py::getitem(&c_ENZYMES()?, &v_enzyme.clone())?).items()?;
    if t5872.len() != 2 { return Err(PyErr::Value("unpack".into())); }
    v_site = t5872[0].clone();
    v_cut = t5872[1].clone();
    v_pos = { let mut out = vec![]; for t5873 in py::range_iter(&[(py::add(&py::sub(&py::len(&v_s.clone())?, &py::len(&v_site.clone())?)?, &Py::Int(1i64))?).i()?]) { let mut v_i: Py; 
v_i = t5873;
if (py::method(&v_s.clone(), "startswith", &[v_site.clone(), v_i.clone()])?).truthy() { out.push(py::add(&py::add(&v_i.clone(), &v_cut.clone())?, &Py::Int(1i64))?); }} Py::list(out) };
    v_bounds = py::add(&py::add(&Py::list(vec![Py::Int(1i64)]), &v_pos.clone())?, &Py::list(vec![py::add(&py::len(&v_s.clone())?, &Py::Int(1i64))?]))?;
    return Ok(Py::dict(vec![(Py::str("n_sites"), py::len(&v_pos.clone())?), (Py::str("positions"), v_pos.clone()), (Py::str("fragments"), { let mut out = vec![]; for t5874 in py::range_iter(&[(py::sub(&py::len(&v_bounds.clone())?, &Py::Int(1i64))?).i()?]) { let mut v_i: Py; 
v_i = t5874;
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
    let t5998 = (Py::tuple(vec![py::sum(&v_c.clone(), None)?, py::len(&v_c.clone())?])).items()?;
    if t5998.len() != 2 { return Err(PyErr::Value("unpack".into())); }
    v_N = t5998[0].clone();
    v_S = t5998[1].clone();
    let _ = f__need({ let t5999 = v_N.clone(); Py::Bool({ let t6000 = Py::Int(2i64); if py::ge(&t5999, &t6000)? { true } else { false } }) }, Py::str("at least 2 individuals are needed (ln N must be positive)"))?;
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
    let t6100 = (f__parse_newick(v_tree1.clone())?).items()?;
    if t6100.len() != 2 { return Err(PyErr::Value("unpack".into())); }
    v_r1 = t6100[0].clone();
    v_l1 = t6100[1].clone();
    let t6101 = (f__parse_newick(v_tree2.clone())?).items()?;
    if t6101.len() != 2 { return Err(PyErr::Value("unpack".into())); }
    v_r2 = t6101[0].clone();
    v_l2 = t6101[1].clone();
    let _ = f__need({ let t6102 = py::set(&v_l1.clone())?; Py::Bool({ let t6103 = py::set(&v_l2.clone())?; if py::eq(&t6102, &t6103) { true } else { false } }) }, Py::str("both trees must have the same taxa"))?;
    let _ = f__need({ let t6104 = py::len(&v_l1.clone())?; Py::Bool({ let t6105 = Py::Int(4i64); if py::ge(&t6104, &t6105)? { true } else { false } }) }, Py::str("at least 4 taxa are needed"))?;
    let t6113 = (Py::tuple(vec![f__splits(v_r1.clone(), v_l1.clone())?, f__splits(v_r2.clone(), v_l2.clone())?])).items()?;
    if t6113.len() != 2 { return Err(PyErr::Value("unpack".into())); }
    v_s1 = t6113[0].clone();
    v_s2 = t6113[1].clone();
    let _ = f__need({ let t6114 = v_s1.clone(); if t6114.truthy() { t6114 } else { v_s2.clone() } }, Py::str("both trees are unresolved stars"))?;
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
    let _ = f__need({ let t6162 = Py::Bool(matches!(v_C.clone(), Py::List(_))); if !t6162.truthy() { t6162 } else { { let t6161 = Py::Bool(matches!(v_R.clone(), Py::List(_))); if !t6161.truthy() { t6161 } else { { let t6160 = Py::Bool(matches!(v_M.clone(), Py::List(_))); if !t6160.truthy() { t6160 } else { { let t6156 = py::len(&v_C.clone())?; Py::Bool({ let t6157 = py::len(&v_R.clone())?; if py::eq(&t6156, &t6157) { { let t6158 = py::len(&v_M.clone())?; if py::eq(&t6157, &t6158) { { let t6159 = Py::Int(2i64); if py::ge(&t6158, &t6159)? { true } else { false } } } else { false } } } else { false } }) } } } } } } }, Py::str("C, R and M need one value per sampling occasion (at least 2)"))?;
    for t6163 in py::iter(&(py::zip(&[v_C.clone(), v_R.clone(), v_M.clone()])?))? {
        let t6180 = (t6163).items()?;
        if t6180.len() != 3 { return Err(PyErr::Value("unpack".into())); }
        v_c = t6180[0].clone();
        v_r = t6180[1].clone();
        v_m = t6180[2].clone();
        let _ = f__need({ let t6168 = { let t6164 = py::float(&v_c.clone())?; Py::Bool({ let t6165 = py::int(&v_c.clone())?; if py::eq(&t6164, &t6165) { true } else { false } }) }; if !t6168.truthy() { t6168 } else { { let t6166 = v_c.clone(); Py::Bool({ let t6167 = Py::Int(1i64); if py::ge(&t6166, &t6167)? { true } else { false } }) } } }, Py::str("catches must be whole numbers >= 1"))?;
        let _ = f__need({ let t6174 = { let t6169 = py::float(&v_r.clone())?; Py::Bool({ let t6170 = py::int(&v_r.clone())?; if py::eq(&t6169, &t6170) { true } else { false } }) }; if !t6174.truthy() { t6174 } else { { let t6171 = Py::Int(0i64); Py::Bool({ let t6172 = v_r.clone(); if py::le(&t6171, &t6172)? { { let t6173 = v_c.clone(); if py::le(&t6172, &t6173)? { true } else { false } } } else { false } }) } } }, Py::str("recaptures must be whole numbers from 0 to the catch"))?;
        let _ = f__need({ let t6179 = { let t6175 = py::float(&v_m.clone())?; Py::Bool({ let t6176 = py::int(&v_m.clone())?; if py::eq(&t6175, &t6176) { true } else { false } }) }; if !t6179.truthy() { t6179 } else { { let t6177 = v_m.clone(); Py::Bool({ let t6178 = Py::Int(0i64); if py::ge(&t6177, &t6178)? { true } else { false } }) } } }, Py::str("marked animals at large must be whole numbers >= 0"))?;
    }
    v_srm = py::sum(&{ let mut out = vec![]; for t6181 in py::iter(&(py::zip(&[v_R.clone(), v_M.clone()])?))? { let mut v_m: Py; let mut v_r: Py; 
let t6182 = (t6181).items()?;
if t6182.len() != 2 { return Err(PyErr::Value("unpack".into())); }
v_r = t6182[0].clone();
v_m = t6182[1].clone();
out.push(py::fraction(&[py::mul(&py::int(&v_r.clone())?, &py::int(&v_m.clone())?)?])?); } Py::list(out) }, None)?;
    let _ = f__need({ let t6183 = v_srm.clone(); Py::Bool({ let t6184 = Py::Int(0i64); if py::gt(&t6183, &t6184)? { true } else { false } }) }, Py::str("no marked animal was ever recaptured: the population cannot be estimated"))?;
    if (v_schnabel.clone()).truthy() {
        return Ok(Py::dict(vec![(Py::str("N"), py::float(&py::truediv(&py::sum(&{ let mut out = vec![]; for t6185 in py::iter(&(py::zip(&[v_C.clone(), v_M.clone()])?))? { let mut v_c: Py; let mut v_m: Py; 
let t6186 = (t6185).items()?;
if t6186.len() != 2 { return Err(PyErr::Value("unpack".into())); }
v_c = t6186[0].clone();
v_m = t6186[1].clone();
out.push(py::fraction(&[py::mul(&py::int(&v_c.clone())?, &py::int(&v_m.clone())?)?])?); } Py::list(out) }, None)?, &py::sum(&{ let mut out = vec![]; for t6187 in py::iter(&(v_R.clone()))? { let mut v_r: Py; 
v_r = t6187;
out.push(py::fraction(&[py::int(&v_r.clone())?])?); } Py::list(out) }, None)?)?)?)]));
    }
    return Ok(Py::dict(vec![(Py::str("N"), py::float(&py::truediv(&py::sum(&{ let mut out = vec![]; for t6188 in py::iter(&(py::zip(&[v_C.clone(), v_M.clone()])?))? { let mut v_c: Py; let mut v_m: Py; 
let t6189 = (t6188).items()?;
if t6189.len() != 2 { return Err(PyErr::Value("unpack".into())); }
v_c = t6189[0].clone();
v_m = t6189[1].clone();
out.push(py::fraction(&[py::mul(&py::int(&v_c.clone())?, &py::pow(&py::int(&v_m.clone())?, &Py::Int(2i64))?)?])?); } Py::list(out) }, None)?, &v_srm.clone())?)?)]));
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn f_semiglobal(mut v_seq1: Py, mut v_seq2: Py, mut v_match: Py, mut v_mismatch: Py, mut v_gap: Py) -> R {
    let mut v_a: Py = Py::None;
    let mut v_b: Py = Py::None;
    let t6321 = (Py::tuple(vec![f__seq(v_seq1.clone(), Py::str("seq1"), Py::str("ANY"))?, f__seq(v_seq2.clone(), Py::str("seq2"), Py::str("ANY"))?])).items()?;
    if t6321.len() != 2 { return Err(PyErr::Value("unpack".into())); }
    v_a = t6321[0].clone();
    v_b = t6321[1].clone();
    let _ = f__check_align(v_a.clone(), v_b.clone(), Py::dict(vec![(Py::str("match"), v_match.clone()), (Py::str("mismatch"), v_mismatch.clone()), (Py::str("gap"), v_gap.clone())]))?;
    let _ = f__need({ let t6329 = v_gap.clone(); Py::Bool({ let t6330 = Py::Int(0i64); if py::le(&t6329, &t6330)? { true } else { false } }) }, Py::str("the gap score must be zero or negative"))?;
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
    let _ = f__need({ let t6435 = v_base.clone(); Py::Bool({ let t6436 = c__LOGS()?; if py::contains(&t6436, &t6435)? { true } else { false } }) }, Py::str("base must be e, 2 or 10"))?;
    v_N = py::sum(&v_c.clone(), None)?;
    v_H = py::neg(&py::fsum(&{ let mut out = vec![]; for t6437 in py::iter(&(v_c.clone()))? { let mut v_x: Py; 
v_x = t6437;
out.push(py::mul(&py::truediv(&v_x.clone(), &v_N.clone())?, &py::math1("log", &py::truediv(&v_x.clone(), &v_N.clone())?)?)?); } Py::list(out) })?)?;
    v_S = py::len(&v_c.clone())?;
    v_Hb = py::truediv(&v_H.clone(), &py::math1("log", &py::getitem(&c__LOGS()?, &v_base.clone())?)?)?;
    v_J = (if ({ let t6438 = v_S.clone(); Py::Bool({ let t6439 = Py::Int(1i64); if py::eq(&t6438, &t6439) { true } else { false } }) }).truthy() { Py::Float(1.0) } else { (if (v_wrong.clone()).truthy() { py::truediv(&v_H.clone(), &v_S.clone())? } else { py::truediv(&v_H.clone(), &py::math1("log", &v_S.clone())?)? }) });
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
    let _ = f__need({ let t6513 = v_estimator.clone(); Py::Bool({ let t6514 = Py::tuple(vec![Py::str("infinite"), Py::str("finite")]); if py::contains(&t6514, &t6513)? { true } else { false } }) }, Py::str("estimator must be infinite or finite"))?;
    v_est = (if (v_swap.clone()).truthy() { py::getitem(&Py::dict(vec![(Py::str("infinite"), Py::str("finite")), (Py::str("finite"), Py::str("infinite"))]), &v_estimator.clone())? } else { v_estimator.clone() });
    v_N = py::sum(&v_c.clone(), None)?;
    if ({ let t6515 = v_est.clone(); Py::Bool({ let t6516 = Py::str("finite"); if py::eq(&t6515, &t6516) { true } else { false } }) }).truthy() {
        let _ = f__need({ let t6517 = v_N.clone(); Py::Bool({ let t6518 = Py::Int(2i64); if py::ge(&t6517, &t6518)? { true } else { false } }) }, Py::str("the finite-sample estimator needs at least 2 individuals"))?;
        v_D = py::fraction(&[py::sum(&{ let mut out = vec![]; for t6519 in py::iter(&(v_c.clone()))? { let mut v_x: Py; 
v_x = t6519;
out.push(py::mul(&v_x.clone(), &py::sub(&v_x.clone(), &Py::Int(1i64))?)?); } Py::list(out) }, None)?, py::mul(&v_N.clone(), &py::sub(&v_N.clone(), &Py::Int(1i64))?)?])?;
        let _ = f__need({ let t6520 = v_D.clone(); Py::Bool({ let t6521 = Py::Int(0i64); if py::gt(&t6520, &t6521)? { true } else { false } }) }, Py::str("every species is a singleton: the finite-sample dominance is 0 and its inverse undefined"))?;
    } else {
        v_D = py::fraction(&[py::sum(&{ let mut out = vec![]; for t6522 in py::iter(&(v_c.clone()))? { let mut v_x: Py; 
v_x = t6522;
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
    let _ = f__need({ let t6568 = Py::Bool(matches!(v_A.clone(), Py::List(_))); if !t6568.truthy() { t6568 } else { { let t6567 = Py::Bool(matches!(v_S.clone(), Py::List(_))); if !t6567.truthy() { t6567 } else { { let t6564 = py::len(&v_A.clone())?; Py::Bool({ let t6565 = py::len(&v_S.clone())?; if py::eq(&t6564, &t6565) { { let t6566 = Py::Int(2i64); if py::ge(&t6565, &t6566)? { true } else { false } } } else { false } }) } } } } }, Py::str("A and S must be lists of equal length"))?;
    for t6569 in py::iter(&(py::add(&v_A.clone(), &v_S.clone())?))? {
        v_v = t6569;
        let _ = f__pos(Py::str("areas and species counts"), Py::tuple(vec![v_v.clone()]))?;
    }
    let t6587 = (f__ols({ let mut out = vec![]; for t6585 in py::iter(&(v_A.clone()))? { let mut v_x: Py; 
v_x = t6585;
out.push(py::math1("log", &v_x.clone())?); } Py::list(out) }, { let mut out = vec![]; for t6586 in py::iter(&(v_S.clone()))? { let mut v_y: Py; 
v_y = t6586;
out.push(py::math1("log", &v_y.clone())?); } Py::list(out) })?).items()?;
    if t6587.len() != 3 { return Err(PyErr::Value("unpack".into())); }
    v_b = t6587[0].clone();
    v_a = t6587[1].clone();
    v_r2 = t6587[2].clone();
    return Ok(f__out(Py::dict(vec![(Py::str("c"), (if (v_no_exp.clone()).truthy() { v_a.clone() } else { py::math1("exp", &v_a.clone())? })), (Py::str("z"), v_b.clone()), (Py::str("r2"), v_r2.clone())]))?);
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn f_ss_fraction(mut v_seq: Py, mut v_order: Py) -> R {
    let mut v_f: Py = Py::None;
    let mut v_s: Py = Py::None;
    v_s = f__seq(v_seq.clone(), Py::str("seq"), Py::str("PROTEIN"))?;
    v_f = { let mut out = vec![]; for t6254 in py::iter(&(py::zip(&[v_order.clone(), Py::tuple(vec![c__HELIX()?, c__TURN()?, c__SHEET()?])])?))? { let mut v_g: Py; let mut v_k: Py; 
let t6255 = (t6254).items()?;
if t6255.len() != 2 { return Err(PyErr::Value("unpack".into())); }
v_k = t6255[0].clone();
v_g = t6255[1].clone();
out.push((v_k.clone(), py::truediv(&py::sum(&{ let mut out = vec![]; for t6256 in py::iter(&(v_g.clone()))? { let mut v_a: Py; 
v_a = t6256;
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
    v_S = py::sum(&{ let mut out = vec![]; for t6652 in py::iter(&(py::zip(&{ let mut v: Vec<Py> = vec![]; v.extend((v_rows.clone()).items()?); v })?))? { let mut v_col: Py; 
v_col = t6652;
if ({ let t6653 = py::len(&py::set(&v_col.clone())?)?; Py::Bool({ let t6654 = Py::Int(1i64); if py::gt(&t6653, &t6654)? { true } else { false } }) }).truthy() { out.push(Py::Int(1i64)); }} Py::list(out) }, None)?;
    let _ = f__need({ let t6655 = v_S.clone(); Py::Bool({ let t6656 = Py::Int(1i64); if py::ge(&t6655, &t6656)? { true } else { false } }) }, Py::str("Tajima's D needs at least one segregating site"))?;
    v_k = py::float(&f__pairwise_k(v_rows.clone())?)?;
    v_a1 = py::fsum(&{ let mut out = vec![]; for t6663 in py::range_iter(&[(Py::Int(1i64)).i()?, (v_n.clone()).i()?]) { let mut v_i: Py; 
v_i = t6663;
out.push(py::truediv(&Py::Float(1.0), &v_i.clone())?); } Py::list(out) })?;
    v_a2 = py::fsum(&{ let mut out = vec![]; for t6664 in py::range_iter(&[(Py::Int(1i64)).i()?, (v_n.clone()).i()?]) { let mut v_i: Py; 
v_i = t6664;
out.push(py::truediv(&Py::Float(1.0), &py::mul(&v_i.clone(), &v_i.clone())?)?); } Py::list(out) })?;
    v_b1 = py::truediv(&py::add(&v_n.clone(), &Py::Int(1i64))?, &py::mul(&Py::Int(3i64), &py::sub(&v_n.clone(), &Py::Int(1i64))?)?)?;
    v_b2 = py::truediv(&py::mul(&Py::Int(2i64), &py::add(&py::add(&py::mul(&v_n.clone(), &v_n.clone())?, &v_n.clone())?, &Py::Int(3i64))?)?, &py::mul(&py::mul(&Py::Int(9i64), &v_n.clone())?, &py::sub(&v_n.clone(), &Py::Int(1i64))?)?)?;
    v_c1 = py::sub(&v_b1.clone(), &py::truediv(&Py::Int(1i64), &v_a1.clone())?)?;
    v_c2 = py::add(&py::sub(&v_b2.clone(), &py::truediv(&py::add(&v_n.clone(), &Py::Int(2i64))?, &py::mul(&v_a1.clone(), &v_n.clone())?)?)?, &py::truediv(&v_a2.clone(), &py::mul(&v_a1.clone(), &v_a1.clone())?)?)?;
    let t6665 = (Py::tuple(vec![py::truediv(&v_c1.clone(), &v_a1.clone())?, py::truediv(&v_c2.clone(), &py::add(&py::mul(&v_a1.clone(), &v_a1.clone())?, &v_a2.clone())?)?])).items()?;
    if t6665.len() != 2 { return Err(PyErr::Value("unpack".into())); }
    v_e1 = t6665[0].clone();
    v_e2 = t6665[1].clone();
    v_var = py::add(&py::mul(&v_e1.clone(), &v_S.clone())?, &py::mul(&py::mul(&v_e2.clone(), &v_S.clone())?, &py::sub(&v_S.clone(), &Py::Int(1i64))?)?)?;
    return Ok(f__out(Py::dict(vec![(Py::str("D"), py::truediv(&py::sub(&v_k.clone(), &py::truediv(&v_S.clone(), &v_a1.clone())?)?, &(if (v_sqrt_denominator.clone()).truthy() { py::math1("sqrt", &v_var.clone())? } else { v_var.clone() }))?), (Py::str("k"), v_k.clone())]))?);
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn f_testcross_recombination(mut v_parental: Py, mut v_recombinant: Py) -> R {
    let mut v_n: Py = Py::None;
    let mut v_rf: Py = Py::None;
    let mut v_word: Py = Py::None;
    let _ = f__whole(Py::str("parental"), v_parental.clone(), Py::Int(0i64))?;
    let _ = f__whole(Py::str("recombinant"), v_recombinant.clone(), Py::Int(0i64))?;
    v_n = py::add(&v_parental.clone(), &v_recombinant.clone())?;
    let _ = f__need({ let t6729 = v_n.clone(); Py::Bool({ let t6730 = Py::Int(1i64); if py::ge(&t6729, &t6730)? { true } else { false } }) }, Py::str("at least one offspring is needed"))?;
    let _ = f__need({ let t6731 = v_recombinant.clone(); Py::Bool({ let t6732 = v_parental.clone(); if py::le(&t6731, &t6732)? { true } else { false } }) }, Py::str("the recombinant (less frequent) classes cannot outnumber the parental ones: recombination frequency is at most 0.5"))?;
    v_rf = py::truediv(&v_recombinant.clone(), &v_n.clone())?;
    v_word = (if ({ let t6733 = v_rf.clone(); Py::Bool({ let t6734 = Py::Float(0.5); if py::lt(&t6733, &t6734)? { true } else { false } }) }).truthy() { Py::str("linked") } else { Py::str("unlinked (independent assortment)") });
    return Ok(f__out(Py::dict(vec![(Py::str("RF"), v_rf.clone()), (Py::str("RF_percent"), py::mul(&Py::Float(100.0), &v_rf.clone())?), (Py::str("d_cM"), py::mul(&Py::Float(100.0), &v_rf.clone())?), (Py::str("SE"), py::math1("sqrt", &py::truediv(&py::mul(&v_rf.clone(), &py::sub(&Py::Float(1.0), &v_rf.clone())?)?, &v_n.clone())?)?), (Py::str("linkage"), v_word.clone())]))?);
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn f_three_point_testcross(mut v_classes: Py, mut v_counts: Py) -> R {
    let mut v_D: Py = Py::None;
    let mut v_P: Py = Py::None;
    let mut v_a: Py = Py::None;
    let mut v_b: Py = Py::None;
    let mut v_cnt: Py = Py::None;
    let mut v_coc: Py = Py::None;
    let mut v_dco: Py = Py::None;
    let mut v_diff: Py = Py::None;
    let mut v_k: Py = Py::None;
    let mut v_letters: Py = Py::None;
    let mut v_mid: Py = Py::None;
    let mut v_n: Py = Py::None;
    let mut v_names: Py = Py::None;
    let mut v_o: Py = Py::None;
    let mut v_outer: Py = Py::None;
    let mut v_rl: Py = Py::None;
    let mut v_ro: Py = Py::None;
    let mut v_rr: Py = Py::None;
    let t6832 = (f__tp_classes(v_classes.clone(), v_counts.clone())?).items()?;
    if t6832.len() != 4 { return Err(PyErr::Value("unpack".into())); }
    v_names = t6832[0].clone();
    v_cnt = t6832[1].clone();
    v_letters = t6832[2].clone();
    v_n = t6832[3].clone();
    v_o = f__rank_desc(v_cnt.clone())?;
    let _ = f__need({ let t6839 = py::getitem(&v_cnt.clone(), &py::getitem(&v_o.clone(), &Py::Int(1i64))?)?; Py::Bool({ let t6840 = py::getitem(&v_cnt.clone(), &py::getitem(&v_o.clone(), &Py::Int(2i64))?)?; if py::gt(&t6839, &t6840)? { true } else { false } }) }, Py::str("the two parental (most numerous) classes are not clearly the largest"))?;
    let _ = f__need({ let t6841 = py::getitem(&v_cnt.clone(), &py::getitem(&v_o.clone(), &Py::Int(6i64))?)?; Py::Bool({ let t6842 = py::getitem(&v_cnt.clone(), &py::getitem(&v_o.clone(), &Py::Int(5i64))?)?; if py::lt(&t6841, &t6842)? { true } else { false } }) }, Py::str("the two double-crossover (least numerous) classes are not clearly the smallest"))?;
    v_P = py::getitem(&v_names.clone(), &py::getitem(&v_o.clone(), &Py::Int(0i64))?)?;
    let _ = f__need({ let t6843 = py::getitem(&v_names.clone(), &py::getitem(&v_o.clone(), &Py::Int(1i64))?)?; Py::Bool({ let t6844 = f__flip_case(v_P.clone())?; if py::eq(&t6843, &t6844) { true } else { false } }) }, Py::str("the two most numerous classes must be reciprocal (the parental classes)"))?;
    v_D = py::getitem(&v_names.clone(), &py::getitem(&v_o.clone(), &Py::Int(7i64))?)?;
    let _ = f__need({ let t6846 = py::getitem(&v_names.clone(), &py::getitem(&v_o.clone(), &Py::Int(6i64))?)?; Py::Bool({ let t6847 = f__flip_case(v_D.clone())?; if py::eq(&t6846, &t6847) { true } else { false } }) }, Py::str("the two least numerous classes must be reciprocal (the double crossovers)"))?;
    v_diff = Py::Int(0i64);
    for t6848 in py::range_iter(&[(Py::Int(3i64)).i()?]) {
        v_k = t6848;
        if ({ let t6849 = py::getitem(&v_D.clone(), &v_k.clone())?; Py::Bool({ let t6850 = py::getitem(&v_P.clone(), &v_k.clone())?; if py::ne(&t6849, &t6850) { true } else { false } }) }).truthy() {
            v_diff = py::iadd(&v_diff, &Py::Int(1i64))?;
        }
    }
    v_mid = py::neg(&Py::Int(1i64))?;
    for t6851 in py::range_iter(&[(Py::Int(3i64)).i()?]) {
        v_k = t6851;
        if ({ let t6862 = { let t6856 = { let t6852 = v_diff.clone(); Py::Bool({ let t6853 = Py::Int(1i64); if py::eq(&t6852, &t6853) { true } else { false } }) }; if !t6856.truthy() { t6856 } else { { let t6854 = py::getitem(&v_D.clone(), &v_k.clone())?; Py::Bool({ let t6855 = py::getitem(&v_P.clone(), &v_k.clone())?; if py::ne(&t6854, &t6855) { true } else { false } }) } } }; if t6862.truthy() { t6862 } else { { let t6861 = { let t6857 = v_diff.clone(); Py::Bool({ let t6858 = Py::Int(2i64); if py::eq(&t6857, &t6858) { true } else { false } }) }; if !t6861.truthy() { t6861 } else { { let t6859 = py::getitem(&v_D.clone(), &v_k.clone())?; Py::Bool({ let t6860 = py::getitem(&v_P.clone(), &v_k.clone())?; if py::eq(&t6859, &t6860) { true } else { false } }) } } } } }).truthy() {
            v_mid = v_k.clone();
        }
    }
    v_outer = Py::list(vec![]);
    for t6863 in py::range_iter(&[(Py::Int(3i64)).i()?]) {
        v_k = t6863;
        if ({ let t6864 = v_k.clone(); Py::Bool({ let t6865 = v_mid.clone(); if py::ne(&t6864, &t6865) { true } else { false } }) }).truthy() {
            let _ = py::method(&v_outer.clone(), "append", &[v_k.clone()])?;
        }
    }
    let t6866 = (Py::tuple(vec![py::getitem(&v_outer.clone(), &Py::Int(0i64))?, py::getitem(&v_outer.clone(), &Py::Int(1i64))?])).items()?;
    if t6866.len() != 2 { return Err(PyErr::Value("unpack".into())); }
    v_a = t6866[0].clone();
    v_b = t6866[1].clone();
    v_rl = py::truediv(&f__tp_rec(v_names.clone(), v_cnt.clone(), v_P.clone(), v_a.clone(), v_mid.clone())?, &v_n.clone())?;
    v_rr = py::truediv(&f__tp_rec(v_names.clone(), v_cnt.clone(), v_P.clone(), v_mid.clone(), v_b.clone())?, &v_n.clone())?;
    v_ro = py::truediv(&f__tp_rec(v_names.clone(), v_cnt.clone(), v_P.clone(), v_a.clone(), v_b.clone())?, &v_n.clone())?;
    let _ = f__need({ let t6878 = { let t6874 = v_rl.clone(); Py::Bool({ let t6875 = Py::Int(0i64); if py::gt(&t6874, &t6875)? { true } else { false } }) }; if !t6878.truthy() { t6878 } else { { let t6876 = v_rr.clone(); Py::Bool({ let t6877 = Py::Int(0i64); if py::gt(&t6876, &t6877)? { true } else { false } }) } } }, Py::str("each interval needs some recombinants"))?;
    v_dco = py::truediv(&py::add(&py::getitem(&v_cnt.clone(), &py::getitem(&v_o.clone(), &Py::Int(6i64))?)?, &py::getitem(&v_cnt.clone(), &py::getitem(&v_o.clone(), &Py::Int(7i64))?)?)?, &v_n.clone())?;
    v_coc = py::truediv(&v_dco.clone(), &py::mul(&v_rl.clone(), &v_rr.clone())?)?;
    return Ok(f__out(Py::dict(vec![(Py::str("gene_order"), py::add(&py::add(&py::add(&py::add(&py::getitem(&v_letters.clone(), &v_a.clone())?, &Py::str("-"))?, &py::getitem(&v_letters.clone(), &v_mid.clone())?)?, &Py::str("-"))?, &py::getitem(&v_letters.clone(), &v_b.clone())?)?), (Py::str("middle_gene"), py::getitem(&v_letters.clone(), &v_mid.clone())?), (Py::str("parental_classes"), py::add(&py::add(&py::getitem(&v_names.clone(), &py::getitem(&v_o.clone(), &Py::Int(0i64))?)?, &Py::str(" "))?, &py::getitem(&v_names.clone(), &py::getitem(&v_o.clone(), &Py::Int(1i64))?)?)?), (Py::str("dco_classes"), py::add(&py::add(&py::getitem(&v_names.clone(), &py::getitem(&v_o.clone(), &Py::Int(6i64))?)?, &Py::str(" "))?, &py::getitem(&v_names.clone(), &py::getitem(&v_o.clone(), &Py::Int(7i64))?)?)?), (Py::str("RF_left"), v_rl.clone()), (Py::str("RF_right"), v_rr.clone()), (Py::str("d_left_cM"), py::mul(&Py::Float(100.0), &v_rl.clone())?), (Py::str("d_right_cM"), py::mul(&Py::Float(100.0), &v_rr.clone())?), (Py::str("d_total_cM"), py::mul(&Py::Float(100.0), &py::add(&v_rl.clone(), &v_rr.clone())?)?), (Py::str("RF_outer_direct"), v_ro.clone()), (Py::str("coincidence"), v_coc.clone()), (Py::str("interference"), py::sub(&Py::Float(1.0), &v_coc.clone())?)]))?);
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn f_tm_nn(mut v_seq: Py, mut v_Na: Py, mut v_dnac1: Py, mut v_dnac2: Py) -> R {
    let mut v_s: Py = Py::None;
    v_s = f__seq(v_seq.clone(), Py::str("seq"), Py::str("DNA"))?;
    let _ = f__need({ let t6957 = py::len(&v_s.clone())?; Py::Bool({ let t6958 = Py::Int(2i64); if py::ge(&t6957, &t6958)? { true } else { false } }) }, Py::str("the sequence needs at least 2 bases"))?;
    let _ = f__pos(Py::str("Na"), Py::tuple(vec![v_Na.clone()]))?;
    let _ = f__pos(Py::str("dnac1"), Py::tuple(vec![v_dnac1.clone()]))?;
    let _ = f__fin(Py::str("dnac2"), Py::tuple(vec![v_dnac2.clone()]))?;
    let _ = f__need({ let t6968 = v_dnac2.clone(); Py::Bool({ let t6969 = Py::Int(0i64); if py::ge(&t6968, &t6969)? { true } else { false } }) }, Py::str("dnac2 must be zero or positive"))?;
    let _ = f__need({ let t6974 = { let t6970 = v_s.clone(); Py::Bool({ let t6971 = f__revcomp(v_s.clone())?; if py::eq(&t6970, &t6971) { true } else { false } }) }; if t6974.truthy() { t6974 } else { { let t6972 = v_dnac1.clone(); Py::Bool({ let t6973 = py::truediv(&v_dnac2.clone(), &Py::Int(2i64))?; if py::gt(&t6972, &t6973)? { true } else { false } }) } } }, Py::str("dnac1 must exceed dnac2/2 (k = c1 - c2/2 must be positive)"))?;
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
    let _ = f__need({ let t7128 = v_strand.clone(); Py::Bool({ let t7129 = Py::tuple(vec![Py::str("coding"), Py::str("template")]); if py::contains(&t7129, &t7128)? { true } else { false } }) }, Py::str("strand must be coding or template"))?;
    v_m = (if ({ let t7130 = v_strand.clone(); Py::Bool({ let t7131 = Py::str("coding"); if py::eq(&t7130, &t7131) { true } else { false } }) }).truthy() { v_s.clone() } else { f__revcomp(v_s.clone())? });
    return Ok(Py::dict(vec![(Py::str("mrna"), py::method(&v_m.clone(), "replace", &[Py::str("T"), Py::str("U")])?)]));
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn f_translation(mut v_seq: Py, mut v_stops: Py) -> R {
    let mut v_p: Py = Py::None;
    let mut v_s: Py = Py::None;
    v_s = py::method(&f__seq(v_seq.clone(), Py::str("seq"), Py::str("NUC"))?, "replace", &[Py::str("U"), Py::str("T")])?;
    let _ = f__need({ let t7196 = py::modulo(&py::len(&v_s.clone())?, &Py::Int(3i64))?; Py::Bool({ let t7197 = Py::Int(0i64); if py::eq(&t7196, &t7197) { true } else { false } }) }, Py::str("the sequence length must be a multiple of 3 (whole codons read in frame 1)"))?;
    let _ = f__need({ let t7198 = v_stops.clone(); Py::Bool({ let t7199 = Py::tuple(vec![Py::str("show"), Py::str("truncate")]); if py::contains(&t7199, &t7198)? { true } else { false } }) }, Py::str("stops must be show or truncate"))?;
    v_p = f__translate(v_s.clone(), c_CODON()?, { let t7204 = v_stops.clone(); Py::Bool({ let t7205 = Py::str("truncate"); if py::eq(&t7204, &t7205) { true } else { false } }) })?;
    return Ok(Py::dict(vec![(Py::str("protein"), v_p.clone()), (Py::str("length"), py::len(&v_p.clone())?)]));
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn f_tree_count(mut v_n: Py, mut v_rooted: Py) -> R {
    let mut v_lo: Py = Py::None;
    let _ = f__need({ let t7247 = v_rooted.clone(); Py::Bool({ let t7248 = Py::tuple(vec![Py::str("rooted"), Py::str("unrooted")]); if py::contains(&t7248, &t7247)? { true } else { false } }) }, Py::str("rooted must be rooted or unrooted"))?;
    v_lo = (if ({ let t7249 = v_rooted.clone(); Py::Bool({ let t7250 = Py::str("rooted"); if py::eq(&t7249, &t7250) { true } else { false } }) }).truthy() { Py::Int(2i64) } else { Py::Int(3i64) });
    let _ = f__need({ let t7254 = Py::Bool(matches!(v_n.clone(), Py::Int(_) | Py::Big(_) | Py::Bool(_))); if !t7254.truthy() { t7254 } else { { let t7251 = v_lo.clone(); Py::Bool({ let t7252 = v_n.clone(); if py::le(&t7251, &t7252)? { { let t7253 = Py::Int(100i64); if py::le(&t7252, &t7253)? { true } else { false } } } else { false } }) } } }, Py::str(&["n must be a whole number of taxa from ".to_string(), py::format(&v_lo.clone(), "")?, " to 100".to_string()].concat()))?;
    return Ok(Py::dict(vec![(Py::str("count"), (if ({ let t7255 = v_rooted.clone(); Py::Bool({ let t7256 = Py::str("rooted"); if py::eq(&t7255, &t7256) { true } else { false } }) }).truthy() { f__dfact(py::sub(&py::mul(&Py::Int(2i64), &v_n.clone())?, &Py::Int(3i64))?)? } else { f__dfact(py::sub(&py::mul(&Py::Int(2i64), &v_n.clone())?, &Py::Int(5i64))?)? }))]));
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
    v_size = { let mut out = vec![]; for t7317 in py::range_iter(&[(v_n.clone()).i()?]) { let mut v_i: Py; 
v_i = t7317;
out.push((v_i.clone(), Py::Int(1i64))); } Py::dict(out) };
    v_d = { let mut out = vec![]; for t7318 in py::range_iter(&[(v_n.clone()).i()?]) { let mut v_i: Py; 
v_i = t7318;
for t7319 in py::range_iter(&[(py::add(&v_i.clone(), &Py::Int(1i64))?).i()?, (v_n.clone()).i()?]) { let mut v_j: Py; 
v_j = t7319;
out.push((Py::tuple(vec![v_i.clone(), v_j.clone()]), py::getitem(&py::getitem(&v_M.clone(), &v_i.clone())?, &v_j.clone())?)); }} Py::dict(out) };
    v_scale = py::minmax(&[py::method(&v_d.clone(), "values", &[])?], true)?;
    let t7320 = (Py::tuple(vec![Py::list(vec![]), v_n.clone()])).items()?;
    if t7320.len() != 2 { return Err(PyErr::Value("unpack".into())); }
    v_heights = t7320[0].clone();
    v_nxt = t7320[1].clone();
    while ({ let t7344 = py::len(&v_size.clone())?; Py::Bool({ let t7345 = Py::Int(1i64); if py::gt(&t7344, &t7345)? { true } else { false } }) }).truthy() {
        v_best = py::minmax(&[py::method(&v_d.clone(), "values", &[])?], false)?;
        let _ = f__need(py::not(&f__tie(Py::list((py::method(&v_d.clone(), "values", &[])?).items()?), v_best.clone(), v_scale.clone())?), Py::str("two pairs tie at the closest distance: the tree is not unique"))?;
        let t7330 = (py::minmax(&[{ let mut out = vec![]; for t7326 in py::iter(&(py::method(&v_d.clone(), "items", &[])?))? { let mut v_k: Py; let mut v_v: Py; 
let t7327 = (t7326).items()?;
if t7327.len() != 2 { return Err(PyErr::Value("unpack".into())); }
v_k = t7327[0].clone();
v_v = t7327[1].clone();
if ({ let t7328 = v_v.clone(); Py::Bool({ let t7329 = v_best.clone(); if py::eq(&t7328, &t7329) { true } else { false } }) }).truthy() { out.push(v_k.clone()); }} Py::list(out) }], false)?).items()?;
        if t7330.len() != 2 { return Err(PyErr::Value("unpack".into())); }
        v_i = t7330[0].clone();
        v_j = t7330[1].clone();
        let _ = py::method(&v_heights.clone(), "append", &[py::truediv(&v_best.clone(), &Py::Int(2i64))?])?;
        v_new = v_nxt.clone();
        v_nxt = py::iadd(&v_nxt, &Py::Int(1i64))?;
        for t7331 in py::iter(&(Py::list((v_size.clone()).items()?)))? {
            v_k = t7331;
            if ({ let t7332 = v_k.clone(); Py::Bool({ let t7333 = Py::tuple(vec![v_i.clone(), v_j.clone()]); if py::contains(&t7333, &t7332)? { true } else { false } }) }).truthy() {
                continue;
            }
            let t7334 = (Py::tuple(vec![py::getitem(&v_d.clone(), &Py::tuple((py::sorted(&Py::tuple(vec![v_i.clone(), v_k.clone()]), (Py::Bool(false)).truthy())?).items()?))?, py::getitem(&v_d.clone(), &Py::tuple((py::sorted(&Py::tuple(vec![v_j.clone(), v_k.clone()]), (Py::Bool(false)).truthy())?).items()?))?])).items()?;
            if t7334.len() != 2 { return Err(PyErr::Value("unpack".into())); }
            v_dik = t7334[0].clone();
            v_djk = t7334[1].clone();
            if ({ let t7335 = v_linkage.clone(); Py::Bool({ let t7336 = Py::str("average"); if py::eq(&t7335, &t7336) { true } else { false } }) }).truthy() {
                py::setitem(&v_d.clone(), &Py::tuple(vec![v_k.clone(), v_new.clone()]), py::truediv(&py::add(&py::mul(&py::getitem(&v_size.clone(), &v_i.clone())?, &v_dik.clone())?, &py::mul(&py::getitem(&v_size.clone(), &v_j.clone())?, &v_djk.clone())?)?, &py::add(&py::getitem(&v_size.clone(), &v_i.clone())?, &py::getitem(&v_size.clone(), &v_j.clone())?)?)?)?;
            } else {
                py::setitem(&v_d.clone(), &Py::tuple(vec![v_k.clone(), v_new.clone()]), py::minmax(&[v_dik.clone(), v_djk.clone()], false)?)?;
            }
        }
        py::setitem(&v_size.clone(), &v_new.clone(), py::add(&py::method(&v_size.clone(), "pop", &[v_i.clone()])?, &py::method(&v_size.clone(), "pop", &[v_j.clone()])?)?)?;
        v_d = { let mut out = vec![]; for t7337 in py::iter(&(py::method(&v_d.clone(), "items", &[])?))? { let mut v_k: Py; let mut v_v: Py; 
let t7338 = (t7337).items()?;
if t7338.len() != 2 { return Err(PyErr::Value("unpack".into())); }
v_k = t7338[0].clone();
v_v = t7338[1].clone();
if ({ let t7343 = { let t7339 = v_i.clone(); Py::Bool({ let t7340 = v_k.clone(); if !py::contains(&t7340, &t7339)? { true } else { false } }) }; if !t7343.truthy() { t7343 } else { { let t7341 = v_j.clone(); Py::Bool({ let t7342 = v_k.clone(); if !py::contains(&t7342, &t7341)? { true } else { false } }) } } }).truthy() { out.push((v_k.clone(), v_v.clone())); }} Py::dict(out) };
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
    let t7410 = (Py::tuple(vec![py::len(&v_rows.clone())?, py::len(&py::getitem(&v_rows.clone(), &Py::Int(0i64))?)?])).items()?;
    if t7410.len() != 2 { return Err(PyErr::Value("unpack".into())); }
    v_n = t7410[0].clone();
    v_L = t7410[1].clone();
    v_S = py::sum(&{ let mut out = vec![]; for t7411 in py::iter(&(py::zip(&{ let mut v: Vec<Py> = vec![]; v.extend((v_rows.clone()).items()?); v })?))? { let mut v_col: Py; 
v_col = t7411;
if ({ let t7412 = py::len(&py::set(&v_col.clone())?)?; Py::Bool({ let t7413 = Py::Int(1i64); if py::gt(&t7412, &t7413)? { true } else { false } }) }).truthy() { out.push(Py::Int(1i64)); }} Py::list(out) }, None)?;
    v_a = py::fsum(&{ let mut out = vec![]; for t7414 in py::range_iter(&[(Py::Int(1i64)).i()?, (py::add(&v_n.clone(), &(if (v_harmonic_to_n.clone()).truthy() { Py::Int(1i64) } else { Py::Int(0i64) }))?).i()?]) { let mut v_i: Py; 
v_i = t7414;
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
fn s_2496(mut v_inputs: Py, mut v_outputs: Py) -> R {
    let mut v_k: Py = Py::None;
    let mut v_v: Py = Py::None;
    let mut v_want: Py = Py::None;
    let mut v_worst: Py = Py::None;
    let t2498: Result<Option<Py>, PyErr> = (|| -> Result<Option<Py>, PyErr> {
            v_want = { let t2497 = v_inputs.clone(); f_hardy_weinberg_frequencies(py::method(&t2497, "get", &[Py::str("p"), Py::None])?, py::method(&t2497, "get", &[Py::str("q"), Py::None])?, py::method(&t2497, "get", &[Py::str("q2"), Py::None])?, py::method(&t2497, "get", &[Py::str("N"), Py::None])?)? };
            #[allow(unreachable_code)]
            Ok(None)
        })();
    match t2498 {
        Ok(Some(r)) => { return Ok(r); }
        Ok(None) => {}
        Err(e) if (!matches!(e, PyErr::Unsupported(_))) => {
            return Ok(Py::Float(f64::INFINITY));
        }
        Err(e) => return Err(e),
    }
    v_worst = Py::Float(0.0);
    for t2499 in py::iter(&(py::method(&v_want.clone(), "items", &[])?))? {
        let t2536 = (t2499).items()?;
        if t2536.len() != 2 { return Err(PyErr::Value("unpack".into())); }
        v_k = t2536[0].clone();
        v_v = t2536[1].clone();
        if ({ let t2500 = v_k.clone(); Py::Bool({ let t2501 = v_outputs.clone(); if !py::contains(&t2501, &t2500)? { true } else { false } }) }).truthy() {
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
fn s_2564(mut v_inputs: Py, mut v_outputs: Py) -> R {
    let mut v_k: Py = Py::None;
    let mut v_v: Py = Py::None;
    let mut v_want: Py = Py::None;
    let mut v_worst: Py = Py::None;
    let t2566: Result<Option<Py>, PyErr> = (|| -> Result<Option<Py>, PyErr> {
            v_want = { let t2565 = v_inputs.clone(); f_hardy_weinberg_test(py::getitem(&t2565, &Py::str("n_AA"))?, py::getitem(&t2565, &Py::str("n_Aa"))?, py::getitem(&t2565, &Py::str("n_aa"))?)? };
            #[allow(unreachable_code)]
            Ok(None)
        })();
    match t2566 {
        Ok(Some(r)) => { return Ok(r); }
        Ok(None) => {}
        Err(e) if (!matches!(e, PyErr::Unsupported(_))) => {
            return Ok(Py::Float(f64::INFINITY));
        }
        Err(e) => return Err(e),
    }
    v_worst = Py::Float(0.0);
    for t2567 in py::iter(&(py::method(&v_want.clone(), "items", &[])?))? {
        let t2604 = (t2567).items()?;
        if t2604.len() != 2 { return Err(PyErr::Value("unpack".into())); }
        v_k = t2604[0].clone();
        v_v = t2604[1].clone();
        if ({ let t2568 = v_k.clone(); Py::Bool({ let t2569 = v_outputs.clone(); if !py::contains(&t2569, &t2568)? { true } else { false } }) }).truthy() {
            return Ok(Py::Float(f64::INFINITY));
        }
        v_worst = py::minmax(&[v_worst.clone(), f__rel(py::getitem(&v_outputs.clone(), &v_k.clone())?, v_v.clone())?], true)?;
    }
    return Ok(v_worst.clone());
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn s_2630(mut v_inputs: Py, mut v_outputs: Py) -> R {
    let mut v_k: Py = Py::None;
    let mut v_v: Py = Py::None;
    let mut v_want: Py = Py::None;
    let mut v_worst: Py = Py::None;
    let t2632: Result<Option<Py>, PyErr> = (|| -> Result<Option<Py>, PyErr> {
            v_want = { let t2631 = v_inputs.clone(); f_het_advantage(py::getitem(&t2631, &Py::str("w_AA"))?, py::getitem(&t2631, &Py::str("w_Aa"))?, py::getitem(&t2631, &Py::str("w_aa"))?)? };
            #[allow(unreachable_code)]
            Ok(None)
        })();
    match t2632 {
        Ok(Some(r)) => { return Ok(r); }
        Ok(None) => {}
        Err(e) if (!matches!(e, PyErr::Unsupported(_))) => {
            return Ok(Py::Float(f64::INFINITY));
        }
        Err(e) => return Err(e),
    }
    v_worst = Py::Float(0.0);
    for t2633 in py::iter(&(py::method(&v_want.clone(), "items", &[])?))? {
        let t2670 = (t2633).items()?;
        if t2670.len() != 2 { return Err(PyErr::Value("unpack".into())); }
        v_k = t2670[0].clone();
        v_v = t2670[1].clone();
        if ({ let t2634 = v_k.clone(); Py::Bool({ let t2635 = v_outputs.clone(); if !py::contains(&t2635, &t2634)? { true } else { false } }) }).truthy() {
            return Ok(Py::Float(f64::INFINITY));
        }
        v_worst = py::minmax(&[v_worst.clone(), f__rel(py::getitem(&v_outputs.clone(), &v_k.clone())?, v_v.clone())?], true)?;
    }
    return Ok(v_worst.clone());
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn s_2846(mut v_inputs: Py, mut v_outputs: Py) -> R {
    let mut v_k: Py = Py::None;
    let mut v_v: Py = Py::None;
    let mut v_want: Py = Py::None;
    let mut v_worst: Py = Py::None;
    let t2848: Result<Option<Py>, PyErr> = (|| -> Result<Option<Py>, PyErr> {
            v_want = { let t2847 = v_inputs.clone(); f_hill_fit(py::getitem(&t2847, &Py::str("S"))?, py::getitem(&t2847, &Py::str("v"))?, py::method(&t2847, "get", &[Py::str("fixed_n"), Py::Bool(false)])?)? };
            #[allow(unreachable_code)]
            Ok(None)
        })();
    match t2848 {
        Ok(Some(r)) => { return Ok(r); }
        Ok(None) => {}
        Err(e) if (!matches!(e, PyErr::Unsupported(_))) => {
            return Ok(Py::Float(f64::INFINITY));
        }
        Err(e) => return Err(e),
    }
    v_worst = Py::Float(0.0);
    for t2849 in py::iter(&(py::method(&v_want.clone(), "items", &[])?))? {
        let t2886 = (t2849).items()?;
        if t2886.len() != 2 { return Err(PyErr::Value("unpack".into())); }
        v_k = t2886[0].clone();
        v_v = t2886[1].clone();
        if ({ let t2850 = v_k.clone(); Py::Bool({ let t2851 = v_outputs.clone(); if !py::contains(&t2851, &t2850)? { true } else { false } }) }).truthy() {
            return Ok(Py::Float(f64::INFINITY));
        }
        v_worst = py::minmax(&[v_worst.clone(), f__rel(py::getitem(&v_outputs.clone(), &v_k.clone())?, v_v.clone())?], true)?;
    }
    return Ok(v_worst.clone());
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn s_2942(mut v_inputs: Py, mut v_outputs: Py) -> R {
    let mut v_k: Py = Py::None;
    let mut v_v: Py = Py::None;
    let mut v_want: Py = Py::None;
    let mut v_worst: Py = Py::None;
    let t2944: Result<Option<Py>, PyErr> = (|| -> Result<Option<Py>, PyErr> {
            v_want = { let t2943 = v_inputs.clone(); f_hill_plot(py::getitem(&t2943, &Py::str("S"))?, py::getitem(&t2943, &Py::str("Y"))?)? };
            #[allow(unreachable_code)]
            Ok(None)
        })();
    match t2944 {
        Ok(Some(r)) => { return Ok(r); }
        Ok(None) => {}
        Err(e) if (!matches!(e, PyErr::Unsupported(_))) => {
            return Ok(Py::Float(f64::INFINITY));
        }
        Err(e) => return Err(e),
    }
    v_worst = Py::Float(0.0);
    for t2945 in py::iter(&(py::method(&v_want.clone(), "items", &[])?))? {
        let t2982 = (t2945).items()?;
        if t2982.len() != 2 { return Err(PyErr::Value("unpack".into())); }
        v_k = t2982[0].clone();
        v_v = t2982[1].clone();
        if ({ let t2946 = v_k.clone(); Py::Bool({ let t2947 = v_outputs.clone(); if !py::contains(&t2947, &t2946)? { true } else { false } }) }).truthy() {
            return Ok(Py::Float(f64::INFINITY));
        }
        v_worst = py::minmax(&[v_worst.clone(), f__rel(py::getitem(&v_outputs.clone(), &v_k.clone())?, v_v.clone())?], true)?;
    }
    return Ok(v_worst.clone());
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn s_3023(mut v_inputs: Py, mut v_outputs: Py) -> R {
    let mut v_k: Py = Py::None;
    let mut v_v: Py = Py::None;
    let mut v_want: Py = Py::None;
    let mut v_worst: Py = Py::None;
    let t3025: Result<Option<Py>, PyErr> = (|| -> Result<Option<Py>, PyErr> {
            v_want = { let t3024 = v_inputs.clone(); f_instability(py::getitem(&t3024, &Py::str("seq"))?)? };
            #[allow(unreachable_code)]
            Ok(None)
        })();
    match t3025 {
        Ok(Some(r)) => { return Ok(r); }
        Ok(None) => {}
        Err(e) if (!matches!(e, PyErr::Unsupported(_))) => {
            return Ok(Py::Float(f64::INFINITY));
        }
        Err(e) => return Err(e),
    }
    v_worst = Py::Float(0.0);
    for t3026 in py::iter(&(py::method(&v_want.clone(), "items", &[])?))? {
        let t3063 = (t3026).items()?;
        if t3063.len() != 2 { return Err(PyErr::Value("unpack".into())); }
        v_k = t3063[0].clone();
        v_v = t3063[1].clone();
        if ({ let t3027 = v_k.clone(); Py::Bool({ let t3028 = v_outputs.clone(); if !py::contains(&t3028, &t3027)? { true } else { false } }) }).truthy() {
            return Ok(Py::Float(f64::INFINITY));
        }
        v_worst = py::minmax(&[v_worst.clone(), f__rel(py::getitem(&v_outputs.clone(), &v_k.clone())?, v_v.clone())?], true)?;
    }
    return Ok(v_worst.clone());
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn s_3138(mut v_inputs: Py, mut v_outputs: Py) -> R {
    let mut v_k: Py = Py::None;
    let mut v_v: Py = Py::None;
    let mut v_want: Py = Py::None;
    let mut v_worst: Py = Py::None;
    let t3140: Result<Option<Py>, PyErr> = (|| -> Result<Option<Py>, PyErr> {
            v_want = { let t3139 = v_inputs.clone(); f_isoelectric_point(py::getitem(&t3139, &Py::str("seq"))?)? };
            #[allow(unreachable_code)]
            Ok(None)
        })();
    match t3140 {
        Ok(Some(r)) => { return Ok(r); }
        Ok(None) => {}
        Err(e) if (!matches!(e, PyErr::Unsupported(_))) => {
            return Ok(Py::Float(f64::INFINITY));
        }
        Err(e) => return Err(e),
    }
    v_worst = Py::Float(0.0);
    for t3141 in py::iter(&(py::method(&v_want.clone(), "items", &[])?))? {
        let t3178 = (t3141).items()?;
        if t3178.len() != 2 { return Err(PyErr::Value("unpack".into())); }
        v_k = t3178[0].clone();
        v_v = t3178[1].clone();
        if ({ let t3142 = v_k.clone(); Py::Bool({ let t3143 = v_outputs.clone(); if !py::contains(&t3143, &t3142)? { true } else { false } }) }).truthy() {
            return Ok(Py::Float(f64::INFINITY));
        }
        v_worst = py::minmax(&[v_worst.clone(), f__rel(py::getitem(&v_outputs.clone(), &v_k.clone())?, v_v.clone())?], true)?;
    }
    return Ok(v_worst.clone());
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn s_3226(mut v_inputs: Py, mut v_outputs: Py) -> R {
    let mut v_k: Py = Py::None;
    let mut v_v: Py = Py::None;
    let mut v_want: Py = Py::None;
    let mut v_worst: Py = Py::None;
    let t3228: Result<Option<Py>, PyErr> = (|| -> Result<Option<Py>, PyErr> {
            v_want = { let t3227 = v_inputs.clone(); f_jukes_cantor(py::getitem(&t3227, &Py::str("seq1"))?, py::getitem(&t3227, &Py::str("seq2"))?)? };
            #[allow(unreachable_code)]
            Ok(None)
        })();
    match t3228 {
        Ok(Some(r)) => { return Ok(r); }
        Ok(None) => {}
        Err(e) if (!matches!(e, PyErr::Unsupported(_))) => {
            return Ok(Py::Float(f64::INFINITY));
        }
        Err(e) => return Err(e),
    }
    v_worst = Py::Float(0.0);
    for t3229 in py::iter(&(py::method(&v_want.clone(), "items", &[])?))? {
        let t3266 = (t3229).items()?;
        if t3266.len() != 2 { return Err(PyErr::Value("unpack".into())); }
        v_k = t3266[0].clone();
        v_v = t3266[1].clone();
        if ({ let t3230 = v_k.clone(); Py::Bool({ let t3231 = v_outputs.clone(); if !py::contains(&t3231, &t3230)? { true } else { false } }) }).truthy() {
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
fn s_3333(mut v_inputs: Py, mut v_outputs: Py) -> R {
    let mut v_k: Py = Py::None;
    let mut v_v: Py = Py::None;
    let mut v_want: Py = Py::None;
    let mut v_worst: Py = Py::None;
    let t3335: Result<Option<Py>, PyErr> = (|| -> Result<Option<Py>, PyErr> {
            v_want = { let t3334 = v_inputs.clone(); f_kimura(py::getitem(&t3334, &Py::str("seq1"))?, py::getitem(&t3334, &Py::str("seq2"))?)? };
            #[allow(unreachable_code)]
            Ok(None)
        })();
    match t3335 {
        Ok(Some(r)) => { return Ok(r); }
        Ok(None) => {}
        Err(e) if (!matches!(e, PyErr::Unsupported(_))) => {
            return Ok(Py::Float(f64::INFINITY));
        }
        Err(e) => return Err(e),
    }
    v_worst = Py::Float(0.0);
    for t3336 in py::iter(&(py::method(&v_want.clone(), "items", &[])?))? {
        let t3373 = (t3336).items()?;
        if t3373.len() != 2 { return Err(PyErr::Value("unpack".into())); }
        v_k = t3373[0].clone();
        v_v = t3373[1].clone();
        if ({ let t3337 = v_k.clone(); Py::Bool({ let t3338 = v_outputs.clone(); if !py::contains(&t3338, &t3337)? { true } else { false } }) }).truthy() {
            return Ok(Py::Float(f64::INFINITY));
        }
        v_worst = py::minmax(&[v_worst.clone(), f__rel(py::getitem(&v_outputs.clone(), &v_k.clone())?, v_v.clone())?], true)?;
    }
    return Ok(v_worst.clone());
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn s_3406(mut v_inputs: Py, mut v_outputs: Py) -> R {
    let mut v_k: Py = Py::None;
    let mut v_v: Py = Py::None;
    let mut v_want: Py = Py::None;
    let mut v_worst: Py = Py::None;
    let t3408: Result<Option<Py>, PyErr> = (|| -> Result<Option<Py>, PyErr> {
            v_want = { let t3407 = v_inputs.clone(); f_kmer_counts(py::getitem(&t3407, &Py::str("seq"))?, py::getitem(&t3407, &Py::str("k"))?)? };
            #[allow(unreachable_code)]
            Ok(None)
        })();
    match t3408 {
        Ok(Some(r)) => { return Ok(r); }
        Ok(None) => {}
        Err(e) if (!matches!(e, PyErr::Unsupported(_))) => {
            return Ok(Py::Float(f64::INFINITY));
        }
        Err(e) => return Err(e),
    }
    v_worst = Py::Float(0.0);
    for t3409 in py::iter(&(py::method(&v_want.clone(), "items", &[])?))? {
        let t3446 = (t3409).items()?;
        if t3446.len() != 2 { return Err(PyErr::Value("unpack".into())); }
        v_k = t3446[0].clone();
        v_v = t3446[1].clone();
        if ({ let t3410 = v_k.clone(); Py::Bool({ let t3411 = v_outputs.clone(); if !py::contains(&t3411, &t3410)? { true } else { false } }) }).truthy() {
            return Ok(Py::Float(f64::INFINITY));
        }
        v_worst = py::minmax(&[v_worst.clone(), f__rel(py::getitem(&v_outputs.clone(), &v_k.clone())?, v_v.clone())?], true)?;
    }
    return Ok(v_worst.clone());
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn s_3520(mut v_inputs: Py, mut v_outputs: Py) -> R {
    let mut v_k: Py = Py::None;
    let mut v_v: Py = Py::None;
    let mut v_want: Py = Py::None;
    let mut v_worst: Py = Py::None;
    let t3522: Result<Option<Py>, PyErr> = (|| -> Result<Option<Py>, PyErr> {
            v_want = { let t3521 = v_inputs.clone(); f_leslie_growth(py::getitem(&t3521, &Py::str("F"))?, py::getitem(&t3521, &Py::str("S"))?, py::method(&t3521, "get", &[Py::str("use_r0"), Py::Bool(false)])?)? };
            #[allow(unreachable_code)]
            Ok(None)
        })();
    match t3522 {
        Ok(Some(r)) => { return Ok(r); }
        Ok(None) => {}
        Err(e) if (!matches!(e, PyErr::Unsupported(_))) => {
            return Ok(Py::Float(f64::INFINITY));
        }
        Err(e) => return Err(e),
    }
    v_worst = Py::Float(0.0);
    for t3523 in py::iter(&(py::method(&v_want.clone(), "items", &[])?))? {
        let t3560 = (t3523).items()?;
        if t3560.len() != 2 { return Err(PyErr::Value("unpack".into())); }
        v_k = t3560[0].clone();
        v_v = t3560[1].clone();
        if ({ let t3524 = v_k.clone(); Py::Bool({ let t3525 = v_outputs.clone(); if !py::contains(&t3525, &t3524)? { true } else { false } }) }).truthy() {
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
fn s_3614(mut v_inputs: Py, mut v_outputs: Py) -> R {
    let mut v_k: Py = Py::None;
    let mut v_v: Py = Py::None;
    let mut v_want: Py = Py::None;
    let mut v_worst: Py = Py::None;
    let t3616: Result<Option<Py>, PyErr> = (|| -> Result<Option<Py>, PyErr> {
            v_want = { let t3615 = v_inputs.clone(); f_leslie_projection(py::getitem(&t3615, &Py::str("F"))?, py::getitem(&t3615, &Py::str("S"))?, py::getitem(&t3615, &Py::str("n0"))?, py::getitem(&t3615, &Py::str("t"))?, py::method(&t3615, "get", &[Py::str("transpose"), Py::Bool(false)])?)? };
            #[allow(unreachable_code)]
            Ok(None)
        })();
    match t3616 {
        Ok(Some(r)) => { return Ok(r); }
        Ok(None) => {}
        Err(e) if (!matches!(e, PyErr::Unsupported(_))) => {
            return Ok(Py::Float(f64::INFINITY));
        }
        Err(e) => return Err(e),
    }
    v_worst = Py::Float(0.0);
    for t3617 in py::iter(&(py::method(&v_want.clone(), "items", &[])?))? {
        let t3654 = (t3617).items()?;
        if t3654.len() != 2 { return Err(PyErr::Value("unpack".into())); }
        v_k = t3654[0].clone();
        v_v = t3654[1].clone();
        if ({ let t3618 = v_k.clone(); Py::Bool({ let t3619 = v_outputs.clone(); if !py::contains(&t3619, &t3618)? { true } else { false } }) }).truthy() {
            return Ok(Py::Float(f64::INFINITY));
        }
        v_worst = py::minmax(&[v_worst.clone(), f__rel(py::getitem(&v_outputs.clone(), &v_k.clone())?, v_v.clone())?], true)?;
    }
    return Ok(v_worst.clone());
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn s_3703(mut v_inputs: Py, mut v_outputs: Py) -> R {
    let mut v_k: Py = Py::None;
    let mut v_v: Py = Py::None;
    let mut v_want: Py = Py::None;
    let mut v_worst: Py = Py::None;
    let t3705: Result<Option<Py>, PyErr> = (|| -> Result<Option<Py>, PyErr> {
            v_want = { let t3704 = v_inputs.clone(); f_lineweaver_burk(py::getitem(&t3704, &Py::str("S"))?, py::getitem(&t3704, &Py::str("v"))?)? };
            #[allow(unreachable_code)]
            Ok(None)
        })();
    match t3705 {
        Ok(Some(r)) => { return Ok(r); }
        Ok(None) => {}
        Err(e) if (!matches!(e, PyErr::Unsupported(_))) => {
            return Ok(Py::Float(f64::INFINITY));
        }
        Err(e) => return Err(e),
    }
    v_worst = Py::Float(0.0);
    for t3706 in py::iter(&(py::method(&v_want.clone(), "items", &[])?))? {
        let t3743 = (t3706).items()?;
        if t3743.len() != 2 { return Err(PyErr::Value("unpack".into())); }
        v_k = t3743[0].clone();
        v_v = t3743[1].clone();
        if ({ let t3707 = v_k.clone(); Py::Bool({ let t3708 = v_outputs.clone(); if !py::contains(&t3708, &t3707)? { true } else { false } }) }).truthy() {
            return Ok(Py::Float(f64::INFINITY));
        }
        v_worst = py::minmax(&[v_worst.clone(), f__rel(py::getitem(&v_outputs.clone(), &v_k.clone())?, v_v.clone())?], true)?;
    }
    return Ok(v_worst.clone());
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn s_3762(mut v_inputs: Py, mut v_outputs: Py) -> R {
    let mut v_k: Py = Py::None;
    let mut v_v: Py = Py::None;
    let mut v_want: Py = Py::None;
    let mut v_worst: Py = Py::None;
    let t3764: Result<Option<Py>, PyErr> = (|| -> Result<Option<Py>, PyErr> {
            v_want = { let t3763 = v_inputs.clone(); f_linkage_disequilibrium(py::getitem(&t3763, &Py::str("n_AB"))?, py::getitem(&t3763, &Py::str("n_Ab"))?, py::getitem(&t3763, &Py::str("n_aB"))?, py::getitem(&t3763, &Py::str("n_ab"))?)? };
            #[allow(unreachable_code)]
            Ok(None)
        })();
    match t3764 {
        Ok(Some(r)) => { return Ok(r); }
        Ok(None) => {}
        Err(e) if (!matches!(e, PyErr::Unsupported(_))) => {
            return Ok(Py::Float(f64::INFINITY));
        }
        Err(e) => return Err(e),
    }
    v_worst = Py::Float(0.0);
    for t3765 in py::iter(&(py::method(&v_want.clone(), "items", &[])?))? {
        let t3802 = (t3765).items()?;
        if t3802.len() != 2 { return Err(PyErr::Value("unpack".into())); }
        v_k = t3802[0].clone();
        v_v = t3802[1].clone();
        if ({ let t3766 = v_k.clone(); Py::Bool({ let t3767 = v_outputs.clone(); if !py::contains(&t3767, &t3766)? { true } else { false } }) }).truthy() {
            return Ok(Py::Float(f64::INFINITY));
        }
        v_worst = py::minmax(&[v_worst.clone(), f__rel(py::getitem(&v_outputs.clone(), &v_k.clone())?, v_v.clone())?], true)?;
    }
    return Ok(v_worst.clone());
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn s_3886(mut v_inputs: Py, mut v_outputs: Py) -> R {
    let mut v_k: Py = Py::None;
    let mut v_v: Py = Py::None;
    let mut v_want: Py = Py::None;
    let mut v_worst: Py = Py::None;
    let t3888: Result<Option<Py>, PyErr> = (|| -> Result<Option<Py>, PyErr> {
            v_want = { let t3887 = v_inputs.clone(); f_local_affine(py::getitem(&t3887, &Py::str("seq1"))?, py::getitem(&t3887, &Py::str("seq2"))?, py::getitem(&t3887, &Py::str("match"))?, py::getitem(&t3887, &Py::str("mismatch"))?, py::getitem(&t3887, &Py::str("open"))?, py::getitem(&t3887, &Py::str("extend"))?)? };
            #[allow(unreachable_code)]
            Ok(None)
        })();
    match t3888 {
        Ok(Some(r)) => { return Ok(r); }
        Ok(None) => {}
        Err(e) if (!matches!(e, PyErr::Unsupported(_))) => {
            return Ok(Py::Float(f64::INFINITY));
        }
        Err(e) => return Err(e),
    }
    v_worst = Py::Float(0.0);
    for t3889 in py::iter(&(py::method(&v_want.clone(), "items", &[])?))? {
        let t3926 = (t3889).items()?;
        if t3926.len() != 2 { return Err(PyErr::Value("unpack".into())); }
        v_k = t3926[0].clone();
        v_v = t3926[1].clone();
        if ({ let t3890 = v_k.clone(); Py::Bool({ let t3891 = v_outputs.clone(); if !py::contains(&t3891, &t3890)? { true } else { false } }) }).truthy() {
            return Ok(Py::Float(f64::INFINITY));
        }
        v_worst = py::minmax(&[v_worst.clone(), f__rel(py::getitem(&v_outputs.clone(), &v_k.clone())?, v_v.clone())?], true)?;
    }
    return Ok(v_worst.clone());
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn s_4007(mut v_inputs: Py, mut v_outputs: Py) -> R {
    let mut v_k: Py = Py::None;
    let mut v_v: Py = Py::None;
    let mut v_want: Py = Py::None;
    let mut v_worst: Py = Py::None;
    let t4009: Result<Option<Py>, PyErr> = (|| -> Result<Option<Py>, PyErr> {
            v_want = { let t4008 = v_inputs.clone(); f_local_linear(py::getitem(&t4008, &Py::str("seq1"))?, py::getitem(&t4008, &Py::str("seq2"))?, py::getitem(&t4008, &Py::str("match"))?, py::getitem(&t4008, &Py::str("mismatch"))?, py::getitem(&t4008, &Py::str("gap"))?)? };
            #[allow(unreachable_code)]
            Ok(None)
        })();
    match t4009 {
        Ok(Some(r)) => { return Ok(r); }
        Ok(None) => {}
        Err(e) if (!matches!(e, PyErr::Unsupported(_))) => {
            return Ok(Py::Float(f64::INFINITY));
        }
        Err(e) => return Err(e),
    }
    v_worst = Py::Float(0.0);
    for t4010 in py::iter(&(py::method(&v_want.clone(), "items", &[])?))? {
        let t4047 = (t4010).items()?;
        if t4047.len() != 2 { return Err(PyErr::Value("unpack".into())); }
        v_k = t4047[0].clone();
        v_v = t4047[1].clone();
        if ({ let t4011 = v_k.clone(); Py::Bool({ let t4012 = v_outputs.clone(); if !py::contains(&t4012, &t4011)? { true } else { false } }) }).truthy() {
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
fn s_4079(mut v_inputs: Py, mut v_outputs: Py) -> R {
    let mut v_k: Py = Py::None;
    let mut v_v: Py = Py::None;
    let mut v_want: Py = Py::None;
    let mut v_worst: Py = Py::None;
    let t4081: Result<Option<Py>, PyErr> = (|| -> Result<Option<Py>, PyErr> {
            v_want = { let t4080 = v_inputs.clone(); f_longest_orf(py::getitem(&t4080, &Py::str("seq"))?)? };
            #[allow(unreachable_code)]
            Ok(None)
        })();
    match t4081 {
        Ok(Some(r)) => { return Ok(r); }
        Ok(None) => {}
        Err(e) if (!matches!(e, PyErr::Unsupported(_))) => {
            return Ok(Py::Float(f64::INFINITY));
        }
        Err(e) => return Err(e),
    }
    v_worst = Py::Float(0.0);
    for t4082 in py::iter(&(py::method(&v_want.clone(), "items", &[])?))? {
        let t4119 = (t4082).items()?;
        if t4119.len() != 2 { return Err(PyErr::Value("unpack".into())); }
        v_k = t4119[0].clone();
        v_v = t4119[1].clone();
        if ({ let t4083 = v_k.clone(); Py::Bool({ let t4084 = v_outputs.clone(); if !py::contains(&t4084, &t4083)? { true } else { false } }) }).truthy() {
            return Ok(Py::Float(f64::INFINITY));
        }
        v_worst = py::minmax(&[v_worst.clone(), f__rel(py::getitem(&v_outputs.clone(), &v_k.clone())?, v_v.clone())?], true)?;
    }
    return Ok(v_worst.clone());
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn s_4144(mut v_inputs: Py, mut v_outputs: Py) -> R {
    let mut v_k: Py = Py::None;
    let mut v_v: Py = Py::None;
    let mut v_want: Py = Py::None;
    let mut v_worst: Py = Py::None;
    let t4146: Result<Option<Py>, PyErr> = (|| -> Result<Option<Py>, PyErr> {
            v_want = { let t4145 = v_inputs.clone(); f_map_distance(py::getitem(&t4145, &Py::str("r"))?, py::getitem(&t4145, &Py::str("map"))?)? };
            #[allow(unreachable_code)]
            Ok(None)
        })();
    match t4146 {
        Ok(Some(r)) => { return Ok(r); }
        Ok(None) => {}
        Err(e) if (!matches!(e, PyErr::Unsupported(_))) => {
            return Ok(Py::Float(f64::INFINITY));
        }
        Err(e) => return Err(e),
    }
    v_worst = Py::Float(0.0);
    for t4147 in py::iter(&(py::method(&v_want.clone(), "items", &[])?))? {
        let t4184 = (t4147).items()?;
        if t4184.len() != 2 { return Err(PyErr::Value("unpack".into())); }
        v_k = t4184[0].clone();
        v_v = t4184[1].clone();
        if ({ let t4148 = v_k.clone(); Py::Bool({ let t4149 = v_outputs.clone(); if !py::contains(&t4149, &t4148)? { true } else { false } }) }).truthy() {
            return Ok(Py::Float(f64::INFINITY));
        }
        v_worst = py::minmax(&[v_worst.clone(), f__rel(py::getitem(&v_outputs.clone(), &v_k.clone())?, v_v.clone())?], true)?;
    }
    return Ok(v_worst.clone());
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn s_4261(mut v_inputs: Py, mut v_outputs: Py) -> R {
    let mut v_k: Py = Py::None;
    let mut v_v: Py = Py::None;
    let mut v_want: Py = Py::None;
    let mut v_worst: Py = Py::None;
    let t4263: Result<Option<Py>, PyErr> = (|| -> Result<Option<Py>, PyErr> {
            v_want = { let t4262 = v_inputs.clone(); f_mm_fit(py::getitem(&t4262, &Py::str("S"))?, py::getitem(&t4262, &Py::str("v"))?)? };
            #[allow(unreachable_code)]
            Ok(None)
        })();
    match t4263 {
        Ok(Some(r)) => { return Ok(r); }
        Ok(None) => {}
        Err(e) if (!matches!(e, PyErr::Unsupported(_))) => {
            return Ok(Py::Float(f64::INFINITY));
        }
        Err(e) => return Err(e),
    }
    v_worst = Py::Float(0.0);
    for t4264 in py::iter(&(py::method(&v_want.clone(), "items", &[])?))? {
        let t4301 = (t4264).items()?;
        if t4301.len() != 2 { return Err(PyErr::Value("unpack".into())); }
        v_k = t4301[0].clone();
        v_v = t4301[1].clone();
        if ({ let t4265 = v_k.clone(); Py::Bool({ let t4266 = v_outputs.clone(); if !py::contains(&t4266, &t4265)? { true } else { false } }) }).truthy() {
            return Ok(Py::Float(f64::INFINITY));
        }
        v_worst = py::minmax(&[v_worst.clone(), f__rel(py::getitem(&v_outputs.clone(), &v_k.clone())?, v_v.clone())?], true)?;
    }
    return Ok(v_worst.clone());
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn s_4359(mut v_inputs: Py, mut v_outputs: Py) -> R {
    let mut v_k: Py = Py::None;
    let mut v_v: Py = Py::None;
    let mut v_want: Py = Py::None;
    let mut v_worst: Py = Py::None;
    let t4361: Result<Option<Py>, PyErr> = (|| -> Result<Option<Py>, PyErr> {
            v_want = { let t4360 = v_inputs.clone(); f_mm_progress(py::getitem(&t4360, &Py::str("S0"))?, py::getitem(&t4360, &Py::str("Vmax"))?, py::getitem(&t4360, &Py::str("Km"))?, py::getitem(&t4360, &Py::str("t"))?, py::method(&t4360, "get", &[Py::str("first_order"), Py::Bool(false)])?)? };
            #[allow(unreachable_code)]
            Ok(None)
        })();
    match t4361 {
        Ok(Some(r)) => { return Ok(r); }
        Ok(None) => {}
        Err(e) if (!matches!(e, PyErr::Unsupported(_))) => {
            return Ok(Py::Float(f64::INFINITY));
        }
        Err(e) => return Err(e),
    }
    v_worst = Py::Float(0.0);
    for t4362 in py::iter(&(py::method(&v_want.clone(), "items", &[])?))? {
        let t4399 = (t4362).items()?;
        if t4399.len() != 2 { return Err(PyErr::Value("unpack".into())); }
        v_k = t4399[0].clone();
        v_v = t4399[1].clone();
        if ({ let t4363 = v_k.clone(); Py::Bool({ let t4364 = v_outputs.clone(); if !py::contains(&t4364, &t4363)? { true } else { false } }) }).truthy() {
            return Ok(Py::Float(f64::INFINITY));
        }
        v_worst = py::minmax(&[v_worst.clone(), f__rel(py::getitem(&v_outputs.clone(), &v_k.clone())?, v_v.clone())?], true)?;
    }
    return Ok(v_worst.clone());
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn s_4426(mut v_inputs: Py, mut v_outputs: Py) -> R {
    let mut v_k: Py = Py::None;
    let mut v_v: Py = Py::None;
    let mut v_want: Py = Py::None;
    let mut v_worst: Py = Py::None;
    let t4428: Result<Option<Py>, PyErr> = (|| -> Result<Option<Py>, PyErr> {
            v_want = { let t4427 = v_inputs.clone(); f_mm_time(py::getitem(&t4427, &Py::str("S0"))?, py::getitem(&t4427, &Py::str("Vmax"))?, py::getitem(&t4427, &Py::str("Km"))?, py::getitem(&t4427, &Py::str("X"))?, py::method(&t4427, "get", &[Py::str("zero_order"), Py::Bool(false)])?)? };
            #[allow(unreachable_code)]
            Ok(None)
        })();
    match t4428 {
        Ok(Some(r)) => { return Ok(r); }
        Ok(None) => {}
        Err(e) if (!matches!(e, PyErr::Unsupported(_))) => {
            return Ok(Py::Float(f64::INFINITY));
        }
        Err(e) => return Err(e),
    }
    v_worst = Py::Float(0.0);
    for t4429 in py::iter(&(py::method(&v_want.clone(), "items", &[])?))? {
        let t4466 = (t4429).items()?;
        if t4466.len() != 2 { return Err(PyErr::Value("unpack".into())); }
        v_k = t4466[0].clone();
        v_v = t4466[1].clone();
        if ({ let t4430 = v_k.clone(); Py::Bool({ let t4431 = v_outputs.clone(); if !py::contains(&t4431, &t4430)? { true } else { false } }) }).truthy() {
            return Ok(Py::Float(f64::INFINITY));
        }
        v_worst = py::minmax(&[v_worst.clone(), f__rel(py::getitem(&v_outputs.clone(), &v_k.clone())?, v_v.clone())?], true)?;
    }
    return Ok(v_worst.clone());
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn s_4501(mut v_inputs: Py, mut v_outputs: Py) -> R {
    let mut v_k: Py = Py::None;
    let mut v_v: Py = Py::None;
    let mut v_want: Py = Py::None;
    let mut v_worst: Py = Py::None;
    let t4503: Result<Option<Py>, PyErr> = (|| -> Result<Option<Py>, PyErr> {
            v_want = { let t4502 = v_inputs.clone(); f_moran_fixation(py::getitem(&t4502, &Py::str("N"))?, py::getitem(&t4502, &Py::str("i"))?, py::getitem(&t4502, &Py::str("r"))?)? };
            #[allow(unreachable_code)]
            Ok(None)
        })();
    match t4503 {
        Ok(Some(r)) => { return Ok(r); }
        Ok(None) => {}
        Err(e) if (!matches!(e, PyErr::Unsupported(_))) => {
            return Ok(Py::Float(f64::INFINITY));
        }
        Err(e) => return Err(e),
    }
    v_worst = Py::Float(0.0);
    for t4504 in py::iter(&(py::method(&v_want.clone(), "items", &[])?))? {
        let t4541 = (t4504).items()?;
        if t4541.len() != 2 { return Err(PyErr::Value("unpack".into())); }
        v_k = t4541[0].clone();
        v_v = t4541[1].clone();
        if ({ let t4505 = v_k.clone(); Py::Bool({ let t4506 = v_outputs.clone(); if !py::contains(&t4506, &t4505)? { true } else { false } }) }).truthy() {
            return Ok(Py::Float(f64::INFINITY));
        }
        v_worst = py::minmax(&[v_worst.clone(), f__rel(py::getitem(&v_outputs.clone(), &v_k.clone())?, v_v.clone())?], true)?;
    }
    return Ok(v_worst.clone());
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn s_4598(mut v_inputs: Py, mut v_outputs: Py) -> R {
    let mut v_k: Py = Py::None;
    let mut v_v: Py = Py::None;
    let mut v_want: Py = Py::None;
    let mut v_worst: Py = Py::None;
    let t4600: Result<Option<Py>, PyErr> = (|| -> Result<Option<Py>, PyErr> {
            v_want = { let t4599 = v_inputs.clone(); f_mutation_selection(py::getitem(&t4599, &Py::str("mu"))?, py::getitem(&t4599, &Py::str("s"))?, py::getitem(&t4599, &Py::str("h"))?)? };
            #[allow(unreachable_code)]
            Ok(None)
        })();
    match t4600 {
        Ok(Some(r)) => { return Ok(r); }
        Ok(None) => {}
        Err(e) if (!matches!(e, PyErr::Unsupported(_))) => {
            return Ok(Py::Float(f64::INFINITY));
        }
        Err(e) => return Err(e),
    }
    v_worst = Py::Float(0.0);
    for t4601 in py::iter(&(py::method(&v_want.clone(), "items", &[])?))? {
        let t4638 = (t4601).items()?;
        if t4638.len() != 2 { return Err(PyErr::Value("unpack".into())); }
        v_k = t4638[0].clone();
        v_v = t4638[1].clone();
        if ({ let t4602 = v_k.clone(); Py::Bool({ let t4603 = v_outputs.clone(); if !py::contains(&t4603, &t4602)? { true } else { false } }) }).truthy() {
            return Ok(Py::Float(f64::INFINITY));
        }
        v_worst = py::minmax(&[v_worst.clone(), f__rel(py::getitem(&v_outputs.clone(), &v_k.clone())?, v_v.clone())?], true)?;
    }
    return Ok(v_worst.clone());
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn s_4667(mut v_inputs: Py, mut v_outputs: Py) -> R {
    let mut v_k: Py = Py::None;
    let mut v_v: Py = Py::None;
    let mut v_want: Py = Py::None;
    let mut v_worst: Py = Py::None;
    let t4669: Result<Option<Py>, PyErr> = (|| -> Result<Option<Py>, PyErr> {
            v_want = { let t4668 = v_inputs.clone(); f_ne_fluctuating(py::getitem(&t4668, &Py::str("sizes"))?)? };
            #[allow(unreachable_code)]
            Ok(None)
        })();
    match t4669 {
        Ok(Some(r)) => { return Ok(r); }
        Ok(None) => {}
        Err(e) if (!matches!(e, PyErr::Unsupported(_))) => {
            return Ok(Py::Float(f64::INFINITY));
        }
        Err(e) => return Err(e),
    }
    v_worst = Py::Float(0.0);
    for t4670 in py::iter(&(py::method(&v_want.clone(), "items", &[])?))? {
        let t4707 = (t4670).items()?;
        if t4707.len() != 2 { return Err(PyErr::Value("unpack".into())); }
        v_k = t4707[0].clone();
        v_v = t4707[1].clone();
        if ({ let t4671 = v_k.clone(); Py::Bool({ let t4672 = v_outputs.clone(); if !py::contains(&t4672, &t4671)? { true } else { false } }) }).truthy() {
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
fn s_4715(mut v_inputs: Py, mut v_outputs: Py) -> R {
    let mut v_k: Py = Py::None;
    let mut v_v: Py = Py::None;
    let mut v_want: Py = Py::None;
    let mut v_worst: Py = Py::None;
    let t4717: Result<Option<Py>, PyErr> = (|| -> Result<Option<Py>, PyErr> {
            v_want = { let t4716 = v_inputs.clone(); f_ne_sex_ratio(py::getitem(&t4716, &Py::str("N_m"))?, py::getitem(&t4716, &Py::str("N_f"))?)? };
            #[allow(unreachable_code)]
            Ok(None)
        })();
    match t4717 {
        Ok(Some(r)) => { return Ok(r); }
        Ok(None) => {}
        Err(e) if (!matches!(e, PyErr::Unsupported(_))) => {
            return Ok(Py::Float(f64::INFINITY));
        }
        Err(e) => return Err(e),
    }
    v_worst = Py::Float(0.0);
    for t4718 in py::iter(&(py::method(&v_want.clone(), "items", &[])?))? {
        let t4755 = (t4718).items()?;
        if t4755.len() != 2 { return Err(PyErr::Value("unpack".into())); }
        v_k = t4755[0].clone();
        v_v = t4755[1].clone();
        if ({ let t4719 = v_k.clone(); Py::Bool({ let t4720 = v_outputs.clone(); if !py::contains(&t4720, &t4719)? { true } else { false } }) }).truthy() {
            return Ok(Py::Float(f64::INFINITY));
        }
        v_worst = py::minmax(&[v_worst.clone(), f__rel(py::getitem(&v_outputs.clone(), &v_k.clone())?, v_v.clone())?], true)?;
    }
    return Ok(v_worst.clone());
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn s_4835(mut v_inputs: Py, mut v_outputs: Py) -> R {
    let mut v_k: Py = Py::None;
    let mut v_v: Py = Py::None;
    let mut v_want: Py = Py::None;
    let mut v_worst: Py = Py::None;
    let t4837: Result<Option<Py>, PyErr> = (|| -> Result<Option<Py>, PyErr> {
            v_want = { let t4836 = v_inputs.clone(); f_neighbor_joining(py::getitem(&t4836, &Py::str("D"))?, py::method(&t4836, "get", &[Py::str("corrected"), Py::Bool(true)])?)? };
            #[allow(unreachable_code)]
            Ok(None)
        })();
    match t4837 {
        Ok(Some(r)) => { return Ok(r); }
        Ok(None) => {}
        Err(e) if (!matches!(e, PyErr::Unsupported(_))) => {
            return Ok(Py::Float(f64::INFINITY));
        }
        Err(e) => return Err(e),
    }
    v_worst = Py::Float(0.0);
    for t4838 in py::iter(&(py::method(&v_want.clone(), "items", &[])?))? {
        let t4875 = (t4838).items()?;
        if t4875.len() != 2 { return Err(PyErr::Value("unpack".into())); }
        v_k = t4875[0].clone();
        v_v = t4875[1].clone();
        if ({ let t4839 = v_k.clone(); Py::Bool({ let t4840 = v_outputs.clone(); if !py::contains(&t4840, &t4839)? { true } else { false } }) }).truthy() {
            return Ok(Py::Float(f64::INFINITY));
        }
        v_worst = py::minmax(&[v_worst.clone(), f__rel(py::getitem(&v_outputs.clone(), &v_k.clone())?, v_v.clone())?], true)?;
    }
    return Ok(v_worst.clone());
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn s_4924(mut v_inputs: Py, mut v_outputs: Py) -> R {
    let mut v_k: Py = Py::None;
    let mut v_v: Py = Py::None;
    let mut v_want: Py = Py::None;
    let mut v_worst: Py = Py::None;
    let t4926: Result<Option<Py>, PyErr> = (|| -> Result<Option<Py>, PyErr> {
            v_want = { let t4925 = v_inputs.clone(); f_nucleic_mass(py::getitem(&t4925, &Py::str("seq"))?, py::getitem(&t4925, &Py::str("molecule"))?, py::getitem(&t4925, &Py::str("strands"))?)? };
            #[allow(unreachable_code)]
            Ok(None)
        })();
    match t4926 {
        Ok(Some(r)) => { return Ok(r); }
        Ok(None) => {}
        Err(e) if (!matches!(e, PyErr::Unsupported(_))) => {
            return Ok(Py::Float(f64::INFINITY));
        }
        Err(e) => return Err(e),
    }
    v_worst = Py::Float(0.0);
    for t4927 in py::iter(&(py::method(&v_want.clone(), "items", &[])?))? {
        let t4964 = (t4927).items()?;
        if t4964.len() != 2 { return Err(PyErr::Value("unpack".into())); }
        v_k = t4964[0].clone();
        v_v = t4964[1].clone();
        if ({ let t4928 = v_k.clone(); Py::Bool({ let t4929 = v_outputs.clone(); if !py::contains(&t4929, &t4928)? { true } else { false } }) }).truthy() {
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
            v_want = { let t4981 = v_inputs.clone(); f_nucleotide_diversity(py::getitem(&t4981, &Py::str("seqs"))?)? };
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
fn s_5062(mut v_inputs: Py, mut v_outputs: Py) -> R {
    let mut v_k: Py = Py::None;
    let mut v_v: Py = Py::None;
    let mut v_want: Py = Py::None;
    let mut v_worst: Py = Py::None;
    let t5064: Result<Option<Py>, PyErr> = (|| -> Result<Option<Py>, PyErr> {
            v_want = { let t5063 = v_inputs.clone(); f_pedigree_inbreeding(py::getitem(&t5063, &Py::str("pedigree"))?, py::getitem(&t5063, &Py::str("individual"))?)? };
            #[allow(unreachable_code)]
            Ok(None)
        })();
    match t5064 {
        Ok(Some(r)) => { return Ok(r); }
        Ok(None) => {}
        Err(e) if (!matches!(e, PyErr::Unsupported(_))) => {
            return Ok(Py::Float(f64::INFINITY));
        }
        Err(e) => return Err(e),
    }
    v_worst = Py::Float(0.0);
    for t5065 in py::iter(&(py::method(&v_want.clone(), "items", &[])?))? {
        let t5102 = (t5065).items()?;
        if t5102.len() != 2 { return Err(PyErr::Value("unpack".into())); }
        v_k = t5102[0].clone();
        v_v = t5102[1].clone();
        if ({ let t5066 = v_k.clone(); Py::Bool({ let t5067 = v_outputs.clone(); if !py::contains(&t5067, &t5066)? { true } else { false } }) }).truthy() {
            return Ok(Py::Float(f64::INFINITY));
        }
        v_worst = py::minmax(&[v_worst.clone(), f__rel(py::getitem(&v_outputs.clone(), &v_k.clone())?, v_v.clone())?], true)?;
    }
    return Ok(v_worst.clone());
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn s_5110(mut v_inputs: Py, mut v_outputs: Py) -> R {
    let mut v_k: Py = Py::None;
    let mut v_v: Py = Py::None;
    let mut v_want: Py = Py::None;
    let mut v_worst: Py = Py::None;
    let t5112: Result<Option<Py>, PyErr> = (|| -> Result<Option<Py>, PyErr> {
            v_want = { let t5111 = v_inputs.clone(); f_petersen_mle(py::getitem(&t5111, &Py::str("M"))?, py::getitem(&t5111, &Py::str("C"))?, py::getitem(&t5111, &Py::str("R"))?, py::method(&t5111, "get", &[Py::str("ceil"), Py::Bool(false)])?)? };
            #[allow(unreachable_code)]
            Ok(None)
        })();
    match t5112 {
        Ok(Some(r)) => { return Ok(r); }
        Ok(None) => {}
        Err(e) if (!matches!(e, PyErr::Unsupported(_))) => {
            return Ok(Py::Float(f64::INFINITY));
        }
        Err(e) => return Err(e),
    }
    v_worst = Py::Float(0.0);
    for t5113 in py::iter(&(py::method(&v_want.clone(), "items", &[])?))? {
        let t5150 = (t5113).items()?;
        if t5150.len() != 2 { return Err(PyErr::Value("unpack".into())); }
        v_k = t5150[0].clone();
        v_v = t5150[1].clone();
        if ({ let t5114 = v_k.clone(); Py::Bool({ let t5115 = v_outputs.clone(); if !py::contains(&t5115, &t5114)? { true } else { false } }) }).truthy() {
            return Ok(Py::Float(f64::INFINITY));
        }
        v_worst = py::minmax(&[v_worst.clone(), f__rel(py::getitem(&v_outputs.clone(), &v_k.clone())?, v_v.clone())?], true)?;
    }
    return Ok(v_worst.clone());
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn s_5189(mut v_inputs: Py, mut v_outputs: Py) -> R {
    let mut v_k: Py = Py::None;
    let mut v_v: Py = Py::None;
    let mut v_want: Py = Py::None;
    let mut v_worst: Py = Py::None;
    let t5191: Result<Option<Py>, PyErr> = (|| -> Result<Option<Py>, PyErr> {
            v_want = { let t5190 = v_inputs.clone(); f_extinction(py::getitem(&t5190, &Py::str("seq"))?)? };
            #[allow(unreachable_code)]
            Ok(None)
        })();
    match t5191 {
        Ok(Some(r)) => { return Ok(r); }
        Ok(None) => {}
        Err(e) if (!matches!(e, PyErr::Unsupported(_))) => {
            return Ok(Py::Float(f64::INFINITY));
        }
        Err(e) => return Err(e),
    }
    v_worst = Py::Float(0.0);
    for t5192 in py::iter(&(py::method(&v_want.clone(), "items", &[])?))? {
        let t5229 = (t5192).items()?;
        if t5229.len() != 2 { return Err(PyErr::Value("unpack".into())); }
        v_k = t5229[0].clone();
        v_v = t5229[1].clone();
        if ({ let t5193 = v_k.clone(); Py::Bool({ let t5194 = v_outputs.clone(); if !py::contains(&t5194, &t5193)? { true } else { false } }) }).truthy() {
            return Ok(Py::Float(f64::INFINITY));
        }
        v_worst = py::minmax(&[v_worst.clone(), f__rel(py::getitem(&v_outputs.clone(), &v_k.clone())?, v_v.clone())?], true)?;
    }
    return Ok(v_worst.clone());
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn s_5309(mut v_inputs: Py, mut v_outputs: Py) -> R {
    let mut v_k: Py = Py::None;
    let mut v_v: Py = Py::None;
    let mut v_want: Py = Py::None;
    let mut v_worst: Py = Py::None;
    let t5311: Result<Option<Py>, PyErr> = (|| -> Result<Option<Py>, PyErr> {
            v_want = { let t5310 = v_inputs.clone(); f_protein_global(py::getitem(&t5310, &Py::str("seq1"))?, py::getitem(&t5310, &Py::str("seq2"))?, py::getitem(&t5310, &Py::str("open"))?, py::getitem(&t5310, &Py::str("extend"))?)? };
            #[allow(unreachable_code)]
            Ok(None)
        })();
    match t5311 {
        Ok(Some(r)) => { return Ok(r); }
        Ok(None) => {}
        Err(e) if (!matches!(e, PyErr::Unsupported(_))) => {
            return Ok(Py::Float(f64::INFINITY));
        }
        Err(e) => return Err(e),
    }
    v_worst = Py::Float(0.0);
    for t5312 in py::iter(&(py::method(&v_want.clone(), "items", &[])?))? {
        let t5349 = (t5312).items()?;
        if t5349.len() != 2 { return Err(PyErr::Value("unpack".into())); }
        v_k = t5349[0].clone();
        v_v = t5349[1].clone();
        if ({ let t5313 = v_k.clone(); Py::Bool({ let t5314 = v_outputs.clone(); if !py::contains(&t5314, &t5313)? { true } else { false } }) }).truthy() {
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
fn s_5429(mut v_inputs: Py, mut v_outputs: Py) -> R {
    let mut v_k: Py = Py::None;
    let mut v_v: Py = Py::None;
    let mut v_want: Py = Py::None;
    let mut v_worst: Py = Py::None;
    let t5431: Result<Option<Py>, PyErr> = (|| -> Result<Option<Py>, PyErr> {
            v_want = { let t5430 = v_inputs.clone(); f_protein_local(py::getitem(&t5430, &Py::str("seq1"))?, py::getitem(&t5430, &Py::str("seq2"))?, py::getitem(&t5430, &Py::str("open"))?, py::getitem(&t5430, &Py::str("extend"))?)? };
            #[allow(unreachable_code)]
            Ok(None)
        })();
    match t5431 {
        Ok(Some(r)) => { return Ok(r); }
        Ok(None) => {}
        Err(e) if (!matches!(e, PyErr::Unsupported(_))) => {
            return Ok(Py::Float(f64::INFINITY));
        }
        Err(e) => return Err(e),
    }
    v_worst = Py::Float(0.0);
    for t5432 in py::iter(&(py::method(&v_want.clone(), "items", &[])?))? {
        let t5469 = (t5432).items()?;
        if t5469.len() != 2 { return Err(PyErr::Value("unpack".into())); }
        v_k = t5469[0].clone();
        v_v = t5469[1].clone();
        if ({ let t5433 = v_k.clone(); Py::Bool({ let t5434 = v_outputs.clone(); if !py::contains(&t5434, &t5433)? { true } else { false } }) }).truthy() {
            return Ok(Py::Float(f64::INFINITY));
        }
        v_worst = py::minmax(&[v_worst.clone(), f__rel(py::getitem(&v_outputs.clone(), &v_k.clone())?, v_v.clone())?], true)?;
    }
    return Ok(v_worst.clone());
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn s_5508(mut v_inputs: Py, mut v_outputs: Py) -> R {
    let mut v_k: Py = Py::None;
    let mut v_v: Py = Py::None;
    let mut v_want: Py = Py::None;
    let mut v_worst: Py = Py::None;
    let t5510: Result<Option<Py>, PyErr> = (|| -> Result<Option<Py>, PyErr> {
            v_want = { let t5509 = v_inputs.clone(); f_protein_mass(py::getitem(&t5509, &Py::str("seq"))?)? };
            #[allow(unreachable_code)]
            Ok(None)
        })();
    match t5510 {
        Ok(Some(r)) => { return Ok(r); }
        Ok(None) => {}
        Err(e) if (!matches!(e, PyErr::Unsupported(_))) => {
            return Ok(Py::Float(f64::INFINITY));
        }
        Err(e) => return Err(e),
    }
    v_worst = Py::Float(0.0);
    for t5511 in py::iter(&(py::method(&v_want.clone(), "items", &[])?))? {
        let t5548 = (t5511).items()?;
        if t5548.len() != 2 { return Err(PyErr::Value("unpack".into())); }
        v_k = t5548[0].clone();
        v_v = t5548[1].clone();
        if ({ let t5512 = v_k.clone(); Py::Bool({ let t5513 = v_outputs.clone(); if !py::contains(&t5513, &t5512)? { true } else { false } }) }).truthy() {
            return Ok(Py::Float(f64::INFINITY));
        }
        v_worst = py::minmax(&[v_worst.clone(), f__rel(py::getitem(&v_outputs.clone(), &v_k.clone())?, v_v.clone())?], true)?;
    }
    return Ok(v_worst.clone());
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn s_5603(mut v_inputs: Py, mut v_outputs: Py) -> R {
    let mut v_k: Py = Py::None;
    let mut v_v: Py = Py::None;
    let mut v_want: Py = Py::None;
    let mut v_worst: Py = Py::None;
    let t5605: Result<Option<Py>, PyErr> = (|| -> Result<Option<Py>, PyErr> {
            v_want = { let t5604 = v_inputs.clone(); f_net_charge(py::getitem(&t5604, &Py::str("seq"))?, py::getitem(&t5604, &Py::str("pH"))?)? };
            #[allow(unreachable_code)]
            Ok(None)
        })();
    match t5605 {
        Ok(Some(r)) => { return Ok(r); }
        Ok(None) => {}
        Err(e) if (!matches!(e, PyErr::Unsupported(_))) => {
            return Ok(Py::Float(f64::INFINITY));
        }
        Err(e) => return Err(e),
    }
    v_worst = Py::Float(0.0);
    for t5606 in py::iter(&(py::method(&v_want.clone(), "items", &[])?))? {
        let t5643 = (t5606).items()?;
        if t5643.len() != 2 { return Err(PyErr::Value("unpack".into())); }
        v_k = t5643[0].clone();
        v_v = t5643[1].clone();
        if ({ let t5607 = v_k.clone(); Py::Bool({ let t5608 = v_outputs.clone(); if !py::contains(&t5608, &t5607)? { true } else { false } }) }).truthy() {
            return Ok(Py::Float(f64::INFINITY));
        }
        v_worst = py::minmax(&[v_worst.clone(), f__rel(py::getitem(&v_outputs.clone(), &v_k.clone())?, v_v.clone())?], true)?;
    }
    return Ok(v_worst.clone());
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn s_5670(mut v_inputs: Py, mut v_outputs: Py) -> R {
    let mut v_k: Py = Py::None;
    let mut v_v: Py = Py::None;
    let mut v_want: Py = Py::None;
    let mut v_worst: Py = Py::None;
    let t5672: Result<Option<Py>, PyErr> = (|| -> Result<Option<Py>, PyErr> {
            v_want = { let t5671 = v_inputs.clone(); f_rarefaction(py::getitem(&t5671, &Py::str("counts"))?, py::getitem(&t5671, &Py::str("n"))?, py::method(&t5671, "get", &[Py::str("with_replacement"), Py::Bool(false)])?)? };
            #[allow(unreachable_code)]
            Ok(None)
        })();
    match t5672 {
        Ok(Some(r)) => { return Ok(r); }
        Ok(None) => {}
        Err(e) if (!matches!(e, PyErr::Unsupported(_))) => {
            return Ok(Py::Float(f64::INFINITY));
        }
        Err(e) => return Err(e),
    }
    v_worst = Py::Float(0.0);
    for t5673 in py::iter(&(py::method(&v_want.clone(), "items", &[])?))? {
        let t5710 = (t5673).items()?;
        if t5710.len() != 2 { return Err(PyErr::Value("unpack".into())); }
        v_k = t5710[0].clone();
        v_v = t5710[1].clone();
        if ({ let t5674 = v_k.clone(); Py::Bool({ let t5675 = v_outputs.clone(); if !py::contains(&t5675, &t5674)? { true } else { false } }) }).truthy() {
            return Ok(Py::Float(f64::INFINITY));
        }
        v_worst = py::minmax(&[v_worst.clone(), f__rel(py::getitem(&v_outputs.clone(), &v_k.clone())?, v_v.clone())?], true)?;
    }
    return Ok(v_worst.clone());
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn s_5742(mut v_inputs: Py, mut v_outputs: Py) -> R {
    let mut v_k: Py = Py::None;
    let mut v_v: Py = Py::None;
    let mut v_want: Py = Py::None;
    let mut v_worst: Py = Py::None;
    let t5744: Result<Option<Py>, PyErr> = (|| -> Result<Option<Py>, PyErr> {
            v_want = { let t5743 = v_inputs.clone(); f_recessive_lethal(py::getitem(&t5743, &Py::str("q0"))?, py::getitem(&t5743, &Py::str("t"))?)? };
            #[allow(unreachable_code)]
            Ok(None)
        })();
    match t5744 {
        Ok(Some(r)) => { return Ok(r); }
        Ok(None) => {}
        Err(e) if (!matches!(e, PyErr::Unsupported(_))) => {
            return Ok(Py::Float(f64::INFINITY));
        }
        Err(e) => return Err(e),
    }
    v_worst = Py::Float(0.0);
    for t5745 in py::iter(&(py::method(&v_want.clone(), "items", &[])?))? {
        let t5782 = (t5745).items()?;
        if t5782.len() != 2 { return Err(PyErr::Value("unpack".into())); }
        v_k = t5782[0].clone();
        v_v = t5782[1].clone();
        if ({ let t5746 = v_k.clone(); Py::Bool({ let t5747 = v_outputs.clone(); if !py::contains(&t5747, &t5746)? { true } else { false } }) }).truthy() {
            return Ok(Py::Float(f64::INFINITY));
        }
        v_worst = py::minmax(&[v_worst.clone(), f__rel(py::getitem(&v_outputs.clone(), &v_k.clone())?, v_v.clone())?], true)?;
    }
    return Ok(v_worst.clone());
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn s_5806(mut v_inputs: Py, mut v_outputs: Py) -> R {
    let mut v_k: Py = Py::None;
    let mut v_v: Py = Py::None;
    let mut v_want: Py = Py::None;
    let mut v_worst: Py = Py::None;
    let t5808: Result<Option<Py>, PyErr> = (|| -> Result<Option<Py>, PyErr> {
            v_want = { let t5807 = v_inputs.clone(); f_recombination_fraction(py::getitem(&t5807, &Py::str("d_cM"))?, py::getitem(&t5807, &Py::str("map"))?)? };
            #[allow(unreachable_code)]
            Ok(None)
        })();
    match t5808 {
        Ok(Some(r)) => { return Ok(r); }
        Ok(None) => {}
        Err(e) if (!matches!(e, PyErr::Unsupported(_))) => {
            return Ok(Py::Float(f64::INFINITY));
        }
        Err(e) => return Err(e),
    }
    v_worst = Py::Float(0.0);
    for t5809 in py::iter(&(py::method(&v_want.clone(), "items", &[])?))? {
        let t5846 = (t5809).items()?;
        if t5846.len() != 2 { return Err(PyErr::Value("unpack".into())); }
        v_k = t5846[0].clone();
        v_v = t5846[1].clone();
        if ({ let t5810 = v_k.clone(); Py::Bool({ let t5811 = v_outputs.clone(); if !py::contains(&t5811, &t5810)? { true } else { false } }) }).truthy() {
            return Ok(Py::Float(f64::INFINITY));
        }
        v_worst = py::minmax(&[v_worst.clone(), f__rel(py::getitem(&v_outputs.clone(), &v_k.clone())?, v_v.clone())?], true)?;
    }
    return Ok(v_worst.clone());
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn s_5875(mut v_inputs: Py, mut v_outputs: Py) -> R {
    let mut v_k: Py = Py::None;
    let mut v_v: Py = Py::None;
    let mut v_want: Py = Py::None;
    let mut v_worst: Py = Py::None;
    let t5877: Result<Option<Py>, PyErr> = (|| -> Result<Option<Py>, PyErr> {
            v_want = { let t5876 = v_inputs.clone(); f_restriction_sites(py::getitem(&t5876, &Py::str("seq"))?, py::getitem(&t5876, &Py::str("enzyme"))?)? };
            #[allow(unreachable_code)]
            Ok(None)
        })();
    match t5877 {
        Ok(Some(r)) => { return Ok(r); }
        Ok(None) => {}
        Err(e) if (!matches!(e, PyErr::Unsupported(_))) => {
            return Ok(Py::Float(f64::INFINITY));
        }
        Err(e) => return Err(e),
    }
    v_worst = Py::Float(0.0);
    for t5878 in py::iter(&(py::method(&v_want.clone(), "items", &[])?))? {
        let t5915 = (t5878).items()?;
        if t5915.len() != 2 { return Err(PyErr::Value("unpack".into())); }
        v_k = t5915[0].clone();
        v_v = t5915[1].clone();
        if ({ let t5879 = v_k.clone(); Py::Bool({ let t5880 = v_outputs.clone(); if !py::contains(&t5880, &t5879)? { true } else { false } }) }).truthy() {
            return Ok(Py::Float(f64::INFINITY));
        }
        v_worst = py::minmax(&[v_worst.clone(), f__rel(py::getitem(&v_outputs.clone(), &v_k.clone())?, v_v.clone())?], true)?;
    }
    return Ok(v_worst.clone());
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn s_5939(mut v_inputs: Py, mut v_outputs: Py) -> R {
    let mut v_k: Py = Py::None;
    let mut v_v: Py = Py::None;
    let mut v_want: Py = Py::None;
    let mut v_worst: Py = Py::None;
    let t5941: Result<Option<Py>, PyErr> = (|| -> Result<Option<Py>, PyErr> {
            v_want = { let t5940 = v_inputs.clone(); f_reverse_complement(py::getitem(&t5940, &Py::str("seq"))?)? };
            #[allow(unreachable_code)]
            Ok(None)
        })();
    match t5941 {
        Ok(Some(r)) => { return Ok(r); }
        Ok(None) => {}
        Err(e) if (!matches!(e, PyErr::Unsupported(_))) => {
            return Ok(Py::Float(f64::INFINITY));
        }
        Err(e) => return Err(e),
    }
    v_worst = Py::Float(0.0);
    for t5942 in py::iter(&(py::method(&v_want.clone(), "items", &[])?))? {
        let t5979 = (t5942).items()?;
        if t5979.len() != 2 { return Err(PyErr::Value("unpack".into())); }
        v_k = t5979[0].clone();
        v_v = t5979[1].clone();
        if ({ let t5943 = v_k.clone(); Py::Bool({ let t5944 = v_outputs.clone(); if !py::contains(&t5944, &t5943)? { true } else { false } }) }).truthy() {
            return Ok(Py::Float(f64::INFINITY));
        }
        v_worst = py::minmax(&[v_worst.clone(), f__rel(py::getitem(&v_outputs.clone(), &v_k.clone())?, v_v.clone())?], true)?;
    }
    return Ok(v_worst.clone());
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn s_6015(mut v_inputs: Py, mut v_outputs: Py) -> R {
    let mut v_k: Py = Py::None;
    let mut v_v: Py = Py::None;
    let mut v_want: Py = Py::None;
    let mut v_worst: Py = Py::None;
    let t6017: Result<Option<Py>, PyErr> = (|| -> Result<Option<Py>, PyErr> {
            v_want = { let t6016 = v_inputs.clone(); f_richness_indices(py::getitem(&t6016, &Py::str("counts"))?, py::method(&t6016, "get", &[Py::str("log10"), Py::Bool(false)])?)? };
            #[allow(unreachable_code)]
            Ok(None)
        })();
    match t6017 {
        Ok(Some(r)) => { return Ok(r); }
        Ok(None) => {}
        Err(e) if (!matches!(e, PyErr::Unsupported(_))) => {
            return Ok(Py::Float(f64::INFINITY));
        }
        Err(e) => return Err(e),
    }
    v_worst = Py::Float(0.0);
    for t6018 in py::iter(&(py::method(&v_want.clone(), "items", &[])?))? {
        let t6055 = (t6018).items()?;
        if t6055.len() != 2 { return Err(PyErr::Value("unpack".into())); }
        v_k = t6055[0].clone();
        v_v = t6055[1].clone();
        if ({ let t6019 = v_k.clone(); Py::Bool({ let t6020 = v_outputs.clone(); if !py::contains(&t6020, &t6019)? { true } else { false } }) }).truthy() {
            return Ok(Py::Float(f64::INFINITY));
        }
        v_worst = py::minmax(&[v_worst.clone(), f__rel(py::getitem(&v_outputs.clone(), &v_k.clone())?, v_v.clone())?], true)?;
    }
    return Ok(v_worst.clone());
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn s_6115(mut v_inputs: Py, mut v_outputs: Py) -> R {
    let mut v_k: Py = Py::None;
    let mut v_v: Py = Py::None;
    let mut v_want: Py = Py::None;
    let mut v_worst: Py = Py::None;
    let t6117: Result<Option<Py>, PyErr> = (|| -> Result<Option<Py>, PyErr> {
            v_want = { let t6116 = v_inputs.clone(); f_robinson_foulds(py::getitem(&t6116, &Py::str("tree1"))?, py::getitem(&t6116, &Py::str("tree2"))?, py::method(&t6116, "get", &[Py::str("one_way"), Py::Bool(false)])?)? };
            #[allow(unreachable_code)]
            Ok(None)
        })();
    match t6117 {
        Ok(Some(r)) => { return Ok(r); }
        Ok(None) => {}
        Err(e) if (!matches!(e, PyErr::Unsupported(_))) => {
            return Ok(Py::Float(f64::INFINITY));
        }
        Err(e) => return Err(e),
    }
    v_worst = Py::Float(0.0);
    for t6118 in py::iter(&(py::method(&v_want.clone(), "items", &[])?))? {
        let t6155 = (t6118).items()?;
        if t6155.len() != 2 { return Err(PyErr::Value("unpack".into())); }
        v_k = t6155[0].clone();
        v_v = t6155[1].clone();
        if ({ let t6119 = v_k.clone(); Py::Bool({ let t6120 = v_outputs.clone(); if !py::contains(&t6120, &t6119)? { true } else { false } }) }).truthy() {
            return Ok(Py::Float(f64::INFINITY));
        }
        v_worst = py::minmax(&[v_worst.clone(), f__rel(py::getitem(&v_outputs.clone(), &v_k.clone())?, v_v.clone())?], true)?;
    }
    return Ok(v_worst.clone());
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn s_6190(mut v_inputs: Py, mut v_outputs: Py) -> R {
    let mut v_k: Py = Py::None;
    let mut v_v: Py = Py::None;
    let mut v_want: Py = Py::None;
    let mut v_worst: Py = Py::None;
    let t6192: Result<Option<Py>, PyErr> = (|| -> Result<Option<Py>, PyErr> {
            v_want = { let t6191 = v_inputs.clone(); f_schumacher(py::getitem(&t6191, &Py::str("C"))?, py::getitem(&t6191, &Py::str("R"))?, py::getitem(&t6191, &Py::str("M"))?, py::method(&t6191, "get", &[Py::str("schnabel"), Py::Bool(false)])?)? };
            #[allow(unreachable_code)]
            Ok(None)
        })();
    match t6192 {
        Ok(Some(r)) => { return Ok(r); }
        Ok(None) => {}
        Err(e) if (!matches!(e, PyErr::Unsupported(_))) => {
            return Ok(Py::Float(f64::INFINITY));
        }
        Err(e) => return Err(e),
    }
    v_worst = Py::Float(0.0);
    for t6193 in py::iter(&(py::method(&v_want.clone(), "items", &[])?))? {
        let t6230 = (t6193).items()?;
        if t6230.len() != 2 { return Err(PyErr::Value("unpack".into())); }
        v_k = t6230[0].clone();
        v_v = t6230[1].clone();
        if ({ let t6194 = v_k.clone(); Py::Bool({ let t6195 = v_outputs.clone(); if !py::contains(&t6195, &t6194)? { true } else { false } }) }).truthy() {
            return Ok(Py::Float(f64::INFINITY));
        }
        v_worst = py::minmax(&[v_worst.clone(), f__rel(py::getitem(&v_outputs.clone(), &v_k.clone())?, v_v.clone())?], true)?;
    }
    return Ok(v_worst.clone());
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn s_6257(mut v_inputs: Py, mut v_outputs: Py) -> R {
    let mut v_k: Py = Py::None;
    let mut v_v: Py = Py::None;
    let mut v_want: Py = Py::None;
    let mut v_worst: Py = Py::None;
    let t6259: Result<Option<Py>, PyErr> = (|| -> Result<Option<Py>, PyErr> {
            v_want = { let t6258 = v_inputs.clone(); f_ss_fraction(py::getitem(&t6258, &Py::str("seq"))?, py::method(&t6258, "get", &[Py::str("order"), Py::tuple(vec![Py::str("helix"), Py::str("turn"), Py::str("sheet")])])?)? };
            #[allow(unreachable_code)]
            Ok(None)
        })();
    match t6259 {
        Ok(Some(r)) => { return Ok(r); }
        Ok(None) => {}
        Err(e) if (!matches!(e, PyErr::Unsupported(_))) => {
            return Ok(Py::Float(f64::INFINITY));
        }
        Err(e) => return Err(e),
    }
    v_worst = Py::Float(0.0);
    for t6260 in py::iter(&(py::method(&v_want.clone(), "items", &[])?))? {
        let t6297 = (t6260).items()?;
        if t6297.len() != 2 { return Err(PyErr::Value("unpack".into())); }
        v_k = t6297[0].clone();
        v_v = t6297[1].clone();
        if ({ let t6261 = v_k.clone(); Py::Bool({ let t6262 = v_outputs.clone(); if !py::contains(&t6262, &t6261)? { true } else { false } }) }).truthy() {
            return Ok(Py::Float(f64::INFINITY));
        }
        v_worst = py::minmax(&[v_worst.clone(), f__rel(py::getitem(&v_outputs.clone(), &v_k.clone())?, v_v.clone())?], true)?;
    }
    return Ok(v_worst.clone());
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn s_6376(mut v_inputs: Py, mut v_outputs: Py) -> R {
    let mut v_k: Py = Py::None;
    let mut v_v: Py = Py::None;
    let mut v_want: Py = Py::None;
    let mut v_worst: Py = Py::None;
    let t6378: Result<Option<Py>, PyErr> = (|| -> Result<Option<Py>, PyErr> {
            v_want = { let t6377 = v_inputs.clone(); f_semiglobal(py::getitem(&t6377, &Py::str("seq1"))?, py::getitem(&t6377, &Py::str("seq2"))?, py::getitem(&t6377, &Py::str("match"))?, py::getitem(&t6377, &Py::str("mismatch"))?, py::getitem(&t6377, &Py::str("gap"))?)? };
            #[allow(unreachable_code)]
            Ok(None)
        })();
    match t6378 {
        Ok(Some(r)) => { return Ok(r); }
        Ok(None) => {}
        Err(e) if (!matches!(e, PyErr::Unsupported(_))) => {
            return Ok(Py::Float(f64::INFINITY));
        }
        Err(e) => return Err(e),
    }
    v_worst = Py::Float(0.0);
    for t6379 in py::iter(&(py::method(&v_want.clone(), "items", &[])?))? {
        let t6416 = (t6379).items()?;
        if t6416.len() != 2 { return Err(PyErr::Value("unpack".into())); }
        v_k = t6416[0].clone();
        v_v = t6416[1].clone();
        if ({ let t6380 = v_k.clone(); Py::Bool({ let t6381 = v_outputs.clone(); if !py::contains(&t6381, &t6380)? { true } else { false } }) }).truthy() {
            return Ok(Py::Float(f64::INFINITY));
        }
        v_worst = py::minmax(&[v_worst.clone(), f__rel(py::getitem(&v_outputs.clone(), &v_k.clone())?, v_v.clone())?], true)?;
    }
    return Ok(v_worst.clone());
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn s_6454(mut v_inputs: Py, mut v_outputs: Py) -> R {
    let mut v_k: Py = Py::None;
    let mut v_v: Py = Py::None;
    let mut v_want: Py = Py::None;
    let mut v_worst: Py = Py::None;
    let t6456: Result<Option<Py>, PyErr> = (|| -> Result<Option<Py>, PyErr> {
            v_want = { let t6455 = v_inputs.clone(); f_shannon(py::getitem(&t6455, &Py::str("counts"))?, py::getitem(&t6455, &Py::str("base"))?, py::method(&t6455, "get", &[Py::str("wrong"), Py::Bool(false)])?)? };
            #[allow(unreachable_code)]
            Ok(None)
        })();
    match t6456 {
        Ok(Some(r)) => { return Ok(r); }
        Ok(None) => {}
        Err(e) if (!matches!(e, PyErr::Unsupported(_))) => {
            return Ok(Py::Float(f64::INFINITY));
        }
        Err(e) => return Err(e),
    }
    v_worst = Py::Float(0.0);
    for t6457 in py::iter(&(py::method(&v_want.clone(), "items", &[])?))? {
        let t6494 = (t6457).items()?;
        if t6494.len() != 2 { return Err(PyErr::Value("unpack".into())); }
        v_k = t6494[0].clone();
        v_v = t6494[1].clone();
        if ({ let t6458 = v_k.clone(); Py::Bool({ let t6459 = v_outputs.clone(); if !py::contains(&t6459, &t6458)? { true } else { false } }) }).truthy() {
            return Ok(Py::Float(f64::INFINITY));
        }
        v_worst = py::minmax(&[v_worst.clone(), f__rel(py::getitem(&v_outputs.clone(), &v_k.clone())?, v_v.clone())?], true)?;
    }
    return Ok(v_worst.clone());
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn s_6523(mut v_inputs: Py, mut v_outputs: Py) -> R {
    let mut v_k: Py = Py::None;
    let mut v_v: Py = Py::None;
    let mut v_want: Py = Py::None;
    let mut v_worst: Py = Py::None;
    let t6525: Result<Option<Py>, PyErr> = (|| -> Result<Option<Py>, PyErr> {
            v_want = { let t6524 = v_inputs.clone(); f_simpson(py::getitem(&t6524, &Py::str("counts"))?, py::getitem(&t6524, &Py::str("estimator"))?, py::method(&t6524, "get", &[Py::str("swap"), Py::Bool(false)])?)? };
            #[allow(unreachable_code)]
            Ok(None)
        })();
    match t6525 {
        Ok(Some(r)) => { return Ok(r); }
        Ok(None) => {}
        Err(e) if (!matches!(e, PyErr::Unsupported(_))) => {
            return Ok(Py::Float(f64::INFINITY));
        }
        Err(e) => return Err(e),
    }
    v_worst = Py::Float(0.0);
    for t6526 in py::iter(&(py::method(&v_want.clone(), "items", &[])?))? {
        let t6563 = (t6526).items()?;
        if t6563.len() != 2 { return Err(PyErr::Value("unpack".into())); }
        v_k = t6563[0].clone();
        v_v = t6563[1].clone();
        if ({ let t6527 = v_k.clone(); Py::Bool({ let t6528 = v_outputs.clone(); if !py::contains(&t6528, &t6527)? { true } else { false } }) }).truthy() {
            return Ok(Py::Float(f64::INFINITY));
        }
        v_worst = py::minmax(&[v_worst.clone(), f__rel(py::getitem(&v_outputs.clone(), &v_k.clone())?, v_v.clone())?], true)?;
    }
    return Ok(v_worst.clone());
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn s_6602(mut v_inputs: Py, mut v_outputs: Py) -> R {
    let mut v_k: Py = Py::None;
    let mut v_v: Py = Py::None;
    let mut v_want: Py = Py::None;
    let mut v_worst: Py = Py::None;
    let t6604: Result<Option<Py>, PyErr> = (|| -> Result<Option<Py>, PyErr> {
            v_want = { let t6603 = v_inputs.clone(); f_species_area(py::getitem(&t6603, &Py::str("A"))?, py::getitem(&t6603, &Py::str("S"))?, py::method(&t6603, "get", &[Py::str("no_exp"), Py::Bool(false)])?)? };
            #[allow(unreachable_code)]
            Ok(None)
        })();
    match t6604 {
        Ok(Some(r)) => { return Ok(r); }
        Ok(None) => {}
        Err(e) if (!matches!(e, PyErr::Unsupported(_))) => {
            return Ok(Py::Float(f64::INFINITY));
        }
        Err(e) => return Err(e),
    }
    v_worst = Py::Float(0.0);
    for t6605 in py::iter(&(py::method(&v_want.clone(), "items", &[])?))? {
        let t6642 = (t6605).items()?;
        if t6642.len() != 2 { return Err(PyErr::Value("unpack".into())); }
        v_k = t6642[0].clone();
        v_v = t6642[1].clone();
        if ({ let t6606 = v_k.clone(); Py::Bool({ let t6607 = v_outputs.clone(); if !py::contains(&t6607, &t6606)? { true } else { false } }) }).truthy() {
            return Ok(Py::Float(f64::INFINITY));
        }
        v_worst = py::minmax(&[v_worst.clone(), f__rel(py::getitem(&v_outputs.clone(), &v_k.clone())?, v_v.clone())?], true)?;
    }
    return Ok(v_worst.clone());
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn s_6680(mut v_inputs: Py, mut v_outputs: Py) -> R {
    let mut v_k: Py = Py::None;
    let mut v_v: Py = Py::None;
    let mut v_want: Py = Py::None;
    let mut v_worst: Py = Py::None;
    let t6682: Result<Option<Py>, PyErr> = (|| -> Result<Option<Py>, PyErr> {
            v_want = { let t6681 = v_inputs.clone(); f_tajima_d(py::getitem(&t6681, &Py::str("seqs"))?, py::method(&t6681, "get", &[Py::str("sqrt_denominator"), Py::Bool(true)])?)? };
            #[allow(unreachable_code)]
            Ok(None)
        })();
    match t6682 {
        Ok(Some(r)) => { return Ok(r); }
        Ok(None) => {}
        Err(e) if (!matches!(e, PyErr::Unsupported(_))) => {
            return Ok(Py::Float(f64::INFINITY));
        }
        Err(e) => return Err(e),
    }
    v_worst = Py::Float(0.0);
    for t6683 in py::iter(&(py::method(&v_want.clone(), "items", &[])?))? {
        let t6720 = (t6683).items()?;
        if t6720.len() != 2 { return Err(PyErr::Value("unpack".into())); }
        v_k = t6720[0].clone();
        v_v = t6720[1].clone();
        if ({ let t6684 = v_k.clone(); Py::Bool({ let t6685 = v_outputs.clone(); if !py::contains(&t6685, &t6684)? { true } else { false } }) }).truthy() {
            return Ok(Py::Float(f64::INFINITY));
        }
        v_worst = py::minmax(&[v_worst.clone(), f__rel(py::getitem(&v_outputs.clone(), &v_k.clone())?, v_v.clone())?], true)?;
    }
    return Ok(v_worst.clone());
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn s_6749(mut v_inputs: Py, mut v_outputs: Py) -> R {
    let mut v_k: Py = Py::None;
    let mut v_v: Py = Py::None;
    let mut v_want: Py = Py::None;
    let mut v_worst: Py = Py::None;
    let t6751: Result<Option<Py>, PyErr> = (|| -> Result<Option<Py>, PyErr> {
            v_want = { let t6750 = v_inputs.clone(); f_testcross_recombination(py::getitem(&t6750, &Py::str("parental"))?, py::getitem(&t6750, &Py::str("recombinant"))?)? };
            #[allow(unreachable_code)]
            Ok(None)
        })();
    match t6751 {
        Ok(Some(r)) => { return Ok(r); }
        Ok(None) => {}
        Err(e) if (!matches!(e, PyErr::Unsupported(_))) => {
            return Ok(Py::Float(f64::INFINITY));
        }
        Err(e) => return Err(e),
    }
    v_worst = Py::Float(0.0);
    for t6752 in py::iter(&(py::method(&v_want.clone(), "items", &[])?))? {
        let t6789 = (t6752).items()?;
        if t6789.len() != 2 { return Err(PyErr::Value("unpack".into())); }
        v_k = t6789[0].clone();
        v_v = t6789[1].clone();
        if ({ let t6753 = v_k.clone(); Py::Bool({ let t6754 = v_outputs.clone(); if !py::contains(&t6754, &t6753)? { true } else { false } }) }).truthy() {
            return Ok(Py::Float(f64::INFINITY));
        }
        v_worst = py::minmax(&[v_worst.clone(), f__rel(py::getitem(&v_outputs.clone(), &v_k.clone())?, v_v.clone())?], true)?;
    }
    return Ok(v_worst.clone());
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn s_6893(mut v_inputs: Py, mut v_outputs: Py) -> R {
    let mut v_k: Py = Py::None;
    let mut v_v: Py = Py::None;
    let mut v_want: Py = Py::None;
    let mut v_worst: Py = Py::None;
    let t6895: Result<Option<Py>, PyErr> = (|| -> Result<Option<Py>, PyErr> {
            v_want = { let t6894 = v_inputs.clone(); f_three_point_testcross(py::getitem(&t6894, &Py::str("classes"))?, py::getitem(&t6894, &Py::str("counts"))?)? };
            #[allow(unreachable_code)]
            Ok(None)
        })();
    match t6895 {
        Ok(Some(r)) => { return Ok(r); }
        Ok(None) => {}
        Err(e) if (!matches!(e, PyErr::Unsupported(_))) => {
            return Ok(Py::Float(f64::INFINITY));
        }
        Err(e) => return Err(e),
    }
    v_worst = Py::Float(0.0);
    for t6896 in py::iter(&(py::method(&v_want.clone(), "items", &[])?))? {
        let t6933 = (t6896).items()?;
        if t6933.len() != 2 { return Err(PyErr::Value("unpack".into())); }
        v_k = t6933[0].clone();
        v_v = t6933[1].clone();
        if ({ let t6897 = v_k.clone(); Py::Bool({ let t6898 = v_outputs.clone(); if !py::contains(&t6898, &t6897)? { true } else { false } }) }).truthy() {
            return Ok(Py::Float(f64::INFINITY));
        }
        v_worst = py::minmax(&[v_worst.clone(), f__rel(py::getitem(&v_outputs.clone(), &v_k.clone())?, v_v.clone())?], true)?;
    }
    return Ok(v_worst.clone());
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn s_7000(mut v_inputs: Py, mut v_outputs: Py) -> R {
    let mut v_k: Py = Py::None;
    let mut v_v: Py = Py::None;
    let mut v_want: Py = Py::None;
    let mut v_worst: Py = Py::None;
    let t7002: Result<Option<Py>, PyErr> = (|| -> Result<Option<Py>, PyErr> {
            v_want = { let t7001 = v_inputs.clone(); f_tm_nn(py::getitem(&t7001, &Py::str("seq"))?, py::getitem(&t7001, &Py::str("Na"))?, py::getitem(&t7001, &Py::str("dnac1"))?, py::getitem(&t7001, &Py::str("dnac2"))?)? };
            #[allow(unreachable_code)]
            Ok(None)
        })();
    match t7002 {
        Ok(Some(r)) => { return Ok(r); }
        Ok(None) => {}
        Err(e) if (!matches!(e, PyErr::Unsupported(_))) => {
            return Ok(Py::Float(f64::INFINITY));
        }
        Err(e) => return Err(e),
    }
    v_worst = Py::Float(0.0);
    for t7003 in py::iter(&(py::method(&v_want.clone(), "items", &[])?))? {
        let t7040 = (t7003).items()?;
        if t7040.len() != 2 { return Err(PyErr::Value("unpack".into())); }
        v_k = t7040[0].clone();
        v_v = t7040[1].clone();
        if ({ let t7004 = v_k.clone(); Py::Bool({ let t7005 = v_outputs.clone(); if !py::contains(&t7005, &t7004)? { true } else { false } }) }).truthy() {
            return Ok(Py::Float(f64::INFINITY));
        }
        v_worst = py::minmax(&[v_worst.clone(), f__rel(py::getitem(&v_outputs.clone(), &v_k.clone())?, v_v.clone())?], true)?;
    }
    return Ok(v_worst.clone());
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn s_7064(mut v_inputs: Py, mut v_outputs: Py) -> R {
    let mut v_k: Py = Py::None;
    let mut v_v: Py = Py::None;
    let mut v_want: Py = Py::None;
    let mut v_worst: Py = Py::None;
    let t7066: Result<Option<Py>, PyErr> = (|| -> Result<Option<Py>, PyErr> {
            v_want = { let t7065 = v_inputs.clone(); f_tm_wallace(py::getitem(&t7065, &Py::str("seq"))?)? };
            #[allow(unreachable_code)]
            Ok(None)
        })();
    match t7066 {
        Ok(Some(r)) => { return Ok(r); }
        Ok(None) => {}
        Err(e) if (!matches!(e, PyErr::Unsupported(_))) => {
            return Ok(Py::Float(f64::INFINITY));
        }
        Err(e) => return Err(e),
    }
    v_worst = Py::Float(0.0);
    for t7067 in py::iter(&(py::method(&v_want.clone(), "items", &[])?))? {
        let t7104 = (t7067).items()?;
        if t7104.len() != 2 { return Err(PyErr::Value("unpack".into())); }
        v_k = t7104[0].clone();
        v_v = t7104[1].clone();
        if ({ let t7068 = v_k.clone(); Py::Bool({ let t7069 = v_outputs.clone(); if !py::contains(&t7069, &t7068)? { true } else { false } }) }).truthy() {
            return Ok(Py::Float(f64::INFINITY));
        }
        v_worst = py::minmax(&[v_worst.clone(), f__rel(py::getitem(&v_outputs.clone(), &v_k.clone())?, v_v.clone())?], true)?;
    }
    return Ok(v_worst.clone());
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn s_7132(mut v_inputs: Py, mut v_outputs: Py) -> R {
    let mut v_k: Py = Py::None;
    let mut v_v: Py = Py::None;
    let mut v_want: Py = Py::None;
    let mut v_worst: Py = Py::None;
    let t7134: Result<Option<Py>, PyErr> = (|| -> Result<Option<Py>, PyErr> {
            v_want = { let t7133 = v_inputs.clone(); f_transcription(py::getitem(&t7133, &Py::str("seq"))?, py::getitem(&t7133, &Py::str("strand"))?)? };
            #[allow(unreachable_code)]
            Ok(None)
        })();
    match t7134 {
        Ok(Some(r)) => { return Ok(r); }
        Ok(None) => {}
        Err(e) if (!matches!(e, PyErr::Unsupported(_))) => {
            return Ok(Py::Float(f64::INFINITY));
        }
        Err(e) => return Err(e),
    }
    v_worst = Py::Float(0.0);
    for t7135 in py::iter(&(py::method(&v_want.clone(), "items", &[])?))? {
        let t7172 = (t7135).items()?;
        if t7172.len() != 2 { return Err(PyErr::Value("unpack".into())); }
        v_k = t7172[0].clone();
        v_v = t7172[1].clone();
        if ({ let t7136 = v_k.clone(); Py::Bool({ let t7137 = v_outputs.clone(); if !py::contains(&t7137, &t7136)? { true } else { false } }) }).truthy() {
            return Ok(Py::Float(f64::INFINITY));
        }
        v_worst = py::minmax(&[v_worst.clone(), f__rel(py::getitem(&v_outputs.clone(), &v_k.clone())?, v_v.clone())?], true)?;
    }
    return Ok(v_worst.clone());
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn s_7206(mut v_inputs: Py, mut v_outputs: Py) -> R {
    let mut v_k: Py = Py::None;
    let mut v_v: Py = Py::None;
    let mut v_want: Py = Py::None;
    let mut v_worst: Py = Py::None;
    let t7208: Result<Option<Py>, PyErr> = (|| -> Result<Option<Py>, PyErr> {
            v_want = { let t7207 = v_inputs.clone(); f_translation(py::getitem(&t7207, &Py::str("seq"))?, py::getitem(&t7207, &Py::str("stops"))?)? };
            #[allow(unreachable_code)]
            Ok(None)
        })();
    match t7208 {
        Ok(Some(r)) => { return Ok(r); }
        Ok(None) => {}
        Err(e) if (!matches!(e, PyErr::Unsupported(_))) => {
            return Ok(Py::Float(f64::INFINITY));
        }
        Err(e) => return Err(e),
    }
    v_worst = Py::Float(0.0);
    for t7209 in py::iter(&(py::method(&v_want.clone(), "items", &[])?))? {
        let t7246 = (t7209).items()?;
        if t7246.len() != 2 { return Err(PyErr::Value("unpack".into())); }
        v_k = t7246[0].clone();
        v_v = t7246[1].clone();
        if ({ let t7210 = v_k.clone(); Py::Bool({ let t7211 = v_outputs.clone(); if !py::contains(&t7211, &t7210)? { true } else { false } }) }).truthy() {
            return Ok(Py::Float(f64::INFINITY));
        }
        v_worst = py::minmax(&[v_worst.clone(), f__rel(py::getitem(&v_outputs.clone(), &v_k.clone())?, v_v.clone())?], true)?;
    }
    return Ok(v_worst.clone());
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn s_7259(mut v_inputs: Py, mut v_outputs: Py) -> R {
    let mut v_k: Py = Py::None;
    let mut v_v: Py = Py::None;
    let mut v_want: Py = Py::None;
    let mut v_worst: Py = Py::None;
    let t7261: Result<Option<Py>, PyErr> = (|| -> Result<Option<Py>, PyErr> {
            v_want = { let t7260 = v_inputs.clone(); f_tree_count(py::getitem(&t7260, &Py::str("n"))?, py::getitem(&t7260, &Py::str("rooted"))?)? };
            #[allow(unreachable_code)]
            Ok(None)
        })();
    match t7261 {
        Ok(Some(r)) => { return Ok(r); }
        Ok(None) => {}
        Err(e) if (!matches!(e, PyErr::Unsupported(_))) => {
            return Ok(Py::Float(f64::INFINITY));
        }
        Err(e) => return Err(e),
    }
    v_worst = Py::Float(0.0);
    for t7262 in py::iter(&(py::method(&v_want.clone(), "items", &[])?))? {
        let t7299 = (t7262).items()?;
        if t7299.len() != 2 { return Err(PyErr::Value("unpack".into())); }
        v_k = t7299[0].clone();
        v_v = t7299[1].clone();
        if ({ let t7263 = v_k.clone(); Py::Bool({ let t7264 = v_outputs.clone(); if !py::contains(&t7264, &t7263)? { true } else { false } }) }).truthy() {
            return Ok(Py::Float(f64::INFINITY));
        }
        v_worst = py::minmax(&[v_worst.clone(), f__rel(py::getitem(&v_outputs.clone(), &v_k.clone())?, v_v.clone())?], true)?;
    }
    return Ok(v_worst.clone());
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn s_7360(mut v_inputs: Py, mut v_outputs: Py) -> R {
    let mut v_k: Py = Py::None;
    let mut v_v: Py = Py::None;
    let mut v_want: Py = Py::None;
    let mut v_worst: Py = Py::None;
    let t7362: Result<Option<Py>, PyErr> = (|| -> Result<Option<Py>, PyErr> {
            v_want = { let t7361 = v_inputs.clone(); f_upgma(py::getitem(&t7361, &Py::str("D"))?, py::method(&t7361, "get", &[Py::str("linkage"), Py::str("average")])?)? };
            #[allow(unreachable_code)]
            Ok(None)
        })();
    match t7362 {
        Ok(Some(r)) => { return Ok(r); }
        Ok(None) => {}
        Err(e) if (!matches!(e, PyErr::Unsupported(_))) => {
            return Ok(Py::Float(f64::INFINITY));
        }
        Err(e) => return Err(e),
    }
    v_worst = Py::Float(0.0);
    for t7363 in py::iter(&(py::method(&v_want.clone(), "items", &[])?))? {
        let t7400 = (t7363).items()?;
        if t7400.len() != 2 { return Err(PyErr::Value("unpack".into())); }
        v_k = t7400[0].clone();
        v_v = t7400[1].clone();
        if ({ let t7364 = v_k.clone(); Py::Bool({ let t7365 = v_outputs.clone(); if !py::contains(&t7365, &t7364)? { true } else { false } }) }).truthy() {
            return Ok(Py::Float(f64::INFINITY));
        }
        v_worst = py::minmax(&[v_worst.clone(), f__rel(py::getitem(&v_outputs.clone(), &v_k.clone())?, v_v.clone())?], true)?;
    }
    return Ok(v_worst.clone());
    #[allow(unreachable_code)]
    Ok(Py::None)
}

#[allow(unused_mut, unused_variables, unused_assignments, non_snake_case, unreachable_code)]
fn s_7429(mut v_inputs: Py, mut v_outputs: Py) -> R {
    let mut v_k: Py = Py::None;
    let mut v_v: Py = Py::None;
    let mut v_want: Py = Py::None;
    let mut v_worst: Py = Py::None;
    let t7431: Result<Option<Py>, PyErr> = (|| -> Result<Option<Py>, PyErr> {
            v_want = { let t7430 = v_inputs.clone(); f_watterson(py::getitem(&t7430, &Py::str("seqs"))?, py::method(&t7430, "get", &[Py::str("harmonic_to_n"), Py::Bool(false)])?)? };
            #[allow(unreachable_code)]
            Ok(None)
        })();
    match t7431 {
        Ok(Some(r)) => { return Ok(r); }
        Ok(None) => {}
        Err(e) if (!matches!(e, PyErr::Unsupported(_))) => {
            return Ok(Py::Float(f64::INFINITY));
        }
        Err(e) => return Err(e),
    }
    v_worst = Py::Float(0.0);
    for t7432 in py::iter(&(py::method(&v_want.clone(), "items", &[])?))? {
        let t7469 = (t7432).items()?;
        if t7469.len() != 2 { return Err(PyErr::Value("unpack".into())); }
        v_k = t7469[0].clone();
        v_v = t7469[1].clone();
        if ({ let t7433 = v_k.clone(); Py::Bool({ let t7434 = v_outputs.clone(); if !py::contains(&t7434, &t7433)? { true } else { false } }) }).truthy() {
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
