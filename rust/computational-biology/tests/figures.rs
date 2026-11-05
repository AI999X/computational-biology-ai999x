//! Every Computational Biology capability's own live figures (src/figures*.rs): its worked example gives each
//! published figure kind, a changed input changes the drawing, and no drawing holds a non-finite number or text outside
//! its canvas.

use ax_engine::figures::{figures, Figure};
use ax_engine::Input;
use computational_biology_ai999x::pack;

/// The published figure kinds of every capability (publish/website/models/computational-biology-ai999x/figures,
/// except the generic "curve"), as tools/figure_parity.py lists them.
const KINDS: &[(&str, &[&str])] = &[
    ("bio_ace", &["bars"]),
    ("bio_amino_acid_composition", &["bars"]),
    ("bio_aromaticity", &["example", "sequence"]),
    ("bio_beta_diversity", &["bars"]),
    ("bio_chao1", &["bars"]),
    ("bio_chemostat_steady_state", &["operating"]),
    ("bio_coalescent_times", &["genealogy"]),
    ("bio_competitive_inhibition_fit", &["fit", "linear"]),
    ("bio_cross_genotype_probability", &["bars", "punnett"]),
    ("bio_cross_phenotype_ratio", &["punnett"]),
    ("bio_drift_heterozygosity", &["heatmap"]),
    ("bio_eadie_hofstee", &["linear"]),
    ("bio_edit_distance", &["table"]),
    ("bio_faith_pd", &["bars", "tree"]),
    ("bio_fitch_parsimony", &["tree"]),
    ("bio_flux_balance_analysis", &["network"]),
    ("bio_fst", &["heatmap"]),
    ("bio_gc_content", &["profile", "sequence"]),
    ("bio_global_alignment_affine", &["table"]),
    ("bio_global_alignment_score", &["table"]),
    ("bio_gravy", &["profile"]),
    ("bio_hamming_distance", &["alignment"]),
    ("bio_hanes_woolf", &["linear"]),
    ("bio_haploid_selection", &["trajectory"]),
    ("bio_hardy_weinberg_exact_test", &["bars"]),
    ("bio_hardy_weinberg_test", &["bars", "triangle"]),
    ("bio_heterozygote_advantage", &["heatmap"]),
    ("bio_hill_fit", &["fit"]),
    ("bio_hill_plot", &["linear"]),
    ("bio_instability_index", &["profile"]),
    ("bio_isoelectric_point", &["titration"]),
    ("bio_kmer_counts", &["sequence"]),
    ("bio_leslie_projection", &["population", "trajectory"]),
    ("bio_lineweaver_burk", &["linear"]),
    ("bio_linkage_disequilibrium", &["table"]),
    ("bio_local_alignment_affine", &["table"]),
    ("bio_local_alignment_score", &["table"]),
    ("bio_longest_orf", &["orf"]),
    ("bio_michaelis_menten_fit", &["fit", "fit-2"]),
    ("bio_michaelis_menten_progress", &["kinetics", "kinetics-2"]),
    ("bio_ne_fluctuating_size", &["bars"]),
    ("bio_ne_sex_ratio", &["field"]),
    ("bio_neighbor_joining", &["matrix", "tree"]),
    ("bio_nucleotide_diversity", &["sample"]),
    ("bio_pedigree_inbreeding", &["bars", "pedigree"]),
    ("bio_protein_extinction", &["heatmap"]),
    ("bio_protein_global_blosum62", &["table"]),
    ("bio_protein_local_blosum62", &["table"]),
    ("bio_protein_mass", &["bars"]),
    ("bio_recessive_lethal", &["trajectory"]),
    ("bio_restriction_sites", &["digest", "gel"]),
    ("bio_reverse_complement", &["duplex"]),
    ("bio_robinson_foulds", &["trees"]),
    ("bio_schumacher_eschmeyer", &["fit"]),
    ("bio_secondary_structure_fraction", &["bars"]),
    ("bio_semiglobal_alignment_score", &["table"]),
    ("bio_species_area", &["fit"]),
    ("bio_tajima_d", &["bars", "sample"]),
    ("bio_transcription", &["sequence"]),
    ("bio_translation", &["codons"]),
    ("bio_upgma", &["matrix", "tree"]),
    ("bio_watterson_theta", &["sample"]),
];

/// Inputs to change for a capability whose inputs are not plain quantities (a second, different example).
const OTHER: &[(&str, &str, &str)] = &[
    ("bio_amino_acid_composition", "seq", "MKKLLA"),
    ("bio_aromaticity", "seq", "MFWAAKY"),
    ("bio_gc_content", "seq", "ATGCATTTGCGA"),
    ("bio_gravy", "seq", "MKTLLILAV"),
    ("bio_instability_index", "seq", "MAEGEITTFTALTEKFNLPPGNYKKPKLLYCSNGGHFLRILPDGTVDGTRDRSDQHIQLQLSAESVGEVYIKST"),
    ("bio_isoelectric_point", "seq", "INGARKK"),
    ("bio_kmer_counts", "seq", "ATATATGG"),
    ("bio_longest_orf", "seq", "GGATGAAACCCTTTTAGCC"),
    ("bio_reverse_complement", "seq", "ATGCCGTTA"),
    ("bio_transcription", "seq", "ATGGCTGAATAA"),
    ("bio_translation", "seq", "ATGGCCATTGTAATGTGA"),
    ("bio_hamming_distance", "seq2", "CATCGTAATGACGGCCA"),
    ("bio_restriction_sites", "seq", "AAGAATTCTTGGATCCGAATTCGG"),
    ("bio_protein_extinction", "seq", "MWYYCAKR"),
    ("bio_protein_mass", "seq", "MAEGEK"),
    ("bio_secondary_structure_fraction", "seq", "EMALKNPGSDV"),
    ("bio_nucleotide_diversity", "seqs", "AAAA AATT ATTT"),
    ("bio_watterson_theta", "seqs", "AAAA AATT ATTT"),
    ("bio_tajima_d", "seqs", "AAAAAA AAAAAT AAATTT AATTTT"),
    ("bio_edit_distance", "seq2", "SITTINGS"),
    ("bio_global_alignment_score", "seq2", "GCATGCA"),
    ("bio_global_alignment_affine", "seq2", "ACGTCG"),
    ("bio_protein_global_blosum62", "seq2", "PAWHEA"),
    ("bio_semiglobal_alignment_score", "seq2", "ACGTA"),
    ("bio_local_alignment_score", "seq2", "GGTTGACT"),
    ("bio_local_alignment_affine", "seq2", "GGACGTTACG"),
    ("bio_protein_local_blosum62", "seq2", "PAWHEAEW"),
    ("bio_faith_pd", "taxa", "A C"),
    ("bio_fitch_parsimony", "states", "A:AACC B:AATG C:GACT D:GGCT"),
    ("bio_neighbor_joining", "D", "[[0, 5, 9, 9, 8], [5, 0, 10, 10, 9], [9, 10, 0, 8, 7], [9, 10, 8, 0, 5], [8, 9, 7, 5, 0]]"),
    ("bio_upgma", "D", "[[0, 2, 6, 6], [2, 0, 6, 6], [6, 6, 0, 4], [6, 6, 4, 0]]"),
    ("bio_robinson_foulds", "tree2", "((A,B),(C,E),(D,F));"),
    ("bio_pedigree_inbreeding", "pedigree", "P-S, Q-S, P-D, R-D, S-X, D-X"),
    ("bio_cross_genotype_probability", "offspring", "AA"),
    ("bio_cross_phenotype_ratio", "parent2", "Aabb"),
    ("bio_flux_balance_analysis", "ub", "[10, 6, 100, 3]"),
    ("bio_chao1", "counts", "[1, 1, 2, 2, 2, 5, 10, 12]"),
    ("bio_ace", "counts", "[1, 2, 2, 3, 5, 8, 20, 60]"),
    ("bio_beta_diversity", "y", "[8, 1, 6, 3]"),
    ("bio_ne_fluctuating_size", "sizes", "[1000, 50, 800]"),
    ("bio_fst", "p", "[0.3, 0.7]"),
    ("bio_leslie_projection", "n0", "[80, 60, 30]"),
    ("bio_michaelis_menten_fit", "v", "[4.8, 6.7, 8.1, 8.8]"),
    ("bio_lineweaver_burk", "v", "[4.8, 6.7, 8.1, 8.8]"),
    ("bio_eadie_hofstee", "v", "[4.8, 6.7, 8.1, 8.8]"),
    ("bio_hanes_woolf", "v", "[4.8, 6.7, 8.1, 8.8]"),
    ("bio_hill_fit", "v", "[0.6, 2.1, 4.9, 8.1, 9.3]"),
    ("bio_hill_plot", "Y", "[0.06, 0.21, 0.49, 0.81]"),
    ("bio_competitive_inhibition_fit", "v", "[5, 6.7, 8, 8.9, 2.4, 4.1, 5.6, 7.3]"),
    ("bio_schumacher_eschmeyer", "R", "[0, 3, 4]"),
    ("bio_species_area", "S", "[10, 22, 41, 85]"),
    ("bio_hardy_weinberg_test", "n_Aa", "400"),
    ("bio_hardy_weinberg_exact_test", "n_Aa", "2"),
    ("bio_linkage_disequilibrium", "n_AB", "30"),
];

/// A capability refused by the model for these inputs gives no figure (capability, inputs).
const REFUSED: (&str, &[(&str, &str)]) = ("bio_gc_content", &[("seq", "")]);


/// The arguments of a structured example: solve('key', a='1 V', b=8, c="(A·B)'").
fn example(key: &str) -> Vec<(String, Input)> {
    let p = pack();
    let ex = &p.procedure(key).expect(key).structured_example;
    let cs: Vec<char> = ex.chars().collect();
    let mut i = ex.find(',').map(|x| ex[..x].chars().count() + 1).unwrap_or(cs.len());
    let mut out = vec![];
    while i < cs.len() {
        while i < cs.len() && (cs[i] == ' ' || cs[i] == ',') {
            i += 1;
        }
        if i >= cs.len() || cs[i] == ')' {
            break;
        }
        let st = i;
        while cs[i] != '=' {
            i += 1;
        }
        let name: String = cs[st..i].iter().collect();
        i += 1;
        let value: String;
        if cs[i] == '\'' || cs[i] == '"' {
            let qch = cs[i];
            let s0 = i + 1;
            i += 1;
            while cs[i] != qch {
                i += 1;
            }
            value = cs[s0..i].iter().collect();
            i += 1;
        } else if cs[i] == '[' {
            let s0 = i;
            let mut depth = 0;
            loop {
                if cs[i] == '[' {
                    depth += 1;
                } else if cs[i] == ']' {
                    depth -= 1;
                    if depth == 0 {
                        break;
                    }
                }
                i += 1;
            }
            i += 1;
            value = cs[s0..i].iter().collect();
        } else {
            let s0 = i;
            while i < cs.len() && cs[i] != ',' && cs[i] != ')' {
                i += 1;
            }
            value = cs[s0..i].iter().collect();
        }
        out.push((name.trim().to_string(), Input::Text(value.trim().to_string())));
    }
    out
}

fn check_texts(f: &Figure) {
    if let Err(e) = ax_engine::figkit::check_figure(f) {
        panic!("{e}");
    }
}

#[test]
fn every_capability_draws_its_published_kinds() {
    let p = pack();
    let mut missing = vec![];
    for (key, kinds) in KINDS {
        let given = example(key);
        assert!(ax_engine::solve(p, key, &given, None).solved(), "{key}: the worked example is not solved");
        let figs = figures(p, key, &given, None, None);
        for kind in *kinds {
            match figs.iter().find(|f| f.kind == *kind) {
                Some(f) => {
                    assert_eq!(f.slug, format!("{}-{kind}", key.replace('_', "-")));
                    assert!(!f.svg.contains("NaN") && !f.svg.contains("inf"), "{key} {kind}: a non-finite number");
                    check_texts(f);
                }
                None => missing.push(format!("{key}: {kind}")),
            }
        }
        if let Some(i) = figs.iter().position(|f| f.slug.contains("-live-")) {
            assert!(figs[i..].iter().all(|f| f.slug.contains("-live-")), "{key}: generic figures before own ones");
        }
    }
    assert!(missing.is_empty(), "missing figures: {missing:?}");
}

/// The numbers of a text multiplied by f ("12 V" -> "16.44 V"; "[3, 6]" -> "[4.11, 8.22]"; "8e-05 m^3/s" ->
/// "0.0001096 m^3/s"): a number that follows a letter or "^" (a unit's power) is kept.
fn scale_text(t: &str, f: f64) -> Option<String> {
    let cs: Vec<char> = t.chars().collect();
    let mut out = String::new();
    let mut changed = false;
    let mut i = 0;
    while i < cs.len() {
        let prev = if i == 0 { ' ' } else { cs[i - 1] };
        let starts = cs[i].is_ascii_digit() || (cs[i] == '-' && i + 1 < cs.len() && cs[i + 1].is_ascii_digit());
        if starts && !(prev.is_alphanumeric() || prev == '^' || prev == '.' || prev == '_') {
            let s0 = i;
            i += 1;
            while i < cs.len() && (cs[i].is_ascii_digit() || cs[i] == '.') {
                i += 1;
            }
            if i + 1 < cs.len() && (cs[i] == 'e' || cs[i] == 'E')
                && (cs[i + 1].is_ascii_digit() || ((cs[i + 1] == '-' || cs[i + 1] == '+') && i + 2 < cs.len())) {
                i += 2;
                while i < cs.len() && cs[i].is_ascii_digit() {
                    i += 1;
                }
            }
            let tok: String = cs[s0..i].iter().collect();
            match tok.parse::<f64>() {
                Ok(x) if x != 0.0 => {
                    out.push_str(&format!("{}", x * f));
                    changed = true;
                }
                _ => out.push_str(&tok),
            }
        } else {
            out.push(cs[i]);
            i += 1;
        }
    }
    changed.then_some(out)
}

#[test]
fn a_changed_input_changes_the_drawing() {
    let p = pack();
    let mut unchanged = vec![];
    for (key, kinds) in KINDS {
        let given = example(key);
        let base = figures(p, key, &given, None, None);
        let mut tries: Vec<Vec<(String, Input)>> = vec![];
        for (_, sym, val) in OTHER.iter().filter(|o| o.0 == *key) {
            let mut g = given.clone();
            for e in g.iter_mut() {
                if e.0 == *sym {
                    e.1 = Input::Text(val.to_string());
                }
            }
            tries.push(g);
        }
        for f in [1.37, 0.73, 2.0] {
            for i in 0..given.len() {
                if let Input::Text(t) = &given[i].1 {
                    if let Some(t2) = scale_text(t, f) {
                        let mut g = given.clone();
                        g[i].1 = Input::Text(t2);
                        tries.push(g);
                    }
                }
            }
        }
        let mut ok = false;
        for g in tries {
            let a = ax_engine::solve(p, key, &g, None);
            if !a.solved() {
                continue;
            }
            let figs = figures(p, key, &g, None, None);
            let differs = kinds.iter().all(|kind| {
                match (base.iter().find(|f| f.kind == *kind), figs.iter().find(|f| f.kind == *kind)) {
                    (Some(a), Some(b)) => a.svg != b.svg,
                    _ => false,
                }
            });
            if differs {
                ok = true;
                for f in figs.iter().filter(|f| !f.slug.contains("-live-")) {
                    assert!(!f.svg.contains("NaN"), "{key}: a non-finite number");
                    check_texts(f);
                }
                break;
            }
        }
        if !ok {
            unchanged.push(*key);
        }
    }
    assert!(unchanged.is_empty(), "drawing did not follow a changed input: {unchanged:?}");
}

#[test]
fn a_refused_input_draws_nothing() {
    let p = pack();
    let (key, inputs) = REFUSED;
    let given: Vec<(String, Input)> = inputs.iter().map(|(k, v)| (k.to_string(), Input::Text(v.to_string()))).collect();
    let a = ax_engine::solve(p, key, &given, None);
    assert!(!a.solved(), "{key}: expected a refusal");
    assert!(figures(p, key, &given, None, None).is_empty());
}

/// Writes every worked example's own figures (and those of the OTHER inputs, prefixed "x-") as SVG files into
/// $FIG_DUMP, for a contact sheet: cargo test --test figures -- --ignored dump
#[test]
#[ignore]
fn dump() {
    let Ok(dir) = std::env::var("FIG_DUMP") else { return };
    let p = pack();
    for (key, _) in KINDS {
        let given = example(key);
        for f in figures(p, key, &given, None, None).iter().filter(|f| !f.slug.contains("-live-")) {
            std::fs::write(format!("{dir}/{}.svg", f.slug), &f.svg).unwrap();
        }
        for (_, sym, val) in OTHER.iter().filter(|o| o.0 == *key) {
            let mut g = given.clone();
            for e in g.iter_mut() {
                if e.0 == *sym {
                    e.1 = Input::Text(val.to_string());
                }
            }
            for f in figures(p, key, &g, None, None).iter().filter(|f| !f.slug.contains("-live-")) {
                std::fs::write(format!("{dir}/x-{}.svg", f.slug), &f.svg).unwrap();
            }
        }
    }
}
