//! Computational Biology: 15 calculations through the Rust library (the same engine as the desktop app, the
//! command line and the browser). Run from the rust/ folder:
//!     cargo run --release --example readme -p computational-biology-ai999x
//! Inputs are text, read like a form field ("2 m", "0.5", "x^3 - 3x"); a JSON list stays a list.

use serde_json::{json, Value};

/// One calculation: the request in, its outcome and every value with its unit out.
fn show(request: Value) {
    let answer: Value = serde_json::from_str(&computational_biology_ai999x::solve_json(&request.to_string())).unwrap();
    println!("{}", answer["outcome"].as_str().unwrap_or(""));
    for v in answer["values"].as_array().into_iter().flatten() {
        let unit = v[2].as_str().unwrap_or("");
        let unit = if unit.is_empty() || unit == "1" { String::new() } else { format!(" {unit}") };
        println!("{} = {}{}", v[0].as_str().unwrap_or(""), v[1], unit);
    }
}

fn main() {
    show(json!({"key": "bio_global_alignment_score", "inputs": {"seq1": "GATTACA", "seq2": "GCATGCU", "match": "1", "mismatch": "-1", "gap": "-1"}}));
    println!();
    show(json!({"key": "bio_ace", "inputs": {"counts": [1, 1, 2, 3, 5, 8, 20, 50], "rare_threshold": "10"}}));
    println!();
    show(json!({"key": "bio_competitive_inhibition_fit", "inputs": {"S": [1, 2, 4, 8, 1, 2, 4, 8], "I": [0, 0, 0, 0, 2, 2, 2, 2], "v": [5, 6.66667, 8, 8.88889, 2.5, 4, 5.71429, 7.27273]}}));
    println!();
    show(json!({"key": "bio_pedigree_inbreeding", "inputs": {"pedigree": "P-S, Q-S, P-D, Q-D, S-X, D-X", "individual": "X"}}));
    println!();
    show(json!({"key": "bio_chemostat_steady_state", "inputs": {"mu_max": "1", "Ks": "1", "Y": "0.5", "D": "0.5", "S0": "10"}}));
    println!();
    show(json!({"key": "bio_jukes_cantor_distance", "inputs": {"seq1": "AAAAAAAAAA", "seq2": "AAAAAAAACC"}}));
    println!();
    show(json!({"key": "bio_leslie_growth_rate", "inputs": {"F": [0, 2], "S": [0.5]}}));
    println!();
    show(json!({"key": "bio_haploid_selection", "inputs": {"p0": "0.5", "s": "1", "t": "1"}}));
    println!();
    show(json!({"key": "bio_amino_acid_composition", "inputs": {"seq": "MAAAK"}}));
    println!();
    show(json!({"key": "bio_flux_balance_analysis", "inputs": {"S": [[1, -1, 0, 0], [0, 1, -1, -1]], "lb": [0, 0, 0, 0], "ub": [10, 8, 100, 3], "c": [0, 0, 1, 0]}}));
    println!();
    show(json!({"key": "bio_global_alignment_affine", "inputs": {"seq1": "ACGTTTACG", "seq2": "ACGACG", "match": "2", "mismatch": "-1", "open": "-5", "extend": "-1"}}));
    println!();
    show(json!({"key": "bio_dominance_indices", "inputs": {"counts": [50, 30, 20]}}));
    println!();
    show(json!({"key": "bio_eadie_hofstee", "inputs": {"S": [1, 2, 4, 8], "v": [5, 6.66667, 8, 8.88889]}}));
    println!();
    show(json!({"key": "bio_map_distance", "inputs": {"r": "0.2", "map": "haldane"}}));
    println!();
    show(json!({"key": "bio_kimura_distance", "inputs": {"seq1": "AAAAAAAAAA", "seq2": "AAAAAAAGGC"}}));
}
