//! Alignment tables: cell (i, j) is the model's score (or distance) of the first i letters of one sequence against
//! the first j letters of the other, each one run of the model; the edge row and column follow from the inputs.

use crate::figs_bio::{g, heat_table, run, ti};
use ax_engine::figkit::Solved;
use ax_engine::figures::Figure;

#[derive(Clone, Copy, PartialEq)]
pub enum Mode {
    Edit,
    GlobalLinear,
    GlobalAffine,
    SemiGlobal,
    Local,
}

pub fn dp_table(s: &Solved, mode: Mode, title: &str, what: &str) -> Option<Figure> {
    let a: Vec<char> = s.txt("seq1")?.trim().chars().collect();
    let b: Vec<char> = s.txt("seq2")?.trim().chars().collect();
    let (n1, n2) = (a.len(), b.len());
    let out = if mode == Mode::Edit { "distance" } else { "score" };
    let total = s.out(out)?;
    if n1 == 0 || n2 == 0 || n1 > 16 || n2 > 16 || n1 * n2 + 2 > s.left() {
        return None;
    }
    let gap = s.inp("gap");
    let edge = |k: usize| -> Option<f64> {
        match mode {
            Mode::Edit => Some(k as f64),
            Mode::GlobalLinear => gap.map(|gp| gp * k as f64),
            Mode::SemiGlobal | Mode::Local => Some(0.0),
            Mode::GlobalAffine => None,
        }
    };
    let with_edge = mode != Mode::GlobalAffine;
    let off = if with_edge { 1 } else { 0 };
    let mut m: Vec<Vec<Option<f64>>> = vec![vec![None; n2 + off]; n1 + off];
    if with_edge {
        for j in 0..=n2 {
            m[0][j] = edge(j);
        }
        for (i, row) in m.iter_mut().enumerate() {
            row[0] = edge(i);
        }
    }
    for i in 1..=n1 {
        for j in 1..=n2 {
            let pa: String = a[..i].iter().collect();
            let pb: String = b[..j].iter().collect();
            m[i - 1 + off][j - 1 + off] = run(s, &[("seq1", ti(&pa)), ("seq2", ti(&pb))], out);
        }
    }
    let last = m[n1 - 1 + off][n2 - 1 + off]?;
    if (last - total).abs() > 1e-9 * total.abs().max(1.0) {
        return None;
    }
    // the cells to outline: a traceback through the recurrence, or where the local optimum is first reached
    let mut outline: Vec<(usize, usize)> = vec![];
    match mode {
        Mode::Edit | Mode::GlobalLinear => {
            let (sm, sx, gp) = match mode {
                Mode::Edit => (0.0, 1.0, 1.0),
                _ => (s.inp("match")?, s.inp("mismatch")?, gap?),
            };
            let (mut i, mut j) = (n1, n2);
            outline.push((i, j));
            while i > 0 || j > 0 {
                let here = m[i][j]?;
                let near = |v: Option<f64>, d: f64| v.map(|x| (x + d - here).abs() < 1e-9).unwrap_or(false);
                if i > 0 && j > 0 && near(m[i - 1][j - 1], if a[i - 1] == b[j - 1] { sm } else { sx }) {
                    i -= 1;
                    j -= 1;
                } else if i > 0 && near(m[i - 1][j], gp) {
                    i -= 1;
                } else if j > 0 && near(m[i][j - 1], gp) {
                    j -= 1;
                } else {
                    outline.clear();
                    break;
                }
                outline.push((i, j));
            }
        }
        Mode::Local => {
            'f: for i in 0..m.len() {
                for j in 0..m[0].len() {
                    if m[i][j].map(|v| (v - total).abs() < 1e-9).unwrap_or(false) {
                        outline.push((i, j));
                        break 'f;
                    }
                }
            }
        }
        _ => {}
    }
    let mut rows: Vec<String> = a.iter().map(|c| c.to_string()).collect();
    let mut cols: Vec<String> = b.iter().map(|c| c.to_string()).collect();
    if with_edge {
        rows.insert(0, "-".into());
        cols.insert(0, "-".into());
    }
    let mark = if mode == Mode::Local { outline.first().copied() } else { Some((n1 - 1 + off, n2 - 1 + off)) };
    let ol: Vec<(usize, usize)> = if mode == Mode::Local { vec![] } else { outline.clone() };
    let (body, h) = heat_table(&m, &rows, &cols, "first sequence (prefix length down)", "second sequence (prefix length across)",
                               &format!("{out} = {}", g(total, 6)), mark, &ol, 4);
    let runs = n1 * n2;
    let tail = match mode {
        Mode::Local => format!("The ringed cell is the first where the returned optimum {} is reached: where the best local \
                                alignment ends.", g(total, 6)),
        Mode::Edit | Mode::GlobalLinear if !outline.is_empty() =>
            format!("The outlined cells trace one optimal path back through cells that satisfy the recurrence; the ringed \
                     bottom-right cell is the returned {out} {}.", g(total, 6)),
        _ => format!("The ringed bottom-right cell is the returned {out} {}.", g(total, 6)),
    };
    Some(s.figure("table", title,
                  &format!("{what} of {} against {} as a table: cell (i, j) is the model's {out} of the first i letters of the \
                            first sequence and the first j letters of the second ({runs} runs){}. {tail}",
                           a.iter().collect::<String>(), b.iter().collect::<String>(),
                           match mode {
                               Mode::Edit => "; the top row and left column are the distances to the empty string",
                               Mode::GlobalLinear => "; the top row and left column are the all-gap scores",
                               Mode::SemiGlobal => "; end gaps are free, so the edge row and column are 0",
                               Mode::Local => "; a local score is never below 0, so the edge row and column are 0",
                               Mode::GlobalAffine => "",
                           }), &body, h))
}
