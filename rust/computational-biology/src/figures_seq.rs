//! Sequence figures: residues and bases drawn letter by letter with what the model returns about them, window
//! profiles along a sequence (each window one run of the model), digests and gels, titration of a peptide.

use crate::figs_bio::{g, other, run, ti, fi};
use crate::figs_common::{bars_body, clip, esc, notes, tick, Chart};
use crate::figs_extra::bars_full;
use ax_engine::figkit::{self as k, Solved};
use ax_engine::figures::Figure;

/// Letters in boxes, wrapped at `per` a row, the row's first position on the left; `class(i, c)` shades a letter.
pub fn seq_rows(seq: &[char], per: usize, y0: f64, class: &dyn Fn(usize, char) -> &'static str) -> (String, f64) {
    let per = per.max(1).min(seq.len().max(1));
    let cell = (440.0 / per as f64).min(30.0);
    let x0 = 60.0;
    let mut out = String::new();
    let mut y = y0;
    let wrapped = seq.len() > per;
    for (r, chunk) in seq.chunks(per).enumerate() {
        if wrapped {
            out.push_str(&format!("<text class=\"axp-t\" x=\"{:.1}\" y=\"{:.1}\" text-anchor=\"end\" font-size=\"12\">{}</text>",
                                  x0 - 8.0, y + cell * 0.66, r * per + 1));
        }
        for (j, c) in chunk.iter().enumerate() {
            let i = r * per + j;
            let x = x0 + cell * j as f64;
            let cl = class(i, *c);
            out.push_str(&format!("<rect class=\"{}\" x=\"{x:.1}\" y=\"{y:.1}\" width=\"{:.1}\" height=\"{:.1}\"/>",
                                  if cl.is_empty() { "axp-f" } else { cl }, cell - 2.0, cell - 2.0));
            let fs = (cell * 0.55).clamp(8.0, 15.0);
            out.push_str(&format!("<text class=\"axp-bt\" x=\"{:.1}\" y=\"{:.1}\" text-anchor=\"middle\" font-size=\"{fs:.0}\">{}</text>",
                                  x + (cell - 2.0) / 2.0, y + cell * 0.66, c));
        }
        y += cell + 4.0;
    }
    (out, y)
}

/// The sequence input as upper-case letters (at most `max`).
fn seq_of(s: &Solved, sym: &str, max: usize) -> Option<Vec<char>> {
    let t = s.txt(sym)?;
    let v: Vec<char> = t.chars().filter(|c| c.is_alphabetic()).map(|c| c.to_ascii_uppercase()).collect();
    (!v.is_empty() && v.len() <= max).then_some(v)
}

fn legend(items: &[(&str, &str)], y: f64) -> String {
    let mut out = String::new();
    let mut x = 60.0;
    for (class, name) in items {
        out.push_str(&format!("<rect class=\"{class}\" x=\"{x:.1}\" y=\"{:.1}\" width=\"12\" height=\"12\"/>\
                               <text class=\"axp-lt\" x=\"{:.1}\" y=\"{y:.1}\">{}</text>", y - 10.0, x + 16.0, esc(name)));
        x += name.chars().count() as f64 * 7.4 + 34.0;
    }
    out
}

// ------------------------------------------------------------------------------------------------ aromaticity
pub fn aromatic_example(s: &Solved) -> Option<Figure> {
    let seq = seq_of(s, "seq", 40)?;
    let (ar, n) = (s.out("aromaticity")?, s.out("aromatic_count")?);
    let (body, y) = seq_rows(&seq, 40, 12.0, &|_, c| if "FWY".contains(c) { "axp-b" } else { "" });
    let (nb, y2) = notes(&[("returned:".into(), format!("aromatic_count = {}, aromaticity = {}", g(n, 4), tick(ar)))], y + 16.0);
    Some(s.figure("example", "Aromatic residues of a peptide",
                  &format!("The sequence {} with its aromatic residues F, W and Y shaded; the model returns aromatic_count \
                            = {} and aromaticity = {}.", seq.iter().collect::<String>(), g(n, 4), tick(ar)),
                  &(body + &nb), y2.max(80.0)))
}

pub fn aromatic_sequence(s: &Solved) -> Option<Figure> {
    let seq = seq_of(s, "seq", 300)?;
    let (ar, n) = (s.out("aromaticity")?, s.out("aromatic_count")?);
    let mut out = String::from("<text class=\"axp-mt\" x=\"16\" y=\"22\">Aromatic residues along the sequence</text>");
    out.push_str(&legend(&[("axp-b", "F phenylalanine"), ("axp-b2", "W tryptophan"), ("axp-b3", "Y tyrosine")], 46.0));
    let per = if seq.len() <= 20 { seq.len() } else { 20 };
    let (body, y) = seq_rows(&seq, per, 58.0, &|_, c| match c { 'F' => "axp-b", 'W' => "axp-b2", 'Y' => "axp-b3", _ => "" });
    out.push_str(&body);
    let count = |x: char| seq.iter().filter(|c| **c == x).count();
    let len = seq.len() as f64;
    let lines: Vec<(String, String)> = ['F', 'W', 'Y'].iter()
        .map(|c| (format!("{c}:"), format!("{} of {} residues ({})", count(*c), seq.len(), tick(count(*c) as f64 / len)))).collect();
    let (nb, y2) = notes(&lines, y + 18.0);
    out.push_str(&nb);
    let (nb2, y3) = notes(&[("returned:".into(), format!("aromatic_count = {}, aromaticity = {}", g(n, 4), tick(ar)))], y2);
    out.push_str(&nb2);
    Some(s.figure("sequence", "Aromatic residues along a protein",
                  &format!("The {}-residue sequence written {per} residues per row (row start numbers on the left), F, W and \
                            Y shaded; below, the count of each (read from the sequence) and the aromaticity the model \
                            returns: aromatic_count = {}, aromaticity = {}.", seq.len(), g(n, 4), tick(ar)), &out, y3 + 4.0))
}

// ------------------------------------------------------------------------------------------------ GC content
pub fn gc_sequence(s: &Solved) -> Option<Figure> {
    let seq = seq_of(s, "seq", 300)?;
    let (f, n) = (s.out("gc_fraction")?, s.out("gc_count")?);
    let per = seq.len().min(30);
    let (body, y) = seq_rows(&seq, per, 12.0, &|_, c| if c == 'G' || c == 'C' { "axp-b" } else { "" });
    let (nb, y2) = notes(&[("returned:".into(), format!("gc_count = {}, gc_fraction = {}", g(n, 4), tick(f)))], y + 16.0);
    Some(s.figure("sequence", "GC content of a DNA sequence, base by base",
                  &format!("The sequence written base by base, G and C shaded. The model counts {} G or C bases and returns \
                            the GC fraction {}.", g(n, 4), tick(f)), &(body + &nb), y2.max(90.0)))
}

/// A window profile: each window of `w` letters one run of the model, its output against the window's centre.
fn windows(s: &Solved, seq: &[char], w: usize, out: &str) -> Vec<(f64, f64)> {
    if seq.len() < w || w == 0 {
        return vec![];
    }
    let n = seq.len() - w + 1;
    let step = n.div_ceil(300).max(1);
    let mut v = vec![];
    let mut i = 0;
    while i < n && s.left() > 4 {
        let win: String = seq[i..i + w].iter().collect();
        if let Some(y) = run(s, &[("seq", ti(&win))], out) {
            v.push((i as f64 + 1.0 + (w as f64 - 1.0) / 2.0, y));
        }
        i += step;
    }
    v
}

fn profile(s: &Solved, out: &str, w: usize, title: &str, ylab: &str, guide: Option<(f64, &str)>, what: &str)
           -> Option<Figure> {
    let seq = seq_of(s, "seq", 100000)?;
    let whole = s.out(out)?;
    let pts = windows(s, &seq, w, out);
    if pts.len() < 2 {
        return None;
    }
    let best = pts.iter().cloned().fold(pts[0], |a, b| if b.1 > a.1 { b } else { a });
    let mut ys: Vec<f64> = pts.iter().map(|p| p.1).collect();
    ys.push(whole);
    if let Some((gv, _)) = guide {
        ys.push(gv);
    }
    let mut c = Chart::new((1.0, seq.len() as f64), Chart::range(&ys), (false, false), (false, true),
                           "position in the sequence (window centre)", ylab, &format!("whole sequence: {}", tick(whole)));
    if let Some((gv, gl)) = guide {
        c.hguide(gv, gl);
    }
    c.line(&pts, "axp-c", "", false);
    if pts.len() <= 60 {
        c.dots(&pts, "axp-d", 3.0);
    }
    c.line(&[(1.0, whole), (seq.len() as f64, whole)], "axp-c2", "whole sequence", true);
    c.mark(best, &format!("{} at {}", tick(best.1), tick(best.0)));
    let body = c.body();
    Some(s.chart("profile", title, &format!("{what} along the {}-residue sequence: each point is the model's value for one \
                                             {w}-letter window ({} windows); line: the whole sequence, {}. Marked: the highest \
                                             window, centred at {}, {}.", seq.len(), pts.len(), tick(whole), tick(best.0),
                                            tick(best.1)), &body))
}

pub fn gc_profile(s: &Solved) -> Option<Figure> {
    let n = seq_of(s, "seq", 100000)?.len();
    let w = (n / 3).clamp(3, 100).min(n);
    profile(s, "gc_fraction", w, "GC fraction in a sliding window along a sequence", "GC fraction", None, "GC fraction")
}

pub fn gravy_profile(s: &Solved) -> Option<Figure> {
    let n = seq_of(s, "seq", 100000)?.len();
    let w = if n < 30 { 3.min(n) } else { 9 };
    profile(s, "gravy", w, "Hydropathy plot: GRAVY along a protein", "GRAVY (hydropathy)", Some((0.0, "0")),
            "Mean hydropathy (GRAVY)")
}

pub fn instability_profile(s: &Solved) -> Option<Figure> {
    let n = seq_of(s, "seq", 100000)?.len();
    let w = if n >= 30 { 20 } else { (n / 2).max(3).min(n) };
    profile(s, "index", w, "Instability index along a protein", "instability index", Some((40.0, "40: unstable above")),
            "Instability index")
}

// ------------------------------------------------------------------------------------------------ k-mers
pub fn kmer_sequence(s: &Solved) -> Option<Figure> {
    let seq = seq_of(s, "seq", 60)?;
    let kk = s.val("k").map(|v| v.num() as usize)?;
    let mf = s.out_text("most_frequent")?.to_ascii_uppercase();
    let (tot, dis, mc) = (s.out("total")?, s.out("distinct")?, s.out("max_count")?);
    let m: Vec<char> = mf.chars().collect();
    let mut hit = vec![false; seq.len()];
    if kk == m.len() && kk >= 1 && seq.len() >= kk {
        for i in 0..=seq.len() - kk {
            if seq[i..i + kk] == m[..] {
                for h in hit.iter_mut().skip(i).take(kk) {
                    *h = true;
                }
            }
        }
    }
    let (body, y) = seq_rows(&seq, seq.len().min(30), 12.0, &|i, _| if hit[i] { "axp-b" } else { "" });
    let (nb, y2) = notes(&[("returned:".into(), format!("total = {}, distinct = {}, most_frequent = {mf}, max_count = {}",
                                                        g(tot, 4), g(dis, 4), g(mc, 4)))], y + 16.0);
    Some(s.figure("sequence", &format!("Counting {kk}-mers along a sequence"),
                  &format!("The {}-letter sequence cut into overlapping {kk}-mers: the model counts {} {kk}-mers, {} of them \
                            distinct, the most frequent being {mf} ({} times, shaded where it occurs).", seq.len(), g(tot, 4),
                           g(dis, 4), g(mc, 4)), &(body + &nb), y2.max(90.0)))
}

// ------------------------------------------------------------------------------------------------ ORF
pub fn longest_orf(s: &Solved) -> Option<Figure> {
    let seq = seq_of(s, "seq", 60)?;
    let (len, plen, cnt) = (s.out("length_nt")? as usize, s.out("protein_length")?, s.out("orf_count")?);
    let stops = ["TAA", "TAG", "TGA"];
    // the ATG whose first in-frame stop closes `len` bases later
    let mut start = None;
    for i in 0..seq.len().saturating_sub(2) {
        if seq[i..i + 3] != ['A', 'T', 'G'] {
            continue;
        }
        let mut j = i;
        while j + 3 <= seq.len() {
            let cod: String = seq[j..j + 3].iter().collect();
            if stops.contains(&cod.as_str()) {
                if j + 3 - i == len {
                    start = Some(i);
                }
                break;
            }
            j += 3;
        }
        if start.is_some() {
            break;
        }
    }
    let st = start?;
    let orf: String = seq[st..st + len].iter().collect();
    let prot = {
        let a = s.solve_other("bio_translation", &[("seq", ti(&orf)), ("stops", ti("truncate"))])?;
        match a.values.iter().find(|v| v.0 == "protein").map(|v| v.1.clone()) {
            Some(serde_json::Value::String(t)) => t,
            _ => return None,
        }
    };
    let per = seq.len();
    let (body, y) = seq_rows(&seq, per, 30.0, &|i, _| {
        if i < st || i >= st + len { "" } else if i < st + 3 { "axp-b2" } else if i >= st + len - 3 { "axp-b3" } else if ((i - st) / 3) % 2 == 0 {
            "axp-b"
        } else {
            "axp-b4"
        }
    });
    let cell = (440.0 / per as f64).min(30.0);
    let mut out = legend(&[("axp-b2", "start codon"), ("axp-b", "codons"), ("axp-b3", "stop codon")], 18.0);
    out.push_str(&body);
    for (k, aa) in prot.chars().enumerate() {
        let x = 60.0 + cell * (st as f64 + 3.0 * k as f64 + 1.5) - 1.0;
        out.push_str(&format!("<text class=\"axp-vt\" x=\"{x:.1}\" y=\"{:.1}\" text-anchor=\"middle\">{aa}</text>", y + 12.0));
    }
    let (nb, y2) = notes(&[("returned:".into(), format!("length_nt = {len}, protein_length = {}, orf_count = {}", g(plen, 4),
                                                        g(cnt, 4)))], y + 36.0);
    out.push_str(&nb);
    Some(s.figure("orf", "Longest open reading frame of a DNA sequence",
                  &format!("The sequence (5' to 3') with its longest open reading frame shaded: the start codon ATG, the \
                            codons and the stop codon; under it, the model's translation of that frame, {prot}. The model \
                            returns length_nt = {len} (stop included), protein_length = {} and orf_count = {}; the frame is \
                            located as the ATG whose first in-frame stop closes {len} bases later.", g(plen, 4), g(cnt, 4)),
                  &out, y2.max(110.0)))
}

// ------------------------------------------------------------------------------------------------ strands
fn comp(c: char) -> Option<char> {
    Some(match c { 'A' => 'T', 'T' => 'A', 'G' => 'C', 'C' => 'G', 'U' => 'A', 'N' => 'N', _ => return None })
}

pub fn duplex(s: &Solved) -> Option<Figure> {
    let seq = seq_of(s, "seq", 40)?;
    let rc = s.out_text("revcomp")?.to_ascii_uppercase();
    let bottom: Vec<char> = seq.iter().map(|c| comp(*c)).collect::<Option<Vec<char>>>()?;
    let rev: String = bottom.iter().rev().collect();
    if rev != rc {
        return None;
    }
    let cl = |c: char| if c == 'A' || c == 'T' || c == 'U' { "axp-b" } else if c == 'G' || c == 'C' { "axp-b2" } else { "" };
    let (top, y) = seq_rows(&seq, seq.len(), 30.0, &|_, c| cl(c));
    let (bot, y2) = seq_rows(&bottom, seq.len(), y + 10.0, &|_, c| cl(c));
    let mut out = String::new();
    out.push_str("<text class=\"axp-lt\" x=\"16\" y=\"22\">given strand 5' to 3'</text>");
    out.push_str(&top);
    out.push_str(&bot);
    out.push_str(&format!("<text class=\"axp-lt\" x=\"16\" y=\"{:.1}\">partner strand, written 3' to 5' under it</text>", y2 + 12.0));
    out.push_str(&legend(&[("axp-b", "A-T pair"), ("axp-b2", "G-C pair")], y2 + 34.0));
    let (nb, y3) = notes(&[("returned:".into(), format!("revcomp = {rc}"))], y2 + 56.0);
    out.push_str(&nb);
    Some(s.figure("duplex", "Reverse complement: the partner strand of a DNA sequence",
                  &format!("The given strand {} (top, 5' to 3') paired with its partner strand (bottom, 3' to 5'); A-T and G-C \
                            pairs in two colours. Read 5' to 3', the partner strand is the model's reverse complement, \
                            revcomp = {rc}.", seq.iter().collect::<String>()), &out, y3.max(120.0)))
}

pub fn transcription(s: &Solved) -> Option<Figure> {
    let seq = seq_of(s, "seq", 40)?;
    let mine = s.txt("strand")?;
    let mut out = String::new();
    let mut y = 16.0;
    let mut said = vec![];
    for strand in ["coding", "template"] {
        let a = s.with_input(&[("strand", ti(strand))])?;
        let m = match a.values.iter().find(|v| v.0 == "mrna").map(|v| v.1.clone()) {
            Some(serde_json::Value::String(t)) => t,
            _ => return None,
        };
        let mc: Vec<char> = m.chars().collect();
        out.push_str(&format!("<text class=\"axp-bt\" x=\"16\" y=\"{:.1}\">{} strand given{}</text>", y + 8.0, strand,
                              if strand.eq_ignore_ascii_case(&mine) { " (the input)" } else { "" }));
        let (b1, y1) = seq_rows(&seq, seq.len().max(mc.len()), y + 16.0, &|_, c| if c == 'T' { "axp-b4" } else { "" });
        let (b2, y2) = seq_rows(&mc, seq.len().max(mc.len()), y1, &|_, c| if c == 'U' { "axp-b" } else { "" });
        out.push_str(&b1);
        out.push_str(&b2);
        out.push_str(&format!("<text class=\"axp-lt\" x=\"52\" y=\"{:.1}\" text-anchor=\"end\">DNA</text><text class=\"axp-lt\" \
                               x=\"52\" y=\"{:.1}\" text-anchor=\"end\">mRNA</text>", y + 36.0, y1 + 20.0));
        said.push(format!("{strand}: {m}"));
        y = y2 + 8.0;
    }
    Some(s.figure("sequence", "Transcription from the coding or the template strand",
                  &format!("The DNA {} (5' to 3') given once as the coding strand and once as the template strand, with the \
                            mRNA the model returns for each (5' to 3'): {}. T and U shaded.", seq.iter().collect::<String>(),
                           said.join("; ")), &out, y + 4.0))
}

pub fn codons(s: &Solved) -> Option<Figure> {
    let seq = seq_of(s, "seq", 60)?;
    let prot = s.out_text("protein")?;
    let pc: Vec<char> = prot.chars().collect();
    let nc = seq.len() / 3;
    if nc == 0 || pc.len() > nc {
        return None;
    }
    let per = (nc * 3).min(60);
    let cell = (440.0 / per as f64).min(30.0);
    let (body, y) = seq_rows(&seq[..nc * 3], per, 30.0, &|i, _| if (i / 3) % 2 == 0 { "axp-b" } else { "axp-b4" });
    let mut out = body;
    for k in 0..nc {
        let x = 60.0 + cell * (3.0 * k as f64 + 1.5) - 1.0;
        out.push_str(&format!("<text class=\"axp-t\" x=\"{x:.1}\" y=\"22\" text-anchor=\"middle\" font-size=\"11\">{}</text>", k + 1));
        if let Some(aa) = pc.get(k) {
            out.push_str(&format!("<text class=\"axp-vt\" x=\"{x:.1}\" y=\"{:.1}\" text-anchor=\"middle\">{aa}</text>", y + 12.0));
        }
    }
    let (nb, y2) = notes(&[("returned:".into(), format!("protein = {prot}, length = {}", g(s.out("length")?, 4)))], y + 36.0);
    out.push_str(&nb);
    Some(s.figure("codons", "Translation: codons and the amino acids they encode",
                  &format!("The coding sequence read in codons (numbered, alternately shaded) with the amino acid the model \
                            assigns under each codon by the standard genetic code (a stop is shown as *). Returned: protein = \
                            {prot}."), &out, y2.max(110.0)))
}

pub fn hamming(s: &Solved) -> Option<Figure> {
    let a = seq_of(s, "seq1", 60)?;
    let b = seq_of(s, "seq2", 60)?;
    if a.len() != b.len() {
        return None;
    }
    let mm = s.out("mismatches")?;
    let diff: Vec<bool> = a.iter().zip(&b).map(|(x, y)| x != y).collect();
    if diff.iter().filter(|d| **d).count() as f64 != mm {
        return None;
    }
    let (t, y1) = seq_rows(&a, a.len(), 12.0, &|i, _| if diff[i] { "axp-b" } else { "" });
    let (bt, y2) = seq_rows(&b, a.len(), y1 + 12.0, &|i, _| if diff[i] { "axp-b" } else { "" });
    let cell = (440.0 / a.len() as f64).min(30.0);
    let mut out = t + &bt;
    for (i, d) in diff.iter().enumerate() {
        if !d {
            let x = 60.0 + cell * i as f64 + (cell - 2.0) / 2.0;
            out.push_str(&format!("<path class=\"axp-k\" d=\"M{x:.1},{:.1} L{x:.1},{:.1}\"/>", y1 - 2.0, y1 + 10.0));
        }
    }
    let (nb, y3) = notes(&[("returned:".into(), format!("mismatches = {}, p_distance = {}, identity = {}", g(mm, 4),
                                                        tick(s.out("p_distance")?), tick(s.out("identity")?)))], y2 + 16.0);
    out.push_str(&nb);
    Some(s.figure("alignment", "Hamming distance: mismatched positions of two sequences",
                  &format!("Two aligned sequences of {} letters: positions where they differ are shaded, equal positions \
                            joined by a bar. The model returns mismatches = {}.", a.len(), g(mm, 4)), &out, y3.max(110.0)))
}

// ------------------------------------------------------------------------------------------------ restriction
pub fn digest(s: &Solved) -> Option<Figure> {
    let seq = seq_of(s, "seq", 100000)?;
    let enz = s.txt("enzyme")?;
    let pos = s.out_list("positions").unwrap_or_default();
    let frags = s.out_list("fragments")?;
    if (frags.iter().sum::<f64>() - seq.len() as f64).abs() > 0.5 {
        return None;
    }
    let n = seq.len() as f64;
    let (x0, w, y) = (40.0, 440.0, 40.0);
    let mut out = format!("<text class=\"axp-mt\" x=\"16\" y=\"22\">{} cuts: {}</text>", esc(&enz), pos.len());
    let classes = ["axp-b", "axp-b2", "axp-b3", "axp-b4"];
    let mut at = 0.0;
    for (k, f) in frags.iter().enumerate() {
        let (a, b) = (x0 + w * at / n, x0 + w * (at + f) / n);
        out.push_str(&format!("<rect class=\"{}\" x=\"{a:.1}\" y=\"{y}\" width=\"{:.1}\" height=\"20\"/>\
                               <text class=\"axp-vt\" x=\"{:.1}\" y=\"{}\" text-anchor=\"middle\">{} bp</text>",
                              classes[k % 4], (b - a).max(1.0), (a + b) / 2.0, y + 38.0, g(*f, 6)));
        at += f;
    }
    for p in &pos {
        let x = x0 + w * (p - 1.0) / n;
        out.push_str(&format!("<path class=\"axp-k\" d=\"M{x:.1},{} L{x:.1},{}\"/><text class=\"axp-t\" x=\"{x:.1}\" y=\"{}\" \
                               text-anchor=\"middle\" font-size=\"11\">{}</text>", y - 6.0, y + 26.0, y - 9.0, g(*p, 6)));
    }
    out.push_str(&format!("<text class=\"axp-t\" x=\"{x0}\" y=\"{}\" font-size=\"11\">1</text><text class=\"axp-t\" x=\"{}\" \
                           y=\"{}\" text-anchor=\"end\" font-size=\"11\">{}</text>", y + 56.0, x0 + w, y + 56.0, seq.len()));
    Some(s.figure("digest", &format!("Restriction map: where {enz} cuts the sequence"),
                  &format!("{enz} digest of the {}-base sequence (5' to 3'): the model returns {} cut positions ({}, the first \
                            base after each cut) and the fragment lengths {} bp, each fragment shaded in its own colour.",
                           seq.len(), pos.len(), pos.iter().map(|p| g(*p, 6)).collect::<Vec<_>>().join(", "),
                           frags.iter().map(|p| g(*p, 6)).collect::<Vec<_>>().join(", ")), &out, 110.0))
}

pub fn gel(s: &Solved) -> Option<Figure> {
    let seq = seq_of(s, "seq", 100000)?;
    let mine = s.txt("enzyme")?;
    let mut lanes: Vec<(String, usize, Vec<f64>)> = vec![];
    for e in crate::figs_common::choices(s, "enzyme") {
        let Some(a) = s.with_input(&[("enzyme", ti(&e))]) else { continue };
        let fr = k::list_of(&a, "fragments").unwrap_or_else(|| vec![seq.len() as f64]);
        let ns = k::num_of(&a, "n_sites").unwrap_or(0.0) as usize;
        lanes.push((e, ns, fr));
    }
    let mi = lanes.iter().position(|l| l.0.eq_ignore_ascii_case(&mine))?;
    let me = lanes.remove(mi);
    lanes.sort_by(|a, b| b.1.cmp(&a.1));
    let mut shown = vec![me];
    shown.extend(lanes.into_iter().take(7));
    let all: Vec<f64> = shown.iter().flat_map(|l| l.2.iter().copied()).filter(|v| *v > 0.0).collect();
    let (lo, hi) = Chart::range_pos(&all);
    let ya = crate::figs_common::Axis::new(lo / 1.5, hi * 1.5, true, 290.0, 50.0, false);
    let mut out = String::from("<text class=\"axp-mt\" x=\"16\" y=\"22\">fragment lengths (bp), one lane per enzyme</text>");
    for t in ya.ticks() {
        let y = ya.px(t);
        out.push_str(&format!("<line class=\"axp-g\" x1=\"70\" y1=\"{y:.1}\" x2=\"504\" y2=\"{y:.1}\"/><text class=\"axp-t\" \
                               x=\"64\" y=\"{:.1}\" text-anchor=\"end\">{}</text>", y + 5.0, tick(t)));
    }
    let lw = 430.0 / shown.len() as f64;
    for (i, (e, ns, fr)) in shown.iter().enumerate() {
        let x = 74.0 + lw * i as f64;
        out.push_str(&format!("<rect class=\"axp-ka\" x=\"{x:.1}\" y=\"46\" width=\"{:.1}\" height=\"248\"/>", lw - 8.0));
        let mut seen: Vec<f64> = vec![];
        for f in fr {
            if seen.iter().any(|v| (v - f).abs() < 0.5) {
                continue;
            }
            seen.push(*f);
            let y = ya.px(*f);
            out.push_str(&format!("<rect class=\"{}\" x=\"{:.1}\" y=\"{:.1}\" width=\"{:.1}\" height=\"5\"/>",
                                  if i == 0 { "axp-b" } else { "axp-b4" }, x + 6.0, y - 2.5, lw - 20.0));
        }
        out.push_str(&format!("<text class=\"axp-bt\" x=\"{:.1}\" y=\"312\" text-anchor=\"middle\" font-size=\"11\">{}</text>\
                               <text class=\"axp-t\" x=\"{:.1}\" y=\"328\" text-anchor=\"middle\" font-size=\"11\">{ns} cut{}</text>",
                              x + (lw - 8.0) / 2.0, esc(&clip(e, 8)), x + (lw - 8.0) / 2.0, if *ns == 1 { "" } else { "s" }));
    }
    Some(s.figure("gel", "Restriction fragments of one sequence under several enzymes",
                  &format!("Fragments of the {}-base sequence cut by {} enzymes, its own ({mine}) first and then those of the \
                            model's enzymes that cut it most, one lane each, on a logarithmic length axis as on a gel: every \
                            band is a fragment length the model returns. Number of sites: {}.", seq.len(), shown.len(),
                           shown.iter().map(|l| format!("{} {}", l.0, l.1)).collect::<Vec<_>>().join(", ")), &out, 340.0))
}

// ------------------------------------------------------------------------------------------------ proteins
pub fn titration(s: &Solved) -> Option<Figure> {
    let seq: String = seq_of(s, "seq", 2000)?.iter().collect();
    let pi = s.out("pI")?;
    let c7 = s.out("charge_at_7")?;
    let peps = [seq.clone(), "DDEEK".to_string(), "KKRHA".to_string()];
    let phs = k::span(0.0, 14.0, 57, false);
    let mut curves: Vec<(usize, String, Vec<(f64, f64)>, Option<f64>)> = vec![];
    for (k_, p) in peps.iter().enumerate() {
        if k_ > 0 && *p == seq {
            continue;
        }
        let pts: Vec<(f64, f64)> = phs.iter().filter_map(|ph| other(s, "bio_protein_net_charge", &[("seq", ti(p)), ("pH", fi(*ph))],
                                                                     "charge").map(|q| (*ph, q))).collect();
        if pts.len() < 2 {
            continue;
        }
        let ppi = if k_ == 0 { Some(pi) } else { other(s, "bio_isoelectric_point", &[("seq", ti(p))], "pI") };
        curves.push((k_, p.clone(), pts, ppi));
    }
    if curves.is_empty() || curves[0].0 != 0 {
        return None;
    }
    let mut all: Vec<f64> = curves.iter().flat_map(|c| c.2.iter().map(|p| p.1)).collect();
    all.push(0.0);
    let (lo, hi) = Chart::range(&all);
    let mut c2 = Chart::new((0.0, 14.0), (lo, hi), (false, false), (false, true), "pH", "net charge",
                            &format!("{}: pI = {}", clip(&seq, 12), tick(pi)));
    c2.hguide(0.0, "");
    for (k_, p, pts, ppi) in &curves {
        c2.line(pts, crate::figs_common::LINES[*k_], &clip(p, 10), false);
        if let Some(v) = ppi {
            c2.dots(&[(*v, 0.0)], crate::figs_common::DOTS[*k_], 4.5);
        }
    }
    c2.mark((7.0, c7), &format!("{} at pH 7", tick(c7)));
    let body = c2.body();
    let said: Vec<String> = curves.iter().map(|c| format!("{} pI = {}", c.1, c.3.map(tick).unwrap_or_default())).collect();
    Some(s.chart("titration", "Isoelectric points: net charge against pH",
                 &format!("Net charge against pH for the given peptide and two references (an acidic DDEEK, a basic KKRHA), \
                           each curve 57 runs of the model's net-charge capability; the dot on the zero line is the isoelectric \
                           point the model returns: {}. Marked: the given peptide carries {} at pH 7 (charge_at_7).",
                          said.join(", "), tick(c7)), &body))
}

pub fn extinction_heat(s: &Solved) -> Option<Figure> {
    let er = s.out("eps_reduced")?;
    let mut m = vec![];
    for a in 0..4 {
        let mut row = vec![];
        for b in 0..5 {
            let p = format!("GGGGGGGG{}{}", "W".repeat(a), "Y".repeat(b));
            row.push(run(s, &[("seq", ti(&p))], "eps_reduced"));
        }
        m.push(row);
    }
    if m.iter().flatten().all(|v| v.is_none()) {
        return None;
    }
    let rows: Vec<String> = (0..4).map(|a| format!("{a} W")).collect();
    let cols: Vec<String> = (0..5).map(|b| format!("{b} Y")).collect();
    let (body, h) = crate::figs_bio::heat_table(&m, &rows, &cols, "tryptophans", "tyrosines",
                                                 &format!("given sequence: eps_reduced = {}", g(er, 6)), None, &[], 6);
    Some(s.figure("heatmap", "Extinction coefficient by tryptophan and tyrosine count",
                  &format!("Molar extinction coefficient at 280 nm (reduced cysteines) of the peptide G8 followed by a \
                            tryptophans and b tyrosines, a = 0 to 3 (rows), b = 0 to 4 (columns): each cell a run of the model. \
                            The given sequence gives eps_reduced = {}, eps_cystines = {}.", g(er, 6), g(s.out("eps_cystines")?, 6)),
                  &body, h))
}

pub fn mass_bars(s: &Solved) -> Option<Figure> {
    let mass = s.out("mass")?;
    let mut v: Vec<(String, f64)> = vec![];
    for aa in "ACDEFGHIKLMNPQRSTVWY".chars() {
        let p = aa.to_string().repeat(10);
        if let Some(m) = run(s, &[("seq", ti(&p))], "mass") {
            v.push((format!("{aa}10"), m));
        }
    }
    if v.len() < 2 {
        return None;
    }
    v.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal));
    let labels: Vec<String> = v.iter().map(|x| x.0.clone()).collect();
    let vals: Vec<f64> = v.iter().map(|x| x.1).collect();
    let n = vals.len();
    let body = bars_body(&labels, &vals, &[0, n - 1], "mass of the chain (Da)", &format!("given protein: {} Da", g(mass, 7)), false)?;
    Some(s.chart("bars", "Protein mass: ten-residue chains of each amino acid",
                 &format!("Mass of a ten-residue chain of each amino acid, each a run of the model, sorted. Marked: the \
                           lightest, {} ({} Da), and the heaviest, {} ({} Da). The given protein weighs {} Da.", labels[0],
                          g(vals[0], 6), labels[n - 1], g(vals[n - 1], 6), g(mass, 7)), &body))
}

pub fn composition_bars(s: &Solved) -> Option<Figure> {
    let m = crate::figs_common::out_map(s, "percent")?;
    let mc = s.out_text("most_common")?;
    let labels: Vec<String> = m.iter().map(|x| x.0.clone()).collect();
    let vals: Vec<f64> = m.iter().map(|x| x.1).collect();
    let hl: Vec<usize> = labels.iter().enumerate().filter(|(_, l)| **l == mc).map(|(i, _)| i).collect();
    let body = bars_body(&labels, &vals, &hl, "percent of residues", &format!("most common: {mc}"), false)?;
    Some(s.chart("bars", "Amino-acid composition of a protein",
                 &format!("Amino-acid composition of the {}-residue protein: one bar per amino acid, the percentage the model \
                           returns. Highlighted: the most common residue, {mc}.", g(s.out("length")?, 6)), &body))
}

pub fn structure_bars(s: &Solved) -> Option<Figure> {
    let seq: String = seq_of(s, "seq", 100000)?.iter().collect();
    let peps = [seq.clone(), "VIVYTFWE".into(), "AELKKLLEEA".into(), "NPGDGSPN".into()];
    let mut labels = vec![];
    let mut groups: Vec<(String, Vec<Option<f64>>)> = vec![("helix".into(), vec![]), ("turn".into(), vec![]),
                                                           ("sheet".into(), vec![])];
    for (i, p) in peps.iter().enumerate() {
        if i > 0 && *p == seq {
            continue;
        }
        let a = s.with_input(&[("seq", ti(p))])?;
        labels.push(if i == 0 { format!("{} (given)", clip(p, 10)) } else { p.clone() });
        for (gname, v) in groups.iter_mut() {
            v.push(k::num_of(&a, gname));
        }
    }
    let body = k::bar_groups(&labels, &groups, "", "fraction of residues", "secondary-structure propensities");
    Some(s.chart("bars", "Secondary-structure propensities of four peptides",
                 &format!("Fractions of helix-, turn- and sheet-forming residues the model returns for the given peptide \
                           (helix {}, turn {}, sheet {}) and three references: a beta-branched VIVYTFWE, a helix-like \
                           AELKKLLEEA and a turn-rich NPGDGSPN. The classes overlap, so the bars of a peptide need not sum to 1.",
                          tick(s.out("helix")?), tick(s.out("turn")?), tick(s.out("sheet")?)), &body))
}

// ------------------------------------------------------------------------------------------------ samples of sequences
fn sample_rows(s: &Solved) -> Option<(Vec<Vec<char>>, Vec<bool>)> {
    let t = s.txt("seqs")?;
    let seqs: Vec<Vec<char>> = t.split(|c: char| c.is_whitespace() || c == ',').filter(|x| !x.is_empty())
        .map(|x| x.chars().map(|c| c.to_ascii_uppercase()).collect()).collect();
    if seqs.len() < 2 || seqs.len() > 12 || seqs.iter().any(|x| x.len() != seqs[0].len()) || seqs[0].len() > 40 {
        return None;
    }
    let segr: Vec<bool> = (0..seqs[0].len()).map(|j| seqs.iter().any(|x| x[j] != seqs[0][j])).collect();
    Some((seqs, segr))
}

pub fn sample(s: &Solved, lines: &[(String, String)], title: &str, what: &str) -> Option<Figure> {
    let (seqs, segr) = sample_rows(s)?;
    let n = seqs[0].len();
    let cell = (440.0 / n as f64).min(30.0);
    let mut out = String::new();
    for (j, sg) in segr.iter().enumerate() {
        if *sg {
            out.push_str(&format!("<text class=\"axp-bt\" x=\"{:.1}\" y=\"22\" text-anchor=\"middle\">s</text>",
                                  60.0 + cell * j as f64 + (cell - 2.0) / 2.0));
        }
    }
    let mut y = 30.0;
    for (r, sq) in seqs.iter().enumerate() {
        out.push_str(&format!("<text class=\"axp-t\" x=\"52\" y=\"{:.1}\" text-anchor=\"end\" font-size=\"12\">{}</text>",
                              y + cell * 0.66, r + 1));
        let (b, y2) = seq_rows(sq, n, y, &|j, c| if c != seqs[0][j] { "axp-b" } else if segr[j] { "axp-b4" } else { "" });
        out.push_str(&b);
        y = y2;
    }
    let ns = segr.iter().filter(|x| **x).count();
    let (nb, y3) = notes(lines, y + 16.0);
    out.push_str(&nb);
    Some(s.figure("sample", title,
                  &format!("A sample of {} aligned sequences; the {ns} segregating sites (marked s, read from the sample) are \
                            shaded where a sequence differs from the first. {what}", seqs.len()), &out, y3.max(120.0)))
}

pub fn nd_sample(s: &Solved) -> Option<Figure> {
    let (kk, pi) = (s.out("k")?, s.out("pi")?);
    sample(s, &[("returned:".into(), format!("k = {}, pi = {}", tick(kk), tick(pi)))],
           "Nucleotide diversity of a sample of sequences",
           &format!("The model returns the mean number of differences between two sequences k = {} and the nucleotide \
                     diversity pi = {}.", tick(kk), tick(pi)))
}

pub fn watterson_sample(s: &Solved) -> Option<Figure> {
    let (sg, th, ts) = (s.out("S")?, s.out("theta")?, s.out("theta_site")?);
    sample(s, &[("returned:".into(), format!("S = {}, theta = {}, theta per site = {}", g(sg, 4), tick(th), tick(ts)))],
           "Watterson's theta from the segregating sites of a sample",
           &format!("The model counts S = {} segregating sites and returns Watterson's theta = {} ({} per site).", g(sg, 4),
                    tick(th), tick(ts)))
}

pub fn tajima_sample(s: &Solved) -> Option<Figure> {
    let (d, kk) = (s.out("D")?, s.out("k")?);
    let seqs = s.txt("seqs")?;
    let th = other(s, "bio_watterson_theta", &[("seqs", ti(&seqs))], "theta")?;
    sample(s, &[("returned:".into(), format!("D = {}, k = {}", tick(d), tick(kk))),
                ("Watterson's theta:".into(), tick(th))],
           "Tajima's D: pairwise differences against segregating sites",
           &format!("Tajima's D compares the mean pairwise difference k = {} (returned with D) with Watterson's theta = {} \
                     (the model's Watterson capability on the same sample); the model returns D = {}.", tick(kk), tick(th),
                    tick(d)))
}

pub fn tajima_bars(s: &Solved) -> Option<Figure> {
    let mine = s.txt("seqs")?;
    let samples = [("rare variants", "AAAAAAAA TAAAAAAA ATAAAAAA AATAAAAA AAATAAAA AAAATAAA"), ("given", mine.as_str()),
                   ("mixed", "AAAAAAAA AAAAAAAT AAAAAATT TAAAAAAA AAAAAATT TTAAAAAA"),
                   ("intermediate", "AAAAAAAA AAAAAAAA AAAAAAAA TTTTAAAA TTTTAAAA TTTTAAAA")];
    let mut labels = vec![];
    let mut vals = vec![];
    for (n, sq) in samples {
        if let Some(d) = run(s, &[("seqs", ti(sq))], "D") {
            labels.push(n.to_string());
            vals.push(d);
        }
    }
    let hl = vec![labels.iter().position(|l| l == "given")?];
    let classes: Vec<&str> = (0..labels.len()).map(|i| if hl.contains(&i) { "axp-b" } else { "axp-b4" }).collect();
    let body = bars_full(&labels, &vals, &classes, None, "Tajima's D", &format!("given sample: D = {}", tick(s.out("D")?)),
                         false, None)?;
    Some(s.chart("bars", "Tajima's D for rare and intermediate-frequency variants",
                 &format!("Tajima's D the model returns for four samples: variants each carried by one sequence (an excess of \
                           rare variants gives D below 0), the given sample, a mixed sample, and two groups differing at four \
                           sites (variants at intermediate frequency give D above 0): {}.",
                          labels.iter().zip(&vals).map(|(l, v)| format!("{l} {}", tick(*v))).collect::<Vec<_>>().join(", ")),
                 &body))
}
