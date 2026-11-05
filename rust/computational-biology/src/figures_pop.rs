//! Genetics, populations and diversity: crosses (Punnett squares and genotype probabilities), allele-frequency
//! trajectories, heat tables over two inputs, Hardy-Weinberg diagrams, Leslie projections and sample bar charts.
//! Every number is a run of the model.

use crate::figs_bio::{cat_chart, chart_bars, chart_xlabels, g, heat_table, ii, other, run, ti, vi, fi};
use crate::figs_common::{self as c, esc, tick, Chart, FieldSpec, DOTS, LINES};
use ax_engine::figkit::{self as k, Solved};
use ax_engine::figures::Figure;

fn list(v: &[f64]) -> String {
    format!("[{}]", v.iter().map(|x| g(*x, 6)).collect::<Vec<_>>().join(", "))
}

// ------------------------------------------------------------------------------------------------ diversity bars
pub fn ace_bars(s: &Solved) -> Option<Figure> {
    let counts = s.vec("counts")?;
    let thr = s.val("rare_threshold").map(|v| v.num())?;
    let ace = s.out("ace")?;
    if counts.is_empty() || counts.len() > 60 || counts.iter().any(|x| !(*x > 0.0)) {
        return None;
    }
    let mut v = counts.clone();
    v.sort_by(|a, b| b.partial_cmp(a).unwrap_or(std::cmp::Ordering::Equal));
    let (lo, hi) = Chart::range_pos(&v);
    let mut ch = cat_chart(v.len(), (lo / 1.5, hi * 1.5), true, "individuals counted (logarithmic)",
                           &format!("ACE = {} species", tick(ace)));
    for (i, x) in v.iter().enumerate() {
        let class = if *x <= thr { "axp-b" } else { "axp-b4" };
        chart_bars(&mut ch, &[(i as f64 + 0.5, lo / 1.5, *x)], (416.0 / v.len() as f64 * 0.7).min(50.0), class);
    }
    ch.hguide(thr, &format!("rare at or below {}", g(thr, 6)));
    let labels: Vec<(f64, String)> = (0..v.len()).filter(|i| v.len() <= 20 || i % 5 == 0)
        .map(|i| (i as f64 + 0.5, format!("{}", i + 1))).collect();
    chart_xlabels(&mut ch, &labels);
    ch.text(292.0, 326.0, "middle", "axp-l", "species, most abundant first");
    let body = ch.body();
    Some(s.chart("bars", "ACE richness from the rare species of a sample",
                 &format!("The sample counts = {} (logarithmic), with the rare-species threshold {}: ACE estimates the richness \
                           from the species at or below it (shaded). {} species are observed; the model returns ace = {}.",
                          list(&counts), g(thr, 6), counts.len(), tick(ace)), &body))
}

pub fn chao1_bars(s: &Solved) -> Option<Figure> {
    let counts = s.vec("counts")?;
    let (ch1, f1, f2) = (s.out("chao1")?, s.out("singletons")?, s.out("doubletons")?);
    let mut v: Vec<f64> = counts.iter().copied().filter(|x| *x > 0.0).collect();
    if v.is_empty() || v.len() > 60 {
        return None;
    }
    v.sort_by(|a, b| b.partial_cmp(a).unwrap_or(std::cmp::Ordering::Equal));
    let labels: Vec<String> = (1..=v.len()).map(|i| i.to_string()).collect();
    let hl: Vec<usize> = v.iter().enumerate().filter(|(_, x)| **x <= 2.0).map(|(i, _)| i).collect();
    let body = c::bars_body(&labels, &v, &hl, "individuals counted", &format!("Chao1 = {}", tick(ch1)), false)?;
    Some(s.chart("bars", "Chao1: unseen species estimated from the rarest counts",
                 &format!("The sample counts = {} as a rank-abundance bar chart: {} species observed, {} of them once and {} \
                           twice (highlighted). From these rare counts the model's Chao1 ({}) estimates the species present, \
                           seen or not: chao1 = {}.", list(&counts), v.len(), g(f1, 6), g(f2, 6), s.txt("form").unwrap_or_default(),
                          tick(ch1)), &body))
}

pub fn beta_bars(s: &Solved) -> Option<Figure> {
    let (x, y) = (s.vec("x")?, s.vec("y")?);
    if x.len() != y.len() || x.len() > 24 {
        return None;
    }
    let (bc, ja, so) = (s.out("bray_curtis")?, s.out("jaccard")?, s.out("sorensen")?);
    let labels: Vec<String> = (1..=x.len()).map(|i| format!("sp {i}")).collect();
    let groups = vec![("sample x".to_string(), x.iter().map(|v| Some(*v)).collect()),
                      ("sample y".to_string(), y.iter().map(|v| Some(*v)).collect())];
    let mut body = k::bar_groups(&labels, &groups, "species", "individuals", "two samples of the same species");
    let (nb, _) = c::notes(&[("returned:".into(), format!("Bray-Curtis {}, Jaccard {}, Sorensen {}", tick(bc),
                                                                          tick(ja), tick(so)))], 360.0);
    body.push_str(&nb);
    Some(s.figure("bars", "Beta diversity of two samples",
                  &format!("Two community samples, x = {} and y = {} (individuals of the same species), and the dissimilarities \
                            the model returns: Bray-Curtis {} (from the abundances), Jaccard {} and Sorensen {} (from presence \
                            and absence only).", list(&x), list(&y), tick(bc), tick(ja), tick(so)), &body, 374.0))
}

pub fn ne_fluct(s: &Solved) -> Option<Figure> {
    let sizes = s.vec("sizes")?;
    let (ne, ret) = (s.out("Ne")?, s.out("retained")?);
    if sizes.len() > 40 || sizes.iter().any(|x| !(*x > 0.0)) {
        return None;
    }
    let mean = sizes.iter().sum::<f64>() / sizes.len() as f64;
    let (lo, hi) = Chart::range_pos(&sizes.iter().copied().chain([ne]).collect::<Vec<_>>());
    let mut ch = cat_chart(sizes.len(), (lo / 2.0, hi * 2.0), true, "census size (logarithmic)", &format!("Ne = {}", tick(ne)));
    for (i, v) in sizes.iter().enumerate() {
        chart_bars(&mut ch, &[(i as f64 + 0.5, lo / 2.0, *v)], (416.0 / sizes.len() as f64 * 0.6).min(60.0), "axp-b4");
    }
    ch.line(&[(0.0, mean), (sizes.len() as f64, mean)], "axp-c4", "arithmetic mean", true);
    ch.line(&[(0.0, ne), (sizes.len() as f64, ne)], "axp-c", "Ne (harmonic mean)", false);
    let labels: Vec<(f64, String)> = (0..sizes.len()).filter(|i| sizes.len() <= 20 || i % 5 == 0)
        .map(|i| (i as f64 + 0.5, format!("{}", i + 1))).collect();
    chart_xlabels(&mut ch, &labels);
    ch.text(292.0, 326.0, "middle", "axp-l", "generation");
    let body = ch.body();
    Some(s.chart("bars", "Effective size of a fluctuating population",
                 &format!("Census sizes {} in successive generations (logarithmic axis): the model returns the effective size \
                           Ne = {}, the harmonic mean, which a bottleneck drags below the arithmetic mean {} (computed from the \
                           sizes for comparison); heterozygosity retained {}.", list(&sizes), tick(ne), tick(mean), tick(ret)),
                 &body))
}

pub fn ne_sex_field(s: &Solved) -> Option<Figure> {
    c::field_map(s, &FieldSpec { z: "Ne", x: "N_m", y: "N_f", log: (false, false), zunit: "",
                                 title: "Effective population size over the numbers of breeding males and females" })
}

// ------------------------------------------------------------------------------------------------ crosses
/// The genes of a genotype written as letter pairs ("AaBb" -> [('A','a'), ('B','b')]).
fn genes(t: &str) -> Option<Vec<(char, char)>> {
    let cs: Vec<char> = t.chars().filter(|c| c.is_alphabetic()).collect();
    if cs.is_empty() || cs.len() % 2 != 0 || cs.len() > 6 {
        return None;
    }
    let v: Vec<(char, char)> = cs.chunks(2).map(|p| (p[0], p[1])).collect();
    v.iter().all(|(a, b)| a.eq_ignore_ascii_case(b)).then_some(v)
}

/// A genotype written in order: the dominant allele first in each gene.
fn canon(alleles: &[(char, char)]) -> String {
    alleles.iter().map(|(a, b)| if b.is_uppercase() && a.is_lowercase() { format!("{b}{a}") } else { format!("{a}{b}") })
        .collect()
}

/// The gametes of a parent: one allele of each gene.
fn gametes(gs: &[(char, char)]) -> Vec<Vec<char>> {
    let mut out: Vec<Vec<char>> = vec![vec![]];
    for (a, b) in gs {
        let mut next = vec![];
        for g_ in &out {
            for al in [a, b] {
                let mut x = g_.clone();
                x.push(*al);
                next.push(x);
            }
        }
        out = next;
    }
    out
}

fn square(p1: &[(char, char)], p2: &[(char, char)]) -> Option<(Vec<Vec<char>>, Vec<Vec<char>>, Vec<Vec<String>>)> {
    if p1.len() != p2.len() || p1.iter().zip(p2).any(|(a, b)| !a.0.eq_ignore_ascii_case(&b.0)) {
        return None;
    }
    let (g1, g2) = (gametes(p1), gametes(p2));
    let cells: Vec<Vec<String>> = g1.iter().map(|a| g2.iter().map(|b| {
        canon(&a.iter().zip(b).map(|(x, y)| (*x, *y)).collect::<Vec<_>>())
    }).collect()).collect();
    Some((g1, g2, cells))
}

fn draw_square(g1: &[Vec<char>], g2: &[Vec<char>], cells: &[Vec<String>], class: &dyn Fn(&str) -> &'static str, y0: f64)
               -> (String, f64) {
    let n = g1.len().max(g2.len());
    let cell = (380.0 / (n as f64 + 1.0)).min(64.0);
    let x0 = (520.0 - cell * (g2.len() as f64 + 1.0)) / 2.0;
    let mut out = String::new();
    for (j, gm) in g2.iter().enumerate() {
        out.push_str(&format!("<text class=\"axp-bt\" x=\"{:.1}\" y=\"{:.1}\" text-anchor=\"middle\">{}</text>",
                              x0 + cell * (j as f64 + 1.5), y0 + cell * 0.6, gm.iter().collect::<String>()));
    }
    for (i, gm) in g1.iter().enumerate() {
        let y = y0 + cell * (i as f64 + 1.0);
        out.push_str(&format!("<text class=\"axp-bt\" x=\"{:.1}\" y=\"{:.1}\" text-anchor=\"middle\">{}</text>",
                              x0 + cell * 0.5, y + cell * 0.6, gm.iter().collect::<String>()));
        for (j, gt) in cells[i].iter().enumerate() {
            let x = x0 + cell * (j as f64 + 1.0);
            let cl = class(gt);
            out.push_str(&format!("<rect class=\"{}\" x=\"{x:.1}\" y=\"{y:.1}\" width=\"{:.1}\" height=\"{:.1}\"/>\
                                   <rect class=\"axp-f\" x=\"{x:.1}\" y=\"{y:.1}\" width=\"{:.1}\" height=\"{:.1}\"/>\
                                   <text class=\"axp-lt\" x=\"{:.1}\" y=\"{:.1}\" text-anchor=\"middle\" font-size=\"{:.0}\">{}</text>",
                                  if cl.is_empty() { "axp-ka" } else { cl }, cell - 2.0, cell - 2.0, cell - 2.0, cell - 2.0,
                                  x + cell / 2.0 - 1.0, y + cell * 0.58, (cell * 0.28).clamp(9.0, 14.0), esc(gt)));
        }
    }
    (out, y0 + cell * (g1.len() as f64 + 1.0) + 10.0)
}

pub fn genotype_punnett(s: &Solved) -> Option<Figure> {
    let (p1, p2, off) = (genes(&s.txt("parent1")?)?, genes(&s.txt("parent2")?)?, genes(&s.txt("offspring")?)?);
    let prob = s.out("probability")?;
    let (g1, g2, cells) = square(&p1, &p2)?;
    let want = canon(&off);
    let hits = cells.iter().flatten().filter(|c| **c == want).count() as f64;
    let total = (g1.len() * g2.len()) as f64;
    if (hits / total - prob).abs() > 1e-9 {
        return None;
    }
    let (body, y) = draw_square(&g1, &g2, &cells, &|gt| if gt == want { "axp-b" } else { "" }, 34.0);
    let mut out = format!("<text class=\"axp-mt\" x=\"16\" y=\"22\">P({want}) = {}</text>", tick(prob));
    out.push_str(&body);
    let frac = s.out_text("fraction").unwrap_or_default();
    Some(s.figure("punnett", &format!("Punnett square: probability of {want}"),
                  &format!("Punnett square of the cross {} x {}: the gametes of each parent along the edges (one allele of each \
                            gene, equally likely), the offspring genotype in each cell; cells of the asked genotype {want} are \
                            shaded ({} of {}). The model returns probability = {} ({frac}).", canon(&p1), canon(&p2), hits,
                           total, tick(prob)), &out, y.max(150.0)))
}

pub fn genotype_bars(s: &Solved) -> Option<Figure> {
    let (p1, p2, off) = (genes(&s.txt("parent1")?)?, genes(&s.txt("parent2")?)?, genes(&s.txt("offspring")?)?);
    let (_, _, cells) = square(&p1, &p2)?;
    let mut kinds: Vec<String> = vec![];
    for c_ in cells.iter().flatten() {
        if !kinds.contains(c_) {
            kinds.push(c_.clone());
        }
    }
    kinds.sort();
    if kinds.len() > 27 {
        return None;
    }
    let want = canon(&off);
    let mut labels = vec![];
    let mut vals = vec![];
    for kd in &kinds {
        if let Some(p) = run(s, &[("offspring", ti(kd))], "probability") {
            labels.push(kd.clone());
            vals.push(p);
        }
    }
    let hl: Vec<usize> = labels.iter().enumerate().filter(|(_, l)| **l == want).map(|(i, _)| i).collect();
    let body = c::bars_body(&labels, &vals, &hl, "probability", &format!("sum = {}", tick(vals.iter().sum())), false)?;
    Some(s.chart("bars", &format!("Offspring genotype probabilities of {} x {}", canon(&p1), canon(&p2)),
                 &format!("Probability of every offspring genotype of the cross {} x {}, each bar a run of the model: {}; \
                           together they sum to {}. Highlighted: the asked genotype, {want}.", canon(&p1), canon(&p2),
                          labels.iter().zip(&vals).map(|(l, v)| format!("{l} {}", tick(*v))).collect::<Vec<_>>().join(", "),
                          tick(vals.iter().sum())), &body))
}

pub fn phenotype_punnett(s: &Solved) -> Option<Figure> {
    let (p1, p2) = (genes(&s.txt("parent1")?)?, genes(&s.txt("parent2")?)?);
    let ratio = s.out_text("ratio")?;
    let (g1, g2, cells) = square(&p1, &p2)?;
    if g1.len() > 8 {
        return None;
    }
    // a cell's phenotype: per gene, dominant (an upper-case allele) or recessive
    let pheno = |gt: &str| -> String {
        let cs: Vec<char> = gt.chars().collect();
        cs.chunks(2).map(|p| if p.iter().any(|c| c.is_uppercase()) { format!("{}_", p[0].to_ascii_uppercase()) } else {
            format!("{}{}", p[0].to_ascii_lowercase(), p[0].to_ascii_lowercase())
        }).collect::<Vec<_>>().join(" ")
    };
    let mut classes: Vec<(String, usize)> = vec![];
    for c_ in cells.iter().flatten() {
        let p = pheno(c_);
        match classes.iter_mut().find(|x| x.0 == p) {
            Some(x) => x.1 += 1,
            None => classes.push((p, 1)),
        }
    }
    classes.sort_by(|a, b| b.1.cmp(&a.1).then(a.0.cmp(&b.0)));
    let mut counts: Vec<usize> = classes.iter().map(|c_| c_.1).collect();
    let gcd = counts.iter().fold(0usize, |a, b| { let (mut x, mut y) = (a, *b); while y != 0 { let t = x % y; x = y; y = t; } x });
    if gcd > 0 {
        counts.iter_mut().for_each(|c_| *c_ /= gcd);
    }
    let mine = counts.iter().map(|c_| c_.to_string()).collect::<Vec<_>>().join(":");
    if mine != ratio.trim() {
        return None;
    }
    let fills = ["axp-b", "axp-b2", "axp-b3", "axp-b4", "axp-s", "axp-s2", "axp-s3", "axp-s4"];
    let cls = |gt: &str| -> &'static str {
        let p = pheno(gt);
        let i = classes.iter().position(|c_| c_.0 == p).unwrap_or(0);
        fills[i % 8]
    };
    let mut out = format!("<text class=\"axp-mt\" x=\"16\" y=\"22\">phenotype ratio {}</text>", esc(&ratio));
    let (body, y) = draw_square(&g1, &g2, &cells, &cls, 34.0);
    out.push_str(&body);
    let mut x = 16.0;
    for (i, (p, n)) in classes.iter().enumerate() {
        out.push_str(&format!("<rect class=\"{}\" x=\"{x:.1}\" y=\"{:.1}\" width=\"12\" height=\"12\"/><text class=\"axp-lt\" \
                               x=\"{:.1}\" y=\"{:.1}\">{} ({n})</text>", fills[i % 8], y + 2.0, x + 16.0, y + 12.0, esc(p)));
        x += (p.chars().count() + 5) as f64 * 7.4 + 24.0;
    }
    let lines = vec![("returned:".to_string(), format!("ratio {ratio}, {} genotypes, {} phenotypes", g(s.out("n_genotypes")?, 4),
                                                       g(s.out("n_phenotypes")?, 4)))];
    let (nb, y2) = c::notes(&lines, y + 40.0);
    out.push_str(&nb);
    Some(s.figure("punnett", &format!("Punnett square of {} x {}: phenotype classes", canon(&p1), canon(&p2)),
                  &format!("Punnett square of {} x {}, each offspring cell shaded by its phenotype class (A_ means at least one \
                            dominant A, aa recessive; classes and their cell counts read from the square). The model returns \
                            ratio = {ratio}, p_all_dominant = {}, p_all_recessive = {}.", canon(&p1), canon(&p2),
                           tick(s.out("p_all_dominant")?), tick(s.out("p_all_recessive")?)), &out, y2.max(150.0)))
}

// ------------------------------------------------------------------------------------------------ heat tables
fn grid_heat(s: &Solved, kind: &str, title: &str, rows: &[f64], cols: &[f64], f: &dyn Fn(f64, f64) -> Option<f64>,
             rlab: &dyn Fn(f64) -> String, clab: &dyn Fn(f64) -> String, rt: &str, ct: &str, header: &str,
             mine: Option<(f64, f64)>, caption: &str) -> Option<Figure> {
    let m: Vec<Vec<Option<f64>>> = rows.iter().map(|r| cols.iter().map(|c_| f(*r, *c_)).collect()).collect();
    if m.iter().flatten().filter(|v| v.is_some()).count() < 2 {
        return None;
    }
    let mark = mine.and_then(|(r, c_)| {
        let i = rows.iter().position(|x| (x - r).abs() < 1e-9)?;
        let j = cols.iter().position(|x| (x - c_).abs() < 1e-9)?;
        Some((i, j))
    });
    let rl: Vec<String> = rows.iter().map(|r| rlab(*r)).collect();
    let cl: Vec<String> = cols.iter().map(|x| clab(*x)).collect();
    let (body, h) = heat_table(&m, &rl, &cl, rt, ct, header, mark, &[], 3);
    Some(s.figure(kind, title, caption, &body, h))
}

fn with_mine(base: &[f64], v: f64) -> Vec<f64> {
    let mut out = base.to_vec();
    if !out.iter().any(|x| (x - v).abs() < 1e-9) {
        out.push(v);
    }
    out.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    out
}

pub fn drift_heat(s: &Solved) -> Option<Figure> {
    let (n, t) = (s.inp("N").or_else(|| s.val("N").map(|v| v.num()))?, s.val("t").map(|v| v.num())?);
    let ret = s.out("retained")?;
    let ns = with_mine(&[5.0, 10.0, 25.0, 50.0, 100.0, 250.0], n);
    let ts = with_mine(&[10.0, 25.0, 50.0, 100.0, 200.0], t);
    grid_heat(s, "heatmap", "Heterozygosity retained under drift", &ns, &ts,
              &|nn, tt| run(s, &[("N", ii(nn as i64)), ("i", ii(nn as i64)), ("t", ii(tt as i64))], "retained"),
              &|v| format!("N {}", g(v, 6)), &|v| format!("t {}", g(v, 6)), "population size N", "generations t",
              &format!("N = {}, t = {}: retained {}", g(n, 6), g(t, 6), tick(ret)), Some((n, t)),
              &format!("Fraction of heterozygosity retained under random genetic drift for population sizes {} (rows) after {} \
                        generations (columns), each cell a run of the model rounded to 3 digits. Ringed: the inputs given, N = \
                        {}, t = {}, retained = {}.", ns.iter().map(|v| g(*v, 6)).collect::<Vec<_>>().join(", "),
                       ts.iter().map(|v| g(*v, 6)).collect::<Vec<_>>().join(", "), g(n, 6), g(t, 6), tick(ret)))
}

pub fn fst_heat(s: &Solved) -> Option<Figure> {
    let p = s.vec("p")?;
    let fst = s.out("FST")?;
    let mine = (p.len() == 2).then(|| (p[0], p[1]));
    let base = [0.1, 0.3, 0.5, 0.7, 0.9];
    let (rows, cols) = match mine {
        Some((a, b)) if (0.0..=1.0).contains(&a) && (0.0..=1.0).contains(&b) => (with_mine(&base, a), with_mine(&base, b)),
        _ => (base.to_vec(), base.to_vec()),
    };
    grid_heat(s, "heatmap", "Fixation index FST between two subpopulations", &rows, &cols,
              &|a, b| run(s, &[("p", vi(&[a, b]))], "FST"), &|v| tick(v), &|v| tick(v), "allele frequency, subpopulation 1",
              "allele frequency, subpopulation 2", &format!("given p = {}: FST = {}", list(&p), tick(fst)), mine,
              &format!("Fixation index between two subpopulations whose allele frequencies are 0.1 to 0.9 (rows: the first, \
                        columns: the second, the given frequencies added), each cell a run of the model rounded to 3 digits: \
                        zero on the diagonal, largest for opposite frequencies. The given p = {} gives HS = {}, HT = {}, FST = {}.", list(&p),
                       tick(s.out("HS")?), tick(s.out("HT")?), tick(fst)))
}

pub fn het_adv_heat(s: &Solved) -> Option<Figure> {
    let (wa, wh, wb) = (s.inp("w_AA")?, s.inp("w_Aa")?, s.inp("w_aa")?);
    let pe = s.out("p_eq")?;
    let top = wh;
    let grid: Vec<f64> = (0..5).map(|i| top * (0.5 + 0.1 * i as f64)).chain([top * 0.95]).collect();
    let (rows, cols) = (with_mine(&grid, wa), with_mine(&grid, wb));
    grid_heat(s, "heatmap", "Equilibrium allele frequency under heterozygote advantage", &rows, &cols,
              &|a, b| run(s, &[("w_AA", fi(a)), ("w_aa", fi(b))], "p_eq"), &|v| tick(v), &|v| tick(v), "fitness w_AA",
              "fitness w_aa", &format!("w_AA = {}, w_aa = {}: p_eq = {}", tick(wa), tick(wb), tick(pe)), Some((wa, wb)),
              &format!("Equilibrium frequency of allele A under heterozygote advantage (w_Aa = {}) for homozygote fitnesses \
                        w_AA (rows) and w_aa (columns) from half to 0.95 of w_Aa and the given values, each cell a run of the model rounded to 3 \
                        digits: A is common where AA is the fitter homozygote. The inputs given, w_AA = {}, w_aa = {}, give \
                        p_eq = {}.", tick(wh), tick(wa), tick(wb), tick(pe)))
}

// ------------------------------------------------------------------------------------------------ trajectories
pub fn haploid_traj(s: &Solved) -> Option<Figure> {
    let (p0, sel, t) = (s.inp("p0")?, s.inp("s")?, s.val("t").map(|v| v.num())?);
    let pt = s.out("p_t")?;
    let mut ss = with_mine(&[0.25, 0.5, 1.0, 2.0], sel);
    ss.truncate(5);
    let tmax = 40.0f64.max(t);
    let ts: Vec<f64> = (1..=40).map(|i| (tmax * i as f64 / 40.0).round()).collect();
    let mut ch = Chart::new((0.0, tmax), (0.0, 1.0), (false, false), (false, false), "generations t", "frequency of A",
                            &format!("p0 = {}, s = {}", tick(p0), tick(sel)));
    for (k_, sv) in ss.iter().enumerate() {
        let mut pts = vec![(0.0, p0)];
        for tt in &ts {
            if let Some(v) = run(s, &[("s", fi(*sv)), ("t", ii(*tt as i64))], "p_t") {
                pts.push((*tt, v));
            }
        }
        let mine = (sv - sel).abs() < 1e-12;
        ch.line(&pts, if mine { "axp-c" } else { LINES[1 + k_ % 3] }, &format!("s = {}{}", tick(*sv), if mine { " (given)" } else { "" }), false);
    }
    ch.mark((t, pt), &format!("t = {}, p = {}", g(t, 6), tick(pt)));
    let body = ch.body();
    Some(s.chart("trajectory", "Haploid selection: allele frequency over generations",
                 &format!("Spread of a favoured allele under haploid selection from p0 = {} for selection coefficients {} (the \
                           given s = {} among them): each curve is p_t the model returns at 40 times up to {} generations (the \
                           start p0 drawn at t = 0). Marked: the inputs given, t = {}, p = {}.", tick(p0),
                          ss.iter().map(|v| tick(*v)).collect::<Vec<_>>().join(", "), tick(sel), g(tmax, 6), g(t, 6), tick(pt)),
                 &body))
}

pub fn lethal_traj(s: &Solved) -> Option<Figure> {
    let (q0, t) = (s.inp("q0")?, s.val("t").map(|v| v.num())?);
    let qt = s.out("q_t")?;
    let mut qs = vec![q0];
    for q in [0.2, 0.05] {
        if (q - q0).abs() > 1e-9 {
            qs.push(q);
        }
    }
    let tmax = 60.0f64.max(t);
    let ts: Vec<f64> = (1..=60).map(|i| (tmax * i as f64 / 60.0).round()).collect();
    let mut ch = Chart::new((0.0, tmax), (0.0, qs.iter().cloned().fold(0.0, f64::max)), (false, false), (false, true),
                            "generations t", "frequency of the lethal allele q", &format!("q0 = {}", tick(q0)));
    for (k_, q) in qs.iter().enumerate() {
        let mut pts = vec![(0.0, *q)];
        for tt in &ts {
            if let Some(v) = run(s, &[("q0", fi(*q)), ("t", ii(*tt as i64))], "q_t") {
                pts.push((*tt, v));
            }
        }
        ch.line(&pts, LINES[k_], &format!("q0 = {}{}", tick(*q), if k_ == 0 { " (given)" } else { "" }), false);
    }
    ch.mark((t, qt), &format!("t = {}, q = {}", g(t, 6), tick(qt)));
    let body = ch.body();
    Some(s.chart("trajectory", "Recessive lethal allele: frequency over generations",
                 &format!("Decline of a recessive lethal allele, whose homozygotes die before breeding, from q0 = {}: each curve \
                           is q_t the model returns for t = 1 to {} (q0 drawn at t = 0). Rare alleles hide in heterozygotes and \
                           fall ever more slowly. Marked: the inputs given, q0 = {}, t = {}, q_t = {}.",
                          qs.iter().map(|v| tick(*v)).collect::<Vec<_>>().join(", "), g(tmax, 6), tick(q0), g(t, 6), tick(qt)),
                 &body))
}

// ------------------------------------------------------------------------------------------------ Hardy-Weinberg
pub fn hw_bars(s: &Solved) -> Option<Figure> {
    let obs = [s.inp("n_AA").or_else(|| s.val("n_AA").map(|v| v.num()))?, s.val("n_Aa").map(|v| v.num())?,
               s.val("n_aa").map(|v| v.num())?];
    let exp = [s.out("exp_AA")?, s.out("exp_Aa")?, s.out("exp_aa")?];
    let labels: Vec<String> = ["AA", "Aa", "aa"].iter().map(|x| x.to_string()).collect();
    let groups = vec![("observed".to_string(), obs.iter().map(|v| Some(*v)).collect()),
                      ("expected".to_string(), exp.iter().map(|v| Some(*v)).collect())];
    let body = k::bar_groups(&labels, &groups, "genotype", "count",
                             &format!("chi2 = {}, p-value {}", tick(s.out("chi2")?), tick(s.out("p_value")?)));
    Some(s.chart("bars", "Hardy-Weinberg: observed against expected genotypes",
                 &format!("Observed genotype counts AA = {}, Aa = {}, aa = {} beside the counts expected under Hardy-Weinberg \
                           equilibrium, computed by the model from the allele frequencies p = {}, q = {}: {}, {}, {}. chi2 = {}, \
                           p-value {}; observed heterozygosity {} against expected {} (F = {}).", g(obs[0], 6), g(obs[1], 6),
                          g(obs[2], 6), tick(s.out("p")?), tick(s.out("q")?), tick(exp[0]), tick(exp[1]), tick(exp[2]),
                          tick(s.out("chi2")?), tick(s.out("p_value")?), tick(s.out("Ho")?), tick(s.out("He")?),
                          tick(s.out("F")?)), &body))
}

pub fn hw_triangle(s: &Solved) -> Option<Figure> {
    let obs = [s.val("n_AA").map(|v| v.num())?, s.val("n_Aa").map(|v| v.num())?, s.val("n_aa").map(|v| v.num())?];
    let exp = [s.out("exp_AA")?, s.out("exp_Aa")?, s.out("exp_aa")?];
    let n: f64 = obs.iter().sum();
    if !(n > 0.0) {
        return None;
    }
    // corners: AA bottom left, aa bottom right, Aa top
    let (a, b, cc) = ((90.0, 360.0), (430.0, 360.0), (260.0, 360.0 - 340.0 * 0.866));
    let pt = |f: [f64; 3]| (f[0] * a.0 + f[1] * cc.0 + f[2] * b.0, f[0] * a.1 + f[1] * cc.1 + f[2] * b.1);
    let mut curve = vec![];
    for i in 0..=50 {
        let naa = (1000.0 * i as f64 / 50.0).round() as i64;
        let set = [("n_AA", ii(naa)), ("n_Aa", ii(0)), ("n_aa", ii(1000 - naa))];
        if let Some(an) = s.with_input(&set) {
            if let (Some(x), Some(y), Some(z)) = (k::num_of(&an, "exp_AA"), k::num_of(&an, "exp_Aa"), k::num_of(&an, "exp_aa")) {
                let t = x + y + z;
                curve.push(pt([x / t, y / t, z / t]));
            }
        }
    }
    if curve.len() < 3 {
        return None;
    }
    let mut out = format!("<text class=\"axp-mt\" x=\"16\" y=\"22\">chi2 = {}, p-value {}</text>", tick(s.out("chi2")?),
                          tick(s.out("p_value")?));
    out.push_str(&format!("<path class=\"axp-k\" d=\"M{:.1},{:.1} L{:.1},{:.1} L{:.1},{:.1} Z\"/>", a.0, a.1, b.0, b.1, cc.0, cc.1));
    let d: String = curve.iter().enumerate().map(|(i, p)| format!("{}{:.1},{:.1} ", if i == 0 { "M" } else { "L" }, p.0, p.1)).collect();
    out.push_str(&format!("<path class=\"axp-c\" d=\"{}\"/>", d.trim_end()));
    let o = pt([obs[0] / n, obs[1] / n, obs[2] / n]);
    let et: f64 = exp.iter().sum();
    let e = pt([exp[0] / et, exp[1] / et, exp[2] / et]);
    out.push_str(&format!("<circle class=\"axp-d3\" cx=\"{:.1}\" cy=\"{:.1}\" r=\"5\"/><circle class=\"axp-m\" cx=\"{:.1}\" \
                           cy=\"{:.1}\" r=\"6\"/>", e.0, e.1, o.0, o.1));
    out.push_str(&format!("<text class=\"axp-bt\" x=\"{:.1}\" y=\"{:.1}\" text-anchor=\"end\">AA</text><text class=\"axp-bt\" \
                           x=\"{:.1}\" y=\"{:.1}\">aa</text><text class=\"axp-bt\" x=\"{:.1}\" y=\"{:.1}\" text-anchor=\"middle\">Aa</text>",
                          a.0 - 6.0, a.1 + 4.0, b.0 + 6.0, b.1 + 4.0, cc.0, cc.1 - 8.0));
    out.push_str(&format!("<circle class=\"axp-m\" cx=\"24\" cy=\"390\" r=\"6\"/><text class=\"axp-lt\" x=\"36\" y=\"395\">observed \
                           ({}, {}, {})</text><circle class=\"axp-d3\" cx=\"270\" cy=\"390\" r=\"5\"/><text class=\"axp-lt\" \
                           x=\"282\" y=\"395\">expected under Hardy-Weinberg</text>", g(obs[0], 6), g(obs[1], 6), g(obs[2], 6)));
    out.push_str("<text class=\"axp-lt\" x=\"420\" y=\"60\" text-anchor=\"middle\">curve: Hardy-Weinberg</text>");
    Some(s.figure("triangle", "De Finetti diagram: genotypes against Hardy-Weinberg",
                  &format!("De Finetti diagram: a population is a point of the triangle whose distances to the sides are its \
                            genotype frequencies (AA, Aa, aa at the corners). The curve joins the Hardy-Weinberg proportions the \
                            model returns as expected counts (populations of 1000 at {} allele frequencies). Marked: the given \
                            counts AA = {}, Aa = {}, aa = {} (observed frequencies) beside their expected point ({}, {}, {}).",
                           curve.len(), g(obs[0], 6), g(obs[1], 6), g(obs[2], 6), tick(exp[0]), tick(exp[1]), tick(exp[2])),
                  &out, 410.0))
}

pub fn hw_exact_bars(s: &Solved) -> Option<Figure> {
    let (aa, ab, bb) = (s.val("n_AA").map(|v| v.num() as i64)?, s.val("n_Aa").map(|v| v.num() as i64)?,
                        s.val("n_aa").map(|v| v.num() as i64)?);
    let n = aa + ab + bb;
    let na = 2 * aa + ab;
    let hmax = na.min(2 * n - na);
    let hs: Vec<i64> = (0..=hmax).filter(|h| (na - h) % 2 == 0).collect();
    if hs.len() < 2 || hs.len() > 60 || n == 0 {
        return None;
    }
    let mut bars = vec![];
    let mut pv = vec![];
    for h in &hs {
        let x_aa = (na - h) / 2;
        let x_bb = n - h - x_aa;
        if x_bb < 0 {
            continue;
        }
        if let Some(a) = s.with_input(&[("n_AA", ii(x_aa)), ("n_Aa", ii(*h)), ("n_aa", ii(x_bb))]) {
            if let (Some(p), Some(q)) = (k::num_of(&a, "p_observed"), k::num_of(&a, "p_value")) {
                bars.push((*h as f64, p));
                pv.push((*h as f64, q));
            }
        }
    }
    if bars.len() < 2 {
        return None;
    }
    let best = bars.iter().cloned().fold(bars[0], |a, b| if b.1 > a.1 { b } else { a });
    let xr = (bars[0].0 - 1.0, bars[bars.len() - 1].0 + 1.0);
    let mut ch = Chart::new(xr, (0.0, 1.0), (false, false), (false, false), "number of heterozygotes",
                            "probability", &format!("given Aa = {ab}: p-value {}", tick(s.out("p_value")?)));
    let w = (380.0 / bars.len() as f64 * 0.6).min(36.0);
    for (h, p) in &bars {
        chart_bars(&mut ch, &[(*h, 0.0, *p)], w, if *h as i64 == ab { "axp-b" } else { "axp-b4" });
    }
    ch.line(&pv, "axp-c2", "two-sided p-value", false);
    ch.dots(&pv, "axp-d2", 3.5);
    ch.legend.push(("axp-c4".into(), "probability (bars)".into()));
    ch.mark(best, &format!("{}: {}", g(best.0, 6), tick(best.1)));
    let body = ch.body();
    Some(s.chart("bars", "Exact Hardy-Weinberg test: every heterozygote count",
                 &format!("Exact Hardy-Weinberg test for {n} individuals carrying {na} copies of allele A: for every possible \
                           number of heterozygotes, the model's probability of that configuration (bars) and its two-sided \
                           exact p-value (dots). Too few or too many heterozygotes are both improbable. Marked: the most probable \
                           count, {}, {}. The given counts (AA {aa}, Aa {ab}, aa {bb}) give p-value {} (highlighted bar).",
                          g(best.0, 6), tick(best.1), tick(s.out("p_value")?)), &body))
}

// ------------------------------------------------------------------------------------------------ Leslie
pub fn leslie_population(s: &Solved) -> Option<Figure> {
    let (f, sv, n0) = (s.vec("F")?, s.vec("S")?, s.vec("n0")?);
    if n0.len() > 8 {
        return None;
    }
    let steps = 12usize;
    let mut cols: Vec<Vec<f64>> = vec![n0.clone()];
    for t in 1..=steps {
        let a = s.with_input(&[("t", ii(t as i64))])?;
        let v = k::list_of(&a, "n_t")?;
        if v.len() != n0.len() {
            return None;
        }
        cols.push(v);
    }
    let lam = other(s, "bio_leslie_growth_rate", &[("F", vi(&f)), ("S", vi(&sv))], "lambda1");
    let tot: Vec<f64> = cols.iter().map(|c_| c_.iter().sum()).collect();
    let hi = tot.iter().cloned().fold(0.0, f64::max);
    let mut ch = Chart::new((-0.6, steps as f64 + 0.6), (0.0, hi), (false, false), (false, true), "step t",
                            "population by age class", &format!("t = {steps}: total {}", tick(tot[steps])));
    let fills = ["axp-b", "axp-b2", "axp-b3", "axp-b4"];
    for (t, col) in cols.iter().enumerate() {
        let mut base = 0.0;
        for (i, v) in col.iter().enumerate() {
            chart_bars(&mut ch, &[(t as f64, base, base + v)], 22.0, fills[i % 4]);
            base += v;
        }
    }
    let names: Vec<String> = (0..n0.len().min(4)).map(|i| format!("class {}", i + 1)).collect();
    let items: Vec<(&str, &str)> = names.iter().enumerate().map(|(i, n)| (fills[i], n.as_str())).collect();
    ch.dots(&[(steps as f64, tot[steps])], "axp-m", 5.0);
    let body = ch.body() + &crate::figs_extra::fill_legend(&items);
    Some(s.chart("population", "Age-structured population growth (Leslie matrix)",
                 &format!("Leslie-matrix projection of an age-structured population (fecundities F = {}, survival S = {}, start \
                           n0 = {}): each bar is the population the model returns at that step, split by age class{}; at t = \
                           {steps} the classes are {}. Marked: total = {}.", list(&f), list(&sv), list(&n0),
                          lam.map(|l| format!("; growth settles to lambda = {} per step (the model's growth-rate capability)", tick(l)))
                              .unwrap_or_default(), list(&cols[steps]), tick(tot[steps])), &body))
}

pub fn leslie_trajectory(s: &Solved) -> Option<Figure> {
    let (f, sv, n0) = (s.vec("F")?, s.vec("S")?, s.vec("n0")?);
    let m = n0.len();
    if m < 2 || m > 8 {
        return None;
    }
    let total: f64 = n0.iter().sum();
    let mut young = vec![0.0; m];
    young[0] = total;
    let mut old = vec![0.0; m];
    old[m - 1] = total;
    let starts = [("given n0", n0.clone()), ("all newborn", young), ("all in the oldest class", old)];
    let steps = 20i64;
    let mut ys = vec![];
    let mut series = vec![];
    for (name, st) in &starts {
        let mut pts = vec![(0.0, total)];
        for t in 1..=steps {
            if let Some(v) = run(s, &[("n0", vi(st)), ("t", ii(t))], "total") {
                if v > 0.0 {
                    pts.push((t as f64, v));
                }
            }
        }
        ys.extend(pts.iter().map(|p| p.1));
        series.push((name.to_string(), pts));
    }
    let lam = other(s, "bio_leslie_growth_rate", &[("F", vi(&f)), ("S", vi(&sv))], "lambda1");
    let mut ch = Chart::new((0.0, steps as f64), Chart::range_pos(&ys), (false, true), (false, true), "step t",
                            "total population (logarithmic)", &lam.map(|l| format!("lambda = {}", tick(l))).unwrap_or_default());
    for (k_, (n, p)) in series.iter().enumerate() {
        ch.line(p, LINES[k_], n, false);
        ch.dots(p, DOTS[k_], 2.5);
    }
    let body = ch.body();
    Some(s.chart("trajectory", "Leslie projections converge to one growth rate",
                 &format!("Total population projected by the Leslie matrix F = {}, S = {} from three starts of {} individuals \
                           (the given n0 = {}, all newborn, all in the oldest class) on a logarithmic axis, each step a run of \
                           the model: after the starting oscillations all three grow in parallel at the dominant rate{}.",
                          list(&f), list(&sv), tick(total), list(&n0),
                          lam.map(|l| format!(" lambda = {} per step (the model's growth-rate capability)", tick(l))).unwrap_or_default()),
                 &body))
}

// ------------------------------------------------------------------------------------------------ linkage
pub fn ld_table(s: &Solved) -> Option<Figure> {
    let v = [s.val("n_AB")?.num(), s.val("n_Ab")?.num(), s.val("n_aB")?.num(), s.val("n_ab")?.num()];
    let (d, dp, r2) = (s.out("D")?, s.out("D_prime")?, s.out("r2")?);
    let m = vec![vec![Some(v[0]), Some(v[1])], vec![Some(v[2]), Some(v[3])]];
    let (body, h) = heat_table(&m, &["A".into(), "a".into()], &["B".into(), "b".into()], "first locus", "second locus",
                               &format!("D = {}, D' = {}, r2 = {}", tick(d), tick(dp), tick(r2)), None, &[], 6);
    Some(s.figure("table", "Haplotype counts of two loci",
                  &format!("Haplotype counts of two loci as a 2 x 2 table (rows: allele A or a, columns: B or b): AB = {}, Ab = \
                            {}, aB = {}, ab = {}. The excess of AB and ab over what independent loci would give is the \
                            disequilibrium the model returns: D = {}, D' = {}, r2 = {}.", g(v[0], 6), g(v[1], 6), g(v[2], 6),
                           g(v[3], 6), tick(d), tick(dp), tick(r2)), &body, h))
}
