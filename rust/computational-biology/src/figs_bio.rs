//! Computational Biology's drawing helpers beyond figs_common.rs: inputs for further solves, a table of numbers shaded
//! by value, trees read from Newick text (drawn to scale), and a dynamic-programming table of letters.
#![allow(dead_code)]

use crate::figs_common::{clip, esc, tick};
use ax_engine::figkit::{self as k, Solved};
use ax_engine::{Answer, Input};

pub fn ti(t: &str) -> Input {
    Input::Text(t.to_string())
}
pub fn fi(x: f64) -> Input {
    Input::Float(x)
}
pub fn ii(x: i64) -> Input {
    Input::Int(x)
}
pub fn vi(v: &[f64]) -> Input {
    Input::Json(serde_json::json!(v))
}
pub fn mi(m: &[Vec<f64>]) -> Input {
    Input::Json(serde_json::json!(m))
}

/// The same capability with some inputs replaced, one number output read.
pub fn run(s: &Solved, set: &[(&str, Input)], out: &str) -> Option<f64> {
    let a = s.with_input(set)?;
    k::num_of(&a, out)
}

/// Another capability of the model, one number output read.
pub fn other(s: &Solved, key: &str, set: &[(&str, Input)], out: &str) -> Option<f64> {
    let a: Answer = s.solve_other(key, set)?;
    k::num_of(&a, out)
}

/// A number written with up to `d` significant digits (integers as integers).
pub fn g(v: f64, d: usize) -> String {
    if v.is_finite() && v == v.round() && v.abs() < 1e9 { format!("{}", v as i64) } else { ax_engine::num::py_g(v, d) }
}

/// Bars on an x-y chart: (x centre, bottom, top) in data units, `w` pixels wide.
pub fn chart_bars(c: &mut crate::figs_common::Chart, bars: &[(f64, f64, f64)], w: f64, class: &str) {
    for (x, b, t) in bars {
        let (px, y1, y2) = (c.xa.px(*x), c.ya.px(*b), c.ya.px(*t));
        if !(px.is_finite() && y1.is_finite() && y2.is_finite()) {
            continue;
        }
        let (y1, y2) = (y1.clamp(crate::figs_common::CY0, crate::figs_common::CY1), y2.clamp(crate::figs_common::CY0, crate::figs_common::CY1));
        c.s.push_str(&format!("<rect class=\"{class}\" x=\"{:.1}\" y=\"{:.1}\" width=\"{w:.1}\" height=\"{:.1}\"/>",
                              px - w / 2.0, y1.min(y2), (y1 - y2).abs().max(0.5)));
        let mut yy = y1.min(y2);
        while yy < y1.max(y2) {
            c.pix.push((px, yy));
            yy += 12.0;
        }
    }
}

/// A chart for n categories across (category i centred at x = i + 0.5) and a value axis up: grid, ticks, frame.
pub fn cat_chart(n: usize, yr: (f64, f64), log: bool, ylab: &str, header: &str) -> crate::figs_common::Chart {
    use crate::figs_common::{Axis, Chart, CX0, CX1, CY0, CY1};
    let xa = Axis::new(0.0, n.max(1) as f64, false, CX0, CX1, false);
    let ya = Axis::new(yr.0, yr.1, log, CY1, CY0, true);
    let mut c = Chart { xa, ya, s: String::new(), legend: vec![], pix: vec![], h: 340.0 };
    for t in ya.ticks() {
        let y = ya.px(t);
        c.s.push_str(&format!("<line class=\"axp-g\" x1=\"{CX0}\" y1=\"{y:.1}\" x2=\"{CX1}\" y2=\"{y:.1}\"/>\
                               <text class=\"axp-t\" x=\"{}\" y=\"{:.1}\" text-anchor=\"end\">{}</text>", CX0 - 6.0, y + 5.0,
                              esc(&tick(t))));
    }
    c.s.push_str(&format!("<rect class=\"axp-f\" x=\"{CX0}\" y=\"{CY0}\" width=\"{}\" height=\"{}\"/>", CX1 - CX0, CY1 - CY0));
    c.s.push_str(&format!("<clipPath id=\"CLIPID\"><rect x=\"{CX0}\" y=\"{CY0}\" width=\"{}\" height=\"{}\"/></clipPath>",
                          CX1 - CX0, CY1 - CY0));
    if !ylab.is_empty() {
        c.s.push_str(&format!("<text class=\"axp-l\" transform=\"translate(16,{:.1}) rotate(-90)\" text-anchor=\"middle\">{}</text>",
                              (CY0 + CY1) / 2.0, esc(&clip(ylab, 32))));
    }
    if !header.is_empty() {
        c.s.push_str(&format!("<text class=\"axp-mt\" x=\"{CX0}\" y=\"22\">{}</text>", esc(&clip(header, 46))));
    }
    c
}

/// Category labels under a chart's x positions.
pub fn chart_xlabels(c: &mut crate::figs_common::Chart, labels: &[(f64, String)]) {
    for (x, l) in labels {
        let px = c.xa.px(*x);
        c.s.push_str(&format!("<text class=\"axp-t\" x=\"{px:.1}\" y=\"{}\" text-anchor=\"middle\" font-size=\"12\">{}</text>",
                              crate::figs_common::CY1 + 16.0, esc(&clip(l, 10))));
    }
}

// ------------------------------------------------------------------------------------------------ tables
/// A table of numbers, each cell shaded by its value (accent positive, second colour negative), its value written;
/// row and column labels; `mark`: one cell ringed; `outline`: cells outlined. Returns the body and its height.
#[allow(clippy::too_many_arguments)]
pub fn heat_table(m: &[Vec<Option<f64>>], rows: &[String], cols: &[String], row_title: &str, col_title: &str,
                  header: &str, mark: Option<(usize, usize)>, outline: &[(usize, usize)], digits: usize) -> (String, f64) {
    let r = m.len();
    let c = m.first().map(|x| x.len()).unwrap_or(0);
    let mut out = String::new();
    if r == 0 || c == 0 {
        return (out, 100.0);
    }
    let top = if header.is_empty() { 20.0 } else { 34.0 } + if col_title.is_empty() { 0.0 } else { 18.0 } + 18.0;
    let left = 34.0 + rows.iter().map(|x| x.chars().count()).max().unwrap_or(1) as f64 * 8.0
        + if row_title.is_empty() { 0.0 } else { 18.0 };
    let cell = ((500.0 - left) / c as f64).min(52.0).min(560.0 / r as f64).max(14.0);
    let x0 = left + ((500.0 - left) - cell * c as f64) / 2.0;
    let y0 = top;
    if !header.is_empty() {
        out.push_str(&format!("<text class=\"axp-mt\" x=\"16\" y=\"22\">{}</text>", esc(&clip(header, 54))));
    }
    if !col_title.is_empty() {
        out.push_str(&format!("<text class=\"axp-lt\" x=\"{:.1}\" y=\"{:.1}\" text-anchor=\"middle\">{}</text>",
                              x0 + cell * c as f64 / 2.0, top - 22.0, esc(col_title)));
    }
    if !row_title.is_empty() {
        out.push_str(&format!("<text class=\"axp-lt\" transform=\"translate(18,{:.1}) rotate(-90)\" text-anchor=\"middle\">{}</text>",
                              y0 + cell * r as f64 / 2.0, esc(row_title)));
    }
    let mx = m.iter().flatten().flatten().fold(0.0f64, |a, v| a.max(v.abs()));
    let font = (cell * 0.36).clamp(8.0, 13.0);
    let chars = ((cell - 2.0) / (font * 0.62)).floor().max(2.0) as usize;
    for (i, row) in m.iter().enumerate() {
        for (j, v) in row.iter().enumerate() {
            let (x, y) = (x0 + cell * j as f64, y0 + cell * i as f64);
            match v {
                Some(v) => {
                    let op = if mx > 0.0 { 0.06 + 0.7 * v.abs() / mx } else { 0.06 };
                    out.push_str(&format!("<rect class=\"{}\" x=\"{x:.1}\" y=\"{y:.1}\" width=\"{:.1}\" height=\"{:.1}\" \
                                           fill-opacity=\"{op:.2}\"/>", if *v < 0.0 { "axp-hn" } else { "axp-hp" }, cell - 1.0,
                                          cell - 1.0));
                    let t = g(*v, digits);
                    let t = if t.chars().count() > chars { ax_engine::num::py_g(*v, 2) } else { t };
                    out.push_str(&format!("<text class=\"axp-ht\" x=\"{:.1}\" y=\"{:.1}\" text-anchor=\"middle\" font-size=\"{font:.0}\">{}</text>",
                                          x + cell / 2.0, y + cell / 2.0 + font * 0.35, esc(&t)));
                }
                // a column without a label is a gap between two panels
                None if cols.get(j).map(|l| l.is_empty()).unwrap_or(false) => {}
                None => out.push_str(&format!("<rect class=\"axp-f\" x=\"{x:.1}\" y=\"{y:.1}\" width=\"{:.1}\" height=\"{:.1}\"/>",
                                              cell - 1.0, cell - 1.0)),
            }
        }
    }
    for (i, j) in outline {
        out.push_str(&format!("<rect class=\"axp-k\" x=\"{:.1}\" y=\"{:.1}\" width=\"{:.1}\" height=\"{:.1}\"/>",
                              x0 + cell * *j as f64 + 0.5, y0 + cell * *i as f64 + 0.5, cell - 2.0, cell - 2.0));
    }
    if let Some((i, j)) = mark {
        out.push_str(&format!("<rect class=\"axp-c2\" x=\"{:.1}\" y=\"{:.1}\" width=\"{:.1}\" height=\"{:.1}\" style=\"stroke-width:3\"/>",
                              x0 + cell * j as f64 - 1.0, y0 + cell * i as f64 - 1.0, cell + 1.0, cell + 1.0));
    }
    let lf = (cell * 0.5).clamp(9.0, 13.0);
    for (i, l) in rows.iter().enumerate() {
        out.push_str(&format!("<text class=\"axp-bt\" x=\"{:.1}\" y=\"{:.1}\" text-anchor=\"end\" font-size=\"{lf:.0}\">{}</text>",
                              x0 - 6.0, y0 + cell * (i as f64 + 0.5) + 4.0, esc(l)));
    }
    for (j, l) in cols.iter().enumerate() {
        out.push_str(&format!("<text class=\"axp-bt\" x=\"{:.1}\" y=\"{:.1}\" text-anchor=\"middle\" font-size=\"{lf:.0}\">{}</text>",
                              x0 + cell * (j as f64 + 0.5), y0 - 6.0, esc(l)));
    }
    (out, y0 + cell * r as f64 + 16.0)
}

// ------------------------------------------------------------------------------------------------ trees
/// A tree read from Newick text: nodes with a name, a branch length (to the parent) and children; node 0 is the root.
#[derive(Clone, Debug, Default)]
pub struct Tree {
    pub name: Vec<String>,
    pub len: Vec<Option<f64>>,
    pub kids: Vec<Vec<usize>>,
    pub parent: Vec<Option<usize>>,
}

impl Tree {
    pub fn add(&mut self, parent: Option<usize>) -> usize {
        self.name.push(String::new());
        self.len.push(None);
        self.kids.push(vec![]);
        self.parent.push(parent);
        let i = self.name.len() - 1;
        if let Some(p) = parent {
            self.kids[p].push(i);
        }
        i
    }
    pub fn leaves(&self) -> Vec<usize> {
        let mut v = vec![];
        self.walk(0, &mut |i| if self.kids[i].is_empty() { v.push(i) });
        v
    }
    fn walk(&self, i: usize, f: &mut dyn FnMut(usize)) {
        f(i);
        for &c in &self.kids[i] {
            self.walk(c, f);
        }
    }
    /// The leaves under a node.
    pub fn under(&self, i: usize) -> Vec<usize> {
        let mut v = vec![];
        self.walk(i, &mut |j| if self.kids[j].is_empty() { v.push(j) });
        v
    }
    /// Distance of each node from the root (branch lengths, 1 per branch when a length is missing and `unit`).
    pub fn depth(&self, unit: bool) -> Vec<f64> {
        let mut d = vec![0.0; self.name.len()];
        let mut stack = vec![0usize];
        while let Some(i) = stack.pop() {
            for &c in &self.kids[i] {
                d[c] = d[i] + self.len[c].unwrap_or(if unit { 1.0 } else { 0.0 });
                stack.push(c);
            }
        }
        d
    }
    pub fn leaf(&self, name: &str) -> Option<usize> {
        (0..self.name.len()).find(|&i| self.kids[i].is_empty() && self.name[i] == name)
    }
}

/// A Newick tree ("((A:1,B:2):3,C);"); None when it does not parse.
pub fn newick(t: &str) -> Option<Tree> {
    let cs: Vec<char> = t.trim().trim_end_matches(';').chars().filter(|c| !c.is_whitespace()).collect();
    let mut tr = Tree::default();
    let mut cur = tr.add(None);
    let mut i = 0;
    let read_label = |i: &mut usize| -> String {
        let s0 = *i;
        while *i < cs.len() && !"(),:;".contains(cs[*i]) {
            *i += 1;
        }
        cs[s0..*i].iter().collect()
    };
    while i < cs.len() {
        match cs[i] {
            '(' => {
                cur = tr.add(Some(cur));
                i += 1;
            }
            ',' => {
                let p = tr.parent[cur]?;
                cur = tr.add(Some(p));
                i += 1;
            }
            ')' => {
                cur = tr.parent[cur]?;
                i += 1;
                let l = read_label(&mut i);
                if !l.is_empty() {
                    tr.name[cur] = l;
                }
            }
            ':' => {
                i += 1;
                let s0 = i;
                while i < cs.len() && !"(),:;".contains(cs[i]) {
                    i += 1;
                }
                let v: String = cs[s0..i].iter().collect();
                tr.len[cur] = Some(v.parse().ok()?);
            }
            _ => {
                let l = read_label(&mut i);
                if l.is_empty() {
                    return None;
                }
                tr.name[cur] = l;
            }
        }
    }
    if tr.leaves().iter().any(|&l| tr.name[l].is_empty()) || tr.leaves().len() < 2 {
        return None;
    }
    Some(tr)
}

/// A drawn tree: rectangular, the root on the left, to scale (`scale`) or by branch count; each branch's class from
/// `class(node)` (thick when "thick"); leaf names right of the tips; branch lengths over the branches when `lengths`.
/// Returns the drawing and each leaf's tip position.
#[allow(clippy::too_many_arguments)]
pub fn draw_tree(t: &Tree, x0: f64, y0: f64, w: f64, h: f64, scale: bool, lengths: bool,
                 class: &dyn Fn(usize) -> &'static str) -> (String, Vec<(usize, f64, f64)>) {
    let d = t.depth(!scale);
    let leaves = t.leaves();
    let maxd = d.iter().cloned().fold(0.0f64, f64::max).max(1e-12);
    let mut y = vec![0.0; t.name.len()];
    let n = leaves.len().max(2) as f64;
    for (k, &l) in leaves.iter().enumerate() {
        y[l] = y0 + h * k as f64 / (n - 1.0);
    }
    // internal nodes: the middle of their children, children first
    let mut order = vec![];
    t.walk(0, &mut |i| order.push(i));
    for &i in order.iter().rev() {
        if !t.kids[i].is_empty() {
            let ys: Vec<f64> = t.kids[i].iter().map(|c| y[*c]).collect();
            y[i] = (ys.iter().cloned().fold(f64::INFINITY, f64::min) + ys.iter().cloned().fold(f64::NEG_INFINITY, f64::max)) / 2.0;
        }
    }
    let px = |v: f64| x0 + w * v / maxd;
    let mut out = String::new();
    let mut tips = vec![];
    for &i in &order {
        let x = px(d[i]);
        if let Some(p) = t.parent[i] {
            let xp = px(d[p]);
            let cl = class(i);
            let style = if cl == "thick" { "class=\"axp-c\" style=\"stroke-width:4\"" } else if cl == "second" {
                "class=\"axp-c2\""
            } else {
                "class=\"axp-k\""
            };
            out.push_str(&format!("<path {style} d=\"M{xp:.1},{:.1} L{xp:.1},{:.1} L{x:.1},{:.1}\"/>", y[p], y[i], y[i]));
            if lengths {
                if let Some(l) = t.len[i] {
                    out.push_str(&format!("<text class=\"axp-vt\" x=\"{:.1}\" y=\"{:.1}\" text-anchor=\"middle\" font-size=\"11\">{}</text>",
                                          (xp + x) / 2.0, y[i] - 5.0, esc(&g(l, 4))));
                }
            }
        }
        if t.kids[i].is_empty() {
            out.push_str(&format!("<text class=\"axp-bt\" x=\"{:.1}\" y=\"{:.1}\">{}</text>", x + 6.0, y[i] + 5.0,
                                  esc(&clip(&t.name[i], 10))));
            tips.push((i, x, y[i]));
        } else {
            out.push_str(&format!("<circle class=\"axp-kd\" cx=\"{x:.1}\" cy=\"{:.1}\" r=\"2.5\"/>", y[i]));
            tips.push((i, x, y[i]));
        }
    }
    (out, tips)
}

/// The splits of a tree's internal branches: each the sorted leaf names on one side, normalised to the side without
/// the alphabetically first leaf; with the node each split belongs to.
pub fn splits(t: &Tree) -> Vec<(Vec<String>, usize)> {
    let mut all: Vec<String> = t.leaves().iter().map(|&l| t.name[l].clone()).collect();
    all.sort();
    let mut out = vec![];
    for i in 1..t.name.len() {
        if t.kids[i].is_empty() {
            continue;
        }
        let mut side: Vec<String> = t.under(i).iter().map(|&l| t.name[l].clone()).collect();
        side.sort();
        if side.len() < 2 || side.len() + 2 > all.len() {
            continue;
        }
        if side.contains(&all[0]) {
            side = all.iter().filter(|n| !side.contains(n)).cloned().collect();
        }
        if !out.iter().any(|(s, _): &(Vec<String>, usize)| *s == side) {
            out.push((side, i));
        }
    }
    out
}
