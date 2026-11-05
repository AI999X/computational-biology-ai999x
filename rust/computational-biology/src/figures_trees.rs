//! Trees, pedigrees, genealogies and networks: drawn from the inputs as read, each drawn only when what it shows
//! reproduces the numbers the model returns (a tree's lengths, heights, splits; a parsimony score site by site).

use crate::figs_bio::{draw_tree, g, heat_table, newick, run, splits, ti, Tree};
use crate::figs_common::{self as c, clip, esc, tick};
use crate::figs_extra::bars_full;
use ax_engine::figkit::{self as k, Solved};
use ax_engine::figures::Figure;
use ax_engine::Val;

fn names(t: &str) -> Vec<String> {
    t.split(|c: char| c.is_whitespace() || c == ',').filter(|x| !x.is_empty()).map(String::from).collect()
}

// ------------------------------------------------------------------------------------------------ Faith's PD
/// The branches joining some leaves to the root.
fn path_branches(t: &Tree, leaves: &[usize]) -> Vec<usize> {
    let mut v: Vec<usize> = vec![];
    for &l in leaves {
        let mut i = l;
        while let Some(p) = t.parent[i] {
            if !v.contains(&i) {
                v.push(i);
            }
            i = p;
        }
    }
    v
}

pub fn faith_tree(s: &Solved) -> Option<Figure> {
    let src = s.txt("tree")?;
    let t = newick(&src)?;
    let taxa = names(&s.txt("taxa")?);
    let pd = s.out("pd")?;
    let ls: Vec<usize> = taxa.iter().map(|n| t.leaf(n)).collect::<Option<Vec<usize>>>()?;
    let br = path_branches(&t, &ls);
    let sum: f64 = br.iter().map(|i| t.len[*i].unwrap_or(0.0)).sum();
    if (sum - pd).abs() > 1e-9 * pd.abs().max(1.0) || t.leaves().len() > 16 {
        return None;
    }
    let h = 22.0 * t.leaves().len() as f64;
    let (body, _) = draw_tree(&t, 40.0, 44.0, 380.0, h, true, true, &|i| if br.contains(&i) { "thick" } else { "" });
    let mut out = format!("<text class=\"axp-mt\" x=\"16\" y=\"22\">PD of {} = {}</text>", esc(&clip(&taxa.join(" "), 20)), g(pd, 6));
    out.push_str(&body);
    Some(s.figure("tree", "Faith's phylogenetic diversity on a tree",
                  &format!("The tree {} drawn to scale (branch lengths written on the branches). The branches joining the taxa \
                            {} to the root are drawn thick: their summed length is Faith's phylogenetic diversity, which the \
                            model returns as pd = {}.", src, taxa.join(", "), g(pd, 6)), &out, 44.0 + h + 30.0))
}

pub fn faith_bars(s: &Solved) -> Option<Figure> {
    let src = s.txt("tree")?;
    let t = newick(&src)?;
    let leaves: Vec<String> = t.leaves().iter().map(|l| t.name[*l].clone()).collect();
    if leaves.len() > 8 {
        return None;
    }
    let mine = names(&s.txt("taxa")?);
    let mut labels = vec![];
    let mut vals = vec![];
    let mut hl = vec![];
    for i in 0..leaves.len() {
        for j in i + 1..leaves.len() {
            if let Some(v) = run(s, &[("taxa", ti(&format!("{} {}", leaves[i], leaves[j])))], "pd") {
                if mine.len() == 2 && mine.contains(&leaves[i]) && mine.contains(&leaves[j]) {
                    hl.push(labels.len());
                }
                labels.push(format!("{}{}", leaves[i], leaves[j]));
                vals.push(v);
            }
        }
    }
    let all = run(s, &[("taxa", ti(&leaves.join(" ")))], "pd")?;
    if labels.len() < 2 {
        return None;
    }
    let classes: Vec<&str> = (0..labels.len()).map(|i| if hl.contains(&i) { "axp-b" } else { "axp-b4" }).collect();
    let mut body = bars_full(&labels, &vals, &classes, None, "phylogenetic diversity", &format!("all {} taxa: {}", leaves.len(), g(all, 6)),
                             false, None)?;
    body.push_str(&format!("<text class=\"axp-lt\" x=\"292\" y=\"326\" text-anchor=\"middle\">pair of taxa</text>"));
    Some(s.chart("bars", "Phylogenetic diversity kept by each pair of taxa",
                 &format!("Faith's phylogenetic diversity of every pair of taxa of the tree {}, each bar a run of the model: a \
                           pair from different clades keeps more of the tree than two sister taxa. All {} taxa together: {}.{}",
                          src, leaves.len(), g(all, 6), if hl.is_empty() { String::new() } else { " Highlighted: the taxa given.".into() }),
                 &body))
}

// ------------------------------------------------------------------------------------------------ Fitch
pub fn fitch_tree(s: &Solved) -> Option<Figure> {
    let src = s.txt("tree")?;
    let t = newick(&src)?;
    let st = s.txt("states")?;
    let score = s.out("score")?;
    let mut states: Vec<(String, Vec<char>)> = vec![];
    for item in st.split_whitespace() {
        let (n, v) = item.split_once(':')?;
        states.push((n.to_string(), v.chars().collect()));
    }
    let nsite = states.first()?.1.len();
    if nsite == 0 || nsite > 20 || states.iter().any(|x| x.1.len() != nsite) || t.leaves().len() > 12 {
        return None;
    }
    let mut per = vec![];
    for i in 0..nsite {
        let one: Vec<String> = states.iter().map(|(n, v)| format!("{n}:{}", v[i])).collect();
        per.push(run(s, &[("states", ti(&one.join(" ")))], "score")?);
    }
    if (per.iter().sum::<f64>() - score).abs() > 1e-9 {
        return None;
    }
    let nl = t.leaves().len();
    let h = 26.0 * (nl as f64 - 1.0);
    let (mut body, tips) = draw_tree(&t, 30.0, 46.0, 180.0, h, false, false, &|_| "");
    let cell = (240.0 / nsite as f64).min(26.0);
    let x0 = 260.0;
    let cls = |c: char| match c { 'A' => "axp-b", 'C' => "axp-b2", 'G' => "axp-b3", 'T' | 'U' => "axp-b4", _ => "axp-ka" };
    for (leaf, _x, y) in tips.into_iter().filter(|p| t.kids[p.0].is_empty()) {
        let v = &states.iter().find(|x| x.0 == t.name[leaf])?.1;
        for (j, ch) in v.iter().enumerate() {
            let x = x0 + cell * j as f64;
            body.push_str(&format!("<rect class=\"{}\" x=\"{x:.1}\" y=\"{:.1}\" width=\"{:.1}\" height=\"{:.1}\"/><text class=\"axp-bt\" \
                                    x=\"{:.1}\" y=\"{:.1}\" text-anchor=\"middle\">{ch}</text>", cls(*ch), y - cell / 2.0 + 1.0,
                                   cell - 2.0, cell - 2.0, x + cell / 2.0 - 1.0, y + 5.0));
        }
    }
    let yb = 46.0 + h + 30.0;
    body.push_str(&format!("<text class=\"axp-lt\" x=\"{:.1}\" y=\"{yb:.1}\" text-anchor=\"end\">changes per site</text>", x0 - 8.0));
    for (j, p) in per.iter().enumerate() {
        body.push_str(&format!("<text class=\"axp-vt\" x=\"{:.1}\" y=\"{yb:.1}\" text-anchor=\"middle\">{}</text>",
                               x0 + cell * j as f64 + cell / 2.0 - 1.0, g(*p, 4)));
    }
    let mut out = format!("<text class=\"axp-mt\" x=\"16\" y=\"22\">parsimony score = {}</text>", g(score, 6));
    out.push_str(&body);
    Some(s.figure("tree", "Fitch parsimony: changes needed at each site",
                  &format!("The tree {} with the character states of each taxon at {nsite} sites ({st}); letters of a site are \
                            coloured by state. Under each site, the minimum number of changes the model returns for that site \
                            alone ({}); their sum is the Fitch parsimony score of all sites, score = {}.", src,
                           per.iter().map(|p| g(*p, 4)).collect::<Vec<_>>().join(", "), g(score, 6)), &out, yb + 16.0))
}

// ------------------------------------------------------------------------------------------------ distance trees
struct Built {
    kids: Vec<Vec<(usize, f64)>>,
    name: Vec<String>,
}

impl Built {
    fn to_tree(&self, root: usize) -> Tree {
        let mut t = Tree::default();
        let r = t.add(None);
        t.name[r] = self.name[root].clone();
        let mut stack = vec![(root, r)];
        while let Some((b, ti_)) = stack.pop() {
            for (c_, l) in &self.kids[b] {
                let n = t.add(Some(ti_));
                t.name[n] = self.name[*c_].clone();
                t.len[n] = Some(*l);
                stack.push((*c_, n));
            }
        }
        t
    }
}

fn taxon(i: usize) -> String {
    if i < 26 { ((b'A' + i as u8) as char).to_string() } else { format!("T{}", i + 1) }
}

/// Neighbour joining: the tree (rooted at the last joining point) and each taxon's terminal branch.
fn neighbor_join(d: &[Vec<f64>]) -> Option<(Built, usize, Vec<f64>)> {
    let n = d.len();
    if n < 3 || n > 12 || d.iter().any(|r| r.len() != n) {
        return None;
    }
    let mut b = Built { kids: vec![vec![]; n], name: (0..n).map(taxon).collect() };
    let mut dist: Vec<Vec<f64>> = d.to_vec();
    let mut active: Vec<usize> = (0..n).collect();
    let mut term = vec![0.0; n];
    while active.len() > 3 {
        let m = active.len();
        let r: Vec<f64> = active.iter().map(|&i| active.iter().map(|&j| dist[i][j]).sum()).collect();
        let mut best = (f64::INFINITY, 0, 1);
        for a in 0..m {
            for bb in a + 1..m {
                let q = (m as f64 - 2.0) * dist[active[a]][active[bb]] - r[a] - r[bb];
                if q < best.0 - 1e-12 {
                    best = (q, a, bb);
                }
            }
        }
        let (_, a, bb) = best;
        let (i, j) = (active[a], active[bb]);
        let li = dist[i][j] / 2.0 + (r[a] - r[bb]) / (2.0 * (m as f64 - 2.0));
        let lj = dist[i][j] - li;
        let u = b.kids.len();
        b.kids.push(vec![(i, li), (j, lj)]);
        b.name.push(String::new());
        for row in dist.iter_mut() {
            row.push(0.0);
        }
        dist.push(vec![0.0; u + 1]);
        for &k_ in &active {
            if k_ != i && k_ != j {
                let v = (dist[i][k_] + dist[j][k_] - dist[i][j]) / 2.0;
                dist[u][k_] = v;
                dist[k_][u] = v;
            }
        }
        if i < n {
            term[i] = li;
        }
        if j < n {
            term[j] = lj;
        }
        active.retain(|x| *x != i && *x != j);
        active.push(u);
    }
    let (x, y, z) = (active[0], active[1], active[2]);
    let lx = (dist[x][y] + dist[x][z] - dist[y][z]) / 2.0;
    let ly = dist[x][y] - lx;
    let lz = dist[x][z] - lx;
    let root = b.kids.len();
    b.kids.push(vec![(x, lx), (y, ly), (z, lz)]);
    b.name.push(String::new());
    for (k_, l) in [(x, lx), (y, ly), (z, lz)] {
        if k_ < n {
            term[k_] = l;
        }
    }
    Some((b, root, term))
}

pub fn nj_tree(s: &Solved) -> Option<Figure> {
    let d = s.mat("D")?;
    let br = s.out_list("branch")?;
    let tl = s.out("total_length")?;
    let (b, root, term) = neighbor_join(&d)?;
    let tree = b.to_tree(root);
    let total: f64 = tree.len.iter().flatten().sum();
    if br.len() != term.len() || br.iter().zip(&term).any(|(a, x)| (a - x).abs() > 1e-6 * a.abs().max(1.0))
        || (total - tl).abs() > 1e-6 * tl.abs().max(1.0) || tree.len.iter().flatten().any(|l| *l < 0.0) {
        return None;
    }
    let nl = d.len();
    let h = 26.0 * (nl as f64 - 1.0);
    let (body, _) = draw_tree(&tree, 50.0, 46.0, 380.0, h, true, true, &|i| if tree.kids[i].is_empty() { "" } else { "second" });
    let mut out = format!("<text class=\"axp-mt\" x=\"16\" y=\"22\">total length = {}</text>", g(tl, 6));
    out.push_str(&body);
    Some(s.figure("tree", "Neighbour-joining tree with its branch lengths",
                  &format!("Neighbour-joining tree of {nl} taxa (A, B, ... the rows of D), drawn to scale from the last joining \
                            point (an unrooted tree has no root; the dot at the left is where its three last branches meet). \
                            The model returns the terminal branch lengths {} and the total length {}; the topology and the \
                            internal branches (second colour) are the neighbour-joining pairing read from D, which reproduces \
                            both.", br.iter().map(|x| g(*x, 6)).collect::<Vec<_>>().join(", "), g(tl, 6)), &out, 46.0 + h + 30.0))
}

fn matrix_fig(s: &Solved, title: &str, extra: Option<(&str, Vec<f64>)>, note: &str) -> Option<Figure> {
    let d = s.mat("D")?;
    let n = d.len();
    if n > 12 {
        return None;
    }
    let mut m: Vec<Vec<Option<f64>>> = d.iter().map(|r| r.iter().map(|v| Some(*v)).collect()).collect();
    let mut cols: Vec<String> = (1..=n).map(|i| i.to_string()).collect();
    if let Some((name, v)) = &extra {
        if v.len() != n {
            return None;
        }
        for (i, row) in m.iter_mut().enumerate() {
            row.push(None);
            row.push(Some(v[i]));
        }
        cols.push(String::new());
        cols.push(name.to_string());
    }
    let rows: Vec<String> = (1..=n).map(|i| i.to_string()).collect();
    let (body, h) = heat_table(&m, &rows, &cols, "taxon", "", note, None, &[], 4);
    Some(s.figure("matrix", title,
                  &format!("The matrices of the inputs given as a heat map (rows and columns numbered; each cell shaded by its \
                            magnitude): D (input){}. {note}.", extra.as_ref().map(|e| format!(" and {} (returned, the last \
                            column)", e.0)).unwrap_or_default()), &body, h))
}

pub fn nj_matrix(s: &Solved) -> Option<Figure> {
    let br = s.out_list("branch")?;
    matrix_fig(s, "Neighbour-joining branch lengths: the matrices", Some(("branch", br)),
               &format!("Returned: total_length = {}", g(s.out("total_length")?, 6)))
}

pub fn upgma_matrix(s: &Solved) -> Option<Figure> {
    matrix_fig(s, "UPGMA tree heights: the matrix", None, &format!("Returned: root_height = {}", g(s.out("root_height")?, 6)))
}

/// UPGMA: the tree and its node heights (ascending).
fn upgma(d: &[Vec<f64>]) -> Option<(Built, usize, Vec<f64>)> {
    let n = d.len();
    if n < 2 || n > 12 || d.iter().any(|r| r.len() != n) {
        return None;
    }
    let mut b = Built { kids: vec![vec![]; n], name: (0..n).map(taxon).collect() };
    let mut height = vec![0.0; n];
    let mut size = vec![1.0; n];
    let mut dist: Vec<Vec<f64>> = d.to_vec();
    let mut active: Vec<usize> = (0..n).collect();
    let mut hs = vec![];
    while active.len() > 1 {
        let mut best = (f64::INFINITY, 0, 1);
        for a in 0..active.len() {
            for bb in a + 1..active.len() {
                let v = dist[active[a]][active[bb]];
                if v < best.0 - 1e-12 {
                    best = (v, a, bb);
                }
            }
        }
        let (v, a, bb) = best;
        let (i, j) = (active[a], active[bb]);
        let hu = v / 2.0;
        let u = b.kids.len();
        b.kids.push(vec![(i, hu - height[i]), (j, hu - height[j])]);
        b.name.push(String::new());
        height.push(hu);
        size.push(size[i] + size[j]);
        hs.push(hu);
        for row in dist.iter_mut() {
            row.push(0.0);
        }
        dist.push(vec![0.0; u + 1]);
        for &k_ in &active {
            if k_ != i && k_ != j {
                let x = (dist[i][k_] * size[i] + dist[j][k_] * size[j]) / (size[i] + size[j]);
                dist[u][k_] = x;
                dist[k_][u] = x;
            }
        }
        active.retain(|x| *x != i && *x != j);
        active.push(u);
    }
    hs.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let root = b.kids.len() - 1;
    Some((b, root, hs))
}

pub fn upgma_tree(s: &Solved) -> Option<Figure> {
    let d = s.mat("D")?;
    let hs = s.out_list("heights")?;
    let rh = s.out("root_height")?;
    let (b, root, mine) = upgma(&d)?;
    if hs.len() != mine.len() || hs.iter().zip(&mine).any(|(a, x)| (a - x).abs() > 1e-6 * a.abs().max(1.0)) {
        return None;
    }
    let tree = b.to_tree(root);
    let nl = d.len();
    let h = 26.0 * (nl as f64 - 1.0);
    let (body, pos) = draw_tree(&tree, 50.0, 46.0, 380.0, h, true, false, &|_| "");
    let mut out = format!("<text class=\"axp-mt\" x=\"16\" y=\"22\">root height = {}</text>", g(rh, 6));
    out.push_str(&body);
    // each node's height (root height minus its distance from the root), written at the node
    let depth = tree.depth(false);
    for (i, x, y) in pos.iter().filter(|p| !tree.kids[p.0].is_empty()) {
        out.push_str(&format!("<circle class=\"axp-d\" cx=\"{x:.1}\" cy=\"{y:.1}\" r=\"4\"/><text class=\"axp-vt\" x=\"{:.1}\" \
                               y=\"{:.1}\" text-anchor=\"end\">{}</text>", x - 6.0, y - 6.0, g(rh - depth[*i], 6)));
    }
    let yb = 46.0 + h + 34.0;
    out.push_str(&format!("<text class=\"axp-lt\" x=\"16\" y=\"{yb:.1}\">node heights returned: {}</text>",
                          esc(&hs.iter().map(|x| g(*x, 6)).collect::<Vec<_>>().join(", "))));
    Some(s.figure("tree", "UPGMA tree with its node heights",
                  &format!("UPGMA tree of {nl} taxa (named A, B, ... in the order of the rows of D), drawn to scale (each \
                            node at its height, written beside it; the tips level). The model returns the node heights {} and the root height \
                            {}; which clusters join at each node is the UPGMA pairing read from D, whose heights match the \
                            returned ones.", hs.iter().map(|x| g(*x, 6)).collect::<Vec<_>>().join(", "), g(rh, 6)), &out, yb + 14.0))
}

// ------------------------------------------------------------------------------------------------ Robinson-Foulds
pub fn rf_trees(s: &Solved) -> Option<Figure> {
    let (s1, s2) = (s.txt("tree1")?, s.txt("tree2")?);
    let (t1, t2) = (newick(&s1)?, newick(&s2)?);
    let rf = s.out("rf")?;
    let (a, b) = (splits(&t1), splits(&t2));
    let only1: Vec<usize> = a.iter().filter(|x| !b.iter().any(|y| y.0 == x.0)).map(|x| x.1).collect();
    let only2: Vec<usize> = b.iter().filter(|x| !a.iter().any(|y| y.0 == x.0)).map(|x| x.1).collect();
    if (only1.len() + only2.len()) as f64 != rf || t1.leaves().len() > 14 {
        return None;
    }
    let n = t1.leaves().len().max(t2.leaves().len());
    let h = 24.0 * (n as f64 - 1.0);
    let (b1, _) = draw_tree(&t1, 20.0, 64.0, 190.0, h, false, false,
                            &|i| if only1.contains(&i) { "thick" } else if t1.kids[i].is_empty() { "" } else { "second" });
    let (b2, _) = draw_tree(&t2, 280.0, 64.0, 190.0, h, false, false,
                            &|i| if only2.contains(&i) { "thick" } else if t2.kids[i].is_empty() { "" } else { "second" });
    let mut out = format!("<text class=\"axp-mt\" x=\"16\" y=\"22\">rf = {}, normalised {}</text>", g(rf, 6),
                          tick(s.out("rf_normalized")?));
    out.push_str("<text class=\"axp-bt\" x=\"20\" y=\"48\">tree 1</text><text class=\"axp-bt\" x=\"280\" y=\"48\">tree 2</text>");
    out.push_str(&b1);
    out.push_str(&b2);
    let yb = 64.0 + h + 34.0;
    out.push_str(&format!("<path class=\"axp-c\" style=\"stroke-width:4\" d=\"M20,{:.1} L44,{:.1}\"/><text class=\"axp-lt\" \
                           x=\"50\" y=\"{yb:.1}\">split missing from the other tree</text><path class=\"axp-c2\" d=\"M290,{:.1} \
                           L314,{:.1}\"/><text class=\"axp-lt\" x=\"320\" y=\"{yb:.1}\">shared split</text>", yb - 4.0, yb - 4.0,
                          yb - 4.0, yb - 4.0));
    Some(s.figure("trees", "Robinson-Foulds distance between two trees",
                  &format!("Tree 1 = {s1} and tree 2 = {s2} drawn side by side. Each internal branch splits the taxa in two; \
                            branches whose split the other tree lacks are drawn thick (read from the two trees), shared ones in \
                            the second colour. The model counts them: rf = {}, normalised {}.", g(rf, 6),
                           tick(s.out("rf_normalized")?)), &out, yb + 14.0))
}

// ------------------------------------------------------------------------------------------------ pedigree
fn edges(s: &Solved) -> Option<Vec<(String, String)>> {
    match s.val("pedigree")? {
        Val::Edges(e) => Some(e.into_iter().map(|(a, b, _)| (a, b)).collect()),
        Val::Text(t) => Some(t.split(',').filter_map(|x| x.trim().split_once('-').map(|(a, b)| (a.trim().to_string(), b.trim().to_string())))
            .collect()),
        _ => None,
    }
}

fn ancestors(e: &[(String, String)], x: &str) -> Vec<String> {
    let mut out: Vec<String> = vec![];
    let mut stack = vec![x.to_string()];
    while let Some(n) = stack.pop() {
        for (p, c_) in e {
            if *c_ == n && !out.contains(p) {
                out.push(p.clone());
                stack.push(p.clone());
            }
        }
    }
    out
}

pub fn pedigree(s: &Solved) -> Option<Figure> {
    let e = edges(s)?;
    let x = s.txt("individual")?;
    let f = s.out("F")?;
    let mut nodes: Vec<String> = vec![];
    for (a, b) in &e {
        for n in [a, b] {
            if !nodes.contains(n) {
                nodes.push(n.clone());
            }
        }
    }
    if nodes.len() > 18 || !nodes.contains(&x) {
        return None;
    }
    // generation: the longest chain of parents above a node
    let mut gen = vec![0usize; nodes.len()];
    for _ in 0..nodes.len() {
        for (a, b) in &e {
            let (ia, ib) = (nodes.iter().position(|n| n == a)?, nodes.iter().position(|n| n == b)?);
            if gen[ib] < gen[ia] + 1 {
                gen[ib] = gen[ia] + 1;
            }
        }
    }
    let ng = gen.iter().max().copied().unwrap_or(0) + 1;
    if ng > 6 {
        return None;
    }
    let parents: Vec<String> = e.iter().filter(|(_, c_)| *c_ == x).map(|(p, _)| p.clone()).collect();
    let shared: Vec<String> = if parents.len() == 2 {
        let a1 = ancestors(&e, &parents[0]);
        let a2 = ancestors(&e, &parents[1]);
        a1.into_iter().filter(|a| a2.contains(a)).collect()
    } else {
        vec![]
    };
    let mut pos = vec![(0.0, 0.0); nodes.len()];
    for gi in 0..ng {
        let row: Vec<usize> = (0..nodes.len()).filter(|i| gen[*i] == gi).collect();
        for (k_, i) in row.iter().enumerate() {
            pos[*i] = (520.0 * (k_ as f64 + 1.0) / (row.len() as f64 + 1.0), 60.0 + 200.0 * gi as f64 / (ng as f64 - 1.0).max(1.0));
        }
    }
    let mut out = format!("<text class=\"axp-mt\" x=\"16\" y=\"22\">F({}) = {}</text>", esc(&x), tick(f));
    for (a, b) in &e {
        let (ia, ib) = (nodes.iter().position(|n| n == a)?, nodes.iter().position(|n| n == b)?);
        out.push_str(&format!("<path class=\"axp-k\" d=\"M{:.1},{:.1} L{:.1},{:.1}\"/>", pos[ia].0, pos[ia].1 + 12.0, pos[ib].0,
                              pos[ib].1 - 12.0));
    }
    for (i, n) in nodes.iter().enumerate() {
        let cl = if *n == x { "axp-b" } else if shared.contains(n) { "axp-b2" } else { "axp-ka" };
        let w = (n.chars().count() as f64 * 8.0 + 16.0).max(28.0);
        out.push_str(&format!("<rect class=\"{cl}\" x=\"{:.1}\" y=\"{:.1}\" width=\"{w:.1}\" height=\"24\" rx=\"5\"/>\
                               <rect class=\"axp-k\" x=\"{:.1}\" y=\"{:.1}\" width=\"{w:.1}\" height=\"24\" rx=\"5\"/>\
                               <text class=\"axp-bt\" x=\"{:.1}\" y=\"{:.1}\" text-anchor=\"middle\">{}</text>", pos[i].0 - w / 2.0,
                              pos[i].1 - 12.0, pos[i].0 - w / 2.0, pos[i].1 - 12.0, pos[i].0, pos[i].1 + 5.0, esc(n)));
    }
    out.push_str("<rect class=\"axp-b\" x=\"16\" y=\"288\" width=\"12\" height=\"12\"/><text class=\"axp-lt\" x=\"32\" y=\"299\">the \
                  individual</text><rect class=\"axp-b2\" x=\"140\" y=\"288\" width=\"12\" height=\"12\"/><text class=\"axp-lt\" \
                  x=\"156\" y=\"299\">ancestors shared by its two parents</text>");
    Some(s.figure("pedigree", &format!("Pedigree of {x} with its inbreeding coefficient"),
                  &format!("The pedigree {} drawn by generation, each line joining a parent to its child; the individual {x} is \
                            shaded, and the ancestors shared by its two parents (read from the pedigree: {}) in a second colour. \
                            The model returns the inbreeding coefficient F = {}.",
                           e.iter().map(|(a, b)| format!("{a}-{b}")).collect::<Vec<_>>().join(", "),
                           if shared.is_empty() { "none".into() } else { shared.join(", ") }, tick(f)), &out, 310.0))
}

pub fn pedigree_bars(s: &Solved) -> Option<Figure> {
    let e = edges(s)?;
    let x = s.txt("individual")?;
    let mine: Vec<String> = e.iter().map(|(a, b)| format!("{a}-{b}")).collect();
    let cases = [("parent offspring", "P-S, Q-S, P-X, S-X"), ("full sibs", "P-S, Q-S, P-D, Q-D, S-X, D-X"),
                 ("half sibs", "P-S, Q-S, P-D, R-D, S-X, D-X"), ("uncle niece", "P-A, Q-A, P-B, Q-B, B-C, M-C, A-X, C-X"),
                 ("first cousins", "P-A, Q-A, P-B, Q-B, A-C, M-C, B-D, N-D, C-X, D-X")];
    let mut labels = vec![];
    let mut vals = vec![];
    let mut hl = vec![];
    let norm = |t: &str| -> Vec<String> {
        let mut v: Vec<String> = t.split(',').map(|x| x.trim().replace(' ', "")).collect();
        v.sort();
        v
    };
    let mut mine_n = mine.clone();
    mine_n.sort();
    for (name, ped) in cases {
        if let Some(v) = run(s, &[("pedigree", ti(ped)), ("individual", ti("X"))], "F") {
            if x == "X" && norm(ped) == mine_n {
                hl.push(labels.len());
            }
            labels.push(name.replacen(' ', "\n", 1));
            vals.push(v);
        }
    }
    if hl.is_empty() {
        hl.push(labels.len());
        labels.push("given".into());
        vals.push(s.out("F")?);
    }
    let classes: Vec<&str> = (0..labels.len()).map(|i| if hl.contains(&i) { "axp-b" } else { "axp-b4" }).collect();
    let body = bars_full(&labels, &vals, &classes, None, "inbreeding coefficient F", &format!("given pedigree: F = {}",
                         tick(s.out("F")?)), false, None)?;
    Some(s.chart("bars", "Inbreeding coefficient for kinds of related parents",
                 &format!("Inbreeding coefficient of an individual X whose parents are related in five ways, each a run of the \
                           model on the pedigree written: {}. Highlighted: the given pedigree ({}).",
                          cases.iter().map(|c_| format!("{}: {}", c_.0, c_.1)).collect::<Vec<_>>().join("; "), mine.join(", ")),
                 &body))
}

// ------------------------------------------------------------------------------------------------ coalescent
pub fn genealogy(s: &Solved) -> Option<Figure> {
    let ne = s.inp("Ne")?;
    let tm = s.out("t_mrca")?;
    let n = 8usize;
    let mut tk = vec![0.0; n + 1];
    for (kk, slot) in tk.iter_mut().enumerate().skip(2) {
        *slot = run(s, &[("n", crate::figs_bio::ii(kk as i64))], "t_mrca")?;
    }
    let total = tk[n];
    if !(total > 0.0) {
        return None;
    }
    let (y0, y1) = (290.0, 50.0);
    let ty = |t: f64| y0 + (y1 - y0) * t / total;
    let mut out = format!("<text class=\"axp-mt\" x=\"16\" y=\"22\">Ne = {}: T(8) = {} generations</text>", tick(ne), tick(total));
    // lineages: x positions, merging two neighbours at each level
    let mut lin: Vec<(f64, f64)> = (0..n).map(|i| (110.0 + 340.0 * i as f64 / (n as f64 - 1.0), 0.0)).collect();
    let mut ev = 0usize;
    let mut last_label = f64::INFINITY;
    for j in (2..=n).rev() {
        let t = total - tk[j - 1];
        let y = ty(t);
        out.push_str(&format!("<path class=\"axp-c4\" stroke-dasharray=\"5 4\" d=\"M90,{y:.1} L470,{y:.1}\"/>"));
        // a level too close to the one labelled below it keeps its line only
        if (last_label - y).abs() >= 13.0 {
            out.push_str(&format!("<text class=\"axp-t\" x=\"86\" y=\"{:.1}\" text-anchor=\"end\" font-size=\"11\">{}</text>",
                                  y + 4.0, tick(t)));
            last_label = y;
        }
        let m = lin.len();
        let k_ = (ev * 2 + 1) % (m - 1);
        ev += 1;
        let (a, b) = (lin[k_], lin[k_ + 1]);
        for (x, t0) in [a, b] {
            out.push_str(&format!("<path class=\"axp-k\" d=\"M{x:.1},{:.1} L{x:.1},{y:.1}\"/>", ty(t0)));
        }
        out.push_str(&format!("<path class=\"axp-k\" d=\"M{:.1},{y:.1} L{:.1},{y:.1}\"/>", a.0, b.0));
        lin.remove(k_ + 1);
        lin[k_] = ((a.0 + b.0) / 2.0, t);
    }
    let root = lin[0];
    out.push_str(&format!("<path class=\"axp-k\" d=\"M{:.1},{:.1} L{:.1},{:.1}\"/><circle class=\"axp-m\" cx=\"{:.1}\" cy=\"{:.1}\" \
                           r=\"6\"/><text class=\"axp-mt\" x=\"{:.1}\" y=\"{:.1}\">MRCA</text>", root.0, ty(root.1), root.0,
                          ty(root.1) - 10.0, root.0, ty(root.1), root.0 + 10.0, ty(root.1) - 4.0));
    for i in 0..n {
        out.push_str(&format!("<circle class=\"axp-kd\" cx=\"{:.1}\" cy=\"{y0}\" r=\"3\"/>", 110.0 + 340.0 * i as f64 / (n as f64 - 1.0)));
    }
    out.push_str("<text class=\"axp-lt\" x=\"280\" y=\"312\" text-anchor=\"middle\">8 gene copies sampled now</text>\
                  <text class=\"axp-lt\" x=\"16\" y=\"40\">generations back</text>");
    let (nb, _) = c::notes(&[("given:".into(), format!("n = {}, t_mrca = {} generations", g(s.val("n")?.num(), 6), tick(tm)))], 334.0);
    out.push_str(&nb);
    Some(s.figure("genealogy", "Expected coalescent genealogy of eight gene copies",
                  &format!("Expected genealogy of 8 gene copies in a population of effective size Ne = {}: the dashed levels are \
                            the expected times of the successive coalescences, T(8) - T(j - 1), from the expected times to the \
                            most recent common ancestor T(k) the model returns for k = 2 to 8 (T(1) = 0). Which lineages join at \
                            each level is drawn arbitrarily. The given n has t_mrca = {}.", tick(ne), tick(tm)), &out, 344.0))
}

// ------------------------------------------------------------------------------------------------ flux network
pub fn flux_network(s: &Solved) -> Option<Figure> {
    let sm = s.mat("S")?;
    let (lb, ub, cv) = (s.vec("lb")?, s.vec("ub")?, s.vec("c")?);
    let obj = s.out("objective")?;
    let (m, r) = (sm.len(), sm[0].len());
    if m == 0 || m > 6 || r > 8 || lb.len() != r || ub.len() != r || cv.len() != r {
        return None;
    }
    let mx = |i: usize| 80.0 + 360.0 * (i as f64 + 0.5) / m as f64;
    let ym = 120.0;
    let mut out = format!("<text class=\"axp-mt\" x=\"16\" y=\"22\">optimal objective = {}</text>", tick(obj));
    for i in 0..m {
        out.push_str(&format!("<circle class=\"axp-ka\" cx=\"{:.1}\" cy=\"{ym}\" r=\"18\"/><circle class=\"axp-k\" cx=\"{:.1}\" \
                               cy=\"{ym}\" r=\"18\"/><text class=\"axp-bt\" x=\"{:.1}\" y=\"{}\" text-anchor=\"middle\">M{}</text>",
                              mx(i), mx(i), mx(i), ym + 5.0, i + 1));
    }
    let arrow = |out: &mut String, a: (f64, f64), b: (f64, f64), cl: &str| {
        let (dx, dy) = (b.0 - a.0, b.1 - a.1);
        let d = (dx * dx + dy * dy).sqrt().max(1e-9);
        let (ux, uy) = (dx / d, dy / d);
        let h = (b.0 - ux * 8.0, b.1 - uy * 8.0);
        out.push_str(&format!("<path class=\"{cl}\" d=\"M{:.1},{:.1} L{:.1},{:.1}\"/><path class=\"axp-kd\" d=\"M{:.1},{:.1} \
                               L{:.1},{:.1} L{:.1},{:.1} Z\"/>", a.0, a.1, h.0, h.1, b.0, b.1, h.0 - uy * 5.0, h.1 + ux * 5.0,
                              h.0 + uy * 5.0, h.1 - ux * 5.0));
    };
    let mut pair_seen: Vec<(usize, usize)> = vec![];
    let (mut uptake, mut output) = (vec![0usize; m], vec![0usize; m]);
    for j in 0..r {
        let cons: Vec<usize> = (0..m).filter(|i| sm[*i][j] < 0.0).collect();
        let prod: Vec<usize> = (0..m).filter(|i| sm[*i][j] > 0.0).collect();
        let cl = if cv[j] != 0.0 { "axp-c" } else { "axp-k" };
        let label = format!("{}R{} [{}, {}]", if cv[j] != 0.0 { "* " } else { "" }, j + 1, g(lb[j], 4), g(ub[j], 4));
        let (lx, ly);
        match (cons.first(), prod.first()) {
            (Some(&a), Some(&b)) => {
                let rep = pair_seen.iter().filter(|p| **p == (a.min(b), a.max(b))).count() as f64;
                pair_seen.push((a.min(b), a.max(b)));
                let off = 40.0 + 30.0 * rep;
                let (x1, x2) = (mx(a), mx(b));
                out.push_str(&format!("<path class=\"{cl}\" d=\"M{:.1},{:.1} Q{:.1},{:.1} {:.1},{:.1}\" fill=\"none\"/>", x1, ym - 18.0,
                                      (x1 + x2) / 2.0, ym - 18.0 - off * 1.6, x2 - 4.0, ym - 20.0));
                arrow(&mut out, (x2 - 10.0, ym - 30.0), (x2 - 3.0, ym - 19.0), cl);
                lx = (x1 + x2) / 2.0;
                ly = ym - 18.0 - off * 0.8 - 4.0;
            }
            (None, Some(&b)) => {
                let k_ = uptake[b] as f64;
                uptake[b] += 1;
                let x = mx(b) - 20.0 - 70.0 * k_;
                arrow(&mut out, (x, ym + 80.0 + 14.0 * k_), (mx(b) - 8.0, ym + 16.0), cl);
                lx = x;
                ly = ym + 96.0 + 14.0 * k_;
            }
            (Some(&a), None) => {
                let k_ = output[a] as f64;
                output[a] += 1;
                let x = mx(a) + 20.0 + 70.0 * k_;
                arrow(&mut out, (mx(a) + 8.0, ym + 16.0), (x, ym + 80.0 + 14.0 * k_), cl);
                lx = x;
                ly = ym + 96.0 + 14.0 * k_;
            }
            (None, None) => continue,
        }
        out.push_str(&format!("<text class=\"axp-vt\" x=\"{lx:.1}\" y=\"{ly:.1}\" text-anchor=\"middle\" font-size=\"11\">{}</text>",
                              esc(&label)));
    }
    out.push_str("<text class=\"axp-lt\" x=\"16\" y=\"250\">* the reaction the objective maximises; [lower, upper] flux bounds</text>");
    Some(s.figure("network", "Flux balance: a small network and its optimal flux",
                  &format!("The metabolic network of the stoichiometric matrix S ({m} metabolites M1 to M{m}, {r} reactions R1 to \
                            R{r}, each arrow running from what a reaction consumes to what it makes; an arrow from or to the \
                            border is an uptake or an output), with the flux bounds [lb, ub]. The starred reaction is the one the \
                            objective c maximises; under steady state the model's flux balance returns objective = {}.",
                           tick(obj)), &out, 262.0))
}

// silence unused imports when a figure family is trimmed
#[allow(dead_code)]
fn _unused(_s: &Solved) -> Option<f64> {
    k::num_of(_s.a, "x")
}
