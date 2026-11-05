//! Enzyme kinetics, fits and growth: measured points with the line or curve of the parameters the model's fit
//! returns (drawn with those returned values), progress curves and steady states from further runs of the model.

use crate::figs_bio::{fi, other, run};
use crate::figs_common::{tick, Chart, DOTS, LINES};
use ax_engine::figkit::{self as k, Solved};
use ax_engine::figures::Figure;

fn pairs(s: &Solved, x: &str, y: &str) -> Option<Vec<(f64, f64)>> {
    let (a, b) = (s.vec(x)?, s.vec(y)?);
    (a.len() == b.len() && a.len() >= 2).then(|| a.into_iter().zip(b).collect())
}

fn list(v: &[f64]) -> String {
    format!("[{}]", v.iter().map(|x| tick(*x)).collect::<Vec<_>>().join(", "))
}

/// The Michaelis-Menten curve v = Vmax S/(Km + S) over S, drawn with returned parameters.
fn mm(vmax: f64, km: f64, xs: &[f64]) -> Vec<(f64, f64)> {
    xs.iter().map(|x| (*x, vmax * x / (km + x))).collect()
}

pub fn mm_fit(s: &Solved) -> Option<Figure> {
    let p = pairs(s, "S", "v")?;
    let (vmax, km) = (s.out("Vmax")?, s.out("Km")?);
    if !(km > 0.0 && vmax > 0.0) {
        return None;
    }
    let smax = p.iter().map(|q| q.0).fold(km, f64::max) * 1.2;
    let xs = k::span(0.0, smax, 81, false);
    let cur = mm(vmax, km, &xs);
    let ys: Vec<f64> = p.iter().map(|q| q.1).chain([0.0, vmax]).collect();
    let mut c = Chart::new((0.0, smax), Chart::range(&ys), (false, false), (false, true), "substrate S", "rate v",
                           &format!("Vmax = {}, Km = {}", tick(vmax), tick(km)));
    c.hguide(vmax, "Vmax");
    c.line(&cur, "axp-c", "", false);
    c.dots(&p, "axp-d2", 5.0);
    c.mark((km, vmax / 2.0), &format!("S = Km = {}", tick(km)));
    let body = c.body();
    Some(s.chart("fit", "Michaelis-Menten kinetics: measured rates and fitted curve",
                 &format!("Michaelis-Menten fit to {} measured rates (S = {}): the model returns Vmax = {} and Km = {}; the \
                           curve drawn is v = Vmax S/(Km + S) with those two values. Marked: S = Km, v = Vmax/2 = {}.",
                          p.len(), list(&s.vec("S")?), tick(vmax), tick(km), tick(vmax / 2.0)), &body))
}

pub fn mm_direct(s: &Solved) -> Option<Figure> {
    let p = pairs(s, "S", "v")?;
    let (vmax, km) = (s.out("Vmax")?, s.out("Km")?);
    if p.iter().any(|q| !(q.0 > 0.0)) {
        return None;
    }
    let smax = p.iter().map(|q| q.0).fold(0.0, f64::max);
    let (x0, x1) = (-smax * 1.1, km.max(smax * 0.2) * 2.0);
    let mut c = Chart::new((x0, x1), (0.0, vmax * 2.0), (false, false), (false, false), "Km", "Vmax",
                           &format!("lines meet at Km = {}, Vmax = {}", tick(km), tick(vmax)));
    for (i, (sv, v)) in p.iter().enumerate() {
        // the line through (-S, 0) and (0, v): Vmax = v + (v/S) Km
        let l = vec![(-sv, 0.0), (x1, v + v / sv * x1)];
        c.line(&l, LINES[i % 3], "", false);
        c.dots(&[(-sv, 0.0), (0.0, *v)], DOTS[i % 3], 3.5);
    }
    c.vguide(0.0, "");
    c.mark((km, vmax), &format!("Km = {}, Vmax = {}", tick(km), tick(vmax)));
    let body = c.body();
    Some(s.chart("fit-2", "Direct linear plot: every measurement as a line in (Km, Vmax)",
                 &format!("Direct linear plot of the {} measured rates (S = {}, v = {}): each measurement is the line through \
                           (-S, 0) and (0, v), the (Km, Vmax) pairs that would explain it exactly. The lines meet near the \
                           parameters the model's fit returns. Marked: Km = {}, Vmax = {}.", p.len(), list(&s.vec("S")?),
                          list(&s.vec("v")?), tick(km), tick(vmax)), &body))
}

/// Points and the straight line y = a + b x drawn with returned parameters; the intercepts marked.
#[allow(clippy::too_many_arguments)]
fn linear(s: &Solved, kind: &str, title: &str, pts: &[(f64, f64)], a: f64, b: f64, xlab: &str, ylab: &str, header: &str,
          marks: &[((f64, f64), String)], caption: &str) -> Option<Figure> {
    let mut xs: Vec<f64> = pts.iter().map(|p| p.0).chain(marks.iter().map(|m| m.0 .0)).collect();
    xs.push(0.0);
    let (xl, xh) = Chart::range(&xs);
    let line = vec![(xl, a + b * xl), (xh, a + b * xh)];
    let ys: Vec<f64> = pts.iter().map(|p| p.1).chain(line.iter().map(|p| p.1)).chain(marks.iter().map(|m| m.0 .1)).collect();
    let mut c = Chart::new((xl, xh), Chart::range(&ys), (false, false), (true, true), xlab, ylab, header);
    c.vguide(0.0, "");
    c.hguide(0.0, "");
    c.line(&line, "axp-c", "", false);
    c.dots(pts, "axp-d2", 5.0);
    // the first intercept is the marked point; the others a dot with the label below it
    for (i, (p, l)) in marks.iter().enumerate() {
        if i == 0 {
            c.mark(*p, l);
        } else if c.xa.inside(p.0) && c.ya.inside(p.1) {
            c.dots(&[*p], "axp-d3", 5.0);
            let (x, y) = (c.xa.px(p.0), c.ya.px(p.1));
            let anchor = if x > 300.0 { "end" } else { "start" };
            c.text(x + if anchor == "end" { -8.0 } else { 8.0 }, (y + 20.0).min(crate::figs_common::CY1 - 4.0), anchor,
                   "axp-mt", l);
        }
    }
    let body = c.body();
    Some(s.chart(kind, title, caption, &body))
}

pub fn lineweaver(s: &Solved) -> Option<Figure> {
    let p = pairs(s, "S", "v")?;
    let (vmax, km) = (s.out("Vmax")?, s.out("Km")?);
    let pts: Vec<(f64, f64)> = p.iter().filter(|q| q.0 != 0.0 && q.1 != 0.0).map(|q| (1.0 / q.0, 1.0 / q.1)).collect();
    linear(s, "linear", "Lineweaver-Burk plot with its fitted line", &pts, 1.0 / vmax, km / vmax, "1/S", "1/v",
           &format!("Vmax = {}, Km = {}", tick(vmax), tick(km)),
           &[((0.0, 1.0 / vmax), format!("1/Vmax = {}", tick(1.0 / vmax))), ((-1.0 / km, 0.0), format!("-1/Km = {}", tick(-1.0 / km)))],
           &format!("Lineweaver-Burk (double-reciprocal) plot of S = {}, v = {}: points at (1/S, 1/v) and the line 1/v = \
                     (Km/Vmax)(1/S) + 1/Vmax drawn with the model's fit Vmax = {}, Km = {} (r2 = {}). Marked: the intercepts \
                     1/Vmax and -1/Km, derived from the returned values.", list(&s.vec("S")?), list(&s.vec("v")?), tick(vmax),
                    tick(km), tick(s.out("r2").unwrap_or(f64::NAN))))
}

pub fn eadie(s: &Solved) -> Option<Figure> {
    let p = pairs(s, "S", "v")?;
    let (vmax, km) = (s.out("Vmax")?, s.out("Km")?);
    let pts: Vec<(f64, f64)> = p.iter().filter(|q| q.0 != 0.0).map(|q| (q.1 / q.0, q.1)).collect();
    linear(s, "linear", "Eadie-Hofstee plot with its fitted line", &pts, vmax, -km, "v/S", "v",
           &format!("Vmax = {}, Km = {}", tick(vmax), tick(km)),
           &[((0.0, vmax), format!("Vmax = {}", tick(vmax))), ((vmax / km, 0.0), format!("Vmax/Km = {}", tick(vmax / km)))],
           &format!("Eadie-Hofstee plot of S = {}, v = {}: points at (v/S, v) and the line v = Vmax - Km (v/S) with the \
                     model's fit Vmax = {}, Km = {}. Marked: the intercepts Vmax and Vmax/Km (derived).", list(&s.vec("S")?),
                    list(&s.vec("v")?), tick(vmax), tick(km)))
}

pub fn hanes(s: &Solved) -> Option<Figure> {
    let p = pairs(s, "S", "v")?;
    let (vmax, km) = (s.out("Vmax")?, s.out("Km")?);
    let pts: Vec<(f64, f64)> = p.iter().filter(|q| q.1 != 0.0).map(|q| (q.0, q.0 / q.1)).collect();
    linear(s, "linear", "Hanes-Woolf plot with its fitted line", &pts, km / vmax, 1.0 / vmax, "S", "S/v",
           &format!("Vmax = {}, Km = {}", tick(vmax), tick(km)),
           &[((0.0, km / vmax), format!("Km/Vmax = {}", tick(km / vmax))), ((-km, 0.0), format!("-Km = {}", tick(-km)))],
           &format!("Hanes-Woolf plot of S = {}, v = {}: points at (S, S/v) and the line S/v = S/Vmax + Km/Vmax with the \
                     model's fit Vmax = {}, Km = {}. Marked: the intercepts -Km and Km/Vmax, derived from the returned values.",
                    list(&s.vec("S")?), list(&s.vec("v")?), tick(vmax), tick(km)))
}

pub fn hill_plot(s: &Solved) -> Option<Figure> {
    let p = pairs(s, "S", "Y")?;
    let (n, kk) = (s.out("n")?, s.out("K")?);
    if !(kk > 0.0) {
        return None;
    }
    let pts: Vec<(f64, f64)> = p.iter().filter(|q| q.0 > 0.0 && q.1 > 0.0 && q.1 < 1.0)
        .map(|q| (q.0.log10(), (q.1 / (1.0 - q.1)).log10())).collect();
    let lk = kk.log10();
    linear(s, "linear", "Hill plot: cooperativity from the slope", &pts, -n * lk, n, "log10 S", "log10(Y/(1 - Y))",
           &format!("n = {}, K = {}", tick(n), tick(kk)), &[((lk, 0.0), format!("S = K = {}", tick(kk)))],
           &format!("Hill plot of the fractional saturations Y = {} at S = {}: points at (log S, log(Y/(1 - Y))) and the line \
                     of slope n through log K, with the model's fit n = {}, K = {}; a slope above 1 means cooperative \
                     binding. Marked: half saturation at S = K.", list(&s.vec("Y")?), list(&s.vec("S")?), tick(n), tick(kk)))
}

pub fn hill_fit(s: &Solved) -> Option<Figure> {
    let p = pairs(s, "S", "v")?;
    let (vmax, kk, n) = (s.out("Vmax")?, s.out("K")?, s.out("n")?);
    if p.iter().any(|q| !(q.0 > 0.0)) || !(kk > 0.0) {
        return None;
    }
    let (lo, hi) = Chart::range_pos(&p.iter().map(|q| q.0).chain([kk]).collect::<Vec<_>>());
    let xs = k::span(lo / 2.0, hi * 2.0, 81, true);
    let h = |m: f64| -> Vec<(f64, f64)> { xs.iter().map(|x| (*x, vmax * x.powf(m) / (kk.powf(m) + x.powf(m)))).collect() };
    let ys: Vec<f64> = p.iter().map(|q| q.1).chain([0.0, vmax]).collect();
    let mut c = Chart::new((lo / 2.0, hi * 2.0), Chart::range(&ys), (true, false), (false, true), "S (logarithmic)", "rate v",
                           &format!("Vmax = {}, K = {}, n = {}", tick(vmax), tick(kk), tick(n)));
    c.line(&h(n), "axp-c", &format!("fit, n = {}", tick(n)), false);
    c.line(&h(1.0), "axp-c4", "n = 1", true);
    c.dots(&p, "axp-d2", 5.0);
    c.mark((kk, vmax / 2.0), &format!("S = K = {}", tick(kk)));
    let body = c.body();
    Some(s.chart("fit", "Hill fit: a cooperative, sigmoid response",
                 &format!("Hill-equation fit to v = {} at S = {} (logarithmic S): the curve v = Vmax S^n/(K^n + S^n) drawn \
                           with the model's fit Vmax = {}, K = {}, n = {}, beside the non-cooperative curve with the same Vmax \
                           and K (n = 1, dashed). Marked: S = K, v = Vmax/2.", list(&s.vec("v")?), list(&s.vec("S")?),
                          tick(vmax), tick(kk), tick(n)), &body))
}

pub fn competitive_fit(s: &Solved) -> Option<Figure> {
    let (sv, iv, vv) = (s.vec("S")?, s.vec("I")?, s.vec("v")?);
    if sv.len() != iv.len() || sv.len() != vv.len() {
        return None;
    }
    let (vmax, km, ki) = (s.out("Vmax")?, s.out("Km")?, s.out("Ki")?);
    let mut levels: Vec<f64> = vec![];
    for i in &iv {
        if !levels.iter().any(|l| (l - i).abs() < 1e-12) {
            levels.push(*i);
        }
    }
    if levels.len() > 4 {
        return None;
    }
    let smax = sv.iter().cloned().fold(km, f64::max) * 1.2;
    let xs = k::span(0.0, smax, 81, false);
    let ys: Vec<f64> = vv.iter().copied().chain([0.0, vmax]).collect();
    let mut c = Chart::new((0.0, smax), Chart::range(&ys), (false, false), (false, true), "substrate S", "rate v",
                           &format!("Vmax = {}, Km = {}, Ki = {}", tick(vmax), tick(km), tick(ki)));
    c.hguide(vmax, "Vmax");
    for (k_, l) in levels.iter().enumerate() {
        let kapp = km * (1.0 + l / ki);
        c.line(&mm(vmax, kapp, &xs), LINES[k_], &format!("I = {}", tick(*l)), false);
        let pts: Vec<(f64, f64)> = sv.iter().zip(&iv).zip(&vv).filter(|((_, i), _)| (*i - l).abs() < 1e-12)
            .map(|((s_, _), v)| (*s_, *v)).collect();
        c.dots(&pts, DOTS[k_], 5.0);
    }
    let body = c.body();
    Some(s.chart("fit", "Competitive inhibition: rates with and without inhibitor",
                 &format!("Rates measured at S = {} with inhibitor I = {} (dots), and the curves v = Vmax S/(Km(1 + I/Ki) + S) \
                           drawn with the model's fit Vmax = {}, Km = {}, Ki = {}. A competitive inhibitor raises the apparent \
                           Km and leaves Vmax unchanged.", list(&sv), list(&iv), tick(vmax), tick(km), tick(ki)), &body))
}

pub fn competitive_linear(s: &Solved) -> Option<Figure> {
    let (sv, iv, vv) = (s.vec("S")?, s.vec("I")?, s.vec("v")?);
    if sv.len() != iv.len() || sv.len() != vv.len() || sv.iter().chain(&vv).any(|x| *x == 0.0) {
        return None;
    }
    let (vmax, km, ki) = (s.out("Vmax")?, s.out("Km")?, s.out("Ki")?);
    let mut levels: Vec<f64> = vec![];
    for i in &iv {
        if !levels.iter().any(|l| (l - i).abs() < 1e-12) {
            levels.push(*i);
        }
    }
    if levels.len() > 4 {
        return None;
    }
    let xmax = sv.iter().map(|x| 1.0 / x).fold(0.0, f64::max) * 1.1;
    let mut ys: Vec<f64> = vv.iter().map(|v| 1.0 / v).collect();
    ys.push(0.0);
    for l in &levels {
        ys.push(1.0 / vmax + km * (1.0 + l / ki) / vmax * xmax);
    }
    let mut c = Chart::new((0.0, xmax), Chart::range(&ys), (false, false), (false, true), "1/S", "1/v",
                           &format!("common intercept 1/Vmax = {}", tick(1.0 / vmax)));
    for (k_, l) in levels.iter().enumerate() {
        let slope = km * (1.0 + l / ki) / vmax;
        c.line(&[(0.0, 1.0 / vmax), (xmax, 1.0 / vmax + slope * xmax)], LINES[k_], &format!("I = {}", tick(*l)), false);
        let pts: Vec<(f64, f64)> = sv.iter().zip(&iv).zip(&vv).filter(|((_, i), _)| (*i - l).abs() < 1e-12)
            .map(|((s_, _), v)| (1.0 / s_, 1.0 / v)).collect();
        c.dots(&pts, DOTS[k_], 5.0);
    }
    c.mark((0.0, 1.0 / vmax), &format!("1/Vmax = {}", tick(1.0 / vmax)));
    let body = c.body();
    Some(s.chart("linear", "Competitive inhibition in a double-reciprocal plot",
                 &format!("The same measurements in double-reciprocal form (points at 1/S, 1/v) with the lines of the model's \
                           fit Vmax = {}, Km = {}, Ki = {}: a competitive inhibitor turns the line about the common intercept \
                           1/Vmax (marked, derived).", tick(vmax), tick(km), tick(ki)), &body))
}

pub fn schumacher(s: &Solved) -> Option<Figure> {
    let (cv, rv, mv) = (s.vec("C")?, s.vec("R")?, s.vec("M")?);
    if cv.len() != rv.len() || cv.len() != mv.len() || cv.iter().any(|c| *c == 0.0) {
        return None;
    }
    let n = s.out("N")?;
    if !(n > 0.0) {
        return None;
    }
    let pts: Vec<(f64, f64)> = mv.iter().zip(rv.iter().zip(&cv)).map(|(m, (r, c))| (*m, r / c)).collect();
    let xmax = mv.iter().cloned().fold(0.0, f64::max).max(1.0) * 1.15;
    let ys: Vec<f64> = pts.iter().map(|p| p.1).chain([0.0, xmax / n]).collect();
    let mut c = Chart::new((0.0, xmax), Chart::range(&ys), (false, false), (false, true), "marked animals at large M",
                           "share of the catch marked R/C", &format!("N = {}", tick(n)));
    c.line(&[(0.0, 0.0), (xmax, xmax / n)], "axp-c", "", false);
    c.dots(&pts, "axp-d2", 5.0);
    let body = c.body();
    Some(s.chart("fit", "Schumacher-Eschmeyer estimate of a population size",
                 &format!("Schumacher-Eschmeyer mark-recapture: on each occasion the share of the catch already marked, R/C \
                           (from C = {}, R = {}), against the marked animals at large M = {}. In a closed population R/C is \
                           about M/N; the line through the origin has slope 1/N with the model's estimate N = {}.",
                          list(&cv), list(&rv), list(&mv), tick(n)), &body))
}

pub fn species_area(s: &Solved) -> Option<Figure> {
    let p = pairs(s, "A", "S")?;
    let (cc, z) = (s.out("c")?, s.out("z")?);
    if p.iter().any(|q| !(q.0 > 0.0 && q.1 > 0.0)) {
        return None;
    }
    let (lo, hi) = Chart::range_pos(&p.iter().map(|q| q.0).chain([1.0]).collect::<Vec<_>>());
    let xs = k::span(lo, hi, 41, true);
    let cur: Vec<(f64, f64)> = xs.iter().map(|a| (*a, cc * a.powf(z))).collect();
    let ys: Vec<f64> = p.iter().map(|q| q.1).chain(cur.iter().map(|q| q.1)).collect();
    let mut c = Chart::new((lo, hi), Chart::range_pos(&ys), (true, true), (true, true), "area A (logarithmic)",
                           "species S (logarithmic)", &format!("S = {} A^{}", tick(cc), tick(z)));
    c.line(&cur, "axp-c", "", false);
    c.dots(&p, "axp-d2", 5.0);
    c.mark((1.0, cc), &format!("c = {}", tick(cc)));
    let body = c.body();
    Some(s.chart("fit", "Species-area relation and its power-law fit",
                 &format!("Species-area relation on logarithmic axes: the counts S = {} on areas A = {} (dots) and the power law \
                           S = c A^z drawn with the model's log-log fit c = {}, z = {} (r2 = {}). Marked: c, the species \
                           expected on a unit area.", list(&s.vec("S")?), list(&s.vec("A")?), tick(cc), tick(z),
                          tick(s.out("r2").unwrap_or(f64::NAN))), &body))
}

// ------------------------------------------------------------------------------------------------ progress curves
pub fn mm_progress(s: &Solved) -> Option<Figure> {
    let (s0, vmax, km) = (s.inp("S0")?, s.inp("Vmax")?, s.inp("Km")?);
    let (t0, sv) = (s.inp("t")?, s.out("S")?);
    if !(s0 > 0.0 && vmax > 0.0 && km > 0.0) {
        return None;
    }
    let tmax = ((s0 * 0.98 + km * 50f64.ln()) / vmax).max(t0 * 1.1);
    let ts = k::span(0.0, tmax, 101, false);
    let mut sp = vec![];
    let mut pp = vec![];
    for t in &ts {
        if let Some(a) = s.with(&[("t", *t)]) {
            if let (Some(x), Some(y)) = (k::num_of(&a, "S"), k::num_of(&a, "P")) {
                sp.push((*t, x));
                pp.push((*t, y));
            }
        }
    }
    if sp.len() < 2 {
        return None;
    }
    let th = other(s, "bio_michaelis_menten_time", &[("S0", fi(s0)), ("Vmax", fi(vmax)), ("Km", fi(km)), ("X", fi(0.5))], "t");
    let mut c = Chart::new((0.0, tmax), (0.0, s0), (false, false), (false, true), "time t", "concentration",
                           &format!("S0 = {}, Vmax = {}, Km = {}", tick(s0), tick(vmax), tick(km)));
    c.line(&sp, "axp-c", "substrate S", false);
    c.line(&pp, "axp-c2", "product P", false);
    if let Some(t) = th {
        c.vguide(t, &format!("half converted at t = {}", tick(t)));
    }
    c.mark((t0, sv), &format!("t = {}, S = {}", tick(t0), tick(sv)));
    let body = c.body();
    Some(s.chart("kinetics", "Enzyme progress curve: substrate and product over time",
                 &format!("Substrate and product over time for an enzyme reaction with S0 = {}, Vmax = {} and Km = {} \
                           (integrated Michaelis-Menten law), computed by the model at {} times{}. Marked: t = {}, S = {}.",
                          tick(s0), tick(vmax), tick(km), sp.len(),
                          th.map(|t| format!("; half the substrate is converted at t = {} (the model's time capability)", tick(t)))
                              .unwrap_or_default(), tick(t0), tick(sv)), &body))
}

pub fn mm_progress_family(s: &Solved) -> Option<Figure> {
    let (vmax, km) = (s.inp("Vmax")?, s.inp("Km")?);
    if !(vmax > 0.0 && km > 0.0) {
        return None;
    }
    let starts = [km / 5.0, km * 2.0, km * 20.0];
    let tmax = (starts[2] * 0.9 + km * 10f64.ln()) / vmax;
    let ts = k::span(0.0, tmax, 60, false);
    let mut c = Chart::new((0.0, tmax), (0.0, 1.0), (false, false), (false, true), "time t", "substrate left S/S0",
                           &format!("Vmax = {}, Km = {}", tick(vmax), tick(km)));
    let mut n = 0;
    for (k_, s0) in starts.iter().enumerate() {
        let pts: Vec<(f64, f64)> = ts.iter().filter_map(|t| run(s, &[("S0", fi(*s0)), ("t", fi(*t))], "S").map(|v| (*t, v / s0)))
            .collect();
        if pts.len() >= 2 {
            c.line(&pts, LINES[k_], &format!("S0 = {}", tick(*s0)), false);
            n += 1;
        }
    }
    if n < 2 {
        return None;
    }
    let body = c.body();
    Some(s.chart("kinetics-2", "Progress curves from below to far above Km",
                 &format!("Substrate left as a fraction of the start, S/S0, against time for an enzyme with Vmax = {} and Km = \
                           {} starting from S0 = {}, {} and {} (a fifth, twice and twenty times Km): S is the model's integrated \
                           Michaelis-Menten value at 60 times, divided by S0. Far below Km the decay is exponential; far above \
                           it the enzyme is saturated and S falls at a nearly constant rate.", tick(vmax), tick(km),
                          tick(starts[0]), tick(starts[1]), tick(starts[2])), &body))
}

pub fn chemostat(s: &Solved) -> Option<Figure> {
    let (mu, d0) = (s.inp("mu_max")?, s.inp("D")?);
    let (sv, xv) = (s.out("S")?, s.out("X")?);
    if !(mu > 0.0) {
        return None;
    }
    let ds = k::span(mu * 0.005, mu, 90, false);
    let mut sp = vec![];
    let mut xp = vec![];
    let mut pr = vec![];
    for d in &ds {
        if let Some(a) = s.with(&[("D", *d)]) {
            if let (Some(a1), Some(a2)) = (k::num_of(&a, "S"), k::num_of(&a, "X")) {
                sp.push((*d, a1));
                xp.push((*d, a2));
                pr.push((*d, d * a2));
            }
        }
    }
    if sp.len() < 2 {
        return None;
    }
    let best = pr.iter().cloned().fold(pr[0], |a, b| if b.1 > a.1 { b } else { a });
    let ys: Vec<f64> = sp.iter().chain(&xp).chain(&pr).map(|p| p.1).chain([0.0]).collect();
    let mut c = Chart::new((0.0, mu), Chart::range(&ys), (false, false), (false, true), "dilution rate D",
                           "steady-state value", &format!("D = {}: X = {}, S = {}", tick(d0), tick(xv), tick(sv)));
    c.line(&xp, "axp-c", "biomass X", false);
    c.line(&sp, "axp-c2", "substrate S", false);
    c.line(&pr, "axp-c3", "productivity D X", true);
    c.mark((d0, xv), &format!("D = {}", tick(d0)));
    let body = c.body();
    Some(s.chart("operating", "Chemostat steady states against dilution rate",
                 &format!("Chemostat operating diagram: steady-state substrate and biomass computed by the model at {} dilution \
                           rates up to washout (rates the model refuses are left out); biomass productivity D X (derived) \
                           peaks near D = {} ({}). Marked: D = {}, X = {}, S = {}.", sp.len(), tick(best.0), tick(best.1),
                          tick(d0), tick(xv), tick(sv)), &body))
}
