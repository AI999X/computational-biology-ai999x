"""Computational biology procedures: deterministic capabilities with typed inputs and several outputs (domain code,
loaded by the pack; enabled only once tools/certify_procedures.py has certified them).

Topics (the biology view): DNA, RNA and protein sequences (reverse complement, transcription, translation with the
standard genetic code, GC content, k-mers, p-distance and edit distance, open reading frames, restriction sites,
molecular masses, extinction coefficients, net charge and isoelectric point, hydropathy, instability index, melting
temperatures), pairwise alignment scores (global, local and semi-global; linear and affine gaps; BLOSUM62),
Mendelian and population genetics (crosses, Hardy-Weinberg tests, selection, mutation-selection balance, drift,
fixation, effective population size, Fst, linkage disequilibrium, map functions, pedigree inbreeding), phylogenetics
(Jukes-Cantor and Kimura distances, UPGMA, neighbour joining, tree counts, Robinson-Foulds distance, Fitch parsimony,
coalescent times, nucleotide diversity, Watterson's theta, Tajima's D), enzyme kinetics (Michaelis-Menten fits and
linearisations, progress curves, Hill and competitive-inhibition fits) and population biology and ecology (Leslie
matrices, diversity and richness indices, rarefaction, beta diversity, species-area fits, mark-recapture estimators,
chemostat steady state, flux balance analysis, Faith's phylogenetic diversity).

Each procedure states its inputs, outputs, conditions, method, LaTeX, the documented question form it parses, worked
examples, an INDEPENDENT reference computation (used only by certification; the checking libraries are imported inside
the reference functions) and a deliberately wrong variant that the checks must reject. The shipped code (run, verify,
parse) is plain Python (math, re, fractions, itertools).

How the references differ from the shipped code:
- sequences: Biopython's own implementations (Seq.reverse_complement/transcribe/translate, gc_fraction,
  molecular_weight, ProtParam, IsoelectricPoint charge (root found by scipy brentq), MeltingTemp.Tm_Wallace and Tm_NN,
  Bio.Restriction search/catalyse), scikit-bio k-mer frequencies, scipy Hamming distance, Biopython's PairwiseAligner
  (unit-cost edit distance), ORFs read from Biopython translations of the six frames;
- alignment: Biopython PairwiseAligner (scores only, never alignments);
- genetics: exhaustive gamete enumeration (Punnett squares) for crosses, scikit-allel allele counts and
  heterozygosities with random union of gametes and scipy chi-square for Hardy-Weinberg, brute-force enumeration of
  allele placements for the exact test, closed forms against iteration for selection, forward iteration of the
  recursion against a root for equilibria, Wright-Fisher and Moran Markov chains (numpy) for drift and fixation, exact
  enumeration of parental gene choices for the sex-ratio effective size, scipy harmonic mean, Wright's variance form
  of Fst, scipy chi-square and a linear programme (scipy linprog) for linkage disequilibrium, a Poisson crossover
  process and quadrature/ODE integration for map functions, the tabular (additive-relationship) method for pedigrees;
- phylogenetics: matrix exponentials of the substitution rate matrices (scipy expm) solved for the distance, scipy
  average linkage, scikit-bio neighbour joining and Robinson-Foulds distance, explicit enumeration of trees, Biopython
  ParsimonyScorer (Fitch), absorbing Markov chains for coalescent times, scikit-allel diversity statistics;
- enzyme kinetics: scipy curve_fit and linregress for fits, ODE integration (scipy DOP853) of progress curves,
  quadrature of the integrated rate law;
- ecology: numpy matrix powers and eigen-decompositions for Leslie matrices, scikit-bio alpha-diversity functions,
  scipy hypergeometric probabilities, scipy distances, scipy linregress, brute-force maximisation of the scipy
  hypergeometric likelihood, numpy weighted least squares, ODE integration to steady state, scipy linprog (HiGHS)
  against an exact rational simplex, scikit-bio Faith's PD.

Conventions (fixed, so run and reference agree exactly):
- Sequences are written 5'->3' as letters only (upper or lower case, read as upper case; any other character is
  refused, never skipped). DNA uses A C G T, RNA uses A C G U, proteins use the 20 standard one-letter codes.
- Genetic code: the standard code (NCBI table 1); stop codons TAA, TAG, TGA. Translation reads frame 1 of a sequence
  whose length is a multiple of 3; stops are shown as '*' or the protein is truncated at the first stop, as stated.
- An open reading frame starts at ATG and ends at the first in-frame stop codon (stop included in its length); all
  six frames are searched; an ORF without a stop codon is not counted.
- Restriction sites: linear DNA; a cut position is the 1-based position of the first base after the cut on the top
  strand (Bio.Restriction's convention); recognition sites and cut offsets are REBASE facts.
- Masses (Biopython's tables, Biopython licence, BSD-style): average nucleotide-monophosphate masses with a 5'
  phosphate (DNA A 331.2218, C 307.1971, G 347.2212, T 322.2085; RNA A 347.2212, C 323.1965, G 363.2206,
  U 324.1813), average amino-acid masses (IUPAC), water 18.0153 removed per bond; g/mol.
- Extinction coefficients at 280 nm (Pace et al. / ProtParam): 5500 per Trp, 1490 per Tyr, 125 per cystine (pair of
  Cys); M^-1 cm^-1; absorbance of a 1 g/L solution in a 1 cm cell = coefficient / mass.
- Net charge and isoelectric point: Bjellqvist pK set as Biopython defines it (N-terminus 7.5, K 10.0, R 12.0,
  H 5.98, C-terminus 3.55, D 4.05, E 4.45, C 9.0, Y 10.0; residue-specific terminal values for N-terminal A, M, S, P,
  T, V, E and C-terminal D, E); Henderson-Hasselbalch partial charges.
- Hydropathy: Kyte-Doolittle scale; instability index: Guruprasad et al. 1990 dipeptide table (DIWV, as in Biopython).
- Melting temperature: Wallace rule 2(A+T) + 4(G+C) degC; nearest-neighbour Tm with the Allawi-SantaLucia 1997 table
  (Biopython DNA_NN3), salt correction 0.368 (N-1) ln[Na+] on the entropy (SantaLucia 1998), k = c1 - c2/2 (strands
  in nM), a self-complementary sequence (equal to its reverse complement) with the symmetry correction and k = c1;
  R = 1.987 cal/(mol K); perfect complement, no mismatches or dangling ends.
- Alignment scores: optimal score only. A gap of length L scores open + (L-1)*extend (linear gaps: open = extend =
  gap). Global alignments score end gaps like internal ones; semi-global alignments leave end gaps unscored; local
  alignments (Smith-Waterman) are never below 0. BLOSUM62 is the NCBI matrix (public domain).
- Genotypes are written locus by locus as letter pairs ('AaBb'): upper case dominant, lower case recessive, complete
  dominance, independent assortment, the same loci in the same order in every genotype.
- Population genetics: discrete, non-overlapping generations; random mating; fitness relative values; Wright-Fisher
  populations of N diploids hold 2N gene copies; heterozygosity is the expected heterozygosity 2pq of the gene pool.
- Map distances in centimorgans; recombination fractions below 0.5; Haldane or Kosambi map as stated.
- Phylogenetics: aligned sequences of equal length without gaps; UPGMA heights are half the cluster distance; neighbour
  joining keeps negative branch lengths as computed; ties in a clustering step are refused (the tree is not unique).
- Diversity: counts per species; zero counts are not species of the sample; the logarithm base is stated.
- Results built from exp, log or roots are rounded to 12 significant digits (identical on every platform and Python
  version); values that are zero by construction are written as exactly 0.

Question form (one per procedure): the capability phrase, a colon, then `name = value` pairs separated by commas
(`Reverse complement: seq = ATGCCGTA`). Numbers are pure numbers (lists in brackets, matrices as [[..], [..]]);
sequences, genotypes and trees are written as they are; choices are words. Any other number, an unknown name, a name
given twice or free prose is refused, never guessed. The reading guard (tools/ask_sweep.py, 2026-10-23):
  1. every stated number is bound: the question is split into 'name = value' pairs and any other text or number
     refuses the whole question;
  2. the capability phrase must open the question and be followed by the colon, so wording that names a different
     capability (global vs local, Hamming vs edit distance, UPGMA vs neighbour joining, Shannon vs Simpson, Haldane
     vs Kosambi) never reaches this reader;
  3. a convention is an input (log base, Simpson estimator, Chao1 form, map function, strand, stop handling, enzyme,
     molecule type): required, never assumed;
  4. a sequence is a text input read exactly as written (letters only); it is never read partially.
"""
from __future__ import annotations

import math
import re
from fractions import Fraction
from itertools import combinations, product as _iproduct


# ================================================================================================ shared helpers
def _r(x):
    """12 significant digits, no negative zero (identical results on every platform and Python version)."""
    if isinstance(x, bool) or isinstance(x, str) or x is None:
        return x
    if isinstance(x, int):
        return x
    if isinstance(x, dict):
        return {k: _r(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [_r(v) for v in x]
    if x == 0 or not math.isfinite(x):
        return 0.0 if x == 0 else x
    return float(f"{x:.12g}") + 0.0


def _out(**kw):
    return {k: _r(v) for k, v in kw.items()}


def _rel(a, b):
    if isinstance(a, dict) or isinstance(b, dict):
        if not (isinstance(a, dict) and isinstance(b, dict)) or a.keys() != b.keys():
            return math.inf
        return max([_rel(a[k], b[k]) for k in a] + [0.0])
    if isinstance(a, (list, tuple)) or isinstance(b, (list, tuple)):
        if not (isinstance(a, (list, tuple)) and isinstance(b, (list, tuple))) or len(a) != len(b):
            return math.inf
        return max([_rel(x, y) for x, y in zip(a, b)] + [0.0])
    if isinstance(a, (str, bool)) or isinstance(b, (str, bool)) or a is None or b is None:
        return 0.0 if a == b else math.inf
    try:
        a, b = float(a), float(b)
    except (TypeError, ValueError):
        return math.inf
    if not (math.isfinite(a) and math.isfinite(b)):
        return 0.0 if a == b else math.inf
    d = abs(a - b)
    if d == 0:
        return 0.0
    if abs(a) < 1e-15 and abs(b) < 1e-15:
        return 0.0
    return d / max(abs(a), abs(b))


def _recheck(fn):
    """A verifier: every output the procedure gives is recomputed and compared (relative); a missing output fails."""
    def verify(inputs, outputs):
        try:
            want = fn(**inputs)
        except Exception:
            return math.inf
        worst = 0.0
        for k, v in want.items():
            if k not in outputs:
                return math.inf
            worst = max(worst, _rel(outputs[k], v))
        return worst
    return verify


def _safe(fn):
    """Guard a wrong variant against exceptions (a variant that raises cannot be compared)."""
    def g(**kw):
        try:
            return fn(**kw)
        except Exception:
            return {"_failed": True}
    return g


def _need(ok, msg):
    if not ok:
        raise ValueError(msg)


def _fin(name, *vals):
    for v in vals:
        _need(isinstance(v, (int, float)) and not isinstance(v, bool) and math.isfinite(v), f"{name} must be a finite number")


def _pos(name, *vals):
    for v in vals:
        _need(isinstance(v, (int, float)) and not isinstance(v, bool) and math.isfinite(v) and v > 0,
              f"{name} must be a positive number")


def _frac(name, v, lo_open=False, hi_open=False):
    _fin(name, v)
    _need((v > 0 if lo_open else v >= 0) and (v < 1 if hi_open else v <= 1), f"{name} must be a fraction between 0 and 1")


def _whole(name, v, lo=0):
    _need(isinstance(v, (int, float)) and not isinstance(v, bool) and math.isfinite(v) and float(v) == int(v) and v >= lo,
          f"{name} must be a whole number of at least {lo}")
    return int(v)


def _counts(name, v, lo=0, positive=False):
    _need(isinstance(v, (list, tuple)) and len(v) >= 1, f"{name} must be a list of counts")
    out = [_whole(name, x, lo) for x in v]
    if positive:
        _need(all(x > 0 for x in out), f"{name}: every count must be positive")
    return out


def _bisect(f, a, b, n=400):
    """Root of f on [a, b] (f(a), f(b) of opposite signs) by bisection, to full double precision."""
    fa, fb = f(a), f(b)
    if fa == 0:
        return a
    if fb == 0:
        return b
    _need((fa > 0) != (fb > 0), "no root in the interval")
    for _ in range(n):
        m = 0.5 * (a + b)
        if m == a or m == b:
            break
        fm = f(m)
        if fm == 0:
            return m
        if (fm > 0) == (fa > 0):
            a, fa = m, fm
        else:
            b = m
    return 0.5 * (a + b)


def _brentq(f, a, b):
    from scipy.optimize import brentq
    return float(brentq(f, a, b, xtol=1e-300, rtol=8.9e-16, maxiter=1000))


def _ivp(fun, t_end, y0, rtol=1e-13):
    from scipy.integrate import solve_ivp
    atol = [1e-16 * max(1.0, abs(v)) for v in y0]
    sol = solve_ivp(fun, (0.0, t_end), list(y0), method="DOP853", rtol=rtol, atol=atol)
    if not sol.success:
        raise RuntimeError(sol.message)
    return [float(v) for v in sol.y[:, -1]]


def _quad(f, a, b):
    from scipy.integrate import quad
    v, _ = quad(f, a, b, epsabs=0.0, epsrel=1e-13, limit=400)
    return float(v)


def _lu(rng, a, b):
    return math.exp(rng.uniform(math.log(a), math.log(b)))


def _sig(x, n=4):
    return float(f"{x:.{n}g}")


# ================================================================================================ the question reader
_NUM = r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?"
_PAT = {
    "seq": r"[A-Za-z]+",
    "seqs": r"[A-Za-z]+(?:\s+[A-Za-z]+)+",
    "names": r"[A-Za-z][A-Za-z0-9_]*(?:\s+[A-Za-z][A-Za-z0-9_]*)*",
    "name": r"[A-Za-z][A-Za-z0-9_]*",
    "newick": r"\([A-Za-z0-9_.:,()\s+\-eE]*\)\s*;?",
    "states": r"[A-Za-z][A-Za-z0-9_]*:[A-Za-z]+(?:\s+[A-Za-z][A-Za-z0-9_]*:[A-Za-z]+)+",
}


def _split(rest):
    """Split on commas outside brackets and parentheses."""
    parts, depth, cur = [], 0, ""
    for ch in rest:
        if ch in "[(":
            depth += 1
        elif ch in "])":
            depth -= 1
        if ch == "," and depth == 0:
            parts.append(cur)
            cur = ""
        else:
            cur += ch
    parts.append(cur)
    return parts


def _value(spec, text):
    """One written value for an input of this type, or None when it cannot be read exactly."""
    t = spec.get("type", "number")
    text = text.strip()
    if t == "number":
        if spec.get("pct") and re.fullmatch(_NUM + r"\s*%", text):
            return float(text.rstrip("% ")) / 100.0
        return float(text) if re.fullmatch(_NUM, text) else None
    if t == "integer":
        return int(text) if re.fullmatch(r"[-+]?\d+", text) else None
    if t == "vector":
        m = re.fullmatch(r"\[([^\[\]]*)\]", text)
        if not m:
            return None
        items = [p for p in re.split(r"[,\s;]+", m.group(1).strip()) if p]
        if not items or not all(re.fullmatch(_NUM, p) for p in items):
            return None
        return ", ".join(items)
    if t == "matrix":
        m = re.fullmatch(r"\[\s*(\[.*\])\s*\]", text, re.S)
        if not m:
            return None
        rows = re.findall(r"\[([^\[\]]*)\]", m.group(1))
        if not rows or re.sub(r"\[[^\[\]]*\]", "", m.group(1)).strip(" ,;\n\t"):
            return None
        out = []
        for row in rows:
            items = [p for p in re.split(r"[,\s]+", row.strip()) if p]
            if not items or not all(re.fullmatch(_NUM, p) for p in items):
                return None
            out.append(" ".join(items))
        return "; ".join(out)
    if t == "edges":
        m = re.fullmatch(r"\[([^\[\]]*)\]", text)
        return m.group(1).strip() if m else None
    if t == "text":
        return text if re.fullmatch(_PAT[spec.get("pattern", "seq")], text) else None
    if t == "choice":
        return text.lower() if re.fullmatch(r"[A-Za-z0-9_\-]+", text) else None
    if t == "boolean":
        return text.lower() if text.lower() in ("true", "false", "yes", "no") else None
    return None


def _reader(phrases, inputs):
    """The one documented question form, '<phrase>: name = value, ...', for a procedure with these inputs."""
    rx = re.compile(r"\s*(?:" + "|".join(re.escape(p) for p in phrases) + r")\s*:\s*", re.I)
    by = {i["symbol"]: i for i in inputs}
    names = {s: s for s in by}
    for s in by:                                  # a lower-case spelling only where it names no other input
        names.setdefault(s.lower(), s)

    def parse(question):
        if not isinstance(question, str):
            return None
        m = rx.match(question)
        if not m:
            return None
        rest = re.sub(r"[.?!]\s*$", "", question[m.end():].strip()).strip()
        if not rest:
            return {}
        out = {}
        for seg in _split(rest):
            mm = re.fullmatch(r"\s*(?:and\s+)?([A-Za-z][A-Za-z_0-9]*)\s*=\s*(.+?)\s*", seg, re.S)
            if not mm:
                return None
            key = names.get(mm.group(1)) or names.get(mm.group(1).lower())
            if key is None or key in out:
                return None
            v = _value(by[key], mm.group(2))
            if v is None:
                return None
            out[key] = v
        return out
    return parse


def _written(spec, v):
    """An input value as a user writes it in the documented form."""
    t = spec.get("type", "number")
    if isinstance(v, str):
        return "[" + v + "]" if t in ("edges", "vector", "matrix") and not v.startswith("[") else v
    if t == "vector":
        return "[" + ", ".join(_num_text(x) for x in v) + "]"
    if t == "matrix":
        return "[" + ", ".join("[" + ", ".join(_num_text(x) for x in row) + "]" for row in v) + "]"
    if t == "edges":
        return "[" + ", ".join("-".join(str(x) for x in e[:2]) for e in v) + "]"
    if t in ("choice", "text", "boolean"):
        return str(v)
    if t == "integer":
        return str(int(v))
    return _num_text(v)


def _num_text(x):
    x = float(x)
    return str(int(x)) if x == int(x) and abs(x) < 1e15 else repr(x)


# ================================================================================================ input/output specs
def _numin(sym, meaning, required=True, pct=False):
    d = {"symbol": sym, "kind": "number", "meaning": meaning, "required": required}
    if pct:
        d["pct"] = True
    return d


def _intin(sym, meaning, lo=None, hi=None, required=True):
    d = {"symbol": sym, "type": "integer", "meaning": meaning, "required": required}
    if lo is not None:
        d["min"] = lo
    if hi is not None:
        d["max"] = hi
    return d


def _vecin(sym, meaning, lo=1, hi=200, required=True):
    return {"symbol": sym, "type": "vector", "kind": "number", "meaning": meaning, "required": required,
            "min_length": lo, "max_length": hi}


def _matin(sym, meaning, required=True):
    return {"symbol": sym, "type": "matrix", "kind": "number", "meaning": meaning, "required": required}


def _choice(sym, choices, meaning, required=True):
    return {"symbol": sym, "type": "choice", "choices": list(choices), "meaning": meaning, "required": required}


def _textin(sym, meaning, pattern="seq", required=True, hi=100000):
    return {"symbol": sym, "type": "text", "meaning": meaning, "required": required, "pattern": pattern,
            "min_length": 1, "max_length": hi}


def _edgesin(sym, meaning, required=True):
    return {"symbol": sym, "type": "edges", "meaning": meaning, "required": required, "min_length": 1, "max_length": 200}


def _o(sym, meaning, kind="number"):
    return {"symbol": sym, "kind": kind, "meaning": meaning}


PROCEDURES: dict = {}
LEVELS: dict = {}


def _add(key, *, name, topic, purpose, inputs, outputs, primary, method, latex, run, reference, wrong, rand,
         phrases, examples, level, conditions=(), extra_q=(), verify=None):
    """Register a procedure: the documented question form is built from its phrases and inputs; templates and
    questions are written from the worked examples exactly as a user writes them."""
    parse = _reader(phrases, inputs)
    by = {i["symbol"]: i for i in inputs}

    def ask_text(ex_inputs):
        return f"{phrases[0]}: " + ", ".join(f"{k} = {_written(by[k], v)}" for k, v in ex_inputs.items())
    templates = [ask_text(e["inputs"]) for e in examples[:2]]
    questions = [{"q": ask_text(e["inputs"]), "outputs": e["outputs"]} for e in examples]
    req = [i["symbol"] for i in inputs if i.get("required")]
    if req and examples:
        partial = {k: v for k, v in examples[0]["inputs"].items() if k != req[0]}
        questions.append({"q": ask_text(partial), "outcome": "NEEDS_INPUT"})
        questions.append({"q": ask_text(examples[0]["inputs"]) + ", zz = 1", "outcome": "UNRECOGNIZED_FORMAT"})
        questions.append({"q": ask_text(examples[0]["inputs"]) + " and 7 more", "outcome": "UNRECOGNIZED_FORMAT"})
    questions += list(extra_q)
    PROCEDURES[key] = {
        "name": name, "topic": topic, "purpose": purpose, "inputs": inputs, "outputs": outputs, "primary": primary,
        "conditions": list(conditions), "method": method, "latex": latex, "run": run,
        "verify": verify or _recheck(run), "parse": parse, "reference": reference, "wrong": _safe(wrong),
        "random_inputs": rand, "templates": templates, "examples": list(examples), "questions": questions}
    LEVELS[key] = level


def _ex(inputs, outputs):
    return {"inputs": inputs, "outputs": outputs}


# ================================================================================================ sequence helpers
_BASES = "TCAG"
_AA_STD = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
CODON = {a + b + c: _AA_STD[16 * i + 4 * j + k]
         for i, a in enumerate(_BASES) for j, b in enumerate(_BASES) for k, c in enumerate(_BASES)}
STOPS = ("TAA", "TAG", "TGA")
_COMP = str.maketrans("ACGTU", "TGCAA")
_COMP_RNA = str.maketrans("ACGU", "UGCA")
AMINO = "ACDEFGHIKLMNPQRSTVWY"


def _seq(s, name, kind):
    """A sequence read exactly as written: letters only, upper-cased; kind is 'DNA', 'RNA', 'NUC' (DNA or RNA),
    'PROTEIN' or 'ANY' (any letters)."""
    _need(isinstance(s, str) and s != "", f"{name} must be a sequence of letters")
    u = s.upper()
    _need(len(u) <= 100000, f"{name}: at most 100000 letters")
    if kind == "DNA":
        _need(re.fullmatch(r"[ACGT]+", u) is not None, f"{name} must be a DNA sequence of A, C, G and T only")
    elif kind == "RNA":
        _need(re.fullmatch(r"[ACGU]+", u) is not None, f"{name} must be an RNA sequence of A, C, G and U only")
    elif kind == "NUC":
        _need(re.fullmatch(r"[ACGT]+|[ACGU]+", u) is not None,
              f"{name} must be a DNA (A, C, G, T) or RNA (A, C, G, U) sequence, not a mixture")
    elif kind == "PROTEIN":
        _need(re.fullmatch("[" + AMINO + "]+", u) is not None,
              f"{name} must be a protein sequence of the 20 standard one-letter amino-acid codes")
    else:
        _need(re.fullmatch(r"[A-Z]+", u) is not None, f"{name} must be letters only")
    return u


def _revcomp(s):
    return s.translate(_COMP)[::-1]


def _rand_dna(rng, lo, hi, alphabet="ACGT"):
    return "".join(rng.choice(alphabet) for _ in range(rng.randint(lo, hi)))


def _rand_protein(rng, lo, hi):
    return "".join(rng.choice(AMINO) for _ in range(rng.randint(lo, hi)))


# ================================================================================================ SEQUENCES
def reverse_complement(seq):
    s = _seq(seq, "seq", "DNA")
    return {"revcomp": _revcomp(s)}


def _ref_revcomp(seq):
    from Bio.Seq import Seq
    return {"revcomp": str(Seq(seq.upper()).reverse_complement())}


_add("bio_reverse_complement", name="Reverse complement of a DNA sequence", topic="sequences / DNA",
     purpose="The reverse complement of a DNA strand (the partner strand read 5'->3').",
     inputs=[_textin("seq", "DNA sequence 5'->3' (A, C, G, T)")],
     outputs=[_o("revcomp", "reverse complement 5'->3'", "text")], primary="revcomp",
     conditions=["DNA letters A, C, G, T only"], method="complement each base (A-T, C-G) and reverse the order",
     latex=r"5'\text{-}x_1\dots x_n\text{-}3' \mapsto 5'\text{-}\bar x_n\dots \bar x_1\text{-}3'",
     run=reverse_complement, reference=_ref_revcomp,
     wrong=lambda seq: {"revcomp": seq.upper().translate(_COMP)},
     rand=lambda rng: {"seq": _rand_dna(rng, 6, 60)},
     phrases=["Reverse complement", "Reverse complement of DNA"],
     examples=[_ex({"seq": "ATGCCGTA"}, {"revcomp": "TACGGCAT"}),
               _ex({"seq": "GAATTC"}, {"revcomp": "GAATTC"}),
               _ex({"seq": "aaccgt"}, {"revcomp": "ACGGTT"})],
     extra_q=[{"q": "Reverse complement: seq = ATG1CG", "outcome": "UNRECOGNIZED_FORMAT"},
              {"q": "Reverse complement: seq = ATGXCG", "outcome": "INVALID_INPUT"},
              {"q": "Reverse complement: seq = AUGC", "outcome": "INVALID_INPUT"}],
     level="beginner")


def transcription(seq, strand):
    s = _seq(seq, "seq", "DNA")
    _need(strand in ("coding", "template"), "strand must be coding or template")
    m = s if strand == "coding" else _revcomp(s)
    return {"mrna": m.replace("T", "U")}


def _ref_transcription(seq, strand):
    from Bio.Seq import Seq
    s = Seq(seq.upper())
    return {"mrna": str(s.transcribe() if strand == "coding" else s.reverse_complement().transcribe())}


def _wrong_transcription(seq, strand):
    s = seq.upper()
    m = _revcomp(s) if strand == "coding" else s.translate(_COMP)      # the strands mixed up
    return {"mrna": m.replace("T", "U")}


_add("bio_transcription", name="Transcription of DNA into mRNA", topic="sequences / gene expression",
     purpose="The mRNA transcribed from a DNA strand given 5'->3', stated as the coding (sense) or template "
             "(antisense) strand.",
     inputs=[_textin("seq", "DNA strand 5'->3' (A, C, G, T)"),
             _choice("strand", ["coding", "template"], "which strand is given: coding (sense) or template (antisense)")],
     outputs=[_o("mrna", "mRNA 5'->3'", "text")], primary="mrna",
     conditions=["DNA letters only", "the strand is stated"],
     method="coding strand: T -> U; template strand: reverse complement, then T -> U",
     latex=r"\text{mRNA} = \text{coding strand with } T\to U = \overline{\text{template}}^{\,R}",
     run=transcription, reference=_ref_transcription, wrong=_wrong_transcription,
     rand=lambda rng: {"seq": _rand_dna(rng, 6, 60), "strand": rng.choice(["coding", "template"])},
     phrases=["Transcription", "Transcribe DNA"],
     examples=[_ex({"seq": "ATGGCTTAA", "strand": "coding"}, {"mrna": "AUGGCUUAA"}),
               _ex({"seq": "TTAAGCCAT", "strand": "template"}, {"mrna": "AUGGCUUAA"})],
     extra_q=[{"q": "Transcription: seq = ATGGCTTAA", "outcome": "NEEDS_INPUT"}],
     level="beginner")


_MITO = dict(CODON, AGA="*", AGG="*", ATA="M", TGA="W")


def _translate(s, table, truncate):
    p = "".join(table[s[i:i + 3]] for i in range(0, len(s), 3))
    if truncate and "*" in p:
        p = p[:p.index("*")]
    return p


def translation(seq, stops):
    s = _seq(seq, "seq", "NUC").replace("U", "T")
    _need(len(s) % 3 == 0, "the sequence length must be a multiple of 3 (whole codons read in frame 1)")
    _need(stops in ("show", "truncate"), "stops must be show or truncate")
    p = _translate(s, CODON, stops == "truncate")
    return {"protein": p, "length": len(p)}


def _ref_translation(seq, stops):
    from Bio.Seq import Seq
    p = str(Seq(seq.upper()).translate(table=1, to_stop=(stops == "truncate")))
    return {"protein": p, "length": len(p)}


def _wrong_translation(seq, stops):
    p = _translate(seq.upper().replace("U", "T"), _MITO, stops == "truncate")     # the mitochondrial code
    return {"protein": p, "length": len(p)}


def _rand_translation(rng):
    n = rng.randint(20, 40)
    codons = [rng.choice([c for c in CODON if c not in STOPS]) for _ in range(n)]
    for _ in range(rng.randint(1, 4)):                                   # codons that read differently elsewhere
        codons[rng.randrange(n)] = rng.choice(["AGA", "AGG", "ATA", "TGA"])
    s = "".join(codons)
    if rng.random() < 0.3:
        s = s.replace("T", "U")
    return {"seq": s, "stops": rng.choice(["show", "truncate"])}


_add("bio_translation", name="Translation with the standard genetic code", topic="sequences / gene expression",
     purpose="The protein encoded by a DNA or mRNA sequence read in frame 1 with the standard genetic code; stop codons "
             "shown as '*' or the protein truncated at the first stop.",
     inputs=[_textin("seq", "DNA or mRNA sequence 5'->3', length a multiple of 3"),
             _choice("stops", ["show", "truncate"], "show stop codons as '*' or truncate at the first stop")],
     outputs=[_o("protein", "amino-acid sequence (one-letter codes)", "text"), _o("length", "number of residues shown")],
     primary="protein", conditions=["length a multiple of 3", "DNA or RNA letters, not a mixture"],
     method="read codons in frame 1 and look each up in the standard genetic code (NCBI table 1)",
     latex=r"\text{codon } x_{3i+1}x_{3i+2}x_{3i+3} \mapsto \text{amino acid}_i",
     run=translation, reference=_ref_translation, wrong=_wrong_translation, rand=_rand_translation,
     phrases=["Translation", "Translate"],
     examples=[_ex({"seq": "ATGGCCATTGTAATGGGCCGCTGA", "stops": "show"}, {"protein": "MAIVMGR*", "length": 8}),
               _ex({"seq": "AUGUUUUGGUAAGGC", "stops": "truncate"}, {"protein": "MFW", "length": 3}),
               _ex({"seq": "ATGAGATGA", "stops": "show"}, {"protein": "MR*", "length": 3})],
     extra_q=[{"q": "Translation: seq = ATGGC, stops = show", "outcome": "INVALID_INPUT"}],
     level="beginner")


def gc_content(seq):
    s = _seq(seq, "seq", "NUC")
    gc = s.count("G") + s.count("C")
    return {"gc_fraction": gc / len(s), "gc_count": gc}


def _ref_gc(seq):
    from Bio.Seq import Seq
    from Bio.SeqUtils import gc_fraction
    s = Seq(seq.upper())
    return {"gc_fraction": float(gc_fraction(s)), "gc_count": s.count("G") + s.count("C")}


_add("bio_gc_content", name="GC content of a nucleic acid", topic="sequences / composition",
     purpose="The fraction and number of G and C bases in a DNA or RNA sequence.",
     inputs=[_textin("seq", "DNA or RNA sequence")],
     outputs=[_o("gc_fraction", "fraction of G + C"), _o("gc_count", "number of G and C bases")],
     primary="gc_fraction", conditions=["unambiguous DNA or RNA letters"],
     method="count G and C and divide by the length", latex=r"\mathrm{GC} = \frac{n_G + n_C}{n}",
     run=gc_content, reference=_ref_gc,
     wrong=lambda seq: {"gc_fraction": seq.upper().count("G") / len(seq), "gc_count": seq.upper().count("G")},
     rand=lambda rng: {"seq": _rand_dna(rng, 5, 80, rng.choice(["ACGT", "ACGU"]))},
     phrases=["GC content"],
     examples=[_ex({"seq": "ATGCGCGA"}, {"gc_fraction": 0.625, "gc_count": 5}),
               _ex({"seq": "GGAUCUUCGGAUCU"}, {"gc_fraction": 0.5, "gc_count": 7})],
     level="beginner")


def kmer_counts(seq, k):
    s = _seq(seq, "seq", "NUC")
    _need(isinstance(k, int) and 1 <= k <= len(s), "k must be a whole number between 1 and the sequence length")
    counts = {}
    for i in range(len(s) - k + 1):
        w = s[i:i + k]
        counts[w] = counts.get(w, 0) + 1
    top = max(counts.values())
    return {"total": len(s) - k + 1, "distinct": len(counts), "max_count": top,
            "most_frequent": min(w for w, c in counts.items() if c == top)}


def _ref_kmer(seq, k):
    from skbio import Sequence
    f = Sequence(seq.upper()).kmer_frequencies(k, overlap=True)
    top = max(f.values())
    return {"total": int(sum(f.values())), "distinct": len(f), "max_count": int(top),
            "most_frequent": min(w for w, c in f.items() if c == top)}


def _wrong_kmer(seq, k):
    s = seq.upper()
    counts = {}
    for i in range(0, len(s) - k + 1, k):                     # non-overlapping windows
        counts[s[i:i + k]] = counts.get(s[i:i + k], 0) + 1
    top = max(counts.values())
    return {"total": sum(counts.values()), "distinct": len(counts), "max_count": top,
            "most_frequent": min(w for w, c in counts.items() if c == top)}


_add("bio_kmer_counts", name="k-mer counts of a sequence", topic="sequences / composition",
     purpose="Overlapping k-mer counts: the number of k-mers, of distinct k-mers, the highest count and the most "
             "frequent k-mer (alphabetically first among ties).",
     inputs=[_textin("seq", "DNA or RNA sequence"), _intin("k", "k-mer length", lo=1)],
     outputs=[_o("total", "number of overlapping k-mers (n - k + 1)"), _o("distinct", "number of distinct k-mers"),
              _o("max_count", "count of the most frequent k-mer"),
              _o("most_frequent", "most frequent k-mer (alphabetically first among ties)", "text")],
     primary="distinct", conditions=["1 <= k <= sequence length"],
     method="slide a window of length k one base at a time and count each word",
     latex=r"N_k = n - k + 1", run=kmer_counts, reference=_ref_kmer, wrong=_wrong_kmer,
     rand=lambda rng: {"seq": _rand_dna(rng, 12, 80, "ACGT" if rng.random() < 0.8 else "AC"), "k": rng.randint(2, 4)},
     phrases=["k-mer counts", "Kmer counts"],
     examples=[_ex({"seq": "ATATAT", "k": 2}, {"total": 5, "distinct": 2, "max_count": 3, "most_frequent": "AT"}),
               _ex({"seq": "GATTACA", "k": 3}, {"total": 5, "distinct": 5, "max_count": 1, "most_frequent": "ACA"})],
     level="intermediate")


def hamming(seq1, seq2):
    a, b = _seq(seq1, "seq1", "ANY"), _seq(seq2, "seq2", "ANY")
    _need(len(a) == len(b), "the two aligned sequences must have the same length")
    d = sum(1 for x, y in zip(a, b) if x != y)
    return {"mismatches": d, "p_distance": d / len(a), "identity": (len(a) - d) / len(a)}


def _ref_hamming(seq1, seq2):
    from scipy.spatial.distance import hamming as h
    p = float(h(list(seq1.upper()), list(seq2.upper())))
    return {"mismatches": int(round(p * len(seq1))), "p_distance": p, "identity": 1.0 - p}


def _wrong_hamming(seq1, seq2):
    m = sum(1 for x, y in zip(seq1.upper(), seq2.upper()) if x == y)           # matches counted
    return {"mismatches": m, "p_distance": m / len(seq1), "identity": 1 - m / len(seq1)}


def _rand_pair_same_length(rng, alphabet="ACGT"):
    a = _rand_dna(rng, 8, 60, alphabet)
    b = "".join(c if rng.random() < 0.75 else rng.choice(alphabet) for c in a)
    return a, b


_add("bio_hamming_distance", name="p-distance (Hamming) between aligned sequences", topic="sequences / comparison",
     purpose="Mismatches, p-distance and identity of two aligned sequences of equal length (no indels): the "
             "proportion of sites that differ.",
     inputs=[_textin("seq1", "first aligned sequence"), _textin("seq2", "second aligned sequence (same length)")],
     outputs=[_o("mismatches", "number of differing sites"), _o("p_distance", "proportion of differing sites"),
              _o("identity", "proportion of identical sites")],
     primary="mismatches", conditions=["equal lengths (aligned, no gaps)"],
     method="compare the sequences site by site", latex=r"p = \frac{1}{n}\sum_i [x_i \ne y_i]",
     run=hamming, reference=_ref_hamming, wrong=_wrong_hamming,
     rand=lambda rng: dict(zip(("seq1", "seq2"), _rand_pair_same_length(rng, rng.choice(["ACGT", AMINO])))),
     phrases=["p-distance", "Hamming distance of aligned sequences", "Hamming distance"],
     examples=[_ex({"seq1": "GAGCCTACTAACGGGAT", "seq2": "CATCGTAATGACGGCCT"},
                   {"mismatches": 7, "p_distance": 7 / 17}),
               _ex({"seq1": "ACGT", "seq2": "ACGT"}, {"mismatches": 0, "p_distance": 0.0, "identity": 1.0})],
     extra_q=[{"q": "Hamming distance: seq1 = ACGT, seq2 = ACG", "outcome": "INVALID_INPUT"}],
     level="beginner")


def _levenshtein(a, b):
    prev = list(range(len(b) + 1))
    for i in range(1, len(a) + 1):
        cur = [i] + [0] * len(b)
        for j in range(1, len(b) + 1):
            cur[j] = min(prev[j] + 1, cur[j - 1] + 1, prev[j - 1] + (a[i - 1] != b[j - 1]))
        prev = cur
    return prev[-1]


def edit_distance(seq1, seq2):
    a, b = _seq(seq1, "seq1", "ANY"), _seq(seq2, "seq2", "ANY")
    _need(len(a) * len(b) <= 4_000_000, "sequences too long (length product at most 4 000 000)")
    return {"distance": _levenshtein(a, b)}


def _ref_edit(seq1, seq2):
    from Bio.Align import PairwiseAligner
    al = PairwiseAligner()
    al.mode = "global"
    al.match_score, al.mismatch_score = 0.0, -1.0
    al.open_gap_score = al.extend_gap_score = -1.0
    return {"distance": int(round(-al.score(seq1.upper(), seq2.upper())))}


def _wrong_edit(seq1, seq2):
    a, b = seq1.upper(), seq2.upper()                      # insertions and deletions only (no substitutions)
    prev = [0] * (len(b) + 1)
    for i in range(1, len(a) + 1):
        cur = [0] * (len(b) + 1)
        for j in range(1, len(b) + 1):
            cur[j] = prev[j - 1] + 1 if a[i - 1] == b[j - 1] else max(prev[j], cur[j - 1])
        prev = cur
    return {"distance": len(a) + len(b) - 2 * prev[-1]}


def _mutate(rng, s, alphabet, p_sub=0.15, p_indel=0.08):
    out = []
    for c in s:
        r = rng.random()
        if r < p_indel / 2:
            continue                                           # deletion
        out.append(c if rng.random() >= p_sub else rng.choice(alphabet))
        if rng.random() < p_indel / 2:
            out.append(rng.choice(alphabet))                   # insertion
    return "".join(out) or rng.choice(alphabet)


def _rand_edit(rng):
    alpha = rng.choice(["ACGT", AMINO])
    a = _rand_dna(rng, 6, 40, alpha)
    return {"seq1": a, "seq2": _mutate(rng, a, alpha, 0.3, 0.15)}


_add("bio_edit_distance", name="Edit (Levenshtein) distance between sequences", topic="sequences / comparison",
     purpose="The minimum number of substitutions, insertions and deletions turning one sequence into the other "
             "(unit costs; an unaligned comparison, unlike the p-distance).",
     inputs=[_textin("seq1", "first sequence"), _textin("seq2", "second sequence")],
     outputs=[_o("distance", "edit distance (unit costs)")], primary="distance",
     conditions=["letters only"], method="dynamic programming over prefixes (Wagner-Fischer)",
     latex=r"D_{i,j} = \min(D_{i-1,j}+1,\ D_{i,j-1}+1,\ D_{i-1,j-1}+[x_i\ne y_j])",
     run=edit_distance, reference=_ref_edit, wrong=_wrong_edit, rand=_rand_edit,
     phrases=["Edit distance", "Levenshtein distance"],
     examples=[_ex({"seq1": "KITTEN", "seq2": "SITTING"}, {"distance": 3}),
               _ex({"seq1": "GATTACA", "seq2": "GCATGCT"}, {"distance": 4})],
     level="intermediate")


def _orf_lengths(s):
    res = []
    for f in range(3):
        start = None
        for i in range(f, len(s) - 2, 3):
            c = s[i:i + 3]
            if start is None:
                if c == "ATG":
                    start = i
            elif c in STOPS:
                res.append(i + 3 - start)
                start = None
    return res


def longest_orf(seq):
    s = _seq(seq, "seq", "DNA")
    res = _orf_lengths(s) + _orf_lengths(_revcomp(s))
    L = max(res) if res else 0
    return {"length_nt": L, "protein_length": L // 3 - 1 if L else 0, "orf_count": len(res)}


def _ref_orf(seq):
    from Bio.Seq import Seq
    res = []
    for strand in (Seq(seq.upper()), Seq(seq.upper()).reverse_complement()):
        for f in range(3):
            sub = strand[f:]
            sub = sub[:len(sub) - len(sub) % 3]
            p = str(sub.translate(table=1))
            prev = -1
            for j, a in enumerate(p):
                if a == "*":
                    seg = p[prev + 1:j]
                    if "M" in seg:
                        res.append((j - (prev + 1 + seg.index("M")) + 1) * 3)
                    prev = j
    L = max(res) if res else 0
    return {"length_nt": L, "protein_length": L // 3 - 1 if L else 0, "orf_count": len(res)}


def _wrong_orf(seq):
    o = longest_orf(seq)
    L = max(o["length_nt"] - 3, 0)                                 # the stop codon left out of the length
    return {"length_nt": L, "protein_length": L // 3 - 1 if L else 0, "orf_count": o["orf_count"]}


def _rand_orf(rng):
    body = "".join(rng.choice([c for c in CODON if c not in STOPS]) for _ in range(rng.randint(5, 40)))
    orf = "ATG" + body + rng.choice(STOPS)
    if rng.random() < 0.4:
        orf = _revcomp(orf)
    return {"seq": _rand_dna(rng, 0, 60) + orf + _rand_dna(rng, 0, 60)}


_add("bio_longest_orf", name="Longest open reading frame", topic="sequences / gene finding",
     purpose="The longest open reading frame (ATG to the first in-frame stop, stop included) over all six frames of "
             "a DNA sequence, its protein length and the number of ORFs.",
     inputs=[_textin("seq", "DNA sequence 5'->3'")],
     outputs=[_o("length_nt", "length of the longest ORF in nucleotides, stop codon included"),
              _o("protein_length", "amino acids encoded by the longest ORF (stop excluded)"),
              _o("orf_count", "number of ORFs (one per stop codon reached from an upstream in-frame ATG)")],
     primary="length_nt", conditions=["DNA letters only", "standard genetic code, start ATG, stops TAA TAG TGA"],
     method="scan the three frames of both strands codon by codon: open at the first ATG, close at the next stop",
     latex=r"\ell = 3(j - i + 1)\ \text{for ATG at codon } i,\ \text{first stop at codon } j",
     run=longest_orf, reference=_ref_orf, wrong=_wrong_orf, rand=_rand_orf,
     phrases=["Longest ORF", "Longest open reading frame"],
     examples=[_ex({"seq": "CCATGAAATTTGGGTAACC"}, {"length_nt": 15, "protein_length": 4, "orf_count": 1}),
               _ex({"seq": "TTACATCCCAT"}, {"length_nt": 6, "protein_length": 1, "orf_count": 1}),
               _ex({"seq": "ACGTACGT"}, {"length_nt": 0, "protein_length": 0, "orf_count": 0})],
     level="intermediate")


# recognition site and cut offset on the top strand (bases before the cut), REBASE; all palindromic sites
ENZYMES = {"EcoRI": ("GAATTC", 1), "BamHI": ("GGATCC", 1), "HindIII": ("AAGCTT", 1), "PstI": ("CTGCAG", 5),
           "SmaI": ("CCCGGG", 3), "KpnI": ("GGTACC", 5), "SacI": ("GAGCTC", 5), "XbaI": ("TCTAGA", 1),
           "SalI": ("GTCGAC", 1), "XhoI": ("CTCGAG", 1), "NotI": ("GCGGCCGC", 2), "EcoRV": ("GATATC", 3),
           "NcoI": ("CCATGG", 1), "NdeI": ("CATATG", 2), "HaeIII": ("GGCC", 2), "AluI": ("AGCT", 2),
           "TaqI": ("TCGA", 1), "MspI": ("CCGG", 1), "ClaI": ("ATCGAT", 2), "BglII": ("AGATCT", 1),
           "SpeI": ("ACTAGT", 1), "NheI": ("GCTAGC", 1), "ApaI": ("GGGCCC", 5), "SphI": ("GCATGC", 5),
           "MluI": ("ACGCGT", 1), "ScaI": ("AGTACT", 3), "PvuII": ("CAGCTG", 3)}


def restriction_sites(seq, enzyme):
    s = _seq(seq, "seq", "DNA")
    _need(enzyme in ENZYMES, "enzyme must be one of " + ", ".join(ENZYMES))
    site, cut = ENZYMES[enzyme]
    pos = [i + cut + 1 for i in range(len(s) - len(site) + 1) if s.startswith(site, i)]
    bounds = [1] + pos + [len(s) + 1]
    return {"n_sites": len(pos), "positions": pos, "fragments": [bounds[i + 1] - bounds[i] for i in range(len(bounds) - 1)]}


def _ref_restriction(seq, enzyme):
    import Bio.Restriction as R
    from Bio.Seq import Seq
    e = getattr(R, enzyme)
    pos = sorted(int(p) for p in e.search(Seq(seq.upper()), linear=True))
    frags = [len(x) for x in e.catalyse(Seq(seq.upper()), linear=True)]
    return {"n_sites": len(pos), "positions": pos, "fragments": frags}


def _wrong_restriction(seq, enzyme):
    s = seq.upper()
    site, cut = ENZYMES[enzyme]
    pos = [i + 1 for i in range(len(s) - len(site) + 1) if s.startswith(site, i)]   # the site start, not the cut
    bounds = [1] + pos + [len(s) + 1]
    return {"n_sites": len(pos), "positions": pos, "fragments": [bounds[i + 1] - bounds[i] for i in range(len(bounds) - 1)]}


def _rand_restriction(rng):
    enz = rng.choice(sorted(ENZYMES))
    parts = [_rand_dna(rng, 2, 30)]
    for _ in range(rng.randint(1, 3)):
        parts += [ENZYMES[enz][0], _rand_dna(rng, 2, 30)]
    return {"seq": "".join(parts), "enzyme": enz}


_add("bio_restriction_sites", name="Restriction sites and fragments", topic="sequences / molecular cloning",
     purpose="Where a restriction enzyme cuts a linear DNA sequence: the number of sites, the cut positions and the "
             "fragment lengths.",
     inputs=[_textin("seq", "linear DNA sequence 5'->3'"),
             _choice("enzyme", sorted(ENZYMES), "restriction enzyme (name)")],
     outputs=[_o("n_sites", "number of cuts"),
              _o("positions", "cut positions (1-based first base after each cut, top strand)"),
              _o("fragments", "fragment lengths in bp, 5'->3' order")],
     primary="n_sites", conditions=["linear DNA", "a listed enzyme"],
     method="find every occurrence of the recognition site (overlaps included) and cut at the enzyme's offset",
     latex=r"\text{cut at } i + c \text{ for each site start } i", run=restriction_sites, reference=_ref_restriction,
     wrong=_wrong_restriction, rand=_rand_restriction,
     phrases=["Restriction sites", "Restriction digest"],
     examples=[_ex({"seq": "AAGAATTCTTGAATTCGGATCC", "enzyme": "EcoRI"},
                   {"n_sites": 2, "positions": [4, 12], "fragments": [3, 8, 11]}),
               _ex({"seq": "AAACTGCAGTTT", "enzyme": "PstI"}, {"n_sites": 1, "positions": [9], "fragments": [8, 4]})],
     extra_q=[{"q": "Restriction sites: seq = AAGAATTC, enzyme = FooI", "outcome": "INVALID_INPUT"}],
     level="intermediate")


DNA_W = {"A": 331.2218, "C": 307.1971, "G": 347.2212, "T": 322.2085}
RNA_W = {"A": 347.2212, "C": 323.1965, "G": 363.2206, "U": 324.1813}
PROTEIN_W = {"A": 89.0932, "C": 121.1582, "D": 133.1027, "E": 147.1293, "F": 165.1891, "G": 75.0666, "H": 155.1546,
             "I": 131.1729, "K": 146.1876, "L": 131.1729, "M": 149.2113, "N": 132.1179, "P": 115.1305, "Q": 146.1445,
             "R": 174.201, "S": 105.0926, "T": 119.1192, "V": 117.1463, "W": 204.2252, "Y": 181.1885}
WATER = 18.0153


def _chain_mass(s, table):
    return sum(table[x] for x in s) - (len(s) - 1) * WATER


def nucleic_mass(seq, molecule, strands):
    _need(molecule in ("DNA", "RNA"), "molecule must be DNA or RNA")
    _need(strands in ("single", "double"), "strands must be single or double")
    s = _seq(seq, "seq", molecule)
    table = DNA_W if molecule == "DNA" else RNA_W
    w = _chain_mass(s, table)
    if strands == "double":
        w += _chain_mass(s.translate(_COMP if molecule == "DNA" else _COMP_RNA), table)
    return _out(mass=w)


def _ref_nucleic_mass(seq, molecule, strands):
    from Bio.SeqUtils import molecular_weight
    return {"mass": float(molecular_weight(seq.upper(), seq_type=molecule, double_stranded=(strands == "double")))}


def _wrong_nucleic_mass(seq, molecule, strands):
    table = DNA_W if molecule == "DNA" else RNA_W
    w = sum(table[x] for x in seq.upper())                                   # no water removed per bond
    return {"mass": 2 * w if strands == "double" else w}


def _rand_nucleic(rng):
    mol = rng.choice(["DNA", "RNA"])
    return {"seq": _rand_dna(rng, 2, 60, "ACGT" if mol == "DNA" else "ACGU"), "molecule": mol,
            "strands": rng.choice(["single", "double"])}


_add("bio_nucleic_acid_mass", name="Molecular mass of DNA or RNA", topic="sequences / molecular mass",
     purpose="The average molecular mass of a single- or double-stranded linear DNA or RNA with a 5' phosphate.",
     inputs=[_textin("seq", "DNA or RNA sequence"), _choice("molecule", ["DNA", "RNA"], "molecule type"),
             _choice("strands", ["single", "double"], "single- or double-stranded")],
     outputs=[_o("mass", "molecular mass in g/mol (Da)")], primary="mass",
     conditions=["the sequence letters match the stated molecule"],
     method="sum of nucleotide-monophosphate masses minus one water per phosphodiester bond (each strand)",
     latex=r"M = \sum_i m(x_i) - (n-1)\,m_{\mathrm{H_2O}}",
     run=nucleic_mass, reference=_ref_nucleic_mass, wrong=_wrong_nucleic_mass, rand=_rand_nucleic,
     phrases=["Nucleic acid mass", "Molecular weight of DNA or RNA"],
     examples=[_ex({"seq": "AGC", "molecule": "DNA", "strands": "single"}, {"mass": 949.6095}),
               _ex({"seq": "AGC", "molecule": "RNA", "strands": "single"}, {"mass": 997.6077}),
               _ex({"seq": "ATGC", "molecule": "DNA", "strands": "double"}, {"mass": 2507.6054})],
     extra_q=[{"q": "Nucleic acid mass: seq = AGU, molecule = DNA, strands = single", "outcome": "INVALID_INPUT"}],
     level="beginner")


def protein_mass(seq):
    s = _seq(seq, "seq", "PROTEIN")
    return _out(mass=_chain_mass(s, PROTEIN_W))


def _ref_protein_mass(seq):
    from Bio.SeqUtils.ProtParam import ProteinAnalysis
    return {"mass": float(ProteinAnalysis(seq.upper()).molecular_weight())}


_PROT_EX = ("MAEGEITTFTALTEKFNLPPGNYKKPKLLYCSNGGHFLRILPDGTVDGTRDRSDQHIQLQLSAESVGEVYIKSTETGQYLAMDTSGLLYGSQTPSEEC"
            "LFLERLEENHYNTYTSKKHAEKNWFVGLKKNGSCKRGPRTHYGQKAILFLPLPV")

_add("bio_protein_mass", name="Molecular mass of a protein", topic="sequences / molecular mass",
     purpose="The average molecular mass of a protein from its amino-acid sequence.",
     inputs=[_textin("seq", "protein sequence (20 standard one-letter codes)")],
     outputs=[_o("mass", "molecular mass in g/mol (Da)")], primary="mass",
     conditions=["standard amino acids only"],
     method="sum of amino-acid masses minus one water per peptide bond",
     latex=r"M = \sum_i m(a_i) - (n-1)\,m_{\mathrm{H_2O}}",
     run=protein_mass, reference=_ref_protein_mass,
     wrong=lambda seq: {"mass": sum(PROTEIN_W[x] for x in seq.upper())},
     rand=lambda rng: {"seq": _rand_protein(rng, 2, 120)},
     phrases=["Protein mass", "Protein molecular weight"],
     examples=[_ex({"seq": _PROT_EX}, {"mass": 17103.1617}), _ex({"seq": "G"}, {"mass": 75.0666})],
     extra_q=[{"q": "Protein mass: seq = MAXB", "outcome": "INVALID_INPUT"}],
     level="beginner")


def extinction(seq):
    s = _seq(seq, "seq", "PROTEIN")
    red = 5500 * s.count("W") + 1490 * s.count("Y")
    ox = red + (s.count("C") // 2) * 125
    m = _chain_mass(s, PROTEIN_W)
    return _out(eps_reduced=red, eps_cystines=ox, abs_reduced=red / m, abs_cystines=ox / m)


def _ref_extinction(seq):
    from Bio.SeqUtils.ProtParam import ProteinAnalysis
    pa = ProteinAnalysis(seq.upper())
    red, ox = pa.molar_extinction_coefficient()
    m = pa.molecular_weight()
    return {"eps_reduced": red, "eps_cystines": ox, "abs_reduced": red / m, "abs_cystines": ox / m}


def _wrong_extinction(seq):
    s = seq.upper()
    red = 1490 * s.count("W") + 5500 * s.count("Y")                         # Trp and Tyr values swapped
    ox = red + (s.count("C") // 2) * 125
    m = _chain_mass(s, PROTEIN_W)
    return {"eps_reduced": red, "eps_cystines": ox, "abs_reduced": red / m, "abs_cystines": ox / m}


def _rand_ext(rng):
    s = list(_rand_protein(rng, 10, 120))
    s[rng.randrange(len(s))] = "W"
    return {"seq": "".join(s)}


_add("bio_protein_extinction", name="Protein extinction coefficient at 280 nm", topic="sequences / protein properties",
     purpose="The molar extinction coefficient at 280 nm (cysteines reduced, and all cysteines as cystines) and the "
             "absorbance of a 1 g/L solution.",
     inputs=[_textin("seq", "protein sequence")],
     outputs=[_o("eps_reduced", "extinction coefficient, cysteines reduced (M^-1 cm^-1)"),
              _o("eps_cystines", "extinction coefficient, cysteines paired as cystines (M^-1 cm^-1)"),
              _o("abs_reduced", "absorbance of 1 g/L, 1 cm, cysteines reduced"),
              _o("abs_cystines", "absorbance of 1 g/L, 1 cm, cystines")],
     primary="eps_reduced", conditions=["standard amino acids only"],
     method="5500 per Trp + 1490 per Tyr (+ 125 per cystine); divide by the molecular mass for 1 g/L",
     latex=r"\varepsilon_{280} = 5500\,n_W + 1490\,n_Y + 125\,n_{\text{cystine}}",
     run=extinction, reference=_ref_extinction, wrong=_wrong_extinction, rand=_rand_ext,
     phrases=["Protein extinction coefficient", "Extinction coefficient"],
     examples=[_ex({"seq": _PROT_EX}, {"eps_reduced": 17420, "eps_cystines": 17545}),
               _ex({"seq": "MWYCC"}, {"eps_reduced": 6990, "eps_cystines": 7115})],
     level="intermediate")


PK_POS = {"Nterm": 7.5, "K": 10.0, "R": 12.0, "H": 5.98}
PK_NEG = {"Cterm": 3.55, "D": 4.05, "E": 4.45, "C": 9.0, "Y": 10.0}
PK_CTERM = {"D": 4.55, "E": 4.75}
PK_NTERM = {"A": 7.59, "M": 7.0, "S": 6.93, "P": 8.36, "T": 6.82, "V": 7.44, "E": 7.7}


def _charge(s, pH, termini=True):
    pos, neg = dict(PK_POS), dict(PK_NEG)
    if s[0] in PK_NTERM:
        pos["Nterm"] = PK_NTERM[s[0]]
    if s[-1] in PK_CTERM:
        neg["Cterm"] = PK_CTERM[s[-1]]
    n = {aa: float(s.count(aa)) for aa in "KRHDECY"}
    n["Nterm"] = n["Cterm"] = 1.0 if termini else 0.0
    plus = 0.0
    for aa, pk in pos.items():
        plus += n[aa] * (1.0 / (10 ** (pH - pk) + 1.0))
    minus = 0.0
    for aa, pk in neg.items():
        minus += n[aa] * (1.0 / (10 ** (pk - pH) + 1.0))
    return plus - minus


def net_charge(seq, pH):
    s = _seq(seq, "seq", "PROTEIN")
    _fin("pH", pH)
    _need(0 <= pH <= 14, "pH must lie between 0 and 14")
    return _out(charge=_charge(s, pH))


def _ref_net_charge(seq, pH):
    from Bio.SeqUtils.ProtParam import ProteinAnalysis
    return {"charge": float(ProteinAnalysis(seq.upper()).charge_at_pH(pH))}


_add("bio_protein_net_charge", name="Net charge of a protein at a pH", topic="sequences / protein properties",
     purpose="The net charge of a protein at a stated pH from its ionisable groups (Bjellqvist pK set).",
     inputs=[_textin("seq", "protein sequence"), _numin("pH", "pH of the solution (0 to 14)")],
     outputs=[_o("charge", "net charge (elementary charges)")], primary="charge",
     conditions=["standard amino acids only", "0 <= pH <= 14"],
     method="Henderson-Hasselbalch partial charges of the termini and the K, R, H, D, E, C, Y side chains",
     latex=r"Z = \sum_{\text{basic}} \frac{n_i}{1+10^{\,pH-pK_i}} - \sum_{\text{acidic}} \frac{n_j}{1+10^{\,pK_j-pH}}",
     run=net_charge, reference=_ref_net_charge,
     wrong=lambda seq, pH: {"charge": _charge(seq.upper(), pH, termini=False)},
     rand=lambda rng: {"seq": _rand_protein(rng, 3, 80), "pH": round(rng.uniform(2, 12), 2)},
     phrases=["Protein net charge", "Net charge of a protein"],
     examples=[_ex({"seq": "INGAR", "pH": 7}, {"charge": 0.760091614289}),
               _ex({"seq": "PETER", "pH": 4.53}, {"charge": 0.00279562756338})],
     level="intermediate")


def isoelectric_point(seq):
    s = _seq(seq, "seq", "PROTEIN")
    pi = _bisect(lambda x: _charge(s, x), 0.0, 14.0)
    return _out(pI=pi, charge_at_7=_charge(s, 7.0))


def _ref_pi(seq):
    from Bio.SeqUtils.ProtParam import ProteinAnalysis
    pa = ProteinAnalysis(seq.upper())
    return {"pI": _brentq(lambda x: pa.charge_at_pH(x), 0.0, 14.0), "charge_at_7": float(pa.charge_at_pH(7.0))}


def _wrong_pi(seq):
    s = seq.upper()
    return {"pI": _bisect(lambda x: _charge(s, x, termini=False), 0.0, 14.0), "charge_at_7": _charge(s, 7.0, False)}


def _rand_pi(rng):
    s = list(_rand_protein(rng, 4, 80))
    s[rng.randrange(len(s))] = rng.choice("KRH")
    s[rng.randrange(len(s))] = rng.choice("DE")
    return {"seq": "".join(s)}


_add("bio_isoelectric_point", name="Isoelectric point of a protein", topic="sequences / protein properties",
     purpose="The isoelectric point (pH of zero net charge) of a protein, with its net charge at pH 7 (Bjellqvist pK "
             "set).",
     inputs=[_textin("seq", "protein sequence")],
     outputs=[_o("pI", "isoelectric point"), _o("charge_at_7", "net charge at pH 7")], primary="pI",
     conditions=["standard amino acids only"],
     method="the net charge decreases with pH: bisection for the pH where it is zero, to full precision",
     latex=r"Z(\mathrm{pI}) = 0", run=isoelectric_point, reference=_ref_pi, wrong=_wrong_pi, rand=_rand_pi,
     phrases=["Isoelectric point", "pI of a protein"],
     examples=[_ex({"seq": "INGAR"}, {"pI": 9.7500246377, "charge_at_7": 0.760091614289}), _ex({"seq": "PETER"}, {"pI": 4.5320880784, "charge_at_7": -1.03586020838})],
     level="advanced")


KD = {"A": 1.8, "R": -4.5, "N": -3.5, "D": -3.5, "C": 2.5, "Q": -3.5, "E": -3.5, "G": -0.4, "H": -3.2, "I": 4.5,
      "L": 3.8, "K": -3.9, "M": 1.9, "F": 2.8, "P": -1.6, "S": -0.8, "T": -0.7, "W": -0.9, "Y": -1.3, "V": 4.2}


def gravy(seq):
    s = _seq(seq, "seq", "PROTEIN")
    return _out(gravy=sum(KD[a] for a in s) / len(s))


def _ref_gravy(seq):
    from Bio.SeqUtils.ProtParam import ProteinAnalysis
    return {"gravy": float(ProteinAnalysis(seq.upper()).gravy())}


_add("bio_gravy", name="Grand average of hydropathy (GRAVY)", topic="sequences / protein properties",
     purpose="The Kyte-Doolittle grand average of hydropathy of a protein (positive = hydrophobic).",
     inputs=[_textin("seq", "protein sequence")], outputs=[_o("gravy", "GRAVY score")], primary="gravy",
     conditions=["standard amino acids only"], method="mean Kyte-Doolittle hydropathy over the residues",
     latex=r"\mathrm{GRAVY} = \frac{1}{n}\sum_i h(a_i)", run=gravy, reference=_ref_gravy,
     wrong=lambda seq: {"gravy": sum(KD[a] for a in seq.upper())},
     rand=lambda rng: {"seq": _rand_protein(rng, 2, 100)},
     phrases=["GRAVY", "Grand average of hydropathy"],
     examples=[_ex({"seq": "MKTLLIL"}, {"gravy": 1.88571428571}), _ex({"seq": "DDEEKR"}, {"gravy": -3.73333333333})],
     level="intermediate")


DIWV_ROWS = {
    "A": "A1 C44.94 E1 D-7.49 G1 F1 I1 H-7.49 K1 M1 L1 N1 Q1 P20.26 S1 R1 T1 W1 V1 Y1",
    "C": "A1 C1 E1 D20.26 G1 F1 I1 H33.6 K1 M33.6 L20.26 N1 Q-6.54 P20.26 S1 R1 T33.6 W24.68 V-6.54 Y1",
    "E": "A1 C44.94 E33.6 D20.26 G1 F1 I20.26 H-6.54 K1 M1 L1 N1 Q20.26 P20.26 S20.26 R1 T1 W-14.03 V1 Y1",
    "D": "A1 C1 E1 D1 G1 F-6.54 I1 H1 K-7.49 M1 L1 N1 Q1 P1 S20.26 R-6.54 T-14.03 W1 V1 Y1",
    "G": "A-7.49 C1 E-6.54 D1 G13.34 F1 I-7.49 H1 K-7.49 M1 L1 N-7.49 Q1 P1 S1 R1 T-7.49 W13.34 V1 Y-7.49",
    "F": "A1 C1 E1 D13.34 G1 F1 I1 H1 K-14.03 M1 L1 N1 Q1 P20.26 S1 R1 T1 W1 V1 Y33.601",
    "I": "A1 C1 E44.94 D1 G1 F1 I1 H13.34 K-7.49 M1 L20.26 N1 Q1 P-1.88 S1 R1 T1 W1 V-7.49 Y1",
    "H": "A1 C1 E1 D1 G-9.37 F-9.37 I44.94 H1 K24.68 M1 L1 N24.68 Q1 P-1.88 S1 R1 T-6.54 W-1.88 V1 Y44.94",
    "K": "A1 C1 E1 D1 G-7.49 F1 I-7.49 H1 K1 M33.6 L-7.49 N1 Q24.64 P-6.54 S1 R33.6 T1 W1 V-7.49 Y1",
    "M": "A13.34 C1 E1 D1 G1 F1 I1 H58.28 K1 M-1.88 L1 N1 Q-6.54 P44.94 S44.94 R-6.54 T-1.88 W1 V1 Y24.68",
    "L": "A1 C1 E1 D1 G1 F1 I1 H1 K-7.49 M1 L1 N1 Q33.6 P20.26 S1 R20.26 T1 W24.68 V1 Y1",
    "N": "A1 C-1.88 E1 D1 G-14.03 F-14.03 I44.94 H1 K24.68 M1 L1 N1 Q-6.54 P-1.88 S1 R1 T-7.49 W-9.37 V1 Y1",
    "Q": "A1 C-6.54 E20.26 D20.26 G1 F-6.54 I1 H1 K1 M1 L1 N1 Q20.26 P20.26 S44.94 R1 T1 W1 V-6.54 Y-6.54",
    "P": "A20.26 C-6.54 E18.38 D-6.54 G1 F20.26 I1 H1 K1 M-6.54 L1 N1 Q20.26 P20.26 S20.26 R-6.54 T1 W-1.88 V20.26 Y1",
    "S": "A1 C33.6 E20.26 D1 G1 F1 I1 H1 K1 M1 L1 N1 Q20.26 P44.94 S20.26 R20.26 T1 W1 V1 Y1",
    "R": "A1 C1 E1 D1 G-7.49 F1 I1 H20.26 K1 M1 L1 N13.34 Q20.26 P20.26 S44.94 R58.28 T1 W58.28 V1 Y-6.54",
    "T": "A1 C1 E20.26 D1 G-7.49 F13.34 I1 H1 K1 M1 L1 N-14.03 Q-6.54 P1 S1 R1 T1 W-14.03 V1 Y1",
    "W": "A-14.03 C1 E1 D1 G-9.37 F1 I1 H24.68 K1 M24.68 L13.34 N13.34 Q1 P1 S1 R1 T-14.03 W1 V-7.49 Y1",
    "V": "A1 C1 E1 D-14.03 G-7.49 F1 I1 H1 K-1.88 M1 L1 N1 Q1 P20.26 S1 R1 T-7.49 W1 V1 Y-6.54",
    "Y": "A24.68 C1 E-6.54 D24.68 G-7.49 F1 I1 H13.34 K1 M44.94 L1 N1 Q1 P13.34 S1 R-15.91 T-7.49 W-9.37 V1 Y13.34"}
DIWV = {a: {t[0]: float(t[1:]) for t in row.split()} for a, row in DIWV_ROWS.items()}


def instability(seq):
    s = _seq(seq, "seq", "PROTEIN")
    _need(len(s) >= 2, "the protein needs at least 2 residues")
    score = 0.0
    for i in range(len(s) - 1):
        score += DIWV[s[i]][s[i + 1]]
    return _out(index=(10.0 / len(s)) * score)


def _ref_instability(seq):
    from Bio.SeqUtils.ProtParam import ProteinAnalysis
    return {"index": float(ProteinAnalysis(seq.upper()).instability_index())}


def _wrong_instability(seq):
    s = seq.upper()
    return {"index": 10.0 / (len(s) - 1) * sum(DIWV[s[i]][s[i + 1]] for i in range(len(s) - 1))}


_add("bio_instability_index", name="Protein instability index", topic="sequences / protein properties",
     purpose="The Guruprasad instability index of a protein from its dipeptides (above 40 suggests an unstable "
             "protein in vivo).",
     inputs=[_textin("seq", "protein sequence (at least 2 residues)")],
     outputs=[_o("index", "instability index")], primary="index", conditions=["standard amino acids, length >= 2"],
     method="sum of the dipeptide instability weights (DIWV) times 10/L",
     latex=r"II = \frac{10}{L}\sum_{i=1}^{L-1} \mathrm{DIWV}(a_i, a_{i+1})",
     run=instability, reference=_ref_instability, wrong=_wrong_instability,
     rand=lambda rng: {"seq": _rand_protein(rng, 3, 120)},
     phrases=["Instability index", "Protein instability index"],
     examples=[_ex({"seq": _PROT_EX}, {"index": 41.9802631579}), _ex({"seq": "AC"}, {"index": 224.7})],
     level="intermediate")


def tm_wallace(seq):
    s = _seq(seq, "seq", "DNA")
    return {"tm": 2 * (s.count("A") + s.count("T")) + 4 * (s.count("G") + s.count("C"))}


def _ref_tm_wallace(seq):
    from Bio.SeqUtils import MeltingTemp as mt
    return {"tm": float(mt.Tm_Wallace(seq.upper()))}


_add("bio_tm_wallace", name="Melting temperature by the Wallace rule", topic="sequences / primers",
     purpose="The approximate melting temperature of a short DNA primer (14-20 nt) by the Wallace rule.",
     inputs=[_textin("seq", "DNA primer sequence")], outputs=[_o("tm", "melting temperature in degC")], primary="tm",
     conditions=["unambiguous DNA letters"], method="2 degC per A or T plus 4 degC per G or C",
     latex=r"T_m = 2(n_A+n_T) + 4(n_G+n_C)\ ^\circ\mathrm{C}", run=tm_wallace, reference=_ref_tm_wallace,
     wrong=lambda seq: {"tm": 4 * (seq.upper().count("A") + seq.upper().count("T"))
                        + 2 * (seq.upper().count("G") + seq.upper().count("C"))},
     rand=lambda rng: {"seq": _rand_dna(rng, 8, 30)},
     phrases=["Wallace melting temperature", "Melting temperature by the Wallace rule"],
     examples=[_ex({"seq": "ACGTTGCAATGCCGTA"}, {"tm": 48}), _ex({"seq": "AAAAATTTTTGGGGG"}, {"tm": 40})],
     level="beginner")


DNA_NN3 = {"init": (0, 0), "init_A/T": (2.3, 4.1), "init_G/C": (0.1, -2.8), "init_oneG/C": (0, 0),
           "init_allA/T": (0, 0), "init_5T/A": (0, 0), "sym": (0, -1.4),
           "AA/TT": (-7.9, -22.2), "AT/TA": (-7.2, -20.4), "TA/AT": (-7.2, -21.3), "CA/GT": (-8.5, -22.7),
           "GT/CA": (-8.4, -22.4), "CT/GA": (-7.8, -21.0), "GA/CT": (-8.2, -22.2), "CG/GC": (-10.6, -27.2),
           "GC/CG": (-9.8, -24.4), "GG/CC": (-8.0, -19.9)}


def _tm_nn(s, Na, dnac1, dnac2, saltcorr=True):
    c = s.translate(_COMP)
    dh, ds = 0.0, 0.0
    t = DNA_NN3
    dh += t["init"][0]
    ds += t["init"][1]
    key = "init_oneG/C" if ("G" in s or "C" in s) else "init_allA/T"
    dh += t[key][0]
    ds += t[key][1]
    if s.startswith("T"):
        dh += t["init_5T/A"][0]
        ds += t["init_5T/A"][1]
    if s.endswith("A"):
        dh += t["init_5T/A"][0]
        ds += t["init_5T/A"][1]
    ends = s[0] + s[-1]
    at = ends.count("A") + ends.count("T")
    gc = ends.count("G") + ends.count("C")
    dh += t["init_A/T"][0] * at
    ds += t["init_A/T"][1] * at
    dh += t["init_G/C"][0] * gc
    ds += t["init_G/C"][1] * gc
    for i in range(len(s) - 1):
        nb = s[i:i + 2] + "/" + c[i:i + 2]
        v = t[nb] if nb in t else t[nb[::-1]]
        dh += v[0]
        ds += v[1]
    selfcomp = s == _revcomp(s)
    k = (dnac1 - dnac2 / 2.0) * 1e-9
    if selfcomp:
        k = dnac1 * 1e-9
        dh += t["sym"][0]
        ds += t["sym"][1]
    if saltcorr:
        ds += 0.368 * (len(s) - 1) * math.log(Na * 1e-3)
    return (1000 * dh) / (ds + 1.987 * math.log(k)) - 273.15


def tm_nn(seq, Na, dnac1, dnac2):
    s = _seq(seq, "seq", "DNA")
    _need(len(s) >= 2, "the sequence needs at least 2 bases")
    _pos("Na", Na)
    _pos("dnac1", dnac1)
    _fin("dnac2", dnac2)
    _need(dnac2 >= 0, "dnac2 must be zero or positive")
    _need(s == _revcomp(s) or dnac1 > dnac2 / 2, "dnac1 must exceed dnac2/2 (k = c1 - c2/2 must be positive)")
    return _out(tm=_tm_nn(s, Na, dnac1, dnac2))


def _ref_tm_nn(seq, Na, dnac1, dnac2):
    from Bio.SeqUtils import MeltingTemp as mt
    s = seq.upper()
    return {"tm": float(mt.Tm_NN(s, Na=Na, dnac1=dnac1, dnac2=dnac2, saltcorr=5, selfcomp=(s == _revcomp(s))))}


def _rand_tm_nn(rng):
    s = _rand_dna(rng, 10, 30)
    if rng.random() < 0.15:
        half = _rand_dna(rng, 5, 10)
        s = half + _revcomp(half)
    d1 = round(rng.uniform(50, 1000))
    return {"seq": s, "Na": round(rng.uniform(10, 300)), "dnac1": d1, "dnac2": round(rng.uniform(0, 1.5 * d1))}


_add("bio_tm_nearest_neighbor", name="Melting temperature by nearest-neighbour thermodynamics",
     topic="sequences / primers",
     purpose="The melting temperature of a DNA duplex (perfect complement) from nearest-neighbour thermodynamics "
             "(Allawi-SantaLucia 1997) with the SantaLucia 1998 sodium correction.",
     inputs=[_textin("seq", "DNA sequence 5'->3'"), _numin("Na", "sodium concentration in mM"),
             _numin("dnac1", "concentration of the more concentrated strand in nM"),
             _numin("dnac2", "concentration of the less concentrated strand in nM")],
     outputs=[_o("tm", "melting temperature in degC")], primary="tm",
     conditions=["unambiguous DNA, length >= 2", "Na > 0", "dnac1 > dnac2/2"],
     method="sum nearest-neighbour enthalpies and entropies with initiation terms; Tm = 1000 dH / (dS + R ln k) - "
            "273.15 with dS corrected by 0.368 (N-1) ln[Na+]",
     latex=r"T_m = \frac{1000\,\Delta H}{\Delta S + 0.368(N-1)\ln[\mathrm{Na^+}] + R\ln k} - 273.15",
     run=tm_nn, reference=_ref_tm_nn,
     wrong=lambda seq, Na, dnac1, dnac2: {"tm": _tm_nn(seq.upper(), Na, dnac1, dnac2, saltcorr=False)},
     rand=_rand_tm_nn,
     phrases=["Nearest-neighbour melting temperature", "Nearest-neighbor melting temperature"],
     examples=[_ex({"seq": "CGTTCCAAAGATGTGGGCATGAGCTTAC", "Na": 50, "dnac1": 25, "dnac2": 25}, {"tm": 60.3209194992}),
               _ex({"seq": "GCGAATTCGC", "Na": 100, "dnac1": 500, "dnac2": 0}, {"tm": 39.5272675469})],
     level="advanced")


# ================================================================================================ ALIGNMENT
_BLOSUM62_TEXT = """
   A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V
A  4 -1 -2 -2  0 -1 -1  0 -2 -1 -1 -1 -1 -2 -1  1  0 -3 -2  0
R -1  5  0 -2 -3  1  0 -2  0 -3 -2  2 -1 -3 -2 -1 -1 -3 -2 -3
N -2  0  6  1 -3  0  0  0  1 -3 -3  0 -2 -3 -2  1  0 -4 -2 -3
D -2 -2  1  6 -3  0  2 -1 -1 -3 -4 -1 -3 -3 -1  0 -1 -4 -3 -3
C  0 -3 -3 -3  9 -3 -4 -3 -3 -1 -1 -3 -1 -2 -3 -1 -1 -2 -2 -1
Q -1  1  0  0 -3  5  2 -2  0 -3 -2  1  0 -3 -1  0 -1 -2 -1 -2
E -1  0  0  2 -4  2  5 -2  0 -3 -3  1 -2 -3 -1  0 -1 -3 -2 -2
G  0 -2  0 -1 -3 -2 -2  6 -2 -4 -4 -2 -3 -3 -2  0 -2 -2 -3 -3
H -2  0  1 -1 -3  0  0 -2  8 -3 -3 -1 -2 -1 -2 -1 -2 -2  2 -3
I -1 -3 -3 -3 -1 -3 -3 -4 -3  4  2 -3  1  0 -3 -2 -1 -3 -1  3
L -1 -2 -3 -4 -1 -2 -3 -4 -3  2  4 -2  2  0 -3 -2 -1 -2 -1  1
K -1  2  0 -1 -3  1  1 -2 -1 -3 -2  5 -1 -3 -1  0 -1 -3 -2 -2
M -1 -1 -2 -3 -1  0 -2 -3 -2  1  2 -1  5  0 -2 -1 -1 -1 -1  1
F -2 -3 -3 -3 -2 -3 -3 -3 -1  0  0 -3  0  6 -4 -2 -2  1  3 -1
P -1 -2 -2 -1 -3 -1 -1 -2 -2 -3 -3 -1 -2 -4  7 -1 -1 -4 -3 -2
S  1 -1  1  0 -1  0  0  0 -1 -2 -2  0 -1 -2 -1  4  1 -3 -2 -2
T  0 -1  0 -1 -1 -1 -1 -2 -2 -1 -1 -1 -1 -2 -1  1  5 -2 -2  0
W -3 -3 -4 -4 -2 -2 -3 -2 -2 -3 -2 -3 -1  1 -4 -3 -2 11  2 -3
Y -2 -2 -2 -3 -2 -1 -2 -3  2 -1 -1 -2 -1  3 -3 -2 -2  2  7 -1
V  0 -3 -3 -3 -1 -2 -2 -3 -3  3  1 -2  1 -1 -2 -2  0 -3 -1  4
"""


def _blosum():
    lines = [ln.split() for ln in _BLOSUM62_TEXT.strip().splitlines()]
    cols = lines[0]
    return {(r[0], c): float(v) for r in lines[1:] for c, v in zip(cols, r[1:])}


BLOSUM62 = _blosum()
_NEG = float("-inf")


def _align(a, b, sc, op, ex, mode, wrong_gaps=False):
    """Optimal alignment score: mode 'global', 'local' or 'semiglobal' (end gaps free). A gap of length L scores
    op + (L-1)*ex (Gotoh's three-state recursion)."""
    n, m = len(a), len(b)
    o2 = op + ex if wrong_gaps else op
    H = [[0.0] * (m + 1) for _ in range(n + 1)]
    E = [[_NEG] * (m + 1) for _ in range(n + 1)]
    F = [[_NEG] * (m + 1) for _ in range(n + 1)]
    for i in range(1, n + 1):
        if mode == "global":
            H[i][0] = F[i][0] = o2 + (i - 1) * ex
    for j in range(1, m + 1):
        if mode == "global":
            H[0][j] = E[0][j] = o2 + (j - 1) * ex
    best = 0.0
    for i in range(1, n + 1):
        Hi, Hp, Ei, Fi, Fp = H[i], H[i - 1], E[i], F[i], F[i - 1]
        ai = a[i - 1]
        for j in range(1, m + 1):
            e = max(Ei[j - 1] + ex, Hi[j - 1] + o2)
            f = max(Fp[j] + ex, Hp[j] + o2)
            h = max(Hp[j - 1] + sc(ai, b[j - 1]), e, f)
            if mode == "local" and h < 0:
                h = 0.0
            Ei[j], Fi[j], Hi[j] = e, f, h
            if h > best:
                best = h
    if mode == "global":
        return H[n][m]
    if mode == "local":
        return best
    return max([H[n][j] for j in range(m + 1)] + [H[i][m] for i in range(n + 1)])


def _check_align(a, b, scores):
    _need(len(a) * len(b) <= 1_000_000, "sequences too long (length product at most 1 000 000)")
    for name, v in scores.items():
        _fin(name, v)


def _mm(match, mismatch):
    return lambda x, y: match if x == y else mismatch


def _aligner(mode, match=None, mismatch=None, op=None, ex=None, matrix=False, end_free=False):
    from Bio.Align import PairwiseAligner, substitution_matrices
    al = PairwiseAligner()
    al.mode = mode
    if matrix:
        al.substitution_matrix = substitution_matrices.load("BLOSUM62")
    else:
        al.match_score, al.mismatch_score = float(match), float(mismatch)
    al.open_gap_score, al.extend_gap_score = float(op), float(ex)
    if end_free:
        al.end_gap_score = 0.0
    return al


def _rand_align(rng, alphabet="ACGT", affine=False, blosum=False):
    a = _rand_dna(rng, 8, 30, alphabet)
    b = _mutate(rng, a, alphabet, 0.2, 0.2)
    if abs(len(a) - len(b)) < 2:
        b = b + "".join(rng.choice(alphabet) for _ in range(rng.randint(2, 4)))
    d = {"seq1": a, "seq2": b}
    if not blosum:
        d.update(match=rng.choice([1, 2, 3, 5]), mismatch=rng.choice([-1, -2, -3, -1.5]))
    if affine:
        d.update(open=rng.choice([-5, -6, -8, -10, -11, -4.5]), extend=rng.choice([-0.5, -1, -2]))
    else:
        d.update(gap=rng.choice([-1, -2, -3, -4, -2.5]))
    return d


def global_linear(seq1, seq2, match, mismatch, gap):
    a, b = _seq(seq1, "seq1", "ANY"), _seq(seq2, "seq2", "ANY")
    _check_align(a, b, {"match": match, "mismatch": mismatch, "gap": gap})
    _need(gap <= 0, "the gap score must be zero or negative")
    return _out(score=_align(a, b, _mm(match, mismatch), gap, gap, "global"))


_ALN_IN = [_textin("seq1", "first sequence"), _textin("seq2", "second sequence"),
           _numin("match", "score of a match"), _numin("mismatch", "score of a mismatch")]

_add("bio_global_alignment_score", name="Global alignment score (Needleman-Wunsch, linear gaps)",
     topic="alignment / pairwise",
     purpose="The optimal global (Needleman-Wunsch) alignment score of two sequences with stated match, mismatch and "
             "linear gap scores.",
     inputs=_ALN_IN + [_numin("gap", "score of each gap position (<= 0)")],
     outputs=[_o("score", "optimal global alignment score")], primary="score",
     conditions=["gap <= 0", "end gaps scored like internal gaps"],
     method="dynamic programming over prefixes (Needleman-Wunsch)",
     latex=r"S_{i,j} = \max(S_{i-1,j-1}+s(x_i,y_j),\ S_{i-1,j}+g,\ S_{i,j-1}+g)",
     run=global_linear, reference=lambda seq1, seq2, match, mismatch, gap: {
         "score": float(_aligner("global", match, mismatch, gap, gap).score(seq1.upper(), seq2.upper()))},
     wrong=lambda seq1, seq2, match, mismatch, gap: {
         "score": _align(seq1.upper(), seq2.upper(), _mm(match, mismatch), gap, gap, "semiglobal")},
     rand=lambda rng: _rand_align(rng, rng.choice(["ACGT", "ACGT", AMINO])),
     phrases=["Global alignment score", "Needleman-Wunsch score"],
     examples=[_ex({"seq1": "GATTACA", "seq2": "GCATGCU", "match": 1, "mismatch": -1, "gap": -1}, {"score": 0}),
               _ex({"seq1": "ACGT", "seq2": "ACT", "match": 2, "mismatch": -1, "gap": -2}, {"score": 4})],
     level="intermediate")


def global_affine(seq1, seq2, match, mismatch, open, extend):
    a, b = _seq(seq1, "seq1", "ANY"), _seq(seq2, "seq2", "ANY")
    _check_align(a, b, {"match": match, "mismatch": mismatch, "open": open, "extend": extend})
    _need(open <= 0 and extend <= 0, "gap scores must be zero or negative")
    return _out(score=_align(a, b, _mm(match, mismatch), open, extend, "global"))


_AFF_IN = [_numin("open", "score of the first position of a gap (<= 0)"),
           _numin("extend", "score of each further gap position (<= 0)")]

_add("bio_global_alignment_affine", name="Global alignment score with affine gaps (Gotoh)", topic="alignment / pairwise",
     purpose="The optimal global alignment score with affine gaps: a gap of length L scores open + (L-1)*extend.",
     inputs=_ALN_IN + _AFF_IN, outputs=[_o("score", "optimal global alignment score")], primary="score",
     conditions=["open, extend <= 0", "a gap of length L scores open + (L-1) extend", "end gaps scored"],
     method="Gotoh's three-state dynamic programming",
     latex=r"g(L) = o + (L-1)\,e", run=global_affine,
     reference=lambda seq1, seq2, match, mismatch, open, extend: {
         "score": float(_aligner("global", match, mismatch, open, extend).score(seq1.upper(), seq2.upper()))},
     wrong=lambda seq1, seq2, match, mismatch, open, extend: {
         "score": _align(seq1.upper(), seq2.upper(), _mm(match, mismatch), open, extend, "global", wrong_gaps=True)},
     rand=lambda rng: _rand_align(rng, "ACGT", affine=True),
     phrases=["Global alignment score with affine gaps", "Affine global alignment score"],
     examples=[_ex({"seq1": "ACGTTTACG", "seq2": "ACGACG", "match": 2, "mismatch": -1, "open": -5, "extend": -1},
                   {"score": 5}),
               _ex({"seq1": "AAAA", "seq2": "AA", "match": 1, "mismatch": -1, "open": -2, "extend": -0.5},
                   {"score": -0.5})],
     level="advanced")


def local_linear(seq1, seq2, match, mismatch, gap):
    a, b = _seq(seq1, "seq1", "ANY"), _seq(seq2, "seq2", "ANY")
    _check_align(a, b, {"match": match, "mismatch": mismatch, "gap": gap})
    _need(gap <= 0, "the gap score must be zero or negative")
    _need(match > 0, "the match score must be positive for a local alignment")
    return _out(score=_align(a, b, _mm(match, mismatch), gap, gap, "local"))


def _wrong_local(seq1, seq2, match, mismatch, gap):
    a, b = seq1.upper(), seq2.upper()                          # no floor at 0: the best prefix-to-anywhere score
    n, m = len(a), len(b)
    H = [[0.0] * (m + 1) for _ in range(n + 1)]
    for i in range(1, n + 1):
        H[i][0] = i * gap
    for j in range(1, m + 1):
        H[0][j] = j * gap
    best = _NEG
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            H[i][j] = max(H[i - 1][j - 1] + (match if a[i - 1] == b[j - 1] else mismatch), H[i - 1][j] + gap,
                          H[i][j - 1] + gap)
            best = max(best, H[i][j])
    return {"score": best}


def _rand_local(rng, affine=False, protein=False):
    alpha = AMINO if protein else "ACGT"
    core = _rand_dna(rng, 8, 20, alpha)
    a = _rand_dna(rng, 3, 15, alpha) + core + _rand_dna(rng, 0, 10, alpha)
    b = _rand_dna(rng, 0, 10, alpha) + _mutate(rng, core, alpha, 0.15, 0.2) + _rand_dna(rng, 3, 15, alpha)
    d = {"seq1": a, "seq2": b}
    if not protein:
        d.update(match=rng.choice([1, 2, 3, 5]), mismatch=rng.choice([-1, -2, -3, -1.5]))
    if affine:
        d.update(open=rng.choice([-2, -3, -4, -5, -2.5]), extend=rng.choice([-0.5, -1]))
    else:
        d.update(gap=rng.choice([-1, -2, -3, -1.5]))
    return d


_add("bio_local_alignment_score", name="Local alignment score (Smith-Waterman, linear gaps)",
     topic="alignment / pairwise",
     purpose="The optimal local (Smith-Waterman) alignment score of two sequences with stated match, mismatch and "
             "linear gap scores.",
     inputs=_ALN_IN + [_numin("gap", "score of each gap position (<= 0)")],
     outputs=[_o("score", "optimal local alignment score (>= 0)")], primary="score",
     conditions=["match > 0", "gap <= 0"], method="dynamic programming with a floor at 0 (Smith-Waterman)",
     latex=r"S_{i,j} = \max(0,\ S_{i-1,j-1}+s(x_i,y_j),\ S_{i-1,j}+g,\ S_{i,j-1}+g)",
     run=local_linear, reference=lambda seq1, seq2, match, mismatch, gap: {
         "score": float(_aligner("local", match, mismatch, gap, gap).score(seq1.upper(), seq2.upper()))},
     wrong=_wrong_local, rand=lambda rng: _rand_local(rng),
     phrases=["Local alignment score", "Smith-Waterman score"],
     examples=[_ex({"seq1": "TGTTACGG", "seq2": "GGTTGACTA", "match": 3, "mismatch": -3, "gap": -2}, {"score": 13}),
               _ex({"seq1": "AAAA", "seq2": "CCCC", "match": 1, "mismatch": -1, "gap": -1}, {"score": 0})],
     level="intermediate")


def local_affine(seq1, seq2, match, mismatch, open, extend):
    a, b = _seq(seq1, "seq1", "ANY"), _seq(seq2, "seq2", "ANY")
    _check_align(a, b, {"match": match, "mismatch": mismatch, "open": open, "extend": extend})
    _need(open <= 0 and extend <= 0, "gap scores must be zero or negative")
    _need(match > 0, "the match score must be positive for a local alignment")
    return _out(score=_align(a, b, _mm(match, mismatch), open, extend, "local"))


_add("bio_local_alignment_affine", name="Local alignment score with affine gaps", topic="alignment / pairwise",
     purpose="The optimal local alignment score with affine gaps (a gap of length L scores open + (L-1)*extend).",
     inputs=_ALN_IN + _AFF_IN, outputs=[_o("score", "optimal local alignment score (>= 0)")], primary="score",
     conditions=["match > 0", "open, extend <= 0"], method="Gotoh's recursion with a floor at 0",
     latex=r"g(L) = o + (L-1)\,e,\quad S \ge 0", run=local_affine,
     reference=lambda seq1, seq2, match, mismatch, open, extend: {
         "score": float(_aligner("local", match, mismatch, open, extend).score(seq1.upper(), seq2.upper()))},
     wrong=lambda seq1, seq2, match, mismatch, open, extend: {
         "score": _align(seq1.upper(), seq2.upper(), _mm(match, mismatch), open, extend, "local", wrong_gaps=True)},
     rand=lambda rng: _rand_local(rng, affine=True),
     phrases=["Local alignment score with affine gaps", "Affine local alignment score"],
     examples=[_ex({"seq1": "TTACGTACGTTT", "seq2": "GGACGTTACGTGG", "match": 2, "mismatch": -1, "open": -3,
                    "extend": -1}, {"score": 13}),
               _ex({"seq1": "AAAA", "seq2": "CCCC", "match": 1, "mismatch": -1, "open": -2, "extend": -1},
                   {"score": 0})],
     level="advanced")


def _blo(x, y):
    return BLOSUM62[(x, y)]


def protein_global(seq1, seq2, open, extend):
    a, b = _seq(seq1, "seq1", "PROTEIN"), _seq(seq2, "seq2", "PROTEIN")
    _check_align(a, b, {"open": open, "extend": extend})
    _need(open <= 0 and extend <= 0, "gap scores must be zero or negative")
    return _out(score=_align(a, b, _blo, open, extend, "global"))


_PROT_IN = [_textin("seq1", "first protein sequence"), _textin("seq2", "second protein sequence")]

_add("bio_protein_global_blosum62", name="Global protein alignment score (BLOSUM62, affine gaps)",
     topic="alignment / protein",
     purpose="The optimal global alignment score of two proteins with BLOSUM62 and affine gaps.",
     inputs=_PROT_IN + _AFF_IN, outputs=[_o("score", "optimal global alignment score")], primary="score",
     conditions=["standard amino acids", "open, extend <= 0"],
     method="Gotoh's dynamic programming with BLOSUM62 substitution scores",
     latex=r"S = \max_{\text{alignments}} \sum B_{62}(x,y) + \sum_{\text{gaps}} [o + (L-1)e]",
     run=protein_global, reference=lambda seq1, seq2, open, extend: {
         "score": float(_aligner("global", op=open, ex=extend, matrix=True).score(seq1.upper(), seq2.upper()))},
     wrong=lambda seq1, seq2, open, extend: {
         "score": _align(seq1.upper(), seq2.upper(), _blo, open, extend, "global", wrong_gaps=True)},
     rand=lambda rng: _rand_align(rng, AMINO, affine=True, blosum=True),
     phrases=["Global protein alignment score", "BLOSUM62 global alignment score"],
     examples=[_ex({"seq1": "HEAGAWGHEE", "seq2": "PAWHEAE", "open": -10, "extend": -1}, {"score": 3}),
               _ex({"seq1": "MKTAYIAK", "seq2": "MKTAYIAK", "open": -11, "extend": -1}, {"score": 39})],
     level="advanced")


def protein_local(seq1, seq2, open, extend):
    a, b = _seq(seq1, "seq1", "PROTEIN"), _seq(seq2, "seq2", "PROTEIN")
    _check_align(a, b, {"open": open, "extend": extend})
    _need(open <= 0 and extend <= 0, "gap scores must be zero or negative")
    return _out(score=_align(a, b, _blo, open, extend, "local"))


_add("bio_protein_local_blosum62", name="Local protein alignment score (BLOSUM62, affine gaps)",
     topic="alignment / protein",
     purpose="The optimal local (Smith-Waterman) alignment score of two proteins with BLOSUM62 and affine gaps.",
     inputs=_PROT_IN + _AFF_IN, outputs=[_o("score", "optimal local alignment score (>= 0)")], primary="score",
     conditions=["standard amino acids", "open, extend <= 0"],
     method="Gotoh's recursion with BLOSUM62 scores and a floor at 0",
     latex=r"S = \max\Big(0, \max_{\text{local alignments}} \textstyle\sum B_{62} + \sum [o + (L-1)e]\Big)",
     run=protein_local, reference=lambda seq1, seq2, open, extend: {
         "score": float(_aligner("local", op=open, ex=extend, matrix=True).score(seq1.upper(), seq2.upper()))},
     wrong=lambda seq1, seq2, open, extend: {
         "score": _align(seq1.upper(), seq2.upper(), _blo, open, extend, "local", wrong_gaps=True)},
     rand=lambda rng: _rand_local(rng, affine=True, protein=True),
     phrases=["Local protein alignment score", "BLOSUM62 local alignment score"],
     examples=[_ex({"seq1": "HEAGAWGHEE", "seq2": "PAWHEAE", "open": -8, "extend": -8}, {"score": 20}),
               _ex({"seq1": "WWWW", "seq2": "WW", "open": -11, "extend": -1}, {"score": 22})],
     level="expert")


def semiglobal(seq1, seq2, match, mismatch, gap):
    a, b = _seq(seq1, "seq1", "ANY"), _seq(seq2, "seq2", "ANY")
    _check_align(a, b, {"match": match, "mismatch": mismatch, "gap": gap})
    _need(gap <= 0, "the gap score must be zero or negative")
    return _out(score=_align(a, b, _mm(match, mismatch), gap, gap, "semiglobal"))


def _rand_semi(rng):
    alpha = "ACGT"
    core = _rand_dna(rng, 6, 20, alpha)
    d = {"seq1": _rand_dna(rng, 2, 10, alpha) + core + _rand_dna(rng, 2, 10, alpha),
         "seq2": _mutate(rng, core, alpha, 0.15, 0.1)}
    d.update(match=rng.choice([1, 2, 3]), mismatch=rng.choice([-1, -2, -3]), gap=rng.choice([-1, -2, -3]))
    return d


_add("bio_semiglobal_alignment_score", name="Semi-global (end-gap-free) alignment score", topic="alignment / pairwise",
     purpose="The optimal alignment score when gaps at either end of either sequence are free (fitting a read inside a "
             "longer sequence, overlaps), with linear internal gaps.",
     inputs=_ALN_IN + [_numin("gap", "score of each internal gap position (<= 0)")],
     outputs=[_o("score", "optimal semi-global alignment score")], primary="score",
     conditions=["gap <= 0", "end gaps unscored"],
     method="Needleman-Wunsch with a zero first row and column and the best score on the last row or column",
     latex=r"S = \max\big(\max_j S_{n,j},\ \max_i S_{i,m}\big),\ S_{i,0}=S_{0,j}=0",
     run=semiglobal, reference=lambda seq1, seq2, match, mismatch, gap: {
         "score": float(_aligner("global", match, mismatch, gap, gap, end_free=True).score(seq1.upper(), seq2.upper()))},
     wrong=lambda seq1, seq2, match, mismatch, gap: {
         "score": _align(seq1.upper(), seq2.upper(), _mm(match, mismatch), gap, gap, "global")},
     rand=_rand_semi,
     phrases=["Semi-global alignment score", "End-gap-free alignment score"],
     examples=[_ex({"seq1": "TTTTACGTACGGGG", "seq2": "ACGTAC", "match": 1, "mismatch": -1, "gap": -1}, {"score": 6}),
               _ex({"seq1": "ACGT", "seq2": "TTTT", "match": 1, "mismatch": -1, "gap": -2}, {"score": 1})],
     level="advanced")


# ================================================================================================ GENETICS
def _genotype(g, name):
    """'AaBb' -> [('A', 'a'), ('B', 'b')]: letter pairs, upper case first; each locus a different letter."""
    _need(isinstance(g, str) and re.fullmatch(r"[A-Za-z]+", g or "") is not None and len(g) % 2 == 0,
          f"{name} must be letter pairs, one pair per locus (e.g. AaBb)")
    loci = []
    for i in range(0, len(g), 2):
        x, y = g[i], g[i + 1]
        _need(x.upper() == y.upper(), f"{name}: '{x}{y}' is not one locus (both letters must name the same gene)")
        loci.append(tuple(sorted((x, y), key=lambda c: (c.islower(), c))))
    letters = [p[0].upper() for p in loci]
    _need(len(set(letters)) == len(letters), f"{name}: each locus must use a different letter")
    _need(len(loci) <= 6, f"{name}: at most 6 loci")
    return loci


def _same_loci(*gs):
    ref = [p[0].upper() for p in gs[0]]
    for g in gs[1:]:
        _need([p[0].upper() for p in g] == ref, "every genotype must list the same loci in the same order")


def _fr_text(f):
    return f"{f.numerator}/{f.denominator}"


def cross_genotype(parent1, parent2, offspring):
    a, b, o = _genotype(parent1, "parent1"), _genotype(parent2, "parent2"), _genotype(offspring, "offspring")
    _same_loci(a, b, o)
    p = Fraction(1)
    for (x1, y1), (x2, y2), want in zip(a, b, o):
        hits = sum(1 for u in (x1, y1) for v in (x2, y2) if tuple(sorted((u, v), key=lambda c: (c.islower(), c))) == want)
        p *= Fraction(hits, 4)
    return {"probability": float(p), "fraction": _fr_text(p)}


def _gametes(loci):
    return ["".join(t) for t in _iproduct(*loci)]


def _offspring_key(g1, g2):
    return "".join("".join(sorted((u, v), key=lambda c: (c.islower(), c))) for u, v in zip(g1, g2))


def _ref_cross_genotype(parent1, parent2, offspring):
    """Every gamete of each parent (2^n, equally likely) against every gamete of the other: a full Punnett square."""
    want = "".join("".join(p) for p in _genotype(offspring, "o"))
    g1, g2 = _gametes(_genotype(parent1, "a")), _gametes(_genotype(parent2, "b"))
    hits = sum(1 for x in g1 for y in g2 if _offspring_key(x, y) == want)
    f = Fraction(hits, len(g1) * len(g2))
    return {"probability": float(f), "fraction": _fr_text(f)}


def _wrong_cross_genotype(parent1, parent2, offspring):
    a, b, o = _genotype(parent1, "a"), _genotype(parent2, "b"), _genotype(offspring, "o")
    p = Fraction(1)
    for (x1, y1), (x2, y2), want in zip(a, b, o):                 # only one order of the two gametes counted
        p *= Fraction(sum(1 for u in (x1, y1) for v in (x2, y2) if (u, v) == want), 4)
    return {"probability": float(p), "fraction": _fr_text(p)}


def _rand_locus(rng, letter, het=0.7):
    if rng.random() < het:
        return letter.upper() + letter.lower()
    return rng.choice([letter.upper() * 2, letter.lower() * 2])


def _rand_cross(rng, with_target=True):
    n = rng.randint(1, 4)
    letters = rng.sample("ABCDEFGHRST", n)
    p1 = "".join(_rand_locus(rng, c) for c in letters)
    p2 = "".join(_rand_locus(rng, c) for c in letters)
    d = {"parent1": p1, "parent2": p2}
    if with_target:
        o = ""
        for i in range(n):
            u, v = rng.choice(p1[2 * i:2 * i + 2]), rng.choice(p2[2 * i:2 * i + 2])
            o += "".join(sorted((u, v), key=lambda c: (c.islower(), c)))
        d["offspring"] = o
    return d


_add("bio_cross_genotype_probability", name="Probability of an offspring genotype in a cross",
     topic="genetics / Mendelian inheritance",
     purpose="The probability that an offspring of two parents has a stated genotype (monohybrid, dihybrid or more "
             "loci; independent assortment).",
     inputs=[_textin("parent1", "genotype of the first parent, letter pairs (e.g. AaBb)", "seq"),
             _textin("parent2", "genotype of the second parent", "seq"),
             _textin("offspring", "offspring genotype asked about", "seq")],
     outputs=[_o("probability", "probability of the offspring genotype"), _o("fraction", "the same as a fraction", "text")],
     primary="probability", conditions=["the same loci in the same order", "independent assortment", "at most 6 loci"],
     method="product over loci of the Punnett-square probability of each locus genotype",
     latex=r"P = \prod_{\ell} P_\ell(g_\ell)", run=cross_genotype, reference=_ref_cross_genotype,
     wrong=_wrong_cross_genotype, rand=_rand_cross,
     phrases=["Cross genotype probability", "Offspring genotype probability"],
     examples=[_ex({"parent1": "Aa", "parent2": "Aa", "offspring": "Aa"}, {"probability": 0.5, "fraction": "1/2"}),
               _ex({"parent1": "AaBb", "parent2": "AaBb", "offspring": "AaBB"}, {"probability": 0.125, "fraction": "1/8"}),
               _ex({"parent1": "AaBbCc", "parent2": "aaBbCC", "offspring": "aabbCc"},
                   {"probability": 0.0625, "fraction": "1/16"})],
     extra_q=[{"q": "Cross genotype probability: parent1 = AaBb, parent2 = Aa, offspring = Aa", "outcome": "INVALID_INPUT"}],
     level="beginner")


def _locus_dist(x, y):
    d = {}
    for u in x:
        for v in y:
            k = tuple(sorted((u, v), key=lambda c: (c.islower(), c)))
            d[k] = d.get(k, 0) + Fraction(1, 4)
    return d


def _ratio_text(probs):
    probs = [p for p in probs if p > 0]
    den = 1
    for p in probs:
        den = den * p.denominator // math.gcd(den, p.denominator)
    ints = [int(p * den) for p in probs]
    g = 0
    for v in ints:
        g = math.gcd(g, v)
    return ":".join(str(v // g) for v in sorted(ints, reverse=True))


def cross_phenotypes(parent1, parent2, complete=True):
    a, b = _genotype(parent1, "parent1"), _genotype(parent2, "parent2")
    _same_loci(a, b)
    per = []
    n_geno = 1
    for (x1, y1), (x2, y2) in zip(a, b):
        d = _locus_dist((x1, y1), (x2, y2))
        n_geno *= len(d)
        if complete:
            dom = sum(p for k, p in d.items() if k[0].isupper())
            classes = [p for p in (dom, 1 - dom) if p > 0]
            per.append((dom, classes))
        else:
            per.append((d.get((x1.upper(), x1.upper()), Fraction(0)), list(d.values())))
    p_dom = Fraction(1)
    for dom, _c in per:
        p_dom *= dom
    probs = [Fraction(1)]
    for _dom, classes in per:
        probs = [p * c for p in probs for c in classes]
    p_rec = Fraction(1)
    for (x1, y1), (x2, y2) in zip(a, b):
        p_rec *= _locus_dist((x1, y1), (x2, y2)).get((x1.lower(), x1.lower()), Fraction(0))
    return {"p_all_dominant": float(p_dom), "p_all_recessive": float(p_rec), "n_genotypes": n_geno,
            "n_phenotypes": len(probs), "ratio": _ratio_text(probs)}


def _ref_cross_phenotypes(parent1, parent2):
    """A full Punnett square: every pair of gametes, phenotypes read from each offspring."""
    g1, g2 = _gametes(_genotype(parent1, "a")), _gametes(_genotype(parent2, "b"))
    genos, phenos = {}, {}
    for x in g1:
        for y in g2:
            k = _offspring_key(x, y)
            genos[k] = genos.get(k, 0) + 1
            ph = "".join("D" if k[2 * i].isupper() else "r" for i in range(len(k) // 2))
            phenos[ph] = phenos.get(ph, 0) + 1
    tot = len(g1) * len(g2)
    n = len(g1[0])
    g = 0
    for v in phenos.values():
        g = math.gcd(g, v)
    return {"p_all_dominant": phenos.get("D" * n, 0) / tot, "p_all_recessive": phenos.get("r" * n, 0) / tot,
            "n_genotypes": len(genos), "n_phenotypes": len(phenos),
            "ratio": ":".join(str(v // g) for v in sorted(phenos.values(), reverse=True))}


_add("bio_cross_phenotype_ratio", name="Phenotype ratio of a cross (complete dominance)",
     topic="genetics / Mendelian inheritance",
     purpose="The offspring phenotype ratio of a cross with complete dominance, the probabilities of the all-dominant "
             "and all-recessive phenotypes and the numbers of genotypes and phenotypes.",
     inputs=[_textin("parent1", "genotype of the first parent, letter pairs (e.g. AaBb)", "seq"),
             _textin("parent2", "genotype of the second parent", "seq")],
     outputs=[_o("p_all_dominant", "probability of showing every dominant trait"),
              _o("p_all_recessive", "probability of showing every recessive trait"),
              _o("n_genotypes", "number of distinct offspring genotypes"),
              _o("n_phenotypes", "number of distinct offspring phenotypes"),
              _o("ratio", "phenotype ratio, largest class first (e.g. 9:3:3:1)", "text")],
     primary="p_all_dominant", conditions=["complete dominance at every locus", "independent assortment"],
     method="per-locus phenotype probabilities (3/4 : 1/4 for Aa x Aa) multiplied across loci",
     latex=r"(\tfrac34 D + \tfrac14 r)^n \Rightarrow 9:3:3:1 \text{ for two loci}", run=cross_phenotypes,
     reference=_ref_cross_phenotypes, wrong=lambda parent1, parent2: cross_phenotypes(parent1, parent2, complete=False),
     rand=lambda rng: _rand_cross(rng, with_target=False),
     phrases=["Cross phenotype ratio", "Phenotype ratio of a cross"],
     examples=[_ex({"parent1": "AaBb", "parent2": "AaBb"},
                   {"p_all_dominant": 0.5625, "p_all_recessive": 0.0625, "n_genotypes": 9, "n_phenotypes": 4,
                    "ratio": "9:3:3:1"}),
               _ex({"parent1": "Aa", "parent2": "aa"},
                   {"p_all_dominant": 0.5, "p_all_recessive": 0.5, "n_genotypes": 2, "n_phenotypes": 2, "ratio": "1:1"}),
               _ex({"parent1": "AaBbCc", "parent2": "AaBbCc"}, {"p_all_dominant": 27 / 64, "ratio": "27:9:9:9:3:3:3:1"})],
     level="beginner")


def hardy_weinberg_test(n_AA, n_Aa, n_aa):
    for nm, v in (("n_AA", n_AA), ("n_Aa", n_Aa), ("n_aa", n_aa)):
        _need(isinstance(v, int) and v >= 0, f"{nm} must be a whole number >= 0")
    N = n_AA + n_Aa + n_aa
    _need(N >= 1, "at least one individual is needed")
    p = Fraction(2 * n_AA + n_Aa, 2 * N)
    _need(0 < p < 1, "both alleles must be present (otherwise the expected counts are zero)")
    q = 1 - p
    exp = [N * p * p, 2 * N * p * q, N * q * q]
    obs = [n_AA, n_Aa, n_aa]
    chi2 = sum((o - e) ** 2 / e for o, e in zip(obs, exp))
    ho, he = Fraction(n_Aa, N), 2 * p * q
    return _out(p=float(p), q=float(q), exp_AA=float(exp[0]), exp_Aa=float(exp[1]), exp_aa=float(exp[2]),
                chi2=float(chi2), p_value=math.erfc(math.sqrt(float(chi2) / 2)), Ho=float(ho), He=float(he),
                F=float(1 - ho / he))


def _ref_hwe(n_AA, n_Aa, n_aa):
    import allel
    import numpy as np
    from scipy.stats import chisquare
    g = allel.GenotypeArray([[[0, 0]] * n_AA + [[0, 1]] * n_Aa + [[1, 1]] * n_aa])
    af = g.count_alleles().to_frequencies()
    N = n_AA + n_Aa + n_aa
    alleles = np.array([0] * (2 * n_AA + n_Aa) + [1] * (n_Aa + 2 * n_aa))
    pairs = np.add.outer(alleles, alleles)                   # random union of gametes: every ordered pair
    exp = [N * float(np.mean(pairs == 0)), N * float(np.mean(pairs == 1)), N * float(np.mean(pairs == 2))]
    chi2, pv = chisquare([n_AA, n_Aa, n_aa], f_exp=exp, ddof=1)
    return {"p": float(af[0][0]), "q": float(af[0][1]), "exp_AA": exp[0], "exp_Aa": exp[1], "exp_aa": exp[2],
            "chi2": float(chi2), "p_value": float(pv), "Ho": float(allel.heterozygosity_observed(g)[0]),
            "He": float(allel.heterozygosity_expected(af, ploidy=2)[0]), "F": float(allel.inbreeding_coefficient(g)[0])}


def _wrong_hwe(n_AA, n_Aa, n_aa):
    N = n_AA + n_Aa + n_aa
    p = (2 * n_AA + n_Aa) / (2 * N)
    q = 1 - p
    exp = [N * p * p, N * p * q, N * q * q]                    # heterozygotes without the factor 2
    chi2 = sum((o - e) ** 2 / e for o, e in zip((n_AA, n_Aa, n_aa), exp))
    return {"p": p, "q": q, "exp_AA": exp[0], "exp_Aa": exp[1], "exp_aa": exp[2], "chi2": chi2,
            "p_value": math.erfc(math.sqrt(chi2 / 2)), "Ho": n_Aa / N, "He": p * q, "F": 1 - (n_Aa / N) / (p * q)}


def _rand_hwe(rng):
    N = rng.randint(30, 300)
    p = rng.uniform(0.2, 0.8)
    f = rng.uniform(-0.15, 0.25)
    pr = [p * p + f * p * (1 - p), 2 * p * (1 - p) * (1 - f), (1 - p) ** 2 + f * p * (1 - p)]
    c = [0, 0, 0]
    for _ in range(N):
        u = rng.random()
        c[0 if u < pr[0] else (1 if u < pr[0] + pr[1] else 2)] += 1
    return {"n_AA": c[0], "n_Aa": c[1], "n_aa": c[2]}


_HWE_IN = [_intin("n_AA", "number of AA individuals", lo=0), _intin("n_Aa", "number of Aa individuals", lo=0),
           _intin("n_aa", "number of aa individuals", lo=0)]

_add("bio_hardy_weinberg_test", name="Hardy-Weinberg chi-square test", topic="population genetics / Hardy-Weinberg",
     purpose="Allele frequencies, Hardy-Weinberg expected genotype counts, the chi-square test of equilibrium "
             "(1 degree of freedom), observed and expected heterozygosity and the inbreeding coefficient F.",
     inputs=_HWE_IN,
     outputs=[_o("p", "frequency of allele A"), _o("q", "frequency of allele a"),
              _o("exp_AA", "expected AA count (N p^2)"), _o("exp_Aa", "expected Aa count (2 N p q)"),
              _o("exp_aa", "expected aa count (N q^2)"), _o("chi2", "chi-square statistic"),
              _o("p_value", "p-value (chi-square, 1 degree of freedom)"), _o("Ho", "observed heterozygosity"),
              _o("He", "expected heterozygosity 2pq"), _o("F", "inbreeding coefficient 1 - Ho/He")],
     primary="chi2", conditions=["both alleles present", "biallelic locus", "no continuity correction"],
     method="allele counting, expected counts N p^2, 2Npq, N q^2, chi-square with 3 - 1 - 1 = 1 degree of freedom",
     latex=r"\chi^2 = \sum \frac{(O-E)^2}{E},\quad E = (Np^2,\ 2Npq,\ Nq^2)",
     run=hardy_weinberg_test, reference=_ref_hwe, wrong=_wrong_hwe, rand=_rand_hwe,
     phrases=["Hardy-Weinberg test", "Hardy-Weinberg chi-square test"],
     examples=[_ex({"n_AA": 298, "n_Aa": 489, "n_aa": 213}, {"p": 0.5425, "exp_Aa": 496.3875}),
               _ex({"n_AA": 50, "n_Aa": 0, "n_aa": 50}, {"chi2": 100.0, "Ho": 0.0, "He": 0.5, "F": 1.0})],
     extra_q=[{"q": "Hardy-Weinberg test: n_AA = 10, n_Aa = 0, n_aa = 0", "outcome": "INVALID_INPUT"}],
     level="intermediate")


def _hwe_probs(nA, N):
    """Exact probabilities of each heterozygote count given nA copies of A among N diploids (rational)."""
    nB = 2 * N - nA
    total = math.comb(2 * N, nA)
    out = {}
    for k in range(nA % 2, min(nA, nB) + 1, 2):
        nAA, nBB = (nA - k) // 2, (nB - k) // 2
        ways = math.factorial(N) // (math.factorial(nAA) * math.factorial(k) * math.factorial(nBB)) * 2 ** k
        out[k] = Fraction(ways, total)
    return out


def hwe_exact(n_AA, n_Aa, n_aa):
    for nm, v in (("n_AA", n_AA), ("n_Aa", n_Aa), ("n_aa", n_aa)):
        _need(isinstance(v, int) and v >= 0, f"{nm} must be a whole number >= 0")
    N = n_AA + n_Aa + n_aa
    _need(1 <= N <= 1000, "between 1 and 1000 individuals")
    nA = 2 * n_AA + n_Aa
    _need(0 < nA < 2 * N, "both alleles must be present")
    pr = _hwe_probs(nA, N)
    po = pr[n_Aa]
    return {"p_value": float(sum(v for v in pr.values() if v <= po)), "p_observed": float(po)}


def _ref_hwe_exact(n_AA, n_Aa, n_aa):
    """Brute force: every placement of the nA copies of A among the 2N gene slots is equally likely."""
    N = n_AA + n_Aa + n_aa
    nA = 2 * n_AA + n_Aa
    counts = {}
    for c in combinations(range(2 * N), nA):
        s = set(c)
        k = sum(1 for i in range(N) if (2 * i in s) != (2 * i + 1 in s))
        counts[k] = counts.get(k, 0) + 1
    tot = sum(counts.values())
    po = Fraction(counts[n_Aa], tot)
    return {"p_value": float(sum(Fraction(v, tot) for v in counts.values() if Fraction(v, tot) <= po)),
            "p_observed": float(po)}


def _wrong_hwe_exact(n_AA, n_Aa, n_aa):
    N = n_AA + n_Aa + n_aa
    pr = _hwe_probs(2 * n_AA + n_Aa, N)
    return {"p_value": float(sum(v for k, v in pr.items() if k <= n_Aa)), "p_observed": float(pr[n_Aa])}


def _rand_hwe_exact(rng):
    while True:
        N = rng.randint(3, 7)
        c = [0, 0, 0]
        for _ in range(N):
            c[rng.choice([0, 1, 1, 2])] += 1
        if 0 < 2 * c[0] + c[1] < 2 * N:
            return {"n_AA": c[0], "n_Aa": c[1], "n_aa": c[2]}


_add("bio_hardy_weinberg_exact_test", name="Hardy-Weinberg exact test", topic="population genetics / Hardy-Weinberg",
     purpose="The exact (Wigginton) test of Hardy-Weinberg equilibrium for a biallelic locus: the two-sided p-value "
             "and the probability of the observed number of heterozygotes given the allele counts.",
     inputs=_HWE_IN,
     outputs=[_o("p_value", "exact two-sided p-value (sum of probabilities not above the observed one)"),
              _o("p_observed", "probability of the observed heterozygote count")],
     primary="p_value", conditions=["both alleles present", "N <= 1000"],
     method="exact conditional distribution of heterozygote counts given the allele counts (rational arithmetic)",
     latex=r"P(n_{AB}\mid n_A) = \frac{N!}{n_{AA}!\,n_{AB}!\,n_{BB}!}\,\frac{2^{n_{AB}}\,n_A!\,n_B!}{(2N)!}",
     run=hwe_exact, reference=_ref_hwe_exact, wrong=_wrong_hwe_exact, rand=_rand_hwe_exact,
     phrases=["Hardy-Weinberg exact test", "Exact test of Hardy-Weinberg equilibrium"],
     examples=[_ex({"n_AA": 2, "n_Aa": 0, "n_aa": 2}, {"p_value": 6 / 70, "p_observed": 6 / 70}),
               _ex({"n_AA": 1, "n_Aa": 2, "n_aa": 1}, {"p_value": 1.0, "p_observed": 48 / 70})],
     level="expert")


def haploid_selection(p0, s, t):
    _frac("p0", p0, True, True)
    _fin("s", s)
    _need(s > -1, "s must be greater than -1 (fitness 1 + s must be positive)")
    _need(isinstance(t, int) and 0 <= t <= 100000, "t must be a whole number of generations, 0 to 100000")
    p = p0
    for _ in range(t):
        p = p * (1 + s) / (1 + s * p)
    return _out(p_t=p)


def _ref_haploid(p0, s, t):
    lo = math.log(p0 / (1 - p0)) + t * math.log1p(s)          # closed form: the odds grow by (1+s) per generation
    return {"p_t": 1 / (1 + math.exp(-lo)) if lo > -700 else math.exp(lo)}


def _wrong_haploid(p0, s, t):
    lo = math.log(p0 / (1 - p0)) + (t - 1) * math.log1p(s)
    return {"p_t": 1 / (1 + math.exp(-lo))}


_add("bio_haploid_selection", name="Allele frequency under haploid (genic) selection",
     topic="population genetics / natural selection",
     purpose="The frequency of allele A after t generations of haploid selection (relative fitness 1 + s against 1).",
     inputs=[_numin("p0", "initial frequency of A"), _numin("s", "selection coefficient (fitness of A = 1 + s)"),
             _intin("t", "number of generations", lo=0)],
     outputs=[_o("p_t", "frequency of A after t generations")], primary="p_t",
     conditions=["0 < p0 < 1", "s > -1", "discrete generations"],
     method="iterate p' = p(1+s)/(1+sp) generation by generation",
     latex=r"p' = \frac{p(1+s)}{1+sp}", run=haploid_selection, reference=_ref_haploid, wrong=_wrong_haploid,
     rand=lambda rng: {"p0": round(rng.uniform(0.05, 0.95), 3), "s": round(rng.uniform(-0.2, 0.2), 3),
                       "t": rng.randint(1, 40)},
     phrases=["Haploid selection", "Genic selection allele frequency"],
     examples=[_ex({"p0": 0.5, "s": 1, "t": 1}, {"p_t": 2 / 3}),
               _ex({"p0": 0.1, "s": 0.1, "t": 10}, {"p_t": 0.223719171702})],
     level="intermediate")


def recessive_lethal(q0, t):
    _frac("q0", q0, True, False)
    _need(isinstance(t, int) and 0 <= t <= 100000, "t must be a whole number of generations, 0 to 100000")
    q = q0
    for _ in range(t):
        q = q / (1 + q)
    return _out(q_t=q)


_add("bio_recessive_lethal", name="Decline of a recessive lethal allele", topic="population genetics / natural selection",
     purpose="The frequency of a recessive lethal allele after t generations of selection (aa never reproduce).",
     inputs=[_numin("q0", "initial frequency of the lethal allele"), _intin("t", "number of generations", lo=0)],
     outputs=[_o("q_t", "frequency after t generations")], primary="q_t",
     conditions=["0 < q0 <= 1", "random mating, aa fitness 0"], method="iterate q' = q/(1+q)",
     latex=r"q' = \frac{q}{1+q}", run=recessive_lethal,
     reference=lambda q0, t: {"q_t": q0 / (1 + t * q0)},
     wrong=lambda q0, t: {"q_t": q0 / (1 + (t + 1) * q0)},
     rand=lambda rng: {"q0": round(rng.uniform(0.01, 0.9), 3), "t": rng.randint(1, 500)},
     phrases=["Recessive lethal allele frequency", "Recessive lethal"],
     examples=[_ex({"q0": 0.5, "t": 1}, {"q_t": 1 / 3}), _ex({"q0": 0.1, "t": 10}, {"q_t": 0.05})],
     level="intermediate")


def het_advantage(w_AA, w_Aa, w_aa):
    _pos("fitness", w_AA, w_Aa, w_aa)
    _need(w_Aa > w_AA and w_Aa > w_aa, "heterozygote advantage needs w_Aa greater than both homozygote fitnesses")
    return _out(p_eq=(w_Aa - w_aa) / (2 * w_Aa - w_AA - w_aa))


def _ref_het_advantage(w_AA, w_Aa, w_aa):
    p = 0.5                                                   # forward iteration of the selection recursion
    for _ in range(2_000_000):
        wbar = p * p * w_AA + 2 * p * (1 - p) * w_Aa + (1 - p) ** 2 * w_aa
        nxt = (p * p * w_AA + p * (1 - p) * w_Aa) / wbar
        if abs(nxt - p) <= 2e-16 * p:                         # a fixed point to rounding level
            p = nxt
            break
        p = nxt
    return {"p_eq": p}


_add("bio_heterozygote_advantage", name="Equilibrium under heterozygote advantage",
     topic="population genetics / natural selection",
     purpose="The stable equilibrium frequency of allele A when the heterozygote is fittest (overdominance).",
     inputs=[_numin("w_AA", "relative fitness of AA"), _numin("w_Aa", "relative fitness of Aa"),
             _numin("w_aa", "relative fitness of aa")],
     outputs=[_o("p_eq", "equilibrium frequency of A")], primary="p_eq",
     conditions=["w_Aa > w_AA and w_Aa > w_aa", "fitnesses positive"],
     method="the interior fixed point of the selection recursion",
     latex=r"\hat p = \frac{w_{Aa} - w_{aa}}{2w_{Aa} - w_{AA} - w_{aa}} = \frac{t}{s+t}",
     run=het_advantage, reference=_ref_het_advantage,
     wrong=lambda w_AA, w_Aa, w_aa: {"p_eq": (w_Aa - w_AA) / (2 * w_Aa - w_AA - w_aa)},
     rand=lambda rng: {"w_AA": round(1 - rng.uniform(0.05, 0.9), 3), "w_Aa": 1.0, "w_aa": round(1 - rng.uniform(0.05, 0.9), 3)},
     phrases=["Heterozygote advantage equilibrium", "Overdominance equilibrium"],
     examples=[_ex({"w_AA": 0.8, "w_Aa": 1, "w_aa": 0.6}, {"p_eq": 2 / 3}),
               _ex({"w_AA": 0.89, "w_Aa": 1, "w_aa": 0.2}, {"p_eq": 0.879120879121})],
     extra_q=[{"q": "Heterozygote advantage equilibrium: w_AA = 1, w_Aa = 0.9, w_aa = 0.8", "outcome": "INVALID_INPUT"}],
     level="advanced")


def _mutsel_map(q, mu, s, h):
    p = 1 - q
    wbar = p * p + 2 * p * q * (1 - h * s) + q * q * (1 - s)
    qs = (p * q * (1 - h * s) + q * q * (1 - s)) / wbar
    return qs + mu * (1 - qs)


def mutation_selection(mu, s, h):
    _fin("mu", mu)
    _need(0 < mu <= 0.01, "mu must lie in (0, 0.01]")
    _need(isinstance(s, (int, float)) and 0.01 <= s <= 1, "s must lie in [0.01, 1]")
    _need(isinstance(h, (int, float)) and 0 <= h <= 1, "h must lie in [0, 1]")
    f = lambda q: _mutsel_map(q, mu, s, h) - q
    _need(f(0.5) < 0, "mutation is too strong relative to selection: no low-frequency balance")
    return _out(q_eq=_bisect(f, 0.0, 0.5))


def _ref_mutsel(mu, s, h):
    q = 0.0                                                   # forward iteration from q = 0 to the fixed point
    for _ in range(5_000_000):
        nxt = _mutsel_map(q, mu, s, h)
        if abs(nxt - q) <= 1e-17 * max(q, 1e-300) or nxt == q:
            q = nxt
            break
        q = nxt
    return {"q_eq": q}


_add("bio_mutation_selection_balance", name="Mutation-selection balance", topic="population genetics / mutation",
     purpose="The equilibrium frequency of a deleterious allele maintained by recurrent mutation against selection "
             "(fitness 1, 1 - hs, 1 - s), from the exact recursion.",
     inputs=[_numin("mu", "mutation rate A -> a per generation"), _numin("s", "selection coefficient against aa"),
             _numin("h", "dominance coefficient (0 recessive, 0.5 additive)")],
     outputs=[_o("q_eq", "equilibrium frequency of the deleterious allele")], primary="q_eq",
     conditions=["0 < mu <= 0.01", "0.01 <= s <= 1", "0 <= h <= 1", "selection then mutation, no back mutation"],
     method="the root of q' - q = 0 for selection followed by mutation, by bisection (approximately sqrt(mu/s) "
            "for a recessive and mu/(hs) otherwise)",
     latex=r"\hat q \approx \sqrt{\mu/s}\ (h=0),\quad \hat q \approx \frac{\mu}{hs}\ (h>0)",
     run=mutation_selection, reference=_ref_mutsel,
     wrong=lambda mu, s, h: {"q_eq": math.sqrt(mu / s) if h == 0 else mu / (h * s)},
     rand=lambda rng: {"mu": _sig(_lu(rng, 1e-6, 1e-3), 3), "s": round(rng.uniform(0.05, 1), 3),
                       "h": 0.0 if rng.random() < 0.3 else round(rng.uniform(0.05, 0.5), 3)},
     phrases=["Mutation-selection balance", "Mutation selection equilibrium"],
     examples=[_ex({"mu": 0.0001, "s": 1, "h": 0}, {"q_eq": 0.01}),
               _ex({"mu": 0.00001, "s": 0.1, "h": 0.5}, {"q_eq": 0.00019999800002})],
     level="expert")


def drift_heterozygosity(N, i, t):
    _need(isinstance(N, int) and 1 <= N <= 100000, "N must be a whole number of diploid individuals (1 to 100000)")
    _need(isinstance(i, int) and 0 < i < 2 * N, "i must be a whole number of copies between 1 and 2N - 1")
    _need(isinstance(t, int) and 0 <= t <= 100000, "t must be a whole number of generations")
    p = Fraction(i, 2 * N)
    h0 = float(2 * p * (1 - p))
    ret = (1 - 1 / (2 * N)) ** t
    return _out(H0=h0, H_t=h0 * ret, retained=ret)


def _ref_drift(N, i, t):
    """The Wright-Fisher Markov chain: binomial sampling of 2N gene copies each generation."""
    import numpy as np
    from scipy.stats import binom
    M = 2 * N
    k = np.arange(M + 1)
    P = np.array([binom.pmf(k, M, j / M) for j in range(M + 1)])
    v = np.zeros(M + 1)
    v[i] = 1.0
    for _ in range(t):
        v = v @ P
    h = 2 * (k / M) * (1 - k / M)
    h0 = float(h[i])
    ht = float(v @ h)
    return {"H0": h0, "H_t": ht, "retained": ht / h0}


def _wrong_drift(N, i, t):
    p = i / (2 * N)
    h0 = 2 * p * (1 - p)
    ret = (1 - 1 / N) ** t if N > 1 else 0.0                   # N gene copies instead of 2N
    return {"H0": h0, "H_t": h0 * ret, "retained": ret}


def _rand_drift(rng):
    N = rng.randint(2, 50)
    return {"N": N, "i": rng.randint(1, 2 * N - 1), "t": rng.randint(1, 120)}


_add("bio_drift_heterozygosity", name="Loss of heterozygosity by genetic drift", topic="population genetics / drift",
     purpose="The expected heterozygosity after t generations of genetic drift in a Wright-Fisher population of N "
             "diploids, starting from i copies of an allele.",
     inputs=[_intin("N", "number of diploid individuals", lo=1), _intin("i", "initial copies of allele A (of 2N)", lo=1),
             _intin("t", "number of generations", lo=0)],
     outputs=[_o("H0", "initial heterozygosity 2pq"), _o("H_t", "expected heterozygosity after t generations"),
              _o("retained", "fraction of heterozygosity retained, (1 - 1/(2N))^t")],
     primary="H_t", conditions=["0 < i < 2N", "constant N, random mating, no mutation or selection"],
     method="expected heterozygosity falls by the factor 1 - 1/(2N) each generation",
     latex=r"H_t = H_0\left(1 - \frac{1}{2N}\right)^t", run=drift_heterozygosity, reference=_ref_drift,
     wrong=_wrong_drift, rand=_rand_drift,
     phrases=["Drift heterozygosity", "Heterozygosity loss by drift"],
     examples=[_ex({"N": 10, "i": 10, "t": 1}, {"H0": 0.5, "H_t": 0.475, "retained": 0.95}),
               _ex({"N": 50, "i": 20, "t": 100}, {"H0": 0.32, "H_t": 0.117130349207, "retained": 0.366032341273})],
     extra_q=[{"q": "Drift heterozygosity: N = 10, i = 20, t = 5", "outcome": "INVALID_INPUT"}],
     level="intermediate")


def _log_abs_expm1(z):
    """log|e^z - 1| without overflow."""
    if z > 0:
        return z + math.log1p(-math.exp(-z)) if z > 30 else math.log(math.expm1(z))
    return math.log(-math.expm1(z))


def moran_fixation(N, i, r):
    _need(isinstance(N, int) and 2 <= N <= 100000, "N must be a whole number of individuals (2 to 100000)")
    _need(isinstance(i, int) and 0 < i < N, "i must be a whole number of mutants between 1 and N - 1")
    _pos("r", r)
    if r == 1:
        return _out(rho=i / N)
    x = math.log(r)
    return _out(rho=math.exp(_log_abs_expm1(-i * x) - _log_abs_expm1(-N * x)))


def _ref_moran(N, i, r):
    """The Moran birth-death chain: absorption probabilities from a linear solve (numpy)."""
    import numpy as np
    A = np.zeros((N + 1, N + 1))
    b = np.zeros(N + 1)
    A[0, 0] = A[N, N] = 1.0
    b[N] = 1.0
    for j in range(1, N):
        tot = r * j + (N - j)
        up = r * j / tot * (N - j) / N
        dn = (N - j) / tot * j / N
        A[j, j - 1], A[j, j], A[j, j + 1] = -dn, up + dn, -up
    return {"rho": float(np.linalg.solve(A, b)[i])}


def _wrong_moran(N, i, r):
    if r == 1:
        return {"rho": (i + 1) / (N + 1)}
    return {"rho": (1 - r ** i) / (1 - r ** N)}                 # the fitness ratio inverted


def _rand_moran(rng):
    N = rng.randint(3, 40)
    return {"N": N, "i": rng.randint(1, N - 1), "r": 1.0 if rng.random() < 0.1 else round(rng.uniform(0.8, 1.5), 3)}


_add("bio_moran_fixation", name="Fixation probability in the Moran process", topic="population genetics / drift",
     purpose="The probability that a mutant with relative fitness r, starting from i copies in a Moran population of "
             "size N, eventually takes over the population.",
     inputs=[_intin("N", "population size", lo=2), _intin("i", "initial number of mutants", lo=1),
             _numin("r", "relative fitness of the mutant")],
     outputs=[_o("rho", "fixation probability")], primary="rho", conditions=["0 < i < N", "r > 0"],
     method="gambler's-ruin solution of the Moran birth-death chain (i/N when r = 1)",
     latex=r"\rho_i = \frac{1 - r^{-i}}{1 - r^{-N}}", run=moran_fixation, reference=_ref_moran, wrong=_wrong_moran,
     rand=_rand_moran, phrases=["Moran fixation probability", "Fixation probability in the Moran process"],
     examples=[_ex({"N": 10, "i": 1, "r": 1}, {"rho": 0.1}),
               _ex({"N": 10, "i": 1, "r": 2}, {"rho": 0.500488758553})],
     level="expert")


def ne_sex_ratio(N_m, N_f):
    _need(isinstance(N_m, int) and isinstance(N_f, int) and N_m >= 1 and N_f >= 1,
          "N_m and N_f must be whole numbers of breeding males and females (at least 1)")
    ne = Fraction(4 * N_m * N_f, N_m + N_f)
    return {"Ne": float(ne), "Ne_over_N": float(ne / (N_m + N_f))}


def _ref_ne_sex(N_m, N_f):
    """Enumerate the two-sex Wright-Fisher model: two genes from different offspring, each from a random parent of a
    random sex (1/2 each) and a random one of that parent's two copies; 1/(2 Ne) is the chance they are one gene."""
    same = Fraction(0)
    for nsex in (N_m, N_f):                                    # both genes from parents of this sex
        for a in range(nsex):
            for b in range(nsex):
                if a == b:
                    for c1 in range(2):
                        for c2 in range(2):
                            if c1 == c2:
                                same += Fraction(1, 4) * Fraction(1, nsex * nsex) * Fraction(1, 4)
    ne = 1 / (2 * same)
    return {"Ne": float(ne), "Ne_over_N": float(ne / (N_m + N_f))}


_add("bio_ne_sex_ratio", name="Effective population size with unequal sex ratio", topic="population genetics / drift",
     purpose="The effective population size of a population with N_m breeding males and N_f breeding females.",
     inputs=[_intin("N_m", "number of breeding males", lo=1), _intin("N_f", "number of breeding females", lo=1)],
     outputs=[_o("Ne", "effective population size"), _o("Ne_over_N", "Ne divided by the census N_m + N_f")],
     primary="Ne", conditions=["random mating within the breeding population", "Wright-Fisher reproduction"],
     method="Wright's formula", latex=r"N_e = \frac{4N_mN_f}{N_m+N_f}", run=ne_sex_ratio, reference=_ref_ne_sex,
     wrong=lambda N_m, N_f: {"Ne": 2 * N_m * N_f / (N_m + N_f), "Ne_over_N": 2 * N_m * N_f / (N_m + N_f) ** 2},
     rand=lambda rng: {"N_m": rng.randint(1, 40), "N_f": rng.randint(1, 40)},
     phrases=["Effective population size with unequal sex ratio", "Sex-ratio effective population size"],
     examples=[_ex({"N_m": 10, "N_f": 90}, {"Ne": 36.0, "Ne_over_N": 0.36}),
               _ex({"N_m": 50, "N_f": 50}, {"Ne": 100.0, "Ne_over_N": 1.0})],
     level="intermediate")


def ne_fluctuating(sizes):
    _need(isinstance(sizes, list) and len(sizes) >= 2, "give the sizes of at least 2 generations")
    ns = [_whole("sizes", v, 1) for v in sizes]
    ne = Fraction(len(ns)) / sum(Fraction(1, n) for n in ns)
    ret = 1.0
    for n in ns:
        ret *= 1 - 1 / (2 * n)
    return _out(Ne=float(ne), retained=ret)


def _ref_ne_fluct(sizes):
    import numpy as np
    from scipy.stats import binom, hmean
    ns = [int(v) for v in sizes]
    m0 = 20                                                     # any starting population: the ratio does not depend on it
    v = np.zeros(m0 + 1)
    v[m0 // 2] = 1.0
    h0 = 0.5
    M = m0
    for n in ns:                                               # Wright-Fisher sampling into a generation of n diploids
        M2 = 2 * n
        P = np.array([binom.pmf(np.arange(M2 + 1), M2, j / M) for j in range(M + 1)])
        v = v @ P
        M = M2
    k = np.arange(M + 1)
    ht = float(v @ (2 * (k / M) * (1 - k / M)))
    return {"Ne": float(hmean(ns)), "retained": ht / h0}


def _wrong_ne_fluct(sizes):
    m = sum(sizes) / len(sizes)
    return {"Ne": m, "retained": (1 - 1 / (2 * m)) ** len(sizes)}


_add("bio_ne_fluctuating_size", name="Effective size of a fluctuating population", topic="population genetics / drift",
     purpose="The effective population size over generations of changing size (harmonic mean) and the fraction of "
             "heterozygosity retained over those generations.",
     inputs=[_vecin("sizes", "numbers of diploid individuals in successive generations", lo=2, hi=200)],
     outputs=[_o("Ne", "effective size (harmonic mean of the sizes)"),
              _o("retained", "heterozygosity retained, product of (1 - 1/(2N_i))")],
     primary="Ne", conditions=["whole positive sizes", "Wright-Fisher reproduction each generation"],
     method="harmonic mean of the generation sizes; heterozygosity falls by 1 - 1/(2N_i) in generation i",
     latex=r"\frac{1}{N_e} = \frac{1}{t}\sum_{i=1}^t \frac{1}{N_i},\quad \frac{H_t}{H_0} = \prod_i \Big(1-\frac{1}{2N_i}\Big)",
     run=ne_fluctuating, reference=_ref_ne_fluct, wrong=_wrong_ne_fluct,
     rand=lambda rng: {"sizes": [rng.randint(2, 30) for _ in range(rng.randint(2, 6))]},
     phrases=["Effective size of a fluctuating population", "Fluctuating population effective size"],
     examples=[_ex({"sizes": [1000, 10, 1000]}, {"Ne": 3000 / 102, "retained": 0.9995 * 0.95 * 0.9995}),
               _ex({"sizes": [20, 20]}, {"Ne": 20.0, "retained": 0.950625})],
     level="intermediate")


def fst(p):
    _need(isinstance(p, list) and len(p) >= 2, "give the allele frequency in at least 2 subpopulations")
    for v in p:
        _frac("p", v)
    pb = sum(p) / len(p)
    _need(0 < pb < 1, "the allele must be polymorphic overall (0 < mean p < 1)")
    hs = sum(2 * v * (1 - v) for v in p) / len(p)
    ht = 2 * pb * (1 - pb)
    return _out(HS=hs, HT=ht, FST=(ht - hs) / ht)


def _ref_fst(p):
    """Wright's variance form: HT - HS = 2 Var(p), FST = Var(p) / (pbar qbar); HT from scikit-allel."""
    import allel
    import numpy as np
    a = np.array(p, dtype=float)
    pb = float(a.mean())
    var = float(a.var())
    ht = float(allel.heterozygosity_expected(np.array([[pb, 1 - pb]]), ploidy=2)[0])
    return {"HS": ht - 2 * var, "HT": ht, "FST": var / (pb * (1 - pb))}


_add("bio_fst", name="Fixation index FST from subpopulation allele frequencies",
     topic="population genetics / population structure",
     purpose="Mean within-subpopulation heterozygosity HS, total heterozygosity HT and FST = (HT - HS)/HT for a "
             "biallelic locus in equally weighted subpopulations.",
     inputs=[_vecin("p", "frequency of one allele in each subpopulation", lo=2, hi=100)],
     outputs=[_o("HS", "mean expected heterozygosity within subpopulations"),
              _o("HT", "expected heterozygosity of the pooled population"), _o("FST", "fixation index (HT - HS)/HT")],
     primary="FST", conditions=["equal subpopulation weights", "biallelic locus", "0 < mean p < 1"],
     method="HS = mean of 2p(1-p), HT = 2 pbar (1 - pbar), FST = (HT - HS)/HT",
     latex=r"F_{ST} = \frac{H_T - H_S}{H_T}", run=fst, reference=_ref_fst,
     wrong=lambda p: (lambda r: {"HS": r["HS"], "HT": r["HT"], "FST": (r["HT"] - r["HS"]) / r["HS"]})(fst(p)),
     rand=lambda rng: {"p": [round(rng.uniform(0.05, 0.95), 3) for _ in range(rng.randint(2, 6))]},
     phrases=["FST", "Fixation index FST"],
     examples=[_ex({"p": [0.2, 0.8]}, {"HS": 0.32, "HT": 0.5, "FST": 0.36}),
               _ex({"p": [0.5, 0.5, 0.5]}, {"FST": 0.0})],
     level="advanced")


def linkage_disequilibrium(n_AB, n_Ab, n_aB, n_ab):
    c = [n_AB, n_Ab, n_aB, n_ab]
    for v in c:
        _need(isinstance(v, int) and v >= 0, "haplotype counts must be whole numbers >= 0")
    n = sum(c)
    _need(n >= 1, "at least one haplotype is needed")
    pA, pB = Fraction(n_AB + n_Ab, n), Fraction(n_AB + n_aB, n)
    _need(0 < pA < 1 and 0 < pB < 1, "both loci must be polymorphic")
    D = Fraction(n_AB, n) - pA * pB
    if D > 0:
        dmax = min(pA * (1 - pB), (1 - pA) * pB)
    else:
        dmax = min(pA * pB, (1 - pA) * (1 - pB))
    return {"D": float(D), "D_prime": float(D / dmax) if D != 0 else 0.0,
            "r2": float(D * D / (pA * (1 - pA) * pB * (1 - pB)))}


def _ref_ld(n_AB, n_Ab, n_aB, n_ab):
    """D from the haplotype cross-products; r^2 = chi^2 / n (scipy); D' from Dmax found by a linear programme over all
    haplotype tables with the same allele frequencies (scipy linprog)."""
    from scipy.optimize import linprog
    from scipy.stats import chi2_contingency
    n = n_AB + n_Ab + n_aB + n_ab
    f = [x / n for x in (n_AB, n_Ab, n_aB, n_ab)]
    D = f[0] * f[3] - f[1] * f[2]
    chi2 = float(chi2_contingency([[n_AB, n_Ab], [n_aB, n_ab]], correction=False)[0])
    pA, pB = f[0] + f[1], f[0] + f[2]
    A_eq = [[1, 1, 0, 0], [1, 0, 1, 0], [1, 1, 1, 1]]
    b_eq = [pA, pB, 1.0]
    if D > 0:
        res = linprog([-1, 0, 0, 0], A_eq=A_eq, b_eq=b_eq, bounds=[(0, None)] * 4, method="highs")
        dmax = res.x[0] - pA * pB
    else:
        res = linprog([1, 0, 0, 0], A_eq=A_eq, b_eq=b_eq, bounds=[(0, None)] * 4, method="highs")
        dmax = pA * pB - res.x[0]
    return {"D": D, "D_prime": D / dmax if D != 0 else 0.0, "r2": chi2 / n}


def _wrong_ld(n_AB, n_Ab, n_aB, n_ab):
    r = linkage_disequilibrium(n_AB, n_Ab, n_aB, n_ab)
    n = n_AB + n_Ab + n_aB + n_ab
    pA, pB = (n_AB + n_Ab) / n, (n_AB + n_aB) / n
    return dict(r, r2=r["D"] ** 2 / (pA * pB))                # the denominator without qA qB


def _rand_ld(rng):
    while True:
        c = [rng.randint(1, 120) for _ in range(4)]
        if c[0] * c[3] != c[1] * c[2]:
            return dict(zip(("n_AB", "n_Ab", "n_aB", "n_ab"), c))


_add("bio_linkage_disequilibrium", name="Linkage disequilibrium D, D' and r^2",
     topic="population genetics / linkage disequilibrium",
     purpose="The linkage disequilibrium between two biallelic loci from haplotype counts: D, Lewontin's D' and r^2.",
     inputs=[_intin("n_AB", "number of AB haplotypes", lo=0), _intin("n_Ab", "number of Ab haplotypes", lo=0),
             _intin("n_aB", "number of aB haplotypes", lo=0), _intin("n_ab", "number of ab haplotypes", lo=0)],
     outputs=[_o("D", "disequilibrium D = pAB - pA pB"), _o("D_prime", "Lewontin's D' = D / Dmax (signed)"),
              _o("r2", "squared correlation r^2")],
     primary="r2", conditions=["both loci polymorphic", "phased haplotype counts"],
     method="haplotype and allele frequencies (exact fractions); Dmax from the allele frequencies",
     latex=r"D = p_{AB} - p_Ap_B,\quad r^2 = \frac{D^2}{p_Aq_Ap_Bq_B},\quad D' = \frac{D}{D_{\max}}",
     run=linkage_disequilibrium, reference=_ref_ld, wrong=_wrong_ld, rand=_rand_ld,
     phrases=["Linkage disequilibrium", "LD from haplotype counts"],
     examples=[_ex({"n_AB": 40, "n_Ab": 10, "n_aB": 10, "n_ab": 40}, {"D": 0.15, "D_prime": 0.6, "r2": 0.36}),
               _ex({"n_AB": 25, "n_Ab": 25, "n_aB": 25, "n_ab": 25}, {"D": 0.0, "D_prime": 0.0, "r2": 0.0})],
     level="advanced")


def map_distance(r, map):
    _fin("r", r)
    _need(0 <= r < 0.5, "the recombination fraction must lie in [0, 0.5)")
    _need(map in ("haldane", "kosambi"), "map must be haldane or kosambi")
    m = -0.5 * math.log1p(-2 * r) if map == "haldane" else 0.25 * math.log((1 + 2 * r) / (1 - 2 * r))
    return _out(d_cM=100 * m)


def _poisson_odd(m):
    from scipy.stats import poisson
    import numpy as np
    k = np.arange(1, 400, 2)
    return float(poisson.pmf(k, m).sum())


def _ref_map_distance(r, map):
    """Haldane: crossovers a Poisson process, r = P(odd number) inverted by brentq; Kosambi: quadrature of
    dd/dr = 1/(1 - 4r^2), its defining differential equation."""
    if r == 0:
        return {"d_cM": 0.0}
    if map == "haldane":
        return {"d_cM": 100 * _brentq(lambda m: _poisson_odd(m) - r, 0.0, 20.0)}
    return {"d_cM": 100 * _quad(lambda x: 1.0 / (1.0 - 4.0 * x * x), 0.0, r)}


_MAP = _choice("map", ["haldane", "kosambi"], "map function: haldane (no interference) or kosambi")

_add("bio_map_distance", name="Map distance from a recombination fraction", topic="genetics / linkage mapping",
     purpose="The genetic map distance (cM) corresponding to a recombination fraction, by the Haldane or Kosambi map "
             "function.",
     inputs=[_numin("r", "recombination fraction (0 to 0.5)"), _MAP],
     outputs=[_o("d_cM", "map distance in centimorgans")], primary="d_cM",
     conditions=["0 <= r < 0.5", "map function stated"],
     method="Haldane d = -ln(1-2r)/2; Kosambi d = ln((1+2r)/(1-2r))/4 (Morgans, x100 for cM)",
     latex=r"d_H = -\tfrac12\ln(1-2r),\quad d_K = \tfrac14\ln\frac{1+2r}{1-2r}", run=map_distance,
     reference=_ref_map_distance, wrong=lambda r, map: {"d_cM": 100 * r},
     rand=lambda rng: {"r": round(rng.uniform(0.005, 0.45), 4), "map": rng.choice(["haldane", "kosambi"])},
     phrases=["Map distance", "Map distance from recombination fraction"],
     examples=[_ex({"r": 0.2, "map": "haldane"}, {"d_cM": 25.5412811883}),
               _ex({"r": 0.2, "map": "kosambi"}, {"d_cM": 21.1824465097})],
     extra_q=[{"q": "Map distance: r = 0.2", "outcome": "NEEDS_INPUT"},
              {"q": "Map distance: r = 0.6, map = haldane", "outcome": "INVALID_INPUT"}],
     level="intermediate")


def recombination_fraction(d_cM, map):
    _fin("d_cM", d_cM)
    _need(d_cM >= 0, "the map distance must be zero or positive")
    _need(map in ("haldane", "kosambi"), "map must be haldane or kosambi")
    m = d_cM / 100
    return _out(r=-0.5 * math.expm1(-2 * m) if map == "haldane" else 0.5 * math.tanh(2 * m))


def _ref_recomb(d_cM, map):
    m = d_cM / 100
    if m == 0:
        return {"r": 0.0}
    if map == "haldane":
        return {"r": _poisson_odd(m)}
    return {"r": _ivp(lambda x, y: [1 - 4 * y[0] ** 2], m, [0.0])[0]}


_add("bio_recombination_fraction", name="Recombination fraction from a map distance", topic="genetics / linkage mapping",
     purpose="The recombination fraction between two loci a stated map distance apart, by the Haldane or Kosambi map "
             "function.",
     inputs=[_numin("d_cM", "map distance in centimorgans"), _MAP],
     outputs=[_o("r", "recombination fraction")], primary="r", conditions=["d >= 0", "map function stated"],
     method="Haldane r = (1 - e^{-2d})/2; Kosambi r = tanh(2d)/2 (d in Morgans)",
     latex=r"r_H = \tfrac12(1 - e^{-2d}),\quad r_K = \tfrac12\tanh 2d", run=recombination_fraction,
     reference=_ref_recomb, wrong=lambda d_cM, map: {"r": min(d_cM / 100, 0.5)},
     rand=lambda rng: {"d_cM": round(rng.uniform(0.5, 200), 2), "map": rng.choice(["haldane", "kosambi"])},
     phrases=["Recombination fraction", "Recombination fraction from map distance"],
     examples=[_ex({"d_cM": 50, "map": "haldane"}, {"r": 0.316060279414}),
               _ex({"d_cM": 50, "map": "kosambi"}, {"r": 0.380797077978})],
     level="intermediate")


def _pedigree(pedigree):
    parents = {}
    for e in pedigree:
        _need(len(e) == 2, "pedigree edges are parent-child pairs without weights")
        p, c = e
        _need(p != c, "an individual cannot be its own parent")
        parents.setdefault(c, [])
        _need(p not in parents[c], f"{p} is listed twice as a parent of {c}")
        parents[c].append(p)
    for c, ps in parents.items():
        _need(len(ps) == 2, f"{c} must have exactly two parents listed (founders have none)")
    seen, done = set(), set()

    def visit(x):                                             # no individual may be its own ancestor
        if x in done:
            return
        _need(x not in seen, "the pedigree has a cycle")
        seen.add(x)
        for p in parents.get(x, []):
            visit(p)
        done.add(x)
    for c in list(parents):
        visit(c)
    return parents


def pedigree_inbreeding(pedigree, individual):
    half = Fraction(1, 2)
    parents = _pedigree(pedigree)
    names = set(parents) | {p for ps in parents.values() for p in ps}
    _need(individual in names, f"{individual} is not in the pedigree")
    memo = {}

    def up_paths(x):
        out = [[x]]
        for p in parents.get(x, []):
            out += [[x] + q for q in up_paths(p)]
        return out

    def F(x):
        if x not in parents:
            return Fraction(0)
        if x in memo:
            return memo[x]
        s, d = parents[x]
        tot = Fraction(0)
        ps, pd = up_paths(s), up_paths(d)
        for a in ps:
            for b in pd:
                if a[-1] == b[-1] and set(a) & set(b) == {a[-1]}:
                    tot += half ** (len(a) + len(b) - 1) * (1 + F(a[-1]))
        memo[x] = tot
        return tot
    return {"F": float(F(individual))}


def _ref_pedigree(pedigree, individual):
    """The tabular method: the additive relationship matrix built in pedigree order; F = A_xx - 1."""
    import numpy as np
    parents = {}
    for p, c in pedigree:
        parents.setdefault(c, []).append(p)
    order, done = [], set()

    def visit(x):
        if x in done:
            return
        for p in parents.get(x, []):
            visit(p)
        done.add(x)
        order.append(x)
    for x in sorted(set(parents) | {p for ps in parents.values() for p in ps}):
        visit(x)
    idx = {x: i for i, x in enumerate(order)}
    A = np.zeros((len(order), len(order)))
    for x in order:
        i = idx[x]
        if x in parents:
            s, d = (idx[p] for p in parents[x])
            for j in range(i):
                A[i, j] = A[j, i] = 0.5 * (A[j, s] + A[j, d])
            A[i, i] = 1 + 0.5 * A[s, d]
        else:
            A[i, i] = 1.0
    return {"F": float(A[idx[individual], idx[individual]] - 1)}


def _rand_pedigree(rng):
    names = [f"I{k}" for k in range(1, 12)]
    nf = rng.randint(2, 3)
    pool = names[:nf]
    edges = []
    for k in range(nf, nf + rng.randint(3, 6)):
        x = names[k]
        recent = pool[-4:] if len(pool) >= 4 else pool
        a, b = rng.sample(recent, 2)
        edges += [(a, x), (b, x)]
        pool.append(x)
    return {"pedigree": edges, "individual": pool[-1]}


_add("bio_pedigree_inbreeding", name="Inbreeding coefficient from a pedigree", topic="genetics / inbreeding",
     purpose="Wright's inbreeding coefficient F of an individual from its pedigree (founders unrelated and not inbred).",
     inputs=[_edgesin("pedigree", "parent-child pairs, e.g. [P-S, Q-S, P-D, Q-D, S-X, D-X]"),
             _textin("individual", "name of the individual", "name", hi=40)],
     outputs=[_o("F", "inbreeding coefficient")], primary="F",
     conditions=["every non-founder has exactly two parents listed", "no cycles", "founders unrelated, F = 0"],
     method="Wright's path counting: sum over common ancestors and non-overlapping paths of (1/2)^(n1+n2+1) (1 + F_A)",
     latex=r"F_X = \sum_{\text{paths}} \left(\tfrac12\right)^{n_1+n_2+1}(1 + F_A)",
     run=pedigree_inbreeding, reference=_ref_pedigree,
     wrong=lambda pedigree, individual: {"F": 2 * pedigree_inbreeding(pedigree, individual)["F"]},
     rand=_rand_pedigree,
     phrases=["Pedigree inbreeding coefficient", "Inbreeding coefficient from a pedigree"],
     examples=[_ex({"pedigree": "P-S, Q-S, P-D, Q-D, S-X, D-X", "individual": "X"}, {"F": 0.25}),
               _ex({"pedigree": "A-C, B-C, A-D, C-D", "individual": "D"}, {"F": 0.25}),
               _ex({"pedigree": "P-S, Q-S, P-D, R-D, S-X, D-X", "individual": "X"}, {"F": 0.125})],
     extra_q=[{"q": "Pedigree inbreeding coefficient: pedigree = [P-S, S-X, D-X], individual = X",
               "outcome": "INVALID_INPUT"}],
     level="expert")


# ================================================================================================ PHYLOGENETICS
def _aligned_dna(seq1, seq2):
    a, b = _seq(seq1, "seq1", "DNA"), _seq(seq2, "seq2", "DNA")
    _need(len(a) == len(b), "the two aligned sequences must have the same length (no gaps)")
    return a, b


def jukes_cantor(seq1, seq2):
    a, b = _aligned_dna(seq1, seq2)
    p = sum(1 for x, y in zip(a, b) if x != y) / len(a)
    _need(p < 0.75, "p >= 0.75: the sequences are saturated and the Jukes-Cantor distance is undefined")
    return _out(p=p, d=-0.75 * math.log1p(-4.0 * p / 3.0))


def _ref_jc(seq1, seq2):
    """The JC69 rate matrix (total rate 1 per site): find t with P(different at t) = p from scipy's expm."""
    import numpy as np
    from scipy.linalg import expm
    from scipy.spatial.distance import hamming as h
    p = float(h(list(seq1.upper()), list(seq2.upper())))
    Q = np.full((4, 4), 1 / 3.0)
    np.fill_diagonal(Q, -1.0)
    if p == 0:
        return {"p": 0.0, "d": 0.0}
    return {"p": p, "d": _brentq(lambda t: 1.0 - float(expm(Q * t)[0, 0]) - p, 0.0, 30.0)}


def _rand_jc(rng):
    a = _rand_dna(rng, 20, 80)
    rate = rng.uniform(0.05, 0.45)
    b = "".join(c if rng.random() > rate else rng.choice([x for x in "ACGT" if x != c]) for c in a)
    return {"seq1": a, "seq2": b}


_add("bio_jukes_cantor_distance", name="Jukes-Cantor distance", topic="phylogenetics / substitution models",
     purpose="The Jukes-Cantor (JC69) evolutionary distance between two aligned DNA sequences: the expected number of "
             "substitutions per site, corrected for multiple hits.",
     inputs=[_textin("seq1", "first aligned DNA sequence"), _textin("seq2", "second aligned DNA sequence")],
     outputs=[_o("p", "proportion of differing sites"), _o("d", "Jukes-Cantor distance (substitutions per site)")],
     primary="d", conditions=["equal lengths, no gaps", "p < 0.75"],
     method="d = -(3/4) ln(1 - 4p/3)", latex=r"d = -\tfrac34 \ln\!\left(1 - \tfrac43 p\right)",
     run=jukes_cantor, reference=_ref_jc,
     wrong=lambda seq1, seq2: (lambda r: {"p": r["p"], "d": -math.log1p(-4 * r["p"] / 3)})(jukes_cantor(seq1, seq2)),
     rand=_rand_jc, phrases=["Jukes-Cantor distance", "JC69 distance"],
     examples=[_ex({"seq1": "AAAAAAAAAA", "seq2": "AAAAAAAACC"}, {"p": 0.2, "d": 0.232616196228}),
               _ex({"seq1": "ACGT", "seq2": "ACGT"}, {"p": 0.0, "d": 0.0})],
     level="intermediate")


_PUR = set("AG")


def kimura(seq1, seq2):
    a, b = _aligned_dna(seq1, seq2)
    n = len(a)
    ts = sum(1 for x, y in zip(a, b) if x != y and ((x in _PUR) == (y in _PUR)))
    tv = sum(1 for x, y in zip(a, b) if (x in _PUR) != (y in _PUR))
    P, Q = ts / n, tv / n
    _need(1 - 2 * P - Q > 0 and 1 - 2 * Q > 0, "too many differences: the Kimura distance is undefined")
    return _out(d=-0.5 * math.log(1 - 2 * P - Q) - 0.25 * math.log(1 - 2 * Q))


def _ref_kimura(seq1, seq2):
    """The K80 rate matrix exponential (scipy expm) solved for (alpha t, beta t) with scipy fsolve; d = (alpha + 2 beta) t."""
    import numpy as np
    from scipy.linalg import expm
    from scipy.optimize import fsolve
    a, b = seq1.upper(), seq2.upper()
    n = len(a)
    P = sum(1 for x, y in zip(a, b) if x != y and ((x in _PUR) == (y in _PUR))) / n
    Q = sum(1 for x, y in zip(a, b) if (x in _PUR) != (y in _PUR)) / n
    if P == 0 and Q == 0:
        return {"d": 0.0}

    def probs(v):
        al, be = v                                     # order A, G, C, T: A<->G and C<->T are transitions
        R = np.array([[0, al, be, be], [al, 0, be, be], [be, be, 0, al], [be, be, al, 0]], dtype=float)
        np.fill_diagonal(R, -R.sum(axis=1))
        M = expm(R)
        return M[0, 1], M[0, 2] + M[0, 3]

    sol, info, ier, msg = fsolve(lambda v: [probs(v)[0] - P, probs(v)[1] - Q], [max(P, 1e-3), max(Q / 2, 1e-3)],
                                 xtol=1e-15, full_output=True)
    if max(abs(x) for x in info["fvec"]) > 1e-14:
        raise RuntimeError("fsolve did not converge: " + msg)
    return {"d": float(sol[0] + 2 * sol[1])}


def _rand_kimura(rng):
    a = _rand_dna(rng, 30, 90)
    swap = {"A": "G", "G": "A", "C": "T", "T": "C"}
    out = []
    for c in a:
        u = rng.random()
        if u < 0.12:
            out.append(swap[c])
        elif u < 0.2:
            out.append(rng.choice(["C", "T"] if c in _PUR else ["A", "G"]))
        else:
            out.append(c)
    if out == list(a) or all((x in _PUR) == (y in _PUR) for x, y in zip(a, out)):
        out[0] = "C" if a[0] in _PUR else "A"
    return {"seq1": a, "seq2": "".join(out)}


_add("bio_kimura_distance", name="Kimura two-parameter distance", topic="phylogenetics / substitution models",
     purpose="The Kimura two-parameter (K80) distance between two aligned DNA sequences, separating transitions "
             "(A<->G, C<->T) from transversions.",
     inputs=[_textin("seq1", "first aligned DNA sequence"), _textin("seq2", "second aligned DNA sequence")],
     outputs=[_o("d", "K80 distance (substitutions per site)")], primary="d",
     conditions=["equal lengths, no gaps", "1 - 2P - Q > 0 and 1 - 2Q > 0"],
     method="P = transition proportion, Q = transversion proportion, d = -ln(1-2P-Q)/2 - ln(1-2Q)/4",
     latex=r"d = -\tfrac12\ln(1-2P-Q) - \tfrac14\ln(1-2Q)", run=kimura, reference=_ref_kimura,
     wrong=lambda seq1, seq2: {"d": jukes_cantor(seq1, seq2)["d"]}, rand=_rand_kimura,
     phrases=["Kimura distance", "Kimura two-parameter distance", "K80 distance"],
     examples=[_ex({"seq1": "AAAAAAAAAA", "seq2": "AAAAAAAGGC"}, {"d": 0.402359478109}),
               _ex({"seq1": "ACGTACGTAC", "seq2": "GCGTACGTAC"}, {"d": 0.111571775657})],
     level="advanced")


def _distance_matrix(D, name="D", lo=2, hi=30):
    n = len(D)
    _need(lo <= n <= hi and all(len(r) == n for r in D), f"{name} must be a square matrix of {lo} to {hi} taxa")
    for i in range(n):
        _need(D[i][i] == 0, f"{name}: the diagonal must be 0")
        for j in range(i + 1, n):
            _need(D[i][j] == D[j][i], f"{name} must be symmetric")
            _need(D[i][j] > 0, f"{name}: distances between different taxa must be positive")
    return [[float(x) for x in r] for r in D]


def _tie(vals, best, scale):
    return sum(1 for v in vals if abs(v - best) <= 1e-12 * scale) > 1


def upgma(D, linkage="average"):
    M = _distance_matrix(D)
    n = len(M)
    size = {i: 1 for i in range(n)}
    d = {(i, j): M[i][j] for i in range(n) for j in range(i + 1, n)}
    scale = max(d.values())
    heights, nxt = [], n
    while len(size) > 1:
        best = min(d.values())
        _need(not _tie(list(d.values()), best, scale), "two pairs tie at the closest distance: the tree is not unique")
        i, j = min(k for k, v in d.items() if v == best)
        heights.append(best / 2)
        new = nxt
        nxt += 1
        for k in list(size):
            if k in (i, j):
                continue
            dik, djk = d[tuple(sorted((i, k)))], d[tuple(sorted((j, k)))]
            if linkage == "average":
                d[(k, new)] = (size[i] * dik + size[j] * djk) / (size[i] + size[j])
            else:
                d[(k, new)] = min(dik, djk)
        size[new] = size.pop(i) + size.pop(j)
        d = {k: v for k, v in d.items() if i not in k and j not in k}
    return _out(root_height=heights[-1], heights=sorted(heights))


def _ref_upgma(D):
    import numpy as np
    from scipy.cluster.hierarchy import linkage
    from scipy.spatial.distance import squareform
    Z = linkage(squareform(np.array(D, dtype=float)), method="average")
    h = sorted(float(x) / 2 for x in Z[:, 2])
    return {"root_height": h[-1], "heights": h}


def _rand_dm(rng, n=None, additive=False):
    n = n or rng.randint(3, 8)
    if additive:
        nodes = {i: [] for i in range(n)}                       # a random tree, then path lengths plus noise
        active = list(range(n))
        nxt = n
        parent = {}
        while len(active) > 1:
            a, b = rng.sample(active, 2)
            parent[a] = (nxt, rng.uniform(0.5, 3))
            parent[b] = (nxt, rng.uniform(0.5, 3))
            active = [x for x in active if x not in (a, b)] + [nxt]
            nxt += 1

        def path(x):
            out, tot = {x: 0.0}, 0.0
            while x in parent:
                x, w = parent[x][0], parent[x][1]
                tot += w
                out[x] = tot
            return out
        P = [path(i) for i in range(n)]
        D = [[0.0] * n for _ in range(n)]
        for i in range(n):
            for j in range(i + 1, n):
                common = min((P[i][x] + P[j][x]) for x in P[i] if x in P[j])
                D[i][j] = D[j][i] = round(common * rng.uniform(0.95, 1.05), 3)
        return D
    D = [[0.0] * n for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            D[i][j] = D[j][i] = round(rng.uniform(1, 10), 3)
    return D


_add("bio_upgma", name="UPGMA tree heights", topic="phylogenetics / distance methods",
     purpose="The node heights of the UPGMA (average-linkage) tree of a distance matrix: the root height and every "
             "merge height (half the distance between the merged clusters).",
     inputs=[_matin("D", "symmetric distance matrix, zero diagonal")],
     outputs=[_o("root_height", "height of the root"), _o("heights", "all node heights, ascending")],
     primary="root_height", conditions=["symmetric, zero diagonal, positive distances", "no ties in a merge step"],
     method="repeatedly merge the closest clusters; distances to a merged cluster are size-weighted averages",
     latex=r"d_{(ij),k} = \frac{n_i d_{ik} + n_j d_{jk}}{n_i + n_j},\quad h = d_{ij}/2",
     run=upgma, reference=_ref_upgma, wrong=lambda D: upgma(D, linkage="single"), rand=lambda rng: {"D": _rand_dm(rng)},
     phrases=["UPGMA", "UPGMA tree heights"],
     examples=[_ex({"D": [[0, 2, 4, 6], [2, 0, 4, 6], [4, 4, 0, 6], [6, 6, 6, 0]]},
                   {"root_height": 3.0, "heights": [1.0, 2.0, 3.0]}),
               _ex({"D": [[0, 17, 21, 31, 23], [17, 0, 30, 34, 21], [21, 30, 0, 28, 39], [31, 34, 28, 0, 43],
                          [23, 21, 39, 43, 0]]}, {"root_height": 16.5, "heights": [8.5, 11.0, 14.0, 16.5]})],
     extra_q=[{"q": "UPGMA: D = [[0, 2], [3, 0]]", "outcome": "INVALID_INPUT"}],
     level="intermediate")


def neighbor_joining(D, corrected=True):
    M = _distance_matrix(D, lo=3)
    n = len(M)
    d = {(i, j): M[i][j] for i in range(n) for j in range(n) if i != j}
    nodes = list(range(n))
    branch = {}
    edges = []
    nxt = n
    scale = max(d.values())
    while len(nodes) > 2:
        r = len(nodes)
        R = {i: sum(d[(i, k)] for k in nodes if k != i) for i in nodes}
        Qv = {(i, j): (r - 2) * d[(i, j)] - R[i] - R[j] for a, i in enumerate(nodes) for j in nodes[a + 1:]}
        best = min(Qv.values())
        i, j = min(k for k, v in Qv.items() if v == best)
        tied = {k for k, v in Qv.items() if abs(v - best) <= 1e-12 * scale * r}
        # with 4 nodes left, a pair and its complement always tie and give the same tree; any other tie is refused
        harmless = {(i, j), tuple(sorted(k for k in nodes if k not in (i, j)))} if r == 4 else {(i, j)}
        _need(r == 3 or tied <= harmless, "two pairs tie in the neighbour-joining criterion: the tree is not unique")
        if corrected:
            li = d[(i, j)] / 2 + (R[i] - R[j]) / (2 * (r - 2))
        else:
            li = d[(i, j)] / 2
        lj = d[(i, j)] - li
        for x, lx in ((i, li), (j, lj)):
            edges.append(lx)
            if x < n:
                branch[x] = lx
        u = nxt
        nxt += 1
        for k in nodes:
            if k not in (i, j):
                d[(u, k)] = d[(k, u)] = (d[(i, k)] + d[(j, k)] - d[(i, j)]) / 2
        nodes = [k for k in nodes if k not in (i, j)] + [u]
    a, b = nodes
    edges.append(d[(a, b)])
    for x in (a, b):
        if x < n:
            branch[x] = d[(a, b)]
    return _out(total_length=sum(edges), branch=[branch[i] for i in range(n)])


def _ref_nj(D):
    from skbio import DistanceMatrix
    from skbio.tree import nj
    ids = [f"t{i}" for i in range(len(D))]
    tr = nj(DistanceMatrix([[float(x) for x in r] for r in D], ids), neg_as_zero=False)
    tips = {t.name: float(t.length) for t in tr.tips()}
    return {"total_length": float(sum(float(x.length) for x in tr.traverse() if x is not tr and x.length is not None)),
            "branch": [tips[i] for i in ids]}


_add("bio_neighbor_joining", name="Neighbour-joining branch lengths", topic="phylogenetics / distance methods",
     purpose="The total tree length and the terminal branch lengths of the neighbour-joining tree of a distance matrix.",
     inputs=[_matin("D", "symmetric distance matrix, zero diagonal (at least 3 taxa)")],
     outputs=[_o("total_length", "sum of all branch lengths"),
              _o("branch", "terminal branch length of each taxon, in input order")],
     primary="total_length", conditions=["at least 3 taxa", "no ties in the joining criterion",
                                         "negative branch lengths kept as computed"],
     method="Saitou-Nei: join the pair minimising Q = (r-2) d_ij - R_i - R_j; branch lengths from the row sums",
     latex=r"Q_{ij} = (r-2)d_{ij} - R_i - R_j,\quad \ell_i = \tfrac12 d_{ij} + \frac{R_i - R_j}{2(r-2)}",
     run=neighbor_joining, reference=_ref_nj, wrong=lambda D: neighbor_joining(D, corrected=False),
     rand=lambda rng: {"D": _rand_dm(rng, rng.randint(4, 8), additive=True)},
     phrases=["Neighbour joining", "Neighbor joining", "Neighbour-joining branch lengths"],
     examples=[_ex({"D": [[0, 5, 9, 9, 8], [5, 0, 10, 10, 9], [9, 10, 0, 8, 7], [9, 10, 8, 0, 3], [8, 9, 7, 3, 0]]},
                   {"total_length": 17.0, "branch": [2.0, 3.0, 4.0, 2.0, 1.0]}),
               _ex({"D": [[0, 3, 4], [3, 0, 5], [4, 5, 0]]}, {"total_length": 6.0, "branch": [1.0, 2.0, 3.0]})],
     level="advanced")


def _dfact(k):
    out = 1
    while k > 1:
        out *= k
        k -= 2
    return out


def tree_count(n, rooted):
    _need(rooted in ("rooted", "unrooted"), "rooted must be rooted or unrooted")
    lo = 2 if rooted == "rooted" else 3
    _need(isinstance(n, int) and lo <= n <= 100, f"n must be a whole number of taxa from {lo} to 100")
    return {"count": _dfact(2 * n - 3) if rooted == "rooted" else _dfact(2 * n - 5)}


def _ref_tree_count(n, rooted):
    """Explicit enumeration: rooted trees by recursive bipartitions of the taxon set; unrooted trees by attaching taxa
    to every edge, deduplicated by their split sets."""
    taxa = list(range(n))
    if rooted == "rooted":
        from functools import lru_cache

        @lru_cache(maxsize=None)
        def trees(s):
            s = tuple(s)
            if len(s) == 1:
                return [s[0]]
            out = set()
            first, rest = s[0], s[1:]
            for k in range(len(rest) + 1):
                for c in combinations(rest, k):
                    left = (first,) + c
                    right = tuple(x for x in rest if x not in c)
                    if not right:
                        continue
                    for a in trees(left):
                        for b in trees(right):
                            out.add(frozenset([a, b]))
            return list(out)
        return {"count": len(trees(tuple(taxa)))}
    seen = set()

    def splits(edges):
        adj = {}
        for u, v in edges:
            adj.setdefault(u, []).append(v)
            adj.setdefault(v, []).append(u)
        res = set()
        for u, v in edges:
            stack, side = [v], {v}
            while stack:
                x = stack.pop()
                for y in adj[x]:
                    if y != u and y not in side:
                        side.add(y)
                        stack.append(y)
            leaves = frozenset(x for x in side if isinstance(x, int))
            if 0 in leaves:
                leaves = frozenset(taxa) - leaves
            if 1 < len(leaves) < n - 1:
                res.add(leaves)
        return frozenset(res)
    level = [[(0, "c"), (1, "c"), (2, "c")]]
    for k in range(3, n):
        nxt = []
        for edges in level:
            for idx, (u, v) in enumerate(edges):
                w = f"w{k}"
                nxt.append(edges[:idx] + edges[idx + 1:] + [(u, w), (w, v), (k, w)])
        level = nxt
    for edges in level:
        seen.add(splits(edges))
    return {"count": len(seen)}


_add("bio_tree_count", name="Number of possible phylogenetic trees", topic="phylogenetics / tree space",
     purpose="The number of distinct labelled bifurcating trees for n taxa, rooted or unrooted.",
     inputs=[_intin("n", "number of taxa", lo=2, hi=100), _choice("rooted", ["rooted", "unrooted"], "rooted or unrooted trees")],
     outputs=[_o("count", "number of distinct bifurcating trees")], primary="count",
     conditions=["n >= 2 (rooted) or n >= 3 (unrooted)", "n <= 100"],
     method="double factorials: (2n-3)!! rooted, (2n-5)!! unrooted",
     latex=r"N_{\text{rooted}} = (2n-3)!!,\quad N_{\text{unrooted}} = (2n-5)!!", run=tree_count,
     reference=_ref_tree_count,
     wrong=lambda n, rooted: {"count": _dfact(2 * n - 5) if rooted == "rooted" else _dfact(2 * n - 3)},
     rand=lambda rng: (lambda r: {"n": rng.randint(4 if r == "unrooted" else 3, 7), "rooted": r})(
         rng.choice(["rooted", "unrooted"])),
     phrases=["Number of trees", "Number of phylogenetic trees"],
     examples=[_ex({"n": 5, "rooted": "unrooted"}, {"count": 15}), _ex({"n": 7, "rooted": "rooted"}, {"count": 10395})],
     level="intermediate")


def _parse_newick(text):
    """A Newick string -> (name, length, children) nodes; whitespace ignored; leaf names required and unique."""
    _need(isinstance(text, str), "the tree must be a Newick string")
    s = re.sub(r"\s+", "", text).rstrip(";")
    pos = [0]

    def num():
        m = re.match(r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?", s[pos[0]:])
        _need(m is not None, "a branch length must be a number")
        pos[0] += m.end()
        return float(m.group(0))

    def node():
        kids = []
        if pos[0] < len(s) and s[pos[0]] == "(":
            pos[0] += 1
            kids.append(node())
            while pos[0] < len(s) and s[pos[0]] == ",":
                pos[0] += 1
                kids.append(node())
            _need(pos[0] < len(s) and s[pos[0]] == ")", "unbalanced parentheses in the tree")
            pos[0] += 1
        m = re.match(r"[A-Za-z0-9_.]*", s[pos[0]:])
        name = m.group(0)
        pos[0] += m.end()
        length = None
        if pos[0] < len(s) and s[pos[0]] == ":":
            pos[0] += 1
            length = num()
        if not kids:
            _need(re.fullmatch(r"[A-Za-z][A-Za-z0-9_]*", name or "") is not None, "every leaf needs a name")
        return (name, length, kids)
    root = node()
    _need(pos[0] == len(s), "unreadable text after the tree")
    leaves = []

    def walk(t):
        if not t[2]:
            leaves.append(t[0])
        for k in t[2]:
            walk(k)
    walk(root)
    _need(len(set(leaves)) == len(leaves), "leaf names must be unique")
    return root, leaves


def _splits(root, leaves):
    allx = frozenset(leaves)
    first = min(leaves)
    out = set()

    def walk(t):
        if not t[2]:
            return frozenset([t[0]])
        s = frozenset().union(*[walk(k) for k in t[2]])
        side = allx - s if first in s else s
        if 2 <= len(side) <= len(allx) - 2:
            out.add(side)
        return s
    walk(root)
    return out


def robinson_foulds(tree1, tree2, one_way=False):
    r1, l1 = _parse_newick(tree1)
    r2, l2 = _parse_newick(tree2)
    _need(set(l1) == set(l2), "both trees must have the same taxa")
    _need(len(l1) >= 4, "at least 4 taxa are needed")
    s1, s2 = _splits(r1, l1), _splits(r2, l2)
    _need(s1 or s2, "both trees are unresolved stars")
    rf = len(s1 - s2) + (0 if one_way else len(s2 - s1))
    return {"rf": rf, "rf_normalized": rf / (len(s1) + len(s2))}


def _ref_rf(tree1, tree2):
    import io
    from skbio import TreeNode
    t1, t2 = TreeNode.read(io.StringIO(tree1.strip().rstrip(";") + ";")), TreeNode.read(io.StringIO(tree2.strip().rstrip(";") + ";"))
    return {"rf": int(round(t1.compare_rfd(t2, rooted=False))),
            "rf_normalized": float(t1.compare_rfd(t2, proportion=True, rooted=False))}


def _rand_newick(rng, taxa, lengths=False):
    nodes = [t for t in taxa]
    while len(nodes) > 1:
        a, b = rng.sample(nodes, 2)
        nodes = [x for x in nodes if x not in (a, b)]
        if lengths:
            nodes.append(f"({a}:{round(rng.uniform(0.1, 2), 2)},{b}:{round(rng.uniform(0.1, 2), 2)})")
        else:
            nodes.append(f"({a},{b})")
    return nodes[0] + ";"


def _rand_rf(rng):
    taxa = [f"T{i}" for i in range(rng.randint(5, 9))]
    return {"tree1": _rand_newick(rng, taxa), "tree2": _rand_newick(rng, taxa)}


_add("bio_robinson_foulds", name="Robinson-Foulds distance between trees", topic="phylogenetics / tree comparison",
     purpose="The Robinson-Foulds (symmetric-difference) distance between two trees on the same taxa, comparing their "
             "bipartitions (trees treated as unrooted), and the distance as a proportion.",
     inputs=[_textin("tree1", "first tree in Newick format", "newick", hi=20000),
             _textin("tree2", "second tree in Newick format", "newick", hi=20000)],
     outputs=[_o("rf", "number of bipartitions found in only one tree"),
              _o("rf_normalized", "RF divided by the total number of non-trivial bipartitions of both trees")],
     primary="rf", conditions=["the same taxa", "at least 4 taxa", "trees compared unrooted; branch lengths ignored"],
     method="non-trivial bipartitions of each tree; count those not shared",
     latex=r"\mathrm{RF} = |S_1 \setminus S_2| + |S_2 \setminus S_1|", run=robinson_foulds, reference=_ref_rf,
     wrong=lambda tree1, tree2: robinson_foulds(tree1, tree2, one_way=True), rand=_rand_rf,
     phrases=["Robinson-Foulds distance", "RF distance between trees"],
     examples=[_ex({"tree1": "((A,B),(C,D),(E,F));", "tree2": "((A,C),(B,D),(E,F));"}, {"rf": 4, "rf_normalized": 4 / 6}),
               _ex({"tree1": "(((A,B),C),(D,E));", "tree2": "((A,B),(C,(D,E)));"}, {"rf": 0, "rf_normalized": 0.0})],
     level="advanced")


def _states(text, need):
    _need(isinstance(text, str), "states must be written as name:SEQUENCE pairs")
    out = {}
    for tok in text.split():
        m = re.fullmatch(r"([A-Za-z][A-Za-z0-9_]*):([A-Za-z]+)", tok)
        _need(m is not None, f"cannot read '{tok}': write name:SEQUENCE")
        _need(m.group(1) not in out, f"{m.group(1)} is given twice")
        out[m.group(1)] = m.group(2).upper()
    _need(set(out) == set(need), "give one sequence for every leaf of the tree, and no others")
    lens = {len(v) for v in out.values()}
    _need(len(lens) == 1, "all sequences must have the same length (aligned)")
    return out


def fitch_parsimony(tree, states, lower_bound=False):
    root, leaves = _parse_newick(tree)

    def binary(t, top):
        if t[2]:
            _need(len(t[2]) == 2, "the tree must be rooted and strictly bifurcating (two children per internal node)")
            for k in t[2]:
                binary(k, False)
    binary(root, True)
    seqs = _states(states, leaves)
    L = len(next(iter(seqs.values())))
    total = 0
    for site in range(L):
        col = {n: s[site] for n, s in seqs.items()}
        alphabet = sorted(set(col.values()))
        if lower_bound:
            total += len(alphabet) - 1
            continue

        def cost(t):                                   # Sankoff with unit costs
            if not t[2]:
                return {a: (0 if a == col[t[0]] else math.inf) for a in alphabet}
            cs = [cost(k) for k in t[2]]
            return {a: sum(min(c[b] + (a != b) for b in alphabet) for c in cs) for a in alphabet}
        total += int(min(cost(root).values()))
    return {"score": total}


def _ref_fitch(tree, states):
    import io
    from Bio import Phylo
    from Bio.Align import MultipleSeqAlignment
    from Bio.Phylo.TreeConstruction import ParsimonyScorer
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    t = Phylo.read(io.StringIO(tree.strip().rstrip(";") + ";"), "newick")
    t.rooted = True
    recs = [SeqRecord(Seq(s.split(":")[1].upper()), id=s.split(":")[0]) for s in states.split()]
    return {"score": int(ParsimonyScorer().get_score(t, MultipleSeqAlignment(recs)))}


def _rand_fitch(rng):
    taxa = [f"T{i}" for i in range(rng.randint(4, 8))]
    L = rng.randint(4, 12)
    base = _rand_dna(rng, L, L)
    seqs = {t: "".join(c if rng.random() > 0.35 else rng.choice("ACGT") for c in base) for t in taxa}
    return {"tree": _rand_newick(rng, taxa), "states": " ".join(f"{t}:{s}" for t, s in seqs.items())}


_add("bio_fitch_parsimony", name="Parsimony score of a tree (Fitch)", topic="phylogenetics / parsimony",
     purpose="The minimum number of character changes (small parsimony score) needed on a given rooted bifurcating "
             "tree to explain aligned sequences at its leaves.",
     inputs=[_textin("tree", "rooted bifurcating tree in Newick format", "newick", hi=20000),
             _textin("states", "aligned sequence of each leaf, written name:SEQUENCE separated by spaces", "states", hi=100000)],
     outputs=[_o("score", "parsimony score (number of changes)")], primary="score",
     conditions=["strictly bifurcating rooted tree", "one aligned sequence per leaf", "unit cost per change"],
     method="dynamic programming over the tree site by site (equivalent to Fitch's algorithm with unit costs)",
     latex=r"c_v(a) = \sum_{w\in\mathrm{children}(v)} \min_b \big(c_w(b) + [a\ne b]\big)", run=fitch_parsimony,
     reference=_ref_fitch, wrong=lambda tree, states: fitch_parsimony(tree, states, lower_bound=True), rand=_rand_fitch,
     phrases=["Parsimony score", "Fitch parsimony score"],
     examples=[_ex({"tree": "((A,B),(C,D));", "states": "A:AACG B:AATG C:GACT D:GGCT"}, {"score": 4}),
               _ex({"tree": "((A,C),(B,D));", "states": "A:A B:A C:G D:G"}, {"score": 2})],
     extra_q=[{"q": "Parsimony score: tree = (A,B,C);, states = A:A B:C C:G", "outcome": "INVALID_INPUT"}],
     level="expert")


def coalescent_times(n, Ne):
    _need(isinstance(n, int) and 2 <= n <= 100000, "n must be a whole number of sampled gene copies (2 to 100000)")
    _pos("Ne", Ne)
    a = math.fsum(1.0 / i for i in range(1, n))
    return _out(t_mrca=4 * Ne * (1 - 1 / n), total_length=4 * Ne * a)


def _ref_coalescent(n, Ne):
    """Absorbing continuous-time Markov chain over the number of lineages (k -> k-1 at rate C(k,2)/(2Ne)): expected
    absorption time and expected total branch length (reward k per unit time) from a numpy linear solve."""
    import numpy as np
    states = list(range(n, 1, -1))
    m = len(states)
    Qm = np.zeros((m, m))
    for a, k in enumerate(states):
        rate = k * (k - 1) / 2 / (2 * Ne)
        Qm[a, a] = -rate
        if a + 1 < m:
            Qm[a, a + 1] = rate
    T = np.linalg.solve(-Qm, np.ones(m))
    Lr = np.linalg.solve(-Qm, np.array(states, dtype=float))
    return {"t_mrca": float(T[0]), "total_length": float(Lr[0])}


_add("bio_coalescent_times", name="Expected coalescence times", topic="population genetics / coalescent",
     purpose="The expected time to the most recent common ancestor and the expected total branch length of the "
             "genealogy of n gene copies from a diploid population of effective size Ne (in generations).",
     inputs=[_intin("n", "number of sampled gene copies", lo=2), _numin("Ne", "diploid effective population size")],
     outputs=[_o("t_mrca", "expected TMRCA in generations"), _o("total_length", "expected total branch length in generations")],
     primary="t_mrca", conditions=["neutral Kingman coalescent", "diploid: a pair coalesces at rate 1/(2Ne) per generation"],
     method="sum of exponential waiting times while k lineages remain, mean 4Ne/(k(k-1))",
     latex=r"E[T_{MRCA}] = 4N_e\left(1 - \frac1n\right),\quad E[L] = 4N_e\sum_{i=1}^{n-1}\frac1i",
     run=coalescent_times, reference=_ref_coalescent,
     wrong=lambda n, Ne: {"t_mrca": 2 * Ne * (1 - 1 / n), "total_length": 2 * Ne * sum(1 / i for i in range(1, n))},
     rand=lambda rng: {"n": rng.randint(2, 60), "Ne": float(rng.choice([100, 500, 1000, 2500, 10000]) * rng.randint(1, 9))},
     phrases=["Coalescent times", "Expected coalescence times"],
     examples=[_ex({"n": 2, "Ne": 1000}, {"t_mrca": 2000.0, "total_length": 4000.0}),
               _ex({"n": 10, "Ne": 10000}, {"t_mrca": 36000.0, "total_length": 113158.730159})],
     level="advanced")


def _alignment(seqs, lo=2):
    _need(isinstance(seqs, str), "give the aligned sequences separated by spaces")
    rows = [s.upper() for s in seqs.split()]
    _need(len(rows) >= lo, f"at least {lo} sequences are needed")
    for r in rows:
        _need(re.fullmatch(r"[ACGT]+", r) is not None, "sequences must be DNA (A, C, G, T) without gaps")
    _need(len({len(r) for r in rows}) == 1, "the sequences must be aligned (equal lengths)")
    return rows


def _pairwise_k(rows):
    n = len(rows)
    tot = sum(sum(1 for x, y in zip(rows[i], rows[j]) if x != y) for i in range(n) for j in range(i + 1, n))
    return Fraction(tot, n * (n - 1) // 2)


def nucleotide_diversity(seqs):
    rows = _alignment(seqs)
    k = _pairwise_k(rows)
    return {"k": float(k), "pi": float(k / len(rows[0]))}


def _allele_counts(rows):
    import numpy as np
    idx = {"A": 0, "C": 1, "G": 2, "T": 3}
    ac = np.zeros((len(rows[0]), 4), dtype=int)
    for r in rows:
        for j, c in enumerate(r):
            ac[j, idx[c]] += 1
    return ac


def _ref_diversity(seqs):
    import allel
    import numpy as np
    rows = [s.upper() for s in seqs.split()]
    ac = _allele_counts(rows)
    L = len(rows[0])
    pos = np.arange(1, L + 1)
    return {"k": float(np.sum(allel.mean_pairwise_difference(ac, fill=0))),
            "pi": float(allel.sequence_diversity(pos, ac, start=1, stop=L))}


def _wrong_diversity(seqs):
    rows = [s.upper() for s in seqs.split()]
    k = _pairwise_k(rows)
    S = sum(1 for col in zip(*rows) if len(set(col)) > 1)
    return {"k": float(k), "pi": float(k / max(S, 1))}               # per segregating site instead of per site


def _rand_alignment(rng, min_seg=1):
    n = rng.randint(3, 10)
    L = rng.randint(20, 60)
    base = _rand_dna(rng, L, L)
    while True:
        rows = ["".join(c if rng.random() > 0.08 else rng.choice("ACGT") for c in base) for _ in range(n)]
        if sum(1 for col in zip(*rows) if len(set(col)) > 1) >= min_seg:
            return {"seqs": " ".join(rows)}


_add("bio_nucleotide_diversity", name="Nucleotide diversity", topic="population genetics / sequence diversity",
     purpose="The average number of pairwise differences k and the nucleotide diversity pi (per site) of a set of "
             "aligned DNA sequences.",
     inputs=[_textin("seqs", "aligned DNA sequences separated by spaces", "seqs", hi=200000)],
     outputs=[_o("k", "mean number of pairwise differences"), _o("pi", "nucleotide diversity per site")],
     primary="pi", conditions=["at least 2 aligned DNA sequences of equal length, no gaps"],
     method="average the differences over all pairs of sequences; divide by the alignment length",
     latex=r"\pi = \frac{1}{L}\binom{n}{2}^{-1}\sum_{i<j} d_{ij}", run=nucleotide_diversity, reference=_ref_diversity,
     wrong=_wrong_diversity, rand=_rand_alignment,
     phrases=["Nucleotide diversity", "Nucleotide diversity pi"],
     examples=[_ex({"seqs": "AAAA AAAT AATT"}, {"k": 4 / 3, "pi": 1 / 3}),
               _ex({"seqs": "ACGTACGT ACGTACGA ACGAACGA ACGTACGT"}, {"k": 7 / 6, "pi": 7 / 48})],
     level="intermediate")


def watterson(seqs, harmonic_to_n=False):
    rows = _alignment(seqs)
    n, L = len(rows), len(rows[0])
    S = sum(1 for col in zip(*rows) if len(set(col)) > 1)
    a = math.fsum(1.0 / i for i in range(1, n + (1 if harmonic_to_n else 0)))
    return _out(S=S, theta=S / a, theta_site=S / a / L)


def _ref_watterson(seqs):
    import allel
    import numpy as np
    rows = [s.upper() for s in seqs.split()]
    ac = allel.AlleleCountsArray(_allele_counts(rows))
    L = len(rows[0])
    th = float(allel.watterson_theta(np.arange(1, L + 1), ac, start=1, stop=L))
    return {"S": int(ac.count_segregating()), "theta": th * L, "theta_site": th}


_add("bio_watterson_theta", name="Watterson's estimator of theta", topic="population genetics / sequence diversity",
     purpose="The number of segregating sites S and Watterson's estimator of the population mutation rate theta, per "
             "sequence and per site.",
     inputs=[_textin("seqs", "aligned DNA sequences separated by spaces", "seqs", hi=200000)],
     outputs=[_o("S", "number of segregating sites"), _o("theta", "Watterson's theta per sequence (S / a_n)"),
              _o("theta_site", "Watterson's theta per site")],
     primary="theta", conditions=["at least 2 aligned DNA sequences of equal length, no gaps"],
     method="count segregating sites and divide by the harmonic number a_n = sum_{i=1}^{n-1} 1/i",
     latex=r"\hat\theta_W = \frac{S}{a_n},\quad a_n = \sum_{i=1}^{n-1}\frac1i", run=watterson, reference=_ref_watterson,
     wrong=lambda seqs: watterson(seqs, harmonic_to_n=True), rand=_rand_alignment,
     phrases=["Watterson theta", "Watterson's theta"],
     examples=[_ex({"seqs": "AAAA AAAT AATT"}, {"S": 2, "theta": 4 / 3, "theta_site": 1 / 3}),
               _ex({"seqs": "ACGT ACGA"}, {"S": 1, "theta": 1.0, "theta_site": 0.25})],
     level="advanced")


def tajima_d(seqs, sqrt_denominator=True):
    rows = _alignment(seqs, lo=4)
    n = len(rows)
    S = sum(1 for col in zip(*rows) if len(set(col)) > 1)
    _need(S >= 1, "Tajima's D needs at least one segregating site")
    k = float(_pairwise_k(rows))
    a1 = math.fsum(1.0 / i for i in range(1, n))
    a2 = math.fsum(1.0 / (i * i) for i in range(1, n))
    b1 = (n + 1) / (3 * (n - 1))
    b2 = 2 * (n * n + n + 3) / (9 * n * (n - 1))
    c1 = b1 - 1 / a1
    c2 = b2 - (n + 2) / (a1 * n) + a2 / (a1 * a1)
    e1, e2 = c1 / a1, c2 / (a1 * a1 + a2)
    var = e1 * S + e2 * S * (S - 1)
    return _out(D=(k - S / a1) / (math.sqrt(var) if sqrt_denominator else var), k=k)


def _ref_tajima(seqs):
    import allel
    import numpy as np
    rows = [s.upper() for s in seqs.split()]
    ac = _allele_counts(rows)
    return {"D": float(allel.tajima_d(ac, min_sites=1)), "k": float(np.sum(allel.mean_pairwise_difference(ac, fill=0)))}


_add("bio_tajima_d", name="Tajima's D", topic="population genetics / neutrality tests",
     purpose="Tajima's D neutrality statistic of aligned DNA sequences (the scaled difference between pi and "
             "Watterson's theta).",
     inputs=[_textin("seqs", "aligned DNA sequences separated by spaces (at least 4)", "seqs", hi=200000)],
     outputs=[_o("D", "Tajima's D"), _o("k", "mean number of pairwise differences")], primary="D",
     conditions=["at least 4 aligned DNA sequences", "at least one segregating site"],
     method="D = (k - S/a1) / sqrt(e1 S + e2 S (S-1)) with Tajima's constants",
     latex=r"D = \frac{\hat\theta_\pi - \hat\theta_W}{\sqrt{e_1S + e_2S(S-1)}}", run=tajima_d, reference=_ref_tajima,
     wrong=lambda seqs: tajima_d(seqs, sqrt_denominator=False),
     rand=lambda rng: _rand_alignment(rng, min_seg=3),
     phrases=["Tajima's D", "Tajima D"],
     examples=[_ex({"seqs": "AAAAAA AAAAAT AAAATT AAATTT"}, {"D": 0.167655795034, "k": 5 / 3}),
               _ex({"seqs": "ACGTACGTAA ACGTACGTAA ACGTACGTAA TCGTACGTAA"}, {"D": -0.612372435696})],
     level="expert")


# ================================================================================================ ENZYME KINETICS
def _kin_data(S, v, lo=3):
    _need(isinstance(S, list) and isinstance(v, list) and len(S) == len(v) and len(S) >= lo,
          f"S and v must be lists of equal length (at least {lo} points)")
    for x in S + v:
        _pos("every concentration and rate", x)
    return [float(x) for x in S], [float(x) for x in v]


def _ols(x, y):
    n = len(x)
    mx, my = math.fsum(x) / n, math.fsum(y) / n
    sxx = math.fsum((a - mx) ** 2 for a in x)
    syy = math.fsum((b - my) ** 2 for b in y)
    sxy = math.fsum((a - mx) * (b - my) for a, b in zip(x, y))
    _need(sxx > 0, "the x values must not all be equal")
    b = sxy / sxx
    r2 = sxy * sxy / (sxx * syy) if syy > 0 else 1.0
    return b, my - b * mx, r2


def _linregress(x, y):
    from scipy.stats import linregress
    r = linregress(x, y)
    return float(r.slope), float(r.intercept), float(r.rvalue) ** 2


def lineweaver_burk(S, v):
    S, v = _kin_data(S, v, lo=2)
    b, a, r2 = _ols([1 / x for x in S], [1 / y for y in v])
    _need(a > 0 and b > 0, "the double-reciprocal line must have a positive slope and intercept")
    return _out(Vmax=1 / a, Km=b / a, r2=r2)


def _ref_lb(S, v):
    b, a, r2 = _linregress([1 / x for x in S], [1 / y for y in v])
    return {"Vmax": 1 / a, "Km": b / a, "r2": r2}


def _rand_mm_data(rng, n=None):
    Vm, Km = rng.uniform(5, 100), rng.uniform(0.5, 50)
    mult = sorted(rng.sample([0.2, 0.3, 0.5, 0.75, 1, 1.5, 2, 3, 5, 8, 10], n or rng.randint(5, 8)))
    S = [_sig(Km * m * rng.uniform(0.95, 1.05)) for m in mult]
    v = [_sig(Vm * s / (Km + s) * (1 + rng.uniform(-0.04, 0.04))) for s in S]
    return {"S": S, "v": v}


_KIN_IN = [_vecin("S", "substrate concentrations (any one unit)", lo=2), _vecin("v", "initial rates (any one unit)", lo=2)]
_KIN_OUT = [_o("Vmax", "maximum rate (unit of v)"), _o("Km", "Michaelis constant (unit of S)"),
            _o("r2", "coefficient of determination of the linear fit")]

_add("bio_lineweaver_burk", name="Lineweaver-Burk fit of Michaelis-Menten data", topic="enzyme kinetics / linearised fits",
     purpose="Vmax and Km from a least-squares line through the double-reciprocal (Lineweaver-Burk) plot of 1/v "
             "against 1/S.",
     inputs=_KIN_IN, outputs=_KIN_OUT, primary="Km", conditions=["positive S and v", "positive slope and intercept"],
     method="ordinary least squares of 1/v on 1/S: intercept 1/Vmax, slope Km/Vmax",
     latex=r"\frac1v = \frac{K_m}{V_{\max}}\frac1S + \frac1{V_{\max}}", run=lineweaver_burk, reference=_ref_lb,
     wrong=lambda S, v: (lambda r: dict(r, Km=r["Km"] / r["Vmax"]))(lineweaver_burk(S, v)), rand=_rand_mm_data,
     phrases=["Lineweaver-Burk fit", "Lineweaver-Burk plot"],
     examples=[_ex({"S": [1, 2, 4, 8], "v": [5, 6.6666666667, 8, 8.8888888889]}, {"Vmax": 10.0, "Km": 1.0}),
               _ex({"S": [0.5, 1, 2, 5], "v": [2.2, 3.5, 5.1, 6.6]}, {"Vmax": 8.69516077399, "Km": 1.47542149819})],
     level="intermediate")


def hanes_woolf(S, v):
    S, v = _kin_data(S, v, lo=2)
    b, a, r2 = _ols(S, [x / y for x, y in zip(S, v)])
    _need(a > 0 and b > 0, "the Hanes-Woolf line must have a positive slope and intercept")
    return _out(Vmax=1 / b, Km=a / b, r2=r2)


def _ref_hw(S, v):
    b, a, r2 = _linregress(S, [x / y for x, y in zip(S, v)])
    return {"Vmax": 1 / b, "Km": a / b, "r2": r2}


_add("bio_hanes_woolf", name="Hanes-Woolf fit of Michaelis-Menten data", topic="enzyme kinetics / linearised fits",
     purpose="Vmax and Km from a least-squares line through the Hanes-Woolf plot of S/v against S.",
     inputs=_KIN_IN, outputs=_KIN_OUT, primary="Km", conditions=["positive S and v", "positive slope and intercept"],
     method="ordinary least squares of S/v on S: slope 1/Vmax, intercept Km/Vmax",
     latex=r"\frac Sv = \frac{1}{V_{\max}}S + \frac{K_m}{V_{\max}}", run=hanes_woolf, reference=_ref_hw,
     wrong=lambda S, v: (lambda r: dict(r, Km=r["Km"] / r["Vmax"]))(hanes_woolf(S, v)), rand=_rand_mm_data,
     phrases=["Hanes-Woolf fit", "Hanes-Woolf plot"],
     examples=[_ex({"S": [1, 2, 4, 8], "v": [5, 6.6666666667, 8, 8.8888888889]}, {"Vmax": 10.0, "Km": 1.0}),
               _ex({"S": [0.5, 1, 2, 5], "v": [2.2, 3.5, 5.1, 6.6]}, {"Vmax": 8.47375759387, "Km": 1.3973707798})],
     level="intermediate")


def eadie_hofstee(S, v):
    S, v = _kin_data(S, v, lo=2)
    b, a, r2 = _ols([y / x for x, y in zip(S, v)], v)
    _need(b < 0 and a > 0, "the Eadie-Hofstee line must have a negative slope and a positive intercept")
    return _out(Vmax=a, Km=-b, r2=r2)


def _ref_eh(S, v):
    b, a, r2 = _linregress([y / x for x, y in zip(S, v)], v)
    return {"Vmax": a, "Km": -b, "r2": r2}


_add("bio_eadie_hofstee", name="Eadie-Hofstee fit of Michaelis-Menten data", topic="enzyme kinetics / linearised fits",
     purpose="Vmax and Km from a least-squares line through the Eadie-Hofstee plot of v against v/S.",
     inputs=_KIN_IN, outputs=_KIN_OUT, primary="Km", conditions=["positive S and v", "negative slope, positive intercept"],
     method="ordinary least squares of v on v/S: intercept Vmax, slope -Km",
     latex=r"v = V_{\max} - K_m\frac{v}{S}", run=eadie_hofstee, reference=_ref_eh,
     wrong=lambda S, v: (lambda r: dict(r, Km=r["Km"] / r["Vmax"]))(eadie_hofstee(S, v)), rand=_rand_mm_data,
     phrases=["Eadie-Hofstee fit", "Eadie-Hofstee plot"],
     examples=[_ex({"S": [1, 2, 4, 8], "v": [5, 6.6666666667, 8, 8.8888888889]}, {"Vmax": 10.0, "Km": 1.0}),
               _ex({"S": [0.5, 1, 2, 5], "v": [2.2, 3.5, 5.1, 6.6]}, {"Vmax": 8.60595833712, "Km": 1.44637496589})],
     level="intermediate")


def _mm_varpro(S, v, K):
    x = [s / (K + s) for s in S]
    xd = [-s / (K + s) ** 2 for s in S]
    svx = math.fsum(a * b for a, b in zip(v, x))
    sxx = math.fsum(a * a for a in x)
    svxd = math.fsum(a * b for a, b in zip(v, xd))
    sxxd = math.fsum(a * b for a, b in zip(x, xd))
    return svx, sxx, svxd * sxx - svx * sxxd


def mm_fit(S, v):
    S, v = _kin_data(S, v)
    _need(len(set(S)) >= 3, "at least 3 different substrate concentrations are needed")
    lo, hi = min(S) * 1e-3, max(S) * 1e3
    grid = [lo * (hi / lo) ** (i / 300) for i in range(301)]
    F = []
    for K in grid:
        svx, sxx, _ = _mm_varpro(S, v, K)
        F.append(svx * svx / sxx)
    i = max(range(len(F)), key=lambda k: F[k])
    _need(0 < i < len(F) - 1, "the data do not determine Km (the best fit lies at the edge of the search range)")
    K = _bisect(lambda k: _mm_varpro(S, v, k)[2], grid[i - 1], grid[i + 1])
    svx, sxx, _ = _mm_varpro(S, v, K)
    return _out(Vmax=svx / sxx, Km=K)


def _ref_mm_fit(S, v):
    import numpy as np
    from scipy.optimize import curve_fit
    S, v = np.array(S, float), np.array(v, float)
    p, _ = curve_fit(lambda s, Vm, Km: Vm * s / (Km + s), S, v,
                     p0=[max(v), float(np.median(S))], method="lm", xtol=1e-15, ftol=1e-15, gtol=1e-15, maxfev=100000)
    p = _gn_polish(p, lambda q: (q[0] * S / (q[1] + S), np.column_stack([S / (q[1] + S), -q[0] * S / (q[1] + S) ** 2])), v)
    return {"Vmax": float(p[0]), "Km": float(p[1])}


_add("bio_michaelis_menten_fit", name="Michaelis-Menten nonlinear least-squares fit", topic="enzyme kinetics / Michaelis-Menten",
     purpose="Vmax and Km of the Michaelis-Menten equation fitted to rate data by nonlinear least squares (unweighted).",
     inputs=[_vecin("S", "substrate concentrations (any one unit)", lo=3), _vecin("v", "initial rates (any one unit)", lo=3)],
     outputs=[_o("Vmax", "maximum rate (unit of v)"), _o("Km", "Michaelis constant (unit of S)")], primary="Km",
     conditions=["at least 3 distinct positive S", "positive v", "an interior least-squares optimum"],
     method="variable projection: Vmax is linear given Km; the best Km zeroes the derivative of the projected fit",
     latex=r"\min_{V_{\max},K_m} \sum_i \Big(v_i - \frac{V_{\max}S_i}{K_m+S_i}\Big)^2", run=mm_fit,
     reference=_ref_mm_fit, wrong=lambda S, v: {k: lineweaver_burk(S, v)[k] for k in ("Vmax", "Km")},
     rand=_rand_mm_data, phrases=["Michaelis-Menten fit", "Michaelis-Menten nonlinear fit"],
     examples=[_ex({"S": [1, 2, 4, 8], "v": [5, 6.6666666667, 8, 8.8888888889]}, {"Vmax": 10.0, "Km": 1.0}),
               _ex({"S": [0.5, 1, 2, 5], "v": [2.2, 3.5, 5.1, 6.6]}, {"Vmax": 8.50005520631, "Km": 1.40001436689})],
     level="advanced")


def _gn_polish(p, model, y, steps=30):
    """Reference helper: Gauss-Newton steps with numpy's least-squares solver from a library fit's optimum (a fit
    that stops on function-value changes is only good to about sqrt(machine epsilon))."""
    import numpy as np
    p = np.array(p, dtype=float)
    last = np.inf
    for _ in range(steps):
        f, J = model(p)
        d = np.linalg.lstsq(J, y - f, rcond=None)[0]
        p = p + d
        size = float(np.max(np.abs(d) / np.maximum(np.abs(p), 1e-300)))
        if size < 1e-16 or size >= last:
            break
        last = size
    return p


def _mm_check(S0, Vmax, Km):
    _pos("S0", S0)
    _pos("Vmax", Vmax)
    _pos("Km", Km)


def mm_progress(S0, Vmax, Km, t, first_order=False):
    _mm_check(S0, Vmax, Km)
    _fin("t", t)
    _need(t >= 0, "t must be zero or positive")
    vt = Vmax * t
    _need(vt / Km <= 600, "the substrate is exhausted beyond representable precision (Vmax t / Km > 600)")
    if first_order:
        S = S0 * math.exp(-vt / Km)
        return {"S": S, "P": S0 - S}
    if vt == 0:
        return {"S": float(S0), "P": 0.0}
    fP = lambda P: -Km * math.log1p(-P / S0) + P - vt          # integrated rate law in the product
    Ph = S0 / 2
    if fP(Ph) >= 0:                                           # less than half converted: solve for P
        P = _bisect(fP, 0.0, Ph)
        return _out(S=S0 - P, P=P)
    fx = lambda x: Km * (math.log(S0) - x) + S0 - math.exp(x) - vt      # more than half: solve for ln S
    x = _bisect(fx, math.log(S0) - vt / Km - 1e-9, math.log(S0 / 2))
    S = math.exp(x)
    return _out(S=S, P=S0 - S)


def _ref_mm_progress(S0, Vmax, Km, t):
    """ODE integration (scipy DOP853) of dS/dt = -Vmax S/(Km + S), with the product carried as its own state."""
    if t == 0:
        return {"S": float(S0), "P": 0.0}
    S, P = _ivp(lambda _t, y: [-Vmax * y[0] / (Km + y[0]), Vmax * y[0] / (Km + y[0])], t, [S0, 0.0])
    return {"S": S, "P": P}


def _rand_mm_progress(rng):
    S0, Km, Vm = round(rng.uniform(1, 100), 2), round(rng.uniform(1, 50), 2), round(rng.uniform(0.1, 10), 3)
    X = rng.uniform(0.05, 0.95)
    t = (Km * -math.log(1 - X) + S0 * X) / Vm
    return {"S0": S0, "Vmax": Vm, "Km": Km, "t": _sig(t)}


_MMP_IN = [_numin("S0", "initial substrate concentration"), _numin("Vmax", "maximum rate (concentration per time)"),
           _numin("Km", "Michaelis constant (concentration)")]

_add("bio_michaelis_menten_progress", name="Substrate remaining in a Michaelis-Menten reaction",
     topic="enzyme kinetics / integrated rate law",
     purpose="The substrate and product concentrations after a time t of an irreversible Michaelis-Menten reaction "
             "(integrated rate law; constant enzyme, no product inhibition).",
     inputs=_MMP_IN + [_numin("t", "elapsed time (time unit of Vmax)")],
     outputs=[_o("S", "substrate concentration at t"), _o("P", "product concentration at t")], primary="S",
     conditions=["positive S0, Vmax, Km", "t >= 0", "consistent units"],
     method="solve Km ln(S0/S) + (S0 - S) = Vmax t by bisection (for P below half conversion, for ln S above)",
     latex=r"K_m\ln\frac{S_0}{S} + (S_0 - S) = V_{\max}t", run=mm_progress, reference=_ref_mm_progress,
     wrong=lambda S0, Vmax, Km, t: mm_progress(S0, Vmax, Km, t, first_order=True), rand=_rand_mm_progress,
     phrases=["Michaelis-Menten progress curve", "Substrate remaining in a Michaelis-Menten reaction"],
     examples=[_ex({"S0": 10, "Vmax": 1, "Km": 2, "t": 5}, {"S": 6.01624392294, "P": 3.98375607706}),
               _ex({"S0": 1, "Vmax": 0.5, "Km": 10, "t": 2}, {"S": 0.912765271609, "P": 0.0872347283914})],
     level="advanced")


def mm_time(S0, Vmax, Km, X, zero_order=False):
    _mm_check(S0, Vmax, Km)
    _fin("X", X)
    _need(0 < X < 1, "the fraction converted must lie strictly between 0 and 1")
    if zero_order:
        return {"t": S0 * X / Vmax}
    return _out(t=(-Km * math.log1p(-X) + S0 * X) / Vmax)


def _ref_mm_time(S0, Vmax, Km, X):
    """Quadrature (scipy quad) of dt = (Km + S)/(Vmax S) dS from S0(1-X) to S0."""
    return {"t": _quad(lambda s: (Km + s) / (Vmax * s), S0 * (1 - X), S0)}


_add("bio_michaelis_menten_time", name="Time to reach a conversion in a Michaelis-Menten reaction",
     topic="enzyme kinetics / integrated rate law",
     purpose="The time an irreversible Michaelis-Menten reaction needs to convert a stated fraction of its substrate.",
     inputs=_MMP_IN + [_numin("X", "fraction of substrate converted (0 to 1, or a percentage)", pct=True)],
     outputs=[_o("t", "time (time unit of Vmax)")], primary="t", conditions=["0 < X < 1", "positive S0, Vmax, Km"],
     method="integrated Michaelis-Menten equation", latex=r"t = \frac{K_m\ln\frac{1}{1-X} + S_0X}{V_{\max}}",
     run=mm_time, reference=_ref_mm_time,
     wrong=lambda S0, Vmax, Km, X: mm_time(S0, Vmax, Km, X, zero_order=True),
     rand=lambda rng: {"S0": round(rng.uniform(1, 100), 2), "Vmax": round(rng.uniform(0.1, 10), 3),
                       "Km": round(rng.uniform(0.5, 50), 2), "X": round(rng.uniform(0.05, 0.95), 3)},
     phrases=["Michaelis-Menten time to conversion", "Time to reach a conversion in a Michaelis-Menten reaction"],
     examples=[_ex({"S0": 10, "Vmax": 1, "Km": 2, "X": 0.5}, {"t": 6.38629436112}),
               _ex({"S0": 100, "Vmax": 5, "Km": 1, "X": 0.9}, {"t": 18.4605170186})],
     extra_q=[{"q": "Michaelis-Menten time to conversion: S0 = 100, Vmax = 5, Km = 1, X = 90%", "outputs": {"t": 18.4605170186}},
              {"q": "Michaelis-Menten time to conversion: S0 = 100, Vmax = 5, Km = 1, X = 90", "outcome": "INVALID_INPUT"}],
     level="advanced")


def hill_plot(S, Y):
    _need(isinstance(S, list) and isinstance(Y, list) and len(S) == len(Y) and len(S) >= 2,
          "S and Y must be lists of equal length (at least 2 points)")
    for s in S:
        _pos("S", s)
    for y in Y:
        _frac("Y", y, True, True)
    b, a, r2 = _ols([math.log10(s) for s in S], [math.log10(y / (1 - y)) for y in Y])
    _need(b > 0, "the Hill plot must rise (positive Hill coefficient)")
    return _out(n=b, K=10 ** (-a / b), r2=r2)


def _ref_hill_plot(S, Y):
    import numpy as np
    b, a, r2 = _linregress(np.log10(S), np.log10(np.array(Y) / (1 - np.array(Y))))
    return {"n": b, "K": 10 ** (-a / b), "r2": r2}


def _rand_hill_plot(rng):
    K, n = rng.uniform(1, 50), rng.uniform(0.8, 3.5)
    S = sorted({_sig(K * m) for m in rng.sample([0.2, 0.3, 0.5, 0.7, 1, 1.4, 2, 3, 5], rng.randint(4, 7))})
    Y = [min(max(_sig(s ** n / (K ** n + s ** n) * (1 + rng.uniform(-0.03, 0.03))), 0.001), 0.999) for s in S]
    return {"S": S, "Y": Y}


_add("bio_hill_plot", name="Hill coefficient from a Hill plot", topic="enzyme kinetics / cooperativity",
     purpose="The Hill coefficient n and half-saturation constant K from a least-squares line through the Hill plot "
             "of log10(Y/(1-Y)) against log10 S (Y = fractional saturation).",
     inputs=[_vecin("S", "ligand or substrate concentrations", lo=2), _vecin("Y", "fractional saturations (0 < Y < 1)", lo=2)],
     outputs=[_o("n", "Hill coefficient"), _o("K", "half-saturation concentration (unit of S)"),
              _o("r2", "coefficient of determination of the line")],
     primary="n", conditions=["0 < Y < 1", "positive S", "rising Hill plot"],
     method="ordinary least squares of log10(Y/(1-Y)) on log10 S: slope n, intercept -n log10 K",
     latex=r"\log\frac{Y}{1-Y} = n\log S - n\log K", run=hill_plot, reference=_ref_hill_plot,
     wrong=lambda S, Y: (lambda r: dict(r, K=r["K"] ** r["n"]))(hill_plot(S, Y)), rand=_rand_hill_plot,
     phrases=["Hill plot", "Hill coefficient from a Hill plot"],
     examples=[_ex({"S": [1, 2, 4, 8], "Y": [0.0588235294118, 0.2, 0.5, 0.8]}, {"n": 2.0, "K": 4.0}),
               _ex({"S": [10, 20, 40], "Y": [0.1, 0.4, 0.8]}, {"n": 2.58496250072, "K": 23.3965269138})],
     level="intermediate")


def _solve3(A, b):
    M = [list(A[i]) + [b[i]] for i in range(3)]
    for c in range(3):
        p = max(range(c, 3), key=lambda r: abs(M[r][c]))
        if abs(M[p][c]) < 1e-300:
            raise ValueError("singular normal equations")
        M[c], M[p] = M[p], M[c]
        for r in range(c + 1, 3):
            f = M[r][c] / M[c][c]
            for k in range(c, 4):
                M[r][k] -= f * M[c][k]
    x = [0.0] * 3
    for r in (2, 1, 0):
        x[r] = (M[r][3] - sum(M[r][k] * x[k] for k in range(r + 1, 3))) / M[r][r]
    return x


def _lm(model, p, ys, ok):
    """Levenberg-Marquardt for 3 parameters: model(p) -> (predictions, jacobian rows)."""
    def cost(pp):
        f, _ = model(pp)
        return math.fsum((y - a) ** 2 for y, a in zip(ys, f))
    lam, c = 1e-3, cost(p)
    for _ in range(2000):
        f, J = model(p)
        r = [y - a for y, a in zip(ys, f)]
        A = [[math.fsum(row[i] * row[j] for row in J) for j in range(3)] for i in range(3)]
        g = [math.fsum(row[i] * ri for row, ri in zip(J, r)) for i in range(3)]
        while True:
            Ad = [[A[i][j] + (lam * A[i][i] if i == j else 0.0) for j in range(3)] for i in range(3)]
            try:
                d = _solve3(Ad, g)
            except ValueError:
                d = None
            q = [p[i] + d[i] for i in range(3)] if d else None
            if q and ok(q):
                cq = cost(q)
                if cq < c:
                    small = max(abs(d[i]) / max(abs(p[i]), 1e-300) for i in range(3)) < 1e-15
                    p, c = q, cq
                    lam = max(lam / 10, 1e-15)
                    if small:
                        return _gn(model, p, ys, ok)
                    break
            lam *= 10
            if lam > 1e20:
                return _gn(model, p, ys, ok)
    return _gn(model, p, ys, ok)


def _gn(model, p, ys, ok):
    """Gauss-Newton polish of a least-squares optimum: steps from the normal equations while they keep shrinking (the
    gradient, not the cost, decides; cost changes vanish below rounding near the optimum)."""
    last = math.inf
    for _ in range(60):
        f, J = model(p)
        r = [y - a for y, a in zip(ys, f)]
        A = [[math.fsum(row[i] * row[j] for row in J) for j in range(3)] for i in range(3)]
        g = [math.fsum(row[i] * ri for row, ri in zip(J, r)) for i in range(3)]
        try:
            d = _solve3(A, g)
        except ValueError:
            return p
        q = [p[i] + d[i] for i in range(3)]
        if not ok(q):
            return p
        size = max(abs(d[i]) / max(abs(q[i]), 1e-300) for i in range(3))
        p = q
        if size < 1e-16 or size >= last:
            return p
        last = size
    return p


def _hill_model(S):
    def m(p):
        Vm, K, n = p
        f, J = [], []
        for s in S:
            u = (s / K) ** n
            x = u / (1 + u)
            f.append(Vm * x)
            dx = x * (1 - x)
            J.append([x, -Vm * dx * n / K, Vm * dx * math.log(s / K)])
        return f, J
    return m


def hill_fit(S, v, fixed_n=False):
    S, v = _kin_data(S, v, lo=4)
    _need(len(set(S)) >= 4, "at least 4 different concentrations are needed")
    if fixed_n:
        r = mm_fit(S, v)
        return {"Vmax": r["Vmax"], "K": r["Km"], "n": 1.0}
    vmax0 = max(v) * 1.1
    Y = [min(max(y / vmax0, 1e-3), 1 - 1e-3) for y in v]
    b, a, _ = _ols([math.log(s) for s in S], [math.log(y / (1 - y)) for y in Y])
    p0 = [vmax0, math.exp(-a / b) if b > 0 else sorted(S)[len(S) // 2], max(b, 0.5)]
    p = _lm(_hill_model(S), p0, v, lambda q: q[0] > 0 and q[1] > 0 and 0.05 < q[2] < 20)
    f, J = _hill_model(S)(p)
    r = [y - a for y, a in zip(v, f)]
    g = [abs(math.fsum(row[i] * ri for row, ri in zip(J, r))) for i in range(3)]
    scale = math.fsum(y * y for y in v)
    _need(max(gi * abs(pi) for gi, pi in zip(g, p)) <= 1e-9 * scale, "the Hill fit did not converge")
    return _out(Vmax=p[0], K=p[1], n=p[2])


def _ref_hill_fit(S, v):
    import numpy as np
    from scipy.optimize import curve_fit
    S, v = np.array(S, float), np.array(v, float)
    p, _ = curve_fit(lambda s, Vm, K, n: Vm * (s / K) ** n / (1 + (s / K) ** n), S, v,
                     p0=[1.1 * v.max(), float(np.median(S)), 1.0], method="lm", xtol=1e-15, ftol=1e-15, gtol=1e-15,
                     maxfev=200000)

    def model(q):
        u = (S / q[1]) ** q[2]
        x = u / (1 + u)
        dx = x * (1 - x)
        return q[0] * x, np.column_stack([x, -q[0] * dx * q[2] / q[1], q[0] * dx * np.log(S / q[1])])
    p = _gn_polish(p, model, v)
    return {"Vmax": float(p[0]), "K": float(p[1]), "n": float(p[2])}


def _rand_hill(rng):
    Vm, K, n = rng.uniform(5, 100), rng.uniform(1, 50), rng.uniform(0.7, 3.5)
    mult = sorted(rng.sample([0.1, 0.2, 0.35, 0.5, 0.75, 1, 1.5, 2, 3, 5, 8, 12], rng.randint(7, 10)))
    S = [_sig(K * m) for m in mult]
    v = [_sig(Vm * (s / K) ** n / (1 + (s / K) ** n) * (1 + rng.uniform(-0.02, 0.02))) for s in S]
    return {"S": S, "v": v}


_add("bio_hill_fit", name="Hill equation nonlinear least-squares fit", topic="enzyme kinetics / cooperativity",
     purpose="Vmax, the half-saturation constant K and the Hill coefficient n fitted to rate or binding data by "
             "nonlinear least squares (unweighted).",
     inputs=[_vecin("S", "concentrations (any one unit)", lo=4), _vecin("v", "rates or responses (any one unit)", lo=4)],
     outputs=[_o("Vmax", "maximum rate"), _o("K", "half-saturation concentration (unit of S)"), _o("n", "Hill coefficient")],
     primary="n", conditions=["at least 4 distinct positive concentrations", "positive responses", "convergence"],
     method="Levenberg-Marquardt on v = Vmax S^n / (K^n + S^n), started from the Hill-plot line",
     latex=r"v = \frac{V_{\max}S^n}{K^n + S^n}", run=hill_fit, reference=_ref_hill_fit,
     wrong=lambda S, v: hill_fit(S, v, fixed_n=True), rand=_rand_hill,
     phrases=["Hill fit", "Hill equation fit"],
     examples=[_ex({"S": [1, 2, 4, 8, 16], "v": [0.588235294118, 2, 5, 8, 9.41176470588]},
                   {"Vmax": 10.0, "K": 4.0, "n": 2.0}),
               _ex({"S": [0.5, 1, 2, 4, 8, 16], "v": [0.3, 1.1, 3.2, 6.1, 8.3, 9.2]}, {"Vmax": 9.58247185219, "K": 2.93749695471, "n": 1.86334907827})],
     level="expert")


def _ci_model(S, I):
    def m(p):
        Vm, Km, Ki = p
        f, J = [], []
        for s, i in zip(S, I):
            den = Km * (1 + i / Ki) + s
            f.append(Vm * s / den)
            J.append([s / den, -Vm * s * (1 + i / Ki) / den ** 2, Vm * s * Km * i / (Ki * Ki) / den ** 2])
        return f, J
    return m


def _nc_model(S, I):
    def m(p):
        Vm, Km, Ki = p
        f, J = [], []
        for s, i in zip(S, I):
            a = 1 + i / Ki
            f.append(Vm * s / ((Km + s) * a))
            J.append([s / ((Km + s) * a), -Vm * s / ((Km + s) ** 2 * a), Vm * s * i / (Ki * Ki) / ((Km + s) * a * a)])
        return f, J
    return m


def competitive_fit(S, I, v, model="competitive"):
    _need(isinstance(S, list) and isinstance(I, list) and isinstance(v, list) and len(S) == len(I) == len(v) >= 4,
          "S, I and v must be lists of equal length (at least 4 points)")
    for x in S + v:
        _pos("every substrate concentration and rate", x)
    for x in I:
        _fin("I", x)
        _need(x >= 0, "inhibitor concentrations must be zero or positive")
    _need(any(x == 0 for x in I) and any(x > 0 for x in I), "include rates without inhibitor and with inhibitor")
    S, I, v = [float(x) for x in S], [float(x) for x in I], [float(x) for x in v]
    base = [(s, y) for s, i, y in zip(S, I, v) if i == 0]
    _need(len({s for s, _ in base}) >= 3, "at least 3 different S without inhibitor are needed")
    r0 = mm_fit([s for s, _ in base], [y for _, y in base])
    mdl = _ci_model(S, I) if model == "competitive" else _nc_model(S, I)
    p = _lm(mdl, [r0["Vmax"], r0["Km"], max(I)], v, lambda q: q[0] > 0 and q[1] > 0 and q[2] > 0)
    f, J = mdl(p)
    r = [y - a for y, a in zip(v, f)]
    g = [abs(math.fsum(row[i] * ri for row, ri in zip(J, r))) for i in range(3)]
    _need(max(gi * abs(pi) for gi, pi in zip(g, p)) <= 1e-9 * math.fsum(y * y for y in v),
          "the inhibition fit did not converge")
    return _out(Vmax=p[0], Km=p[1], Ki=p[2])


def _ref_competitive(S, I, v):
    import numpy as np
    from scipy.optimize import curve_fit
    X = np.vstack([np.array(S, float), np.array(I, float)])
    p, _ = curve_fit(lambda x, Vm, Km, Ki: Vm * x[0] / (Km * (1 + x[1] / Ki) + x[0]), X, np.array(v, float),
                     p0=[max(v), float(np.median(S)), float(max(I))], method="lm", xtol=1e-15, ftol=1e-15, gtol=1e-15,
                     maxfev=200000)
    s, i = X

    def model(q):
        den = q[1] * (1 + i / q[2]) + s
        return q[0] * s / den, np.column_stack([s / den, -q[0] * s * (1 + i / q[2]) / den ** 2,
                                                 q[0] * s * q[1] * i / q[2] ** 2 / den ** 2])
    p = _gn_polish(p, model, np.array(v, float))
    return {"Vmax": float(p[0]), "Km": float(p[1]), "Ki": float(p[2])}


def _rand_competitive(rng):
    Vm, Km, Ki = rng.uniform(5, 100), rng.uniform(1, 30), rng.uniform(0.5, 20)
    S, I, v = [], [], []
    mults = sorted(rng.sample([0.25, 0.5, 1, 2, 4, 8], 5))
    for ii in (0.0, _sig(Ki), _sig(3 * Ki)):
        for m in mults:
            s = _sig(Km * m)
            S.append(s)
            I.append(ii)
            v.append(_sig(Vm * s / (Km * (1 + ii / Ki) + s) * (1 + rng.uniform(-0.02, 0.02))))
    return {"S": S, "I": I, "v": v}


_add("bio_competitive_inhibition_fit", name="Competitive inhibition global fit", topic="enzyme kinetics / inhibition",
     purpose="Vmax, Km and the inhibition constant Ki of a competitive inhibitor fitted globally to rates measured at "
             "several substrate and inhibitor concentrations (nonlinear least squares, unweighted).",
     inputs=[_vecin("S", "substrate concentrations (one per measurement)", lo=4),
             _vecin("I", "inhibitor concentrations (one per measurement, 0 for none)", lo=4),
             _vecin("v", "measured rates", lo=4)],
     outputs=[_o("Vmax", "maximum rate"), _o("Km", "Michaelis constant (unit of S)"),
              _o("Ki", "inhibition constant (unit of I)")],
     primary="Ki", conditions=["rates without inhibitor at 3 or more S", "rates with inhibitor", "convergence"],
     method="Levenberg-Marquardt on v = Vmax S / (Km (1 + I/Ki) + S), started from the uninhibited fit",
     latex=r"v = \frac{V_{\max}S}{K_m\left(1 + \frac{[I]}{K_i}\right) + S}", run=competitive_fit,
     reference=_ref_competitive, wrong=lambda S, I, v: competitive_fit(S, I, v, model="noncompetitive"),
     rand=_rand_competitive, phrases=["Competitive inhibition fit", "Competitive inhibition global fit"],
     examples=[_ex({"S": [1, 2, 4, 8, 1, 2, 4, 8], "I": [0, 0, 0, 0, 2, 2, 2, 2],
                    "v": [5, 6.66666666667, 8, 8.88888888889, 2.5, 4, 5.71428571429, 7.27272727273]},
                   {"Vmax": 10.0, "Km": 1.0, "Ki": 1.0})],
     level="expert")


# ================================================================================================ PROTEIN COMPOSITION
def aromaticity(seq):
    s = _seq(seq, "seq", "PROTEIN")
    n = s.count("F") + s.count("W") + s.count("Y")
    return {"aromaticity": n / len(s), "aromatic_count": n}


def _ref_aromaticity(seq):
    from Bio.SeqUtils.ProtParam import ProteinAnalysis
    pa = ProteinAnalysis(seq.upper())
    c = pa.count_amino_acids()
    return {"aromaticity": float(pa.aromaticity()), "aromatic_count": int(c["F"] + c["W"] + c["Y"])}


_add("bio_aromaticity", name="Aromaticity of a protein", topic="sequences / protein properties",
     purpose="The aromaticity (relative frequency of Phe + Trp + Tyr, Lobry 1994) of a protein and its number of "
             "aromatic residues.",
     inputs=[_textin("seq", "protein sequence")],
     outputs=[_o("aromaticity", "fraction of F, W and Y residues"), _o("aromatic_count", "number of F, W and Y residues")],
     primary="aromaticity", conditions=["standard amino acids only"], method="count F, W and Y and divide by the length",
     latex=r"\mathrm{Ar} = \frac{n_F + n_W + n_Y}{n}", run=aromaticity, reference=_ref_aromaticity,
     wrong=lambda seq: (lambda s: {"aromaticity": (s.count("F") + s.count("W") + s.count("Y") + s.count("H")) / len(s),
                                   "aromatic_count": s.count("F") + s.count("W") + s.count("Y") + s.count("H")})(seq.upper()),
     rand=lambda rng: {"seq": _rand_protein(rng, 20, 120)},
     phrases=["Aromaticity", "Protein aromaticity"],
     examples=[_ex({"seq": "MFWYAK"}, {"aromaticity": 0.5, "aromatic_count": 3}),
               _ex({"seq": "GGGG"}, {"aromaticity": 0.0, "aromatic_count": 0})],
     level="beginner")


def aa_composition(seq):
    s = _seq(seq, "seq", "PROTEIN")
    pct = {a: 100.0 * s.count(a) / len(s) for a in AMINO}
    top = max(pct.values())
    return _out(percent=pct, most_common=min(a for a in AMINO if pct[a] == top), length=len(s))


def _ref_aa_composition(seq):
    from Bio.SeqUtils.ProtParam import ProteinAnalysis
    pa = ProteinAnalysis(seq.upper())
    pct = {a: float(pa.amino_acids_percent[a]) for a in AMINO}
    top = max(pct.values())
    return {"percent": pct, "most_common": min(a for a in AMINO if pct[a] == top), "length": int(pa.length)}


def _wrong_aa_comp(seq):
    s = seq.upper()
    pct = {a: 100.0 * s.count(a) / (len(s) + 1) for a in AMINO}            # length off by one
    top = max(pct.values())
    return {"percent": pct, "most_common": min(a for a in AMINO if pct[a] == top), "length": len(s) + 1}


_add("bio_amino_acid_composition", name="Amino-acid composition of a protein", topic="sequences / protein properties",
     purpose="The percentage of each of the 20 standard amino acids in a protein, its most common residue "
             "(alphabetically first among ties) and its length.",
     inputs=[_textin("seq", "protein sequence")],
     outputs=[_o("percent", "percentage of each amino acid (one-letter code -> %)"),
              _o("most_common", "most common amino acid (alphabetically first among ties)", "text"),
              _o("length", "number of residues")],
     primary="length", conditions=["standard amino acids only"], method="count each residue and divide by the length",
     latex=r"\%a = 100\,\frac{n_a}{n}", run=aa_composition, reference=_ref_aa_composition, wrong=_wrong_aa_comp,
     rand=lambda rng: {"seq": _rand_protein(rng, 5, 120)},
     phrases=["Amino acid composition", "Amino-acid composition"],
     examples=[_ex({"seq": "MAAAK"}, {"most_common": "A", "length": 5}),
               _ex({"seq": _PROT_EX}, {"length": 152, "most_common": "L"})],
     level="beginner")


_HELIX, _TURN, _SHEET = "EMALK", "NPGSD", "VIYFWLT"


def ss_fraction(seq, order=("helix", "turn", "sheet")):
    s = _seq(seq, "seq", "PROTEIN")
    f = {k: sum(s.count(a) for a in g) / len(s) for k, g in zip(order, (_HELIX, _TURN, _SHEET))}
    return f


def _ref_ss(seq):
    from Bio.SeqUtils.ProtParam import ProteinAnalysis
    h, t, sh = ProteinAnalysis(seq.upper()).secondary_structure_fraction()
    return {"helix": float(h), "turn": float(t), "sheet": float(sh)}


_add("bio_secondary_structure_fraction", name="Secondary-structure propensity fractions",
     topic="sequences / protein properties",
     purpose="The fractions of residues that tend to be in helices (E M A L K), turns (N P G S D) and sheets "
             "(V I Y F W L T), as ProtParam defines them (L counts for helix and sheet).",
     inputs=[_textin("seq", "protein sequence")],
     outputs=[_o("helix", "fraction of helix-favouring residues"), _o("turn", "fraction of turn-favouring residues"),
              _o("sheet", "fraction of sheet-favouring residues")],
     primary="helix", conditions=["standard amino acids only"], method="count the residues of each propensity group",
     latex=r"f_{\text{helix}} = \frac{n_E+n_M+n_A+n_L+n_K}{n}", run=ss_fraction, reference=_ref_ss,
     wrong=lambda seq: ss_fraction(seq, order=("sheet", "turn", "helix")),
     rand=lambda rng: {"seq": _rand_protein(rng, 10, 120)},
     phrases=["Secondary structure fraction", "Secondary-structure fractions"],
     examples=[_ex({"seq": "EMALKNPGSD"}, {"helix": 0.5, "turn": 0.5, "sheet": 0.1}),
               _ex({"seq": "VVVV"}, {"helix": 0.0, "turn": 0.0, "sheet": 1.0})],
     level="intermediate")


# ================================================================================================ POPULATION MODELS
def _leslie_in(F, S):
    _need(isinstance(F, list) and isinstance(S, list) and len(F) >= 2 and len(S) == len(F) - 1,
          "F needs one fecundity per age class (at least 2) and S one survival fewer")
    for f in F:
        _fin("F", f)
        _need(f >= 0, "fecundities must be zero or positive")
    for s in S:
        _fin("S", s)
        _need(0 < s <= 1, "survival probabilities must lie in (0, 1]")
    _need(any(f > 0 for f in F), "at least one age class must reproduce")


def leslie_projection(F, S, n0, t, transpose=False):
    _leslie_in(F, S)
    _need(isinstance(n0, list) and len(n0) == len(F), "n0 needs one number per age class")
    for x in n0:
        _fin("n0", x)
        _need(x >= 0, "population numbers must be zero or positive")
    _need(isinstance(t, int) and 0 <= t <= 500, "t must be a whole number of time steps (0 to 500)")
    k = len(F)
    L = [[0.0] * k for _ in range(k)]
    L[0] = [float(f) for f in F]
    for i, s in enumerate(S):
        L[i + 1][i] = float(s)
    if transpose:
        L = [list(r) for r in zip(*L)]
    n = [float(x) for x in n0]
    for _ in range(t):
        n = [math.fsum(L[i][j] * n[j] for j in range(k)) for i in range(k)]
    return _out(n_t=n, total=math.fsum(n))


def _ref_leslie_projection(F, S, n0, t):
    import numpy as np
    k = len(F)
    L = np.zeros((k, k))
    L[0, :] = F
    for i, s in enumerate(S):
        L[i + 1, i] = s
    n = np.linalg.matrix_power(L, t) @ np.array(n0, dtype=float)
    return {"n_t": [float(x) for x in n], "total": float(n.sum())}


def _rand_leslie(rng, with_n0=True):
    k = rng.randint(2, 5)
    F = [0.0] + [round(rng.uniform(0.2, 3), 2) for _ in range(k - 1)]
    if rng.random() < 0.3:
        F[0] = round(rng.uniform(0, 1), 2)
    S = [round(rng.uniform(0.2, 0.95), 2) for _ in range(k - 1)]
    d = {"F": F, "S": S}
    if with_n0:
        d.update(n0=[float(rng.randint(0, 100)) for _ in range(k)], t=rng.randint(1, 25))
        if not any(d["n0"]):
            d["n0"][0] = 10.0
    return d


_add("bio_leslie_projection", name="Leslie matrix population projection", topic="population biology / age structure",
     purpose="The age-class numbers and total population after t time steps of an age-structured (Leslie) model.",
     inputs=[_vecin("F", "fecundity of each age class", lo=2, hi=20),
             _vecin("S", "survival from each age class to the next (one fewer than F)", lo=1, hi=19),
             _vecin("n0", "initial number in each age class", lo=2, hi=20), _intin("t", "number of time steps", lo=0)],
     outputs=[_o("n_t", "number in each age class after t steps"), _o("total", "total population after t steps")],
     primary="total", conditions=["0 < S_i <= 1", "F_i >= 0, at least one positive", "t <= 500"],
     method="multiply the age vector by the Leslie matrix t times",
     latex=r"\mathbf n_{t+1} = L\,\mathbf n_t,\quad L = \begin{pmatrix}F_1 & F_2 & \cdots\\ S_1 & 0 & \\ & \ddots & \end{pmatrix}",
     run=leslie_projection, reference=_ref_leslie_projection,
     wrong=lambda F, S, n0, t: leslie_projection(F, S, n0, t, transpose=True), rand=_rand_leslie,
     phrases=["Leslie projection", "Leslie matrix projection"],
     examples=[_ex({"F": [0, 1.5, 1], "S": [0.5, 0.8], "n0": [100, 50, 20], "t": 1},
                   {"n_t": [95.0, 50.0, 40.0], "total": 185.0}),
               _ex({"F": [0, 2], "S": [0.5], "n0": [10, 0], "t": 2}, {"n_t": [10.0, 0.0], "total": 10.0})],
     level="intermediate")


def _euler_lotka(F, S):
    l = [1.0]
    for s in S:
        l.append(l[-1] * s)
    return l


def leslie_growth(F, S, use_r0=False):
    _leslie_in(F, S)
    l = _euler_lotka(F, S)
    R0 = math.fsum(f * li for f, li in zip(F, l))
    if use_r0:
        lam = R0
    else:
        g = lambda x: math.fsum(f * li * math.exp(-(i + 1) * x) for i, (f, li) in enumerate(zip(F, l))) - 1.0
        hi = math.log(max(R0, 1.0)) + 1.0
        lo = -1.0
        while g(lo) <= 0:
            lo -= 1.0
        while g(hi) >= 0:
            hi += 1.0
        lam = math.exp(_bisect(g, lo, hi))
    w = [li * lam ** (-i) for i, li in enumerate(l)]
    tot = math.fsum(w)
    return _out(lambda1=lam, r=math.log(lam), R0=R0, stable=[x / tot for x in w])


def _ref_leslie_growth(F, S):
    """numpy eigen-decomposition of the Leslie matrix; R0 as the spectral radius of the next-generation matrix
    F (I - T)^-1."""
    import numpy as np
    k = len(F)
    L = np.zeros((k, k))
    L[0, :] = F
    for i, s in enumerate(S):
        L[i + 1, i] = s
    vals, vecs = np.linalg.eig(L)
    i = int(np.argmax(vals.real))
    lam = float(vals[i].real)
    w = np.abs(vecs[:, i].real)
    Fm = np.zeros((k, k))
    Fm[0, :] = F
    T = L - Fm
    R0 = float(max(abs(np.linalg.eigvals(Fm @ np.linalg.inv(np.eye(k) - T)))))
    return {"lambda1": lam, "r": float(np.log(lam)), "R0": R0, "stable": [float(x) for x in w / w.sum()]}


_add("bio_leslie_growth_rate", name="Leslie matrix growth rate and stable age distribution",
     topic="population biology / age structure",
     purpose="The asymptotic growth rate (dominant eigenvalue lambda and r = ln lambda), the net reproductive rate R0 "
             "and the stable age distribution of a Leslie model.",
     inputs=[_vecin("F", "fecundity of each age class", lo=2, hi=20),
             _vecin("S", "survival from each age class to the next (one fewer than F)", lo=1, hi=19)],
     outputs=[_o("lambda1", "dominant eigenvalue (growth factor per step)"), _o("r", "intrinsic rate ln(lambda)"),
              _o("R0", "net reproductive rate (expected offspring per newborn)"),
              _o("stable", "stable age distribution (fractions summing to 1)")],
     primary="lambda1", conditions=["0 < S_i <= 1", "F_i >= 0, at least one positive"],
     method="the positive root of the Euler-Lotka equation sum F_i l_i lambda^-i = 1; stable distribution "
            "proportional to l_i lambda^-(i-1)",
     latex=r"\sum_{i} F_i\,l_i\,\lambda^{-i} = 1,\quad R_0 = \sum_i F_i l_i", run=leslie_growth,
     reference=_ref_leslie_growth, wrong=lambda F, S: leslie_growth(F, S, use_r0=True),
     rand=lambda rng: _rand_leslie(rng, with_n0=False),
     phrases=["Leslie growth rate", "Leslie matrix growth rate"],
     examples=[_ex({"F": [0, 2], "S": [0.5]}, {"lambda1": 1.0, "R0": 1.0, "stable": [2 / 3, 1 / 3]}),
               _ex({"F": [0, 1.5, 1], "S": [0.5, 0.8]}, {"lambda1": 1.06151735158, "R0": 1.15, "stable": [0.547643499894, 0.257953155019, 0.194403345087]})],
     level="advanced")


# ================================================================================================ DIVERSITY
def _abund(counts, name="counts"):
    c = _counts(name, counts, 0)
    pos = [x for x in c if x > 0]
    _need(pos, f"{name}: at least one species must have a positive count")
    return pos


def _rand_counts(rng, lo=3, hi=25, single_rich=True):
    k = rng.randint(lo, hi)
    c = [int(_lu(rng, 1, 200)) for _ in range(k)]
    if single_rich:
        for _ in range(rng.randint(1, 4)):
            c[rng.randrange(k)] = rng.choice([1, 1, 2])
    return c


_LOGS = {"e": math.e, "2": 2.0, "10": 10.0}


def shannon(counts, base, wrong=False):
    c = _abund(counts)
    _need(base in _LOGS, "base must be e, 2 or 10")
    N = sum(c)
    H = -math.fsum((x / N) * math.log(x / N) for x in c)
    S = len(c)
    Hb = H / math.log(_LOGS[base])
    J = (1.0 if S == 1 else (H / S if wrong else H / math.log(S)))
    return _out(H=Hb, evenness=J, effective_number=math.exp(H))


def _ref_shannon(counts, base):
    from skbio.diversity.alpha import pielou_e, shannon as sh
    c = [int(x) for x in counts if x > 0]
    b = None if base == "e" else float(base)
    return {"H": float(sh(c, base=b)), "evenness": float(pielou_e(c, base=b)), "effective_number": float(sh(c, exp=True))}


_COUNTS_IN = _vecin("counts", "number of individuals of each species", lo=1, hi=500)

_add("bio_shannon_diversity", name="Shannon diversity index", topic="ecology / diversity indices",
     purpose="The Shannon diversity index H' of a community in the stated logarithm base, Pielou's evenness J' and "
             "the effective number of species exp(H').",
     inputs=[_COUNTS_IN, _choice("base", ["e", "2", "10"], "logarithm base: e (nats), 2 (bits) or 10")],
     outputs=[_o("H", "Shannon index H' in the stated base"), _o("evenness", "Pielou's evenness J' = H'/ln S"),
              _o("effective_number", "effective number of species exp(H') (Hill number of order 1)")],
     primary="H", conditions=["whole non-negative counts", "zero counts are not species of the sample"],
     method="H' = -sum p_i log p_i with p_i = n_i / N", latex=r"H' = -\sum_i p_i \log_b p_i,\quad J' = \frac{H'}{\log_b S}",
     run=shannon, reference=_ref_shannon, wrong=lambda counts, base: shannon(counts, base, wrong=True),
     rand=lambda rng: {"counts": _rand_counts(rng, 2, 20), "base": rng.choice(["e", "2", "10"])},
     phrases=["Shannon diversity", "Shannon index"],
     examples=[_ex({"counts": [10, 10, 10, 10], "base": "e"}, {"H": 1.38629436112, "evenness": 1.0, "effective_number": 4.0}),
               _ex({"counts": [50, 30, 20], "base": "2"}, {"H": 1.48547529723})],
     extra_q=[{"q": "Shannon diversity: counts = [10, 10, 10, 10]", "outcome": "NEEDS_INPUT"}],
     level="beginner")


def simpson(counts, estimator, swap=False):
    c = _abund(counts)
    _need(estimator in ("infinite", "finite"), "estimator must be infinite or finite")
    est = {"infinite": "finite", "finite": "infinite"}[estimator] if swap else estimator
    N = sum(c)
    if est == "finite":
        _need(N >= 2, "the finite-sample estimator needs at least 2 individuals")
        D = Fraction(sum(x * (x - 1) for x in c), N * (N - 1))
        _need(D > 0, "every species is a singleton: the finite-sample dominance is 0 and its inverse undefined")
    else:
        D = Fraction(sum(x * x for x in c), N * N)
    return {"D": float(D), "simpson": float(1 - D), "inverse": float(1 / D)}


def _ref_simpson(counts, estimator):
    from skbio.diversity.alpha import dominance, enspie, simpson as si
    c = [int(x) for x in counts if x > 0]
    f = estimator == "finite"
    return {"D": float(dominance(c, finite=f)), "simpson": float(si(c, finite=f)), "inverse": float(enspie(c, finite=f))}


_add("bio_simpson_diversity", name="Simpson diversity indices", topic="ecology / diversity indices",
     purpose="Simpson's dominance D, Simpson's diversity 1 - D and the inverse Simpson index 1/D, with the stated "
             "estimator (infinite population: sum p^2; finite sample: sum n(n-1)/(N(N-1))).",
     inputs=[_COUNTS_IN, _choice("estimator", ["infinite", "finite"], "infinite (sum p^2) or finite-sample correction")],
     outputs=[_o("D", "Simpson's dominance D"), _o("simpson", "Simpson's diversity 1 - D"), _o("inverse", "inverse Simpson 1/D")],
     primary="simpson", conditions=["whole non-negative counts", "estimator stated"],
     method="D = sum p_i^2 or sum n_i(n_i - 1)/(N(N - 1))",
     latex=r"D = \sum_i p_i^2 \quad\text{or}\quad D = \frac{\sum n_i(n_i-1)}{N(N-1)}", run=simpson,
     reference=_ref_simpson, wrong=lambda counts, estimator: simpson(counts, estimator, swap=True),
     rand=lambda rng: {"counts": _rand_counts(rng, 2, 20), "estimator": rng.choice(["infinite", "finite"])},
     phrases=["Simpson diversity", "Simpson index"],
     examples=[_ex({"counts": [10, 10, 10, 10], "estimator": "infinite"}, {"D": 0.25, "simpson": 0.75, "inverse": 4.0}),
               _ex({"counts": [3, 2, 1], "estimator": "finite"}, {"D": 4 / 15, "simpson": 11 / 15, "inverse": 3.75})],
     level="beginner")


def _osd(c):
    return len(c), sum(1 for x in c if x == 1), sum(1 for x in c if x == 2)


def chao1(counts, form, swap=False):
    c = _abund(counts)
    _need(form in ("bias_corrected", "classic"), "form must be bias_corrected or classic")
    if swap:
        form = "classic" if form == "bias_corrected" else "bias_corrected"
    o, s, d = _osd(c)
    if form == "classic" and s and d:
        est = o + Fraction(s * s, 2 * d)
    else:
        est = o + Fraction(s * (s - 1), 2 * (d + 1))
    return {"chao1": float(est), "s_obs": o, "singletons": s, "doubletons": d}


def _ref_chao1(counts, form):
    from skbio.diversity.alpha import chao1 as ch, osd
    c = [int(x) for x in counts if x > 0]
    o, s, d = osd(c)
    return {"chao1": float(ch(c, bias_corrected=(form == "bias_corrected"))), "s_obs": int(o), "singletons": int(s),
            "doubletons": int(d)}


def _rand_chao(rng):
    c = _rand_counts(rng, 5, 25)
    c += [1, 1, 2]
    rng.shuffle(c)
    return {"counts": c, "form": rng.choice(["bias_corrected", "classic"])}


_add("bio_chao1", name="Chao1 species richness estimator", topic="ecology / richness estimation",
     purpose="The Chao1 estimate of total species richness from abundance counts (bias-corrected or classic form), "
             "with the observed richness and the numbers of singletons and doubletons.",
     inputs=[_COUNTS_IN, _choice("form", ["bias_corrected", "classic"], "bias_corrected or classic Chao1")],
     outputs=[_o("chao1", "estimated species richness"), _o("s_obs", "observed number of species"),
              _o("singletons", "species seen once (F1)"), _o("doubletons", "species seen twice (F2)")],
     primary="chao1", conditions=["whole counts", "the classic form needs F1 > 0 and F2 > 0, otherwise the "
                                                   "bias-corrected form is used (as scikit-bio does)"],
     method="S_obs + F1^2/(2 F2) (classic) or S_obs + F1(F1 - 1)/(2(F2 + 1)) (bias-corrected)",
     latex=r"\hat S = S_{obs} + \frac{F_1(F_1-1)}{2(F_2+1)}", run=chao1, reference=_ref_chao1,
     wrong=lambda counts, form: chao1(counts, form, swap=True), rand=_rand_chao,
     phrases=["Chao1", "Chao1 richness estimator"],
     examples=[_ex({"counts": [1, 1, 1, 2, 2, 5, 10], "form": "classic"},
                   {"chao1": 9.25, "s_obs": 7, "singletons": 3, "doubletons": 2}),
               _ex({"counts": [1, 1, 1, 2, 2, 5, 10], "form": "bias_corrected"}, {"chao1": 8.0})],
     level="intermediate")


def ace(counts, rare_threshold, no_gamma=False):
    c = _abund(counts)
    _need(isinstance(rare_threshold, int) and rare_threshold >= 1, "rare_threshold must be a whole number >= 1")
    rare = [x for x in c if x <= rare_threshold]
    s_rare = len(rare)
    singles = sum(1 for x in c if x == 1)
    _need(not (singles > 0 and singles == s_rare), "the only rare species are singletons, so ACE is undefined "
                                                     "(use the bias-corrected Chao1)")
    s_abun = len(c) - s_rare
    if s_rare == 0:
        return {"ace": float(s_abun)}
    n_rare = sum(rare)
    c_ace = 1 - Fraction(singles, n_rare)
    top = s_rare * sum(x * (x - 1) for x in rare)
    gamma = Fraction(top) / (c_ace * n_rare * (n_rare - 1)) - 1
    gamma = max(gamma, Fraction(0))
    if no_gamma:
        gamma = Fraction(0)
    return {"ace": float(s_abun + s_rare / c_ace + singles / c_ace * gamma)}


def _ref_ace(counts, rare_threshold):
    from skbio.diversity.alpha import ace as a
    return {"ace": float(a([int(x) for x in counts if x > 0], rare_threshold=rare_threshold))}


def _rand_ace(rng):
    c = _rand_counts(rng, 8, 30) + [1, 2, 3, rng.randint(4, 9)]
    rng.shuffle(c)
    return {"counts": c, "rare_threshold": 10}


_add("bio_ace", name="ACE species richness estimator", topic="ecology / richness estimation",
     purpose="The abundance-based coverage estimator (ACE) of species richness, with a stated rare-species threshold "
             "(10 is the usual value).",
     inputs=[_COUNTS_IN, _intin("rare_threshold", "largest count of a rare species (usually 10)", lo=1)],
     outputs=[_o("ace", "estimated species richness")], primary="ace",
     conditions=["whole counts", "not every rare species a singleton"],
     method="abundant species counted directly; rare species scaled by the sample coverage C = 1 - F1/N_rare with the "
            "coefficient-of-variation term (EstimateS form, as scikit-bio)",
     latex=r"S_{ace} = S_{abund} + \frac{S_{rare}}{C_{ace}} + \frac{F_1}{C_{ace}}\gamma^2", run=ace, reference=_ref_ace,
     wrong=lambda counts, rare_threshold: ace(counts, rare_threshold, no_gamma=True), rand=_rand_ace,
     phrases=["ACE richness estimator", "ACE"],
     examples=[_ex({"counts": [1, 1, 2, 3, 5, 8, 20, 50], "rare_threshold": 10}, {"ace": 9.719298245614034}),
               _ex({"counts": [12, 15, 20], "rare_threshold": 10}, {"ace": 3.0})],
     level="advanced")


def goods_coverage(counts):
    c = _abund(counts)
    return {"coverage": float(1 - Fraction(sum(1 for x in c if x == 1), sum(c)))}


def _ref_goods(counts):
    from skbio.diversity.alpha import goods_coverage as g
    return {"coverage": float(g([int(x) for x in counts if x > 0]))}


_add("bio_goods_coverage", name="Good's coverage of a sample", topic="ecology / sampling completeness",
     purpose="Good's coverage estimator: the estimated fraction of the community's individuals belonging to species "
             "already seen in the sample.",
     inputs=[_COUNTS_IN], outputs=[_o("coverage", "Good's coverage 1 - F1/N")], primary="coverage",
     conditions=["whole counts"], method="one minus the proportion of individuals that are singletons",
     latex=r"C = 1 - \frac{F_1}{N}", run=goods_coverage, reference=_ref_goods,
     wrong=lambda counts: (lambda c: {"coverage": 1 - sum(1 for x in c if x == 1) / len(c)})([x for x in counts if x > 0]),
     rand=lambda rng: {"counts": _rand_counts(rng, 3, 25)},
     phrases=["Good's coverage", "Goods coverage"],
     examples=[_ex({"counts": [1, 1, 2, 6]}, {"coverage": 0.8}), _ex({"counts": [5, 5]}, {"coverage": 1.0})],
     level="beginner")


def richness_indices(counts, log10=False):
    c = _abund(counts)
    N, S = sum(c), len(c)
    _need(N >= 2, "at least 2 individuals are needed (ln N must be positive)")
    return _out(margalef=(S - 1) / (math.log10(N) if log10 else math.log(N)), menhinick=S / math.sqrt(N))


def _ref_richness(counts):
    from skbio.diversity.alpha import margalef, menhinick
    c = [int(x) for x in counts if x > 0]
    return {"margalef": float(margalef(c)), "menhinick": float(menhinick(c))}


_add("bio_richness_indices", name="Margalef and Menhinick richness indices", topic="ecology / diversity indices",
     purpose="Margalef's richness index (S - 1)/ln N and Menhinick's richness index S/sqrt(N) of a sample.",
     inputs=[_COUNTS_IN], outputs=[_o("margalef", "Margalef's index (S - 1)/ln N"), _o("menhinick", "Menhinick's index S/sqrt(N)")],
     primary="margalef", conditions=["whole counts", "N >= 2"], method="species and individuals counted",
     latex=r"D_{Mg} = \frac{S-1}{\ln N},\quad D_{Mn} = \frac{S}{\sqrt N}", run=richness_indices,
     reference=_ref_richness, wrong=lambda counts: richness_indices(counts, log10=True),
     rand=lambda rng: {"counts": _rand_counts(rng, 2, 25)},
     phrases=["Richness indices", "Margalef and Menhinick indices"],
     examples=[_ex({"counts": [10, 20, 30, 40]}, {"margalef": 0.651441722855, "menhinick": 0.4}),
               _ex({"counts": [1, 1, 1, 1]}, {"margalef": 2.16404256133, "menhinick": 2.0})],
     level="beginner")


def dominance_indices(counts, wrong=False):
    c = _abund(counts)
    N = sum(c)
    _need(N >= 2, "at least 2 individuals are needed")
    U = math.sqrt(sum(x * x for x in c))
    return _out(berger_parker=max(c) / N, mcintosh=(N - U) / (N if wrong else (N - math.sqrt(N))))


def _ref_dominance(counts):
    from skbio.diversity.alpha import berger_parker_d, mcintosh_d
    c = [int(x) for x in counts if x > 0]
    return {"berger_parker": float(berger_parker_d(c)), "mcintosh": float(mcintosh_d(c))}


_add("bio_dominance_indices", name="Berger-Parker and McIntosh indices", topic="ecology / diversity indices",
     purpose="The Berger-Parker dominance (share of the most abundant species) and McIntosh's index of a sample.",
     inputs=[_COUNTS_IN], outputs=[_o("berger_parker", "Berger-Parker dominance n_max/N"),
                                   _o("mcintosh", "McIntosh index (N - U)/(N - sqrt N), U = sqrt(sum n_i^2)")],
     primary="berger_parker", conditions=["whole counts", "N >= 2"], method="largest count over N; McIntosh's U statistic",
     latex=r"d = \frac{n_{\max}}{N},\quad D = \frac{N - U}{N - \sqrt N}", run=dominance_indices,
     reference=_ref_dominance, wrong=lambda counts: dominance_indices(counts, wrong=True),
     rand=lambda rng: {"counts": _rand_counts(rng, 2, 25)},
     phrases=["Berger-Parker and McIntosh indices", "Dominance indices"],
     examples=[_ex({"counts": [50, 30, 20]}, {"berger_parker": 0.5, "mcintosh": 0.426176221892}),
               _ex({"counts": [25, 25, 25, 25]}, {"berger_parker": 0.25})],
     level="beginner")


def brillouin(counts, log10=False):
    c = _abund(counts)
    N = sum(c)
    m = math.factorial(N)
    for x in c:
        m //= math.factorial(x)
    return _out(H_B=(math.log10(m) if log10 else math.log(m)) / N)


def _ref_brillouin(counts):
    from skbio.diversity.alpha import brillouin_d
    return {"H_B": float(brillouin_d([int(x) for x in counts if x > 0]))}


_add("bio_brillouin_diversity", name="Brillouin diversity index", topic="ecology / diversity indices",
     purpose="Brillouin's diversity index of a fully censused collection (natural logarithm).",
     inputs=[_vecin("counts", "number of individuals of each species", lo=1, hi=200)],
     outputs=[_o("H_B", "Brillouin index (nats per individual)")], primary="H_B",
     conditions=["whole counts", "N <= 20000"], method="logarithm of the multinomial coefficient over N (exact integer)",
     latex=r"H_B = \frac{\ln N! - \sum_i \ln n_i!}{N}", run=brillouin, reference=_ref_brillouin,
     wrong=lambda counts: brillouin(counts, log10=True),
     rand=lambda rng: {"counts": _rand_counts(rng, 2, 15)},
     phrases=["Brillouin diversity", "Brillouin index"],
     examples=[_ex({"counts": [5, 3, 2]}, {"H_B": 0.783201418051}), _ex({"counts": [1, 1]}, {"H_B": 0.34657359028})],
     level="intermediate")


def rarefaction(counts, n, with_replacement=False):
    c = _abund(counts)
    N = sum(c)
    _need(isinstance(n, int) and 1 <= n <= N, "n must be a whole number from 1 to the number of individuals")
    _need(N <= 100000, "at most 100000 individuals")
    if with_replacement:
        return {"expected_species": math.fsum(1 - (1 - x / N) ** n for x in c)}
    tot = math.comb(N, n)
    e = sum(1 - Fraction(math.comb(N - x, n), tot) for x in c)
    return {"expected_species": float(e)}


def _ref_rarefaction(counts, n):
    from scipy.stats import hypergeom
    c = [int(x) for x in counts if x > 0]
    N = sum(c)
    return {"expected_species": float(sum(1 - hypergeom.pmf(0, N, x, n) for x in c))}


def _rand_rare(rng):
    c = _rand_counts(rng, 3, 20)
    return {"counts": c, "n": rng.randint(2, max(2, sum(c) - 1))}


_add("bio_rarefaction", name="Rarefied species richness", topic="ecology / richness estimation",
     purpose="The expected number of species in a random subsample of n individuals drawn without replacement "
             "(Hurlbert-Sanders rarefaction).",
     inputs=[_COUNTS_IN, _intin("n", "subsample size (individuals)", lo=1)],
     outputs=[_o("expected_species", "expected number of species in the subsample")], primary="expected_species",
     conditions=["whole counts", "1 <= n <= N"], method="one minus the chance that a species is missed, summed (exact)",
     latex=r"E[S_n] = \sum_i \left[1 - \binom{N-N_i}{n}\Big/\binom{N}{n}\right]", run=rarefaction,
     reference=_ref_rarefaction, wrong=lambda counts, n: rarefaction(counts, n, with_replacement=True), rand=_rand_rare,
     phrases=["Rarefaction", "Rarefied species richness"],
     examples=[_ex({"counts": [2, 1, 1], "n": 2}, {"expected_species": 11 / 6}),
               _ex({"counts": [10, 5, 1], "n": 16}, {"expected_species": 3.0})],
     level="advanced")


def beta_diversity(x, y, similarity=False):
    _need(isinstance(x, list) and isinstance(y, list) and len(x) == len(y) and len(x) >= 1,
          "x and y must list the same species in the same order")
    for v in x + y:
        _fin("abundance", v)
        _need(v >= 0, "abundances must be zero or positive")
    _need(any(v > 0 for v in x) and any(v > 0 for v in y), "each sample needs at least one species")
    bc = math.fsum(abs(a - b) for a, b in zip(x, y)) / math.fsum(a + b for a, b in zip(x, y))
    A = {i for i, v in enumerate(x) if v > 0}
    B = {i for i, v in enumerate(y) if v > 0}
    j = Fraction(len(A ^ B), len(A | B))
    s = Fraction(len(A ^ B), len(A) + len(B))
    if similarity:
        j, s = 1 - j, 1 - s
    return {"bray_curtis": bc, "jaccard": float(j), "sorensen": float(s)}


def _ref_beta(x, y):
    import numpy as np
    from scipy.spatial.distance import braycurtis, dice, jaccard
    a, b = np.array(x, float), np.array(y, float)
    return {"bray_curtis": float(braycurtis(a, b)), "jaccard": float(jaccard(a > 0, b > 0)), "sorensen": float(dice(a > 0, b > 0))}


def _rand_beta(rng):
    k = rng.randint(4, 15)
    x = [float(rng.choice([0, 0] + list(range(1, 60)))) for _ in range(k)]
    y = [float(rng.choice([0, 0] + list(range(1, 60)))) for _ in range(k)]
    x[0], y[1] = 3.0, 4.0
    return {"x": x, "y": y}


_add("bio_beta_diversity", name="Beta diversity between two samples", topic="ecology / community comparison",
     purpose="The Bray-Curtis dissimilarity (abundances) and the Jaccard and Sorensen dissimilarities (presence and "
             "absence) between two community samples.",
     inputs=[_vecin("x", "abundances of each species in the first sample", lo=1, hi=500),
             _vecin("y", "abundances of the same species in the second sample", lo=1, hi=500)],
     outputs=[_o("bray_curtis", "Bray-Curtis dissimilarity"), _o("jaccard", "Jaccard dissimilarity (1 - shared/union)"),
              _o("sorensen", "Sorensen dissimilarity (1 - 2 shared/(S_x + S_y))")],
     primary="bray_curtis", conditions=["the same species order in both samples", "each sample non-empty"],
     method="sum |x - y| / sum (x + y); presence-absence set differences",
     latex=r"BC = \frac{\sum|x_i - y_i|}{\sum (x_i + y_i)},\quad J = 1 - \frac{|A\cap B|}{|A\cup B|}",
     run=beta_diversity, reference=_ref_beta, wrong=lambda x, y: beta_diversity(x, y, similarity=True), rand=_rand_beta,
     phrases=["Beta diversity", "Bray-Curtis and Jaccard dissimilarity"],
     examples=[_ex({"x": [6, 7, 4, 0], "y": [10, 0, 6, 3]}, {"bray_curtis": 0.444444444444, "jaccard": 0.5, "sorensen": 1 / 3}),
               _ex({"x": [1, 2], "y": [1, 2]}, {"bray_curtis": 0.0, "jaccard": 0.0, "sorensen": 0.0})],
     level="intermediate")


def species_area(A, S, no_exp=False):
    _need(isinstance(A, list) and isinstance(S, list) and len(A) == len(S) >= 2, "A and S must be lists of equal length")
    for v in A + S:
        _pos("areas and species counts", v)
    b, a, r2 = _ols([math.log(x) for x in A], [math.log(y) for y in S])
    return _out(c=a if no_exp else math.exp(a), z=b, r2=r2)


def _ref_species_area(A, S):
    import numpy as np
    b, a, r2 = _linregress(np.log(A), np.log(S))
    return {"c": float(np.exp(a)), "z": b, "r2": r2}


def _rand_sar(rng):
    c, z = rng.uniform(2, 40), rng.uniform(0.15, 0.4)
    A = sorted({_sig(_lu(rng, 0.5, 1e5), 3) for _ in range(rng.randint(4, 8))})
    return {"A": A, "S": [float(max(1, round(c * a ** z * rng.uniform(0.85, 1.15)))) for a in A]}


_add("bio_species_area", name="Species-area relationship fit", topic="ecology / biogeography",
     purpose="The constants c and z of the Arrhenius species-area relationship S = c A^z fitted by least squares on "
             "log S against log A.",
     inputs=[_vecin("A", "areas (any one unit)", lo=2), _vecin("S", "species counts", lo=2)],
     outputs=[_o("c", "constant c (species in a unit area)"), _o("z", "exponent z"), _o("r2", "r^2 of the log-log fit")],
     primary="z", conditions=["positive areas and counts", "c refers to the area unit used"],
     method="ordinary least squares of ln S on ln A: slope z, intercept ln c",
     latex=r"\ln S = \ln c + z\ln A", run=species_area, reference=_ref_species_area,
     wrong=lambda A, S: species_area(A, S, no_exp=True), rand=_rand_sar,
     phrases=["Species-area relationship", "Species-area fit"],
     examples=[_ex({"A": [1, 10, 100, 1000], "S": [10, 20, 40, 80]}, {"c": 10.0, "z": 0.301029995664, "r2": 1.0}),
               _ex({"A": [2, 8], "S": [5, 10]}, {"z": 0.5, "c": 3.53553390593})],
     level="beginner")


def petersen_mle(M, C, R, ceil=False):
    for nm, v in (("M", M), ("C", C), ("R", R)):
        _need(isinstance(v, int) and v >= 1, f"{nm} must be a whole number >= 1")
    _need(R <= min(M, C), "R cannot exceed M or C")
    return {"N_mle": -((-M * C) // R) if ceil else (M * C) // R}


def _ref_petersen(M, C, R):
    """Brute force: the hypergeometric likelihood of R recaptures maximised over every whole N (scipy)."""
    import numpy as np
    from scipy.stats import hypergeom
    Ns = np.arange(M + C - R, 3 * (M * C) // R + 50)
    ll = hypergeom.logpmf(R, Ns, M, C)
    best = ll.max()
    return {"N_mle": int(Ns[ll >= best - 1e-9].max())}


def _rand_petersen(rng):
    N = rng.randint(50, 3000)
    M = rng.randint(10, max(11, N // 3))
    C = rng.randint(10, max(11, N // 3))
    R = max(1, min(M, C, round(M * C / N * rng.uniform(0.7, 1.3))))
    return {"M": M, "C": C, "R": R}


_add("bio_petersen_estimate", name="Mark-recapture population estimate (Petersen MLE)",
     topic="ecology / population estimation",
     purpose="The maximum-likelihood population size from one mark-recapture experiment: M animals marked, C caught "
             "later, R of them marked (hypergeometric likelihood).",
     inputs=[_intin("M", "number marked and released", lo=1), _intin("C", "number caught in the second sample", lo=1),
             _intin("R", "number of marked animals recaptured", lo=1)],
     outputs=[_o("N_mle", "maximum-likelihood population size (largest maximiser)")], primary="N_mle",
     conditions=["closed population", "R >= 1, R <= min(M, C)", "when MC/R is whole, N = MC/R and MC/R - 1 tie; the "
                                                                "larger is given"],
     method="the hypergeometric likelihood increases while N <= MC/R: the maximiser is floor(MC/R)",
     latex=r"\hat N = \left\lfloor \frac{MC}{R} \right\rfloor", run=petersen_mle, reference=_ref_petersen,
     wrong=lambda M, C, R: petersen_mle(M, C, R, ceil=True), rand=_rand_petersen,
     phrases=["Petersen estimate", "Mark-recapture estimate"],
     examples=[_ex({"M": 100, "C": 80, "R": 12}, {"N_mle": 666}), _ex({"M": 50, "C": 40, "R": 10}, {"N_mle": 200})],
     extra_q=[{"q": "Petersen estimate: M = 10, C = 5, R = 6", "outcome": "INVALID_INPUT"}],
     level="beginner")


def schumacher(C, R, M, schnabel=False):
    _need(isinstance(C, list) and isinstance(R, list) and isinstance(M, list) and len(C) == len(R) == len(M) >= 2,
          "C, R and M need one value per sampling occasion (at least 2)")
    for c, r, m in zip(C, R, M):
        _need(float(c) == int(c) and c >= 1, "catches must be whole numbers >= 1")
        _need(float(r) == int(r) and 0 <= r <= c, "recaptures must be whole numbers from 0 to the catch")
        _need(float(m) == int(m) and m >= 0, "marked animals at large must be whole numbers >= 0")
    srm = sum(Fraction(int(r) * int(m)) for r, m in zip(R, M))
    _need(srm > 0, "no marked animal was ever recaptured: the population cannot be estimated")
    if schnabel:
        return {"N": float(sum(Fraction(int(c) * int(m)) for c, m in zip(C, M)) / sum(Fraction(int(r)) for r in R))}
    return {"N": float(sum(Fraction(int(c) * int(m) ** 2) for c, m in zip(C, M)) / srm)}


def _ref_schumacher(C, R, M):
    """Weighted least-squares line through the origin of R/C on M (weights C), numpy lstsq: slope 1/N."""
    import numpy as np
    C, R, M = (np.array(v, float) for v in (C, R, M))
    w = np.sqrt(C)
    b = np.linalg.lstsq((w * M)[:, None], w * R / C, rcond=None)[0][0]
    return {"N": float(1 / b)}


def _rand_schumacher(rng):
    N = rng.randint(200, 3000)
    k = rng.randint(3, 7)
    marked = 0
    C, R, M = [], [], []
    for _ in range(k):
        c = rng.randint(20, max(21, N // 5))
        r = min(c, round(c * marked / N * rng.uniform(0.8, 1.2)))
        C.append(float(c))
        R.append(float(r))
        M.append(float(marked))
        marked += c - r
    if not any(r * m for r, m in zip(R, M)):
        R[-1] = 1.0
    return {"C": C, "R": R, "M": M}


_add("bio_schumacher_eschmeyer", name="Schumacher-Eschmeyer multiple-census estimate",
     topic="ecology / population estimation",
     purpose="The Schumacher-Eschmeyer estimate of a closed population's size from repeated mark-recapture samples.",
     inputs=[_vecin("C", "number caught on each occasion", lo=2), _vecin("R", "number of those already marked", lo=2),
             _vecin("M", "number of marked animals at large before each occasion", lo=2)],
     outputs=[_o("N", "estimated population size")], primary="N",
     conditions=["closed population", "whole counts, R <= C", "at least one recapture of a marked animal"],
     method="weighted regression of R/C on M through the origin", latex=r"\hat N = \frac{\sum_t C_tM_t^2}{\sum_t R_tM_t}",
     run=schumacher, reference=_ref_schumacher, wrong=lambda C, R, M: schumacher(C, R, M, schnabel=True),
     rand=_rand_schumacher, phrases=["Schumacher-Eschmeyer estimate", "Schumacher-Eschmeyer population estimate"],
     examples=[_ex({"C": [10, 10, 10], "R": [0, 2, 3], "M": [0, 10, 18]}, {"N": 4240 / 74}),
               _ex({"C": [50, 40], "R": [0, 10], "M": [0, 50]}, {"N": 200.0})],
     level="advanced")


def chemostat(mu_max, Ks, Y, D, S0, no_residual=False):
    _pos("mu_max, Ks, Y, D and S0", mu_max, Ks, Y, D, S0)
    _need(D < mu_max, "the dilution rate must be below mu_max (otherwise the culture washes out)")
    Ss = Ks * D / (mu_max - D)
    _need(Ss < S0, "washout: the steady-state substrate would exceed the feed concentration")
    return _out(S=Ss, X=Y * S0 if no_residual else Y * (S0 - Ss))


def _ref_chemostat(mu_max, Ks, Y, D, S0):
    """ODE integration (scipy DOP853) of the Monod chemostat from a small inoculum until the state stops changing."""
    f = lambda _t, y: [D * (S0 - y[0]) - mu_max * y[0] / (Ks + y[0]) * y[1] / Y,
                       (mu_max * y[0] / (Ks + y[0]) - D) * y[1]]
    y = [S0, 0.01 * Y * S0]
    for _ in range(400):
        y2 = _ivp(f, 40.0 / D, y)
        done = all(abs(a - b) <= 1e-14 * max(abs(a), 1e-300) for a, b in zip(y, y2))
        y = y2
        if done:
            break
    return {"S": y[0], "X": y[1]}


def _rand_chemostat(rng):
    mu, Ks, Y, S0 = round(rng.uniform(0.3, 1.5), 3), round(rng.uniform(0.05, 2), 3), round(rng.uniform(0.3, 0.6), 3), round(rng.uniform(5, 20), 2)
    Dc = mu * S0 / (Ks + S0)
    return {"mu_max": mu, "Ks": Ks, "Y": Y, "D": round(Dc * rng.uniform(0.2, 0.85), 4), "S0": S0}


_add("bio_chemostat_steady_state", name="Chemostat steady state (Monod kinetics)", topic="microbial growth / continuous culture",
     purpose="The steady-state substrate and biomass concentrations of a chemostat with Monod growth kinetics.",
     inputs=[_numin("mu_max", "maximum specific growth rate (1/time)"), _numin("Ks", "half-saturation constant (concentration)"),
             _numin("Y", "biomass yield (biomass per substrate)"), _numin("D", "dilution rate (1/time, same unit as mu_max)"),
             _numin("S0", "feed substrate concentration")],
     outputs=[_o("S", "steady-state substrate concentration"), _o("X", "steady-state biomass concentration")],
     primary="X", conditions=["D < mu_max", "no washout (S < S0)", "consistent units"],
     method="growth balances dilution: mu(S) = D gives S; the substrate balance gives X",
     latex=r"\bar S = \frac{K_sD}{\mu_{\max} - D},\quad \bar X = Y(S_0 - \bar S)", run=chemostat,
     reference=_ref_chemostat, wrong=lambda mu_max, Ks, Y, D, S0: chemostat(mu_max, Ks, Y, D, S0, no_residual=True),
     rand=_rand_chemostat, phrases=["Chemostat steady state", "Chemostat"],
     examples=[_ex({"mu_max": 1, "Ks": 1, "Y": 0.5, "D": 0.5, "S0": 10}, {"S": 1.0, "X": 4.5}),
               _ex({"mu_max": 0.6, "Ks": 0.2, "Y": 0.4, "D": 0.3, "S0": 5}, {"S": 0.2, "X": 1.92})],
     extra_q=[{"q": "Chemostat steady state: mu_max = 1, Ks = 1, Y = 0.5, D = 0.95, S0 = 10", "outcome": "INVALID_INPUT"}],
     level="advanced")


def _simplex_max(A, b, c):
    """Maximise c.x subject to A x = b, x >= 0 (exact rationals, two-phase simplex, Bland's rule).
    Returns the optimal value, or raises ValueError when infeasible or unbounded."""
    m, n = len(A), len(A[0])
    A = [[Fraction(x) for x in row] for row in A]
    b = [Fraction(x) for x in b]
    for i in range(m):
        if b[i] < 0:
            A[i] = [-x for x in A[i]]
            b[i] = -b[i]
    T = [A[i] + [Fraction(1 if k == i else 0) for k in range(m)] + [b[i]] for i in range(m)]
    basis = [n + i for i in range(m)]
    N = n + m

    def pivot(r, col):
        pv = T[r][col]
        T[r] = [x / pv for x in T[r]]
        for i in range(m):
            if i != r and T[i][col] != 0:
                f = T[i][col]
                T[i] = [x - f * y for x, y in zip(T[i], T[r])]
        basis[r] = col

    def run(cost, allowed):
        while True:
            red = [cost[j] - sum(cost[basis[i]] * T[i][j] for i in range(m)) for j in range(N)]
            enter = next((j for j in range(N) if allowed(j) and red[j] > 0), None)
            if enter is None:
                return
            rows = [(T[i][-1] / T[i][enter], basis[i], i) for i in range(m) if T[i][enter] > 0]
            if not rows:
                raise ValueError("the objective is unbounded")
            r = min(rows)[2]
            pivot(r, enter)
    run([Fraction(0)] * n + [Fraction(-1)] * m, lambda j: True)
    if sum(T[i][-1] for i in range(m) if basis[i] >= n) != 0:
        raise ValueError("no flux distribution satisfies the steady state and the bounds (infeasible)")
    for i in range(m):                                          # drive artificial variables out of the basis
        if basis[i] >= n:
            col = next((j for j in range(n) if T[i][j] != 0), None)
            if col is not None:
                pivot(i, col)
    cost = [Fraction(x) for x in c] + [Fraction(0)] * m
    run(cost, lambda j: j < n)
    x = [Fraction(0)] * N
    for i in range(m):
        x[basis[i]] = T[i][-1]
    return sum(cf * xv for cf, xv in zip(c, x[:n]))


def fba(S, lb, ub, c, box_only=False):
    _need(isinstance(S, list) and S and all(len(r) == len(S[0]) for r in S), "S must be a stoichiometric matrix")
    m, n = len(S), len(S[0])
    _need(m <= 15 and n <= 20, "at most 15 metabolites and 20 reactions")
    for name, v in (("lb", lb), ("ub", ub), ("c", c)):
        _need(isinstance(v, list) and len(v) == n, f"{name} needs one value per reaction ({n})")
    for row in S:
        for x in row:
            _fin("S", x)
    for a, bnd, w in zip(lb, ub, c):
        _fin("bounds and objective", a, bnd, w)
        _need(a <= bnd, "every lower bound must not exceed its upper bound")
    F = lambda x: Fraction(repr(float(x)))
    if box_only:
        return {"objective": float(sum(max(F(w) * F(a), F(w) * F(bnd)) for a, bnd, w in zip(lb, ub, c)))}
    L, U = [F(a) for a in lb], [F(bnd) for bnd in ub]
    Sf = [[F(x) for x in row] for row in S]
    # v = L + x, 0 <= x <= U - L, slack s: x + s = U - L; S x = -S L
    A = [row + [Fraction(0)] * n for row in Sf]
    bvec = [-sum(row[j] * L[j] for j in range(n)) for row in Sf]
    for j in range(n):
        A.append([Fraction(1 if k == j else 0) for k in range(n)] + [Fraction(1 if k == j else 0) for k in range(n)])
        bvec.append(U[j] - L[j])
    cvec = [F(w) for w in c] + [Fraction(0)] * n
    z = _simplex_max(A, bvec, cvec) + sum(F(w) * l for w, l in zip(c, L))
    return {"objective": float(z)}


def _ref_fba(S, lb, ub, c):
    import numpy as np
    from scipy.optimize import linprog
    res = linprog(-np.array(c, float), A_eq=np.array(S, float), b_eq=np.zeros(len(S)), bounds=list(zip(lb, ub)),
                  method="highs", options={"primal_feasibility_tolerance": 1e-10, "dual_feasibility_tolerance": 1e-10})
    if res.status != 0:
        raise RuntimeError(res.message)
    return {"objective": float(-res.fun)}


def _rand_fba(rng):
    m, n = rng.randint(2, 5), rng.randint(4, 8)
    S = [[float(rng.choice([0, 0, 0, 1, -1, 2, -2])) for _ in range(n)] for _ in range(m)]
    for i in range(m):
        if not any(S[i]):
            S[i][rng.randrange(n)] = 1.0
    lb = [float(rng.choice([-10, -5, 0, 0])) for _ in range(n)]
    ub = [float(rng.choice([5, 10, 20])) for _ in range(n)]
    c = [0.0] * n
    c[rng.randrange(n)] = 1.0
    if rng.random() < 0.3:
        c[rng.randrange(n)] = float(rng.choice([0.5, 2, -1]))
    return {"S": S, "lb": lb, "ub": ub, "c": c}


_add("bio_flux_balance_analysis", name="Flux balance analysis (optimal objective)", topic="systems biology / metabolic networks",
     purpose="The maximum of a linear objective (e.g. a biomass or product flux) over all steady-state flux "
             "distributions of a metabolic network within flux bounds: flux balance analysis.",
     inputs=[_matin("S", "stoichiometric matrix (metabolites x reactions)"),
             _vecin("lb", "lower bound of each reaction flux", lo=1, hi=20),
             _vecin("ub", "upper bound of each reaction flux", lo=1, hi=20),
             _vecin("c", "objective coefficient of each reaction", lo=1, hi=20)],
     outputs=[_o("objective", "optimal objective value")], primary="objective",
     conditions=["steady state S v = 0", "finite bounds lb <= ub", "a feasible flux distribution exists"],
     method="linear programming solved exactly in rational arithmetic (two-phase simplex, Bland's rule); the optimal "
            "value is unique even when the optimal fluxes are not",
     latex=r"\max_v\ c^\top v \ \text{ s.t. } Sv = 0,\ \ l \le v \le u", run=fba, reference=_ref_fba,
     wrong=lambda S, lb, ub, c: fba(S, lb, ub, c, box_only=True), rand=_rand_fba,
     phrases=["Flux balance analysis", "FBA"],
     examples=[_ex({"S": [[1, -1, 0, 0], [0, 1, -1, -1]], "lb": [0, 0, 0, 0], "ub": [10, 8, 100, 3], "c": [0, 0, 1, 0]},
                   {"objective": 8.0}),
               _ex({"S": [[1, -1, -1]], "lb": [0, 0, 0], "ub": [10, 100, 100], "c": [0, 1, 2]}, {"objective": 20.0})],
     level="expert")


def faith_pd(tree, taxa, unrooted=False):
    root, leaves = _parse_newick(tree)
    names = taxa.split() if isinstance(taxa, str) else []
    _need(names and len(set(names)) == len(names), "give the observed taxa as distinct names separated by spaces")
    for x in names:
        _need(x in leaves, f"{x} is not a leaf of the tree")
    want = set(names)
    total = Fraction(0)
    keep = []

    def walk(t, is_root):
        nonlocal total
        if not t[2]:
            hit = t[0] in want
        else:
            hits = [walk(k, False) for k in t[2]]
            hit = any(hits)
            keep.append((t, sum(hits)))
        if hit and not is_root:
            _need(t[1] is not None, "every branch needs a length")
            total += Fraction(repr(t[1]))
        return hit
    walk(root, True)
    if unrooted:                                              # drop the path from the taxa's common ancestor to the root
        node = root
        while node[2]:
            kids = [k for k in node[2] if _has_any(k, want)]
            if len(kids) != 1:
                break
            node = kids[0]
            if node[2] or node[0] not in want:
                total -= Fraction(repr(node[1]))
            else:
                break
    return {"pd": float(total)}


def _has_any(t, want):
    if not t[2]:
        return t[0] in want
    return any(_has_any(k, want) for k in t[2])


def _ref_faith(tree, taxa):
    import io
    from skbio import TreeNode
    from skbio.diversity.alpha import faith_pd as fp
    names = taxa.split()
    t = TreeNode.read(io.StringIO(tree.strip().rstrip(";") + ";"))
    return {"pd": float(fp([1] * len(names), taxa=names, tree=t))}


def _rand_faith(rng):
    taxa = [f"T{i}" for i in range(rng.randint(4, 8))]
    nw = _rand_newick(rng, taxa, lengths=True)
    root, _ = _parse_newick(nw)
    clades = []

    def walk(t, top):
        if not t[2]:
            return [t[0]]
        got = [x for k in t[2] for x in walk(k, False)]
        if not top and len(got) >= 2:
            clades.append(got)
        return got
    walk(root, True)
    pool = rng.choice(clades) if clades and rng.random() < 0.8 else taxa
    return {"tree": nw, "taxa": " ".join(rng.sample(pool, rng.randint(2, len(pool))))}


_add("bio_faith_pd", name="Faith's phylogenetic diversity", topic="ecology / phylogenetic diversity",
     purpose="Faith's phylogenetic diversity of a set of observed taxa: the total branch length of the rooted tree "
             "connecting them to the root.",
     inputs=[_textin("tree", "rooted tree in Newick format with branch lengths", "newick", hi=20000),
             _textin("taxa", "observed taxa (leaf names separated by spaces)", "names", hi=20000)],
     outputs=[_o("pd", "Faith's PD (sum of branch lengths)")], primary="pd",
     conditions=["every branch has a length", "taxa are leaves of the tree", "rooted convention (path to the root included)"],
     method="sum the lengths of every branch on a path from an observed taxon to the root, each once",
     latex=r"PD = \sum_{b \in T(\text{taxa}) \cup R} \ell_b", run=faith_pd, reference=_ref_faith,
     wrong=lambda tree, taxa: faith_pd(tree, taxa, unrooted=True), rand=_rand_faith,
     phrases=["Faith's phylogenetic diversity", "Faith PD"],
     examples=[_ex({"tree": "((A:1,B:2):3,(C:4,D:5):6);", "taxa": "A B"}, {"pd": 6.0}),
               _ex({"tree": "((A:1,B:2):3,(C:4,D:5):6);", "taxa": "A C"}, {"pd": 14.0})],
     level="advanced")


# ==== END OF PROCEDURES
