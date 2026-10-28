# computational-biology-ai999x — offline, verified computational biology solver with step-by-step answers

[![tests](https://github.com/AI999X/computational-biology-ai999x/actions/workflows/ci.yml/badge.svg)](https://github.com/AI999X/computational-biology-ai999x/actions/workflows/ci.yml) [![PyPI](https://img.shields.io/pypi/v/computational-biology-ai999x)](https://pypi.org/project/computational-biology-ai999x/) [![Hugging Face](https://img.shields.io/badge/%F0%9F%A4%97-Hugging%20Face-yellow)](https://huggingface.co/AI999X/computational-biology-ai999x)

**[Install](https://pypi.org/project/computational-biology-ai999x/)** · **[Capability guide](https://huggingface.co/AI999X/computational-biology-ai999x/tree/main/docs)** · **[Evaluation](https://huggingface.co/AI999X/computational-biology-ai999x/tree/main/evaluation)** · **[AI999X on Hugging Face](https://huggingface.co/AI999X)** · **[Website](https://ai999x.com)** · **[Support: Buy Me a Coffee](https://buymeacoffee.com/ibnuljaiffarabi)**

License: Apache-2.0 · Author: Ibnul Jaif Farabi (creator and engineer) · Organization: AI999X · Version 1.0.3

Field: **Life Sciences** · Tasks: Question Answering, Token Classification, Text Classification, Sentence Similarity, Feature Extraction

computational-biology-ai999x is an offline, deterministic Python solver that answers computational biology questions with verified results from 83 certified capabilities, gives 616 documented answers identically on Windows and Linux, and refuses rather than guesses when an answer cannot be verified.

**computational-biology-ai999x** is a compact, deterministic computational biology engine for Python. It reads a question in plain English (or a typed call), picks a certified computational biology capability, computes the answer exactly, checks it, and returns the result with its steps, units, LaTeX and assumptions — or an honest refusal that says what is missing. It runs fully offline on any CPU, needs no GPU, no API key and no internet, and is built to give **zero wrong answers**: when it is not sure, it refuses instead of guessing.

[![Buy me a coffee](https://img.shields.io/badge/Buy%20me%20a%20coffee-support%20development-FFDD00?logo=buymeacoffee&logoColor=black)](https://buymeacoffee.com/ibnuljaiffarabi) ![License](https://img.shields.io/badge/license-Apache--2.0-blue) ![Python](https://img.shields.io/badge/python-3.10%20%7C%203.11%20%7C%203.12%20%7C%203.13%20%7C%203.14-blue) ![Offline](https://img.shields.io/badge/runs-offline%20on%20CPU-green)

## Highlights

- **Verified, not generated.** Every answer is computed by a certified computational biology capability and checked (substituted back, dimension-checked or compared with an independent computation) before it is returned.
- **Zero wrong answers by design.** Seven explicit outcomes — SOLVED, NEEDS_INPUT, AMBIGUOUS, INVALID_INPUT, OUT_OF_SCOPE, UNRECOGNIZED_FORMAT, INTERNAL_ERROR — never a silent guess. Every one of its 616 documented answers is reproduced exactly on every tested platform.
- **Plain-English questions and a guaranteed typed interface.** `ask()` reads documented question forms; `solve(capability_id, **inputs)` always solves a supported task with valid inputs.
- **Explains itself.** Steps, formulas in LaTeX, units, assumptions and the model used, ready for notebooks, tutoring apps and homework checkers.
- **Small and portable.** Pure Python (sympy, mpmath); identical answers on Windows and Linux (Python 3.10, 3.11, 3.12, 3.13, 3.14).
- **Coverage:** DNA and protein sequences, pairwise alignment, molecular masses, Mendelian and population genetics, phylogenetics, enzyme kinetics, population models and ecology.

## Model overview

| Property | Value |
|---|---|
| Package | `computational-biology-ai999x` 1.0.3 (import `computational_biology_ai999x`) |
| Type | Deterministic neuro-symbolic computational biology engine (rule-based reading + certified symbolic capabilities) |
| Capabilities | 83 certified procedures |
| Knowledge pack | version 1.0.1 |
| Golden answers | 616 documented answers, identical on every tested platform |
| Platforms tested | Windows (Python 3.13), Linux x86-64 (Python 3.10, 3.11, 3.12, 3.13, 3.14) |
| Runtime | CPU only, offline; no GPU, no network, no API key |
| License | Apache-2.0 |
| Creator / engineer | Ibnul Jaif Farabi |
| Organization | AI999X |

## Quickstart

```bash
pip install computational-biology-ai999x
```

```python
import computational_biology_ai999x as computational_biology

computational_biology.ask("Reverse complement: seq = ATGCCGTA").show()
# -> SOLVED — revcomp = TACGGCAT  (bio_reverse_complement)
computational_biology.ask("Hardy-Weinberg test: n_AA = 298, n_Aa = 489, n_aa = 213").show()
# -> SOLVED — p = 0.5425, q = 0.4575, exp_AA = 294.306, exp_Aa = 496.387  (bio_hardy_weinberg_test)
computational_biology.ask("Global alignment score: seq1 = GATTACA, seq2 = GCATGCU, match = 1, mismatch = -1, gap = -1").show()
# -> SOLVED — score = 0  (bio_global_alignment_score)
computational_biology.ask("Cross phenotype ratio: parent1 = AaBb, parent2 = AaBb").show()
# -> SOLVED — p_all_dominant = 0.5625, p_all_recessive = 0.0625, n_genotypes = 9, n_phenotypes = 4  (bio_cross_phenotype_ratio)
computational_biology.solve("bio_michaelis_menten_fit", S="0.5, 1, 2, 5", v="2.2, 3.5, 5.1, 6.6").show()
# -> SOLVED — Vmax = 8.50006, Km = 1.40001  (bio_michaelis_menten_fit)
computational_biology.ask("Who created you?").message
# -> computational-biology-ai999x was created and engineered by Ibnul Jaif Farabi at AI999X.
```

Command line: `computational-biology-ai999x "<question>"`, `computational-biology-ai999x --capabilities`. The outputs above were produced by running this exact code on computational-biology-ai999x 1.0.3.

## Find the right capability, and draw it

`computational_biology.search()` finds which certified capability a question is about — across **every** AI999X model, so a question for another model points you to it. It runs **[ai999x-encoder](https://huggingface.co/AI999X/ai999x-encoder)**, AI999X's own capability encoder (4.7M parameters, a 4.9 MB ONNX model, offline; the root `model.safetensors` of this repository): `pip install "computational-biology-ai999x[search]"`. `ask(question, suggest=True)` adds the same pointers to a question it does not answer; an answer itself never changes.

```python
import computational_biology_ai999x as computational_biology

computational_biology.search("Solve x^2 - 5x + 6 = 0", k=2)
# -> Roots of a quadratic equation (math-foundations-ai999x, 0.63); Roots of a cubic equation (math-foundations-ai999x, 0.59)
```

## How every answer is verified

1. **Reading** — the question is turned into typed quantities (value, unit, kind, role) and one asked quantity; every stated number must be used or the question is refused.
2. **Planning** — only certified capabilities may answer; a capability applies only under the conditions the wording states (for example a pendulum, a liquid at rest, boiling).
3. **Computing** — exact symbolic or high-precision numeric evaluation.
4. **Checking** — dimensions, substitution back into the defining equation, and agreement between independent readings; any disagreement becomes a refusal, never a number.

Capabilities are certified before they ship: each must reproduce hand-worked reference cases, reject a deliberately wrong variant (a control), and leave every benchmark answer unchanged.

## Worked examples

Each question was put to `computational-biology-ai999x`; the answer is the one it returned, with the steps it took.

**Reverse complement of a DNA sequence** · _Sequences_

> Reverse complement: seq = ATGCCGTA

→ **revcomp = TACGGCAT**

<details><summary>Steps</summary>

- Reverse complement of a DNA sequence: seq = ATGCCGTA → revcomp = TACGGCAT

</details>

**Hardy-Weinberg chi-square test** · _Population genetics_

> Hardy-Weinberg test: n_AA = 298, n_Aa = 489, n_aa = 213

→ **p = 0.5425, q = 0.4575, exp_AA = 294.30625, exp_Aa = 496.3875, exp_aa = 209.30625, chi2 = 0.2214895876, p_value = 0.6379072638, Ho = 0.489, He = 0.4963875, F = 0.01488252625**

<details><summary>Steps</summary>

- Hardy-Weinberg chi-square test: n_AA = 298, n_Aa = 489, n_aa = 213 → p = 0.5425, q = 0.4575, exp_AA = 294.30625, exp_Aa = 496.3875, exp_aa = 209.30625, chi2 = 0.2214895876, p_value = 0.6379072638, Ho = 0.489, He = 0.4963875, F = 0.01488252625

</details>

**Global alignment score (Needleman-Wunsch, linear gaps)** · _Alignment_

> Global alignment score: seq1 = GATTACA, seq2 = GCATGCU, match = 1, mismatch = -1, gap = -1

→ **score = 0**

<details><summary>Steps</summary>

- Global alignment score (Needleman-Wunsch, linear gaps): seq1 = GATTACA, seq2 = GCATGCU, match = 1, mismatch = -1, gap = -1 → score = 0

</details>

**Phenotype ratio of a cross (complete dominance)** · _Genetics_

> Cross phenotype ratio: parent1 = AaBb, parent2 = AaBb

→ **p_all_dominant = 0.5625, p_all_recessive = 0.0625, n_genotypes = 9, n_phenotypes = 4, ratio = 9:3:3:1**

<details><summary>Steps</summary>

- Phenotype ratio of a cross (complete dominance): parent1 = AaBb, parent2 = AaBb → p_all_dominant = 0.5625, p_all_recessive = 0.0625, n_genotypes = 9, n_phenotypes = 4, ratio = 9:3:3:1

</details>

**Michaelis-Menten nonlinear least-squares fit** · _Enzyme kinetics_

`m.solve("bio_michaelis_menten_fit", S="0.5, 1, 2, 5", v="2.2, 3.5, 5.1, 6.6")`

→ **Vmax = 8.500055206, Km = 1.400014367**

<details><summary>Steps</summary>

- Michaelis-Menten nonlinear least-squares fit: S = [0.5, 1, 2, 5], v = [2.2, 3.5, 5.1, 6.6] → Vmax = 8.500055206, Km = 1.400014367

</details>

**ACE species richness estimator** · _Ecology_

> ACE richness estimator: counts = [1, 1, 2, 3, 5, 8, 20, 50], rare_threshold = 10

→ **ace = 9.719298246**

<details><summary>Steps</summary>

- ACE species richness estimator: counts = [1, 1, 2, 3, 5, 8, 20, 50], rare_threshold = 10 → ace = 9.719298246

</details>

## Figures

Every point was computed by `computational-biology-ai999x`.

**Michaelis–Menten kinetics: measured rates and fitted curve**

![Michaelis–Menten kinetics: measured rates and fitted curve — computed by computational-biology-ai999x](https://huggingface.co/AI999X/computational-biology-ai999x/resolve/main/figures/bio-michaelis-menten-fit-fit.svg)

_Michaelis–Menten fit to 4 measured rates (S = [1, 2, 4, 8]): the package returns Vmax = 10 and Km = 1; the curve drawn is v = Vmax·S/(Km + S) with those two values (the double-reciprocal fit gives Vmax = 10, Km = 1, r² = 1). Marked: S = Km = 1, v = Vmax/2 = 5._

**Enzyme progress curve: substrate and product over time**

![Enzyme progress curve: substrate and product over time — computed by computational-biology-ai999x](https://huggingface.co/AI999X/computational-biology-ai999x/resolve/main/figures/bio-michaelis-menten-progress-kinetics.svg)

_Substrate and product over time for an enzyme reaction with S0 = 10, Vmax = 1 and Km = 2 (integrated Michaelis–Menten law), computed by the package at 101 times. Half the substrate is converted at t = 6.38629; marked: t = 6.38629, S = 5._

**Chemostat steady states against dilution rate**

![Chemostat steady states against dilution rate — computed by computational-biology-ai999x](https://huggingface.co/AI999X/computational-biology-ai999x/resolve/main/figures/bio-chemostat-steady-state-operating.svg)

_Chemostat operating diagram (Monod growth, μmax = 1, Ks = 1, Y = 0.5, feed S0 = 10): steady-state substrate and biomass computed by the package at 90 dilution rates up to washout; biomass productivity D·X peaks near D = 0.7 (2.68333). Marked: D = 0.5, X = 4.5, S = 1._

**Age-structured population growth (Leslie matrix)**

![Age-structured population growth (Leslie matrix) — computed by computational-biology-ai999x](https://huggingface.co/AI999X/computational-biology-ai999x/resolve/main/figures/bio-leslie-projection-population.svg)

_Leslie-matrix projection of an age-structured population (fecundities F = [0, 1.5, 1], survival S = [0.5, 0.8], start n0 = [100, 50, 20]): each bar is the population the package returns at that step, split by age class. Growth settles to λ = 1.06152 per step (r = 0.0596993); at t = 12 the classes are [197.307, 92.9426, 69.9967]; marked: total = 360.247._

**Hardy–Weinberg: observed against expected genotypes**

![Hardy–Weinberg: observed against expected genotypes — computed by computational-biology-ai999x](https://huggingface.co/AI999X/computational-biology-ai999x/resolve/main/figures/bio-hardy-weinberg-test-bars.svg)

_Observed genotype counts AA = 298, Aa = 489, aa = 213 beside the counts expected under Hardy–Weinberg equilibrium, computed by the package from the allele frequencies p = 0.5425, q = 0.4575: 294.306, 496.387, 209.306. χ² = 0.22149, p-value 0.637907; observed heterozygosity 0.489 against expected 0.496387 (F = 0.0148825)._

**Heterozygosity lost by drift: population size against time**

![Heterozygosity lost by drift: population size against time — computed by computational-biology-ai999x](https://huggingface.co/AI999X/computational-biology-ai999x/resolve/main/figures/bio-drift-heterozygosity-heatmap.svg)

_Fraction of heterozygosity retained under random genetic drift, (1 − 1/(2N))^t, for population sizes N = 5, 10, 25, 50, 100, 250 (rows) after 10, 25, 50, 100, 200 generations (columns), each cell computed by the package and rounded to 3 decimals. N = 50 keeps 0.366 after 100 generations; N = 5 keeps 0 after 200._

## Capabilities

83 certified procedures. A few of them:

- **Procedures:** ACE species richness estimator; Allele frequency under haploid (genic) selection; Amino-acid composition of a protein; Aromaticity of a protein; Berger-Parker and McIntosh indices; Beta diversity between two samples; Brillouin diversity index; Chao1 species richness estimator; Chemostat steady state (Monod kinetics); Competitive inhibition global fit; Decline of a recessive lethal allele; Eadie-Hofstee fit of Michaelis-Menten data; Edit (Levenshtein) distance between sequences; Effective population size with unequal sex ratio; Effective size of a fluctuating population; Equilibrium under heterozygote advantage; Expected coalescence times; Faith's phylogenetic diversity; Fixation index FST from subpopulation allele frequencies; Fixation probability in the Moran process; Flux balance analysis (optimal objective); GC content of a nucleic acid; Global alignment score (Needleman-Wunsch, linear gaps); Global alignment score with affine gaps (Gotoh) …

`computational_biology_ai999x.templates()` lists every question form `ask()` reads; the capability guide (PDF, every example executed) is included in this repository.

## Out of scope

- Open-ended explanation, essay writing or tutoring dialogue: it is a solver, not a chatbot.
- Questions outside its certified computational biology capabilities, and any question whose stated numbers it cannot all use — it refuses with the reason instead.
- High-stakes engineering, medical or legal decisions without independent review.

## Limitations

- It answers what its certified capabilities cover (DNA and protein sequences, pairwise alignment, molecular masses, Mendelian and population genetics, phylogenetics, enzyme kinetics, population models and ecology); anything else is refused with a reason.
- `ask()` reads documented question forms and common phrasings; unusual wording may be refused — use `solve()` for guaranteed results.
- Results are as correct as the certified relations and their stated conditions; a certified relation is a checked model, not a claim about every physical situation.

## Citation

```bibtex
@software{computational_biology_ai999x_1_0_3,
  title  = {computational-biology-ai999x: an offline, verified computational biology solver},
  author = {Farabi, Ibnul Jaif}, organization = {AI999X},
  year   = {2026}, version = {1.0.3},
  url    = {https://huggingface.co/AI999X/computational-biology-ai999x}
}
```

## Contact

Questions, bug reports and requests: open a discussion on https://huggingface.co/AI999X/computational-biology-ai999x/discussions or visit https://ai999x.com.

## Support the project

computational-biology-ai999x is free and open source. If it helps you, you can support further development on [Buy Me a Coffee](https://buymeacoffee.com/ibnuljaiffarabi).

Created and engineered by **Ibnul Jaif Farabi** · **[AI999X](https://huggingface.co/AI999X)** · [ai999x.com](https://ai999x.com)

## Source and releases

This repository holds the installable source of each release, exactly as published on [PyPI](https://pypi.org/project/computational-biology-ai999x/) (`pip install .` builds it). Releases are tagged `v<version>`; the wheel is attached to each [GitHub release](https://github.com/AI999X/computational-biology-ai999x/releases). The optional ONNX models (formula graphs and, where measured, a decision model) are on [Hugging Face](https://huggingface.co/AI999X/computational-biology-ai999x/tree/main/onnx).
