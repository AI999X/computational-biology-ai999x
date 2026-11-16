# Computational Biology (computational-biology-ai999x)

[![tests](https://github.com/AI999X/computational-biology-ai999x/actions/workflows/ci.yml/badge.svg)](https://github.com/AI999X/computational-biology-ai999x/actions/workflows/ci.yml) [![PyPI](https://img.shields.io/pypi/v/computational-biology-ai999x)](https://pypi.org/project/computational-biology-ai999x/) [![Hugging Face](https://img.shields.io/badge/%F0%9F%A4%97-Hugging%20Face-yellow)](https://huggingface.co/AI999X/computational-biology-ai999x)

Computational Biology covers DNA and protein sequences, alignment, molecular weights, genetics and population genetics, phylogenetics, enzyme kinetics, metabolic and population models, and ecology: 86 calculations in 10 topics, each with its formula, its steps and its units, and a plain reason instead of a number when an input is missing or out of range. One compiled engine, four ways to use it: the desktop app, the command line, the browser and Python. Everything runs on your own computer, offline.

[![The Computational Biology desktop app: GC content of a nucleic acid, with its worked example and live scene](https://www.ai999x.com/?ax_shot=computational-biology-ai999x/hero.webp)](https://www.ai999x.com/models/computational-biology-ai999x/)

## Four ways to use it

| | |
|---|---|
| **Desktop app** (Windows, Linux) | [Windows](https://github.com/AI999X/ai999x-desktop/releases/download/v1.6.4/ai999x-1.6.4-windows-x86_64.zip) · [Linux](https://github.com/AI999X/ai999x-desktop/releases/download/v1.6.4/ai999x-1.6.4-linux-x86_64.tar.gz) |
| **Command line** `computational-biology` (Windows, Linux) | [Windows](https://github.com/AI999X/ai999x-desktop/releases/download/v1.6.4/ai999x-cli-1.6.4-windows-x86_64.zip) · [Linux](https://github.com/AI999X/ai999x-desktop/releases/download/v1.6.4/ai999x-cli-1.6.4-linux-x86_64.tar.gz) |
| **Browser** (any device) | [Try it on the model page](https://www.ai999x.com/models/computational-biology-ai999x/#try) |
| **Python** 3.10+ | `pip install computational-biology-ai999x`: on Windows and Linux (x86-64) `pip` installs the same compiled engine; on other systems (macOS, ARM) it installs a pure-Python build that gives the same answers |

## Ask in words: 15 examples

Each question below was asked on this release; the comment under it is the line it printed.

```python
import computational_biology_ai999x as computational_biology

a = computational_biology.ask(
    "Global alignment score: seq1 = GATTACA, seq2 = GCATGCU, match = "
    "1, mismatch = -1, gap = -1"
)
print(a.outcome, a.values)
# SOLVED {'score': 0.0}

a = computational_biology.ask(
    "ACE richness estimator: counts = [1, 1, 2, 3, 5, 8, 20, 50], "
    "rare_threshold = 10"
)
print(a.outcome, a.values)
# SOLVED {'ace': 9.719298245614034}

a = computational_biology.ask(
    "Pedigree inbreeding coefficient: pedigree = [P-S, Q-S, P-D, Q-D, "
    "S-X, D-X], individual = X"
)
print(a.outcome, a.values)
# SOLVED {'F': 0.25}

a = computational_biology.ask(
    "Chemostat steady state: mu_max = 1, Ks = 1, Y = 0.5, D = 0.5, S0 "
    "= 10"
)
print(a.outcome, a.values)
# SOLVED {'S': 1.0, 'X': 4.5}

a = computational_biology.ask(
    "Jukes-Cantor distance: seq1 = AAAAAAAAAA, seq2 = AAAAAAAACC"
)
print(a.outcome, a.values)
# SOLVED {'p': 0.2, 'd': 0.232616196228}

a = computational_biology.ask(
    "Haploid selection: p0 = 0.5, s = 1, t = 1"
)
print(a.outcome, a.values)
# SOLVED {'p_t': 0.666666666667}

a = computational_biology.ask(
    "Flux balance analysis: S = [[1, -1, 0, 0], [0, 1, -1, -1]], lb = "
    "[0, 0, 0, 0], ub = [10, 8, 100, 3], c = [0, 0, 1, 0]"
)
print(a.outcome, a.values)
# SOLVED {'objective': 8.0}

a = computational_biology.ask(
    "Global alignment score with affine gaps: seq1 = ACGTTTACG, seq2 = "
    "ACGACG, match = 2, mismatch = -1, open = -5, extend = -1"
)
print(a.outcome, a.values)
# SOLVED {'score': 5.0}

a = computational_biology.ask(
    "Berger-Parker and McIntosh indices: counts = [50, 30, 20]"
)
print(a.outcome, a.values)
# SOLVED {'berger_parker': 0.5, 'mcintosh': 0.426176221892}

a = computational_biology.ask("Map distance: r = 0.2, map = haldane")
print(a.outcome, a.values)
# SOLVED {'d_cM': 25.5412811883}

a = computational_biology.ask(
    "Kimura distance: seq1 = AAAAAAAAAA, seq2 = AAAAAAAGGC"
)
print(a.outcome, a.values)
# SOLVED {'d': 0.402359478109}

a = computational_biology.ask(
    "Recessive lethal allele frequency: q0 = 0.5, t = 1"
)
print(a.outcome, a.values)
# SOLVED {'q_t': 0.333333333333}

a = computational_biology.ask("Aromaticity: seq = MFWYAK")
print(a.outcome, a.values)
# SOLVED {'aromaticity': 0.5, 'aromatic_count': 3}

a = computational_biology.ask(
    "Global protein alignment score: seq1 = HEAGAWGHEE, seq2 = "
    "PAWHEAE, open = -10, extend = -1"
)
print(a.outcome, a.values)
# SOLVED {'score': 3.0}

a = computational_biology.ask(
    "Effective population size with unequal sex ratio: N_m = 10, N_f = "
    "90"
)
print(a.outcome, a.values)
# SOLVED {'Ne': 36.0, 'Ne_over_N': 0.36}
```

## Screenshots

[![The Computational Biology desktop app: GC content of a nucleic acid, with its worked example and live scene](https://www.ai999x.com/?ax_shot=computational-biology-ai999x/hero-dark.webp)](https://www.ai999x.com/models/computational-biology-ai999x/)

[![Decline of a recessive lethal allele: the live scene in the Computational Biology desktop app](https://www.ai999x.com/?ax_shot=computational-biology-ai999x/see-1.webp)](https://www.ai999x.com/models/computational-biology-ai999x/)

*Decline of a recessive lethal allele*

[![ACE species richness estimator: the live scene in the Computational Biology desktop app](https://www.ai999x.com/?ax_shot=computational-biology-ai999x/see-2.webp)](https://www.ai999x.com/models/computational-biology-ai999x/)

*ACE species richness estimator*

[![Lineweaver-Burk fit of Michaelis-Menten data: the live scene in the Computational Biology desktop app](https://www.ai999x.com/?ax_shot=computational-biology-ai999x/see-3.webp)](https://www.ai999x.com/models/computational-biology-ai999x/)

*Lineweaver-Burk fit of Michaelis-Menten data*

Every calculation, with its formula, a worked example and its figures: **[Computational Biology on ai999x.com](https://www.ai999x.com/models/computational-biology-ai999x/)**. Install guides for the app and the command line: [https://www.ai999x.com/desktop-cli/](https://www.ai999x.com/desktop-cli/).

Version 1.6.4 · License: Apache-2.0 · Author: Ibnul Jaif Farabi · Organization: AI999X · Contact: support@ai999x.com · [Support the project](https://buymeacoffee.com/ibnuljaiffarabi)

*The Future Is Self-Improving.*

## Source and releases

This repository holds the installable source of each release, exactly as published on [PyPI](https://pypi.org/project/computational-biology-ai999x/) (`pip install .` builds it). Releases are tagged `v1.6.4`; the [GitHub release](https://github.com/AI999X/computational-biology-ai999x/releases/tag/v1.6.4) carries the wheel, the desktop app and the command-line program for Windows and Linux, and their checksums.

## Build the program from source

`rust/` holds the Rust source of the `computational-biology` command-line program and browser module: the AI999X engine (`rust/ax-engine`), this model (`rust/computational-biology`) and its data (`rust/models/computational_biology/pack.json`). With Rust 1.80 or later:

```bash
cd rust
cargo build --release -p computational-biology-ai999x
./target/release/computational-biology --help
cargo test -p computational-biology-ai999x
```
