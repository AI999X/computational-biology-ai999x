# computational-biology-ai999x: offline, verified computational biology solver with step-by-step answers

[![tests](https://github.com/AI999X/computational-biology-ai999x/actions/workflows/ci.yml/badge.svg)](https://github.com/AI999X/computational-biology-ai999x/actions/workflows/ci.yml) [![PyPI](https://img.shields.io/pypi/v/computational-biology-ai999x)](https://pypi.org/project/computational-biology-ai999x/) [![Hugging Face](https://img.shields.io/badge/%F0%9F%A4%97-Hugging%20Face-yellow)](https://huggingface.co/AI999X/computational-biology-ai999x)

**[Desktop app & command line](https://www.ai999x.com/desktop-cli/)** · **[Compute in your browser](https://www.ai999x.com/models/computational-biology-ai999x/#try)** · **[Install (Python)](https://pypi.org/project/computational-biology-ai999x/)** · **[Guides](https://huggingface.co/AI999X/computational-biology-ai999x/tree/main/docs)** · **[Evaluation](https://huggingface.co/AI999X/computational-biology-ai999x/tree/main/evaluation)** · **[AI999X on Hugging Face](https://huggingface.co/AI999X)** · **[Website](https://ai999x.com)** · **[Support: Buy Me a Coffee](https://buymeacoffee.com/ibnuljaiffarabi)**

License: Apache-2.0 · Author: Ibnul Jaif Farabi (creator and engineer) · Organization: AI999X · Version 1.5.0

Field: **Life Sciences** · Tasks: Question Answering

computational-biology-ai999x is an offline, deterministic computational biology solver written in Rust: a desktop app, a command-line program, a browser module and a Python package. It answers computational biology questions with verified results from 83 certified capabilities, gives 616 documented answers identically on Windows and Linux, and refuses rather than guesses when an answer cannot be verified.

**computational-biology-ai999x** reads a question in plain English (or a typed call), picks a certified computational biology capability, computes the answer exactly, checks it, and returns the result with its steps, units, LaTeX and assumptions, or a refusal that says what is missing. It runs fully offline on any CPU, needs no GPU, no API key and no internet, and is built to give **zero wrong answers**: when it is not sure, it refuses instead of guessing.

[![Buy me a coffee](https://img.shields.io/badge/Buy%20me%20a%20coffee-support%20development-FFDD00?logo=buymeacoffee&logoColor=black)](https://buymeacoffee.com/ibnuljaiffarabi) ![License](https://img.shields.io/badge/license-Apache--2.0-blue) ![Python](https://img.shields.io/badge/python-3.10%20%7C%203.11%20%7C%203.12%20%7C%203.13%20%7C%203.14-blue) ![Offline](https://img.shields.io/badge/runs-offline%20on%20CPU-green) ![Rust](https://img.shields.io/badge/Rust-Windows%20%7C%20Linux%20%7C%20WebAssembly-orange)

**Run it your way.** computational-biology-ai999x comes as the AI999X desktop app and as the command-line program `computational-biology` for Windows and Linux, runs in your browser as WebAssembly, and installs as a Python package. The desktop app, the command line and the browser run the same compiled Rust program on your own device: offline, no Python needed.

## Get it

| Way to run it | For | Get it |
|---|---|---|
| **AI999X Computational Biology desktop app** 1.5.0 | Windows, Linux. This model in its own app: search, ask in words, forms, live figures, explanations, offline | [Windows (x86-64)](https://downloads.ai999x.com/desktop/1.5.0/ai999x-computational-biology-1.5.0-windows-x86_64.zip) · [Linux (x86-64)](https://downloads.ai999x.com/desktop/1.5.0/ai999x-computational-biology-1.5.0-linux-x86_64.tar.gz) |
| **AI999X desktop app** 1.5.0 | Windows, Linux. Every AI999X model in one app: search, ask in words, forms, live figures, explanations, offline | [Windows (x86-64)](https://downloads.ai999x.com/desktop/1.5.0/ai999x-1.5.0-windows-x86_64.zip) · [Linux (x86-64)](https://downloads.ai999x.com/desktop/1.5.0/ai999x-1.5.0-linux-x86_64.tar.gz) |
| **Command line** (`computational-biology`) | Windows, Linux. One program per model, offline, no Python | [Windows (x86-64)](https://downloads.ai999x.com/desktop/1.5.0/ai999x-cli-1.5.0-windows-x86_64.zip) · [Linux (x86-64)](https://downloads.ai999x.com/desktop/1.5.0/ai999x-cli-1.5.0-linux-x86_64.tar.gz) |
| **In your browser** | any device. Nothing to install | [Compute on the model page](https://www.ai999x.com/models/computational-biology-ai999x/#try) |
| **Python** | Python 3.10+ | `pip install computational-biology-ai999x` |

Checksums: [SHA256SUMS](https://downloads.ai999x.com/desktop/1.5.0/SHA256SUMS) · every version: [https://downloads.ai999x.com/desktop/latest/](https://downloads.ai999x.com/desktop/latest/) · install guide: [desktop & command line](https://www.ai999x.com/desktop-cli/)

## Quick start

### Desktop app

Download **AI999X Computational Biology** for your system (above), unpack it and start it; type a question in words or pick a capability, fill in the form and press Compute. Figures and an explanation of the answer appear beside it. The AI999X desktop app has the same, with Computational Biology among all AI999X models.

Linux, from a terminal:

```bash
curl -LO https://downloads.ai999x.com/desktop/1.5.0/ai999x-computational-biology-1.5.0-linux-x86_64.tar.gz
curl -LO https://downloads.ai999x.com/desktop/1.5.0/SHA256SUMS
sha256sum -c --ignore-missing SHA256SUMS
tar xzf ai999x-computational-biology-1.5.0-linux-x86_64.tar.gz
./ai999x-computational-biology-1.5.0-linux-x86_64/ai999x-computational-biology
```

### Command line

Linux:

```bash
curl -LO https://downloads.ai999x.com/desktop/1.5.0/ai999x-cli-1.5.0-linux-x86_64.tar.gz
curl -LO https://downloads.ai999x.com/desktop/1.5.0/SHA256SUMS
sha256sum -c --ignore-missing SHA256SUMS
mkdir -p ~/.local/share && tar xzf ai999x-cli-1.5.0-linux-x86_64.tar.gz -C ~/.local/share
echo 'export PATH="$HOME/.local/share/ai999x-cli-1.5.0-linux-x86_64:$PATH"' >> ~/.bashrc && source ~/.bashrc
computational-biology --version
```

Windows (PowerShell):

```powershell
$F = "ai999x-cli-1.5.0-windows-x86_64.zip"
Invoke-WebRequest "https://downloads.ai999x.com/desktop/1.5.0/$F" -OutFile $F
Invoke-WebRequest "https://downloads.ai999x.com/desktop/1.5.0/SHA256SUMS" -OutFile SHA256SUMS
(Get-FileHash $F -Algorithm SHA256).Hash.ToLower() -eq ((Select-String -SimpleMatch $F SHA256SUMS).Line -split ' ')[0]
Expand-Archive $F -DestinationPath "$env:LOCALAPPDATA\AI999X" -Force
$D = "$env:LOCALAPPDATA\AI999X\ai999x-cli-1.5.0-windows-x86_64"; $env:Path += ";$D"
[Environment]::SetEnvironmentVariable("Path", [Environment]::GetEnvironmentVariable("Path", "User") + ";$D", "User")
computational-biology --version
```

Every command below was run as shown; the text under it is what the program printed.

**Ask in words**: Write the question as you would say it: the program reads the quantities and their units, picks the capability and answers.

```console
$ computational-biology ask "Reverse complement: seq = ATGCCGTA"
Read as: Reverse complement of a DNA sequence
revcomp = TACGGCAT
  by: Reverse complement of a DNA sequence  (complement each base (A-T, C-G) and reverse the order)
  seq = ATGCCGTA
  assumes: DNA letters A, C, G, T only
```

**What a capability needs**: The equation or method, every input and output with its SI unit. Read this before solve: it tells you the input names.

```console
$ computational-biology info bio_hardy_weinberg_test
Hardy-Weinberg chi-square test
  allele counting, expected counts N p^2, 2Npq, N q^2, chi-square with 3 - 1 - 1 = 1 degree of freedom
  in   n_AA     number of AA individuals (integer) [1]
  in   n_Aa     number of Aa individuals (integer) [1]
  in   n_aa     number of aa individuals (integer) [1]
  out  p        frequency of allele A [1]
  out  q        frequency of allele a [1]
  out  exp_AA   expected AA count (N p^2) [1]
  out  exp_Aa   expected Aa count (2 N p q) [1]
  out  exp_aa   expected aa count (N q^2) [1]
  out  chi2     chi-square statistic [1]
  out  p_value  p-value (chi-square, 1 degree of freedom) [1]
  out  Ho       observed heterozygosity [1]
  out  He       expected heterozygosity 2pq [1]
  out  F        inbreeding coefficient 1 - Ho/He [1]
```

**Solve**: Inputs as name=value, each with any unit of the right kind. The answer lists every output, the method used, each input converted to SI, and the conditions the answer assumes.

```console
$ computational-biology solve bio_hardy_weinberg_test n_AA=298 n_Aa=489 n_aa=213
p = 0.5425
q = 0.4575
exp_AA = 294.30625
exp_Aa = 496.3875
exp_aa = 209.30625
chi2 = 0.2214895876
p_value = 0.6379072638
Ho = 0.489
He = 0.4963875
F = 0.01488252625
  by: Hardy-Weinberg chi-square test  (allele counting, expected counts N p^2, 2Npq, N q^2, chi-square with 3 - 1 - 1 = 1 degree of freedom)
  assumes: both alleles present; biallelic locus; no continuity correction
```

**Many questions at once**: One JSON request per line in, one JSON answer per line out, in order: for scripts, notebooks and agents (the program starts once). In PowerShell: Get-Content cases.jsonl | computational-biology batch

`cases.jsonl`:

```json
{"key": "bio_hardy_weinberg_test", "inputs": {"n_AA": "298", "n_Aa": "489", "n_aa": "213"}}
{"key": "bio_reverse_complement", "inputs": {"seq": "ATGCCGTA"}}
```

```console
$ computational-biology batch < cases.jsonl
{"outcome":"SOLVED","message":"","value":0.221489587647,"unit":"","solved_for":"chi2","values":[["p",0.5425,""],["q",0.4575,""],["exp_AA",294.30625,""],["exp_Aa",496.3875,""],["exp_aa",209.30625,""],["chi2",0.221489587647,""],["p_value",0.637907263782,""],["Ho",0.489,""],["He",0.4963875,""],["F",0.0148825262522,""]],"target":"chi2","phrase":"frequency of allele A","assumptions":["both alleles present","biallelic locus","no continuity correction"],"name":"Hardy-Weinberg chi-square test","equation":"allele counting, expected counts N p^2, 2Npq, N q^2, chi-square with 3 - 1 - 1 = 1 degree of freedom"}
{"outcome":"SOLVED","message":"","value":null,"unit":"","solved_for":"revcomp","values":[["revcomp","TACGGCAT",""]],"target":"revcomp","phrase":"reverse complement 5'->3'","assumptions":["DNA letters A, C, G, T only"],"name":"Reverse complement of a DNA sequence","equation":"complement each base (A-T, C-G) and reverse the order"}
```

### In your browser

Open [Computational Biology on ai999x.com](https://www.ai999x.com/models/computational-biology-ai999x/#try), choose a capability and press **Compute**: the model runs on your own device as WebAssembly; nothing is sent anywhere. Change any input (any unit of the right kind) and the answer, its explanation and the figures follow.

In your own page: `computational_biology_ai999x.wasm` and the loader `ax-wasm.js` (https://www.ai999x.com/compute/rust/ax-wasm.js).

### Python

```bash
pip install computational-biology-ai999x
```

## Worked examples in every form

The same worked example from the command line, from Python and in the browser. Each command-line output is what the program printed.

### Global alignment score (Needleman-Wunsch, linear gaps)

The optimal global (Needleman-Wunsch) alignment score of two sequences with stated match, mismatch and linear gap scores.

Command line:

```console
$ computational-biology solve bio_global_alignment_score seq1=GATTACA seq2=GCATGCU match=1 mismatch=-1 gap=-1
score = 0
  by: Global alignment score (Needleman-Wunsch, linear gaps)  (dynamic programming over prefixes (Needleman-Wunsch))
  match = 1
  mismatch = -1
  gap = -1
  seq1 = GATTACA
  seq2 = GCATGCU
  assumes: gap <= 0; end gaps scored like internal gaps
```

Python:

```python
import computational_biology_ai999x as m

a = m.solve("bio_global_alignment_score", seq1="GATTACA", seq2="GCATGCU", match="1", mismatch="-1", gap="-1")
print(a.values or (a.value, a.unit))
```

Browser: [open this capability on the model page](https://www.ai999x.com/models/computational-biology-ai999x/#cap-bio_global_alignment_score) and press **Compute in your browser**.

_What this means:_ This calculation determines the highest possible score for aligning two entire sequences end-to-end. By applying dynamic programming to the prefixes of sequence one and sequence two, it evaluates every possible alignment path. The result, 0, represents the sum of scores for matched characters, mismatched characters, and gap positions, where each gap contributes a fixed penalty defined by the input parameter. This metric quantifies the overall similarity between the two inputs under the specified scoring scheme.

### Global alignment score with affine gaps (Gotoh)

The optimal global alignment score with affine gaps: a gap of length L scores open + (L-1)*extend.

Command line:

```console
$ computational-biology solve bio_global_alignment_affine seq1=ACGTTTACG seq2=ACGACG match=2 mismatch=-1 open=-5 extend=-1
score = 5
  by: Global alignment score with affine gaps (Gotoh)  (Gotoh's three-state dynamic programming)
  match = 2
  mismatch = -1
  open = -5
  extend = -1
  seq1 = ACGTTTACG
  seq2 = ACGACG
  assumes: open, extend <= 0; a gap of length L scores open + (L-1) extend; end gaps scored
```

Python:

```python
import computational_biology_ai999x as m

a = m.solve("bio_global_alignment_affine", seq1="ACGTTTACG", seq2="ACGACG", match="2", mismatch="-1", open="-5", extend="-1")
print(a.values or (a.value, a.unit))
```

Browser: [open this capability on the model page](https://www.ai999x.com/models/computational-biology-ai999x/#cap-bio_global_alignment_affine) and press **Compute in your browser**.

_What this means:_ This calculation determines the highest possible similarity score for aligning two biological sequences globally, treating a new gap opening differently from extending an existing one. By applying Gotoh's dynamic programming method with a match score of 2 and mismatch penalty of -1, the algorithm evaluates all possible alignments while penalizing gaps based on an open cost of -5 and an extension cost of -1. The resulting value of 5 represents the maximum cumulative score achievable under these specific scoring rules, accounting for both substitutions and indels across the entire length of the sequences.

### Global protein alignment score (BLOSUM62, affine gaps)

The optimal global alignment score of two proteins with BLOSUM62 and affine gaps.

Command line:

```console
$ computational-biology solve bio_protein_global_blosum62 seq1=HEAGAWGHEE seq2=PAWHEAE open=-10 extend=-1
score = 3
  by: Global protein alignment score (BLOSUM62, affine gaps)  (Gotoh's dynamic programming with BLOSUM62 substitution scores)
  open = -10
  extend = -1
  seq1 = HEAGAWGHEE
  seq2 = PAWHEAE
  assumes: standard amino acids; open, extend <= 0
```

Python:

```python
import computational_biology_ai999x as m

a = m.solve("bio_protein_global_blosum62", seq1="HEAGAWGHEE", seq2="PAWHEAE", open="-10", extend="-1")
print(a.values or (a.value, a.unit))
```

Browser: [open this capability on the model page](https://www.ai999x.com/models/computational-biology-ai999x/#cap-bio_protein_global_blosum62) and press **Compute in your browser**.

_What this means:_ This calculation determines the best global alignment score between two protein sequences using the BLOSUM62 substitution matrix and affine gap penalties. The result, 3, represents the maximum sum of match scores plus gap costs when aligning HEAGAWGHEE with PAWHEAE from end to end. A higher value indicates greater similarity and a better fit under the chosen scoring scheme.

### Local alignment score (Smith-Waterman, linear gaps)

The optimal local (Smith-Waterman) alignment score of two sequences with stated match, mismatch and linear gap scores.

Command line:

```console
$ computational-biology solve bio_local_alignment_score seq1=TGTTACGG seq2=GGTTGACTA match=3 mismatch=-3 gap=-2
score = 13
  by: Local alignment score (Smith-Waterman, linear gaps)  (dynamic programming with a floor at 0 (Smith-Waterman))
  match = 3
  mismatch = -3
  gap = -2
  seq1 = TGTTACGG
  seq2 = GGTTGACTA
  assumes: match > 0; gap <= 0
```

Python:

```python
import computational_biology_ai999x as m

a = m.solve("bio_local_alignment_score", seq1="TGTTACGG", seq2="GGTTGACTA", match="3", mismatch="-3", gap="-2")
print(a.values or (a.value, a.unit))
```

Browser: [open this capability on the model page](https://www.ai999x.com/models/computational-biology-ai999x/#cap-bio_local_alignment_score) and press **Compute in your browser**.

_What this means:_ This calculation determines the highest possible score for aligning a specific segment from the first sequence with a segment from the second sequence. Using dynamic programming, it evaluates matches, mismatches, and gaps to find the best local similarity. The resulting score of 13 indicates the strength of this local relationship, where higher values suggest more significant biological homology or shared motifs.

## Documentation

- **Rust user guide**: the command line, every capability with a command and what it printed: [PDF](https://huggingface.co/AI999X/computational-biology-ai999x/blob/main/docs/computational_biology_ai999x_rust_guide.pdf) · [Markdown](https://huggingface.co/AI999X/computational-biology-ai999x/blob/main/docs/computational_biology_ai999x_rust_guide.md)
- **Python guide**: every capability with a Python example: [PDF](https://huggingface.co/AI999X/computational-biology-ai999x/blob/main/docs/computational_biology_ai999x_capabilities.pdf) · [Markdown](https://huggingface.co/AI999X/computational-biology-ai999x/blob/main/docs/computational_biology_ai999x_capabilities.md)
- **Desktop app & command line**: downloads, install, ask and search, figures, offline use, checksums: https://www.ai999x.com/desktop-cli/
- **Model page**: compute in your browser, every capability explained: https://www.ai999x.com/models/computational-biology-ai999x/

## Highlights

- **Verified, not generated.** Every answer is computed by a certified computational biology capability and checked (substituted back, dimension-checked or compared with an independent computation) before it is returned.
- **Zero wrong answers by design.** Seven explicit outcomes (SOLVED, NEEDS_INPUT, AMBIGUOUS, INVALID_INPUT, OUT_OF_SCOPE, UNRECOGNIZED_FORMAT, INTERNAL_ERROR), never a silent guess. Every one of its 616 documented answers is reproduced exactly on every tested platform.
- **Plain-English questions and a guaranteed typed interface.** `ask()` reads documented question forms; `solve(capability_id, **inputs)` always solves a supported task with valid inputs.
- **Explains itself.** Steps, formulas in LaTeX, units, assumptions and the model used, ready for notebooks, tutoring apps and homework checkers.
- **Native and in the browser.** Every capability also computes natively in Rust: the desktop app and the `computational-biology` command line for Windows and Linux, and WebAssembly on the model page, no Python, no network.
- **Small and portable.** Pure Python (sympy, mpmath); identical answers on Windows and Linux (Python 3.10, 3.11, 3.12, 3.13, 3.14).
- **Coverage:** DNA and protein sequences, pairwise alignment, molecular masses, Mendelian and population genetics, phylogenetics, enzyme kinetics, population models and ecology.

## Model overview

| Property | Value |
|---|---|
| Package | `computational-biology-ai999x` 1.5.0 (import `computational_biology_ai999x`) |
| Type | Deterministic neuro-symbolic computational biology engine (rule-based reading + certified symbolic capabilities) |
| Capabilities | 83 certified procedures |
| Knowledge pack | version 1.0.1 |
| Golden answers | 616 documented answers, identical on every tested platform |
| Platforms tested | Windows (Python 3.13), Linux x86-64 (Python 3.10, 3.11, 3.12, 3.13, 3.14) |
| Runtime | CPU only, offline; no GPU, no network, no API key |
| License | Apache-2.0 |
| Creator / engineer | Ibnul Jaif Farabi |
| Organization | AI999X |

## Python quick start

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

The Python package also installs a Python command, `computational-biology-ai999x "<question>"`; for the native program use `computational-biology` (above). The outputs above were produced by running this exact code on computational-biology-ai999x 1.5.0.

## How every answer is verified

1. **Reading**: the question is turned into typed quantities (value, unit, kind, role) and one asked quantity; every stated number must be used or the question is refused.
2. **Planning**: only certified capabilities may answer; a capability applies only under the conditions the wording states (for example a pendulum, a liquid at rest, boiling).
3. **Computing**: exact symbolic or high-precision numeric evaluation.
4. **Checking**: dimensions, substitution back into the defining equation, and agreement between independent readings; any disagreement becomes a refusal, never a number.

Capabilities are certified before they ship: each must reproduce hand-worked reference cases, reject a deliberately wrong variant (a control), and leave every benchmark answer unchanged.

## Python worked examples, with their steps

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

![Michaelis–Menten kinetics: measured rates and fitted curve, computed by computational-biology-ai999x](https://huggingface.co/AI999X/computational-biology-ai999x/resolve/main/figures/bio-michaelis-menten-fit-fit.svg)

_Michaelis–Menten fit to 4 measured rates (S = [1, 2, 4, 8]): the package returns Vmax = 10 and Km = 1; the curve drawn is v = Vmax·S/(Km + S) with those two values (the double-reciprocal fit gives Vmax = 10, Km = 1, r² = 1). Marked: S = Km = 1, v = Vmax/2 = 5._

**Enzyme progress curve: substrate and product over time**

![Enzyme progress curve: substrate and product over time, computed by computational-biology-ai999x](https://huggingface.co/AI999X/computational-biology-ai999x/resolve/main/figures/bio-michaelis-menten-progress-kinetics.svg)

_Substrate and product over time for an enzyme reaction with S0 = 10, Vmax = 1 and Km = 2 (integrated Michaelis–Menten law), computed by the package at 101 times. Half the substrate is converted at t = 6.38629; marked: t = 6.38629, S = 5._

**Chemostat steady states against dilution rate**

![Chemostat steady states against dilution rate, computed by computational-biology-ai999x](https://huggingface.co/AI999X/computational-biology-ai999x/resolve/main/figures/bio-chemostat-steady-state-operating.svg)

_Chemostat operating diagram (Monod growth, μmax = 1, Ks = 1, Y = 0.5, feed S0 = 10): steady-state substrate and biomass computed by the package at 90 dilution rates up to washout; biomass productivity D·X peaks near D = 0.7 (2.68333). Marked: D = 0.5, X = 4.5, S = 1._

**Age-structured population growth (Leslie matrix)**

![Age-structured population growth (Leslie matrix), computed by computational-biology-ai999x](https://huggingface.co/AI999X/computational-biology-ai999x/resolve/main/figures/bio-leslie-projection-population.svg)

_Leslie-matrix projection of an age-structured population (fecundities F = [0, 1.5, 1], survival S = [0.5, 0.8], start n0 = [100, 50, 20]): each bar is the population the package returns at that step, split by age class. Growth settles to λ = 1.06152 per step (r = 0.0596993); at t = 12 the classes are [197.307, 92.9426, 69.9967]; marked: total = 360.247._

**Hardy–Weinberg: observed against expected genotypes**

![Hardy–Weinberg: observed against expected genotypes, computed by computational-biology-ai999x](https://huggingface.co/AI999X/computational-biology-ai999x/resolve/main/figures/bio-hardy-weinberg-test-bars.svg)

_Observed genotype counts AA = 298, Aa = 489, aa = 213 beside the counts expected under Hardy–Weinberg equilibrium, computed by the package from the allele frequencies p = 0.5425, q = 0.4575: 294.306, 496.387, 209.306. χ² = 0.22149, p-value 0.637907; observed heterozygosity 0.489 against expected 0.496387 (F = 0.0148825)._

**Heterozygosity lost by drift: population size against time**

![Heterozygosity lost by drift: population size against time, computed by computational-biology-ai999x](https://huggingface.co/AI999X/computational-biology-ai999x/resolve/main/figures/bio-drift-heterozygosity-heatmap.svg)

_Fraction of heterozygosity retained under random genetic drift, (1 − 1/(2N))^t, for population sizes N = 5, 10, 25, 50, 100, 250 (rows) after 10, 25, 50, 100, 200 generations (columns), each cell computed by the package and rounded to 3 decimals. N = 50 keeps 0.366 after 100 generations; N = 5 keeps 0 after 200._

## Capabilities

83 certified procedures. A few of them:

- **Procedures:** ACE species richness estimator; Allele frequency under haploid (genic) selection; Amino-acid composition of a protein; Aromaticity of a protein; Berger-Parker and McIntosh indices; Beta diversity between two samples; Brillouin diversity index; Chao1 species richness estimator; Chemostat steady state (Monod kinetics); Competitive inhibition global fit; Decline of a recessive lethal allele; Eadie-Hofstee fit of Michaelis-Menten data; Edit (Levenshtein) distance between sequences; Effective population size with unequal sex ratio; Effective size of a fluctuating population; Equilibrium under heterozygote advantage; Expected coalescence times; Faith's phylogenetic diversity; Fixation index FST from subpopulation allele frequencies; Fixation probability in the Moran process; Flux balance analysis (optimal objective); GC content of a nucleic acid; Global alignment score (Needleman-Wunsch, linear gaps); Global alignment score with affine gaps (Gotoh) …

`computational_biology_ai999x.templates()` lists every question form `ask()` reads; the capability guide (PDF, every example executed) is included in this repository.

## Out of scope

- Open-ended explanation, essay writing or tutoring dialogue: it is a solver, not a chatbot.
- Questions outside its certified computational biology capabilities, and any question whose stated numbers it cannot all use: it refuses with the reason instead.
- High-stakes engineering, medical or legal decisions without independent review.

## Limitations

- It answers what its certified capabilities cover (DNA and protein sequences, pairwise alignment, molecular masses, Mendelian and population genetics, phylogenetics, enzyme kinetics, population models and ecology); anything else is refused with a reason.
- `ask()` reads documented question forms and common phrasings; unusual wording may be refused; use `solve()` for guaranteed results.
- Results are as correct as the certified relations and their stated conditions; a certified relation is a checked model, not a claim about every physical situation.

## Citation

```bibtex
@software{computational_biology_ai999x_1_5_0,
  title  = {computational-biology-ai999x: an offline, verified computational biology solver},
  author = {Farabi, Ibnul Jaif}, organization = {AI999X},
  year   = {2026}, version = {1.5.0},
  url    = {https://huggingface.co/AI999X/computational-biology-ai999x}
}
```

## Contact

Questions, bug reports and requests: open a discussion on https://huggingface.co/AI999X/computational-biology-ai999x/discussions or visit https://ai999x.com.

## Support the project

computational-biology-ai999x is free and open source. If it helps you, you can support further development on [Buy Me a Coffee](https://buymeacoffee.com/ibnuljaiffarabi).

Created and engineered by **Ibnul Jaif Farabi** · **[AI999X](https://huggingface.co/AI999X)** · [ai999x.com](https://ai999x.com)

## Source and releases

This repository holds the installable source of each release, exactly as published on [PyPI](https://pypi.org/project/computational-biology-ai999x/) (`pip install .` builds it). Releases are tagged `v1.5.0`; the [GitHub release](https://github.com/AI999X/computational-biology-ai999x/releases/tag/v1.5.0) carries the wheel, the desktop app and the command-line program for Windows and Linux, and their checksums.

## Build the program from source

`rust/` holds the Rust source of the `computational-biology` command-line program and browser module: the AI999X engine (`rust/ax-engine`), this model (`rust/computational-biology`) and its data (`rust/models/computational_biology/pack.json`). With Rust 1.80 or later:

```bash
cd rust
cargo build --release -p computational-biology-ai999x
./target/release/computational-biology --help
cargo test -p computational-biology-ai999x
```
