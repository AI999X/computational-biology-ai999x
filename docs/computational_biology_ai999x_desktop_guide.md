# Computational Biology: the desktop app, user guide

**Computational Biology**, desktop app 1.6.5 for Windows and Linux. **Creator / engineer:** Ibnul Jaif Farabi · **Company:** AI999X · support@ai999x.com

Computational Biology covers DNA and protein sequences, alignment, molecular weights, genetics and population genetics, phylogenetics, enzyme kinetics, metabolic and population models, and ecology. Its desktop app holds all 86 calculations in 10 topics: pick one, change its inputs and the answer, a live scene and charts follow at once. It runs on your own computer with no network connection and nothing else to install. The pictures in this guide were taken in this app.

## Contents

- 1. Get it and start it
- 2. The window
- 3. Compute a calculation
- 4. Ask in words
- 5. The live scene and the charts
- 6. Light and dark
- 7. Every calculation
- 8. Other ways to use it

## 1. Get it and start it

| system | download |
|---|---|
| Windows x86-64 | https://github.com/AI999X/ai999x-desktop/releases/download/v1.6.5/ai999x-computational-biology-1.6.5-windows-x86_64.zip |
| Linux x86-64 | https://github.com/AI999X/ai999x-desktop/releases/download/v1.6.5/ai999x-computational-biology-1.6.5-linux-x86_64.tar.gz |

Windows: Download ai999x-computational-biology-1.6.5-windows-x86_64.zip, right-click it → Extract All, then open the folder and start ai999x-computational-biology.exe. If Windows SmartScreen asks, choose More info → Run anyway.

Linux:

```shell
curl -LO https://github.com/AI999X/ai999x-desktop/releases/download/v1.6.5/ai999x-computational-biology-1.6.5-linux-x86_64.tar.gz
curl -LO https://github.com/AI999X/ai999x-desktop/releases/download/v1.6.5/SHA256SUMS
sha256sum -c --ignore-missing SHA256SUMS
tar xzf ai999x-computational-biology-1.6.5-linux-x86_64.tar.gz
./ai999x-computational-biology-1.6.5-linux-x86_64/ai999x-computational-biology
```

Checksums of every download: https://github.com/AI999X/ai999x-desktop/releases/download/v1.6.5/SHA256SUMS. Install steps for every app: https://www.ai999x.com/desktop-cli/

## 2. The window

![The Computational Biology app, opened on Reverse complement of a DNA sequence with its worked example.](https://www.ai999x.com/?ax_shot=computational-biology-ai999x/hero.webp)

*The Computational Biology app, opened on Reverse complement of a DNA sequence with its worked example.*

On the left, a search box and the topics with their calculations. In the middle, the chosen calculation: its answer at the top, then the live scene (or its charts), what the answer means and the formula. On the right, its inputs, each with a slider where that makes sense, and the Compute and Reset buttons. At the bottom, a box for a question in words.

## 3. Compute a calculation

Open Sequences › Reverse complement of a DNA sequence (or type part of its name in the search box). It opens with a worked example filled in: “Reverse complement: seq = ATGCCGTA” The answer: revcomp = TACGGCAT.

Change any input (type a number with its unit, or drag its slider) and the answer, the scene and the charts follow; Compute recomputes after typing, Reset brings back the worked example. An input outside the range the formula holds for is refused with the reason, never answered with a guessed number.

![Reverse complement of a DNA sequence: the live scene of the worked example.](https://www.ai999x.com/?ax_shot=computational-biology-ai999x/scene-1.webp)

*Reverse complement of a DNA sequence: the live scene of the worked example.*

## 4. Ask in words

Type a question as you would say it, for example “Reverse complement: seq = ATGCCGTA”. The app reads the quantities and their units, picks the calculation and answers; the command-line program reads it the same way and printed:

```text
Read as: Reverse complement of a DNA sequence
revcomp = TACGGCAT
  by: Reverse complement of a DNA sequence  (complement each base (A-T, C-G) and reverse the order)
```

When a question names a quantity the app cannot place, it says what is missing instead of guessing.

## 5. The live scene and the charts

Every calculation draws a scene from its answer: it moves as an input changes. The Charts button shows charts that vary one input around your values; Enlarge opens one in a large view. Scenes from other topics:

![Hardy-Weinberg chi-square test](https://www.ai999x.com/?ax_shot=computational-biology-ai999x/scene-2.webp)

*Hardy-Weinberg chi-square test*

![Global alignment score (Needleman-Wunsch, linear gaps)](https://www.ai999x.com/?ax_shot=computational-biology-ai999x/scene-3.webp)

*Global alignment score (Needleman-Wunsch, linear gaps)*

![Phenotype ratio of a cross (complete dominance)](https://www.ai999x.com/?ax_shot=computational-biology-ai999x/scene-4.webp)

*Phenotype ratio of a cross (complete dominance)*

![Michaelis-Menten nonlinear least-squares fit](https://www.ai999x.com/?ax_shot=computational-biology-ai999x/scene-5.webp)

*Michaelis-Menten nonlinear least-squares fit*

![ACE species richness estimator](https://www.ai999x.com/?ax_shot=computational-biology-ai999x/scene-6.webp)

*ACE species richness estimator*

![Chemostat steady state (Monod kinetics)](https://www.ai999x.com/?ax_shot=computational-biology-ai999x/scene-7.webp)

*Chemostat steady state (Monod kinetics)*

![Parsimony score of a tree (Fitch)](https://www.ai999x.com/?ax_shot=computational-biology-ai999x/scene-8.webp)

*Parsimony score of a tree (Fitch)*

![Reverse complement of a DNA sequence: a chart, enlarged.](https://www.ai999x.com/?ax_shot=computational-biology-ai999x/chart-1.webp)

*Reverse complement of a DNA sequence: a chart, enlarged.*

![Hardy-Weinberg chi-square test: a chart, enlarged.](https://www.ai999x.com/?ax_shot=computational-biology-ai999x/chart-2.webp)

*Hardy-Weinberg chi-square test: a chart, enlarged.*

## 6. Light and dark

The app follows a light or a dark theme (the switch is in its top bar).

![The Computational Biology app in the dark theme.](https://www.ai999x.com/?ax_shot=computational-biology-ai999x/hero-dark.webp)

*The Computational Biology app in the dark theme.*

## 7. Every calculation

86 calculations in 10 topics; each has its own entry on the model's page (https://www.ai999x.com/models/computational-biology-ai999x/) with its formula, worked example and charts.

| topic | calculations |
|---|---|
| Alignment | Global alignment score (Needleman-Wunsch, linear gaps), Global alignment score with affine gaps (Gotoh), Global protein alignment score (BLOSUM62, affine gaps), Local alignment score (Smith-Waterman, linear gaps), Local alignment score with affine gaps, Local protein alignment score (BLOSUM62, affine gaps), Semi-global (end-gap-free) alignment score |
| Ecology | ACE species richness estimator, Berger-Parker and McIntosh indices, Beta diversity between two samples, Brillouin diversity index, Chao1 species richness estimator, Faith's phylogenetic diversity, Good's coverage of a sample, Margalef and Menhinick richness indices, Mark-recapture population estimate (Petersen MLE), Rarefied species richness, Schumacher-Eschmeyer multiple-census estimate, Shannon diversity index, Simpson diversity indices, Species-area relationship fit |
| Enzyme kinetics | Competitive inhibition global fit, Eadie-Hofstee fit of Michaelis-Menten data, Hanes-Woolf fit of Michaelis-Menten data, Hill coefficient from a Hill plot, Hill equation nonlinear least-squares fit, Lineweaver-Burk fit of Michaelis-Menten data, Michaelis-Menten nonlinear least-squares fit, Substrate remaining in a Michaelis-Menten reaction, Time to reach a conversion in a Michaelis-Menten reaction |
| Genetics | Inbreeding coefficient from a pedigree, Map distance from a recombination fraction, Phenotype ratio of a cross (complete dominance), Probability of an offspring genotype in a cross, Recombination fraction from a map distance, Recombination frequency from two-point testcross counts, Three-point testcross: gene order, map distances and interference |
| Microbial growth | Chemostat steady state (Monod kinetics) |
| Phylogenetics | Jukes-Cantor distance, Kimura two-parameter distance, Neighbour-joining branch lengths, Number of possible phylogenetic trees, Parsimony score of a tree (Fitch), Robinson-Foulds distance between trees, UPGMA tree heights |
| Population biology | Leslie matrix growth rate and stable age distribution, Leslie matrix population projection |
| Population genetics | Allele frequency under haploid (genic) selection, Decline of a recessive lethal allele, Effective population size with unequal sex ratio, Effective size of a fluctuating population, Equilibrium under heterozygote advantage, Expected coalescence times, Fixation index FST from subpopulation allele frequencies, Fixation probability in the Moran process, Hardy-Weinberg chi-square test, Hardy-Weinberg exact test, Hardy-Weinberg genotype and carrier frequencies, Linkage disequilibrium D, D' and r^2, Loss of heterozygosity by genetic drift, Mutation-selection balance, Nucleotide diversity, Tajima's D, Watterson's estimator of theta |
| Sequences | Amino-acid composition of a protein, Aromaticity of a protein, Edit (Levenshtein) distance between sequences, GC content of a nucleic acid, Grand average of hydropathy (GRAVY), Isoelectric point of a protein, k-mer counts of a sequence, Longest open reading frame, Melting temperature by nearest-neighbour thermodynamics, Melting temperature by the Wallace rule, Molecular mass of a protein, Molecular mass of DNA or RNA, Net charge of a protein at a pH, p-distance (Hamming) between aligned sequences, Protein extinction coefficient at 280 nm, Protein instability index, Restriction sites and fragments, Reverse complement of a DNA sequence, Secondary-structure propensity fractions, Transcription of DNA into mRNA, Translation with the standard genetic code |
| Systems biology | Flux balance analysis (optimal objective) |

## 8. Other ways to use it

The same engine runs on the command line (`computational-biology`, its own guide), in Python (`pip install computational-biology-ai999x`, `import computational_biology_ai999x`, its own guide) and in the browser on https://www.ai999x.com/models/computational-biology-ai999x/.

Licence: Apache-2.0 · AI999X · https://www.ai999x.com/models/computational-biology-ai999x/
