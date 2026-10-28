"""computational-biology-ai999x — offline computational biology answers with their formula, steps and verification.

    import computational_biology_ai999x as m
    m.ask("Reverse complement: seq = ATGCCGTA").show()                               # TACGGCAT
    m.ask("Hardy-Weinberg test: n_AA = 298, n_Aa = 489, n_aa = 213").show()          # p = 0.5425, chi2 = 0.2215
    m.ask("Cross phenotype ratio: parent1 = AaBb, parent2 = AaBb").show()            # 9:3:3:1
    m.solve("bio_michaelis_menten_fit", S="0.5, 1, 2, 5", v="2.2, 3.5, 5.1, 6.6")      # Vmax = 8.5001, Km = 1.4000
    m.capabilities()                   # every certified capability

83 certified capabilities: DNA, RNA and protein sequences (reverse complement, transcription, translation with the
standard genetic code, GC content, k-mer counts, p-distance and edit distance, longest open reading frame, restriction
sites and fragments, nucleic-acid and protein masses, extinction coefficient, net charge, isoelectric point, GRAVY,
instability index, aromaticity, amino-acid composition, secondary-structure propensities, Wallace and nearest-neighbour
melting temperatures), pairwise alignment scores (Needleman-Wunsch, Smith-Waterman and semi-global; linear and affine
gaps; BLOSUM62 protein alignments), Mendelian and population genetics (genotype probabilities and phenotype ratios of
crosses, Hardy-Weinberg chi-square and exact tests, haploid selection, recessive lethals, heterozygote advantage,
mutation-selection balance, drift and heterozygosity loss, Moran fixation probability, effective population size, FST,
linkage disequilibrium D, D' and r^2, Haldane and Kosambi map functions, pedigree inbreeding coefficients),
phylogenetics (Jukes-Cantor and Kimura distances, UPGMA, neighbour joining, tree counts, Robinson-Foulds distance, Fitch
parsimony, coalescent times, nucleotide diversity, Watterson's theta, Tajima's D), enzyme kinetics (Michaelis-Menten
nonlinear fit, Lineweaver-Burk, Hanes-Woolf and Eadie-Hofstee fits, progress curves and conversion times, Hill plot and
Hill fit, competitive-inhibition fit) and population biology and ecology (Leslie projection and growth rate, Shannon,
Simpson, Brillouin, Margalef, Menhinick, Berger-Parker and McIntosh indices, Chao1, ACE, Good's coverage, rarefaction,
Bray-Curtis, Jaccard and Sorensen dissimilarities, species-area fits, Petersen and Schumacher-Eschmeyer mark-recapture
estimates, chemostat steady state, flux balance analysis, Faith's phylogenetic diversity).

Questions use the form "<capability phrase>: name = value, ..." (see templates()): sequences, genotypes and Newick trees
are written as they are, lists in brackets, matrices as [[..], [..]], and conventions (log base, Simpson estimator, Chao1
form, map function, strand, stop handling, enzyme, molecule type) are stated, never assumed. Built by the AI999X model
machine from the same engine as every domain; a question that is not determined is refused rather than guessed.

Creator / engineer: Ibnul Jaif Farabi. Company: AI999X.
"""
from __future__ import annotations

from pathlib import Path

try:                                       # the installed package carries its own copy of the engine
    from ._engine import Answer
    from ._engine.product import Product
except ImportError:                        # the factory's source tree, during development
    from ai999x_factory.runtime import Answer
    from ai999x_factory.runtime.product import Product

__version__ = "1.0.1"
_P = Product("computational_biology", [Path(__file__).resolve().parent / "pack", Path(__file__).resolve().parents[3]])
pack, ask, solve, compile_problem = _P.pack, _P.ask, _P.solve, _P.compile_problem
formula, formulas, quantities, conditions, capabilities = _P.formula, _P.formulas, _P.quantities, _P.conditions, _P.capabilities
examples, family_examples, structured_example, to_si = _P.examples, _P.family_examples, _P.structured_example, _P.to_si
templates, procedure_examples = _P.templates, _P.procedure_examples
prove, derivations, about = _P.prove, _P.derivations, (lambda: {**_P.about(), "version": __version__})
knowledge, lookup, related = _P.knowledge, _P.lookup, _P.related
plot, plot_relation, search = _P.plot, _P.plot_relation, _P.search   # figures (SVG) and capability search ([search])
__all__ = ["ask", "solve", "prove", "derivations", "about", "formula", "formulas", "quantities", "conditions", "capabilities",
           "compile_problem", "pack", "Answer", "plot", "plot_relation", "search", "__version__"]
