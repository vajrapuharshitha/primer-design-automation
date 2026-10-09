# Automated PCR Primer Design & Screening

Python tool that generates candidate primer pairs with Primer3, then independently re-screens each pair and ranks the survivors.
Built for the "Primer Design" case-study theme (Tm calculation, GC filtering, hairpin/dimer screening).

**Template used in the demo:** NC_005816.1, *Yersinia pestis* plasmid pPCP1 (9,609 bp), a public sequence from the Biopython test suite. Any single-record FASTA works.

## Screens applied to every pair
- Tm 57-63 C per primer, Tm difference <= 2 C
- GC content 40-60%
- Hairpin, self-dimer, and cross-dimer free energy >= -6 kcal/mol
- Each primer has exactly one exact-match binding site on the template (both strands)

## Demo result (target region 1000-3500, product 200-500 bp)
Primer3 returned 50 candidate pairs; 19 passed every screen. All 31 rejected pairs failed the secondary-structure screen.
Top pair: TGGCGTTTACTACAGGCAGG / CAGACAAGCTGTGACCGTCT, 438 bp product, Tm 60.0 / 60.0 C, GC 55% / 55%, cross-dimer dG -3.9 kcal/mol.

Limitations: specificity is checked by exact match on the supplied template only, not by genome-wide BLAST; wet-lab validation has not been done.

## Run
```
pip install -r requirements.txt
python design_primers.py --fasta data/NC_005816.fna --start 1000 --end 3500 --out results/primer_pairs.csv
```
