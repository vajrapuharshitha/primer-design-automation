"""
Automated primer design and screening.

Generates candidate PCR primer pairs with Primer3, then independently re-screens each
pair in Python: melting temperature (Tm), GC content, hairpin / self-dimer /
cross-dimer free energy, Tm mismatch between primers, and binding-site uniqueness
on the template. Passing pairs are ranked and written to CSV.

Example:
    python design_primers.py --fasta data/NC_005816.fna --start 1000 --end 3500 \
        --out results/primer_pairs.csv
"""
import argparse
import pandas as pd
import primer3
from Bio import SeqIO
from Bio.Seq import Seq

# Screening thresholds (common PCR design guidelines).
TM_MIN, TM_MAX = 57.0, 63.0
GC_MIN, GC_MAX = 40.0, 60.0
MAX_TM_DIFF = 2.0
MIN_STRUCTURE_DG = -6.0   # kcal/mol; more negative = stronger unwanted structure
MAX_BINDING_SITES = 1


def gc_percent(s):
    return 100.0 * (s.count("G") + s.count("C")) / len(s)


def dg_kcal(result):
    return result.dg / 1000.0  # primer3 returns cal/mol


def binding_sites(template, primer):
    """Exact-match binding sites on both strands of the template."""
    rc = str(Seq(primer).reverse_complement())
    fwd = template.count(primer)
    rev = template.count(rc)
    return fwd + rev


def screen_pair(template, left, right):
    row = {"left_primer": left, "right_primer": right}
    for tag, s in (("L", left), ("R", right)):
        row[f"{tag}_tm"] = round(primer3.calc_tm(s), 2)
        row[f"{tag}_gc"] = round(gc_percent(s), 1)
        row[f"{tag}_hairpin_dg"] = round(dg_kcal(primer3.calc_hairpin(s)), 2)
        row[f"{tag}_homodimer_dg"] = round(dg_kcal(primer3.calc_homodimer(s)), 2)
        row[f"{tag}_sites"] = binding_sites(template, s)
    row["cross_dimer_dg"] = round(dg_kcal(primer3.calc_heterodimer(left, right)), 2)
    row["tm_diff"] = round(abs(row["L_tm"] - row["R_tm"]), 2)
    return row


def passes(r):
    checks = {
        "tm": all(TM_MIN <= r[f"{t}_tm"] <= TM_MAX for t in "LR"),
        "gc": all(GC_MIN <= r[f"{t}_gc"] <= GC_MAX for t in "LR"),
        "tm_diff": r["tm_diff"] <= MAX_TM_DIFF,
        "structure": all(r[k] >= MIN_STRUCTURE_DG for k in
                         ("L_hairpin_dg", "R_hairpin_dg", "L_homodimer_dg",
                          "R_homodimer_dg", "cross_dimer_dg")),
        "unique": all(r[f"{t}_sites"] <= MAX_BINDING_SITES for t in "LR"),
    }
    return checks


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--start", type=int, required=True, help="target region start (0-based)")
    ap.add_argument("--end", type=int, required=True, help="target region end (exclusive)")
    ap.add_argument("--product-min", type=int, default=200)
    ap.add_argument("--product-max", type=int, default=500)
    ap.add_argument("--num-return", type=int, default=50)
    ap.add_argument("--out", default="results/primer_pairs.csv")
    a = ap.parse_args()

    record = next(SeqIO.parse(a.fasta, "fasta"))
    template = str(record.seq).upper()
    print(f"Template: {record.id[:60]}... ({len(template)} bp)")

    res = primer3.bindings.design_primers(
        {"SEQUENCE_ID": record.id[:30], "SEQUENCE_TEMPLATE": template,
         "SEQUENCE_INCLUDED_REGION": [a.start, a.end - a.start]},
        {"PRIMER_OPT_SIZE": 20, "PRIMER_MIN_SIZE": 18, "PRIMER_MAX_SIZE": 25,
         "PRIMER_OPT_TM": 60.0, "PRIMER_MIN_TM": 55.0, "PRIMER_MAX_TM": 65.0,
         "PRIMER_MIN_GC": 35.0, "PRIMER_MAX_GC": 65.0,
         "PRIMER_PRODUCT_SIZE_RANGE": [[a.product_min, a.product_max]],
         "PRIMER_NUM_RETURN": a.num_return})

    n = res["PRIMER_PAIR_NUM_RETURNED"]
    print(f"Primer3 returned {n} candidate pairs")
    rows = []
    for i in range(n):
        left = res[f"PRIMER_LEFT_{i}_SEQUENCE"]
        right = res[f"PRIMER_RIGHT_{i}_SEQUENCE"]
        r = screen_pair(template, left, right)
        r["product_bp"] = res[f"PRIMER_PAIR_{i}_PRODUCT_SIZE"]
        r["left_start"] = res[f"PRIMER_LEFT_{i}"][0]
        r["primer3_penalty"] = round(res[f"PRIMER_PAIR_{i}_PENALTY"], 3)
        chk = passes(r)
        r["passes_all"] = all(chk.values())
        r["failed_checks"] = ",".join(k for k, v in chk.items() if not v)
        rows.append(r)

    df = pd.DataFrame(rows).sort_values(["passes_all", "primer3_penalty"], ascending=[False, True])
    df.to_csv(a.out, index=False)
    ok = int(df["passes_all"].sum())
    print(f"{ok} of {n} pairs passed all screens -> {a.out}")
    cols = ["left_primer", "right_primer", "product_bp", "L_tm", "R_tm", "L_gc", "R_gc", "cross_dimer_dg"]
    print(df[df["passes_all"]].head(5)[cols].to_string(index=False))


if __name__ == "__main__":
    main()
