"""Primer Design & Screening: a small Streamlit app around design_primers.py.
Design candidate primer pairs with Primer3, then re-screen them for Tm, GC, hairpin/dimer structure and binding uniqueness."""
from io import StringIO
import pandas as pd
import primer3
import streamlit as st
from Bio import SeqIO
from design_primers import screen_pair, passes

st.set_page_config(page_title="Primer Design & Screening", layout="wide")
st.title("Primer Design & Screening")
st.caption("Primer3 candidates, independently re-screened in Python. Demonstration tool: always confirm primers in Primer-BLAST and in the lab.")

src = st.sidebar.radio("Template sequence", ["Demo: Yersinia pestis plasmid pPCP1 (9.6 kb)", "Paste or upload FASTA"])
if src.startswith("Demo"):
    record = next(SeqIO.parse("data/NC_005816.fna", "fasta"))
else:
    up = st.sidebar.file_uploader("FASTA file", type=["fa", "fasta", "fna", "txt"])
    txt = up.getvalue().decode() if up else st.sidebar.text_area("...or paste FASTA", height=150)
    try:
        record = next(SeqIO.parse(StringIO(txt), "fasta"))
    except StopIteration:
        st.info("Provide a single-record FASTA sequence in the sidebar.")
        st.stop()
template = str(record.seq).upper()
if set(template) - set("ACGTN") or len(template) < 300:
    st.error("Sequence must be DNA (A, C, G, T, N) and at least 300 bp.")
    st.stop()
st.sidebar.write(f"{record.id[:40]} | {len(template)} bp")

L = len(template)
default_start = min(1000, L // 4)
default_end = min(3500, L)
start, end = st.sidebar.slider(
    "Target region (bp)", 0, L, (default_start, default_end)
)
pmin, pmax = st.sidebar.slider(
    "Product size (bp)", 100, 1000, (200, 500), step=10
)
n = st.sidebar.slider("Candidate pairs to generate", 10, 100, 50, step=10)

if st.sidebar.button("Design primers", type="primary"):
    if end - start < pmin:
        st.error("Target region is shorter than the minimum product size.")
        st.stop()
    res = primer3.bindings.design_primers(
        {"SEQUENCE_ID": record.id[:30], "SEQUENCE_TEMPLATE": template, "SEQUENCE_INCLUDED_REGION": [start, end - start]},
        {"PRIMER_OPT_SIZE": 20, "PRIMER_MIN_SIZE": 18, "PRIMER_MAX_SIZE": 25, "PRIMER_OPT_TM": 60.0,
         "PRIMER_MIN_TM": 55.0, "PRIMER_MAX_TM": 65.0, "PRIMER_MIN_GC": 35.0, "PRIMER_MAX_GC": 65.0,
         "PRIMER_PRODUCT_SIZE_RANGE": [[pmin, pmax]], "PRIMER_NUM_RETURN": n})
    rows = []
    for i in range(res["PRIMER_PAIR_NUM_RETURNED"]):
        r = screen_pair(template, res[f"PRIMER_LEFT_{i}_SEQUENCE"], res[f"PRIMER_RIGHT_{i}_SEQUENCE"])
        r["product_bp"] = res[f"PRIMER_PAIR_{i}_PRODUCT_SIZE"]
        r["primer3_penalty"] = round(res[f"PRIMER_PAIR_{i}_PENALTY"], 3)
        chk = passes(r)
        r["passes_all"] = all(chk.values())
        r["failed_checks"] = ",".join(k for k, v in chk.items() if not v)
        rows.append(r)
    st.session_state["df"] = (pd.DataFrame(rows).sort_values(["passes_all", "primer3_penalty"], ascending=[False, True])
                              if rows else pd.DataFrame())

df = st.session_state.get("df")
if df is None:
    st.info("Choose a region in the sidebar and press **Design primers**.")
    st.stop()
if df.empty:
    st.warning("Primer3 returned no pairs. Try a wider region or a larger product-size range.")
    st.stop()

ok = df[df["passes_all"]]
c1, c2, c3 = st.columns(3)
c1.metric("Candidate pairs", len(df))
c2.metric("Passed all screens", len(ok))
c3.metric("Screens", "Tm, GC, structure, uniqueness")
cols = ["left_primer", "right_primer", "product_bp", "L_tm", "R_tm", "L_gc", "R_gc", "cross_dimer_dg", "primer3_penalty"]
st.subheader("Recommended pairs")
st.dataframe(ok[cols].reset_index(drop=True), width="stretch")
with st.expander("All candidates, including rejected pairs and why"):
    st.dataframe(df[cols + ["failed_checks"]].reset_index(drop=True), width="stretch")
st.download_button("Download all results (CSV)", df.to_csv(index=False), "primer_pairs.csv", "text/csv")
st.markdown("Next step: check specificity of your chosen pair with [NCBI Primer-BLAST](https://www.ncbi.nlm.nih.gov/tools/primer-blast/).")
