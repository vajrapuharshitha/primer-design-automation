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
start, end = st.sidebar.slider("Target region (bp)", 0, L, (min(1000, L //
