"""Biomarker Risk Classifier: a small Streamlit app around the pipeline.
Use the built-in breast cancer dataset or upload your own CSV (one label column plus numeric biomarker columns)."""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import streamlit as st
from sklearn.datasets import load_breast_cancer
from sklearn.model_selection import train_test_split, StratifiedKFold, cross_val_score
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.feature_selection import SelectKBest, mutual_info_classif
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import roc_auc_score, roc_curve
from xgboost import XGBClassifier

SEED = 42
st.set_page_config(page_title="Biomarker Risk Classifier", layout="wide")
st.title("Biomarker Risk Classifier")
st.caption("Ranks biomarkers and classifies risk with leakage-safe cross-validation. Demonstration tool, not for clinical use.")

# ---------- data ----------
src = st.sidebar.radio("Data", ["Built-in: breast cancer (n=569)", "Upload CSV"])
if src.startswith("Built-in"):
    d = load_breast_cancer(as_frame=True)
    df = d.data.copy()
    df["malignant"] = 1 - d.target
    label, pos = "malignant", 1
else:
    f = st.sidebar.file_uploader("CSV file", type="csv")
    if f is None:
        st.info("Upload a CSV with one label column and numeric biomarker columns.")
        st.stop()
    df = pd.read_csv(f)
    label = st.sidebar.selectbox("Label column", df.columns)
    pos = st.sidebar.selectbox("Positive (disease) class", sorted(df[label].dropna().unique().tolist(), key=str))

feats = [c for c in df.select_dtypes("number").columns if c != label]
df = df[feats + [label]].dropna()
y = (df[label] == pos).astype(int)
X = df[feats]
if len(feats) < 2 or y.sum() < 10 or (len(y) - y.sum()) < 10:
    st.error("Need at least 2 numeric features and 10 samples in each class.")
    st.stop()
st.sidebar.write(f"{len(df)} samples, {len(feats)} features, {int(y.sum())} positive")

k = st.sidebar.slider("Top-k biomarkers to keep", 2, len(feats), min(10, len(feats)))
model_name = st.sidebar.selectbox("Model", ["Logistic Regression", "Random Forest", "XGBoost"])

def make(name):
    m = {"Logistic Regression": LogisticRegression(max_iter=2000, random_state=SEED),
         "Random Forest": RandomForestClassifier(n_estimators=300, random_state=SEED),
         "XGBoost": XGBClassifier(n_estimators=200, max_depth=3, learning_rate=0.1,
                                  eval_metric="logloss", random_state=SEED)}[name]
    return Pipeline([("scale", StandardScaler()), ("select", SelectKBest(mutual_info_classif, k=k)), ("clf", m)])

# ---------- run ----------
if st.sidebar.button("Run pipeline", type="primary"):
    Xtr, Xte, ytr, yte = train_test_split(X, y, test_size=0.25, stratify=y, random_state=SEED)
    pipe = make(model_name)
    cv = cross_val_score(pipe, Xtr, ytr, cv=StratifiedKFold(5, shuffle=True, random_state=SEED), scoring="roc_auc")
    pipe.fit(Xtr, ytr)
    ptr, pte = pipe.predict_proba(Xtr)[:, 1], pipe.predict_proba(Xte)[:, 1]
    _, tpr, thr = roc_curve(ytr, ptr)
    t0 = float(thr[np.where(tpr >= 0.95)[0][0]])
    clf = pipe.named_steps["clf"]
    imp = np.abs(clf.coef_[0]) if hasattr(clf, "coef_") else clf.feature_importances_
    rank = pd.DataFrame({"biomarker": X.columns[pipe.named_steps["select"].get_support()], "importance": imp}
                        ).sort_values("importance", ascending=False)
    st.session_state["res"] = dict(cv=cv, yte=yte.values, pte=pte, auc=roc_auc_score(yte, pte),
                                   t0=round(t0, 2), rank=rank, name=model_name, ntr=len(ytr), nte=len(yte))

res = st.session_state.get("res")
if not res:
    st.info("Choose settings in the sidebar and press **Run pipeline**.")
    st.stop()

c1, c2, c3 = st.columns(3)
c1.metric("Model", res["name"])
c2.metric("Cross-validated AUC", f"{res['cv'].mean():.3f} ± {res['cv'].std():.3f}", help=f"5-fold CV on {res['ntr']} training samples")
c3.metric("Hold-out AUC", f"{res['auc']:.3f}", help=f"{res['nte']} samples never used for training or selection")

thr = st.slider("Risk threshold (default = lowest threshold giving ≥95% sensitivity on training data)", 0.0, 1.0, float(min(max(res["t0"], 0), 1)), 0.01)
yte, pte = res["yte"], res["pte"]
flag = pte >= thr
tp, fn = int((flag & (yte == 1)).sum()), int((~flag & (yte == 1)).sum())
fp, tn = int((flag & (yte == 0)).sum()), int((~flag & (yte == 0)).sum())
m1, m2, m3 = st.columns(3)
m1.metric("Sensitivity", f"{tp / max(tp + fn, 1):.1%}")
m2.metric("Specificity", f"{tn / max(tn + fp, 1):.1%}")
m3.write(pd.DataFrame([[tp, fn], [fp, tn]], index=["Disease", "No disease"], columns=["Flagged", "Not flagged"]))

a, b = st.columns(2)
fpr, tpr, _ = roc_curve(yte, pte)
fig, ax = plt.subplots(figsize=(4.5, 4))
ax.plot(fpr, tpr, color="#1F3A5F"); ax.plot([0, 1], [0, 1], "k--", lw=0.8)
ax.scatter([fp / max(fp + tn, 1)], [tp / max(tp + fn, 1)], color="#e0a800", zorder=3, s=60)
ax.set(xlabel="False positive rate", ylabel="True positive rate", title="Hold-out ROC")
a.pyplot(fig)
fig2, ax2 = plt.subplots(figsize=(4.5, 4))
r = res["rank"][::-1]
ax2.barh(r["biomarker"], r["importance"], color="#1F3A5F"); ax2.set(title="Selected biomarkers", xlabel="Importance")
fig2.tight_layout(); b.pyplot(fig2)
st.dataframe(res["rank"].reset_index(drop=True), width="stretch")
