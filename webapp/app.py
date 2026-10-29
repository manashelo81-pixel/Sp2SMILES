"""Spectrum -> Structure web app (Streamlit).

Run with:   python run_webapp.py          (from the smiles_ms folder)
or:         streamlit run webapp/app.py
"""
import io
import json
import os
import shutil
import sys
import tempfile

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import pandas as pd  # noqa: E402
import streamlit as st  # noqa: E402
import streamlit.components.v1 as components  # noqa: E402
from rdkit import Chem  # noqa: E402
from rdkit.Chem import Draw  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import fetch  # noqa: E402
import pipeline as P  # noqa: E402
try:
    import glossary as G  # noqa: E402
except ImportError:                                  # glossary.py not copied: the app still works, without tooltips
    class G:  # noqa: N801
        GLOSSARY, LIBRARY_HELP, RETRIEVAL_HELP = {}, [], []


def gl(key):
    return G.GLOSSARY.get(key, "")

ROOT = os.path.dirname(HERE)
DEFAULT_MODEL = os.environ.get("SMILES_MS_MODEL", os.path.join(ROOT, "models", "model.pt"))
DEFAULT_DATASET = os.environ.get("SMILES_MS_DATASET", os.path.join(ROOT, "models", "dataset.pkl"))

st.set_page_config(page_title="Spectrum → Structure", page_icon="🧪", layout="wide")


# ----------------------------------------------------------------------------- helpers
@st.cache_resource(show_spinner=False)
def get_resources(model_path, dataset_path):
    return P.Resources(model_path, dataset_path)


def show_df(df, **kw):
    """Table whose column headers carry a (?) tooltip with a plain-language explanation."""
    cfg = {c: st.column_config.Column(c, help=gl(c)) for c in df.columns if gl(c)}
    try:
        st.dataframe(df, width="stretch", hide_index=True, column_config=cfg, **kw)
    except TypeError:                      # older Streamlit versions
        st.dataframe(df, use_container_width=True, hide_index=True, column_config=cfg, **kw)


def info(key, label="ℹ️"):
    """An info button that opens a small popover explaining one term."""
    with st.popover(label, help=f"What is '{key}'?"):
        st.markdown(f"**{key}**")
        st.write(gl(key))


def info_many(keys, title="ℹ️ What do these numbers mean?"):
    with st.popover(title):
        for k in keys:
            st.markdown(f"**{k}** - {gl(k)}")


def param(label, value, key=None):
    """One parameter line with its own info button."""
    c1, c2 = st.columns([14, 1])
    c1.markdown(f"**{label}:** {value}")
    if key:
        with c2:
            info(key, "ℹ️")


def mol_image(smiles, size=(260, 200)):
    mol = Chem.MolFromSmiles(smiles)
    return Draw.MolToImage(mol, size=size) if mol is not None else None


def plot_spectrum(f, hit=None):
    fig, ax = plt.subplots(figsize=(8, 3.6))
    p = f["peaks"]
    ax.vlines(p[:, 0], 0, p[:, 1] * 100, color="#1f77b4", lw=1.6, label="measured")
    if hit is not None:
        lp = hit["peaks"]
        ax.vlines(lp[:, 0], 0, -lp[:, 1] * 100, color="#d9822b", lw=1.6, label="library")
    ax.axvline(f["precursor_mz"], ls="--", color="grey", lw=0.8)
    ax.axhline(0, color="black", lw=0.5)
    ax.set_xlabel("m/z")
    ax.set_ylabel("relative intensity (%)")
    ax.legend(loc="upper left", frameon=False)
    fig.tight_layout()
    return fig


def cards(hits, kind, per_row=5):
    for start in range(0, len(hits), per_row):
        cols = st.columns(per_row)
        for col, h in zip(cols, hits[start:start + per_row]):
            with col:
                img = mol_image(h["smiles"])
                if img is not None:
                    st.image(img, width=230)
                if kind == "library":
                    st.caption(f"**#{h['rank']}** cosine {h['score']:.3f} · {h['matched']} peaks · {h['confidence']}  \n"
                               f"{h['formula']} · {h['ppm']:+.1f} ppm  \n`{h['library_id']}`")
                else:
                    st.caption(f"**#{h['rank']}** fingerprint score {h['score']:.3f}  \n"
                               f"{h['formula']} · {h['ppm']:+.1f} ppm")
                st.code(h["smiles"], language=None)


def background(mode):
    """Animated page background (purely cosmetic: any failure is ignored)."""
    try:
        js = open(os.path.join(HERE, "background.js"), encoding="utf-8").read().replace("__MODE__", json.dumps(mode))
        html = f"<script>{js}</script>"
        if hasattr(st, "iframe"):                 # newer Streamlit
            st.iframe(html, height=1)
        else:                                     # older Streamlit (for example 1.51)
            components.html(html, height=0)
    except Exception:
        pass


def to_excel(sheets):
    buf = io.BytesIO()
    with pd.ExcelWriter(buf, engine="openpyxl") as xw:
        for name, df in sheets.items():
            df.to_excel(xw, sheet_name=name[:30], index=False)
    return buf.getvalue()


# ----------------------------------------------------------------------------- sidebar
st.sidebar.title("Settings")
# Model and library locations are set when the app is started (--model / --dataset, or the environment
# variables SMILES_MS_MODEL / SMILES_MS_DATASET), so they are not shown here.


@st.cache_resource(show_spinner="Downloading the model and spectral library (first start only) ...")
def prepare_files(model_path, dataset_path, model_url, dataset_url):
    m, e1 = fetch.ensure(model_path, model_url)
    d, e2 = fetch.ensure(dataset_path, dataset_url)
    return m, d, [e for e in (e1, e2) if e]


# On hosted deployments the files can be fetched from URLs (set SMILES_MS_MODEL_URL / SMILES_MS_DATASET_URL)
_model_path, _dataset_path, _download_errors = prepare_files(
    DEFAULT_MODEL, DEFAULT_DATASET, os.environ.get("SMILES_MS_MODEL_URL", ""), os.environ.get("SMILES_MS_DATASET_URL", ""))
res = get_resources(_model_path, _dataset_path)

BG_CHOICES = {"Electron spirals": "electrons", "Spectrum peaks": "spectrum", "Off": "off"}
_bg_default = {"electrons": 0, "spectrum": 1, "off": 2}.get(os.environ.get("SMILES_MS_BG", "electrons").lower(), 0)
bg_label = st.sidebar.selectbox("Background animation", list(BG_CHOICES), index=_bg_default,
                                help="Electron spirals follow your mouse. Spectrum peaks: move the mouse to scan m/z "
                                     "and click to add ions. Choose Off for a plain background.")

with st.sidebar.expander("Background and grouping"):
    blank_ppm = st.number_input("Blank match tolerance (ppm)", 1.0, 50.0, 10.0)
    do_cluster = st.checkbox("Group scans with the same precursor", True,
                             help="Merges repeated scans of one ion into a single feature. Switch off for peak lists that "
                                  "already contain one spectrum per feature.")
    cluster_ppm = st.number_input("Grouping tolerance (ppm)", 1.0, 50.0, 10.0)
    min_peaks = st.number_input("Minimum peaks per spectrum", 1, 50, 5)
with st.sidebar.expander("Adducts"):
    adduct_mode = st.selectbox("Adduct assumption", ["auto"] + P.C.ADDUCTS_BY_MODE["Positive"] + P.C.ADDUCTS_BY_MODE["Negative"],
                               help="auto: detect co-eluting adduct pairs (e.g. [M+H]+ / [M+NH4]+), otherwise [M+H]+ or [M-H]-")
    adduct_ppm = st.number_input("Adduct-pair tolerance (ppm)", 1.0, 20.0, 5.0)
    rt_tol = st.number_input("Adduct-pair RT tolerance (s)", 1.0, 60.0, 10.0)
with st.sidebar.expander("Library search"):
    do_library = st.checkbox("Enable library search", True, disabled=not res.has_dataset)
    lib_ppm = st.number_input("Precursor tolerance (ppm)", 1.0, 100.0, 20.0)
    frag_tol = st.number_input("Fragment tolerance (Da)", 0.001, 0.5, 0.02, format="%.3f")
    min_matched = st.number_input("Minimum matched peaks", 1, 30, 3)
    min_score = st.slider("Minimum cosine", 0.0, 1.0, 0.4, 0.05)
    top_lib = st.number_input("Hits per feature", 1, 20, 5)
with st.sidebar.expander("Fingerprint retrieval"):
    do_retrieval = st.checkbox("Enable retrieval", True, disabled=not res.has_model)
    ret_ppm = st.number_input("Candidate mass tolerance (ppm)", 1.0, 100.0, 10.0)
    top_ret = st.number_input("Candidates per feature", 1, 30, 10)
    instrument = st.selectbox("Instrument type (if known)", ["UNKNOWN", "ORBITRAP", "QTOF", "QQQ", "ITRAP", "OTHER"])
    ce = st.number_input("Collision energy (0 = unknown)", 0.0, 200.0, 0.0)
formula_ppm = st.sidebar.number_input("Formula tolerance (ppm)", 0.5, 20.0, 5.0)

# ----------------------------------------------------------------------------- main page
background(BG_CHOICES[bg_label])
st.title("🧪 Spectrum → Structure")
st.write("Upload LC-MS/MS data (mzML, mzXML, mzData, or MGF/MSP peak lists). The app groups repeated scans, removes background ions seen in your blanks, "
         "infers adducts, searches a spectral library and ranks candidate structures with a trained model.")

for _err in _download_errors:
    st.warning(_err)
missing = []
if not res.has_dataset:
    missing.append("a spectral library (library search is unavailable)")
if not res.has_model:
    missing.append("a trained model (candidate ranking is unavailable)")
if missing:
    st.info("Not found: " + " and ".join(missing) + ". Start the app with `--dataset` / `--model` to enable "
            "them (see README_WEB.md). You can still upload files to get grouped precursors, background removal, "
            "adduct inference and formulas.")

c1, c2, c3 = st.columns(3)
with c1:
    sample_files = st.file_uploader("Sample file(s): mzML, mzXML, mzData, MGF, MSP", type=["mzML", "mzml", "mzXML", "mzxml", "mzData", "mzdata", "mgf", "msp", "gz"], accept_multiple_files=True)
with c2:
    blank_files = st.file_uploader("Blank file(s) (optional)", type=["mzML", "mzml", "mzXML", "mzxml", "mzData", "mzdata", "mgf", "msp", "gz"], accept_multiple_files=True)
with c3:
    cand_file = st.file_uploader("Your own candidate structures (optional)", type=["txt", "csv", "smi", "tsv"],
                                 help="One SMILES per line. They are ranked in addition to the library molecules.")

if st.button("Run analysis", type="primary", disabled=not sample_files):
    tmp = tempfile.mkdtemp(prefix="smiles_ms_")
    try:
        def save(files):
            out = []
            for uf in files:
                path = os.path.join(tmp, f"{len(out)}_{os.path.basename(uf.name)}")
                with open(path, "wb") as fh:
                    fh.write(uf.getbuffer())
                out.append((uf.name, path))
            return out
        samples, blanks = save(sample_files), save(blank_files or [])
        user_pool, bad = None, 0
        if cand_file is not None:
            user_pool, bad = P.read_candidate_list(cand_file.getvalue().decode("utf-8", errors="replace"))
        params = dict(min_peaks=min_peaks, blank_ppm=blank_ppm, cluster_ppm=cluster_ppm, do_cluster=do_cluster, adduct_mode=adduct_mode,
                      adduct_ppm=adduct_ppm, rt_tol=rt_tol, do_library=do_library, lib_ppm=lib_ppm, frag_tol=frag_tol,
                      min_matched=int(min_matched), min_score=min_score, top_lib=int(top_lib),
                      do_retrieval=do_retrieval, ret_ppm=ret_ppm, top_ret=int(top_ret), formula_ppm=formula_ppm,
                      instrument=instrument, collision_energy=ce or None)
        bar = st.progress(0.0, text="starting")
        try:
            result = P.run_analysis(samples, blanks, params, res, user_pool,
                                    progress=lambda fr, msg: bar.progress(min(max(fr, 0.0), 1.0), text=msg))
            if bad:
                result["warnings"].append(f"{bad} line(s) of your candidate list were not valid SMILES and were skipped.")
            st.session_state["result"] = result
        except Exception as exc:                       # show a readable message instead of a stack trace
            st.session_state.pop("result", None)
            st.error(f"Analysis failed: {exc}")
        bar.empty()
    finally:
        shutil.rmtree(tmp, ignore_errors=True)

result = st.session_state.get("result")
if result:
    feats, summary = result["features"], result["summary"]
    for w in result["warnings"]:
        st.warning(w)
    if not feats:
        st.error("No precursors could be extracted from the uploaded file(s), so there is nothing to analyse. "
                 "See the message above for the reason, or check the terminal for the line "
                 "'<file>: N MS2 spectra found'.")
        st.stop()
    n_bg = sum(f["background"] for f in feats)
    n_lib = sum(1 for f in feats if f["lib_hits"])
    n_ret = sum(1 for f in feats if not f["lib_hits"] and (f["ret_hits"] or f["user_hits"]))
    m = st.columns(5)
    m[0].metric("MS/MS scans", result["n_scans"], help=gl("MS/MS scans"))
    m[1].metric("Unique precursors", len(feats), help=gl("Unique precursors"))
    m[2].metric("Background", n_bg, help=gl("Background"))
    m[3].metric("Library matches", n_lib, help=gl("Library matches"))
    m[4].metric("Retrieval only", n_ret, help=gl("Retrieval only"))

    tab_feat, tab_det, tab_dl, tab_help = st.tabs(["Features", "Feature detail", "Download", "How to read this"])

    with tab_feat:
        statuses = sorted(summary["Status"].unique())
        pick = st.multiselect("Show", statuses, default=[s for s in statuses if not s.startswith("Background")])
        show_df(summary[summary["Status"].isin(pick)] if pick else summary)
        st.caption("Hover over a column header and click the (?) for a plain-language explanation.")
        with st.expander("Glossary: what do the columns mean?"):
            for col in summary.columns:
                if gl(col):
                    st.markdown(f"**{col}** - {gl(col)}")

    with tab_det:
        labels = {f["id"]: f"{f['id']} · m/z {f['precursor_mz']:.4f} · RT {f['rt_range']} s · {f['status']}" for f in feats}
        fid = st.selectbox("Feature", list(labels), format_func=labels.get)
        f = next(x for x in feats if x["id"] == fid)
        a, b = st.columns([3, 2])
        with a:
            hit = None
            if f["lib_hits"]:
                idx = st.radio("Compare with library hit", range(len(f["lib_hits"])), horizontal=True,
                               format_func=lambda i: f"#{i + 1} (cosine {f['lib_hits'][i]['score']:.2f})")
                hit = f["lib_hits"][idx]
            st.pyplot(plot_spectrum(f, hit))
            cap, btn = st.columns([14, 1])
            cap.caption("Blue = measured fragments; orange (mirrored) = library spectrum; dashed line = precursor.")
            with btn:
                info("Spectrum plot")
        with b:
            param("File", f["file"], "File")
            param("Precursor m/z", f"{f['precursor_mz']:.5f} ({f['mode']})", "Precursor m/z")
            param("Ion mode", f["mode"], "Mode")
            param("Scans merged", f"{f['n_scans']} (RT {f['rt_range']} s)", "Scans merged")
            param("Adduct", f["adduct"], "Adduct")
            param("Adduct evidence", f["adduct_evidence"], "Adduct evidence")
            param("Neutral mass", f"{f['neutral_mass']:.5f}", "Neutral mass")
            param("Status", f["status"], "Status")
            if f["background"]:
                st.info("This precursor also appears in your blank, so it is treated as background and was not searched.")
            elif f["formulas"]:
                h1, h2 = st.columns([14, 1])
                h1.markdown("**Possible formulas (CHNOPS)**")
                with h2:
                    info("Possible formulas")
                show_df(pd.DataFrame([{"Formula": fm, "Error (ppm)": round(e, 2)} for fm, e in f["formulas"]]))
            else:
                st.write("No CHNOPS formula fits this mass at the chosen tolerance - the adduct may be wrong, or the "
                         "compound may contain other elements.")
        if f["lib_hits"]:
            h1, h2 = st.columns([5, 2])
            h1.subheader("Library matches")
            with h2:
                info_many(G.LIBRARY_HELP)
            cards(f["lib_hits"], "library")
        elif not f["background"] and result["params"]["do_library"]:
            st.write(f"No library match ({f['lib_candidates']} library spectra at this precursor mass).")
        if f["user_hits"]:
            h1, h2 = st.columns([5, 2])
            h1.subheader("Ranked from your candidate list")
            with h2:
                info_many(G.RETRIEVAL_HELP)
            cards(f["user_hits"], "retrieval")
        if f["ret_hits"]:
            h1, h2 = st.columns([5, 2])
            h1.subheader("Ranked from library molecules (fingerprint model)")
            with h2:
                info_many(G.RETRIEVAL_HELP)
            st.caption(f"{f['ret_candidates']} library molecules fit this mass. The ranking shows which are most similar "
                       "to the model's predicted fingerprint. It is a shortlist, not an identification.")
            cards(f["ret_hits"], "retrieval")

    with tab_dl:
        hits = P.hits_table(feats)
        with st.expander("Preview of all hits (hover over headers for explanations)"):
            show_df(hits)
        st.download_button("Feature summary (CSV)", summary.to_csv(index=False).encode(), "features.csv", "text/csv")
        st.download_button("All hits (CSV)", hits.to_csv(index=False).encode(), "hits.csv", "text/csv")
        st.download_button("Everything (Excel)", to_excel({"Features": summary, "Hits": hits}), "results.xlsx",
                           "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet")

    with tab_help:
        st.markdown("""
**Status values**
- *Background*: the precursor also occurs in your blank(s).
- *Library match (strong / moderate / weak)*: a library spectrum with a similar fragmentation pattern. Strong means
  cosine ≥ 0.7 with ≥ 6 shared peaks, moderate ≥ 0.5 with ≥ 4. A match still needs retention time or a reference
  standard to be confirmed.
- *Retrieval candidates only*: structures whose mass fits, ranked by similarity to the model's predicted fingerprint.
  This is a shortlist and the true structure may not be in it.
- *No candidate*: only the formula list is available.

**Limits.** The library covers known compounds only. Adducts are inferred only when a co-eluting partner exists,
otherwise [M+H]+ / [M-H]- is assumed. Formulas use C, H, N, O, P and S only.
""")
