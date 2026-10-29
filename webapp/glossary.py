"""Plain-language explanations shown in the app (tooltips on table headers and ℹ️ buttons)."""

GLOSSARY = {
    # ---- summary table / feature info ----
    "ID": "A label (F001, F002, ...) given to each detected ion so you can find it again in the 'Feature detail' tab.",
    "File": "The data file this ion was found in.",
    "Precursor m/z": "The mass-to-charge ratio of the ion the instrument picked for fragmentation. Most ions carry one "
                     "charge, so this is roughly the molecule's mass plus the small particle (adduct) it picked up.",
    "RT (s)": "Retention time: when the compound came off the LC column, in seconds. A range (e.g. 98-104) means the "
              "ion was fragmented several times between those moments. 'n/a' means the file has no retention times.",
    "Scans": "How many separate MS/MS scans of this same precursor were merged into this one feature.",
    "Mode": "Ionisation mode. Positive = the instrument detected positively charged ions; negative = negatively charged ions.",
    "Adduct": "The particle attached to the molecule when it became an ion. [M+H]+ means the molecule (M) plus one "
              "hydrogen ion; [M+Na]+ means plus sodium; [M+NH4]+ plus ammonium; [M-H]- means one hydrogen removed. "
              "It decides the molecule's true mass.",
    "Adduct evidence": "How the adduct was chosen. 'assumed' = a default guess ([M+H]+ in positive mode, [M-H]- in "
                       "negative mode). 'co-eluting pair' = two ions appeared at the same retention time and their "
                       "masses fit two adducts of one molecule, which is stronger evidence.",
    "Neutral mass": "The estimated mass of the original, uncharged molecule: precursor m/z minus the adduct's mass. "
                    "It is used to look for formulas and candidate structures.",
    "Status": "What the analysis found for this ion. Background = also present in your blank (not a sample compound). "
              "Library match = its fragment pattern resembles a known spectrum. Retrieval candidates only = structures "
              "with the right mass, ranked by the model. No candidate = only formulas are available.",
    "Library SMILES": "Structure of the best-matching library spectrum. A SMILES is a line of text that encodes a "
                      "molecule's structure; it can be pasted into most chemistry tools.",
    "Library cosine": "How similar your fragment pattern is to the library's, from 0 (nothing in common) to 1 "
                      "(identical). Above about 0.7 with many shared peaks is a good sign.",
    "Matched peaks": "How many fragment peaks sit at the same m/z in both spectra. More is better: with 3 or fewer, "
                     "even a high cosine is weak evidence.",
    "Library ID": "Identifier of the matched library entry (for example a MoNA ID). Use it to look up the compound's "
                  "name and source.",
    "Retrieval top-1 SMILES": "The structure ranked first among candidates with the right mass, judged by how closely it "
                              "fits the fingerprint predicted by the model. A best guess, not an identification.",
    "Retrieval top-1 formula": "Molecular formula of that top-ranked candidate (which atoms, and how many).",
    "Possible formulas": "Combinations of carbon, hydrogen, nitrogen, oxygen, phosphorus and sulfur whose mass matches "
                         "the neutral mass, with the error in ppm. Many are chemically unlikely; treat them as a "
                         "starting point.",
    # ---- hits table ----
    "Source": "Where the suggestion came from: 'library match' compares spectra with known compounds; 'retrieval' ranks "
              "structures of the right mass by a model-predicted fingerprint.",
    "Rank": "Position in the list for that method. 1 = best.",
    "SMILES": "Text code for the molecule's structure.",
    "Formula": "Molecular formula: which atoms the molecule contains, and how many of each.",
    "Score": "For library matches this is the cosine similarity; for retrieval it is the similarity to the model's "
             "predicted fingerprint. Both run from 0 to 1, but they are not comparable with each other.",
    "Mass diff (ppm)": "How far the candidate's mass is from your measured mass, in parts per million. Smaller is "
                       "better; modern instruments are usually within about 5 ppm.",
    "Confidence": "A rough label. For library matches: strong = cosine 0.7+ with 6+ shared peaks, moderate = 0.5+ with "
                  "4+, otherwise weak. 'rank only' = a model-ranked shortlist with no confidence rating.",
    # ---- formulas table ----
    "Error (ppm)": "Difference between the formula's exact mass and your measured mass, in parts per million. "
                   "Closer to 0 is better.",
    # ---- headline numbers ----
    "MS/MS scans": "Fragmentation spectra read from your file(s). Each one is a scan in which an ion was selected and "
                   "broken into pieces.",
    "Unique precursors": "Distinct ions found after merging repeated scans of the same m/z. Each one is a 'feature'.",
    "Background": "Ions that also appear in your blank sample, so they probably come from solvents, tubing or "
                  "contamination and are skipped.",
    "Library matches": "Features whose fragment pattern matched at least one library spectrum above the cutoffs "
                       "set in the sidebar.",
    "Retrieval only": "Features without a library match but with model-ranked candidate structures of the right mass.",
    # ---- detail view ----
    "Scans merged": "The number of MS/MS scans of this precursor that were combined (the one with the most peaks is shown).",
    "Spectrum plot": "Each vertical line is a fragment ion. Its position is the fragment's m/z and its height the "
                     "relative intensity (strongest = 100%). Blue = your measurement. Orange (mirrored downward) = the "
                     "library spectrum. The dashed line marks the precursor m/z. The more lines match up, the better.",
    "Fingerprint score": "How closely a candidate matches the structural fingerprint predicted from your spectrum "
                         "(0 to 1). Only meaningful for ranking candidates against each other.",
}

# Short sentences for the 'What do these numbers mean?' popovers
LIBRARY_HELP = ["Library cosine", "Matched peaks", "Confidence", "Mass diff (ppm)", "Library ID"]
RETRIEVAL_HELP = ["Fingerprint score", "Mass diff (ppm)", "Formula", "Rank"]
