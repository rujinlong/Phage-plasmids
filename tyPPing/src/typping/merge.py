"""Merge MinProteins and Composition branch results into final predictions."""

import numpy as np
import pandas as pd

from .constants import ALL_CUTOFFS, CATEGORIES_SPECIAL, map_category_to_pptype


def merge_predictions_complete(
    mp_all: pd.DataFrame,
    comp_all: pd.DataFrame,
    genome_size: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Merge results and assign confidence levels for complete genomes.

    Returns
    -------
    all_hits_table : pd.DataFrame
        All HMM hits summary table.
    final_prediction : pd.DataFrame
        Final prediction table with confidence levels.
    """
    # Full join on genome_id + PP_category
    merged = mp_all.merge(comp_all, on=["genome_id", "PP_category"], how="outer", suffixes=("_mp", "_comp"))

    # Coalesce size and n_protein
    if "size_mp" in merged.columns and "size_comp" in merged.columns:
        merged["size"] = merged["size_mp"].fillna(merged["size_comp"])
        merged["n_protein"] = merged["n_protein_mp"].fillna(merged["n_protein_comp"])
        merged = merged.drop(columns=["size_mp", "size_comp", "n_protein_mp", "n_protein_comp"])
    elif "size" not in merged.columns:
        # If only one branch produced results
        for suffix in ["_mp", "_comp"]:
            if f"size{suffix}" in merged.columns:
                merged = merged.rename(columns={f"size{suffix}": "size", f"n_protein{suffix}": "n_protein"})

    # Map PP_category to full P-P type name
    merged["PP_category"] = merged["PP_category"].map(map_category_to_pptype)

    # Build All_hmm_hits_table
    all_hits = merged.rename(columns={
        "genome_id": "Genome ID",
        "PP_category": "P-P type",
        "size": "Genome size (bp)",
        "n_protein": "Number of proteins",
        "cutoff_MinProteins": "MinProteins cutoff",
        "n_hits_MinProteins": "MinProteins hits",
        "score_0_MinProteins": "Signature to all genes ratio",
        "score_1_MinProteins": "P-P score sum",
        "score_2_MinProteins": "P-P type score sum",
        "Composition": "Composition",
        "composition_size": "Composition size",
        "tol_thr": "Tolerance to gaps",
        "n_hits_composition": "Composition hits",
    })

    all_hits_cols = [
        "Genome ID", "P-P type", "MinProteins cutoff", "MinProteins hits",
        "Signature to all genes ratio", "P-P score sum", "P-P type score sum",
        "Composition", "Composition size", "Tolerance to gaps", "Composition hits",
        "Genome size (bp)", "Number of proteins",
    ]
    all_hits = all_hits[[c for c in all_hits_cols if c in all_hits.columns]]

    # Filter for predicted candidates
    pred = merged[
        (merged["n_hits_MinProteins"] >= merged["cutoff_MinProteins"]) |
        (merged["Composition"].notna())
    ].copy()

    if pred.empty:
        final_pred = pd.DataFrame(columns=[
            "Genome ID", "Genome size (bp)", "P-P type", "Confidence level",
            "Predicted by", "MinProteins cutoff", "MinProteins hits",
            "Composition", "Composition size", "Composition hits",
        ])
        return all_hits, final_pred

    # Assign Predicted_by
    pred["Predicted_by"] = np.where(
        (pred["n_hits_MinProteins"] >= pred["cutoff_MinProteins"]) & pred["Composition"].notna(),
        "MinProteins + Composition",
        np.where(
            pred["n_hits_MinProteins"] >= pred["cutoff_MinProteins"],
            "MinProteins",
            "Composition",
        ),
    )

    # Join with cutoffs for size ranges
    cutoffs = ALL_CUTOFFS[["PP_type", "size_min", "size_max", "pct10_mean_size"]].copy()
    pred = pred.merge(cutoffs, left_on="PP_category", right_on="PP_type", how="left")

    # Assign confidence levels
    def _assign_confidence(row):
        if (row["size"] <= row["size_max"] and row["size"] >= row["size_min"]
                and row["Predicted_by"] == "MinProteins + Composition"):
            return "High"
        if row["PP_category"] in CATEGORIES_SPECIAL and row["Predicted_by"] == "Composition":
            return "Low"
        if (row["size"] > row["pct10_mean_size"] + row["size_max"]
                or row["size"] < row["size_min"] - row["pct10_mean_size"]):
            return "Low"
        return "Medium"

    pred["Confidence"] = pred.apply(_assign_confidence, axis=1)

    # Scoring for dedup of special categories
    pred["prediction_score"] = pred["Confidence"].map({"High": 3, "Medium": 2, "Low": 1})
    pred["type_score"] = pred["PP_category"].map({"pSLy3": 2, "pKpn": 1}).fillna(3).astype(int)

    # Split into non-special and special
    pred_other = pred[~pred["PP_category"].isin(CATEGORIES_SPECIAL)].copy()

    pred_special = pred[pred["PP_category"].isin(CATEGORIES_SPECIAL)].copy()
    if not pred_special.empty:
        pred_special = pred_special.sort_values(
            ["prediction_score", "type_score"], ascending=[False, False]
        )
        pred_special = pred_special.groupby("genome_id", sort=False).first().reset_index()

    final = pd.concat([pred_other, pred_special], ignore_index=True)

    # Re-evaluate confidence for special categories after dedup
    def _reassign_special(row):
        if (row["PP_category"] in CATEGORIES_SPECIAL
                and row["Predicted_by"] == "Composition"
                and row["size"] < row["pct10_mean_size"] + row["size_max"]
                and row["size"] > row["size_min"] - row["pct10_mean_size"]):
            return "Medium"
        return row["Confidence"]

    final["Confidence"] = final.apply(_reassign_special, axis=1)

    final_pred = final.rename(columns={
        "genome_id": "Genome ID",
        "size": "Genome size (bp)",
        "PP_category": "P-P type",
        "Confidence": "Confidence level",
        "Predicted_by": "Predicted by",
        "cutoff_MinProteins": "MinProteins cutoff",
        "n_hits_MinProteins": "MinProteins hits",
        "composition_size": "Composition size",
        "n_hits_composition": "Composition hits",
    })

    pred_cols = [
        "Genome ID", "Genome size (bp)", "P-P type", "Confidence level",
        "Predicted by", "MinProteins cutoff", "MinProteins hits",
        "Composition", "Composition size", "Composition hits",
    ]
    final_pred = final_pred[[c for c in pred_cols if c in final_pred.columns]]

    return all_hits, final_pred


def merge_predictions_draft(
    mp_all: pd.DataFrame,
    comp_all: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Merge results for draft genomes.

    Returns
    -------
    all_hits_table : pd.DataFrame
    final_prediction : pd.DataFrame
    """
    # Full join
    merged = mp_all.merge(comp_all, on=["genome_id", "PP_category"], how="outer")

    # Map PP_category
    merged["PP_category"] = merged["PP_category"].map(map_category_to_pptype)

    # N15 special cutoff: when no composition, set cutoff to 24
    n15_no_comp = (merged["PP_category"] == "N15") & merged["Composition"].isna()
    merged.loc[n15_no_comp, "cutoff_MinProteins"] = 24

    # Build all_hits_table
    all_hits = merged.rename(columns={
        "genome_id": "Genome ID",
        "PP_category": "P-P type",
        "size_MinProteins": "Genome size MinProteins (bp)",
        "n_protein_MinProteins": "Number of proteins MinProteins",
        "cutoff_MinProteins": "MinProteins cutoff",
        "n_hits_MinProteins": "MinProteins hits",
        "n_hits_MinProteins_list": "MinProteins hits list",
        "contig_id_list_MinProteins": "MinProteins contigs list",
        "score_0_MinProteins": "Signature to all genes ratio",
        "score_1_MinProteins": "P-P score sum",
        "score_2_MinProteins": "P-P type score sum",
        "size_c": "Genome size Composition (bp)",
        "n_protein_c": "Number of proteins Composition",
        "Composition": "Composition",
        "composition_size": "Composition size",
        "tol_thr": "Tolerance to gaps",
        "n_hits_composition": "Composition hits",
        "contig_id_list_p": "Composition contigs list",
    })

    all_hits_cols = [
        "Genome ID", "P-P type",
        "Genome size MinProteins (bp)", "Number of proteins MinProteins",
        "MinProteins cutoff", "MinProteins hits", "MinProteins hits list", "MinProteins contigs list",
        "Signature to all genes ratio", "P-P score sum", "P-P type score sum",
        "Genome size Composition (bp)", "Number of proteins Composition",
        "Composition", "Composition size", "Tolerance to gaps", "Composition hits",
        "Composition contigs list",
    ]
    all_hits = all_hits[[c for c in all_hits_cols if c in all_hits.columns]]

    # Filter for predicted candidates
    pred = merged[
        (merged["n_hits_MinProteins"] >= merged["cutoff_MinProteins"]) |
        (merged["Composition"].notna())
    ].copy()

    if pred.empty:
        final_pred = pd.DataFrame(columns=[
            "Genome ID", "P-P type", "Confidence level", "Predicted by",
            "Genome size MinProteins (bp)", "MinProteins cutoff", "MinProteins hits",
            "MinProteins hits list", "MinProteins contigs list",
            "Genome size Composition (bp)", "Composition", "Composition size",
            "Composition hits", "Composition contigs list",
        ])
        return all_hits, final_pred

    # Predicted_by
    pred["Predicted_by"] = np.where(
        (pred["n_hits_MinProteins"] >= pred["cutoff_MinProteins"]) & pred["Composition"].notna(),
        "MinProteins + Composition",
        np.where(
            pred["n_hits_MinProteins"] >= pred["cutoff_MinProteins"],
            "MinProteins",
            "Composition",
        ),
    )

    # Confidence levels (same logic as complete, but using size_MinProteins)
    cutoffs = ALL_CUTOFFS[["PP_type", "size_min", "size_max", "pct10_mean_size"]].copy()
    pred = pred.merge(cutoffs, left_on="PP_category", right_on="PP_type", how="left")

    def _assign_confidence(row):
        size = row.get("size_MinProteins", np.nan)
        if pd.isna(size):
            size = row.get("size_c", 0)
        if (size <= row["size_max"] and size >= row["size_min"]
                and row["Predicted_by"] == "MinProteins + Composition"):
            return "High"
        if row["PP_category"] in CATEGORIES_SPECIAL and row["Predicted_by"] == "Composition":
            return "Low"
        if (size > row["pct10_mean_size"] + row["size_max"]
                or size < row["size_min"] - row["pct10_mean_size"]):
            return "Low"
        return "Medium"

    pred["Confidence"] = pred.apply(_assign_confidence, axis=1)

    pred["prediction_score"] = pred["Confidence"].map({"High": 3, "Medium": 2, "Low": 1})
    pred["type_score"] = pred["PP_category"].map({"pSLy3": 2, "pKpn": 1}).fillna(3).astype(int)

    pred_other = pred[~pred["PP_category"].isin(CATEGORIES_SPECIAL)].copy()
    pred_special = pred[pred["PP_category"].isin(CATEGORIES_SPECIAL)].copy()
    if not pred_special.empty:
        pred_special = pred_special.sort_values(
            ["prediction_score", "type_score"], ascending=[False, False]
        )
        pred_special = pred_special.groupby("genome_id", sort=False).first().reset_index()

    final = pd.concat([pred_other, pred_special], ignore_index=True)

    def _reassign_special(row):
        size = row.get("size_MinProteins", np.nan)
        if pd.isna(size):
            size = row.get("size_c", 0)
        if (row["PP_category"] in CATEGORIES_SPECIAL
                and row["Predicted_by"] == "Composition"
                and size < row["pct10_mean_size"] + row["size_max"]
                and size > row["size_min"] - row["pct10_mean_size"]):
            return "Medium"
        return row["Confidence"]

    final["Confidence"] = final.apply(_reassign_special, axis=1)

    final_pred = final.rename(columns={
        "genome_id": "Genome ID",
        "PP_category": "P-P type",
        "Confidence": "Confidence level",
        "Predicted_by": "Predicted by",
        "size_MinProteins": "Genome size MinProteins (bp)",
        "cutoff_MinProteins": "MinProteins cutoff",
        "n_hits_MinProteins": "MinProteins hits",
        "n_hits_MinProteins_list": "MinProteins hits list",
        "contig_id_list_MinProteins": "MinProteins contigs list",
        "size_c": "Genome size Composition (bp)",
        "composition_size": "Composition size",
        "n_hits_composition": "Composition hits",
        "contig_id_list_p": "Composition contigs list",
    })

    pred_cols = [
        "Genome ID", "P-P type", "Confidence level", "Predicted by",
        "Genome size MinProteins (bp)", "MinProteins cutoff", "MinProteins hits",
        "MinProteins hits list", "MinProteins contigs list",
        "Composition", "Genome size Composition (bp)", "Composition size",
        "Composition hits", "Composition contigs list",
    ]
    final_pred = final_pred[[c for c in pred_cols if c in final_pred.columns]]

    return all_hits, final_pred
