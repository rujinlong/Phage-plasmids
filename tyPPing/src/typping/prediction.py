"""MinProteins and Composition prediction branches."""

import numpy as np
import pandas as pd

from .constants import ALL_CUTOFFS


def _r_round(x: float, digits: int = 0) -> float:
    """Round using R's round-half-up convention."""
    multiplier = 10 ** digits
    return np.floor(x * multiplier + 0.5) / multiplier


def minproteins_branch_complete(
    category: str,
    hmm_filtered: pd.DataFrame,
    minproteins_score: pd.DataFrame,
    minproteins_cutoff: pd.DataFrame,
    genome_size: pd.DataFrame,
) -> pd.DataFrame:
    """MinProteins branch for complete genomes.

    Parameters
    ----------
    category : str
        PP_category code (e.g. "AB", "P11").
    hmm_filtered : pd.DataFrame
        HMM hits filtered for this category.
    minproteins_score : pd.DataFrame
        Profile scores table.
    minproteins_cutoff : pd.DataFrame
        Sequence score thresholds.
    genome_size : pd.DataFrame
        genome_id, size, n_protein.

    Returns
    -------
    pd.DataFrame
        Per-genome MinProteins results.
    """
    # Filter hits passing sequence score threshold
    hmm_mp = hmm_filtered.merge(
        minproteins_cutoff[["hmm_profile", "sequence_score_thr"]],
        on="hmm_profile", how="left",
    )
    hmm_mp = hmm_mp[hmm_mp["sequence_score"] >= hmm_mp["sequence_score_thr"]]

    if hmm_mp.empty:
        return pd.DataFrame(columns=[
            "genome_id", "size", "n_protein", "PP_category", "cutoff_MinProteins",
            "n_hits_MinProteins", "score_0_MinProteins", "score_1_MinProteins", "score_2_MinProteins",
        ])

    cutoff_row = ALL_CUTOFFS[ALL_CUTOFFS["PP_category"] == category]
    category_cutoff = int(cutoff_row["MinProteins_min_N"].iloc[0])

    # Get category profiles
    cat_scores = minproteins_score[minproteins_score["PP_category"] == category][
        ["hmm_profile", "profile_scores_1", "profile_scores_2"]
    ].copy()

    # Build presence/absence matrix: which profiles hit which genomes
    hit_profiles = hmm_mp[["hmm_profile", "genome_id"]].drop_duplicates()
    hit_profiles["values"] = 1
    presence = hit_profiles.pivot_table(index="hmm_profile", columns="genome_id", values="values", fill_value=0)

    # Align with full profile list for this category (full_join)
    all_profiles = cat_scores[["hmm_profile"]].copy()
    presence = all_profiles.merge(presence.reset_index(), on="hmm_profile", how="outer").fillna(0)
    presence = presence.set_index("hmm_profile")

    genome_cols = [c for c in presence.columns]

    # Merge scores
    score_df = all_profiles.merge(cat_scores, on="hmm_profile", how="left").set_index("hmm_profile")

    # Compute weighted scores per genome
    results = []
    for gid in genome_cols:
        pres_vec = presence[gid].values
        n_hits = int(pres_vec.sum())
        s1 = float((pres_vec * score_df["profile_scores_1"].values).sum())
        s2 = float((pres_vec * score_df["profile_scores_2"].values).sum())
        results.append({"genome_id": gid, "n_hits_MinProteins": n_hits,
                        "score_1_MinProteins": s1, "score_2_MinProteins": s2})

    result_df = pd.DataFrame(results)
    result_df = result_df.merge(genome_size, on="genome_id", how="left")
    result_df["PP_category"] = category
    result_df["cutoff_MinProteins"] = category_cutoff
    result_df["size"] = result_df["size"].astype(float)
    result_df["n_protein"] = result_df["n_protein"].astype(float)
    result_df["score_0_MinProteins"] = result_df["n_hits_MinProteins"] / result_df["n_protein"]

    # Round using R convention
    result_df["score_0_MinProteins"] = result_df["score_0_MinProteins"].apply(lambda x: _r_round(x, 3))
    result_df["score_1_MinProteins"] = result_df["score_1_MinProteins"].apply(lambda x: _r_round(x, 1))
    result_df["score_2_MinProteins"] = result_df["score_2_MinProteins"].apply(lambda x: _r_round(x, 1))

    return result_df[["genome_id", "size", "n_protein", "PP_category", "cutoff_MinProteins",
                       "n_hits_MinProteins", "score_0_MinProteins", "score_1_MinProteins", "score_2_MinProteins"]]


def minproteins_branch_draft(
    category: str,
    hmm_filtered: pd.DataFrame,
    minproteins_score: pd.DataFrame,
    minproteins_cutoff: pd.DataFrame,
    contig_size: pd.DataFrame,
    contig_to_genome: pd.DataFrame,
) -> pd.DataFrame:
    """MinProteins branch for draft genomes.

    Works at contig level first, filters by score_0 >= 0.1, then aggregates to genome level.
    """
    # Filter hits passing sequence score threshold
    hmm_mp = hmm_filtered.merge(
        minproteins_cutoff[["hmm_profile", "sequence_score_thr"]],
        on="hmm_profile", how="left",
    )
    hmm_mp = hmm_mp[hmm_mp["sequence_score"] >= hmm_mp["sequence_score_thr"]]

    if hmm_mp.empty:
        return pd.DataFrame(columns=[
            "genome_id", "size_MinProteins", "n_protein_MinProteins", "PP_category",
            "cutoff_MinProteins", "n_hits_MinProteins", "n_hits_MinProteins_list",
            "contig_id_list_MinProteins", "score_0_MinProteins", "score_1_MinProteins",
            "score_2_MinProteins",
        ])

    cutoff_row = ALL_CUTOFFS[ALL_CUTOFFS["PP_category"] == category]
    category_cutoff = int(cutoff_row["MinProteins_min_N"].iloc[0])

    cat_scores = minproteins_score[minproteins_score["PP_category"] == category][
        ["hmm_profile", "profile_scores_1", "profile_scores_2"]
    ].copy()

    # Build presence/absence matrix per contig (not genome)
    hmm_mp_ungrouped = hmm_mp.copy()
    hit_profiles = hmm_mp_ungrouped[["hmm_profile", "contig_id"]].drop_duplicates()
    hit_profiles["values"] = 1
    presence = hit_profiles.pivot_table(index="hmm_profile", columns="contig_id", values="values", fill_value=0)

    all_profiles = cat_scores[["hmm_profile"]].copy()
    presence = all_profiles.merge(presence.reset_index(), on="hmm_profile", how="outer").fillna(0)
    presence = presence.set_index("hmm_profile")

    contig_cols = [c for c in presence.columns]

    score_df = all_profiles.merge(cat_scores, on="hmm_profile", how="left").set_index("hmm_profile")

    # Compute per-contig scores
    contig_results = []
    for cid in contig_cols:
        pres_vec = presence[cid].values
        n_hits = int(pres_vec.sum())
        s1 = float((pres_vec * score_df["profile_scores_1"].values).sum())
        s2 = float((pres_vec * score_df["profile_scores_2"].values).sum())
        contig_results.append({
            "contig_id": cid, "n_hits_MinProteins": n_hits,
            "score_1_MinProteins": s1, "score_2_MinProteins": s2,
        })

    contig_df = pd.DataFrame(contig_results)
    contig_df = contig_df.merge(contig_size[["contig_id", "size", "n_protein"]], on="contig_id", how="left")
    contig_df["PP_category"] = category
    contig_df["cutoff_MinProteins"] = category_cutoff
    contig_df["size"] = contig_df["size"].astype(float)
    contig_df["n_protein"] = contig_df["n_protein"].astype(float)
    contig_df["score_0_MinProteins"] = contig_df["n_hits_MinProteins"] / contig_df["n_protein"]
    contig_df["score_0_MinProteins"] = contig_df["score_0_MinProteins"].apply(lambda x: _r_round(x, 3))
    contig_df["score_1_MinProteins"] = contig_df["score_1_MinProteins"].apply(lambda x: _r_round(x, 1))
    contig_df["score_2_MinProteins"] = contig_df["score_2_MinProteins"].apply(lambda x: _r_round(x, 1))

    # Filter: score_0 >= 0.1
    contig_df = contig_df[contig_df["score_0_MinProteins"] >= 0.1]

    if contig_df.empty:
        return pd.DataFrame(columns=[
            "genome_id", "size_MinProteins", "n_protein_MinProteins", "PP_category",
            "cutoff_MinProteins", "n_hits_MinProteins", "n_hits_MinProteins_list",
            "contig_id_list_MinProteins", "score_0_MinProteins", "score_1_MinProteins",
            "score_2_MinProteins",
        ])

    # Aggregate to genome level
    contig_df = contig_df.merge(contig_to_genome, on="contig_id", how="left")
    contig_df = contig_df.sort_values("n_hits_MinProteins", ascending=False)

    genome_agg = contig_df.groupby(["genome_id", "PP_category", "cutoff_MinProteins"], sort=False).agg(
        n_hits_MinProteins_list=("n_hits_MinProteins", lambda x: ";".join(str(v) for v in x)),
        contig_id_list_MinProteins=("contig_id", lambda x: ";".join(str(v) for v in x)),
        n_protein=("n_protein", "sum"),
        size=("size", "sum"),
        n_hits_MinProteins=("n_hits_MinProteins", "sum"),
        score_1_MinProteins=("score_1_MinProteins", "sum"),
        score_2_MinProteins=("score_2_MinProteins", "sum"),
    ).reset_index()

    genome_agg["score_0_MinProteins"] = genome_agg["n_hits_MinProteins"] / genome_agg["n_protein"]
    genome_agg["score_0_MinProteins"] = genome_agg["score_0_MinProteins"].apply(lambda x: _r_round(x, 3))
    genome_agg["score_1_MinProteins"] = genome_agg["score_1_MinProteins"].apply(lambda x: _r_round(x, 1))
    genome_agg["score_2_MinProteins"] = genome_agg["score_2_MinProteins"].apply(lambda x: _r_round(x, 1))

    genome_agg = genome_agg.rename(columns={"size": "size_MinProteins", "n_protein": "n_protein_MinProteins"})

    return genome_agg[["genome_id", "size_MinProteins", "n_protein_MinProteins", "PP_category",
                        "cutoff_MinProteins", "n_hits_MinProteins", "n_hits_MinProteins_list",
                        "contig_id_list_MinProteins", "score_0_MinProteins", "score_1_MinProteins",
                        "score_2_MinProteins"]]


def composition_branch_complete(
    category: str,
    hmm_filtered: pd.DataFrame,
    compositions: pd.DataFrame,
    genome_size: pd.DataFrame,
) -> pd.DataFrame:
    """Composition branch for complete genomes."""
    cutoff_row = ALL_CUTOFFS[ALL_CUTOFFS["PP_category"] == category]
    tolerance_threshold = int(cutoff_row["composition_tol_thr"].iloc[0])

    # Filter by coverage >= 50%
    hmm_cov = hmm_filtered[hmm_filtered["cov_profile"] >= 0.5]

    # Join with compositions for this category
    comp_cat = compositions[compositions["PP_category"] == category]
    comp_join = comp_cat.merge(hmm_cov, on="hmm_profile", how="inner")

    if comp_join.empty:
        return pd.DataFrame(columns=[
            "genome_id", "size", "n_protein", "PP_category",
            "Composition", "composition_size", "tol_thr", "n_hits_composition",
        ])

    # Count hits per composition per genome
    comp_join["number_hits"] = comp_join.groupby(["Composition", "genome_id"])["hmm_profile"].transform("count")

    # Filter: exact match OR (hits >= size - tolerance AND hits >= 75% of size)
    mask = (
        (comp_join["number_hits"] == comp_join["composition_size"]) |
        ((comp_join["number_hits"] >= comp_join["composition_size"] - tolerance_threshold) &
         (comp_join["number_hits"] >= 0.75 * comp_join["composition_size"]))
    )
    comp_filt = comp_join[mask].copy()

    if comp_filt.empty:
        return pd.DataFrame(columns=[
            "genome_id", "size", "n_protein", "PP_category",
            "Composition", "composition_size", "tol_thr", "n_hits_composition",
        ])

    # Per genome: take the composition with most hits
    # Secondary sort by composition_size desc to match R tie-breaking behavior
    comp_filt = comp_filt.sort_values(["number_hits", "composition_size"], ascending=[False, False])
    comp_filt = comp_filt.groupby("genome_id", sort=False).first().reset_index()

    comp_filt["tol_thr"] = tolerance_threshold
    comp_filt = comp_filt.merge(genome_size, on="genome_id", how="left")
    comp_filt["size"] = comp_filt["size"].astype(float)
    comp_filt["n_protein"] = comp_filt["n_protein"].astype(float)

    # Use PP_category from compositions (PP_category_x)
    pp_col = "PP_category_x" if "PP_category_x" in comp_filt.columns else "PP_category"

    return comp_filt[["genome_id", "size", "n_protein", pp_col,
                       "Composition", "composition_size", "tol_thr", "number_hits"]].rename(
        columns={pp_col: "PP_category", "number_hits": "n_hits_composition"}
    )


def composition_branch_draft(
    category: str,
    hmm_filtered: pd.DataFrame,
    compositions: pd.DataFrame,
    contig_size: pd.DataFrame,
) -> pd.DataFrame:
    """Composition branch for draft genomes."""
    cutoff_row = ALL_CUTOFFS[ALL_CUTOFFS["PP_category"] == category]
    tolerance_threshold = int(cutoff_row["composition_tol_thr"].iloc[0])

    # Filter by coverage >= 0.5, then filter contigs by score_0 >= 0.1
    hmm_cov = hmm_filtered[hmm_filtered["cov_profile"] >= 0.5].copy()
    if hmm_cov.empty:
        return pd.DataFrame(columns=[
            "genome_id", "size_c", "n_protein_c", "PP_category",
            "Composition", "composition_size", "tol_thr", "n_hits_composition",
            "contig_id_list_p",
        ])

    # Count hits per contig and compute score_0_c
    hmm_cov["number_hits_contig"] = hmm_cov.groupby("contig_id")["hmm_profile"].transform("count")
    hmm_cov = hmm_cov.merge(contig_size[["contig_id", "size", "n_protein"]], on="contig_id", how="left",
                             suffixes=("", "_cs"))
    hmm_cov["score_0_c"] = hmm_cov["number_hits_contig"] / hmm_cov["n_protein"]
    hmm_cov = hmm_cov[hmm_cov["score_0_c"] >= 0.1]
    hmm_cov = hmm_cov.drop(columns=["size", "score_0_c", "number_hits_contig", "n_protein"], errors="ignore")

    # Join with compositions
    comp_cat = compositions[compositions["PP_category"] == category]
    comp_join = comp_cat.merge(hmm_cov, on="hmm_profile", how="inner")

    if comp_join.empty:
        return pd.DataFrame(columns=[
            "genome_id", "size_c", "n_protein_c", "PP_category",
            "Composition", "composition_size", "tol_thr", "n_hits_composition",
            "contig_id_list_p",
        ])

    # Count hits per composition per genome
    comp_join["number_hits"] = comp_join.groupby(["Composition", "genome_id"])["hmm_profile"].transform("count")

    # Filter
    mask = (
        (comp_join["number_hits"] == comp_join["composition_size"]) |
        ((comp_join["number_hits"] >= comp_join["composition_size"] - tolerance_threshold) &
         (comp_join["number_hits"] >= 0.75 * comp_join["composition_size"]))
    )
    comp_filt = comp_join[mask].copy()

    if comp_filt.empty:
        return pd.DataFrame(columns=[
            "genome_id", "size_c", "n_protein_c", "PP_category",
            "Composition", "composition_size", "tol_thr", "n_hits_composition",
            "contig_id_list_p",
        ])

    # Per contig: take best composition (most hits)
    comp_filt = comp_filt.sort_values("number_hits", ascending=False)
    comp_filt = comp_filt.groupby("contig_id", sort=False).first().reset_index()

    # Re-join contig_size for size/n_protein
    comp_filt = comp_filt.merge(contig_size[["contig_id", "size", "n_protein"]], on="contig_id", how="left",
                                 suffixes=("", "_cs"))
    comp_filt = comp_filt.sort_values("size", ascending=False)

    # Aggregate to genome level
    comp_filt["tol_thr"] = tolerance_threshold
    genome_agg = comp_filt.groupby("genome_id", sort=False).agg(
        contig_id_list_p=("contig_id", lambda x: ";".join(str(v) for v in x)),
        total_size=("size", "sum"),
        total_n_protein=("n_protein", "sum"),
    ).reset_index()

    # Take the first row per genome for composition info
    first_per_genome = comp_filt.groupby("genome_id", sort=False).first().reset_index()

    pp_col = "PP_category_x" if "PP_category_x" in first_per_genome.columns else "PP_category"
    result = genome_agg.merge(
        first_per_genome[["genome_id", pp_col, "Composition", "composition_size", "tol_thr", "number_hits"]],
        on="genome_id", how="left", suffixes=("", "_fg"),
    )

    # Use tol_thr from first_per_genome if both exist
    if "tol_thr_fg" in result.columns:
        result["tol_thr"] = result["tol_thr_fg"]
        result = result.drop(columns=["tol_thr_fg"])

    return result.rename(columns={
        "total_size": "size_c",
        "total_n_protein": "n_protein_c",
        pp_col: "PP_category",
        "number_hits": "n_hits_composition",
    })[["genome_id", "size_c", "n_protein_c", "PP_category",
        "Composition", "composition_size", "tol_thr", "n_hits_composition", "contig_id_list_p"]]
