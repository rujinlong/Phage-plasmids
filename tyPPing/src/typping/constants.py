"""Constants and cutoff tables for tyPPing."""

import pandas as pd

# Column names for HMMER domtblout format (23 columns, last is free-text description)
HMM_DOMTBL_COLUMNS = [
    "domain_name", "domain_accession", "domain_len",
    "query_name", "query_accession", "qlen",
    "sequence_evalue", "sequence_score", "sequence_bias",
    "domain_N", "domain_of",
    "domain_cevalue", "domain_ievalue", "domain_score", "domain_bias",
    "hmm_from", "hmm_to", "ali_from", "ali_to",
    "env_from", "env_to", "acc", "description",
]

# Map short PP_category codes extracted from HMM profile names back to full P-P type names
PP_CATEGORY_REVERSE = {
    "AB": "AB_1",
    "P11": "P1_1",
    "P12": "P1_2",
    "SSU5": "SSU5_pHCM2",
}

# Map full P-P type names to short PP_category codes
PP_TYPE_TO_CATEGORY = {
    "AB_1": "AB",
    "P1_1": "P11",
    "P1_2": "P12",
    "SSU5_pHCM2": "SSU5_pHCM2",
    "cp32": "cp32",
    "N15": "N15",
    "pCAV": "pCAV",
    "pKpn": "pKpn",
    "pMT1": "pMT1",
    "pSLy3": "pSLy3",
}

# Special categories that share genome space and need deduplication
CATEGORIES_SPECIAL = {"cp32", "SSU5_pHCM2", "pKpn", "pSLy3"}

# Cutoffs table — one row per P-P type
ALL_CUTOFFS = pd.DataFrame({
    "PP_type": ["AB_1", "cp32", "N15", "P1_1", "P1_2", "pCAV", "pKpn", "pMT1", "pSLy3", "SSU5_pHCM2"],
    "composition_tol_thr": [7, 6, 9, 11, 17, 3, 2, 6, 3, 3],
    "MinProteins_min_N": [55, 18, 7, 49, 46, 30, 38, 50, 47, 71],
    "size_min": [101329, 27653, 39838, 72057, 79071, 102639, 98807, 85066, 85190, 94107],
    "size_max": [122444, 32756, 65948, 125381, 103576, 118316, 125007, 115126, 133911, 125751],
    "pct10_mean_size": [11189, 3020, 5289, 9872, 9132, 11048, 11191, 10010, 10955, 10993],
    "PP_category": ["AB", "cp32", "N15", "P11", "P12", "pCAV", "pKpn", "pMT1", "pSLy3", "SSU5_pHCM2"],
})


def map_category_to_pptype(category: str) -> str:
    """Map a PP_category code to the full P-P type name."""
    return PP_CATEGORY_REVERSE.get(category, category)


def map_pptype_to_category(pptype: str) -> str:
    """Map a full P-P type name to the short PP_category code."""
    return PP_TYPE_TO_CATEGORY.get(pptype, pptype)
