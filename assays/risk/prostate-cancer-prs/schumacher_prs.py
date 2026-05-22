ARTICLE_URL = "https://link.springer.com/article/10.1186/s13059-022-02766-z"
PRS_METHOD_NOTE = (
    "Kim et al. 2022 counted prostate-cancer risk allele doses, used proxy r2 to adjust effect sizes, "
    "summed dose times adjusted beta across loci, and filled missing genotypes with study-site mean risk-allele "
    "counts. This implementation reports an observed-only unstandardized score because study-site mean counts "
    "for baseline imputation are not present in the provided Schumacher proxy tables."
)


def row_value(row, key):
    return row[key] if key in row else ""


def split_list(value):
    if value == "":
        return []
    return [item.strip() for item in value.split("|") if item.strip() != ""]


def float_or_none(value):
    if value == "" or value == "NA":
        return None
    return float(value)


def int_or_zero(value):
    if value == "":
        return 0
    return int(float(value))


def index_observations(observations):
    by_rsid = {}
    for obs in observations:
        for key in ["matched_rsid", "rsid"]:
            rsid = row_value(obs, key)
            if rsid != "":
                by_rsid[rsid] = obs
    return by_rsid


def index_variants(variants):
    by_rsid = {}
    for variant in variants:
        rsid = row_value(variant, "rsid")
        if rsid != "":
            by_rsid[rsid] = variant
        for alias in split_list(row_value(variant, "aliases")):
            by_rsid[alias] = variant
    return by_rsid


def variant_rsids(variant):
    rsids = []
    primary = row_value(variant, "rsid")
    if primary != "":
        rsids.append(primary)
    for alias in split_list(row_value(variant, "aliases")):
        if alias not in rsids:
            rsids.append(alias)
    return rsids


def find_observation_for_variant(observations_by_rsid, variant, preferred_rsid):
    lookup_order = []
    if preferred_rsid != "":
        lookup_order.append(preferred_rsid)
    for rsid in variant_rsids(variant):
        if rsid not in lookup_order:
            lookup_order.append(rsid)
    for rsid in lookup_order:
        obs = observations_by_rsid.get(rsid)
        if obs is not None:
            return obs
    return None


def genotype_tokens(genotype):
    if genotype == "":
        return []
    if "/" in genotype:
        return [token for token in genotype.replace("|", "/").split("/") if token != ""]
    if "|" in genotype:
        return [token for token in genotype.split("|") if token != ""]
    if len(genotype) <= 2:
        return list(genotype)
    return [genotype]


def allele_dose(observation, variant, allele):
    if allele == "":
        return None
    genotype = row_value(observation, "genotype")
    tokens = genotype_tokens(genotype)
    if tokens:
        if len(allele) == 1 and all(len(token) == 1 for token in tokens):
            return sum(1 for token in tokens if token == allele)
        return sum(1 for token in tokens if token == allele)

    ref = row_value(variant, "ref")
    alts = split_list(row_value(variant, "alts"))
    alt_count = row_value(observation, "alt_count")
    if allele in alts and alt_count != "":
        return int_or_zero(alt_count)
    if allele == ref and alt_count != "":
        # Diploid fallback only; haploid/ploidy-aware dosage requires genotype tokens.
        return max(0, 2 - int_or_zero(alt_count))
    return None


def choose_candidates(prs_rows):
    by_primary = {}
    for row in prs_rows:
        if row_value(row, "use_in_default_score") != "true":
            continue
        primary = row_value(row, "primary_rsid")
        if primary == "":
            continue
        by_primary.setdefault(primary, []).append(row)
    return by_primary


def score_prs(observations, variants, prs_rows, participant):
    observations_by_rsid = index_observations(observations)
    variants_by_rsid = index_variants(variants)
    candidates = choose_candidates(prs_rows)

    total = 0.0
    primary_found = []
    proxy_found = []
    missing = []
    matched_rows = []

    for primary in candidates:
        chosen = None
        chosen_obs = None
        chosen_variant = None
        chosen_dose = None
        for row in candidates[primary]:
            observed = row_value(row, "observed_rsid")
            variant = variants_by_rsid.get(observed)
            if variant is None:
                continue
            obs = find_observation_for_variant(observations_by_rsid, variant, observed)
            if obs is None:
                continue
            dose = allele_dose(obs, variant, row_value(row, "observed_effect_allele"))
            if dose is None:
                continue
            chosen = row
            chosen_obs = obs
            chosen_variant = variant
            chosen_dose = dose
            break

        if chosen is None:
            missing.append(primary)
            continue

        adjusted_beta = float_or_none(row_value(chosen, "adjusted_beta"))
        if adjusted_beta is None:
            missing.append(primary)
            continue
        contribution = chosen_dose * adjusted_beta
        total += contribution
        detail = (
            primary + "->" + row_value(chosen, "observed_rsid") + ":" +
            row_value(chosen, "observed_effect_allele") + ":dose=" + str(chosen_dose) +
            ":beta=" + row_value(chosen, "adjusted_beta")
        )
        if row_value(chosen, "is_proxy") == "true":
            proxy_found.append(detail)
        else:
            primary_found.append(detail)
        matched_rows.append({
            "participant_id": participant,
            "primary_rsid": primary,
            "observed_rsid": row_value(chosen, "observed_rsid"),
            "is_proxy": row_value(chosen, "is_proxy"),
            "effect_allele": row_value(chosen, "observed_effect_allele"),
            "dose": str(chosen_dose),
            "adjusted_beta": row_value(chosen, "adjusted_beta"),
            "contribution": str(contribution),
            "genotype": row_value(chosen_obs, "genotype"),
            "variant_id": row_value(chosen_variant, "variant_id"),
        })

    return {
        "primary_rsids_found": primary_found,
        "proxies_found": proxy_found,
        "missing": missing,
        "matched_rows": matched_rows,
        "total": total,
    }


try:
    observations_path = bioscript.context["observations_file"]
except Exception:
    observations_path = observations_file

try:
    output_path = bioscript.context["output_file"]
except Exception:
    output_path = output_file

try:
    current_participant_id = bioscript.context["participant_id"]
except Exception:
    current_participant_id = participant_id

try:
    assets = bioscript.context["assets"]
except Exception:
    try:
        assets = bioscript.context["asset_paths"]
    except Exception:
        assets = asset_paths

observations = bioscript.read_tsv(observations_path)
if len(observations) == 0:
    fallback_observations_path = observations_path.replace("test-output/report/data/", "test-output/").replace(".observations.tsv", ".tsv")
    try:
        observations = bioscript.read_tsv(fallback_observations_path)
    except Exception:
        observations = []
variants = bioscript.read_tsv(assets["variants"])
prs_rows = bioscript.read_tsv(assets["prs_metadata"])
result = score_prs(observations, variants, prs_rows, current_participant_id)
bioscript.write_tsv("prostate_cancer_prs_matched_markers.tsv", result["matched_rows"])

bioscript.write_tsv(output_path, [{
    "participant_id": current_participant_id,
    "primary_rsids_found_count": str(len(result["primary_rsids_found"])),
    "primary_rsids_found": "; ".join(result["primary_rsids_found"]),
    "proxies_found_count": str(len(result["proxies_found"])),
    "proxies_found": "; ".join(result["proxies_found"]),
    "missing_baseline_count": str(len(result["missing"])),
    "missing_baseline_markers": ";".join(result["missing"]),
    "baseline_imputation": "not_applied_no_study_site_mean_risk_allele_counts_in_source_tables",
    "total_score": str(result["total"]),
    "score_units": "observed_only_unstandardized_sum_of_effect_allele_dose_times_beta_times_proxy_r2",
    "report_notes": PRS_METHOD_NOTE + " Source: " + ARTICLE_URL,
}])
