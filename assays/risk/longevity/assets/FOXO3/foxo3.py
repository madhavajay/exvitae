MARKERS = [
    {
        "rsid": "rs2802292",
        "beneficial_allele": "G",
        "label": "FOXO3 rs2802292",
        "study_note": (
            "You carry {copies} of the rs2802292 G allele. Zeng et al. reported rs2802292*G was associated with "
            "longevity in Asians in meta-analysis; Torigoe et al. reported protective cellular-aging associations for "
            "rs2802292 G-allele carriers."
        ),
    },
    {
        "rsid": "rs2802288",
        "beneficial_allele": "A",
        "label": "FOXO3 rs2802288",
        "study_note": (
            "You carry {copies} of the rs2802288 A allele. Zeng et al. reported rs2802288*A was associated with "
            "longevity in Southern Chinese and remained significant after Bonferroni correction; their meta-analysis "
            "also supported rs2802288*A."
        ),
    },
    {
        "rsid": "rs13217795",
        "beneficial_allele": "C",
        "label": "FOXO3 rs13217795",
        "study_note": (
            "You carry {copies} of the rs13217795 C allele. Zeng et al. reported rs13217795*C had a higher minor "
            "allele frequency in longevity subjects and remained significant after Bonferroni correction in their "
            "Southern Chinese case-control study."
        ),
    },
]


def normalize_genotype(gt):
    if gt is None:
        return None
    text = ""
    for ch in gt:
        if ch in "ACGT":
            text = text + ch
    if text == "":
        return None
    return text


def dosage(gt, allele):
    text = normalize_genotype(gt)
    if text is None:
        return None
    return text.count(allele)


def outcome_for_dosage(d):
    if d is None:
        return "unknown"
    if d >= 1:
        return "variant"
    return "normal"


def copies_text(d):
    if d is None:
        return "an unknown number of copies"
    if d == 1:
        return "1 copy"
    return str(d) + " copies"


def overall_for_rows(rows):
    if any(row["foxo3_outcome"] == "variant" for row in rows):
        return "variant"
    if any(row["foxo3_outcome"] == "unknown" for row in rows):
        return "unknown"
    return "normal"


def row_matches_marker(row, marker):
    rsid = row.get("rsid") or row.get("id") or row.get("marker")
    if rsid == marker:
        return True
    variant_key = row.get("variant_key") or row.get("name") or row.get("variant") or row.get("path") or ""
    return variant_key == "FOXO3_" + marker or variant_key.endswith(marker) or marker in variant_key


def row_genotype(row):
    for key in ["genotype_display", "genotype", "call", "gt"]:
        value = row.get(key)
        if value is not None and value != "" and value != "." and value != "./.":
            return value
    return None


def observation_map():
    rows = bioscript.read_tsv(observations_file)
    values = {}
    for marker in MARKERS:
        marker_id = marker["rsid"]
        for row in rows:
            if row_matches_marker(row, marker_id):
                genotype = row_genotype(row)
                if genotype is not None:
                    values[marker_id] = genotype
                    break
    return values


def main():
    observations = observation_map()
    marker_rows = []
    for marker in MARKERS:
        gt = observations.get(marker["rsid"])
        d = dosage(gt, marker["beneficial_allele"])
        outcome = outcome_for_dosage(d)
        marker_rows.append({
            "participant_id": participant_id,
            "foxo3_outcome": outcome,
            "rsid": marker["rsid"],
            "genotype": gt if gt is not None and gt != "." else "missing",
            "notes": marker["study_note"].replace("{copies}", copies_text(d)),
        })

    bioscript.write_tsv(output_file, marker_rows)
    print(overall_for_rows(marker_rows))


if __name__ == "__main__":
    main()
