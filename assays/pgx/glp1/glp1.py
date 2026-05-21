MARKERS = [
    {
        "rsid": "rs10305420",
        "gene": "GLP1R",
        "path_token": "rs10305420",
        "effect_allele": "T",
        "study_finding": "The 2026 23andMe Nature GWAS found that people with the rs10305420 T allele in GLP1R p.Pro7Leu had greater GLP-1 receptor agonist weight-loss response, about 0.76 kg greater weight loss per effect allele in the GLP1RA cohort.",
        "present_result": "You have the rs10305420 T allele that the study associated with greater GLP-1 receptor agonist weight-loss response.",
        "absent_result": "You do not have the rs10305420 T allele that the study associated with greater GLP-1 receptor agonist weight-loss response.",
    },
    {
        "rsid": "rs11760106",
        "gene": "GLP1R",
        "path_token": "rs11760106",
        "effect_allele": "T",
        "study_finding": "The 2026 23andMe Nature GWAS found that people with the rs11760106 T allele at the GLP1R locus had increased moderate-to-severe vomiting risk during GLP-1 therapy in the GLP1RA cohort, with a reported odds ratio of about 1.57.",
        "present_result": "You have the rs11760106 T allele that the study associated with increased moderate-to-severe vomiting risk during GLP-1 therapy.",
        "absent_result": "You do not have the rs11760106 T allele that the study associated with increased moderate-to-severe vomiting risk during GLP-1 therapy.",
    },
    {
        "rsid": "rs9357296",
        "gene": "GLP1R",
        "path_token": "rs9357296",
        "effect_allele": "G",
        "study_finding": "The 2026 23andMe Nature GWAS found that people with the rs9357296 G allele at the GLP1R locus had increased nausea risk during GLP-1 therapy in the GLP1RA cohort, with a reported odds ratio of about 1.36.",
        "present_result": "You have the rs9357296 G allele that the study associated with increased nausea risk during GLP-1 therapy.",
        "absent_result": "You do not have the rs9357296 G allele that the study associated with increased nausea risk during GLP-1 therapy.",
    },
    {
        "rsid": "rs1800437",
        "gene": "GIPR",
        "path_token": "rs1800437",
        "effect_allele": "C",
        "study_finding": "The 2026 23andMe Nature GWAS found that people with the rs1800437 C allele in GIPR p.Glu354Gln had increased vomiting or nausea risk in tirzepatide-treated individuals, with a reported odds ratio of about 1.83. The paper prioritized this missense variant as the probable causal GIPR signal.",
        "present_result": "You have the rs1800437 C allele that the study associated with increased vomiting or nausea risk in tirzepatide-treated individuals.",
        "absent_result": "You do not have the rs1800437 C allele that the study associated with increased vomiting or nausea risk in tirzepatide-treated individuals.",
    },
    {
        "rsid": "rs71338792",
        "gene": "GIPR",
        "path_token": "rs71338792",
        "effect_allele": "T[14]",
        "effect_sequence": "TTTTTTTTTTTTTT",
        "study_finding": "The 2026 23andMe Nature GWAS reported rs71338792 as a GIPR locus tagging marker for the side-effect signal. It is in high linkage disequilibrium with rs1800437, which the paper prioritized as the more likely causal variant.",
        "present_result": "You have the rs71338792 T[14] repeat allele represented as rs71338792:AT in the study, so this GIPR tagging marker is present. The paper still prioritizes rs1800437 as the more likely causal variant at this locus.",
        "absent_result": "You do not have the rs71338792 T[14] repeat allele for this GIPR tagging marker.",
    },
]

GENERAL_NOTE = (
    "This interpretation summarizes cohort-level associations from the 2026 23andMe Nature paper GWAS study. "
    "It is research-use genetic context, not a treatment recommendation."
)


def narrative_info(marker, call):
    genotype = call["genotype"]
    return (
        f"Genotype: {genotype}. "
        f"Finding: {marker['study_finding']} "
        f"Your result: {call['your_result']}"
    )


def observation_for_marker(rows, marker):
    token = marker["path_token"]
    for row in rows:
        path = row.get("path") or ""
        rsid = row.get("matched_rsid") or ""
        name = row.get("name") or ""
        if rsid == marker["rsid"] or token in path or token in name:
            return row
    return None


def normalize_snv_genotype(genotype):
    if genotype is None:
        return ""
    text = ""
    for ch in genotype:
        if ch in "ACGT":
            text = text + ch
    return text


def effect_allele_count(genotype, allele):
    text = normalize_snv_genotype(genotype)
    if text == "":
        return None
    return text.count(allele)


def has_indel_alt(row, genotype):
    try:
        if int(row.get("alt_count") or "0") > 0:
            return True
    except ValueError:
        pass
    evidence = row.get("evidence") or ""
    if "matched_alts:" in evidence:
        return True
    if genotype is None or genotype == "":
        return False
    if "/" in genotype or "|" in genotype:
        parts = genotype.replace("|", "/").split("/")
        if len(parts) > 1:
            first = parts[0]
            for part in parts:
                if part != first:
                    return True
    return False


def indel_effect_allele_count(genotype, effect_sequence):
    if genotype is None or genotype == "":
        return None
    if "/" in genotype or "|" in genotype:
        parts = genotype.replace("|", "/").split("/")
        return sum(1 for part in parts if part == effect_sequence)
    if "I" in genotype or "D" in genotype:
        return None
    return 0


def is_missing_genotype(genotype):
    if genotype is None:
        return True
    if genotype == "" or genotype == "missing" or genotype == "./." or genotype == ".|." or genotype == "??":
        return True
    return False


def call_marker(marker, row):
    if row is None:
        return {
            "genotype": "missing",
            "variant_status": "unknown",
            "glp1_outcome": "unknown",
            "effect_allele_count": "",
            "your_result": "No observation was available for this marker.",
            "interpretation": "Your result for this marker is unknown, so this report cannot say whether you have the study-associated allele.",
        }

    genotype = row.get("genotype") or "missing"
    if marker["rsid"] == "rs71338792":
        if is_missing_genotype(row.get("genotype")):
            return {
                "genotype": genotype,
                "variant_status": "unknown",
                "glp1_outcome": "unknown",
                "effect_allele_count": "",
                "your_result": "No genotype call was available for this indel tagging marker.",
                "interpretation": "Your result for this marker is unknown, so this report cannot say whether you have the GIPR tagging marker.",
            }
        count = indel_effect_allele_count(row.get("genotype"), marker["effect_sequence"])
        if count is None:
            present = has_indel_alt(row, row.get("genotype") or "")
            count = 1 if present else 0
        if count > 0:
            return {
                "genotype": genotype,
                "variant_status": "present",
                "glp1_outcome": "variant",
                "effect_allele_count": str(count),
                "your_result": marker["present_result"],
                "interpretation": marker["present_result"],
            }
        return {
            "genotype": genotype,
            "variant_status": "not_observed",
            "glp1_outcome": "normal",
            "effect_allele_count": "0",
            "your_result": marker["absent_result"],
            "interpretation": marker["absent_result"],
        }

    count = effect_allele_count(row.get("genotype"), marker["effect_allele"])
    if count is None:
        return {
            "genotype": genotype,
            "variant_status": "unknown",
            "glp1_outcome": "unknown",
            "effect_allele_count": "",
            "your_result": "No genotype call was available for this marker.",
            "interpretation": "Your result for this marker is unknown, so this report cannot say whether you have the study-associated allele.",
        }
    if count > 0:
        return {
            "genotype": genotype,
            "variant_status": "present",
            "glp1_outcome": "variant",
            "effect_allele_count": str(count),
            "your_result": marker["present_result"],
            "interpretation": marker["present_result"],
        }
    return {
        "genotype": genotype,
        "variant_status": "not_observed",
        "glp1_outcome": "normal",
        "effect_allele_count": "0",
        "your_result": marker["absent_result"],
        "interpretation": marker["absent_result"],
    }


def main():
    observations_path = bioscript.context["observations_file"]
    rows = bioscript.read_tsv(observations_path)
    output_rows = []
    for marker in MARKERS:
        row = observation_for_marker(rows, marker)
        call = call_marker(marker, row)
        output_rows.append({
            "participant_id": participant_id,
            "rsid": marker["rsid"],
            "info": narrative_info(marker, call),
            "notes": GENERAL_NOTE,
        })
    bioscript.write_tsv(output_file, output_rows)
    print("glp1_nature_findings")


if __name__ == "__main__":
    main()
