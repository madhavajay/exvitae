RS429358 = bioscript.variant(
    rsid=["rs429358", "rs630496", "rs61228756"],
    grch37="19:45411941-45411941",
    grch38="19:44908684-44908684",
    ref="T",
    alt="C",
    kind="snp",
)

RS7412 = bioscript.variant(
    rsid=["rs7412", "rs3200542"],
    grch37="19:45412079-45412079",
    grch38="19:44908822-44908822",
    ref="C",
    alt="T",
    kind="snp",
)

APOE_QUERY_PLAN = bioscript.query_plan([
    RS429358,
    RS7412,
])


def sorted_pair(a, b):
    if a <= b:
        return a + b
    return b + a


def normalize_genotype(gt):
    if gt is None:
        return None
    text = ""
    for ch in gt:
        if ch != "/" and ch != "|" and ch != " " and ch != "-":
            text = text + ch
    if len(text) == 0:
        return None
    if len(text) == 2:
        return sorted_pair(text[0], text[1])
    return text


def display_genotype(gt):
    if gt is None:
        return "missing"
    if gt == "":
        return "missing"
    return gt


def classify_apoe(rs429358_gt, rs7412_gt):
    g429 = normalize_genotype(rs429358_gt)
    g7412 = normalize_genotype(rs7412_gt)

    if g429 is None or g7412 is None:
        return "unresolved_missing_variant"

    # Common unambiguous APOE epsilon-genotype calls from rs429358 + rs7412.
    if g429 == "TT" and g7412 == "TT":
        return "e2/e2"
    if g429 == "TT" and g7412 == "CT":
        return "e2/e3"
    if g429 == "TT" and g7412 == "CC":
        return "e3/e3"
    if g429 == "CT" and g7412 == "CC":
        return "e3/e4"
    if g429 == "CC" and g7412 == "CC":
        return "e4/e4"

    # Double heterozygotes are usually reported as e2/e4 in consumer contexts,
    # but strictly require phase to exclude rare/noncanonical haplotypes.
    if g429 == "CT" and g7412 == "CT":
        return "likely_e2/e4_phase_unconfirmed"

    return "noncanonical_or_requires_phasing"


def apoe_outcome(status):
    if status == "e3/e3":
        return "normal"
    if status == "unresolved_missing_variant":
        return "unknown"
    return "variant"


def report_notes(status, outcome):
    if outcome == "unknown":
        return (
            "APOE epsilon status could not be resolved because rs429358 or rs7412 was missing. "
            "Consult a licensed doctor for advice."
        )
    return (
        'APOE epsilon status is "' + status + '" and is reported as "' + outcome + '" for this assay. '
        "This result is derived from rs429358 and rs7412. Consult a licensed doctor for advice."
    )


def main():
    genotypes = bioscript.load_genotypes(input_file)
    rs429358_gt, rs7412_gt = genotypes.lookup_variants(APOE_QUERY_PLAN)
    apoe_status = classify_apoe(rs429358_gt, rs7412_gt)
    outcome = apoe_outcome(apoe_status)

    rows = [{
        "participant_id": participant_id,
        "apoe_outcome": outcome,
        "apoe_status": apoe_status,
        "rs429358": display_genotype(rs429358_gt),
        "rs7412": display_genotype(rs7412_gt),
        "notes": report_notes(apoe_status, outcome),
    }]
    bioscript.write_tsv(output_file, rows)
    print(apoe_status)


if __name__ == "__main__":
    main()
