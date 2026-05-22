def sorted_pair(a, b):
    if a <= b:
        return a + b
    return b + a


def normalize_genotype(gt):
    if gt is None:
        return None
    text = ""
    for ch in gt:
        if ch != "/" and ch != "|" and ch != " " and ch != "-" and ch != ".":
            text = text + ch
    if len(text) == 0:
        return None
    if len(text) == 2:
        return sorted_pair(text[0], text[1])
    return text


def display_genotype(gt):
    if gt is None:
        return "missing"
    if gt == "" or gt == "." or gt == "./.":
        return "missing"
    return gt


def classify_apoe(rs429358_gt, rs7412_gt):
    g429 = normalize_genotype(rs429358_gt)
    g7412 = normalize_genotype(rs7412_gt)

    if g429 is None or g7412 is None:
        return "unresolved_missing_variant"

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
    if g429 == "CT" and g7412 == "CT":
        return "likely_e2/e4_phase_unconfirmed"

    return "noncanonical_or_requires_phasing"


def apoe_outcome(status):
    if status == "unresolved_missing_variant":
        return "unknown"
    if status == "e3/e3":
        return "normal"
    return "variant"


def apoe_longevity_status(status):
    if status == "e2/e2" or status == "e2/e3":
        return "good"
    if status == "e3/e3":
        return "normal"
    if status == "e3/e4" or status == "e4/e4":
        return "less_good"
    if status == "likely_e2/e4_phase_unconfirmed":
        return "mixed"
    return "unknown"


def apoe_longevity_context():
    return (
        "Hirose et al. (Tokyo Centenarian Study, PMID 9212680) reported that APOE epsilon2 was more frequent and "
        "epsilon4 was less frequent in centenarians than controls, and concluded that epsilon2 was positively and "
        "epsilon4 negatively associated with longevity."
    )


def apoe_longevity_summary(status):
    source = apoe_longevity_context()
    if status == "e2/e2":
        return (
            "APOE status e2/e2 carries two epsilon2 alleles. " + source +
            " This is summarized here as a favorable longevity-associated APOE pattern."
        )
    if status == "e2/e3":
        return (
            "APOE status e2/e3 carries one epsilon2 allele and one common epsilon3 allele. " + source +
            " This is summarized here as a favorable longevity-associated APOE pattern."
        )
    if status == "e3/e3":
        return (
            "APOE status e3/e3 carries two epsilon3 alleles. In the Tokyo Centenarian Study, epsilon3 was the most common "
            "allele in both centenarians and controls; this assay treats e3/e3 as the reference/common APOE longevity pattern."
        )
    if status == "e3/e4":
        return (
            "APOE status e3/e4 carries one epsilon4 allele. " + source +
            " This is summarized here as a less favorable APOE longevity pattern than e3/e3 or epsilon2-containing genotypes."
        )
    if status == "e4/e4":
        return (
            "APOE status e4/e4 carries two epsilon4 alleles. " + source +
            " This is summarized here as a less favorable APOE longevity pattern."
        )
    if status == "likely_e2/e4_phase_unconfirmed":
        return (
            "APOE status is likely e2/e4 from rs429358 and rs7412, but phase is unconfirmed. The genotype contains signals "
            "consistent with both epsilon2 and epsilon4; Hirose et al. reported epsilon2 positively and epsilon4 negatively "
            "associated with longevity, so this assay summarizes the result as mixed."
        )
    if status == "unresolved_missing_variant":
        return "APOE epsilon status could not be resolved because rs429358 or rs7412 was missing."
    return "APOE epsilon status is noncanonical or requires phasing; no longevity summary is assigned."


def report_notes(status, outcome):
    if outcome == "unknown":
        return apoe_longevity_summary(status) + " Consult a licensed doctor for advice."
    return (
        'APOE epsilon status is "' + status + '" and is reported as "' + outcome + '" for this assay. '
        + apoe_longevity_summary(status) + " Consult a licensed doctor for advice."
    )


def row_matches_marker(row, marker):
    rsid = row.get("rsid") or row.get("matched_rsid") or row.get("id") or row.get("marker")
    if rsid == marker:
        return True
    variant_key = row.get("variant_key") or row.get("name") or row.get("variant") or row.get("path") or ""
    return variant_key == "APOE_" + marker or variant_key.endswith(marker) or marker in variant_key


def row_genotype(row):
    for key in ["genotype_display", "genotype", "call", "gt"]:
        value = row.get(key)
        if value is not None and value != "" and value != "." and value != "./.":
            return value
    return None


def observation_genotypes():
    rows = bioscript.read_tsv(observations_file)
    values = {}
    for marker in ["rs429358", "rs7412"]:
        for row in rows:
            if row_matches_marker(row, marker):
                genotype = row_genotype(row)
                if genotype is not None:
                    values[marker] = genotype
                    break
    return values.get("rs429358"), values.get("rs7412")


def main():
    rs429358_gt, rs7412_gt = observation_genotypes()
    apoe_status = classify_apoe(rs429358_gt, rs7412_gt)
    outcome = apoe_outcome(apoe_status)
    longevity_status = apoe_longevity_status(apoe_status)
    longevity_summary = apoe_longevity_summary(apoe_status)
    longevity_context = apoe_longevity_context()

    rows = [{
        "participant_id": participant_id,
        "apoe_outcome": outcome,
        "apoe_status": apoe_status,
        "apoe_longevity_status": longevity_status,
        "apoe_longevity_context": longevity_context,
        "apoe_longevity_summary": longevity_summary,
        "rs429358": display_genotype(rs429358_gt),
        "rs7412": display_genotype(rs7412_gt),
        "notes": report_notes(apoe_status, outcome),
    }]
    bioscript.write_tsv(output_file, rows)
    print(apoe_status)


if __name__ == "__main__":
    main()
