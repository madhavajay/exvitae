MARKERS = [
    {
        "rsid": "rs11591147",
        "name": "R46L / Arg46Leu",
        "path_token": "rs11591147",
        "protective_allele": "T",
        "present": "PCSK9 R46L protective allele observed.",
        "absent": "PCSK9 R46L protective allele not observed.",
    },
    {
        "rsid": "rs67608943",
        "name": "Y142X / Tyr142Ter",
        "path_token": "rs67608943",
        "protective_allele": "G",
        "present": "PCSK9 Y142X protective nonsense allele observed.",
        "absent": "PCSK9 Y142X protective nonsense allele not observed.",
    },
    {
        "rsid": "rs28362286",
        "name": "C679X / Cys679Ter",
        "path_token": "rs28362286",
        "protective_allele": "A",
        "present": "PCSK9 C679X protective nonsense allele observed.",
        "absent": "PCSK9 C679X protective nonsense allele not observed.",
    },
]

GENERAL_NOTE = (
    "These are naturally occurring PCSK9 LDL-lowering alleles. In the Cohen et al. 2006 population study, carrier groups "
    "had lower LDL cholesterol and lower coronary heart disease risk. In this narrow research context they are generally "
    "favorable, but this report is not a diagnosis, a full cardiovascular risk estimate, or medication advice."
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


def protective_allele_count(row, marker):
    if row is None:
        return None
    genotype = row.get("genotype")
    text = normalize_snv_genotype(genotype)
    if text == "":
        return None
    return text.count(marker["protective_allele"])


def genotype_or_missing(row):
    if row is None:
        return "missing"
    genotype = row.get("genotype")
    if genotype is None or genotype == "":
        return "missing"
    return genotype


def classify(total_protective, missing_count):
    if total_protective == 0 and missing_count == len(MARKERS):
        return "unknown"
    if total_protective == 0:
        return "typical"
    if total_protective == 1:
        return "protective_carrier"
    return "protective_homozygous_or_multiple"


def status_label(status):
    if status == "typical":
        return "Typical PCSK9 result"
    if status == "protective_carrier":
        return "PCSK9 LDL-lowering allele carrier"
    if status == "protective_homozygous_or_multiple":
        return "Multiple PCSK9 LDL-lowering allele copies observed"
    return "PCSK9 result unknown"


def interpretation(status, total_protective, details):
    joined = " ".join(details)
    if status == "typical":
        return (
            "No curated PCSK9 LDL-lowering protective allele was observed in the available data. "
            "This is the usual result for these specific markers and does not rule out other lipid-related genetics. "
            + joined
        )
    if status == "unknown":
        return (
            "No usable genotype call was available for these PCSK9 markers, so this report cannot classify carrier status. "
            + joined
        )
    return (
        f"{total_protective} curated PCSK9 LDL-lowering protective allele copy/copies were observed. "
        "For these specific markers, the variant allele is generally the favorable finding because it is associated with "
        "lower LDL-C and lower coronary heart disease risk in the cited population studies. "
        + joined
    )


def main():
    observations_path = bioscript.context["observations_file"]
    rows = bioscript.read_tsv(observations_path)

    total_protective = 0
    missing_count = 0
    calls = {}
    details = []

    for marker in MARKERS:
        row = observation_for_marker(rows, marker)
        genotype = genotype_or_missing(row)
        count = protective_allele_count(row, marker)
        if count is None:
            missing_count = missing_count + 1
            calls[marker["rsid"]] = genotype
            details.append(f"{marker['rsid']} ({marker['name']}): missing.")
        else:
            total_protective = total_protective + count
            calls[marker["rsid"]] = genotype
            result = marker["present"] if count > 0 else marker["absent"]
            details.append(f"{marker['rsid']} ({marker['name']}): genotype {genotype}; {result}")

    status = classify(total_protective, missing_count)
    rows = [{
        "participant_id": participant_id,
        "pcsk9_status": status,
        "status_label": status_label(status),
        "protective_allele_count": str(total_protective) if missing_count < len(MARKERS) else "",
        "rs11591147": calls.get("rs11591147", "missing"),
        "rs67608943": calls.get("rs67608943", "missing"),
        "rs28362286": calls.get("rs28362286", "missing"),
        "interpretation": interpretation(status, total_protective, details),
        "notes": GENERAL_NOTE,
    }]
    bioscript.write_tsv(output_file, rows)
    print(status)


if __name__ == "__main__":
    main()
