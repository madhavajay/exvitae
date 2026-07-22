HLA_PROXY_MARKERS = [
    ("HLA-A*3101", "rs1061235", ""),
    ("HLA-B*1502", "rs3909184", ""),
    ("HLA-B*1502", "rs2844682", ""),
    ("HLA-B*5701", "rs2395029", ""),
    ("HLA-B*5801", "rs9263726", "allopurinol-induced SCAR"),
]
HLA_PROXY_REFERENCE_ALLELES = {
    "rs1061235": "A",
    "rs3909184": "G",
    "rs2844682": "G",
    "rs2395029": "T",
}


def observation_for_rsid(rows, rsid):
    for row in rows:
        matched = row.get("matched_rsid") or row.get("rsid") or ""
        path = row.get("path") or ""
        name = row.get("name") or row.get("variant_key") or ""
        if matched == rsid or rsid in path or rsid in name:
            return row
    return None


def raw_genotype(row):
    if row is None:
        return "missing"
    genotype = row.get("genotype") or ""
    display = row.get("genotype_display") or ""
    if genotype == "" and display != "":
        genotype = display
    if genotype == "" or genotype == "./." or genotype == ".|.":
        return "missing"
    return genotype


def normalize_genotype(genotype):
    if genotype is None:
        return "missing"

    text = genotype.strip()
    if text == "" or text.lower() in ("missing", "unknown", "no_call", "no-call", "nocall"):
        return "missing"

    alleles = []
    for ch in text:
        if ch != "/" and ch != "|" and ch != " " and ch != "-":
            alleles.append(ch)

    if len(alleles) == 0:
        return "missing"

    alleles.sort()
    return "".join(alleles)


def display_genotype(genotype):
    if genotype == "missing":
        return "missing"
    if len(genotype) == 2:
        return genotype[0] + "/" + genotype[1]
    return genotype


def display_rs9263726_genotype(genotype):
    if genotype == "missing":
        return "missing"
    if genotype == "GG":
        return "G/G"
    if genotype == "AG":
        return "G/A"
    if genotype == "AA":
        return "A/A"
    if len(genotype) == 2:
        return genotype[0] + "/" + genotype[1]
    return genotype


def display_marker_genotype(rsid, genotype):
    if rsid == "rs9263726":
        return display_rs9263726_genotype(genotype)
    return display_genotype(genotype)


def rs9263726_proxy_interpretation(genotype):
    if genotype == "GG":
        return "proxy-negative"
    if genotype == "AG":
        return "proxy-positive / likely HLA-B*58:01 carrier"
    if genotype == "AA":
        return "proxy-positive / likely HLA-B*58:01 carrier, possibly homozygous tag"
    return "unknown"


def rs9263726_proxy_status(genotype):
    interpretation = rs9263726_proxy_interpretation(genotype)
    if interpretation == "proxy-negative":
        return "negative"
    if interpretation.startswith("proxy-positive"):
        return "positive"
    return "unknown"


def marker_variant_status(rsid, genotype):
    if genotype == "missing":
        return "missing"
    ref = HLA_PROXY_REFERENCE_ALLELES.get(rsid)
    if ref is None or len(genotype) != 2:
        return "unknown"
    if genotype == ref + ref:
        return "negative"
    return "positive"


def marker_status(rsid, genotype):
    if rsid == "rs9263726":
        return rs9263726_proxy_status(genotype)
    return marker_variant_status(rsid, genotype)


def marker_phenotype(hla_allele, rsid, genotype):
    if rsid == "rs9263726":
        status = rs9263726_proxy_status(genotype)
        if status == "positive":
            return "Positive for allopurinol-induced SCAR proxy"
        if status == "negative":
            return "Negative for allopurinol-induced SCAR proxy"
        return "Unknown allopurinol-induced SCAR proxy"
    status = marker_variant_status(rsid, genotype)
    if status == "positive":
        return "Positive"
    if status == "negative":
        return "Negative"
    return "Unknown"


def marker_interpretation(hla_allele, rsid, genotype):
    if rsid == "rs9263726":
        return rs9263726_proxy_interpretation(genotype)
    return (
        hla_allele + " proxy marker " + rsid + " genotype is " + display_marker_genotype(rsid, genotype) + ". "
        "This row lists the marker genotype only; no aggregate HLA inference is emitted."
    )


def marker_notes(rsid):
    if rsid == "rs9263726":
        return (
            "Research-use HLA-B*58:01 proxy only. G/G is proxy-negative; G/A and A/A are proxy-positive; "
            "C or other alleles are unknown. "
            "This is not definitive HLA typing or clinical medication advice."
        )
    return "Research-use HLA proxy marker listing only; this is not definitive HLA typing."


def main():
    rows = bioscript.read_tsv(bioscript.context["observations_file"])
    calls = []
    for marker in HLA_PROXY_MARKERS:
        rsid = marker[1]
        calls.append(normalize_genotype(raw_genotype(observation_for_rsid(rows, rsid))))

    rows = []
    for index in range(len(HLA_PROXY_MARKERS)):
        marker = HLA_PROXY_MARKERS[index]
        hla_allele = marker[0]
        rsid = marker[1]
        condition = marker[2]
        genotype = calls[index]
        rows.append({
            "participant_id": participant_id,
            "hla_proxy_status": marker_status(rsid, genotype),
            "hla_allele": hla_allele,
            "rsid": rsid,
            "genotype": display_marker_genotype(rsid, genotype),
            "drug": "allopurinol" if rsid == "rs9263726" else "",
            "condition": condition,
            "phenotype": marker_phenotype(hla_allele, rsid, genotype),
            "interpretation": marker_interpretation(hla_allele, rsid, genotype),
            "notes": marker_notes(rsid),
        })

    bioscript.write_tsv(output_file, rows)
    print("hla_proxy_markers")


if __name__ == "__main__":
    main()
