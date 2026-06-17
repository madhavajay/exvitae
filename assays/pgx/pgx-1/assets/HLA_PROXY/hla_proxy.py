RS1061235 = bioscript.variant(
    rsid=["rs1061235", "rs3173409", "rs3823348", "rs41559417", "rs114782388", "rs117115314"],
    grch37="6:29913298-29913298",
    grch38="6:29945521-29945521",
    ref="A",
    alt="T",
    kind="snp",
)

RS3909184 = bioscript.variant(
    rsid=["rs3909184", "rs115153408", "rs118099270"],
    grch37="6:30699384-30699384",
    grch38="6:30731607-30731607",
    ref="G",
    alt=["C", "T"],
    kind="snp",
)

RS2844682 = bioscript.variant(
    rsid=["rs2844682", "rs60077366", "rs115766652", "rs117711722"],
    grch37="6:30946148-30946148",
    grch38="6:30978371-30978371",
    ref="G",
    alt=["A", "C", "T"],
    kind="snp",
)

RS2395029 = bioscript.variant(
    rsid=["rs2395029", "rs3997925", "rs60378661", "rs111645003", "rs114783691"],
    grch37="6:31431780-31431780",
    grch38="6:31464003-31464003",
    ref="T",
    alt="G",
    kind="snp",
)

RS9263726 = bioscript.variant(
    rsid=["rs9263726", "rs52792269", "rs61603795", "rs112940667", "rs114541727", "rs117673022", "rs186562789"],
    grch37="6:31106499-31106499",
    grch38="6:31138722-31138722",
    ref="G",
    alt=["A", "C"],
    kind="snp",
)

HLA_PROXY_MARKERS = [
    ("rs1061235", RS1061235),
    ("rs3909184", RS3909184),
    ("rs2844682", RS2844682),
    ("rs2395029", RS2395029),
    ("rs9263726", RS9263726),
]

HLA_PROXY_QUERY_PLAN = bioscript.query_plan([marker[1] for marker in HLA_PROXY_MARKERS])


def normalize_genotype(genotype):
    if genotype is None:
        return "missing"

    alleles = []
    for ch in genotype:
        if ch != "/" and ch != "|" and ch != " " and ch != "-":
            alleles.append(ch)

    if len(alleles) == 0:
        return "missing"

    alleles.sort()
    return "".join(alleles)


def hla_proxy_status(calls):
    observed = 0
    for call in calls:
        if call != "missing":
            observed = observed + 1

    if observed == len(calls):
        return "complete"
    if observed == 0:
        return "missing"
    return "partial"


def combined_key(marker_ids, calls):
    parts = []
    for index in range(len(marker_ids)):
        parts.append(marker_ids[index])
        parts.append(calls[index])
    return "_".join(parts)


def report_notes(status, key):
    return (
        "HLA proxy marker key " + key + " emitted with status " + status + ". "
        "This is a deterministic combined SNP genotype key for grouping and cohort statistics; it does not infer phase or "
        "definitively classify HLA-A or HLA-B allele carriage."
    )


def main():
    genotypes = bioscript.load_genotypes(input_file)
    raw_calls = genotypes.lookup_variants(HLA_PROXY_QUERY_PLAN)

    marker_ids = [marker[0] for marker in HLA_PROXY_MARKERS]
    calls = []
    for raw_call in raw_calls:
        calls.append(normalize_genotype(raw_call))

    status = hla_proxy_status(calls)
    observed_count = 0
    for call in calls:
        if call != "missing":
            observed_count = observed_count + 1

    key = combined_key(marker_ids, calls)

    rows = [{
        "participant_id": participant_id,
        "hla_proxy_status": status,
        "hla_proxy_key": key,
        "hla_proxy_observed_marker_count": observed_count,
        "rs1061235": calls[0],
        "rs3909184": calls[1],
        "rs2844682": calls[2],
        "rs2395029": calls[3],
        "rs9263726": calls[4],
        "notes": report_notes(status, key),
    }]

    bioscript.write_tsv(output_file, rows)
    print(key)


if __name__ == "__main__":
    main()
