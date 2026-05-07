MTHFR_SOURCE = "https://pmc.ncbi.nlm.nih.gov/articles/PMC6743281/"

RS1801133 = bioscript.variant(
    rsid=["rs1801133", "rs4134713", "rs59514310", "rs386545618"],
    grch37="1:11856378-11856378",
    grch38="1:11796321-11796321",
    ref="G",
    alt="A",
    kind="snp",
)

RS1801131 = bioscript.variant(
    rsid=["rs1801131", "rs4134712", "rs17367365", "rs17857426"],
    grch37="1:11854476-11854476",
    grch38="1:11794419-11794419",
    ref="T",
    alt="G",
    kind="snp",
)

MTHFR_QUERY_PLAN = bioscript.query_plan([
    RS1801133,
    RS1801131,
])


MTHFR_COMBINATIONS = {
    (0, 0): {
        "status": "mthfr_677cc_1298aa_reference",
        "display": "677C>T (CC) 1298A>C (AA)",
        "estimated_activity_pct": "100",
        "activity_note": "reference/common activity",
        "source": MTHFR_SOURCE,
    },
    (1, 0): {
        "status": "mthfr_677ct_1298aa_677_heterozygous",
        "display": "677C>T (CT) 1298A>C (AA)",
        "estimated_activity_pct": "65",
        "activity_note": "C677T heterozygous; the cited review states CT genotype has 65% of normal enzyme activity.",
        "source": MTHFR_SOURCE,
    },
    (2, 0): {
        "status": "mthfr_677tt_1298aa_677_homozygous",
        "display": "677C>T (TT) 1298A>C (AA)",
        "estimated_activity_pct": "<=30",
        "activity_note": "C677T homozygous; the cited review states 677TT has no more than 30% of normal enzyme activity.",
        "source": MTHFR_SOURCE,
    },
    (0, 1): {
        "status": "mthfr_677cc_1298ac_1298_heterozygous",
        "display": "677C>T (CC) 1298A>C (AC)",
        "estimated_activity_pct": "~85",
        "activity_note": "A1298C heterozygous; the cited review states AC carriers have a 15% reduction of enzymatic activity.",
        "source": MTHFR_SOURCE,
    },
    (0, 2): {
        "status": "mthfr_677cc_1298cc_1298_homozygous",
        "display": "677C>T (CC) 1298A>C (CC)",
        "estimated_activity_pct": "~70",
        "activity_note": "A1298C homozygous; the cited review states CC carriers have a 30% reduction of enzymatic activity.",
        "source": MTHFR_SOURCE,
    },
    (1, 1): {
        "status": "mthfr_677ct_1298ac_compound_heterozygous",
        "display": "677C>T (CT) 1298A>C (AC)",
        "estimated_activity_pct": "unknown",
        "activity_note": "compound heterozygous C677T/A1298C; this assay does not assign a percent activity estimate because the cited source does not provide a specific combined estimate.",
        "source": MTHFR_SOURCE,
    },
    (1, 2): {
        "status": "mthfr_677ct_1298cc_rare_combined",
        "display": "677C>T (CT) 1298A>C (CC)",
        "estimated_activity_pct": "unknown",
        "activity_note": "rare combined genotype; this assay does not assign a percent activity estimate because the cited source does not provide a specific combined estimate.",
        "source": MTHFR_SOURCE,
    },
    (2, 1): {
        "status": "mthfr_677tt_1298ac_rare_combined",
        "display": "677C>T (TT) 1298A>C (AC)",
        "estimated_activity_pct": "unknown",
        "activity_note": "rare combined genotype; this assay does not assign a percent activity estimate because the cited source does not provide a specific combined estimate.",
        "source": MTHFR_SOURCE,
    },
    (2, 2): {
        "status": "mthfr_677tt_1298cc_very_rare_combined",
        "display": "677C>T (TT) 1298A>C (CC)",
        "estimated_activity_pct": "unknown",
        "activity_note": "very rare combined genotype; this assay does not assign a percent activity estimate because the cited source does not provide a specific combined estimate.",
        "source": MTHFR_SOURCE,
    },
}


def allele_dosage(genotype, allele):
    if genotype is None:
        return None

    count = 0
    seen = 0

    for ch in genotype:
        if ch != "/" and ch != "|" and ch != " " and ch != "-":
            seen = seen + 1
            if ch == allele:
                count = count + 1

    if seen == 0:
        return None

    return count


def classify_mthfr(rs1801133_gt, rs1801131_gt):
    # rs1801133 alt=A corresponds to legacy 677T.
    c677t_dosage = allele_dosage(rs1801133_gt, "A")

    # rs1801131 alt=G corresponds to legacy 1298C.
    a1298c_dosage = allele_dosage(rs1801131_gt, "G")

    if c677t_dosage is None or a1298c_dosage is None:
        return {
            "status": "mthfr_unresolved_missing_variant",
            "display": "MTHFR unresolved",
            "estimated_activity_pct": "unknown",
            "activity_note": "missing rs1801133 or rs1801131 genotype",
            "source": MTHFR_SOURCE,
            "c677t_dosage": c677t_dosage,
            "a1298c_dosage": a1298c_dosage,
        }

    result = MTHFR_COMBINATIONS[(c677t_dosage, a1298c_dosage)].copy()
    result["c677t_dosage"] = c677t_dosage
    result["a1298c_dosage"] = a1298c_dosage

    return result


def mthfr_outcome(result):
    if result["status"] == "mthfr_677cc_1298aa_reference":
        return "normal"
    if result["status"] == "mthfr_unresolved_missing_variant":
        return "unknown"
    return "variant"


def report_notes(result):
    if result["estimated_activity_pct"] == "unknown":
        return (
            'Individuals with "' + result["status"] + '" have the combined genotype ' + result["display"] + ". "
            "This assay does not assign a percent activity estimate for this combination from the cited source. "
            "Consult a licensed doctor for advice."
        )

    return (
        'Individuals with "' + result["status"] + '" have the combined genotype ' + result["display"]
        + " and should expect an estimated " + result["estimated_activity_pct"]
        + "% MTHFR enzymatic activity based on the cited source. "
        "Consult a licensed doctor for advice."
    )


def main():
    genotypes = bioscript.load_genotypes(input_file)

    rs1801133_gt, rs1801131_gt = genotypes.lookup_variants(
        MTHFR_QUERY_PLAN
    )

    result = classify_mthfr(rs1801133_gt, rs1801131_gt)

    rows = [{
        "participant_id": participant_id,
        "mthfr_outcome": mthfr_outcome(result),
        "mthfr_display": result["display"],
        "rs1801131": rs1801131_gt,
        "rs1801133": rs1801133_gt,
        "mthfr_estimated_activity_pct": result["estimated_activity_pct"],
        "mthfr_activity_note": result["activity_note"],
        "mthfr_status": result["status"],
        "mthfr_c677t_variant_dosage": result["c677t_dosage"],
        "mthfr_a1298c_variant_dosage": result["a1298c_dosage"],
        "mthfr_source": result["source"],
        "notes": report_notes(result),
    }]

    bioscript.write_tsv(output_file, rows)

    print(result["status"])


if __name__ == "__main__":
    main()
