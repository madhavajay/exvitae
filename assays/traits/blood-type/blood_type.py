MARKERS = [
    "rs8176719",
    "rs8176746",
    "rs590787",
    "rs8176058",
    "rs12075",
    "rs2814778",
    "rs1058396",
    "rs601338",
]

GENERAL_NOTE = (
    "Genotype-derived exploratory prediction only; this is not clinical blood typing. "
    "Confirm clinically relevant results with serology or validated molecular blood-group testing."
)


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
    genotype = row.get("genotype_display") or row.get("genotype") or ""
    if genotype == "" or genotype == "./." or genotype == ".|.":
        return "missing"
    return genotype


def alt_count(row):
    if row is None:
        return None
    value = row.get("alt_count")
    if value is None or value == "":
        return None
    try:
        return int(value)
    except ValueError:
        return None


def normalize_genotype(genotype):
    if genotype is None:
        return ""
    text = genotype.upper().strip()
    if text == "" or text == "MISSING":
        return ""

    text = text.replace(" ", "")

    if text in ("ID", "DI"):
        return "DI"
    if text in ("DD", "II"):
        return text

    if "/" in text or "|" in text:
        parts = text.replace("|", "/").split("/")
        compact = ""
        for part in parts:
            if part not in ("", "."):
                compact = compact + part
        text = compact

    if len(text) == 2:
        return "".join(sorted(text))
    return text


def normalize_snv_genotype(genotype):
    text = normalize_genotype(genotype)
    bases = ""
    for ch in text:
        if ch in "ACGT":
            bases = bases + ch
    if len(bases) == 2:
        return "".join(sorted(bases))
    return bases


def insertion_state(row, genotype):
    text = normalize_genotype(genotype)
    count = alt_count(row)

    if text in ("DD", "DI", "II"):
        return text
    if count == 0:
        return "DD"
    if count == 1:
        return "DI"
    if count == 2:
        return "II"

    if "/" in genotype or "|" in genotype:
        parts = genotype.upper().replace("|", "/").split("/")
        insertion_count = 0
        usable_parts = 0
        for part in parts:
            if part in ("", "."):
                continue
            usable_parts = usable_parts + 1
            if part == "TC" or part.endswith("C") and len(part) > 1:
                insertion_count = insertion_count + 1
        if usable_parts == 2:
            if insertion_count == 0:
                return "DD"
            if insertion_count == 1:
                return "DI"
            if insertion_count == 2:
                return "II"

    return ""


def make_result(system, phenotype, genotype_interpretation, confidence, genotypes, limitations, interpretation):
    return {
        "participant_id": participant_id,
        "system": system,
        "phenotype": phenotype,
        "genotype_interpretation": genotype_interpretation,
        "confidence": confidence,
        "genotypes": genotypes,
        "interpretation": interpretation,
        "limitations": limitations,
        "notes": GENERAL_NOTE,
    }


def genotype_summary(calls, rsids):
    parts = []
    for rsid in rsids:
        parts.append(rsid + "=" + calls.get(rsid, "missing"))
    return "; ".join(parts)


def predict_abo(rows, calls):
    row_o = observation_for_rsid(rows, "rs8176719")
    row_b = observation_for_rsid(rows, "rs8176746")
    o_marker = insertion_state(row_o, calls.get("rs8176719", ""))
    b_marker = normalize_snv_genotype(calls.get("rs8176746", ""))
    b_count = b_marker.count("T")
    genotypes = genotype_summary(calls, ["rs8176719", "rs8176746"])
    limitations = (
        "Common-marker ABO prediction only. Requires verified strand/REF/ALT and does not assess rs8176747, phase, "
        "rare ABO alleles, cis/trans ambiguity, Bombay phenotype, or other causes of ABO discrepancy."
    )

    if o_marker == "" or b_marker == "":
        return make_result(
            "ABO",
            "Unresolved",
            "Missing rs8176719 or rs8176746",
            "low",
            genotypes,
            limitations,
            "ABO could not be resolved because one or both required marker calls were missing or unsupported.",
        )

    if o_marker == "DD":
        phenotype = "O"
        interpretation = "O/O"
        reason = "Two common O-deletion marker alleles were observed."
    elif o_marker == "DI" and b_count == 0:
        phenotype = "A"
        interpretation = "A/O"
        reason = "One common O marker allele and one non-B functional ABO allele were observed."
    elif o_marker == "DI" and b_count >= 1:
        phenotype = "B"
        interpretation = "B/O"
        reason = "One common O marker allele and one B-associated marker allele were observed."
    elif o_marker == "II" and b_count == 0:
        phenotype = "A"
        interpretation = "A/A"
        reason = "No common O marker allele and no B-associated marker allele were observed."
    elif o_marker == "II" and b_count == 1:
        phenotype = "AB"
        interpretation = "A/B"
        reason = "One A-like and one B-associated marker allele were observed."
    elif o_marker == "II" and b_count == 2:
        phenotype = "B"
        interpretation = "B/B"
        reason = "Two B-associated marker alleles and no common O marker allele were observed."
    else:
        return make_result(
            "ABO",
            "Unresolved",
            "rs8176719=" + o_marker + "; rs8176746=" + b_marker,
            "low",
            genotypes,
            limitations,
            "The marker combination does not fit the simplified common ABO rule set.",
        )

    return make_result("ABO", phenotype, interpretation, "moderate", genotypes, limitations, reason)


def predict_rhd(calls):
    genotype = normalize_snv_genotype(calls.get("rs590787", ""))
    genotypes = genotype_summary(calls, ["rs590787"])
    limitations = (
        "rs590787 is only a population- and platform-dependent RhD proxy, most informative in European-ancestry contexts. "
        "It does not detect RHD deletion structure, RHD/RHCE hybrids, weak D, partial D, DEL alleles, or mapping artifacts."
    )

    if genotype == "":
        return make_result(
            "RhD",
            "Unresolved",
            "rs590787 missing",
            "low",
            genotypes,
            limitations,
            "No usable rs590787 proxy genotype was available.",
        )
    if "A" in genotype:
        return make_result(
            "RhD",
            "Likely RhD positive proxy",
            "rs590787 contains A",
            "low",
            genotypes,
            limitations,
            "The available proxy marker is consistent with likely RHD presence, but this is not definitive RhD typing.",
        )
    return make_result(
        "RhD",
        "Possibly RhD negative proxy",
        "rs590787 lacks A",
        "low",
        genotypes,
        limitations,
        "The available proxy marker may be consistent with RHD absence in some cohorts, but this is not definitive RhD typing.",
    )


def predict_kell(calls):
    genotype = normalize_snv_genotype(calls.get("rs8176058", ""))
    genotypes = genotype_summary(calls, ["rs8176058"])
    limitations = "Predicts only the common K/k pair. Does not assess rare KEL null, weak, or variant alleles."

    if genotype == "GG":
        return make_result("Kell", "K-k+", "KEL*02/KEL*02", "moderate", genotypes, limitations, "Common k/k marker genotype.")
    if genotype == "AG":
        return make_result("Kell", "K+k+", "KEL*01/KEL*02", "moderate", genotypes, limitations, "One K-associated and one k-associated marker allele.")
    if genotype == "AA":
        return make_result("Kell", "K+k-", "KEL*01/KEL*01", "moderate", genotypes, limitations, "Two K-associated marker alleles.")
    return make_result("Kell", "Unresolved", "Unsupported rs8176058=" + genotype, "low", genotypes, limitations, "Kell K/k could not be resolved from this marker call.")


def predict_duffy(calls):
    antigen = normalize_snv_genotype(calls.get("rs12075", ""))
    promoter = normalize_snv_genotype(calls.get("rs2814778", ""))
    genotypes = genotype_summary(calls, ["rs12075", "rs2814778"])
    limitations = (
        "Prediction concerns red-cell Duffy expression. ACKR1 may remain expressed in non-erythroid tissues. "
        "Weak FY*X, rare null alleles, and phase are not fully assessed."
    )

    if promoter == "CC":
        if antigen == "AA":
            interpretation = "Likely FY*02N.01/FY*02N.01"
        elif antigen == "AG":
            interpretation = "Erythroid-silent promoter genotype; coding alleles require phasing"
        elif antigen == "GG":
            interpretation = "Erythroid-silent promoter genotype with unusual rs12075 background"
        else:
            interpretation = "Homozygous erythroid-silent promoter genotype"
        return make_result(
            "Duffy",
            "Fy(a-b-), erythroid Duffy-null",
            interpretation,
            "high",
            genotypes,
            limitations,
            "rs2814778 CC is consistent with homozygous erythroid silencing of ACKR1 on red cells.",
        )

    if antigen == "GG":
        return make_result("Duffy", "Fy(a+b-)", "FY*A/FY*A", "moderate", genotypes, limitations, "Two common Fy(a) marker alleles.")
    if antigen == "AG":
        return make_result("Duffy", "Fy(a+b+)", "FY*A/FY*B", "moderate", genotypes, limitations, "One common Fy(a) and one Fy(b) marker allele.")
    if antigen == "AA":
        return make_result("Duffy", "Fy(a-b+)", "FY*B/FY*B", "moderate", genotypes, limitations, "Two common Fy(b) marker alleles without homozygous erythroid-silent promoter call.")
    return make_result("Duffy", "Unresolved", "Unsupported rs12075=" + antigen + "; rs2814778=" + promoter, "low", genotypes, limitations, "Duffy could not be resolved from these marker calls.")


def predict_kidd(calls):
    genotype = normalize_snv_genotype(calls.get("rs1058396", ""))
    genotypes = genotype_summary(calls, ["rs1058396"])
    limitations = "Predicts only the common JK*01/JK*02 Jk(a)/Jk(b) marker. Rare weak and JK-null alleles are not assessed."

    if genotype == "GG":
        return make_result("Kidd", "Jk(a+b-)", "JK*01/JK*01", "moderate", genotypes, limitations, "Two common Jk(a) marker alleles.")
    if genotype == "AG":
        return make_result("Kidd", "Jk(a+b+)", "JK*01/JK*02", "moderate", genotypes, limitations, "One Jk(a) and one Jk(b) marker allele.")
    if genotype == "AA":
        return make_result("Kidd", "Jk(a-b+)", "JK*02/JK*02", "moderate", genotypes, limitations, "Two common Jk(b) marker alleles.")
    return make_result("Kidd", "Unresolved", "Unsupported rs1058396=" + genotype, "low", genotypes, limitations, "Kidd could not be resolved from this marker call.")


def predict_fut2(calls):
    genotype = normalize_snv_genotype(calls.get("rs601338", ""))
    genotypes = genotype_summary(calls, ["rs601338"])
    limitations = (
        "rs601338 performs best where FUT2 W143X is the main non-secretor allele, especially many European-ancestry "
        "cohorts. Other ancestry-specific FUT2 loss-of-function alleles are not assessed."
    )

    if genotype == "GG":
        return make_result("FUT2 secretor", "Secretor", "FUT2 functional/functional", "moderate", genotypes, limitations, "Two common functional FUT2 alleles.")
    if genotype == "AG":
        return make_result("FUT2 secretor", "Secretor", "FUT2 functional/common loss-of-function", "moderate", genotypes, limitations, "One functional FUT2 copy is generally sufficient for secretor status.")
    if genotype == "AA":
        return make_result("FUT2 secretor", "Non-secretor", "Common FUT2 loss-of-function/loss-of-function", "moderate", genotypes, limitations, "Two copies of the common FUT2 W143X loss-of-function marker allele.")
    return make_result("FUT2 secretor", "Unresolved", "Unsupported rs601338=" + genotype, "low", genotypes, limitations, "FUT2 secretor status could not be resolved from this marker call.")


def main():
    observations_path = bioscript.context["observations_file"]
    rows = bioscript.read_tsv(observations_path)

    calls = {}
    for rsid in MARKERS:
        calls[rsid] = raw_genotype(observation_for_rsid(rows, rsid))

    output_rows = [
        predict_abo(rows, calls),
        predict_rhd(calls),
        predict_kell(calls),
        predict_duffy(calls),
        predict_kidd(calls),
        predict_fut2(calls),
    ]
    bioscript.write_tsv(output_file, output_rows)
    print("blood_type_prediction")


if __name__ == "__main__":
    main()
