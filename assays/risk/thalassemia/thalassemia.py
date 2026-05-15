ALPHA_CLASSIFICATION_SOURCE = "https://www.ncbi.nlm.nih.gov/books/NBK1435/"
TIF_PREVENTION_SOURCE = "https://www.thalassaemia.org.cy/wp-content/uploads/2019/11/Preventiondiagnosis-of-Hbpathies_BOOKLETNEW-1.pdf"
NOT_INFERRED = "not_inferred_from_snp_array_catalogue"


def row_value(row, key):
    return row[key] if key in row else ""


def split_list(value):
    if value == "":
        return []
    out = []
    for item in value.split("|"):
        item = item.strip()
        if item != "":
            out.append(item)
    return out


def int_or_zero(value):
    if value == "":
        return 0
    return int(value)


def append_unique(values, value):
    if value == "":
        return
    if value not in values:
        values.append(value)


def gene_bucket(gene):
    if gene == "HBA1" or gene == "HBA2":
        return "alpha"
    if gene == "HBB":
        return "beta"
    if gene == "HBD" or gene == "HBG1" or gene == "HBG2" or gene == "HBE1" or gene == "HBZ":
        return "other_globin"
    return "modifier_or_other"


def phenotype_contains(finding, text):
    return text in row_value(finding, "phenotype").lower()


def is_causative(finding):
    return row_value(finding, "functionality").lower() == "causative"


def thalassemia_family(finding, bucket):
    phenotype = row_value(finding, "phenotype").lower()
    if "hpfh" in phenotype:
        return "HPFH"
    if "δβ" in phenotype or "delta-beta" in phenotype:
        return "delta_beta_thalassemia"
    if "εγδβ" in phenotype:
        return "epsilon_gamma_delta_beta_thalassemia"
    if "δ-chain" in phenotype or "delta-chain" in phenotype:
        return "delta_thalassemia_or_delta_chain_variant"
    if "α-thalassaemia" in phenotype or "alpha-thalassaemia" in phenotype:
        return "alpha_thalassemia"
    if "β-thalassaemia" in phenotype or "beta-thalassaemia" in phenotype:
        return "beta_thalassemia"
    if bucket == "alpha" and "α-chain variant" in phenotype:
        return "structural_hemoglobin_variant"
    if bucket == "beta" and "β-chain variant" in phenotype:
        return "structural_hemoglobin_variant"
    return ""


def allele_functional_class(finding, bucket):
    if not is_causative(finding):
        return row_value(finding, "functionality")
    family = thalassemia_family(finding, bucket)
    if family == "alpha_thalassemia":
        return "alpha_thalassemia_allele_unclassified_alpha_plus_alpha0_or_non_deletional"
    if family == "beta_thalassemia":
        return "beta_thalassemia_allele_unclassified_beta0_beta_plus_or_dominant"
    if family == "HPFH":
        return "HPFH_allele"
    if family != "":
        return family + "_allele"
    return "causative_allele_unclassified"


def index_variants(rows):
    by_rsid = {}
    for row in rows:
        rsid = row_value(row, "rsid")
        if rsid != "":
            by_rsid[rsid] = row
        aliases = split_list(row_value(row, "aliases"))
        for alias in aliases:
            by_rsid[alias] = row
    return by_rsid


def index_findings(rows):
    by_variant = {}
    for row in rows:
        variant_id = row_value(row, "variant_id")
        if variant_id not in by_variant:
            by_variant[variant_id] = []
        by_variant[variant_id].append(row)
    return by_variant


def genotype_allele_count(genotype, allele):
    if genotype == "" or allele == "":
        return 0
    count = 0
    for ch in genotype:
        if ch == allele:
            count = count + 1
    return count


def allele_observed(observation, finding):
    alt = row_value(finding, "alt")
    if alt == "*":
        return int_or_zero(row_value(observation, "alt_count")) > 0
    genotype = row_value(observation, "genotype")
    if len(alt) == 1 and genotype != "":
        return genotype_allele_count(genotype, alt) > 0
    return int_or_zero(row_value(observation, "alt_count")) > 0


def classify_counts(alpha_count, beta_count, other_count):
    if alpha_count == 0 and beta_count == 0 and other_count == 0:
        return "no_catalogue_thalassemia_variant_observed"
    labels = []
    if beta_count == 1:
        labels.append("possible_beta_thalassemia_carrier")
    elif beta_count >= 2:
        labels.append("possible_biallelic_HBB_thalassemia_risk")
    if alpha_count == 1:
        labels.append("possible_alpha_thalassemia_silent_carrier_or_trait")
    elif alpha_count == 2:
        labels.append("possible_alpha_thalassemia_trait_or_compound_alpha_globin_carrier")
    elif alpha_count >= 3:
        labels.append("possible_HbH_or_Hb_Bart_alpha_thalassemia_spectrum")
    if other_count > 0:
        labels.append("other_globin_chain_or_modifier_variant_observed")
    return ";".join(labels)


def alpha_state(alpha_count):
    if alpha_count == 0:
        return "no_alpha_thalassemia_causative_variant_observed"
    if alpha_count == 1:
        return "possible_alpha_thalassemia_silent_carrier_or_trait"
    if alpha_count == 2:
        return "possible_alpha_thalassemia_trait_or_compound_alpha_globin_carrier"
    return "possible_HbH_or_Hb_Bart_alpha_thalassemia_spectrum"


def beta_state(beta_count):
    if beta_count == 0:
        return "no_beta_thalassemia_causative_variant_observed"
    if beta_count == 1:
        return "possible_beta_thalassemia_carrier_or_trait"
    return "possible_biallelic_HBB_thalassemia_risk"


def finding_detail(rsid, variant, finding):
    phenotype = row_value(finding, "phenotype")
    if phenotype == "":
        phenotype = "N/A"
    return (
        rsid
        + ":"
        + row_value(variant, "gene")
        + ":"
        + row_value(finding, "alt")
        + ":"
        + row_value(finding, "functionality")
        + ":"
        + phenotype
    )


def classify_thalassemia(observations, variants, findings):
    variants_by_rsid = index_variants(variants)
    findings_by_variant = index_findings(findings)
    alpha_count = 0
    beta_count = 0
    other_count = 0
    observed = []
    observed_alpha = []
    observed_beta = []
    observed_other = []
    disorder_types = []
    functional_classes = []
    matched_rows = []

    for observation in observations:
        rsid = row_value(observation, "matched_rsid")
        if rsid == "":
            rsid = row_value(observation, "rsid")
        if rsid == "":
            continue
        if rsid not in variants_by_rsid:
            continue
        variant = variants_by_rsid[rsid]
        variant_id = row_value(variant, "variant_id")
        gene = row_value(variant, "gene")
        bucket = gene_bucket(row_value(variant, "gene"))
        variant_findings = findings_by_variant[variant_id] if variant_id in findings_by_variant else []
        for finding in variant_findings:
            if not allele_observed(observation, finding):
                continue
            family = thalassemia_family(finding, bucket)
            functional_class = allele_functional_class(finding, bucket)
            detail = finding_detail(rsid, variant, finding)
            if is_causative(finding) and family == "alpha_thalassemia":
                alpha_count = alpha_count + 1
            elif is_causative(finding) and family == "beta_thalassemia":
                beta_count = beta_count + 1
            elif bucket != "alpha" and bucket != "beta":
                other_count = other_count + 1
            append_unique(disorder_types, family)
            append_unique(functional_classes, functional_class)
            observed.append(detail)
            if bucket == "alpha":
                observed_alpha.append(detail)
            elif bucket == "beta":
                observed_beta.append(detail)
            else:
                observed_other.append(detail)
            matched_rows.append(
                {
                    "participant_id": bioscript.context["participant_id"],
                    "rsid": rsid,
                    "variant_id": variant_id,
                    "gene": gene,
                    "bucket": bucket,
                    "alt": row_value(finding, "alt"),
                    "functionality": row_value(finding, "functionality"),
                    "phenotype": row_value(finding, "phenotype"),
                    "globin_disorder_type": family,
                    "allele_functional_class": functional_class,
                    "genotype": row_value(observation, "genotype"),
                    "alt_count": row_value(observation, "alt_count"),
                }
            )

    status = classify_counts(alpha_count, beta_count, other_count)
    notes = (
        "Screening summary only; not diagnostic. Alpha-thalassemia severity depends on the number of inactive "
        "alpha-globin alleles; GeneReviews describes HbH disease as usually three inactive alpha-globin alleles "
        "and Hb Bart syndrome as four inactive alpha-globin alleles. This assay does not determine phase, deletion "
        "structure, Hb fractions, or clinical severity. Sources: "
        + ALPHA_CLASSIFICATION_SOURCE
        + " ; "
        + TIF_PREVENTION_SOURCE
    )
    return {
        "thalassemia_status": status,
        "globin_disorder_types": "|".join(disorder_types),
        "allele_functional_classes": "|".join(functional_classes),
        "alpha_diplotype": NOT_INFERRED,
        "beta_diplotype": NOT_INFERRED,
        "alpha_carrier_or_disease_state": alpha_state(alpha_count),
        "beta_carrier_or_disease_state": beta_state(beta_count),
        "clinical_severity_phenotype": "not_inferred_requires_cbc_hba2_hbf_and_clinical_context",
        "transfusion_dependence_phenotype": "not_inferred_requires_clinical_course",
        "special_named_syndromes": NOT_INFERRED,
        "alpha_variant_findings": str(alpha_count),
        "beta_variant_findings": str(beta_count),
        "other_globin_or_modifier_findings": str(other_count),
        "matched_alpha_findings": "; ".join(observed_alpha),
        "matched_beta_findings": "; ".join(observed_beta),
        "matched_other_findings": "; ".join(observed_other),
        "thalassemia_findings": "; ".join(observed),
        "report_notes": notes,
        "matched_rows": matched_rows,
    }


observations = bioscript.read_tsv(bioscript.context["observations_file"])
variants = bioscript.read_tsv(bioscript.context["assets"]["variants"])
findings = bioscript.read_tsv(bioscript.context["assets"]["findings"])
result = classify_thalassemia(observations, variants, findings)
bioscript.write_tsv("/output/thalassemia_matched_findings.tsv", result["matched_rows"])

bioscript.write_tsv(bioscript.context["output_file"], [
    {
        "participant_id": bioscript.context["participant_id"],
        "thalassemia_status": result["thalassemia_status"],
        "globin_disorder_types": result["globin_disorder_types"],
        "allele_functional_classes": result["allele_functional_classes"],
        "alpha_diplotype": result["alpha_diplotype"],
        "beta_diplotype": result["beta_diplotype"],
        "alpha_carrier_or_disease_state": result["alpha_carrier_or_disease_state"],
        "beta_carrier_or_disease_state": result["beta_carrier_or_disease_state"],
        "clinical_severity_phenotype": result["clinical_severity_phenotype"],
        "transfusion_dependence_phenotype": result["transfusion_dependence_phenotype"],
        "special_named_syndromes": result["special_named_syndromes"],
        "alpha_variant_findings": result["alpha_variant_findings"],
        "beta_variant_findings": result["beta_variant_findings"],
        "other_globin_or_modifier_findings": result["other_globin_or_modifier_findings"],
        "matched_alpha_findings": result["matched_alpha_findings"],
        "matched_beta_findings": result["matched_beta_findings"],
        "matched_other_findings": result["matched_other_findings"],
        "thalassemia_findings": result["thalassemia_findings"],
        "report_notes": result["report_notes"],
    }
])
