def count_char(text, needle):
    if text is None:
        return 0
    total = 0
    for ch in text:
        if ch == needle:
            total = total + 1
    return total


def count_non_ref(text, ref):
    if text is None:
        return 0
    total = 0
    for ch in text:
        if ch != ref and ch != "-":
            total = total + 1
    return total


def classify_apol1(site1, site2, g2):
    if site1 is None and site2 is None and g2 is None:
        return "G-/G-"

    d_count = count_char(g2, "D")
    site1_variants = count_non_ref(site1, "A")
    site2_variants = count_non_ref(site2, "T")

    has_g1 = site1_variants > 0 and site2_variants > 0
    if has_g1:
        g1_total = site1_variants + site2_variants
    else:
        g1_total = 0

    if d_count == 2:
        return "G2/G2"
    if d_count == 1:
        if g1_total >= 2:
            return "G2/G1"
        return "G2/G0"
    if g1_total == 4:
        return "G1/G1"
    if g1_total >= 2:
        return "G1/G0"
    return "G0/G0"


def main():
    status = classify_apol1("AA", None, "II")
    rows = [{
        "participant_id": participant_id,
        "apol1_status": status,
    }]
    bioscript.write_tsv(output_file, rows)
    print(status)


if __name__ == "__main__":
    main()
