def main():
    genotypes = bioscript.load_genotypes(input_file)
    value = genotypes.get("rs73885319")
    rows = [{
        "participant_id": participant_id,
        "rs73885319": value,
    }]
    bioscript.write_tsv(output_file, rows)
    print(value)


if __name__ == "__main__":
    main()
