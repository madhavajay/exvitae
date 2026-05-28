// ExVitae APOL1 risk report flow.

nextflow.enable.dsl=2

def envInt = { name, fallback ->
    def raw = System.getenv(name)
    if (raw && raw.isInteger()) {
        return raw.toInteger()
    }
    return fallback
}

def DEFAULT_REPORT_MAX_FORKS = envInt('EXVITAE_APOL1_MAX_FORKS', 10)
def REPORT_MAX_FORKS = params.nextflow?.report_max_forks ?: (params.nextflow?.max_forks ?: DEFAULT_REPORT_MAX_FORKS)

workflow USER {
    take:
        context
        participants

    main:
        def assetsDir = context.assets_dir
        if (!assetsDir) {
            throw new IllegalStateException("Missing assets directory in context")
        }

        def assayPackage = file("${assetsDir}/APOL1.zip")
        def assayPackageCh = Channel.value(assayPackage)
        def participantItems = participants.map { record ->
            tuple(record.participant_id.toString(), file(record.genotype_file))
        }
        def perParticipantReports = exvitae_apol1_report(assayPackageCh, participantItems)
        def report_outputs = aggregate_reports(perParticipantReports.report_dir.collect())

    emit:
        participant_reports = report_outputs.participant_reports
        observations = report_outputs.observations
        reports = report_outputs.reports
        analysis = report_outputs.analysis
}

process exvitae_apol1_report {
    container 'ghcr.io/openmined/bioscript:0.2.0'
    stageInMode 'copy'
    tag { participant_id }
    errorStrategy { params.nextflow.error_strategy }
    maxRetries { params.nextflow.max_retries }
    maxForks REPORT_MAX_FORKS

    input:
        path assay_package
        tuple val(participant_id), path(input_file)

    output:
        path "${participant_id}", emit: report_dir

    script:
    def assayPackageName = assay_package.name
    def inputFileName = input_file.name
    """
    bs report "\${PWD}/${assayPackageName}" \
      --root "\${PWD}" \
      --input-file "\${PWD}/${inputFileName}" \
      --output-dir "${participant_id}" \
      --detect-sex \
      --html \
      --analysis-max-duration-ms 30000
    """
}

process aggregate_reports {
    container 'debian:bookworm-slim'
    publishDir params.results_dir, mode: 'copy', overwrite: true
    stageInMode 'copy'

    input:
        path report_dirs

    output:
        path "participants", emit: participant_reports
        path "observations.tsv", emit: observations
        path "reports.jsonl", emit: reports
        path "analysis.jsonl", emit: analysis

    script:
    """
    set -euo pipefail

    mkdir -p participants
    : > reports.jsonl
    : > analysis.jsonl
    first_observations=1

    for report_dir in \$(find . -maxdepth 1 -mindepth 1 -type d ! -name participants | sort); do
      name="\$(basename "\${report_dir}")"
      cp -R "\${report_dir}" "participants/\${name}"

      if [ -f "\${report_dir}/observations.tsv" ]; then
        if [ "\${first_observations}" = "1" ]; then
          cat "\${report_dir}/observations.tsv" > observations.tsv
          first_observations=0
        else
          tail -n +2 "\${report_dir}/observations.tsv" >> observations.tsv
        fi
      fi

      if [ -f "\${report_dir}/reports.jsonl" ]; then
        cat "\${report_dir}/reports.jsonl" >> reports.jsonl
      fi

      if [ -f "\${report_dir}/analysis.jsonl" ]; then
        cat "\${report_dir}/analysis.jsonl" >> analysis.jsonl
      fi
    done

    if [ ! -f observations.tsv ]; then
      printf "participant_id\\tassay_id\\tassay_version\\tvariant_key\\trsid\\tassembly\\tchrom\\tpos_start\\tpos_end\\tref\\talt\\tkind\\tmatch_status\\tcoverage_status\\tcall_status\\tgenotype\\tgenotype_display\\tzygosity\\tref_count\\talt_count\\tdepth\\tgenotype_quality\\tallele_balance\\toutcome\\tevidence_type\\tevidence_raw\\tfacets\\n" > observations.tsv
    fi
    """
}
