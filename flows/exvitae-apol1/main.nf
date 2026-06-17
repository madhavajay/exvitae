// ExVitae APOL1 risk report flow.

nextflow.enable.dsl=2

def EXVITAE_CONTAINER = 'ghcr.io/openmined/exvitae:0.2.3'

def countryValue(record) {
    def facets = record.facets ?: [:]
    return (
        record.country ?:
        record.Country ?:
        record.country_code ?:
        record.countryCode ?:
        facets.country ?:
        facets.Country ?:
        facets.country_code ?:
        facets.countryCode ?:
        ''
    ).toString().trim()
}

def shellQuote(value) {
    return "'" + value.toString().replace("'", "'\"'\"'") + "'"
}

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
        def assayPackage = params.assay_package ?: '/opt/exvitae/assays/APOL1.zip'
        def assayPackageCh = Channel.value(assayPackage)
        def participantItems = participants.map { record ->
            tuple(record.participant_id.toString(), countryValue(record), file(record.genotype_file))
        }
        def perParticipantReports = exvitae_apol1_report(assayPackageCh, participantItems)
        def report_outputs = aggregate_reports(perParticipantReports.report_dir.collect())

    emit:
        participant_reports = report_outputs.participant_reports
        country_aggregates = report_outputs.country_aggregates
        observations = report_outputs.observations
        reports = report_outputs.reports
        analysis = report_outputs.analysis
}

process exvitae_apol1_report {
    container EXVITAE_CONTAINER
    stageInMode 'copy'
    tag { participant_id }
    errorStrategy { params.nextflow.error_strategy }
    maxRetries { params.nextflow.max_retries }
    maxForks REPORT_MAX_FORKS

    input:
        val assay_package
        tuple val(participant_id), val(country), path(input_file)

    output:
        path "${participant_id}", emit: report_dir

    script:
    def inputFileName = input_file.name
    """
    bs report "${assay_package}" \
      --root "\${PWD}" \
      --input-file "\${PWD}/${inputFileName}" \
      --output-dir "${participant_id}" \
      --detect-sex \
      --html \
      --analysis-max-duration-ms 30000
    { printf 'participant_id\\tcountry\\n'; printf '%s\\t%s\\n' ${shellQuote(participant_id)} ${shellQuote(country)}; } \
      > "${participant_id}/metadata.tsv"
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
        path "countries", emit: country_aggregates

    script:
    """
    set -euo pipefail

    python3 - <<'PY'
import re
import shutil
from pathlib import Path

OBS_HEADER = "participant_id\\tassay_id\\tassay_version\\tvariant_key\\trsid\\tassembly\\tchrom\\tpos_start\\tpos_end\\tref\\talt\\tkind\\tmatch_status\\tcoverage_status\\tcall_status\\tgenotype\\tgenotype_display\\tzygosity\\tref_count\\talt_count\\tdepth\\tgenotype_quality\\tallele_balance\\toutcome\\tevidence_type\\tevidence_raw\\tfacets\\n"

def safe_id(value):
    safe = re.sub(r"[^A-Za-z0-9_.-]+", "_", value.strip())
    return safe.strip("._") or "unknown"

def read_country(report_dir):
    path = report_dir / "metadata.tsv"
    if not path.is_file():
        return ""
    for line in path.read_text(encoding="utf-8").splitlines()[1:2]:
        parts = line.split("\\t")
        return parts[1].strip() if len(parts) > 1 else ""
    return ""

report_dirs = sorted(
    path for path in Path(".").iterdir()
    if path.is_dir() and path.name not in {"participants", "countries"}
)

Path("participants").mkdir(exist_ok=True)
Path("countries").mkdir(exist_ok=True)
reports = Path("reports.jsonl").open("w", encoding="utf-8")
analysis = Path("analysis.jsonl").open("w", encoding="utf-8")
obs = Path("observations.tsv").open("w", encoding="utf-8")
obs.write(OBS_HEADER)
country_obs_started = {}

for report_dir in report_dirs:
    name = report_dir.name
    target = Path("participants") / name
    if target.exists():
        shutil.rmtree(target)
    shutil.copytree(report_dir, target)

    obs_path = report_dir / "observations.tsv"
    if obs_path.is_file():
        lines = obs_path.read_text(encoding="utf-8").splitlines()
        for line in lines[1:]:
            if line:
                obs.write(line + "\\n")

    for source, handle in ((report_dir / "reports.jsonl", reports), (report_dir / "analysis.jsonl", analysis)):
        if source.is_file():
            text = source.read_text(encoding="utf-8")
            if text:
                handle.write(text)
                if not text.endswith("\\n"):
                    handle.write("\\n")

    country = read_country(report_dir)
    if not country:
        continue
    country_root = Path("countries") / safe_id(country)
    country_participants = country_root / "participants"
    country_participants.mkdir(parents=True, exist_ok=True)
    country_target = country_participants / name
    if country_target.exists():
        shutil.rmtree(country_target)
    shutil.copytree(report_dir, country_target)

    country_obs = country_root / "observations.tsv"
    if not country_obs_started.get(country_root):
        country_obs.write_text(OBS_HEADER, encoding="utf-8")
        country_obs_started[country_root] = True
    if obs_path.is_file():
        with country_obs.open("a", encoding="utf-8") as handle:
            for line in obs_path.read_text(encoding="utf-8").splitlines()[1:]:
                if line:
                    handle.write(line + "\\n")
    for filename in ("reports.jsonl", "analysis.jsonl"):
        source = report_dir / filename
        if source.is_file():
            text = source.read_text(encoding="utf-8")
            if text:
                with (country_root / filename).open("a", encoding="utf-8") as handle:
                    handle.write(text)
                    if not text.endswith("\\n"):
                        handle.write("\\n")

reports.close()
analysis.close()
obs.close()
PY
    """
}
