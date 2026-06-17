FROM rust:1-bookworm AS builder

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        bash \
        ca-certificates \
        clang \
        cmake \
        git \
        libbz2-dev \
        libclang-dev \
        libcurl4-openssl-dev \
        liblzma-dev \
        libssl-dev \
        pkg-config \
        zlib1g-dev \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /app
COPY . .

RUN test -f bioscript/rust/Cargo.toml \
    || (echo "bioscript submodule is missing. Run: git submodule update --init --recursive" >&2; exit 1)

RUN cargo build --release \
    --manifest-path bioscript/rust/Cargo.toml \
    -p bioscript-cli

RUN install -m 0755 bioscript/rust/target/release/bioscript /usr/local/bin/bioscript

FROM debian:bookworm-slim AS runtime

COPY --from=ghcr.io/astral-sh/uv:latest /uv /uvx /usr/local/bin/

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        bash \
        ca-certificates \
        libbz2-1.0 \
        libcurl4 \
        liblzma5 \
        libssl3 \
        python3 \
        zlib1g \
    && rm -rf /var/lib/apt/lists/*

COPY --from=builder /usr/local/bin/bioscript /usr/local/bin/bioscript
RUN ln -s /usr/local/bin/bioscript /usr/local/bin/bs
COPY . /app

ENV BIOSCRIPT_BIN=/usr/local/bin/bioscript
ENV UV_CACHE_DIR=/opt/uv-cache

WORKDIR /app
RUN uv run --with pyyaml python -c "import yaml" \
    && BIOSCRIPT_BIN=/usr/local/bin/bioscript FORMATTER=/usr/bin/true ./build-assay-zips.sh --self-contained \
    && mkdir -p /opt/exvitae/assays \
    && find assays -name '*.zip' -type f -exec cp {} /opt/exvitae/assays/ \;
COPY --chmod=0755 docker/exvitae-report /usr/local/bin/exvitae-report

WORKDIR /workspace

CMD ["bioscript"]
