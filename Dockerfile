# =============================================================================
# Stage 1: Builder — install R packages using PPM binary repos
# =============================================================================
FROM rocker/r-ver:4.3.3 AS builder

# Build-time system deps (headers + libs needed to install R packages)
RUN apt-get update && apt-get install -y --no-install-recommends \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libfontconfig1-dev \
    libfreetype6-dev \
    libpng-dev \
    libtiff-dev \
    libjpeg-turbo8-dev \
    zlib1g-dev \
    && rm -rf /var/lib/apt/lists/*

# Install CRAN packages (rocker/r-ver defaults to PPM binary repo)
# then Bioconductor annotation packages — all in one layer
RUN R -q -e " \
  install.packages(c( \
    'plumber', 'ggplot2', 'ggrepel', 'dplyr', \
    'readr', 'jsonlite', 'BiocManager' \
  )); \
  BiocManager::install( \
    c('AnnotationDbi', 'org.Mm.eg.db'), \
    ask = FALSE, update = FALSE \
  )"

# Strip docs, help pages, vignettes, tests, and compiled object files
# to shrink the library that gets copied to the runtime stage
RUN find /usr/local/lib/R/site-library -type d \
      \( -name "doc" -o -name "html" -o -name "help" -o -name "demo" \
         -o -name "examples" -o -name "tests" -o -name "unitTests" \
         -o -name "testdata" -o -name "tinytest" -o -name "vignettes" \) \
      -exec rm -rf {} + 2>/dev/null; \
    find /usr/local/lib/R/site-library -type f \
      \( -name "*.o" -o -name "*.a" -o -name "*.so.dSYM" \) \
      -delete 2>/dev/null; \
    exit 0

# =============================================================================
# Stage 2: Runtime — minimal image with only what's needed to run
# =============================================================================
FROM rocker/r-ver:4.3.3

LABEL maintainer="cfeng" \
      description="Volcano plot API (plumber + ggplot2)" \
      org.opencontainers.image.source="https://github.com/cfeng204/r1-volcano-plotting"

# Runtime-only system libraries (NO -dev headers)
RUN apt-get update && apt-get install -y --no-install-recommends \
    libcurl4 \
    libxml2 \
    libfontconfig1 \
    libfreetype6 \
    libpng16-16 \
    libtiff5 \
    libjpeg-turbo8 \
    curl \
    ca-certificates \
    && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

# Copy pre-built R packages from builder
COPY --from=builder /usr/local/lib/R/site-library /usr/local/lib/R/site-library

WORKDIR /app

COPY api/ /app/api/
COPY R/   /app/R/
COPY web/ /app/web/

RUN mkdir -p /data
COPY data/example.csv /data/example.csv

ENV RES_FILE=/data/example.csv \
    PORT=8000

EXPOSE 8000

# Run as non-root
RUN groupadd -r appuser && useradd -r -g appuser -d /app appuser \
    && chown -R appuser:appuser /app /data
USER appuser

CMD ["Rscript", "api/run_api.R"]