# Fast builds with binary CRAN packages
FROM rocker/r2u:jammy

RUN apt-get update && apt-get install -y --no-install-recommends \
    libcurl4-openssl-dev libxml2-dev libssl-dev \
    libfontconfig1-dev libfreetype6-dev libpng-dev libtiff-dev libjpeg-dev \
    ca-certificates && \
    rm -rf /var/lib/apt/lists/*

# CRAN + Bioconductor
RUN R -q -e "install.packages(c('plumber','ggplot2','ggrepel','dplyr','readr','AnnotationDbi','jsonlite'), repos='https://cloud.r-project.org')" \
 && R -q -e "install.packages('BiocManager', repos='https://cloud.r-project.org')" \
 && R -q -e "BiocManager::install('org.Mm.eg.db', ask=FALSE, update=FALSE)"

WORKDIR /app

# Copy only the directories the app needs (keeps image small)
COPY api/  /app/api/
COPY R/    /app/R/
COPY web/  /app/web/

# Ship a default dataset inside the image
RUN mkdir -p /data
COPY data/example.csv /data/example.csv

# Prefer the baked-in CSV by default
ENV RES_FILE=/data/example.csv
ENV PORT=8000
EXPOSE 8000

CMD ["Rscript", "api/run_api.R"]