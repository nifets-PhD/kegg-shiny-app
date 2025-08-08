FROM rocker/shiny:latest

# Set environment variables for C++ compilation and R optimization
ENV CXX14=g++
ENV CXX14FLAGS="-O3 -march=native -mtune=native"
ENV CXX14STD="-std=c++14"

# R performance optimizations
ENV R_COMPILE_PKGS=1
ENV R_DISABLE_HTTPD=1
ENV OMP_NUM_THREADS=8
ENV OPENBLAS_NUM_THREADS=8
ENV MKL_NUM_THREADS=8

# Install system dependencies including C++14 support and performance libraries
RUN apt-get update && apt-get install -y \
    gdebi \
    build-essential \
    libcurl4-openssl-dev \
    libssl-dev \
    libfontconfig1-dev \
    libxml2-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libfreetype6-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev \
    libsodium-dev \
    libudunits2-dev \
    libgdal-dev \
    libgeos-dev \
    libproj-dev \
    g++ \
    gcc \
    gfortran \
    libblas-dev \
    liblapack-dev \
    libopenblas-dev \
    libnode-dev \
    libglpk-dev \
    libglpk40 \
    curl \
    htop \
    && rm -rf /var/lib/apt/lists/*

# Set working directory
WORKDIR /srv/shiny-server/kegg-app

# Copy dependency files first for better caching
COPY DESCRIPTION .
COPY global.R .

# Install BiocManager first
RUN R -e "if (!requireNamespace('BiocManager', quietly = TRUE)) install.packages('BiocManager')"

# Install remotes for GitHub package installation (more reliable than pak for GitHub packages)
RUN R -e "install.packages('remotes')"

# Install R packages in specific order to handle dependencies
RUN R -e "install.packages(c('shiny', 'shinydashboard', 'visNetwork', 'DT', 'dplyr', 'plotly', 'shinyWidgets', 'ggplot2', 'scales', 'xml2', 'igraph', 'tidyr', 'stringr'))"

# Install Bioconductor packages
RUN R -e "BiocManager::install(c('KEGGREST', 'KEGGgraph', 'biomaRt', 'clusterProfiler', 'org.Hs.eg.db', 'DOSE'), ask = FALSE, update = FALSE)"

# Install myTAI from GitHub using remotes (more reliable for GitHub packages)
RUN R -e "remotes::install_github('drostlab/myTAI', upgrade = 'never')"

# Copy application files
COPY . .

# Create cache directories with proper permissions
RUN mkdir -p kegg_cache cache data && \
    chown -R shiny:shiny /srv/shiny-server/kegg-app && \
    chmod -R 755 /srv/shiny-server/kegg-app

# Copy custom shiny-server configuration AFTER other files (cache bust: v3)
COPY deployment/shiny-server.conf /etc/shiny-server/shiny-server.conf

EXPOSE 3838

# Add healthcheck with longer start period for KEGG data loading
HEALTHCHECK --interval=60s --timeout=10s --start-period=120s --retries=3 \
    CMD curl -f http://localhost:3838/ || exit 1

CMD ["/usr/bin/shiny-server"]
