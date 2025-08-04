FROM rocker/shiny:latest

# Install system dependencies
RUN apt-get update && apt-get install -y \
    libxml2-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    libgit2-dev \
    && rm -rf /var/lib/apt/lists/*

# Install R packages
RUN R -e "install.packages(c('shinydashboard', 'visNetwork', 'DT', 'xml2', 'dplyr', 'readr'), repos='https://cran.r-project.org/')"

# Install Bioconductor packages
RUN R -e "if (!require('BiocManager', quietly = TRUE)) install.packages('BiocManager'); BiocManager::install(c('KEGGREST', 'KEGGgraph', 'biomaRt'))"

# Install devtools for GitHub installations
RUN R -e "install.packages('devtools', repos='https://cran.r-project.org/')"

# Install myTAI from GitHub
RUN R -e "devtools::install_github('drostlab/myTAI')"

# Copy app files
COPY . /srv/shiny-server/kegg-app/

# Create cache directory with proper permissions
RUN mkdir -p /srv/shiny-server/kegg-app/cache && \
    mkdir -p /srv/shiny-server/kegg-app/kegg_cache && \
    mkdir -p /srv/shiny-server/kegg-app/hsa_cache && \
    chown -R shiny:shiny /srv/shiny-server/kegg-app/

# Copy pre-built cache if available
COPY kegg_cache/ /srv/shiny-server/kegg-app/kegg_cache/
COPY hsa_cache/ /srv/shiny-server/kegg-app/hsa_cache/
COPY data/ /srv/shiny-server/kegg-app/data/

# Set proper permissions
RUN chown -R shiny:shiny /srv/shiny-server/

# Expose port
EXPOSE 3838

# Run Shiny Server
CMD ["/usr/bin/shiny-server.sh"]
