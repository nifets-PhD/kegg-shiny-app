FROM rocker/shiny:latest

ENV CXX14=g++ \
    CXX14FLAGS="-O3 -march=native -mtune=native" \
    CXX14STD="-std=c++14" \
    R_COMPILE_PKGS=1 \
    R_DISABLE_HTTPD=1 \
    OMP_NUM_THREADS=8 \
    OPENBLAS_NUM_THREADS=8 \
    MKL_NUM_THREADS=8

RUN apt-get update && apt-get install -y --no-install-recommends \
    gdebi \
    build-essential \
    libcurl4-openssl-dev libssl-dev libfontconfig1-dev libxml2-dev \
    libharfbuzz-dev libfribidi-dev libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev \
    libsodium-dev libudunits2-dev libgdal-dev libgeos-dev libproj-dev \
    g++ gcc gfortran libblas-dev liblapack-dev libopenblas-dev \
    libnode-dev libglpk-dev libglpk40 curl htop \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /srv/shiny-server/kegg-app

COPY DESCRIPTION .
COPY global.R .

RUN R -e "install.packages('remotes')" && \
    R -e "install.packages(c('shiny','shinydashboard','visNetwork','DT','dplyr','plotly','shinyWidgets','ggplot2','scales','xml2','igraph','tidyr','stringr'))" && \
    R -e "if (!requireNamespace('BiocManager', quietly=TRUE)) install.packages('BiocManager')" && \
    R -e "BiocManager::install(c('KEGGREST','KEGGgraph','biomaRt','clusterProfiler','org.Hs.eg.db','DOSE'), ask=FALSE, update=FALSE)" && \
    R -e "remotes::install_github('drostlab/myTAI', upgrade='never')"

COPY . .
RUN mkdir -p kegg_cache cache data && \
    chown -R shiny:shiny /srv/shiny-server/kegg-app && \
    chmod -R 755 /srv/shiny-server/kegg-app

COPY deployment/shiny-server.conf /etc/shiny-server/shiny-server.conf

EXPOSE 3838
HEALTHCHECK --interval=60s --timeout=10s --start-period=120s --retries=3 \
    CMD curl -f http://localhost:3838/ || exit 1

CMD ["/usr/bin/shiny-server"]
