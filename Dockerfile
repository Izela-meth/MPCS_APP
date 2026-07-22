FROM rocker/shiny:4.3.3

RUN apt-get update && apt-get install -y \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    git \
    && rm -rf /var/lib/apt/lists/*

RUN R -e "install.packages(c('shiny', 'bslib', 'DT', 'dplyr', 'tidyr', 'readxl', 'haven', 'ggplot2', 'igraph', 'RColorBrewer', 'patchwork'), repos='https://cloud.r-project.org/')"

WORKDIR /srv/shiny-server/

COPY app.R .
COPY functions/ /srv/shiny-server/functions/
COPY data/ /srv/shiny-server/data/

EXPOSE 7860

CMD ["R", "-e", "shiny::runApp('/srv/shiny-server/', host='0.0.0.0', port=7860)"]
