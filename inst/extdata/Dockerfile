FROM bioconductor/bioconductor_docker:devel

RUN apt-get update \
	&& apt-get install -y \
	apt-utils \
	libxml2-dev \
	libglpk-dev \
	libnode-dev \
	libgit2-dev

LABEL org.opencontainers.image.licenses="GPL-2.0-or-later" \
      org.opencontainers.image.source="https://github.com/rocker-org/rocker-versioned2" 

COPY install_shiny_server.sh ./install_shiny_server.sh
RUN chmod 755 ./install_shiny_server.sh
RUN ./install_shiny_server.sh

RUN R -e 'BiocManager::install(c("survminer", "survival","fastmatch","reshape2","beeswarm","grDevices","shinycssloaders","shinythemes","networkD3","httr","RColorBrewer","psych","stringr","shiny","shinydashboard","shinyWidgets","shinybusy","matrixStats","flatxml","excelR","shinyjs","shinyFiles","DT","plotly","openxlsx","yaml","curl","sortable","BiocManager","password","ggseqlogo","devtools","Matrix","RSQLite"))'
RUN R -e 'BiocManager::install(c("Biobase", "fgsea", "S4Vectors", "SummarizedExperiment"), update = FALSE)'
RUN R -e 'devtools::install_github("mengchen18/omicsViewer", dependencies = FALSE)'


RUN mkdir /srv/shiny-server/sample-apps/ExpressionSetViewer
COPY app/ui.R /srv/shiny-server/sample-apps/ExpressionSetViewer/
COPY app/server.R /srv/shiny-server/sample-apps/ExpressionSetViewer/

RUN mkdir /media/ExpressionSetViewerData
RUN chmod 775 /media/ExpressionSetViewerData

# expose port
EXPOSE 3838

# run app on container start
CMD ["R", "-e", "shiny::runApp('/srv/shiny-server/sample-apps/ExpressionSetViewer/', host = '0.0.0.0', port = 3838)"]
