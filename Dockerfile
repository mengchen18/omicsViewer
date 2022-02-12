FROM rocker/shiny:4.1.2

RUN apt-get update \
	&& apt-get install -y \
	apt-utils \
	libxml2-dev \
	libglpk-dev \
	libnode-dev \
	libgit2-dev

RUN R -e 'install.packages("survminer")'
RUN R -e 'install.packages("survival")'
RUN R -e 'install.packages("fastmatch")'
RUN R -e 'install.packages("reshape2")'
RUN R -e 'install.packages("beeswarm")'
RUN R -e 'install.packages("grDevices")'
RUN R -e 'install.packages("shinycssloaders")'
RUN R -e 'install.packages("shinythemes")'
RUN R -e 'install.packages("networkD3")'
RUN R -e 'install.packages("httr")'
RUN R -e 'install.packages("RColorBrewer")'
RUN R -e 'install.packages("psych")'
RUN R -e 'install.packages("stringr")'
RUN R -e 'install.packages("shiny")'
RUN R -e 'install.packages("shinydashboard")'
RUN R -e 'install.packages("shinyWidgets")'
RUN R -e 'install.packages("shinybusy")'
RUN R -e 'install.packages("matrixStats")'
RUN R -e 'install.packages("flatxml")'
RUN R -e 'install.packages("excelR")'
RUN R -e 'install.packages("shinyjs")'
RUN R -e 'install.packages("shinyFiles")'
RUN R -e 'install.packages("DT")'
RUN R -e 'install.packages("plotly")'
RUN R -e 'install.packages("openxlsx")'
RUN R -e 'install.packages("yaml")'
RUN R -e 'install.packages("curl",)'
RUN R -e 'install.packages("sortable")'
RUN R -e 'install.packages("BiocManager")'
RUN R -e 'install.packages("password")'
RUN R -e 'install.packages("ggseqlogo")'
RUN R -e 'install.packages("devtools")'
RUN R -e 'install.packages("Matrix")'
RUN R -e 'install.packages("RSQLite")'
RUN R -e 'BiocManager::install("Biobase", update = FALSE)'
RUN R -e 'BiocManager::install("fgsea", update = FALSE)'
RUN R -e 'BiocManager::install("S4Vectors", update = FALSE)'
RUN R -e 'BiocManager::install("SummarizedExperiment", update = FALSE)'
RUN R -e 'devtools::install_github("mengchen18/omicsViewer", dependencies = FALSE)'


RUN mkdir /srv/shiny-server/sample-apps/ExpressionSetViewer
COPY inst/app/ui.R /srv/shiny-server/sample-apps/ExpressionSetViewer/
COPY inst/app/server.R /srv/shiny-server/sample-apps/ExpressionSetViewer/

RUN mkdir /media/ExpressionSetViewerData
RUN chmod 775 /media/ExpressionSetViewerData

# expose port
EXPOSE 3838

# run app on container start
CMD ["R", "-e", "shiny::runApp('/srv/shiny-server/sample-apps/ExpressionSetViewer/', host = '0.0.0.0', port = 3838)"]


