#!/bin/bash
set -e

## System + R runtime setup for the omicsViewer application image.
##
## shiny-server is intentionally NOT installed: the Dockerfile runs the
## app directly via `R -e shiny::runApp(...)` (see its CMD), so the
## shiny-server .deb, its s6 service scripts and the rstudio example apps
## were dead weight - and their download source retired: Posit purged the
## download3.rstudio.org/ubuntu-14.04/x86_64 shiny-server binaries, so the
## `wget .../shiny-server-${VERSION}-amd64.deb` started returning 404
## (wget exit code 8) and broke the CI image build. Only the pieces the
## image actually uses are installed here:
##   - system libraries needed to build/run the app's R packages
##     (libcurl for httr/curl, Cairo/X11 for plotly and base graphics);
##   - shiny + rmarkdown (everything else comes from the
##     BiocManager::install lines in the Dockerfile).

## build ARGs
NCPUS=${NCPUS:--1}

# System libraries required to build and run the app's R packages
apt-get update
apt-get install -y --no-install-recommends \
    libcurl4-openssl-dev \
    libcairo2-dev \
    libxt-dev

# Core runtime packages
install2.r --error --skipinstalled -n "$NCPUS" shiny rmarkdown

# Clean up
rm -rf /var/lib/apt/lists/*
rm -rf /tmp/downloaded_packages
