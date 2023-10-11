FROM rocker/tidyverse:4.2.0

# APT packages
RUN apt-get update \
&& apt-get install -y --no-install-recommends \
  fftw3 \
  libgsl-dev \
  libmagick++-dev \
  libglpk-dev \
  libgdal-dev \
  libproj-dev \
  libudunits2-dev \
  gawk \
  libxt6 \
  libzmq3-dev

# Other software
# taxon-tools
ENV APPS_HOME=/apps
RUN mkdir $APPS_HOME
WORKDIR $APPS_HOME

ENV APP_NAME=taxon-tools
RUN git clone https://github.com/camwebb/$APP_NAME.git && \
cd $APP_NAME && \
git checkout 509cb5241d97d13cf254c01f202dcb3581c3cd53 && \
make check && \
make install

# taxastand
RUN Rscript -e 'remotes::install_github("joelnitta/taxastand")'

# ggtreeExtra
RUN Rscript -e 'remotes::install_github("YuLab-SMU/ggtreeExtra")'

# R packages
RUN install2.r -e \
  rgdal \
  sf \
  raster \
  diversitree \
  ape \
  readr \
  phytools \
  exactextractr \
  foreach \
  doParallel \
  here \
  targets \
  clustermq \
  ggstar \
  ggimage
