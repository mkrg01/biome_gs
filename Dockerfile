FROM rocker/tidyverse:4.3.2

# install by apt-get
RUN apt-get update && apt-get install -y \
	locales \
	fftw3 \
	libgsl-dev \
	libmagick++-dev \
	libglpk-dev \
	libgdal-dev \
	libproj-dev \
	libudunits2-dev \
	gawk \
	libxt6 \
	libzmq3-dev \
 && rm -rf /var/lib/apt/lists/*

# Generate C.UTF-8 locale
RUN locale-gen C.UTF-8

# Set environment variables
ENV LC_ALL C.UTF-8
ENV LANG C.UTF-8
ENV LANGUAGE C.UTF-8

# Install from CRAN
RUN install2.r --deps TRUE --repos 'https://cloud.r-project.org' \
	BiocManager \
	remotes

# Install by BioManager
RUN Rscript -e "BiocManager::install('treeio')"
RUN Rscript -e "BiocManager::install('ggtree')"
RUN Rscript -e "BiocManager::install('ggtreeExtra')"

# Install by remotes
RUN R -e "remotes::install_version('sf', version='1.0.16', repos='https://cloud.r-project.org')"
RUN R -e "remotes::install_version('raster', version='3.6.26', repos='https://cloud.r-project.org')"
RUN R -e "remotes::install_version('diversitree', version='0.9.16', repos='https://cloud.r-project.org')"
RUN R -e "remotes::install_version('ape', version='5.8', repos='https://cloud.r-project.org')"
RUN R -e "remotes::install_version('readr', version='2.1.2', repos='https://cloud.r-project.org')"
RUN R -e "remotes::install_version('phytools', version='1.0.3', repos='https://cloud.r-project.org')"
RUN R -e "remotes::install_version('exactextractr', version='0.8.2', repos='https://cloud.r-project.org')"
RUN R -e "remotes::install_version('foreach', version='1.5.2', repos='https://cloud.r-project.org')"
RUN R -e "remotes::install_version('doParallel', version='1.0.17', repos='https://cloud.r-project.org')"
RUN R -e "remotes::install_version('here', version='1.0.1', repos='https://cloud.r-project.org')"
RUN R -e "remotes::install_version('targets', version='0.12.0', repos='https://cloud.r-project.org')"
RUN R -e "remotes::install_version('clustermq', version='0.8.95.3', repos='https://cloud.r-project.org')"
RUN R -e "remotes::install_version('ggstar', version='1.0.3', repos='https://cloud.r-project.org')"
RUN R -e "remotes::install_version('mvtnorm', version='1.2.5', repos='https://cloud.r-project.org')"
