FROM rocker/tidyverse:4.4.0

RUN apt-get update && apt-get install -y \
    nano \
    neovim \
    git \
    bash-completion \
    openssh-client \
    cmake \
    libglpk-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libfontconfig1-dev \
    libfreetype6-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libpng-dev \
    libjpeg-dev \
    libtiff-dev \
    libwebp-dev \
    gdal-bin \
    libgdal-dev \
    libproj-dev \
    libgeos-dev \
    libudunits2-dev \
    && rm -rf /var/lib/apt/lists/*

# --- CONFIGURATION FIXES ---
# 1. Move library out of project folder (fixes Volume Trap)
ENV RENV_PATHS_LIBRARY=/renv/library
# 2. Disable Symlinks (fixes Root Permission Trap) <--- CRITICAL NEW LINE
ENV RENV_CONFIG_CACHE_SYMLINKS=FALSE

RUN mkdir -p /renv/library && chmod 777 /renv/library

WORKDIR /home/rstudio/project

COPY renv.lock renv.lock
COPY .Rprofile .Rprofile
COPY renv/activate.R renv/activate.R
COPY renv/settings.json renv/settings.json

# Restore (Binaries + No Symlinks)
RUN R -e "options(repos = c(CRAN = 'https://packagemanager.posit.co/cran/__linux__/jammy/latest')); install.packages('renv'); renv::restore()"

COPY . .

# Ensure the actual files are readable by everyone
RUN chmod -R 777 /renv

CMD ["R"]