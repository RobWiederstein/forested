# This image already includes Tidyverse, Quarto, and the Spatial Stack (GDAL/GEOS/PROJ)
FROM rocker/geospatial:latest

# Create a home for the project
WORKDIR /home/rstudio/project

# 1. Copy renv files first (better for Docker caching)
COPY renv.lock renv.lock
COPY .Rprofile .Rprofile
COPY renv/activate.R renv/activate.R
COPY renv/settings.json renv/settings.json

# 2. Restore the R packages (this uses the system libs already in the image)
RUN R -e "install.packages('renv'); renv::restore()"

# 3. Copy your actual code
COPY . .

# Default to starting R
CMD ["R"]
