# Use the official R image
FROM bioconductor/r-ver:3.21-R-4.5.1

COPY ./inst/install_packages.sh /tmp/install_packages.sh
RUN chmod +x /tmp/install_packages.sh

# Install additional R packages
RUN /tmp/install_packages.sh

COPY . /app

# Install CORTAR
RUN Rscript -e "devtools::install_local('/app')"

WORKDIR /app

# By default, drop into R
ENTRYPOINT [ "/app/inst/run_cortar.sh" ]