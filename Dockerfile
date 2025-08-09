# Use the official R image
FROM bioconductor/r-ver:3.21-R-4.5.1

# Install additional R packages
COPY ./inst/install_packages.sh /tmp/install_packages.sh
RUN chmod +x /tmp/install_packages.sh && /tmp/install_packages.sh

WORKDIR /tmp

# Install micromamba (arch-specific)
RUN ARCH=$(uname -m) && \
    if [ "$ARCH" = "x86_64" ]; then \
        MAMBA_URL="https://micro.mamba.pm/api/micromamba/linux-64/latest"; \
    elif [ "$ARCH" = "aarch64" ]; then \
        MAMBA_URL="https://micro.mamba.pm/api/micromamba/linux-aarch64/latest"; \
    else \
        echo "Unsupported architecture: $ARCH" && exit 1; \
    fi && \
    curl -L "$MAMBA_URL" | tar -xvj bin/micromamba && \
    mv bin/micromamba /usr/local/bin/micromamba && \
    rm -rf bin

# Create an env in a fixed path and install samtools
RUN /usr/local/bin/micromamba create -y -p /opt/conda/envs/bioenv -c conda-forge -c bioconda samtools=1.22.1 && \
    /usr/local/bin/micromamba clean --all -y

# Put that env first on PATH so binaries are available without activation
ENV PATH=/opt/conda/envs/bioenv/bin:$PATH

WORKDIR /app
COPY . .

# Install CORTAR
RUN Rscript -e "devtools::install_local( \
      '.', \
      dependencies = FALSE, \
      upgrade      = 'never', \
      args         = c('--no-multiarch','--no-test-load') \
    )"

# By default, drop into R
ENTRYPOINT [ "/app/inst/run_cortar.sh" ]