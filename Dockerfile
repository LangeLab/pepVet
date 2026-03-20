FROM bioconductor/bioconductor_docker:RELEASE_3_22

WORKDIR /opt/pepVet

RUN Rscript -e 'install.packages(c("renv", "pak", "littler"), repos = "https://cloud.r-project.org")'

COPY . .

RUN if [ -f renv.lock ]; then Rscript -e 'renv::restore(prompt = FALSE)'; fi

CMD ["R"]
