# Base image
# Use rocker/tidyverse as the base image for R and tidyverse
FROM rocker/tidyverse:latest
# Copy entrypoint.sh into the container
COPY entrypoint.sh /entrypoint.sh


# Install system dependencies (CLI tools + R/Bioconductor build requirements)
RUN apt-get update && apt-get install -y \
	wget \
	ca-certificates \
	build-essential \
	libcurl4-openssl-dev \
	libssl-dev \
	libxml2-dev \
	zlib1g-dev \
	libbz2-dev \
	liblzma-dev \
	&& rm -rf /var/lib/apt/lists/*

# Download prebuilt fastp binary (v0.23.4) from OpenGene
RUN wget -O /usr/local/bin/fastp "http://opengene.org/fastp/fastp.0.23.4" \
	&& chmod a+x /usr/local/bin/fastp

# Install read-mapping and variant-calling tools
RUN apt-get update && apt-get install -y \
	minimap2 \
	samtools \
	bcftools \
	bedtools \
	&& rm -rf /var/lib/apt/lists/*

# Install Bioconductor dependencies for SNP annotation
RUN Rscript -e "\
	if (!requireNamespace('BiocManager', quietly = TRUE)) install.packages('BiocManager'); \
	BiocManager::install(c('GenomicFeatures', 'GenomicRanges', 'AnnotationDbi', 'VariantAnnotation', 'BSgenome', 'Rsamtools', 'GenomicAlignments', 'rtracklayer'), ask = FALSE, update = FALSE) \
"

# Install PA14-specific annotation packages from GitHub
RUN Rscript -e "\
	if (!requireNamespace('remotes', quietly = TRUE)) install.packages('remotes'); \
	remotes::install_github('utubun/BSgenome.Paeruginosa.NCBI.PA14'); \
	remotes::install_github('utubun/TxDb.Paeruginosa.PA14') \
"