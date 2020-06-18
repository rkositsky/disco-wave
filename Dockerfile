# Base Image
FROM r-base:3.6.0

# Metadata
LABEL base.image="disco-wave:v.20200618"
LABEL version="1"
LABEL software="disco-wave"
LABEL software.version="1.0.0"
LABEL description="A bioinformatics tool to identify translocations using discordant reads"
LABEL tags="Translocation Structural Variant"

# Maintainer
MAINTAINER DaveLab <lab.dave@gmail.com>

# update the OS related packages
RUN apt-get update -y && apt-get install -y \
    build-essential \
    libnss-sss \
    curl \
    vim \
    less \
    wget \
    unzip \
    cmake \
    python3 \
    gawk \
    python-pip \
    zlib1g-dev \
    libncurses5-dev \
    libncursesw5-dev \
    libbz2-dev \
    liblzma-dev \
    bzip2 \
    libcurl4-openssl-dev \
    libssl-dev \
    git \
    autoconf \
    bsdmainutils \
    bedtools

# install Python libraries
WORKDIR /usr/local/bin
RUN pip install argparse
RUN pip install pysam

# install R required dependencies
RUN R --vanilla -e 'install.packages(c("optparse", "tidyverse", "gridExtra", "viridis"), repos="http://cran.us.r-project.org")'

# clone disco-wave repo
ADD https://api.github.com/repos/rkositsky/disco-wave/git/refs/heads/ version.json
RUN git clone https://github.com/rkositsky/disco-wave.git

# add disco-wave repo to SYSPATH
ENV PATH /disco-wave:$PATH

# change the permission of the repo
RUN chmod 777 -R /disco-wave

# make the Main.py as default script to execute when no command is given
CMD ["python3 Main.py -h"]

