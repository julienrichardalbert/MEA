##########################################################################################
# Dockerfile
#
# Software:         ALEA2
# Software Version: alea.2.0
# Description:      A toolbox for Allele-specific Epigenomics Analysis
# Base Image:       java:7
# Build Cmd:        docker build -t taskkoike/alea:2.0 /path/to/dockerfile-directory/
# Note bene:        export ALEA2='docker run -v <path/to/host/alea/data>:/alea-data -i taskkoike/alea:2.0'
# Run Cmd:          $ALEA2 setting
##########################################################################################

FROM java:7

MAINTAINER Tasuku KOIKE <task_koike_2.718281828459@hotmail.co.jp>


##########################################################################################
# Preparation
##########################################################################################

RUN apt-get update && \
  apt-get install -y \
  bc \
  wget \
  gzip \
  bzip2 \
  g++ \
  make \
  pkg-config \
  ncurses-dev \
  zlib1g-dev \
  build-essential \
  gcc-multilib \
  apt-utils \
  vim \
  python \
  python-dev \
  python-setuptools \
  libgomp1


##########################################################################################
# samtools 0.1.19
##########################################################################################

WORKDIR /tmp
RUN wget -q https://github.com/samtools/samtools/archive/0.1.19.tar.gz && \
  tar -xf 0.1.19.tar.gz -C /opt

WORKDIR /opt/samtools-0.1.19
RUN make all

RUN rm -f /tmp/0.1.19.tar.gz


##########################################################################################
# samtools 1.3.1 & tabix & gbzip
##########################################################################################

ENV SAMTOOLS_INSTALL_DIR=/opt/samtools-1.3.1
WORKDIR /tmp
RUN wget -q https://github.com/samtools/samtools/releases/download/1.3.1/samtools-1.3.1.tar.bz2 && \
  tar -xjf samtools-1.3.1.tar.bz2

WORKDIR /tmp/samtools-1.3.1
RUN ./configure --enable-plugins --prefix=$SAMTOOLS_INSTALL_DIR && \
  make all all-htslib && \
  make install install-htslib

RUN rm -f /tmp/samtools-1.3.1.tar.bz2 && \
  find /tmp/samtools-1.3.1 -type f -exec rm -f '{}' \; 
#RUN find /tmp/samtools-1.3.1/test/mpileup -type f -exec rm -vf '{}' \;
#RUN rm -rf /tmp/samtools-1.3.1


##########################################################################################
# bwa 0.7.12
##########################################################################################

WORKDIR /tmp
RUN wget -q https://downloads.sourceforge.net/project/bio-bwa/bwa-0.7.12.tar.bz2 && \
  tar -xjf bwa-0.7.12.tar.bz2 -C /opt

WORKDIR /opt/bwa-0.7.12
# Patch Makefile and Compile
RUN sed -i 's/CFLAGS=\\t\\t-g -Wall -Wno-unused-function -O2/CFLAGS=-g -Wall -Wno-unused-function -O2 -static/' Makefile && \
  make

RUN rm -f /tmp/bwa-0.7.12.tar.bz2


##########################################################################################
# bowtie 1.1.2
##########################################################################################

WORKDIR /tmp
RUN wget -q https://sourceforge.net/projects/bowtie-bio/files/bowtie/1.1.2/bowtie-1.1.2-linux-x86_64.zip && \
  unzip bowtie-1.1.2-linux-x86_64.zip -d /opt && \
  rm -f /tmp/bowtie-1.1.2-linux-x86_64.zip


##########################################################################################
# bowtie2-2.2.9
##########################################################################################

WORKDIR /tmp
RUN wget -q https://github.com/BenLangmead/bowtie2/releases/download/v2.2.9/bowtie2-2.2.9-linux-x86_64.zip && \
  unzip bowtie2-2.2.9-linux-x86_64.zip -d /opt && \
  rm -f /tmp/bowtie2-2.2.9-linux-x86_64.zip


##########################################################################################
# bedGraphToBigWig
##########################################################################################

ENV BEDGRAPHTOBIGWIG_URL=http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/bedGraphToBigWig

WORKDIR /tmp
RUN wget -q ${BEDGRAPHTOBIGWIG_URL} && \
  chmod +x bedGraphToBigWig

RUN mkdir /opt/bedGraphToBigWig && \
  cp -p /tmp/bedGraphToBigWig /opt/bedGraphToBigWig/ && \
  rm -f /tmp/bedGraphToBigWig


##########################################################################################
# shapeit.v2
##########################################################################################

WORKDIR /tmp
RUN wget -q https://mathgen.stats.ox.ac.uk/genetics_software/shapeit/shapeit.v2.r837.GLIBCv2.12.Linux.static.tgz && \
  mkdir /opt/shapeit && \
  tar -xf shapeit.v2.r837.GLIBCv2.12.Linux.static.tgz -C /opt/shapeit && \
  rm -f /tmp/shapeit.v2.r837.GLIBCv2.12.Linux.static.tgz


##########################################################################################
# VCFtools 0.1.14
##########################################################################################

ENV VCFTOOLS_INSTALL_DIR=/opt/vcftools-0.1.14
WORKDIR /tmp
RUN wget -q https://github.com/vcftools/vcftools/releases/download/v0.1.14/vcftools-0.1.14.tar.gz && \
  tar -xf vcftools-0.1.14.tar.gz -C /opt

ENV PERL5LIB=$VCFTOOLS_INSTALL_DIR/src/perl
WORKDIR /opt/vcftools-0.1.14
RUN ./configure --prefix=$VCFTOOLS_INSTALL_DIR && \
  make && \
  make install

RUN rm -f /tmp/vcftools-0.1.14.tar.gz


##########################################################################################
# bedtools v2.26.0
##########################################################################################

WORKDIR /tmp
RUN wget -q https://github.com/arq5x/bedtools2/releases/download/v2.26.0/bedtools-2.26.0.tar.gz && \
  tar -xf bedtools-2.26.0.tar.gz -C /opt
WORKDIR /opt/bedtools2
RUN make && \
  rm -f /tmp/bedtools-2.26.0.tar.gz


##########################################################################################
# bedtools v2.22.0
##########################################################################################

WORKDIR /tmp
RUN wget -q https://github.com/arq5x/bedtools2/releases/download/v2.22.0/bedtools-2.22.0.tar.gz && \
  mkdir /opt/bedtools22/ && \
  tar -xf bedtools-2.22.0.tar.gz -C /opt/bedtools22/
WORKDIR /opt/bedtools22/bedtools2
RUN make && \
  rm -f /tmp/bedtools-2.22.0.tar.gz


##########################################################################################
# bedops  v2.4.27
##########################################################################################

WORKDIR /tmp
RUN mkdir -p /opt/bedops && \
  wget -q https://github.com/bedops/bedops/releases/download/v2.4.27/bedops_linux_x86_64-v2.4.27.tar.bz2 && \
  tar -xf bedops_linux_x86_64-v2.4.27.tar.bz2 -C /opt/bedops
WORKDIR /opt/bedops
RUN rm -f /tmp/bedops_linux_x86_64-v2.4.27.tar.bz2


##########################################################################################
# Bismark 0.16.3
##########################################################################################

WORKDIR /tmp
RUN wget -q http://www.bioinformatics.babraham.ac.uk/projects/bismark/bismark_v0.16.3.tar.gz && \
  tar -xf bismark_v0.16.3.tar.gz -C /opt && \
  rm -f bismark_v0.16.3.tar.gz


##########################################################################################
# tophat2-2.1.1
##########################################################################################

WORKDIR /tmp
RUN wget -q http://ccb.jhu.edu/software/tophat/downloads/tophat-2.1.1.Linux_x86_64.tar.gz && \
  tar -xf tophat-2.1.1.Linux_x86_64.tar.gz -C /opt && \
  rm -f tophat-2.1.1.Linux_x86_64.tar.gz


##########################################################################################
# STAR-2.5.1b
##########################################################################################

WORKDIR /tmp
RUN wget -q https://github.com/alexdobin/STAR/archive/2.5.1b.tar.gz && \
  tar -xf 2.5.1b.tar.gz -C /opt && \
  make STAR -C /opt/STAR-2.5.1b/source && \
  rm -f 2.5.1b.tar.gz


##########################################################################################
# twoBitToFa
##########################################################################################

WORKDIR /opt/2bit
RUN wget -q  http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/twoBitToFa && \
  chmod +x twoBitToFa


##########################################################################################
# hubCheck
##########################################################################################

WORKDIR /opt/hubCheck
RUN wget -q  http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/hubCheck && \
  chmod +x hubCheck


##########################################################################################
# ALEA2 2.0
##########################################################################################

ENV ALEA_INSTALL_DIR=/opt/alea
ENV ALEA_VERSION=2.0
ENV ALEA_FILENAME=alea.${ALEA_VERSION}.tar.gz
#ENV ALEA_URL=https://github.com/hyounesy/ALEA/raw/Task/dist/${ALEA_FILENAME}
#https://github.com/hyounesy/ALEA/raw/master/dist/alea.1.2.2.tar.gz

ENV ALEA_URL=https://github.com/hyounesy/ALEA/raw/master/dist/${ALEA_FILENAME}


WORKDIR /tmp
COPY  ${ALEA_FILENAME} ${ALEA_FILENAME}
RUN wget -q $ALEA_URL
RUN  tar -zxvf ${ALEA_FILENAME} -C /opt && \
  rm -f /tmp/${ALEA_FILENAME}


##########################################################################################
# Create links
##########################################################################################

RUN ln -s /opt/samtools-0.1.19/samtools $ALEA_INSTALL_DIR/bin/ && \
  ln -s $SAMTOOLS_INSTALL_DIR/bin/tabix $ALEA_INSTALL_DIR/bin/tabix && \
  ln -s $SAMTOOLS_INSTALL_DIR/bin/bgzip $ALEA_INSTALL_DIR/bin/bgzip && \
  ln -s /opt/bwa-0.7.12/bwa $ALEA_INSTALL_DIR/bin/ && \
  find /opt/bowtie-1.1.2/ -maxdepth 1 -perm -+x -type f | xargs -i ln -s {} $ALEA_INSTALL_DIR/bin/ && \
  find /opt/bowtie2-2.2.9/ -maxdepth 1 -perm -+x -type f | xargs -i ln -s {} $ALEA_INSTALL_DIR/bin/ && \
  ln -s /opt/bedGraphToBigWig/bedGraphToBigWig $ALEA_INSTALL_DIR/bin/ && \
  ln -s /opt/shapeit/bin/shapeit $ALEA_INSTALL_DIR/bin/ && \
  find $VCFTOOLS_INSTALL_DIR/bin/ -maxdepth 1 -perm -+x -type f | xargs -i ln -s {} $ALEA_INSTALL_DIR/bin/ && \
  find /opt/bedtools2/bin/ -maxdepth 1 -perm -+x -type f | xargs -i ln -s {} $ALEA_INSTALL_DIR/bin/ && \
  find /opt/bedops/bin/ -maxdepth 1 -perm -+x -type f | xargs -i ln -s {} $ALEA_INSTALL_DIR/bin/ && \
  find /opt/bismark_v0.16.3/ -maxdepth 1 -perm -+x -type f | xargs -i ln -s {} $ALEA_INSTALL_DIR/bin/ && \
  find /opt/tophat-2.1.1.Linux_x86_64/ -maxdepth 1 -perm -+x -type f | xargs -i ln -s {} $ALEA_INSTALL_DIR/bin/ && \
  ln -s /opt/STAR-2.5.1b/source/STAR $ALEA_INSTALL_DIR/bin/

#  find /opt/bedtools22/bedtools2/bin/ -maxdepth 1 -perm -+x -type f | xargs -i ln -s {} $ALEA_INSTALL_DIR/bin/bedtools22 && \

RUN  find /opt/bowtie2-2.2.9/ -maxdepth 1 -perm -+x -type f | xargs -i ln -s {} /usr/local/bin/ 

RUN ln -s /opt/2bit/twoBitToFa $ALEA_INSTALL_DIR/bin/ && \
  ln -s /opt/hubCheck/hubCheck $ALEA_INSTALL_DIR/bin/

##########################################################################################
# Clean
##########################################################################################

WORKDIR /
RUN apt-get clean && \
  apt-get remove -y --force-yes --purge \
  bzip2 \
  g++ \
  make \
  pkg-config \
  ncurses-dev \
  build-essential \
  gcc-multilib \
  apt-utils

RUN apt-get -y autoremove


#  wget \
#  zlib1g-dev \

##########################################################################################
# COPY alea.config
##########################################################################################
#
#COPY alea.config /opt/alea/bin/
#
#
##########################################################################################
# ENTRYPOINT
##########################################################################################

ENTRYPOINT ["/opt/alea/bin/alea"]
#ENTRYPOINT ["/bin/bash"]
