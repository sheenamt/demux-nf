FROM ubuntu:18.04
ENV DEBIAN_FRONTEND noninteractive

RUN apt-get update
RUN apt-get -y install git cmake gcc g++ make
#RUN apk --no-cache add git cmake gcc g++ make zlib-dev

RUN git clone https://github.com/Illumina/interop.git /build/interop \
    && mkdir /build/illumina-interop-build/ \
    && mkdir /build/illumina-interop/ \
    && cd /build/illumina-interop-build/ \
    && cmake ../interop -DCMAKE_INSTALL_PREFIX=/build/illumina-interop/ \
    && cmake --build . \
    && cmake --build . --target install

RUN apt-get install zlib1g-dev

# Seqtk
RUN git clone https://github.com/lh3/seqtk.git /build/seqtk \
    && make -C /build/seqtk

# Python and python dependencies
RUN apt-get install -y rpm python3 python3-pip
RUN pip3 install --upgrade pip setuptools
RUN pip3 install sample-sheet==0.11.0 pandas==0.25.2


# bcl2fastq
ADD https://s3-us-west-2.amazonaws.com/management-ncglsoftwarerepobucket-zeblyt3py6dz/bcl2fastq/bcl2fastq2-v2.20.0.422-Linux-x86_64.rpm /bcl2fastq.rpm
RUN rpm --nodeps -i /bcl2fastq.rpm

# AWS tool chain
RUN pip3 --no-cache-dir install awscli


