# see https://hub.docker.com/_/python
# for basic instructions on how to use the python image

FROM python:3.9

# expectation is to map your local data folder using -v to here
WORKDIR /usr/data

COPY requirements.txt ./
RUN pip install --no-cache-dir -r requirements.txt

# install all the deps we can from apt-get in preference
# builds and installs below were required such that our version is the latest
# since the apt package versions are sometimes way out of date
RUN apt-get -y update && apt-get -y install bcftools minimap2 cmake

# build htslib and install it
RUN cd /home && \
    wget https://github.com/samtools/htslib/releases/download/1.19/htslib-1.19.tar.bz2 && \
    tar xfjp htslib-1.19.tar.bz2 && \
    cd htslib-1.19 && \
    ./configure && make && make install

# build spoa and install it
RUN cd /home && \
    git clone https://github.com/rvaser/spoa && \
    cd spoa && \
    cmake -B build -DCMAKE_BUILD_TYPE=Release && \
    make -C build && \
    make -C build install

# build samtools and install it
RUN cd /home && \
    wget https://github.com/samtools/samtools/releases/download/1.19/samtools-1.19.tar.bz2 && \
    tar xfjp samtools-1.19.tar.bz2 && \
    cd samtools-1.19 && \
    ./configure && make && make install

# this is the arm version - unforunately conda has very poor support for arm
# so forcing linux/x86_64 instead
#RUN cd /tmp && \
#    wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-aarch64.sh && \
#    chmod a+x Miniconda3-latest-Linux-aarch64.sh && \
#    ./Miniconda3-latest-Linux-aarch64.sh -b -p /opt/miniconda3 && \
#    rm Miniconda3-latest-Linux-aarch64.sh

# install latest miniconda3
RUN cd /tmp && \
    wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh && \
    mv Miniconda*.sh miniconda.sh && \
    chmod a+x miniconda.sh && \
    ./miniconda.sh -b -p /opt/miniconda3 && \
    rm miniconda.sh
ENV PATH="/opt/miniconda3/bin:${PATH}"
RUN conda init

# build pipeline-umi-amplicon
RUN cd /usr && \
    git clone https://github.com/nanoporetech/pipeline-umi-amplicon.git && \
    cd pipeline-umi-amplicon && \
    conda env create -f environment.yml
# equivalent to "conda activate pipeline-umi-amplicon"
SHELL ["conda", "run", "-n", "pipeline-umi-amplicon", "/bin/bash", "-c"]
RUN cd /usr/pipeline-umi-amplicon/lib && pip install .

# use the root login profile to default to pipeline-umi-amplicon
RUN echo "conda activate pipeline-umi-amplicon" >> /root/.profile
    
CMD [ "bash", "--login" ]
