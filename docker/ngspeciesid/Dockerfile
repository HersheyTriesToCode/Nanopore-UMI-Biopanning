# see https://hub.docker.com/_/python
# for basic instructions on how to use the python image

FROM python:3.9

WORKDIR /usr/data

COPY requirements.txt ./
RUN pip install --no-cache-dir -r requirements.txt

# install all the deps we can from apt-get
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

CMD [ "bash" ]