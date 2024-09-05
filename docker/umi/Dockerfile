# see https://hub.docker.com/_/python
# for basic instructions on how to use the python image

FROM python:3.9

WORKDIR /usr/data

COPY requirements.txt ./
RUN pip install --no-cache-dir -r requirements.txt

RUN apt-get -y update && apt-get -y install bcftools minimap2 cmake less dos2unix

# build spoa and install it
RUN cd /home && \
    git clone https://github.com/rvaser/spoa && \
    cd spoa && \
    cmake -B build -DCMAKE_BUILD_TYPE=Release && \
    make -C build && \
    make -C build install

# build meshclust from HersheyTriesToCode repo
RUN cd /usr/local/src && \
    git clone https://github.com/HersheyTriesToCode/Identity.git && \
    cd Identity/bin && \
    cmake .. && \
    make && \
    cp -pv identity meshclust *.so /usr/local/bin

CMD [ "bash" ]