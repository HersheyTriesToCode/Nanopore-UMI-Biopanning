# see https://hub.docker.com/_/python
# for basic instructions on how to use the python image

FROM python:3.9

WORKDIR /usr/data

COPY requirements.txt ./
RUN pip install --no-cache-dir -r requirements.txt

CMD [ "bash" ]