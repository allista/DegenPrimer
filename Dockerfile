FROM python:2.7.18
LABEL tag=degen-primer
ENV PYTHONUNBUFFERED=1

WORKDIR /DegenPrimer
ADD requirements-external.txt /DegenPrimer
ADD BioUtils/requirements.txt /DegenPrimer/BioUtils/requirements.txt
RUN pip install -r requirements-external.txt
ADD . /DegenPrimer/
RUN pip install -v ./BioUtils
RUN pip install -v .

VOLUME /input /output

WORKDIR /output
ENTRYPOINT ["degen_primer"]
