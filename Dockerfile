FROM python:3.11

WORKDIR /umivar

COPY . ./

RUN apt-get update && apt-get install -y gcc g++ make git

RUN git clone https://github.com/samtools/htslib && \
	cd htslib && \
    git submodule update --init --recursive && \
	make && \
	make install && \
	cd .. && \
	git clone https://github.com/samtools/samtools && \
	cd samtools && \
	make && \
	make install && \
	cd .. && \
	rm -rf htslib samtools 

RUN pip install -r requirements.txt 

RUN apt-get purge -y --auto-remove gcc g++ make git

ENTRYPOINT ["python3", "bin/umivar.py"]