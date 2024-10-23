FROM ubuntu:latest
RUN apt-get update -q -y
RUN apt-get install samtools bedtools default-jre r-base python3 python3-pip curl build-essential -y
RUN python3 -m pip config set global.break-system-packages true
RUN ln -sf /usr/bin/python3 /usr/bin/python; ln -sf /usr/bin/pip3 /usr/bin/pip
RUN pip install MACS3
RUN cd /usr/bin; curl -LJO https://github.com/broadinstitute/picard/releases/download/3.3.0/picard.jar
RUN echo "\npython3 -m pip config set global.break-system-packages true\n" >> /etc/bash.bashrc
RUN echo 'picard() {\n\tjava -Xmx$1 -jar /usr/bin/picard.jar "${@:2}"\n}\n' >> /etc/bash.bashrc