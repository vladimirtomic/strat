FROM ubuntu:22.04

WORKDIR /opt

RUN apt update -y && apt upgrade -y && apt install -y curl && apt install -y bzip2
RUN curl -L https://github.com/bwa-mem2/bwa-mem2/releases/download/v2.2.1/bwa-mem2-2.2.1_x64-linux.tar.bz2 | tar jxf -

ENTRYPOINT ["/bin/bash"]
