FROM ubuntu:22.04

RUN apt-get update -y && apt-get upgrade -y && apt-get install wget -y

WORKDIR /opt

ENV VERSION="0.2.2"
RUN wget https://github.com/hasindu2008/squigulator/releases/download/v${VERSION}/squigulator-v${VERSION}-x86_64-linux-binaries.tar.gz
RUN tar xf squigulator-v${VERSION}-x86_64-linux-binaries.tar.gz

WORKDIR /opt/squigulator-v${VERSION}

ENTRYPOINT [ "./squigulator" ]
