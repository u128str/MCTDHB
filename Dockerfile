## RUN MCTDHB through docker: docker run --rm -it -v $(pwd):/tmp mctdhb
## RESTART: sudo service docker restart
## ADD DNS to sudo vim /etc/default/docker
##  DOCKER_OPTS="--dns 141.51.8.4 --dns 8.8.8.8 --dns 8.8.4.4"

##rm all images ..
## sudo docker ps -a | awk '{print $1}' | grep -v CONTAINER | xargs sudo docker rm
## sudo docker images | grep "<none>" | awk '{print $3}' | xargs sudo docker rmi

## sudo vim /etc/resolv.conf  and put: 
## nameserver 141.51.8.4
## nameserver 8.8.8.8
## nameserver 8.8.4.4

## Build:   docker build --no-cache -t u128str/mctdhb .
## Run:     docker run --rm -it -v $(pwd):/tmp u128str/mctdhb

#  apt-get install make
#  apt-get install openmpi-bin libopenmpi-dev
#  apt-get install fftw3  fftw3-dev
#  apt-get install libblas-dev liblapack-dev
# Version: 0.0.1

#FROM ubuntu:16.04
#MAINTAINER  Alexej I. Streltsov  <u128str@gmail.com>
#RUN apt-get update && apt-get install -y  vim make  openmpi-bin libopenmpi-dev fftw3  fftw3-dev libblas-dev liblapack-dev

## RUN MCTDHB through docker: docker run --rm -it -v $(pwd):/tmp mctdhb
FROM mctdhb/minunix:latest
MAINTAINER  Alexej I. Streltsov  <u128str@gmail.com>
##COPY ./MCTDHB_V3.3.01  /mctdhb
#COPY Makefile  /mctdhb
#WORKDIR /mctdhb/
RUN  make
##COPY ./MCTDHB_V3.3.01  /mctdhb
#WORKDIR /tmp
#CMD  ["/mctdhb/bin/boson_MCTDHB_gnu_FFTW"]
