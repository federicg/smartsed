## Start from the official Ubuntu 
FROM ubuntu:20.04

WORKDIR /smartsed

RUN apt-get update 

# -y means that the answer to the question is yes!
RUN apt-get install -y gdal-bin python3-gdal


RUN apt-get install -y g++

RUN apt-get install make -y

RUN apt-get install -y tzdata
RUN apt-get install -y wget
RUN apt-get install -y gpg

#RUN echo "deb https://cloud.r-project.org/bin/linux/ubuntu xenial-cran35/" | tee /etc/apt/sources.list.d/docker.list

#RUN apt-get -y install --no-install-recommends --no-install-suggests \
#        ca-certificates software-properties-common gnupg2 gnupg1 \
#      && apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9 \
#      && add-apt-repository 'deb https://cloud.r-project.org/bin/linux/ubuntu xenial-cran35/' \
#      && apt-get install r-base

#RUN gpg --keyserver hkp://keyserver.ubuntu.com:80 --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9
#RUN gpg -a --export E298A3A825C0D65DFD57CBB651716619E084DAB9 | apt-key add -
#RUN apt-get update

#RUN deb https://cloud.r-project.org/bin/linux/ubuntu xenial-cran35/
RUN apt-get install -y r-base


RUN apt-get update
 
RUN Rscript -e "packages = c('raster','gstat','compositions', 'dissever', 'fields', 'soiltexture', 'viridis','psych','latex2exp'); for (pck in packages) {if ( ! (pck %in% installed.packages()[,'Package']) ) install.packages(pck) }"


RUN apt-get update
RUN apt-get install -y python3-pip

RUN pip3 install richdem

RUN apt-get install -y libopenmpi-dev


RUN apt-get install -y libgdal-dev

RUN Rscript -e "packages = c('raster','dissever'); for (pck in packages) {if ( ! (pck %in% installed.packages()[,'Package']) ) install.packages(pck) }"

RUN Rscript -e "packages = c('rgdal'); for (pck in packages) {if ( ! (pck %in% installed.packages()[,'Package']) ) install.packages(pck) }"



#RUN wget http://tau.uoregon.edu/tau.tgz
#RUN tar xzvf tau.tgz
#RUN cd tau
#RUN ./configure -cc=gcc -c++=g++ -openmp -mpi
#RUN make -j install
#RUN ./configure -cc=nvc -c++=nvc++ -openmp -mpi

#RUN rm -r *.tar.gz



#RUN apt-get update
#RUN apt-get install -y libboost-all-dev

RUN apt-get install -y vim

ENV OMPI_ALLOW_RUN_AS_ROOT=1
ENV OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1
ENV LD_LIBRARY_PATH=/smartsed/DeterministicProgram/include/lib:$LD_LIBRARY_PATH




# To build the Docker container
#docker build -t smartsed_env -f ./Dockerfile ./
