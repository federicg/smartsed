## Start from the official Ubuntu 18.04 image
FROM ubuntu:19.04

WORKDIR /smartsed

RUN apt-get update 

# -y means that the answer to the question is yes!
RUN apt-get install -y gdal-bin python-gdal python3-gdal


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





RUN wget -P /smartsed "https://cran.r-project.org/src/contrib/tensorA_0.36.1.tar.gz"
RUN R CMD INSTALL tensorA_0.36.1.tar.gz

RUN wget -P /smartsed "https://cran.r-project.org/src/contrib/DEoptimR_1.0-8.tar.gz"
RUN R CMD INSTALL DEoptimR_1.0-8.tar.gz

RUN wget -P /smartsed "https://cran.r-project.org/src/contrib/robustbase_0.93-5.tar.gz"
RUN R CMD INSTALL robustbase_0.93-5.tar.gz

RUN wget -P /smartsed "https://cran.r-project.org/src/contrib/sp_1.3-2.tar.gz"
RUN R CMD INSTALL sp_1.3-2.tar.gz

RUN wget -P /smartsed "https://cran.r-project.org/src/contrib/Rcpp_1.0.3.tar.gz"
RUN R CMD INSTALL Rcpp_1.0.3.tar.gz

RUN wget -P /smartsed "https://cran.r-project.org/src/contrib/raster_3.0-12.tar.gz"
RUN R CMD INSTALL raster_3.0-12.tar.gz

RUN wget -P /smartsed "https://cran.r-project.org/src/contrib/RcppArmadillo_0.9.850.1.0.tar.gz"
RUN R CMD INSTALL RcppArmadillo_0.9.850.1.0.tar.gz

RUN wget -P /smartsed "https://cran.r-project.org/src/contrib/bayesm_3.1-4.tar.gz"
RUN R CMD INSTALL bayesm_3.1-4.tar.gz

RUN wget -P /smartsed "https://cran.r-project.org/src/contrib/compositions_1.40-3.tar.gz"
RUN R CMD INSTALL compositions_1.40-3.tar.gz

RUN apt-get install -y libgdal-dev


RUN wget -P /smartsed "https://cran.r-project.org/src/contrib/rgdal_1.4-8.tar.gz"
RUN R CMD INSTALL rgdal_1.4-8.tar.gz


RUN wget -P /smartsed "https://cran.r-project.org/src/contrib/zoo_1.8-7.tar.gz"
RUN R CMD INSTALL zoo_1.8-7.tar.gz

RUN wget -P /smartsed "https://cran.r-project.org/src/contrib/xts_0.12-0.tar.gz"
RUN R CMD INSTALL xts_0.12-0.tar.gz

RUN wget -P /smartsed "https://cran.r-project.org/src/contrib/intervals_0.15.1.tar.gz"
RUN R CMD INSTALL intervals_0.15.1.tar.gz

RUN wget -P /smartsed "https://cran.r-project.org/src/contrib/spacetime_1.2-3.tar.gz"
RUN R CMD INSTALL spacetime_1.2-3.tar.gz

RUN wget -P /smartsed "https://cran.r-project.org/src/contrib/FNN_1.1.3.tar.gz"
RUN R CMD INSTALL FNN_1.1.3.tar.gz

RUN wget -P /smartsed "https://cran.r-project.org/src/contrib/gstat_2.0-4.tar.gz"
RUN R CMD INSTALL gstat_2.0-4.tar.gz


RUN wget -P /smartsed "https://cran.r-project.org/src/contrib/dotCall64_1.0-0.tar.gz"
RUN R CMD INSTALL dotCall64_1.0-0.tar.gz

RUN wget -P /smartsed "https://cran.r-project.org/src/contrib/spam_2.5-1.tar.gz"
RUN R CMD INSTALL spam_2.5-1.tar.gz

RUN wget -P /smartsed "https://cran.r-project.org/src/contrib/maps_3.3.0.tar.gz"
RUN R CMD INSTALL maps_3.3.0.tar.gz

RUN wget -P /smartsed "https://cran.r-project.org/src/contrib/fields_10.3.tar.gz"
RUN R CMD INSTALL fields_10.3.tar.gz

RUN apt-get install -y python3-pip

RUN pip3 install richdem








RUN rm -r *.tar.gz

# To build the Docker container
#docker build -t smartsed_env -f ./Dockerfile ./
