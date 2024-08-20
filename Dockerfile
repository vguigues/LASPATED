FROM ubuntu:22.04
ARG USE_GUROBI=0
ARG GRB_VERSION=11.0.1
ARG GRB_SHORT_VERSION=11.0
ARG TARGETPLATFORM
ARG PROJECT_DIR=/LASPATED
ENV GUROBI_HOME=/opt/gurobi/linux
ENV PATH=$PATH:$GUROBI_HOME/bin
ENV LD_LIBRARY_PATH=$GUROBI_HOME/lib
COPY . $PROJECT_DIR

# Create a script to determine GRB_PLATFORM based on TARGETPLATFORM
RUN if [ "$TARGETPLATFORM" = "linux/arm64" ]; then \
        echo "armlinux64" > /platform.txt; \
      else \
        echo "linux64" > /platform.txt; \
    fi

# install gurobi package and copy the files
WORKDIR /opt

RUN apt-get update && apt-get install -y python3 python3-pip libboost-all-dev \
        build-essential vim gdal-bin libgdal-dev
RUN pip install --upgrade pip



# update system and certificates
RUN apt-get install --no-install-recommends -y\
       ca-certificates  \
       p7zip-full \
       zip \
    && update-ca-certificates \
    && python3 -m pip install gurobipy==${GRB_VERSION} \
    && rm -rf /var/lib/apt/lists/*

RUN if [ "$USE_GUROBI" = "1" ]; then \
        export GRB_PLATFORM=$(cat /platform.txt) && echo $GRB_PLATFORM \
        && apt-get update \
        && apt-get install --no-install-recommends -y\
        ca-certificates  \
        wget \
        && update-ca-certificates \
        && wget -v https://packages.gurobi.com/${GRB_SHORT_VERSION}/gurobi${GRB_VERSION}_$GRB_PLATFORM.tar.gz \
        && tar -xvf gurobi${GRB_VERSION}_$GRB_PLATFORM.tar.gz  \
        && rm -f gurobi${GRB_VERSION}_$GRB_PLATFORM.tar.gz \
        && mv -f gurobi* gurobi \
        && rm -rf gurobi/$GRB_PLATFORM/docs \
        && mv -f gurobi/$GRB_PLATFORM*  gurobi/linux; \
    fi

    
WORKDIR $PROJECT_DIR
RUN pip3 install ${PROJECT_DIR}/laspated/.

RUN ln -s $(which python3) /usr/local/bin/python

RUN if [ "${USE_GUROBI}" = "1" ]; then \
        make -C Model_Calibration/Cpp USE_GUROBI=1 GUROBI_VER=110 && \
        make -C Model_Calibration/Cpp USE_GUROBI=1 GUROBI_VER=110 test; \
    else \
        make -C Model_Calibration/Cpp USE_GUROBI=0 && \
        make -C Model_Calibration/Cpp USE_GUROBI=0 test; \
    fi

RUN if [ ${USE_GUROBI} = 1 ]; then \ 
        make -C Replication/cpp_tests USE_GUROBI=1 GUROBI_VER=110; \
    else \
        make -C Replication/cpp_tests USE_GUROBI=0; \
    fi