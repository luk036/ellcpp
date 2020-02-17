FROM gitpod/workspace-full

USER root
# Install util tools.
RUN apt-get update \
 && apt-get install -y \
  apt-utils \
  sudo \
  git \
  less \
  wget

RUN mkdir -p /workspace/data \
    && chown -R gitpod:gitpod /workspace/data
  
RUN mkdir /home/gitpod/.conda
# Install conda
RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh && \
    /bin/bash ~/miniconda.sh -b -p /opt/conda && \
    rm ~/miniconda.sh && \
    ln -s /opt/conda/etc/profile.d/conda.sh /etc/profile.d/conda.sh && \
    echo ". /opt/conda/etc/profile.d/conda.sh" >> ~/.bashrc && \
    echo "conda activate base" >> ~/.bashrc
    
RUN chown -R gitpod:gitpod /opt/conda \
    && chmod -R 777 /opt/conda \
    && chown -R gitpod:gitpod /home/gitpod/.conda \
    && chmod -R 777 /home/gitpod/.conda

RUN /opt/conda/bin/conda config --set always_yes yes --set changeps1 no \
    && /opt/conda/bin/conda update -q conda \
    && /opt/conda/bin/conda info -a

RUN /opt/conda/bin/conda install -y \
    ninja \
    lcov

RUN /opt/conda/bin/conda install -y -c conda-forge \
    catch2 \
    fmt \
    lapack \
    libboost \
    openblas \
    spdlog \
    xtensor-fftw=0.2.5 \
    xtensor-blas=0.16.1 \
    xtensor=0.20.10 \
    cppcheck

RUN apt-get clean && rm -rf /var/cache/apt/* && rm -rf /var/lib/apt/lists/* && rm -rf /tmp/*
