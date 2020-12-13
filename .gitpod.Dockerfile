FROM gitpod/workspace-full

USER root
# Install util tools.

RUN apt-get update \
 && apt-get install -y \
  apt-utils \
  sudo \
  git \
  less \
  neofetch \
  neovim \
  asciinema \
  tmux \
  vifm \
  w3m \
  cowsay \
  figlet \
  fortune \
  toilet \
  tty-clock \
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
    
RUN /opt/conda/bin/conda config --set always_yes yes --set changeps1 no \
    && /opt/conda/bin/conda update -q conda \
    && /opt/conda/bin/conda info -a

RUN /opt/conda/bin/conda install -y \
    ninja

RUN /opt/conda/bin/conda install -y -c conda-forge \
    benchmark \
    boost \
    catch2 \
    fmt \
    lapack \
    openblas \
    spdlog \
    xtensor-fftw \
    xtensor-blas \
    xtensor \
    cppcheck

RUN chown -R gitpod:gitpod /opt/conda \
    && chmod -R 777 /opt/conda \
    && chown -R gitpod:gitpod /home/gitpod/.conda \
    && chmod -R 777 /home/gitpod/.conda

RUN apt-get clean && rm -rf /var/cache/apt/* && rm -rf /var/lib/apt/lists/* && rm -rf /tmp/*
