# Dockerfile for building PySCF_front Linux distribution
# This creates a complete build environment with Miniforge, Node.js, and all dependencies

FROM ubuntu:22.04

# Prevent interactive prompts during package installation
ENV DEBIAN_FRONTEND=noninteractive

# Install system dependencies
RUN apt-get update && apt-get install -y \
    curl \
    git \
    build-essential \
    libx11-dev \
    libxext-dev \
    libxi-dev \
    libxrender-dev \
    libxrandr-dev \
    libxcursor-dev \
    libxinerama-dev \
    libgl1-mesa-dev \
    libglu1-mesa-dev \
    libasound2-dev \
    libpulse-dev \
    libudev-dev \
    libdbus-1-dev \
    libglib2.0-dev \
    libgtk-3-dev \
    libnss3-dev \
    libatk1.0-dev \
    libatk-bridge2.0-dev \
    libcups2-dev \
    libdrm-dev \
    libgbm-dev \
    libxkbcommon-dev \
    ca-certificates \
    && rm -rf /var/lib/apt/lists/*

# Install Node.js 20.x
RUN curl -fsSL https://deb.nodesource.com/setup_20.x | bash - && \
    apt-get install -y nodejs && \
    rm -rf /var/lib/apt/lists/*

# Install Miniforge
RUN curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh" && \
    bash Miniforge3-Linux-x86_64.sh -b -p /root/miniforge3 && \
    rm Miniforge3-Linux-x86_64.sh

# Set up conda environment
ENV PATH="/root/miniforge3/bin:${PATH}"
ENV CONDA_DEFAULT_ENV=pyscf-env
ENV HOME=/root

# Create working directory
WORKDIR /app

# Copy environment file first for better caching
COPY .github/environment.yml .github/environment.yml

# Configure conda for better network stability
RUN /root/miniforge3/bin/conda config --set remote_connect_timeout_secs 30.0 && \
    /root/miniforge3/bin/conda config --set remote_read_timeout_secs 120.0 && \
    /root/miniforge3/bin/conda config --set remote_max_retries 5

# Create conda environment with retry logic
RUN for i in 1 2 3; do \
        /root/miniforge3/bin/conda env create -f .github/environment.yml && break || \
        (echo "Attempt $i failed, retrying..." && sleep 5); \
    done && \
    /root/miniforge3/bin/conda clean -afy

# Activate conda environment in shell
SHELL ["/bin/bash", "-c"]
RUN echo "source /root/miniforge3/etc/profile.d/conda.sh && conda activate pyscf-env" >> ~/.bashrc

# Copy package files for dependency installation
COPY package*.json ./

# Install Node.js dependencies
RUN npm ci

# Copy the entire project
COPY . .

# Set default command to build the application
CMD ["bash", "-c", "source /root/miniforge3/etc/profile.d/conda.sh && conda activate pyscf-env && npm run package:linux"]
