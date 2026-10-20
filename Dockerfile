# Dockerfile for building PySCF_front Linux distribution
# This creates a complete build environment with Miniforge, Node.js, and all dependencies

FROM --platform=linux/amd64 ubuntu:22.04

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

# Install Miniforge for the linux/amd64 AppImage target
RUN test "$(uname -m)" = "x86_64" && \
    curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh" && \
    bash Miniforge3-Linux-x86_64.sh -b -p /root/miniforge3 && \
    rm Miniforge3-Linux-x86_64.sh

# Set up conda environment
ENV PATH="/root/miniforge3/bin:${PATH}"
ENV CONDA_DEFAULT_ENV=pyscf-env
ENV HOME=/root

# Create working directory
WORKDIR /app

# Copy lockfile first for better caching (matches CI/release lock-based install)
COPY .github/pyscf-env.conda-lock.yml .github/pyscf-env.conda-lock.yml

# Configure conda for better network stability
RUN /root/miniforge3/bin/conda config --set remote_connect_timeout_secs 30.0 && \
    /root/miniforge3/bin/conda config --set remote_read_timeout_secs 120.0 && \
    /root/miniforge3/bin/conda config --set remote_max_retries 5

# Install conda-lock (pinned to match CI) and create environment from lockfile
RUN /root/miniforge3/bin/conda install -n base -c conda-forge conda-lock=3.0.4 -y && \
    for i in 1 2 3; do \
        conda-lock install --name pyscf-env .github/pyscf-env.conda-lock.yml && break || { \
            if [ "$i" -eq 3 ]; then \
                echo "conda-lock install failed after 3 attempts"; \
                exit 1; \
            fi; \
            echo "Attempt $i failed, retrying..."; \
            sleep 5; \
        }; \
    done && \
    /root/miniforge3/bin/conda run -n pyscf-env python --version && \
    /root/miniforge3/bin/conda clean -afy

# Activate conda environment in shell
SHELL ["/bin/bash", "-c"]
RUN echo "source /root/miniforge3/etc/profile.d/conda.sh && conda activate pyscf-env" >> ~/.bashrc

# Copy package files for dependency installation
COPY package*.json ./

# Install Node.js dependencies
RUN npm config set fetch-retries 5 && \
    npm config set fetch-retry-mintimeout 20000 && \
    npm config set fetch-retry-maxtimeout 120000 && \
    npm config set fetch-timeout 300000 && \
    for i in 1 2 3; do \
        npm ci --legacy-peer-deps && break || { \
            if [ "$i" -eq 3 ]; then \
                echo "npm ci failed after 3 attempts"; \
                exit 1; \
            fi; \
            echo "npm ci failed, retrying (attempt $i/3)..."; \
            sleep 10; \
        }; \
    done

# Copy the entire project
COPY . .

# Set default command to build the application
CMD ["bash", "-c", "source /root/miniforge3/etc/profile.d/conda.sh && conda activate pyscf-env && npm run package:linux"]
