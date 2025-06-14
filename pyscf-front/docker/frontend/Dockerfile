FROM ubuntu:22.04

# Install system dependencies
RUN apt-get update && apt-get install -y \
    curl \
    git \
    unzip \
    xz-utils \
    zip \
    libglu1-mesa \
    cmake \
    ninja-build \
    build-essential \
    libgtk-3-dev \
    liblzma-dev \
    && rm -rf /var/lib/apt/lists/*

# Install Flutter
RUN git clone https://github.com/flutter/flutter.git -b stable /flutter
ENV PATH="/flutter/bin:${PATH}"

# Enable desktop support
RUN flutter config --enable-linux-desktop
RUN flutter config --no-analytics

# Set working directory
WORKDIR /app

# Copy pubspec files first for better caching
COPY frontend/pubspec.yaml frontend/pubspec.lock ./

# Get dependencies
RUN flutter pub get

# Copy the rest of the application
COPY frontend/ .

# Build the application
RUN flutter build linux --release

# Expose port for development server (if needed)
EXPOSE 8080

# Default command
CMD ["./build/linux/x64/release/bundle/pyscf_front"]