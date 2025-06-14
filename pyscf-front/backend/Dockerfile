# PySCF_Front Backend Dockerfile
FROM python:3.10-slim

# Set work directory
WORKDIR /app

# Install system dependencies
RUN apt-get update && apt-get install -y \
    gcc \
    g++ \
    gfortran \
    libopenblas-dev \
    liblapack-dev \
    libhdf5-dev \
    libhdf5-serial-dev \
    cmake \
    git \
    pkg-config \
    curl \
    && rm -rf /var/lib/apt/lists/*

# Create non-root user
RUN useradd --create-home --shell /bin/bash app
RUN chown -R app:app /app
USER app

# Copy requirements first for better caching
COPY --chown=app:app requirements.txt .

# Install Python dependencies
RUN pip install --user --no-cache-dir -r requirements.txt

# Copy source code
COPY --chown=app:app src/ ./src/
COPY --chown=app:app protos/ ./protos/

# Set Python path
ENV PYTHONPATH="/app/src"
ENV PATH="/home/app/.local/bin:$PATH"

# Generate protobuf files
RUN python -m grpc_tools.protoc \
    --proto_path=protos \
    --python_out=src/generated \
    --grpc_python_out=src/generated \
    protos/calculation.proto

# Expose port
EXPOSE 50051

# Health check
HEALTHCHECK --interval=30s --timeout=10s --start-period=5s --retries=3 \
    CMD python -c "import grpc; channel = grpc.insecure_channel('localhost:50051'); grpc.channel_ready_future(channel).result(timeout=10)"

# Run the application
CMD ["python", "src/main.py"]