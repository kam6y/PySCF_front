FROM python:3.13-slim

# Install system dependencies
RUN apt-get update && apt-get install -y \
    build-essential \
    gfortran \
    libopenblas-dev \
    liblapack-dev \
    libhdf5-dev \
    pkg-config \
    cmake \
    git \
    curl \
    && rm -rf /var/lib/apt/lists/*

# Set environment variables
ENV PYTHONDONTWRITEBYTECODE=1
ENV PYTHONUNBUFFERED=1
ENV PYTHONPATH=/app
ENV PYSCF_USE_GPU=false
ENV OMP_NUM_THREADS=4

# Set working directory
WORKDIR /app

# Copy requirements first for better caching
COPY backend/requirements.txt requirements.txt

# Install Python dependencies
RUN pip install --no-cache-dir --upgrade pip && \
    pip install --no-cache-dir -r requirements.txt

# Copy application code
COPY backend/ .
COPY protos/ /protos/

# Generate gRPC stubs
RUN mkdir -p src/grpc_stubs && \
    python -m grpc_tools.protoc \
    -I=/protos \
    --python_out=src/grpc_stubs \
    --grpc_python_out=src/grpc_stubs \
    /protos/calculation.proto && \
    touch src/grpc_stubs/__init__.py

# Fix the generated import in grpc stub
RUN sed -i 's/import calculation_pb2 as calculation__pb2/from . import calculation_pb2 as calculation__pb2/' src/grpc_stubs/calculation_pb2_grpc.py

# Create necessary directories
RUN mkdir -p /app/data /app/logs /app/temp

# Set proper permissions
RUN chmod -R 755 /app

# Expose ports
EXPOSE 8000 50051

# Health check
HEALTHCHECK --interval=30s --timeout=30s --start-period=5s --retries=3 \
    CMD curl -f http://localhost:8000/health || exit 1

# Production command
CMD ["python", "-m", "uvicorn", "src.main:app", "--host", "0.0.0.0", "--port", "8000"]