# containers/Dockerfile
# Reproducible base + system tools (works on amd64/arm64)
FROM ubuntu:22.04

ENV DEBIAN_FRONTEND=noninteractive \
    PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    LC_ALL=C.UTF-8 \
    LANG=C.UTF-8

SHELL ["/bin/bash", "-lc"]

# 1) System deps + bio tools
RUN apt-get update && apt-get install -y --no-install-recommends \
      ca-certificates curl wget git \
      python3 python3-pip python3-dev \
      openjdk-17-jre-headless \
      build-essential \
      ncbi-blast+ \
      prodigal \
      diamond-aligner \
    && rm -rf /var/lib/apt/lists/*

# --- Nextflow baked-in & writable for any UID ---
ARG NXF_VER=25.04.6
ENV NXF_HOME=/opt/.nextflow \
    NXF_TEMP=/tmp \
    HOME=/opt/saphari_home

# Install launcher and prefetch the exact JAR (no network at runtime)
RUN curl -fsSL https://get.nextflow.io -o /usr/local/bin/nextflow \
 && chmod 0755 /usr/local/bin/nextflow \
 && chown root:root /usr/local/bin/nextflow \
 && mkdir -p ${NXF_HOME}/framework/${NXF_VER} ${HOME} \
 && curl -fsSL https://www.nextflow.io/releases/v${NXF_VER}/nextflow-${NXF_VER}-one.jar \
      -o ${NXF_HOME}/framework/${NXF_VER}/nextflow-${NXF_VER}-one.jar \
 && chmod 0644 ${NXF_HOME}/framework/${NXF_VER}/nextflow-${NXF_VER}-one.jar \
 && chmod -R 0777 ${NXF_HOME} ${HOME}

# 3) Install your Python package from src/ (correct for src/ layout)
#    Copy only what’s needed for better Docker layer caching
WORKDIR /opt/saphari
COPY pyproject.toml ./pyproject.toml
COPY src ./src

# If you have runtime data outside src/ that must be included, copy it here.
# (Families.json is already under src/saphari/data/, so no extra copy needed.)

RUN python3 -m pip install --no-cache-dir --upgrade pip \
 && python3 -m pip install --no-cache-dir .

# 4) Create a neutral working dir for mounted volumes
WORKDIR /work

# Keep the container flexible for Nextflow & ad-hoc debugging.
# `saphari` is now on PATH so both of these work:
#   docker run --rm saphari:0.1.0 saphari --help
#   docker run -it --rm saphari:0.1.0 /bin/bash
CMD ["/bin/bash"]

# Optional: a quick health check for CI (uncomment if you like)
# HEALTHCHECK --interval=30s --timeout=10s --retries=3 \
#   CMD saphari --help >/dev/null 2>&1 || exit 1
