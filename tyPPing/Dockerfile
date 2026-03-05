FROM mambaorg/micromamba:latest

USER root

# System dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    python3 \
    python3-pip \
    unzip \
    && rm -rf /var/lib/apt/lists/*

# Install HMMER via micromamba
RUN micromamba install -y -n base -c bioconda -c conda-forge hmmer \
    && micromamba clean --all --yes
ENV PATH="/opt/conda/bin:$PATH"

# Copy tyPPing package and test example for verification
COPY pyproject.toml /opt/tyPPing/
COPY src/ /opt/tyPPing/src/
COPY tyPPing_input_data/ /opt/tyPPing/tyPPing_input_data/
COPY tyPPing_signature_profiles.zip /opt/tyPPing/
COPY test_example/ /opt/tyPPing/test_example/

# Install tyPPing Python package
RUN pip install --break-system-packages /opt/tyPPing/

# Unzip HMM profiles
RUN unzip /opt/tyPPing/tyPPing_signature_profiles.zip -d /opt/tyPPing/ \
    && rm /opt/tyPPing/tyPPing_signature_profiles.zip

# WORKDIR /data

ENTRYPOINT ["typping"]
