FROM mambaorg/micromamba:2.3.2-ubuntu22.04

LABEL org.opencontainers.image.title="GLADE RetroPath local service" \
      org.opencontainers.image.version="3.9.1-r20260212" \
      org.opencontainers.image.source="https://github.com/brsynth/retropath2-wrapper"

ARG MAMBA_DOCKERFILE_ACTIVATE=1

USER root
RUN packages="ca-certificates fontconfig libasound2 libdrm2 libfreetype6 libgbm1 libgl1 libgtk-3-0 libnss3 libx11-xcb1 libxcomposite1 libxdamage1 libxkbcommon0 libxrandr2 libxss1 libxtst6 procps unzip" \
    && for attempt in 1 2 3 4 5; do \
        rm -rf /var/lib/apt/lists/*; \
        apt-get -o Acquire::Retries=3 update && break; \
        if [ "${attempt}" = "5" ]; then exit 1; fi; \
        sleep 15; \
    done \
    && for attempt in 1 2 3 4 5; do \
        DEBIAN_FRONTEND=noninteractive apt-get -o Acquire::Retries=5 install -y --fix-missing --no-install-recommends ${packages} && break; \
        if [ "${attempt}" = "5" ]; then exit 1; fi; \
        sleep 15; \
    done \
    && rm -rf /var/lib/apt/lists/* \
    && mkdir -p /opt/knime /opt/service /opt/retropath/rules /var/lib/retropath/jobs \
    && chown -R "${MAMBA_USER}:${MAMBA_USER}" /opt/knime /opt/service /opt/retropath /var/lib/retropath

COPY --chown=${MAMBA_USER}:${MAMBA_USER} environment.yml /tmp/environment.yml

USER ${MAMBA_USER}
RUN micromamba install --yes --name base --file /tmp/environment.yml \
    && micromamba clean --all --yes

# Install the original KNIME platform and the chemistry/RDKit nodes used by
# the RetroPath2.0 workflow. The result is baked into the image so runtime
# jobs never download executables or plugins.
COPY --chown=${MAMBA_USER}:${MAMBA_USER} install_knime.py /tmp/install_knime.py
RUN --mount=type=cache,target=/var/cache/knime,uid=57439,gid=57439 \
    python /tmp/install_knime.py \
      --cache-dir /var/cache/knime \
      --install-dir /opt/knime

# KNIME 4.7 bundles a legacy Chromium binary whose Linux build uses the
# CentOS OpenSSL 1.0 SONAMEs. Keep those ABI-compatible libraries in a
# separate prefix so the service's Python environment can stay on OpenSSL 3.
RUN micromamba create --yes \
      --prefix /home/mambauser/.openssl10 \
      --channel conda-forge \
      openssl=1.0.2k=0 \
    && ln -s libssl.so.1.0.0 /home/mambauser/.openssl10/lib/libssl.so.10 \
    && ln -s libcrypto.so.1.0.0 /home/mambauser/.openssl10/lib/libcrypto.so.10 \
    && micromamba clean --all --yes

RUN LD_LIBRARY_PATH=/home/mambauser/.openssl10/lib:/opt/conda/lib \
    /opt/knime/knime_4.7.0/knime \
      -nosplash \
      -application org.knime.product.KNIME_BATCH_APPLICATION

WORKDIR /opt/service
COPY --chown=${MAMBA_USER}:${MAMBA_USER} app /opt/service/app
COPY --chown=${MAMBA_USER}:${MAMBA_USER} tests /opt/service/tests

ENV RETROPATH_RULES_PATH=/opt/retropath/rules/rr02.csv \
    RETROPATH_RULES_VERSION=rr02-rp2-hs \
    RETROPATH_RULES_SHA256=e24eb97d3172195d03abed6e7da07a4cfd53965553853d126aaa8a93b4bc552f \
    RETROPATH_KNIME_DIR=/opt/knime \
    RETROPATH_WORKFLOW_VERSION=r20260212 \
    RETROPATH_KNIME_VERSION=4.7.0 \
    RETROPATH_RDKIT_PLUGIN_VERSION=4.9.1 \
    RETROPATH_OPENSSL10_DIR=/home/mambauser/.openssl10/lib \
    RETROPATH_JOBS_DIR=/var/lib/retropath/jobs \
    RETROPATH_MAX_QUEUE=8 \
    RETROPATH_JOB_TIMEOUT_SECONDS=3600 \
    PYTHONPATH=/opt/service \
    PYTHONUNBUFFERED=1 \
    LD_LIBRARY_PATH=/home/mambauser/.openssl10/lib:/opt/conda/lib

EXPOSE 8000

CMD ["uvicorn", "app.main:app", "--host", "0.0.0.0", "--port", "8000", "--workers", "1"]
