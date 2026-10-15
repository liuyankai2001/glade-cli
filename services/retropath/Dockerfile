FROM mambaorg/micromamba:2.3.2-ubuntu22.04

LABEL org.opencontainers.image.title="GLADE RetroPath local service" \
      org.opencontainers.image.version="3.9.1-r20260212-memoryfix2-cpu2" \
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

COPY --chown=${MAMBA_USER}:${MAMBA_USER} services/retropath/environment-linux-64.lock /tmp/environment.lock

USER ${MAMBA_USER}
RUN micromamba install --yes --name base --file /tmp/environment.lock \
    && micromamba clean --all --yes

# Install the original KNIME platform and the chemistry/RDKit nodes used by
# the RetroPath2.0 workflow. The result is baked into the image so runtime
# jobs never download executables or plugins.
COPY --chown=${MAMBA_USER}:${MAMBA_USER} services/retropath/install_knime.py /tmp/install_knime.py
RUN --mount=type=cache,target=/var/cache/knime,uid=57439,gid=57439 \
    python /tmp/install_knime.py \
      --cache-dir /var/cache/knime \
      --install-dir /opt/knime \
    && sed -i -E 's/^-Xmx.*/-Xmx1024m/' /opt/knime/knime_4.7.0/knime.ini \
    && grep -qx -- '-Xmx1024m' /opt/knime/knime_4.7.0/knime.ini

# KNIME 4.7 bundles a legacy Chromium binary whose Linux build uses the
# CentOS OpenSSL 1.0 SONAMEs. Keep those ABI-compatible libraries in a
# separate prefix so the service's Python environment can stay on OpenSSL 3.
COPY --chown=${MAMBA_USER}:${MAMBA_USER} services/retropath/openssl10-linux-64.lock /tmp/openssl10.lock
RUN micromamba create --yes \
      --prefix /home/mambauser/.openssl10 \
      --file /tmp/openssl10.lock \
    && ln -s libssl.so.1.0.0 /home/mambauser/.openssl10/lib/libssl.so.10 \
    && ln -s libcrypto.so.1.0.0 /home/mambauser/.openssl10/lib/libcrypto.so.10 \
    && micromamba clean --all --yes

RUN LD_LIBRARY_PATH=/home/mambauser/.openssl10/lib:/opt/conda/lib \
    /opt/knime/knime_4.7.0/knime \
      -nosplash \
      -application org.knime.product.KNIME_BATCH_APPLICATION

WORKDIR /opt/service
COPY --chown=${MAMBA_USER}:${MAMBA_USER} services/retropath/app /opt/service/app
COPY --chown=${MAMBA_USER}:${MAMBA_USER} services/retropath/tests /opt/service/tests
COPY --chown=${MAMBA_USER}:${MAMBA_USER} services/retropath/memoryfix /opt/retropath/memoryfix
RUN /opt/conda/bin/python /opt/retropath/memoryfix/build.py
COPY --chown=${MAMBA_USER}:${MAMBA_USER} data/retropath/rules/rr02/retrorules_rr02_rp2_flat_retro.csv /opt/retropath/rules/rr02.csv

ENV RETROPATH_RULES_PATH=/opt/retropath/rules/rr02.csv \
    RETROPATH_RULES_VERSION=rr02-rp2-hs \
    RETROPATH_RULES_SHA256=e24eb97d3172195d03abed6e7da07a4cfd53965553853d126aaa8a93b4bc552f \
    RETROPATH_KNIME_DIR=/opt/knime \
    RETROPATH_WORKFLOW_VERSION=r20260212 \
    RETROPATH_KNIME_VERSION=4.7.0 \
    RETROPATH_RDKIT_PLUGIN_VERSION=4.9.1-glade-memoryfix2 \
    RETROPATH_OPENSSL10_DIR=/home/mambauser/.openssl10/lib \
    RETROPATH_JOBS_DIR=/var/lib/retropath/jobs \
    RETROPATH_MAX_QUEUE=8 \
    RETROPATH_JOB_TIMEOUT_SECONDS=1800 \
    RETROPATH_MEMORY_LIMIT_BYTES=8053063680 \
    RETROPATH_MEMORY_LIMIT_CONSECUTIVE_SAMPLES=3 \
    RETROPATH_RESOURCE_POLL_SECONDS=2 \
    RETROPATH_CGROUP_MEMORY_STAT_PATH=/sys/fs/cgroup/memory.stat \
    PYTHONPATH=/opt/service \
    PYTHONUNBUFFERED=1 \
    LD_LIBRARY_PATH=/home/mambauser/.openssl10/lib:/opt/conda/lib \
    _JAVA_OPTIONS="-Xmx1024m -XX:MaxDirectMemorySize=256m" \
    JAVA_TOOL_OPTIONS="-XX:ActiveProcessorCount=2 -Dorg.knime.core.maxThreads=2 -javaagent:/opt/retropath/memoryfix/native-memory-agent.jar" \
    OMP_NUM_THREADS=1 \
    OPENBLAS_NUM_THREADS=1

RUN /opt/conda/bin/python /opt/retropath/memoryfix/verify.py
HEALTHCHECK --interval=30s --timeout=10s --start-period=120s --retries=5 CMD ["/opt/conda/bin/python", "-m", "app.healthcheck"]

EXPOSE 8000

CMD ["uvicorn", "app.main:app", "--host", "0.0.0.0", "--port", "8000", "--workers", "1"]
