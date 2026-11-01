FROM mambaorg/micromamba:2.3.2-ubuntu22.04@sha256:0e87302b8b802b595c947f408e02c436c1d49aa582132ddaa4cd5e1c991a4871

COPY --chown=$MAMBA_USER:$MAMBA_USER environment-linux-64.lock /tmp/environment-linux-64.lock
RUN micromamba install --yes --name base --file /tmp/environment-linux-64.lock \
    && micromamba clean --all --yes
COPY --chown=$MAMBA_USER:$MAMBA_USER requirements-http.txt /tmp/requirements-http.txt
COPY --chown=$MAMBA_USER:$MAMBA_USER verify_environment.py /tmp/verify_environment.py
RUN /opt/conda/bin/python -m pip install --no-cache-dir --no-deps --require-hashes -r /tmp/requirements-http.txt \
    && /opt/conda/bin/python /tmp/verify_environment.py

USER root
RUN mkdir -p /opt/service /var/lib/bactabolize/jobs \
    && chown -R "$MAMBA_USER:$MAMBA_USER" /opt/service /var/lib/bactabolize
WORKDIR /opt/service
COPY --chown=$MAMBA_USER:$MAMBA_USER app ./app
USER $MAMBA_USER

ENV PYTHONUNBUFFERED=1 \
    PYTHONDONTWRITEBYTECODE=1 \
    BACTABOLIZE_JOBS_DIR=/var/lib/bactabolize/jobs \
    BACTABOLIZE_MAX_QUEUE=8 \
    BACTABOLIZE_JOB_TIMEOUT_SECONDS=1800 \
    BACTABOLIZE_UPLOAD_LIMIT_BYTES=67108864 \
    BACTABOLIZE_SHUTDOWN_GRACE_SECONDS=30
EXPOSE 8000
HEALTHCHECK --interval=30s --timeout=10s --start-period=30s --retries=5 \
    CMD ["/opt/conda/bin/python", "-m", "app.healthcheck"]
CMD ["/opt/conda/bin/python", "-m", "uvicorn", "app.main:app", "--host", "0.0.0.0", "--port", "8000", "--workers", "1"]
