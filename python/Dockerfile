# Container for the Python MVT analysis pipeline.
#
# The image bundles only the code and its dependencies; the data and results
# live on the host and are mounted at run time, so nothing large is baked in and
# outputs land directly on the host.
#
# Build (from this python/ directory):
#   docker build -t mvtpy .
#
# Run a stage, mounting the workspace's data (read-only) and results (writable):
#   docker run --rm \
#     -v /abs/path/to/data:/work/data:ro \
#     -v /abs/path/to/results:/work/results \
#     mvtpy samples --day 16
#
# The default MVT_DATA_DIR / MVT_RESULTS_DIR point at those mount points, so no
# --data-dir/--results-dir flags are needed. See docker-compose.yml for a
# ready-made service, and README.md for the full walkthrough.

FROM python:3.12-slim AS base

# h5py/scipy/matplotlib wheels are self-contained on manylinux, so no system
# build tools are needed for a standard install.
ENV PYTHONUNBUFFERED=1 \
    MVT_DATA_DIR=/work/data \
    MVT_RESULTS_DIR=/work/results

WORKDIR /opt/mvtpy

# Install dependencies first (better layer caching), then the package.
COPY pyproject.toml README.md ./
COPY mvtpy ./mvtpy
RUN pip install --no-cache-dir ".[full]"

# The workspace mount points; overridable, but these are the defaults the
# environment variables above expect.
WORKDIR /work
VOLUME ["/work/data", "/work/results"]

ENTRYPOINT ["mvt"]
CMD ["--help"]
