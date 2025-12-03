FROM public.ecr.aws/docker/library/python:3.10-bookworm

WORKDIR /app

# Ensure a writable temp directory exists for pip and Python's tempfile
ENV TMPDIR=/tmp
RUN mkdir -p /tmp /var/tmp /usr/tmp && chmod 1777 /tmp /var/tmp /usr/tmp

# Configure apt
RUN echo "Acquire::http::Pipeline-Depth 0;" > /etc/apt/apt.conf.d/99custom && \
    echo "Acquire::http::No-Cache true;" >> /etc/apt/apt.conf.d/99custom && \
    echo "Acquire::BrokenProxy    true;" >> /etc/apt/apt.conf.d/99custom

# Harden apt: add retries/timeouts and switch sources to HTTPS (supports deb822 and legacy formats)
RUN set -eux; \
    printf 'Acquire::Retries "5";\nAcquire::http::Timeout "30";\nAcquire::https::Timeout "30";\n' > /etc/apt/apt.conf.d/80-retries; \
    if [ -f /etc/apt/sources.list ]; then \
      sed -i -e 's|http://|https://|g' /etc/apt/sources.list; \
    fi; \
    if ls /etc/apt/sources.list.d/* >/dev/null 2>&1; then \
      sed -i -e 's|http://|https://|g' /etc/apt/sources.list.d/*; \
    fi

# Base build tools
RUN apt-get update && DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
    build-essential cmake ffmpeg git curl \
    && rm -rf /var/lib/apt/lists/* && apt-get clean

# Install Python development dependencies
RUN apt-get update && DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
    python3-dev \
    python3-tk \
    && rm -rf /var/lib/apt/lists/* \
    && apt-get clean

# Install only necessary library dependencies
RUN apt-get update && DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
    libpq-dev \
    libffi-dev \
    zlib1g-dev \
    && rm -rf /var/lib/apt/lists/* \
    && apt-get clean

# Install only necessary GUI and rendering dependencies
RUN apt-get update && DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
    libgl1 \
    libxrender1 \
    && rm -rf /var/lib/apt/lists/* \
    && apt-get clean

# Install Playwright system dependencies
RUN apt-get update && DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
    libnss3 \
    libnspr4 \
    libatk1.0-0 \
    libatk-bridge2.0-0 \
    libcups2 \
    libdrm2 \
    libdbus-1-3 \
    libxkbcommon0 \
    libxcomposite1 \
    libxdamage1 \
    libxfixes3 \
    libxrandr2 \
    libgbm1 \
    libasound2 \
    && rm -rf /var/lib/apt/lists/* \
    && apt-get clean

# Set PYTHONPATH to include the parent directory of new_src
ENV PYTHONPATH="/app"

COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

# Copy the source code
COPY . /app/src

WORKDIR /app/src

EXPOSE 8000

CMD exec gunicorn main:app \
    --bind :8000 \
    --workers 2 \
    --timeout 300 \
    --worker-class uvicorn.workers.UvicornWorker \
    --threads 8 \
    --keep-alive 120 \
    --access-logfile - \
    --error-logfile - \
    --log-level debug
