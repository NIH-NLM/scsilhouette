# Use slim Python base image
FROM python:3.10-slim

# Set non-root user (optional for security)
# RUN useradd -ms /bin/bash scuser
# USER scuser

# Set working directory
WORKDIR /app

# Copy only the necessary files first to leverage Docker cache
COPY pyproject.toml poetry.lock* ./

# Install Poetry
RUN pip install --no-cache-dir poetry

# Install dependencies
RUN poetry config virtualenvs.create false \
  && poetry install --no-interaction --no-ansi

# Copy source code
COPY src/ ./src/

# Install CLI entry point
ENV PYTHONPATH=/app/src
ENTRYPOINT ["scsilhouette"]
CMD ["--help"]

