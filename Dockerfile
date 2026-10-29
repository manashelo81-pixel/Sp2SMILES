# Web app container (CPU). Build:  docker build -t smiles-ms .
# Run:    docker run -p 8501:8501 -v /path/to/models:/app/models smiles-ms
FROM python:3.11-slim

RUN apt-get update && apt-get install -y --no-install-recommends libxrender1 libxext6 libsm6 libglib2.0-0 \
    && rm -rf /var/lib/apt/lists/*
WORKDIR /app
COPY webapp/requirements-web.txt webapp/requirements-web.txt
RUN pip install --no-cache-dir torch --index-url https://download.pytorch.org/whl/cpu \
    && pip install --no-cache-dir -r webapp/requirements-web.txt
COPY . /app

EXPOSE 8501
ENV SMILES_MS_MODEL=/app/models/model.pt SMILES_MS_DATASET=/app/models/dataset.pkl SMILES_MS_LOCK_PATHS=1
CMD ["python", "run_webapp.py", "--host", "0.0.0.0", "--no-browser", "--lock-paths"]
