FROM python:3.10 AS builder
RUN pip install -U --no-cache-dir pip setuptools wheel
RUN pip install -U --no-cache-dir pdm

COPY pyproject.toml README.rst pdm.lock /project/
COPY metabolike/ /project/metabolike
COPY tests/ /project/tests

WORKDIR /project
RUN pdm install --prod --no-lock --no-editable

FROM python:3.10

ENV PYTHONPATH=/project/pkgs
COPY --from=builder /project/__pypackages__/3.10/lib /project/pkgs

ENTRYPOINT ["python", "-m", "metabolike.main"]
CMD ["--help"]
