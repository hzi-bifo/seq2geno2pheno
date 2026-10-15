# s2g2p only: the Python steps (k-mers, GPA, compression, splits, evaluation).
# Build from the repository root:  docker build -f containers/s2g2p/Dockerfile -t seq2geno2pheno:dev .
FROM mambaorg/micromamba:1.5.10
COPY --chown=$MAMBA_USER:$MAMBA_USER . /opt/s2g2p
RUN micromamba install -y -n base -c conda-forge python=3.12 numpy pandas scikit-learn pyyaml joblib pip \
 && micromamba run -n base pip install --no-deps /opt/s2g2p \
 && micromamba clean -a -y
ENV PATH=/opt/conda/bin:$PATH
