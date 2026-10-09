# DRAFT -- not built yet (ROADMAP.md, "Images").
# The legacy toolchain of the published method plus s2g2p:
#   * Roary 3.8.2 from Seq2Geno's fork, hzi-bifo/Roary @ 7402d1b, whose Blastp.pm adds the
#     >95% identity / >95% coverage-of-both edge filter. Its install_dependencies.sh builds the
#     versions Seq2Geno ran: BLAST 2.4.0, CD-HIT 4.6.6, MCL 14-137, PRANK, bedtools 2.26, parallel.
#   * Perl from Seq2Geno's perl5_22_env pins (denovo/perl5_22_env.yml).
#   * mafft 7.310 and jvarkit-msa2vcf 1.0 (denovo/indel_cluster_env.yml, indel_env.yml).
# BLAST 2.4.0 needs libnsl.so.1; its absence is what silently disabled Roary's BLAST step in the
# 2022/23 runs behind Hu et al. 2024. The image must provide it, and the smoke test below fails
# if blastp cannot start.
FROM mambaorg/micromamba:1.5.10
USER root
RUN apt-get update -qq && apt-get install -y --no-install-recommends git make gcc g++ zlib1g-dev \
      libnsl2 wget bzip2 ca-certificates && rm -rf /var/lib/apt/lists/*
USER $MAMBA_USER
COPY --chown=$MAMBA_USER:$MAMBA_USER containers/tools/env-legacy.yml /tmp/env-legacy.yml
RUN micromamba create -y -n legacy -f /tmp/env-legacy.yml && micromamba clean -a -y
RUN git clone https://github.com/hzi-bifo/Roary.git /opt/roary \
 && cd /opt/roary && git checkout 7402d1bab868d73b9c1e7cf30e8af8c5b2578313 \
 && micromamba run -n legacy bash ./install_dependencies.sh
COPY --chown=$MAMBA_USER:$MAMBA_USER . /opt/s2g2p
RUN micromamba create -y -n s2g2p -c conda-forge python=3.12 numpy pandas scikit-learn pyyaml joblib pip \
 && micromamba run -n s2g2p pip install --no-deps /opt/s2g2p && micromamba clean -a -y
ENV ROARY_HOME=/opt/roary
ENV PATH=/opt/roary/bin:/opt/roary/build/fasttree:/opt/roary/build/mcl-14-137/src/alien/oxygen/src:/opt/roary/build/mcl-14-137/src/shmcl:/opt/roary/build/ncbi-blast-2.4.0+/bin:/opt/roary/build/prank-msa-master/src:/opt/roary/build/cd-hit-v4.6.6-2016-0711:/opt/roary/build/bedtools2/bin:/opt/roary/build/parallel-20160722/src:/opt/conda/envs/legacy/bin:/opt/conda/envs/s2g2p/bin:$PATH
ENV PERL5LIB=/opt/roary/lib:/opt/roary/build/bedtools2/lib
RUN blastp -version && mafft --version && roary -w && s2g2p --version
