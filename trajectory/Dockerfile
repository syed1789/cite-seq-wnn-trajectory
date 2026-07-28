FROM mambaorg/micromamba:2.0.5
ARG MAMBA_DOCKERFILE_ACTIVATE=1
COPY --chown=$MAMBA_USER:$MAMBA_USER environment.yml /tmp/environment.yml
RUN micromamba install -y -n base -f /tmp/environment.yml && micromamba clean --all --yes
# bioconda's r-monocle3 recipe is pinned to a broken/unsolvable old R 3.6
# dependency chain (r-spelling -> r-hunspell version conflict) -- installing
# from source instead, same pattern as the CellChat install in my
# heart-nextflow-r/Dockerfile.
# Pinned to v1.3.1, not master/HEAD: newer monocle3 (>= 1.4) added BPCells,
# ggdist, and sf as hard Imports, none of which build in this minimal image
# (BPCells especially -- it's fetched from a separate GitHub fork with its
# own native build requirements). v1.3.1 is the last tag before that, its
# full Imports: list is satisfied entirely by the conda packages above
# (confirmed against the actual DESCRIPTION file, not guessed), and it's
# still the version most tutorials/papers cite.
RUN R -q -e "remotes::install_github('cole-trapnell-lab/leidenbase', upgrade='never', dependencies=FALSE)"
RUN R -q -e "remotes::install_github('cole-trapnell-lab/monocle3', ref='v1.3.1', upgrade='never', dependencies=FALSE)"
ENV R_DEFAULT_DEVICE=pdf
WORKDIR /work
