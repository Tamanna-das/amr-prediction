# Run from the repo root with the `amr` env active.
ENVS := $(CONDA_PREFIX)/..
export PATH := $(ENVS)/amr-amrfinder/bin:$(ENVS)/amr-mash/bin:$(PATH)

.PHONY: all test dryrun

all:
	snakemake -s workflow/Snakefile --cores 6 report --config only_present=False

dryrun:
	snakemake -s workflow/Snakefile -n report --config only_present=False

test:
	pytest -q
