# Makefile for the CIRCLES MegaVanderTest pipeline.
#
# Encodes the stage dependency graph so that `make` rebuilds only what is out
# of date and `make -j` runs independent work in parallel. Each recipe shells
# out to MATLAB in batch mode; the MATLAB side performs its own, finer-grained
# staleness check (see Scripts/+mvt/isStale.m), so a target that make decides
# to visit may still legitimately do nothing.
#
# Layout assumption (unchanged from the original pipeline): this repository is
# a sibling of the data/ and results/ folders.
#
#   make                     # everything, all days (same as `make all`)
#   make -j3 slim            # slim data for all days, three MATLAB processes
#   make SHARDS=4 slim-17    # one day, four MATLAB processes over 24 segments
#   make FORCE=1 fields-17   # rebuild even if outputs look fresh
#   make -n all              # dry run: show the plan without running MATLAB
#   make status              # ask MATLAB what is stale and why
#   make clean-figures-17    # narrow, explicit cleanup
#
# (C) 2026 CIRCLES Consortium. BSD-3-Clause.

SHELL := /bin/bash
.DELETE_ON_ERROR:

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
MATLAB       ?= /Applications/MATLAB_R2025b.app/bin/matlab
MATLAB_FLAGS ?= -nodisplay -nosplash -batch
DAYS         ?= 16 17 18
SHARDS       ?= 1
FORCE        ?= 0
CLEAN        ?= 0
VERBOSE      ?= 1
SETTLE       ?=

REPO_ROOT := $(patsubst %/,%,$(dir $(abspath $(lastword $(MAKEFILE_LIST)))))
WORKSPACE := $(patsubst %/,%,$(dir $(REPO_ROOT)))
SCRIPTS   := $(REPO_ROOT)/Scripts
MODELS    := $(REPO_ROOT)/Models
DATA      ?= $(WORKSPACE)/data
# Override to write a second copy of the outputs for comparison, e.g.
#   make RESULTS=$(WORKSPACE)/results_verify slim-16
RESULTS   ?= $(WORKSPACE)/results
STATE     := $(RESULTS)/.mvt
STAMPS    := $(STATE)/stamps

# Options and locations are passed to MATLAB through the environment
# (see mvt.options and mvt.paths).
export MVT_FORCE       := $(FORCE)
export MVT_CLEAN       := $(CLEAN)
export MVT_VERBOSE     := $(VERBOSE)
export MVT_DATA_DIR    := $(DATA)
export MVT_RESULTS_DIR := $(RESULTS)
ifneq ($(strip $(SETTLE)),)
export MVT_SETTLE_SECONDS := $(SETTLE)
endif

# ---------------------------------------------------------------------------
# Per-stage source dependencies (editing these invalidates the stage)
# ---------------------------------------------------------------------------
FUEL_MODELS := $(wildcard $(MODELS)/*.m) $(MODELS)/Eastbound_grade_fit.csv
SRC_gps     := $(SCRIPTS)/assemble_data_GPS.m
SRC_slim    := $(SCRIPTS)/generate_data_mvt_slim.m $(FUEL_MODELS)
SRC_full    := $(SCRIPTS)/generate_data_mvt_full.m $(FUEL_MODELS)
SRC_samples := $(SCRIPTS)/generate_data_samples.m
SRC_fields  := $(SCRIPTS)/generate_macroscopic_fields.m
SRC_macro   := $(SCRIPTS)/plot_macroscopic_fields.m
SRC_micro   := $(SCRIPTS)/plot_microscopic_trajectories.m
SRC_av      := $(SCRIPTS)/plot_AV_analysis.m

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
# $(call raw_files,DAY) - raw MOTION segments for a day
raw_files = $(wildcard $(DATA)/i24motion/2022-11-$(1)/*_0_*.json)
# $(call car_files,DAY) - control-vehicle inputs for a day
car_files = $(wildcard $(DATA)/cars/cars_gps/circles_v2_1_car*.csv) \
            $(DATA)/cars/cars_vins.csv $(DATA)/cars/veh_ping_202211$(1).csv

# $(call run_stage,STAGE,DAY) - one MATLAB process for a stage
define run_stage
cd $(SCRIPTS) && $(MATLAB) $(MATLAB_FLAGS) "mvt.build('$(1)', $(2))"
endef

# $(call run_sharded,STAGE,DAY) - SHARDS MATLAB processes over the 24 segments
define run_sharded
@if [ "$(SHARDS)" -le 1 ]; then \
  cd $(SCRIPTS) && $(MATLAB) $(MATLAB_FLAGS) "mvt.build('$(1)', $(2))"; \
else \
  echo "[make] $(1) 2022-11-$(2): $(SHARDS) shards"; \
  pids=""; \
  for k in $$(seq 1 $(SHARDS)); do \
    ( cd $(SCRIPTS) && MVT_SHARD=$$k/$(SHARDS) \
      $(MATLAB) $(MATLAB_FLAGS) "mvt.build('$(1)', $(2))" ) & \
    pids="$$pids $$!"; \
  done; \
  status=0; \
  for pid in $$pids; do wait $$pid || status=1; done; \
  exit $$status; \
fi
endef

# ---------------------------------------------------------------------------
# Aggregate targets
# ---------------------------------------------------------------------------
.PHONY: all data figures gps slim full samples fields macro micro av \
        status config test help migrate-cache accept

all: figures av

data: gps slim samples fields

figures: macro micro

gps:     $(foreach d,$(DAYS),$(STAMPS)/gps-$(d))
slim:    $(foreach d,$(DAYS),$(STAMPS)/slim-$(d))
full:    $(foreach d,$(DAYS),$(STAMPS)/full-$(d))
samples: $(foreach d,$(DAYS),$(STAMPS)/samples-$(d))
fields:  $(foreach d,$(DAYS),$(STAMPS)/fields-$(d))
macro:   $(foreach d,$(DAYS),$(STAMPS)/macro-$(d))
micro:   $(foreach d,$(DAYS),$(STAMPS)/micro-$(d))
av:      $(STAMPS)/av

# Convenience aliases so `make slim-17` works as well as `make $(STAMPS)/slim-17`
define day_aliases
.PHONY: gps-$(1) slim-$(1) full-$(1) samples-$(1) fields-$(1) macro-$(1) micro-$(1) day-$(1)
gps-$(1):     $(STAMPS)/gps-$(1)
slim-$(1):    $(STAMPS)/slim-$(1)
full-$(1):    $(STAMPS)/full-$(1)
samples-$(1): $(STAMPS)/samples-$(1)
fields-$(1):  $(STAMPS)/fields-$(1)
macro-$(1):   $(STAMPS)/macro-$(1)
micro-$(1):   $(STAMPS)/micro-$(1)
day-$(1):     $(STAMPS)/macro-$(1) $(STAMPS)/micro-$(1)
endef
$(foreach d,16 17 18,$(eval $(call day_aliases,$(d))))

$(STAMPS):
	@mkdir -p $(STAMPS)

# ---------------------------------------------------------------------------
# Stage rules
#
# Each stamp records "make last ran this stage successfully". The real output
# freshness test lives in MATLAB, which is what makes a rebuild after a code
# change correct even if the stamp is missing or stale.
# ---------------------------------------------------------------------------
$(STAMPS)/gps-%: $(SRC_gps) | $(STAMPS)
	$(call run_stage,gps,$*)
	@touch $@

$(STAMPS)/slim-%: $(SRC_slim) $(STAMPS)/gps-% | $(STAMPS)
	$(call run_sharded,slim,$*)
	@touch $@

$(STAMPS)/full-%: $(SRC_full) $(STAMPS)/gps-% | $(STAMPS)
	$(call run_sharded,full,$*)
	@touch $@

$(STAMPS)/samples-%: $(SRC_samples) $(STAMPS)/slim-% | $(STAMPS)
	$(call run_stage,samples,$*)
	@touch $@

$(STAMPS)/fields-%: $(SRC_fields) $(STAMPS)/slim-% | $(STAMPS)
	$(call run_stage,fields,$*)
	@touch $@

$(STAMPS)/macro-%: $(SRC_macro) $(STAMPS)/fields-% | $(STAMPS)
	$(call run_stage,macro,$*)
	@touch $@

$(STAMPS)/micro-%: $(SRC_micro) $(STAMPS)/slim-% | $(STAMPS)
	$(call run_stage,micro,$*)
	@touch $@

# Cross-day: needs every requested day's samples first.
$(STAMPS)/av: $(SRC_av) $(foreach d,$(DAYS),$(STAMPS)/samples-$(d)) | $(STAMPS)
	cd $(SCRIPTS) && MVT_DAYS="$(DAYS)" $(MATLAB) $(MATLAB_FLAGS) "mvt.build('av', [])"
	@touch $@

# ---------------------------------------------------------------------------
# Inspection
# ---------------------------------------------------------------------------
status:
	@cd $(SCRIPTS) && MVT_DAYS="$(DAYS)" $(MATLAB) $(MATLAB_FLAGS) "mvt.status()"

# Mark existing outputs as current, for code changes that provably do not alter
# results. Verify first by rebuilding into a separate tree and comparing:
#   make RESULTS=$(WORKSPACE)/results_verify slim-16 && md5 <old> <new>
accept:
	@cd $(SCRIPTS) && MVT_DAYS="$(DAYS)" $(MATLAB) $(MATLAB_FLAGS) \
	  "files = mvt.accept(); fprintf('accepted %d files\n', numel(files));"
	@mkdir -p $(STAMPS)
	@for d in $(DAYS); do \
	  for s in gps slim full samples fields macro micro; do touch $(STAMPS)/$$s-$$d; done; \
	done
	@echo "[make] stamps refreshed; 'make' will not revisit accepted stages"

config:
	@echo "MATLAB    = $(MATLAB)"
	@echo "REPO_ROOT = $(REPO_ROOT)"
	@echo "WORKSPACE = $(WORKSPACE)"
	@echo "DATA      = $(DATA)"
	@echo "RESULTS   = $(RESULTS)"
	@echo "DAYS      = $(DAYS)"
	@echo "SHARDS    = $(SHARDS)"
	@echo "FORCE     = $(FORCE)"

test:
	@cd $(SCRIPTS) && $(MATLAB) $(MATLAB_FLAGS) \
	  "results = runtests('$(REPO_ROOT)/tests'); disp(results); exit(any([results.Failed]))"

help:
	@sed -n '1,20p' $(lastword $(MAKEFILE_LIST))

# ---------------------------------------------------------------------------
# Cleaning - deliberately narrow. results/ holds ~150 GB; nothing here removes
# a whole tree, and there is no `make clean`.
# ---------------------------------------------------------------------------
.PHONY: clean-stamps clean-cache
clean-stamps:
	rm -f $(STAMPS)/*

clean-cache:
	rm -rf $(STATE)/cache

define clean_aliases
.PHONY: clean-figures-$(1) clean-slim-$(1) clean-full-$(1)
clean-figures-$(1):
	rm -f $(RESULTS)/figures/2022-11-$(1)/fig_*.png $(RESULTS)/figures/2022-11-$(1)/fig_*.fig
	rm -f $(STAMPS)/macro-$(1) $(STAMPS)/micro-$(1)

clean-slim-$(1):
	rm -f $(RESULTS)/slim/2022-11-$(1)/I-24MOTION_*.json
	rm -f $(STAMPS)/slim-$(1)

clean-full-$(1):
	rm -f $(RESULTS)/full/2022-11-$(1)/I-24MOTION_*.json
	rm -f $(STAMPS)/full-$(1)
endef
$(foreach d,16 17 18,$(eval $(call clean_aliases,$(d))))

# One-time migration: earlier versions wrote the *_reduced.mat plotting caches
# into results/slim/<day>/ alongside the released JSON. Move them to the cache
# folder so results/slim contains only published artifacts.
migrate-cache:
	@for d in $(DAYS); do \
	  src="$(RESULTS)/slim/2022-11-$$d"; dst="$(STATE)/cache/2022-11-$$d"; \
	  if compgen -G "$$src/*_reduced.mat" > /dev/null; then \
	    mkdir -p "$$dst"; \
	    echo "[make] moving $$(ls $$src/*_reduced.mat | wc -l | tr -d ' ') cache files for 2022-11-$$d"; \
	    mv "$$src"/*_reduced.mat "$$dst"/; \
	  else \
	    echo "[make] no cache files to move for 2022-11-$$d"; \
	  fi; \
	done
