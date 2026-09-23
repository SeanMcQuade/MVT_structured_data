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
#   make watch               # live progress of a run, from a second terminal
#   make clean-figures-17    # narrow, explicit cleanup
#
# (C) 2026 CIRCLES Consortium. BSD-3-Clause.

SHELL := /bin/bash
.DELETE_ON_ERROR:

# V=1 echoes the raw MATLAB command lines; by default recipes print a short
# labelled header instead, so `make -j` output stays readable.
V ?= 0
Q := $(if $(filter 1,$(V)),,@)

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

# Resolve DATA/RESULTS to absolute paths against the directory make was invoked
# from. Recipes cd into Scripts/ or python/ before running, and MATLAB and the
# Python CLI resolve relative paths against different bases again, so a relative
# override otherwise lands somewhere the caller did not choose.
# `override` is required: a variable set on the command line normally wins over
# every assignment in the makefile, so a plain := here would be ignored for
# exactly the case that needs it.
override DATA    := $(abspath $(DATA))
override RESULTS := $(abspath $(RESULTS))
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
SRC_lanes   := $(SCRIPTS)/generate_orig_dist_lanes.m
SRC_lc      := $(SCRIPTS)/extract_lane_changes_v_dist_to_av.m
SRC_relspeed := $(SCRIPTS)/relative_speed_histogram.m
SRC_lcplot  := $(SCRIPTS)/plotting_LC_analysis.m
SRC_relspeedplot := $(SCRIPTS)/plot_relative_speed.m $(SCRIPTS)/binned_relative_speed.m
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
$(Q)printf '[make] %-8s 2022-11-%s  ->  %s\n' "$(1)" "$(2)" "$(RESULTS)"
$(Q)cd $(SCRIPTS) && $(MATLAB) $(MATLAB_FLAGS) "mvt.build('$(1)', $(2))"
endef

# $(call run_sharded,STAGE,DAY) - SHARDS MATLAB processes over the 24 segments
define run_sharded
@if [ "$(SHARDS)" -le 1 ]; then \
  cd $(SCRIPTS) && $(MATLAB) $(MATLAB_FLAGS) "mvt.build('$(1)', $(2))"; \
else \
  printf '[make] %-8s 2022-11-%s  ->  %s  (%s shards)\n' "$(1)" "$(2)" "$(RESULTS)" "$(SHARDS)"; \
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
.PHONY: all data figures figures-from-mat figures-from-slim \
        gps slim full lanes lc lcplot relspeed relspeedplot samples fields \
        macro micro av \
        status config test help migrate-cache migrate-lanes accept watch \
        watch-once verify

# `make` does the figures that need only the small .mat intermediates first, so
# that a results-only download produces output before anything reaches for the
# 51 GB slim tree. See docs/DOWNLOADS.md.
all: figures-from-mat figures-from-slim

# `lanes` and `lc` are part of the data set but not yet of `all`: the figures
# that consume LC_data (plotting_LC_analysis) are not wired in yet, so nothing
# downstream of them would be built.
data: gps slim lanes lc relspeed samples fields

figures: figures-from-mat figures-from-slim

# Buildable from results/gps plus the .mat intermediates (fields_*, samples_*,
# LC_data_*) - no slim tree required.
figures-from-mat: macro lcplot av relspeedplot

# These read the slim JSON, or (for micro) the reduced plotting caches derived
# from it, so they need one of the larger downloads.
figures-from-slim: micro

gps:     $(foreach d,$(DAYS),$(STAMPS)/gps-$(d))
slim:    $(foreach d,$(DAYS),$(STAMPS)/slim-$(d))
lanes:   $(foreach d,$(DAYS),$(STAMPS)/lanes-$(d))
lc:      $(foreach d,$(DAYS),$(STAMPS)/lc-$(d))
relspeed: $(foreach d,$(DAYS),$(STAMPS)/relspeed-$(d))
relspeedplot: $(foreach d,$(DAYS),$(STAMPS)/relspeedplot-$(d))
lcplot:  $(foreach d,$(DAYS),$(STAMPS)/lcplot-$(d))
full:    $(foreach d,$(DAYS),$(STAMPS)/full-$(d))
samples: $(foreach d,$(DAYS),$(STAMPS)/samples-$(d))
fields:  $(foreach d,$(DAYS),$(STAMPS)/fields-$(d))
macro:   $(foreach d,$(DAYS),$(STAMPS)/macro-$(d))
micro:   $(foreach d,$(DAYS),$(STAMPS)/micro-$(d))
av:      $(STAMPS)/av

# Convenience aliases so `make slim-17` works as well as `make $(STAMPS)/slim-17`
define day_aliases
.PHONY: gps-$(1) slim-$(1) full-$(1) lanes-$(1) lc-$(1) lcplot-$(1) relspeed-$(1) relspeedplot-$(1) samples-$(1) fields-$(1) macro-$(1) micro-$(1) day-$(1)
gps-$(1):     $(STAMPS)/gps-$(1)
slim-$(1):    $(STAMPS)/slim-$(1)
lanes-$(1):   $(STAMPS)/lanes-$(1)
lc-$(1):      $(STAMPS)/lc-$(1)
relspeed-$(1): $(STAMPS)/relspeed-$(1)
relspeedplot-$(1): $(STAMPS)/relspeedplot-$(1)
lcplot-$(1):  $(STAMPS)/lcplot-$(1)
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
# The derived inputs, as files
#
# The figure stages depend on these .mat files directly, not on the stamps of
# the stages that make them. A results-only download already has them, and
# depending on the stamps would rebuild gps and slim to recreate files that are
# already on disk - which is exactly what `make figures-from-mat` used to do in
# a fresh tree. If a file really is missing, the recipe below builds it; the
# stamp is removed first so the stage runs even when a stale stamp claims it is
# already done.
# ---------------------------------------------------------------------------
define derived_inputs
$$(RESULTS)/analysis/2022-11-$(1)/fields_motion_2022-11-$(1).mat:
	@rm -f $$(STAMPS)/fields-$(1)
	$$(MAKE) $$(STAMPS)/fields-$(1)

$$(RESULTS)/analysis/2022-11-$(1)/samples_for_distance_analysis_$(1).mat:
	@rm -f $$(STAMPS)/samples-$(1)
	$$(MAKE) $$(STAMPS)/samples-$(1)

$$(RESULTS)/analysis/2022-11-$(1)/LC_data_$(1).mat:
	@rm -f $$(STAMPS)/lc-$(1)
	$$(MAKE) $$(STAMPS)/lc-$(1)

$$(RESULTS)/analysis/2022-11-$(1)/relspeed_data_$(1).mat:
	@rm -f $$(STAMPS)/relspeed-$(1)
	$$(MAKE) $$(STAMPS)/relspeed-$(1)

$$(RESULTS)/gps/CIRCLES_GPS_10Hz_2022-11-$(1).json:
	@rm -f $$(STAMPS)/gps-$(1)
	$$(MAKE) $$(STAMPS)/gps-$(1)
endef
$(foreach d,16 17 18,$(eval $(call derived_inputs,$(d))))

# What a stage that reads another stage's output should depend on: the files
# themselves when they are already there, and the producing stage's stamp only
# when they are not. A results-only download has the files but can never run the
# producing stage - slim and lanes both read the raw data - so depending on the
# stamp unconditionally made those trees unbuildable.
slim_files  = $(wildcard $(RESULTS)/slim/2022-11-$(1)/I-24MOTION_*.json)
lane_files  = $(wildcard $(RESULTS)/analysis/2022-11-$(1)/I-24MOTION_*_orig_dist_lane.mat)
slim_dep    = $(if $(call slim_files,$(1)),$(call slim_files,$(1)),$(STAMPS)/slim-$(1))
lanes_dep   = $(if $(call lane_files,$(1)),$(call lane_files,$(1)),$(STAMPS)/lanes-$(1))

# Shorthands for the per-day derived inputs
fields_mat  = $(RESULTS)/analysis/2022-11-$(1)/fields_motion_2022-11-$(1).mat
samples_mat = $(RESULTS)/analysis/2022-11-$(1)/samples_for_distance_analysis_$(1).mat
lc_mat      = $(RESULTS)/analysis/2022-11-$(1)/LC_data_$(1).mat
relspeed_mat= $(RESULTS)/analysis/2022-11-$(1)/relspeed_data_$(1).mat
gps_json    = $(RESULTS)/gps/CIRCLES_GPS_10Hz_2022-11-$(1).json

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

# lanes reads the raw segments directly, so it does not wait for slim; lc pairs
# the two by index and needs both.
$(STAMPS)/lanes-%: $(SRC_lanes) | $(STAMPS)
	$(call run_sharded,lanes,$*)
	@touch $@

















# Stages that read the slim trajectories, and lc which also reads the lane
# sidecars. Written per day so each can depend on the files it actually needs.
define slim_reader_rules
$$(STAMPS)/lc-$(1): $$(SRC_lc) $(call slim_dep,$(1)) $(call lanes_dep,$(1)) | $$(STAMPS)
	$$(call run_stage,lc,$(1))
	@touch $$@

$$(STAMPS)/relspeed-$(1): $$(SRC_relspeed) $(call slim_dep,$(1)) | $$(STAMPS)
	$$(call run_stage,relspeed,$(1))
	@touch $$@

$$(STAMPS)/samples-$(1): $$(SRC_samples) $(call slim_dep,$(1)) | $$(STAMPS)
	$$(call run_stage,samples,$(1))
	@touch $$@

$$(STAMPS)/fields-$(1): $$(SRC_fields) $(call slim_dep,$(1)) | $$(STAMPS)
	$$(call run_stage,fields,$(1))
	@touch $$@

$$(STAMPS)/micro-$(1): $$(SRC_micro) $(call slim_dep,$(1)) | $$(STAMPS)
	$$(call run_stage,micro,$(1))
	@touch $$@
endef
$(foreach d,16 17 18,$(eval $(call slim_reader_rules,$(d))))

# Figures built from the derived .mat files. These depend on the files, so a
# results-only tree builds them without touching gps or slim.
define figure_from_mat_rules
$$(STAMPS)/macro-$(1): $$(SRC_macro) $(call fields_mat,$(1)) $(call gps_json,$(1)) | $$(STAMPS)
	$$(call run_stage,macro,$(1))
	@touch $$@

$$(STAMPS)/lcplot-$(1): $$(SRC_lcplot) $(call lc_mat,$(1)) $(call gps_json,$(1)) | $$(STAMPS)
	$$(call run_stage,lcplot,$(1))
	@touch $$@

$$(STAMPS)/relspeedplot-$(1): $$(SRC_relspeedplot) $(call relspeed_mat,$(1)) | $$(STAMPS)
	$$(call run_stage,relspeedplot,$(1))
	@touch $$@
endef
$(foreach d,16 17 18,$(eval $(call figure_from_mat_rules,$(d))))

# Cross-day: needs every requested day's samples first.
$(STAMPS)/av: $(SRC_av) $(foreach d,$(DAYS),$(call samples_mat,$(d))) | $(STAMPS)
	cd $(SCRIPTS) && MVT_DAYS="$(DAYS)" $(MATLAB) $(MATLAB_FLAGS) "mvt.build('av', [])"
	@touch $@

# ---------------------------------------------------------------------------
# Inspection
# ---------------------------------------------------------------------------
status:
	@cd $(SCRIPTS) && MVT_DAYS="$(DAYS)" $(MATLAB) $(MATLAB_FLAGS) "mvt.status()"

# ---------------------------------------------------------------------------
# Output verification
#
# Checks generated files against the expected md5s in python/expected/. The
# manifest is derived from the MATLAB outputs, so a pass is evidence that a
# pipeline reproduced them byte for byte - which is the claim worth making
# about the Python port. PNGs are excluded: renderers do not agree pixel for
# pixel, so a checksum there would fail for reasons unrelated to correctness.
#
#   make verify                       # every day in DAYS
#   make verify RESULTS=/tmp/other    # verify a tree built somewhere else
#   make verify-against RESULTS=/tmp/other   # explain how it differs
#   make verify-update                # re-record expected checksums
# ---------------------------------------------------------------------------
PYTHON ?= python3
.PHONY: verify verify-update verify-against
verify:
	@cd $(REPO_ROOT)/python && $(PYTHON) -m mvtpy verify \
	  $(foreach d,$(DAYS),--day $(d)) --results-dir "$(RESULTS)"

verify-against:
	@cd $(REPO_ROOT)/python && $(PYTHON) -m mvtpy verify \
	  $(foreach d,$(DAYS),--day $(d)) --results-dir "$(RESULTS)" \
	  --reference "$(WORKSPACE)/results"

verify-update:
	@cd $(REPO_ROOT)/python && $(PYTHON) -m mvtpy verify \
	  $(foreach d,$(DAYS),--day $(d)) --results-dir "$(RESULTS)" --update

# Live progress of a run in flight. Start this in a second terminal while
# `make -jN` works: it reads the outputs on disk, so it needs no cooperation
# from the stages and adds no dependency that would make them stale.
# `make watch-once` prints a single snapshot, which is what you want in a log.
.PHONY: watch watch-once
watch:
	@cd $(SCRIPTS) && MVT_DAYS="$(DAYS)" $(MATLAB) $(MATLAB_FLAGS) "mvt.watch()"

watch-once:
	@cd $(SCRIPTS) && MVT_DAYS="$(DAYS)" $(MATLAB) $(MATLAB_FLAGS) \
	  "mvt.watch('Once', true)"

# Check the data are laid out where the pipeline expects them.
.PHONY: check-data
check-data:
	@MVT_DATA_DIR="$(DATA)" MVT_RESULTS_DIR="$(RESULTS)" $(REPO_ROOT)/check_data.sh

# Rebuild everything from the raw inputs, ignoring what is already on disk.
# FORCE=1 makes every stage's own staleness check report stale, so this does not
# depend on timestamps being meaningful - which is the point after a `git pull`
# rewrites every .m mtime.
#
# Intended for a from-scratch run into a fresh tree, which leaves the released
# results untouched:
#   make DATA=/path/to/data RESULTS=/path/to/new_results -j3 rebuild
#
# Covers the released products (gps, slim, lanes, lc, samples, fields) and the
# figures. `full` is opt-in as always: add `make ... FORCE=1 full` for it.
.PHONY: rebuild
rebuild:
	$(MAKE) FORCE=1 data figures av

# Accept existing outputs as current, but only once their checksums prove they
# are the expected bytes: md5 decides, not the clock. This is the safe answer
# when `make status` reports everything stale after a pull.
#
# Scope: `verify` covers the gps and slim JSON only - the .mat products carry
# gzip creation timestamps and the .png figures are renderer-dependent, so
# neither is checksummed (see python/mvtpy/verify.py). Those stages are still
# accepted here; only the JSON is actually proven.
.PHONY: accept-verified
accept-verified:
	@$(MAKE) verify
	@$(MAKE) accept

# Mark existing outputs as current, for code changes that provably do not alter
# results. Verify first by rebuilding into a separate tree and comparing:
#   make RESULTS=$(WORKSPACE)/results_verify slim-16 && md5 <old> <new>
accept:
	@cd $(SCRIPTS) && MVT_DAYS="$(DAYS)" $(MATLAB) $(MATLAB_FLAGS) \
	  "files = mvt.accept(); fprintf('accepted %d files\n', numel(files));"
	@mkdir -p $(STAMPS)
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

# Fast suite: temporary files only, no data tree needed. Seconds to run.
test:
	@cd $(SCRIPTS) && $(MATLAB) $(MATLAB_FLAGS) \
	  "r = runtests('$(REPO_ROOT)/tests'); \
	   fprintf('\n%d passed, %d failed (%.1f s)\n', sum([r.Passed]), sum([r.Failed]), sum([r.Duration])); \
	   if any([r.Failed]), exit(1); end"

# Full check: compare real outputs against the recorded manifests.
.PHONY: verify-full manifests
verify-full:
	@for d in $(DAYS); do \
	  cd $(SCRIPTS) && $(MATLAB) $(MATLAB_FLAGS) \
	    "addpath('$(REPO_ROOT)/tests'); verify_outputs($$d)" || exit 1; \
	done

# Record manifests from a known-good tree, e.g.
#   make RESULTS=$(WORKSPACE)/results_groundtruth manifests
manifests:
	@for d in $(DAYS); do \
	  cd $(SCRIPTS) && $(MATLAB) $(MATLAB_FLAGS) \
	    "addpath('$(REPO_ROOT)/tests'); generate_manifest($$d)" || exit 1; \
	done

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
.PHONY: clean-figures-$(1) clean-slim-$(1) clean-full-$(1) clean-lanes-$(1) clean-lc-$(1) clean-lcplot-$(1) clean-relspeed-$(1) clean-relspeedplot-$(1)
clean-figures-$(1):
	rm -f $(RESULTS)/figures/2022-11-$(1)/fig_*.png $(RESULTS)/figures/2022-11-$(1)/fig_*.fig
	rm -f $(STAMPS)/macro-$(1) $(STAMPS)/micro-$(1)

clean-slim-$(1):
	rm -f $(RESULTS)/slim/2022-11-$(1)/I-24MOTION_*.json
	rm -f $(STAMPS)/slim-$(1)

clean-full-$(1):
	rm -f $(RESULTS)/full/2022-11-$(1)/I-24MOTION_*.json
	rm -f $(STAMPS)/full-$(1)
clean-lanes-$(1):
	rm -f $(RESULTS)/analysis/2022-11-$(1)/I-24MOTION_*_orig_dist_lane.mat
	rm -f $(STAMPS)/lanes-$(1)
clean-lcplot-$(1):
	rm -f $(RESULTS)/figures/2022-11-$(1)/fig_lc_*_202211$(1).png
	rm -f $(STAMPS)/lcplot-$(1)
clean-lc-$(1):
	rm -f $(RESULTS)/analysis/2022-11-$(1)/LC_data_$(1).mat
	rm -f $(STAMPS)/lc-$(1)
clean-relspeed-$(1):
	rm -f $(RESULTS)/analysis/2022-11-$(1)/relspeed_data_$(1).mat
	rm -f $(STAMPS)/relspeed-$(1)
clean-relspeedplot-$(1):
	rm -f "$(RESULTS)/figures/Relative speed histogram, day $(1) for files j = "*.pdf
	rm -f "$(RESULTS)/figures/Relative speeds behind AV, day $(1).pdf"
	rm -f $(STAMPS)/relspeedplot-$(1)
endef
$(foreach d,16 17 18,$(eval $(call clean_aliases,$(d))))

# One-time migration: the lane-change work first wrote its .mat products into
# results/slim/<day>/ alongside the released JSON. Move them to the figures
# folder, where the other analysis .mat files live, so results/slim contains
# only published artifacts. Never clobbers: a file already rebuilt at the
# destination is left alone and the superseded copy stays in slim/ for you to
# delete deliberately, because the pre-move sidecars carry the accumulator bug
# (their length is a running maximum of the segment trajectory counts).
migrate-lanes:
	@for d in $(DAYS); do \
	  src="$(RESULTS)/slim/2022-11-$$d"; dst="$(RESULTS)/figures/2022-11-$$d"; \
	  mkdir -p "$$dst"; \
	  moved=0; kept=0; \
	  for f in "$$src"/*_orig_dist_lane.mat "$$src/LC_data_$$d.mat"; do \
	    [ -f "$$f" ] || continue; \
	    if [ -e "$$dst/$$(basename "$$f")" ]; then \
	      kept=$$((kept+1)); \
	    else \
	      mv "$$f" "$$dst"/; moved=$$((moved+1)); \
	    fi; \
	  done; \
	  echo "[make] 2022-11-$$d: moved $$moved file(s) to $$dst; left $$kept in place (already rebuilt there)"; \
	done

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
