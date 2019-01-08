## Directory information
RAWDATDIR = ./data/raw-data
CLEANDATDIR = ./data/clean-data
CODEDIR = ./code
RESDIR = ./results
MSDIR = ./manuscript
SIMDIR = ./simulations


## Analysis information
NREPS = 50


## Set the manuscript as the overall target
.PHONY: all
all: $(MSDIR)/measles-ews-manuscript.pdf


## Clean the data
DATA_DEPS := $(RAWDATDIR)/district_pops.csv
DATA_DEPS += $(RAWDATDIR)/niger_crude_birth_rates.csv
DATA_DEPS += $(RAWDATDIR)/niger_regional_1995_2005.csv
DATA_DEPS += $(CODEDIR)/fetch-clean-data.R
DATA_TARG = $(CLEANDATDIR)/weekly-measles-incidence-niger-cities-clean.RDS
COVAR_TARG = $(CLEANDATDIR)/annual-demographic-data-niger-cities-clean.RDS

$(DATA_TARG): $(DATA_DEPS)
	Rscript $(CODEDIR)/fetch-clean-data.R

$(COVAR_TARG): $(DATA_DEPS)
	Rscript $(CODEDIR)/fetch-clean-data.R


## Make pomp objects
POMP_DEPS := $(DATA_TARG)
POMP_DEPS += $(COVAR_TARG)
POMP_DEPS += $(CODEDIR)/define-continuous-measles-pomp.R
POMP_TARG = $(CODEDIR)/measles-pomp-object-*.RDS

$(POMP_TARG): $(POMP_DEPS)
	Rscript $(CODEDIR)/define-continuous-measles-pomp.R
	$(warning Pomp objects updated. Do you need to refit on HPC?)


## Fit the POMP models
## NOTE: model fitting was done on a high performance computing cluster,
## using 1000s of cores. Therefore, this Makefile skips over the model
## fitting stage, but the R script global-search-mif.R can be interrogated
## by interested users.
.INTERMEDIATE = $(RESDIR)/initial-mif-lls-*.csv
#FIT_DEPS := $(POMP_TARG)
#FIT_DEPS += $(CODEDIR)/global-search-mif.R
#FIT_TARG := $(RESDIR)/initial-mif-lls-Agadez.csv
#FIT_TARG += $(RESDIR)/initial-mif-lls-Maradi.csv
#FIT_TARG += $(RESDIR)/initial-mif-lls-Niamey.csv
#FIT_TARG += $(RESDIR)/initial-mif-lls-Zinder.csv

#$(FIT_TARG): $(FIT_DEPS)
	#Rscript $(CODEDIR)/global-search-mif.R 1

LOGLIKS = $(RESDIR)/initial-mif-lls-*.csv

## Make one-week-ahead predictions
PRED_DEPS := $(LOGLIKS)
PRED_DEPS += $(CODEDIR)/estimate-transmission-state.R
PRED_TARG = $(RESDIR)/predictive-dist-states-*.RDS

$(PRED_TARG): $(PRED_DEPS)
	Rscript $(CODEDIR)/estimate-transmission-state.R 1
	Rscript $(CODEDIR)/estimate-transmission-state.R 2
	Rscript $(CODEDIR)/estimate-transmission-state.R 3
	Rscript $(CODEDIR)/estimate-transmission-state.R 4
	

## Simulate bootstrap replicates for parameter uncertainty
BOOT_DEPS := $(CODEDIR)/simulate-bootstrap-series.R
BOOT_DEPS += $(CODEDIR)/make-pomp-simulator-function.R
BOOT_DEPS += $(LOGLIKS)
BOOT_TARG = $(SIMDIR)/bootstrap-sims-*.RDS

$(BOOT_TARG): $(BOOT_DEPS)
	Rscript $(CODEDIR)/simulate-bootstrap-series.R


## Simulate emergence replicates
EMERG_DEPS := $(LOGLIKS)
EMERG_DEPS += $(CODEDIR)/make-pomp-simulator-function.R
EMERG_DEPS += $(CODEDIR)/simulate-emergence-grid.R
EMERG_TARG = $(SIMDIR)/emergence-simulations-grid-*.RDS

$(EMERG_TARG): $(EMERG_DEPS)
	Rscript $(CODEDIR)/simulate-emergence-grid.R $(NREPS)


## Simulate elimination replicates
ELIMIN_DEPS := $(LOGLIKS)
ELIMIN_DEPS += $(CODEDIR)/make-pomp-simulator-function.R
ELIMIN_DEPS += $(CODEDIR)/simulate-emergence-grid.R
ELIMIN_TARG = $(SIMDIR)/elimination-simulations-grid-*.RDS

$(ELIMIN_TARG): $(ELIMIN_DEPS)
	Rscript $(CODEDIR)/simulate-elimination-grid.R $(NREPS)


## Analyze the emergence sims (EWS and AUC)
EM_EWS_DEPS := $(SIMDIR)/emergence-simulations-grid-*.RDS
EM_EWS_DEPS += $(CODEDIR)/analyze-emergence-grid-sims.R
EM_EWS_TARG = $(RESDIR)/emergence-grid-aucs.csv

$(EM_EWS_TARG): $(EM_EWS_DEPS)
	Rscript $(CODEDIR)/analyze-emergence-grid-sims.R


## Analyze the elimination sims (EWS and AUC)
EL_EWS_DEPS := $(SIMDIR)/elimination-simulations-grid-*.RDS
EL_EWS_DEPS += $(CODEDIR)/analyze-elimination-grid-sims.R
EL_EWS_TARG = $(RESDIR)/elimination-grid-aucs.csv

$(EL_EWS_TARG): $(EL_EWS_DEPS)
	Rscript $(CODEDIR)/analyze-elimination-grid-sims.R


## Make the manuscript
MS_DEPS := $(MSDIR)/measles-ews-manuscript.Rmd
MS_DEPS += $(POMP_TARG)
MS_DEPS += $(PRED_TARG)
MS_DEPS += $(BOOT_TARG)
MS_DEPS += $(EM_EWS_TARG)
MS_DEPS += $(EL_EWS_TARG)

$(MSDIR)/measles-ews-manuscript.pdf: $(MS_DEPS)
	Rscript -e 'rmarkdown::render("$<")'
