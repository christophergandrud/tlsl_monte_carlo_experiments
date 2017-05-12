RDIR = .
SETUP_OUT = a_setup.Rout
PLOTS_OUT = z_mc_results_plots.Rout

RSOURCE = $(wildcard $(RDIR)/*.R)

OUT_FILES = $(RSOURCE:.R=.Rout)

all: $(OUT_FILES)

$(RDIR)/%.Rout: $(RDIR)/%.R
		R CMD BATCH $<

clean:
	rm -fv $(OUT_FILES)

cleanSet:
	rm -fv $(SETUP_OUT)
	rm -fv $(PLOTS_OUT)
