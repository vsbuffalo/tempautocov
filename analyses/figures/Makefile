all: figures copy

figures: figure-1.pdf cummulative-cov-var.pdf cummulative-cov-var-all.pdf mom-fits-va1.pdf mom-fits-N.pdf sim-pred-covs-varyn.pdf estimate-g.pdf estimate-gp.pdf fluct-sel.pdf sim-pred13-covs-varyl.pdf sim-pred-covs-varyl-va-oom.pdf keynote-cartoons-figure-1.pdf unlinked_linked_contributions.pdf

.PHONY: copy clean

# cartoons
figure-1.tex: cartoons.r
	Rscript --vanilla cartoons.r

figure-1.pdf: figure-1.tex
	pdflatex $<

ALL_FIGS=sim-pred-covs-sl.pdf cummulative-cov-var.pdf expfit-vark-types.pdf sim-pred-covs-L500-ssh.pdf wvar-sl-trajectories.pdf sim-pred-covs-varyl-alt.pdf mom-fits-va1.pdf mom-fits-N.pdf mom-fits-both.pdf supp-ld-neutral.pdf supp-het-neutral.pdf supp-Ne-est-neutral.pdf sim-pred-covs-varyn.pdf estimate-g.pdf estimate-gp.pdf fluct-sel.pdf sim-pred13-covs-varyl.pdf sim-pred-covs-varyl-va-oom.pdf keynote-cartoons-wide.pdf keynote-cartoons-figure-1.pdf

DATA=$(addprefix '../../data/', cumcov_plot_all_df.rda cumcov_plot_df.rda pd5_vard.rda pd5f.rda pd5f_ssh.rda predf.rda predf_ssh.rda sl_pd5f.rda sl_predf.rda sl_wvar_trajs.rda vark_types.rda predf_varyn.rda pd5f_varyn.rda mom_fits.rda, cumd.rda, cumd_abs.rda predf_ssh2.rda)

# figures
$(ALL_FIGS): ms-plots.r ../../R/plots.r fluct-plot.r
	Rscript --vanilla $<

copy: $(ALL_FIGS)
	cp figure-1.pdf ../../manuscript/images/figure-1.pdf
	cp cummulative-cov-var.pdf ../../manuscript/images/
	cp cummulative-cov-var-all.pdf ../../manuscript/images/
	cp va-r-cov.pdf ../../manuscript/images
	cp sim-pred-covs-sl.pdf ../../manuscript/images/
	cp expfit-vark-types.pdf ../../manuscript/images/
	cp sim-pred-covs-L500-ssh.pdf ../../manuscript/images/
	cp wvar-sl-trajectories.pdf ../../manuscript/images/
	cp sim-pred-covs-varyl-alt.pdf ../../manuscript/images/
	cp mom-fits-va1.pdf ../../manuscript/images/
	cp mom-fits-N.pdf ../../manuscript/images/
	cp mom-fits-both.pdf ../../manuscript/images/
	cp supp-ld-neutral.pdf ../../manuscript/images/
	cp supp-het-neutral.pdf ../../manuscript/images/
	cp supp-Ne-est-neutral.pdf ../../manuscript/images/
	cp sim-pred-covs-varyn.pdf ../../manuscript/images/
	cp estimate-g.pdf ../../manuscript/images/
	cp estimate-gp.pdf ../../manuscript/images/
	cp sim-pred13-covs-varyl.pdf ../../manuscript/images/
	cp fluct-sel.pdf ../../manuscript/images/
	cp sim-pred-covs-varyl-va-oom.pdf ../../manuscript/images/
	cp keynote-cartoons-wide.pdf ../../manuscript/images/
	cp mom-fits-both-alt-finite.pdf ../../manuscript/images/
	cp mom-fits-both-alt.pdf ../../manuscript/images/
	cp supp-relative-error.pdf ../../manuscript/images/
	cp unlinked_linked_contributions.pdf ../../manuscript/images/

clean:
	rm -f $(ALL_FIGS)

clean-data: 
	mv $(DATA) ~/.Trash/
