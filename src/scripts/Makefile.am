bin_SCRIPTS = pve_insect discretize_pve pve_curvature

EXTRA_DIST = pve_insect.pl discretize_pve.pl pve_curvature.pl

SUFFIXES = .pl
.pl:
	@rm -f $@
	sed < $< > $@-t \
		-e 's,xPERLx,$(PERL),' \
		-e 's,xINCDIRx,$(perllibdir),g' \
		-e 's,xDATADIRx,$(datadir),g' \
		-e 's,xPKGDATADIRx,$(pkgdatadir),g' &&\
	mv $@-t $@

CLEANFILES = pve_insect discretize_pve pve_curvature

