bin_SCRIPTS = pve_insect discretize_pve

EXTRA_DIST = pve_insect.pl discretize_pve.pl

SUFFIXES = .pl
.pl:
	@rm -f $@
	sed < $< > $@-t \
		-e 's,xPERLx,$(PERL),' \
		-e 's,xINCDIRx,$(perllibdir),g' \
		-e 's,xDATADIRx,$(datadir),g' \
		-e 's,xPKGDATADIRx,$(pkgdatadir),g' &&\
	mv $@-t $@

CLEANFILES = pve_insect discretize_pve

