MAKEFLAGS += --no-builtin-rules
.SUFFIXES:
KDB=../Databases/std/database.kdb
KRAKENHLL_DIR="../krakenuniq/install"

all: hll-report.html

hll-report.html: hll-report.R hll-on-1m-db-kmers.txt random
	R -e "rmarkdown::render('hll-report.R', output_file='hll-report.html')"
	cp -v figures/relative-error-sd.std_x.percent_y.pdf ../Figures/Figure2-hll-relative-error-sd-on-db-kmers.pdf
	cp -v figures/hll-method-comparison-ertl_heule_flajolet.pdf ../Figures/Suppl-Figure1-hll-ertl_heule_flajolet.pdf
	cp -v figures/hll-method-comparison-ertl_heule.pdf ../Figures/Suppl-Figure2-hll-ertl_heule.pdf

hll-on-1m-db-kmers.txt: #$(KDB) $(KRAKENHLL_DIR)/test_hll_on_db
	$(KRAKENHLL_DIR)/test_hll_on_db $(KDB) 1 1010000 > hll-on-1m-dbkmers.txt

random: #$(KRAKENHLL_DIR)/count_unique
	mkdir -p random
	for P in 10 12 14 16 18; do \
        $(KRAKENHLL_DIR)/count_unique -p $P -r 10000000000 -s -t -y -f -e -x 100 > random/hll-p$P.txt; \
    done

