# Test CNVkit.
#
# Dependency: pdfunite (poppler-utils)
# (Otherwise, all-scatters.pdf and all-diagrams.pdf will be empty files.)

rscript_exe=Rscript
cnvkit=cnvkit.py

# ------------------------------------------------------------------------------
# Samples pre-processed with Picard CalculateHsMetrics

picard_normals=$(wildcard picard/p2-*_5.*.csv)
picard_tumor_cnns=$(patsubst picard/%.csv,build/%.cnn,$(wildcard picard/*.csv))

picard_cnrs=$(addprefix build/p2-,$(addsuffix .cnr,5_5 9_2 20_1 20_2 20_3 20_4 20_5))
picard_segs=$(picard_cnrs:.cnr=.cns)

picard_targets= build/reference-picard.cnn $(picard_cnrs) $(picard_segs) \
	     all-scatters.pdf all-diagrams.pdf heatmap-picard.pdf \
	     p2-9_2-breaks.txt p2-9_2-genemetrics.txt gender-picard.txt \
	     p2-all.cdt p2-all-jtv.txt p2-all.seg p2-9_2.nexus p2-9_2.nexus-ogt \
	     p2-all.bed p2-9_2.vcf p2-9_2.theta2.input \
	     p2-9_2-segmetrics.cns p2-5_5-segmetrics.cns

# ------------------------------------------------------------------------------
#  Action!

all: mini $(picard_targets) p2-5_5-metrics.tsv p2-9_2-metrics.tsv p2-20-metrics.tsv


.PHONY: clean
clean:
	# Picard
	rm -rf build/p*targetcoverage.cnn $(picard_targets) p*-scatter.pdf p*-diagram.pdf
	rm -rf build/na12878-chrM-Y-trunc* build/chrM-Y-trunc*

.PHONY: test
test:
	pytest


# ------------------------------------------------------------------------------
# Minimal 'batch' command

.PHONY: mini
mini: build/na12878-chrM-Y-trunc.bintest.cns

build/na12878-chrM-Y-trunc.bintest.cns: formats/na12878-chrM-Y-trunc.bam formats/chrM-Y-trunc.hg19.fa
	$(cnvkit) batch -n -f formats/chrM-Y-trunc.hg19.fa -m wgs formats/na12878-chrM-Y-trunc.bam -d build


# ------------------------------------------------------------------------------
# Standard workflow

# == Build pooled references from normal samples

build/reference-picard.cnn: $(picard_normals)
	$(cnvkit) import-picard $^ -d build/
	$(cnvkit) reference build/p2-*_5.*targetcoverage.cnn -y -o $@


# == Import Picard files (including GC values, used in the reference)

$(picard_tumor_cnns): build/%.cnn: picard/%.csv
	$(cnvkit) import-picard $^ -d build/


# == Build components

$(picard_cnrs): %.cnr: %.targetcoverage.cnn %.antitargetcoverage.cnn build/reference-picard.cnn
	$(cnvkit) fix $^ -o $@

$(picard_segs): %.cns: %.cnr
	$(cnvkit) segment --rscript-path $(rscript_exe) -p 2 --drop-low-coverage -t .01 $< -o $@

$(picard_segs:.cns=.call.cns): build/%.call.cns: build/%.cns
	$(cnvkit) call $^ -y -m clonal --purity 0.65 -o $@


# ------------------------------------------------------------------------------
# Demo other features with Picard outputs

heatmap-picard.pdf: $(picard_segs)
	$(cnvkit) heatmap $^ -y -o $@

# == Plot coverages with segmentation calls

all-scatters.pdf: p2-5_5-scatter.pdf p2-9_2-scatter.pdf p2-9_2-chr1-scatter.pdf p2-9_2-chr21-scatter.pdf p2-9_2-chr9p-scatter.pdf p2-9_2-SMARCA2-PTPRD-scatter.pdf p2-9_2-bed_regions-scatter.pdf
	which pdfunite && pdfunite $^ $@ || touch $@

p2-5_5-scatter.pdf p2-9_2-scatter.pdf: %-scatter.pdf: build/%.cns build/%.cnr
	$(cnvkit) scatter -s $^ -t --y-min=-2.5 -o $@

p2-9_2-chr1-scatter.pdf: build/p2-9_2.cns build/p2-9_2.cnr
	$(cnvkit) scatter -s $^ -c chr1 -t -o $@

p2-9_2-chr21-scatter.pdf: build/p2-9_2.cns build/p2-9_2.cnr
	$(cnvkit) scatter -s $^ -c chr21 -t -o $@

p2-9_2-chr9p-scatter.pdf: build/p2-9_2.cns build/p2-9_2.cnr
	$(cnvkit) scatter -s $^ -c 'chr9:150000-45000000' -o $@

p2-9_2-SMARCA2-PTPRD-scatter.pdf: build/p2-9_2.cns build/p2-9_2.cnr
	$(cnvkit) scatter -s $^ -g SMARCA2,PTPRD -w 4e6 -o $@

p2-9_2-bed_regions-scatter.pdf: build/p2-9_2.cns build/p2-9_2.cnr
	$(cnvkit) scatter -s $^ -l regions.bed -o $@

# == Draw diagrams

all-diagrams.pdf: $(foreach sfx,- -cbs- -both-,$(addsuffix $(sfx)diagram.pdf,p2-5_5 p2-9_2 p2-20_1 p2-20_2))
	which pdfunite && pdfunite $^ $@ || touch $@

p2-5_5-diagram.pdf p2-9_2-diagram.pdf p2-20_1-diagram.pdf p2-20_2-diagram.pdf: %-diagram.pdf: build/%.cnr
	$(cnvkit) diagram -y $< -o $@

p2-5_5-cbs-diagram.pdf p2-9_2-cbs-diagram.pdf p2-20_1-cbs-diagram.pdf p2-20_2-cbs-diagram.pdf: %-cbs-diagram.pdf: build/%.cns
	$(cnvkit) diagram -y --segment=$< -o $@

p2-5_5-both-diagram.pdf p2-9_2-both-diagram.pdf p2-20_1-both-diagram.pdf p2-20_2-both-diagram.pdf: %-both-diagram.pdf: build/%.cns build/%.cnr
	$(cnvkit) diagram -y --segment=$^ -t 0.3 -o $@

# == Text reporting ==

p2-9_2-breaks.txt: build/p2-9_2.cnr build/p2-9_2.cns
	$(cnvkit) breaks $^ -o $@

p2-9_2-genemetrics.txt: build/p2-9_2.cns build/p2-9_2.cnr
	$(cnvkit) genemetrics -y -m 2 -s $^ -o $@

gender-picard.txt: build/p2-5_5.cnr  build/p2-5_5.cns build/p2-9_2.cnr build/p2-9_2.cns build/p2-20_5.cnr build/p2-20_5.cns
	$(cnvkit) sex -y $^ -o $@

p2-5_5-metrics.tsv p2-9_2-metrics.tsv: %-metrics.tsv: build/%.cnr build/%.cns
	$(cnvkit) metrics  $< -s $(lastword $^) -o $@

p2-20-metrics.tsv: $(addprefix build/p2-20_,$(addsuffix .cns,1 2 3 4 5))
	$(cnvkit) metrics $(^:.cns=.cnr) -s $^ -o $@

p2-9_2-segmetrics.cns: build/p2-9_2.cns build/p2-9_2.cnr
	$(cnvkit) segmetrics -s $^ -o $@ \
		--mean --median --mode --t-test \
		--stdev --mad --mse --iqr --bivar \
		--ci --pi --sem --smooth-bootstrap

p2-5_5-segmetrics.cns: build/p2-5_5.cns build/p2-5_5.cnr
	$(cnvkit) segmetrics -s $^ -o $@ \
		--ci -b 50 -a 0.5


# == Export to other formats

p2-all.seg: $(picard_segs)
	$(cnvkit) export seg $^ -o $@

p2-all.bed: $(picard_segs:.cns=.call.cns)
	$(cnvkit) export bed $^ -y --show variant -o $@

p2-9_2.vcf: build/p2-9_2.cnr build/p2-9_2.call.cns
	$(cnvkit) export vcf -o $@ -y --cnr $^ -y -o $@

p2-9_2.theta2.input: build/p2-9_2.cns build/reference-picard.cnn
	$(cnvkit) export theta $< -r $(lastword $^) -o $@

p2-9_2.nexus: build/p2-9_2.cnr
	$(cnvkit) export nexus-basic $^ -o $@

p2-9_2.nexus-ogt: build/p2-9_2.cnr formats/na12878_na12882_mix.vcf
	$(cnvkit) export nexus-ogt $^ -o $@

p2-all.cdt: $(picard_cnrs)
	$(cnvkit) export cdt $^ -o $@

p2-all-jtv.txt: $(picard_cnrs)
	$(cnvkit) export jtv $^ -o $@
