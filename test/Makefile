# Test CNVkit.
#
# Dependency: pdfunite (poppler-utils)
# (Otherwise, all-scatters.pdf and all-diagrams.pdf will be empty files.)

cnvkit=python ../cnvkit.py

# ------------------------------------------------------------------------------
# Samples pre-processed with Picard CalculateHsMetrics

picard_normals=$(wildcard picard/p2-*_5.*.csv)
picard_tumor_cnns=$(patsubst picard/%.csv,build/%.cnn,$(wildcard picard/*.csv))

picard_cnrs=$(addprefix build/p2-,$(addsuffix .cnr,5_5 9_2 20_1 20_2 20_3 20_4 20_5))
picard_segs=$(picard_cnrs:.cnr=.cns)

picard_targets= build/reference-picard.cnn $(picard_cnrs) $(picard_segs) \
	     p2-9_2.cdt p2-all.cdt p2-all-jtv.txt p2-all-cbs.seg \
	     all-scatters.pdf all-diagrams.pdf heatmap-picard.pdf \
	     p2-9_2-breaks.txt p2-9_2-gainloss.txt gender-picard.txt


# ------------------------------------------------------------------------------
# Targeted resequencing samples ("TR")

tr_tcnn=$(wildcard targeted/TR*.targetcoverage.cnn)
tr_cnrs=$(patsubst targeted/%.targetcoverage.cnn,build/%.cnr,$(tr_tcnn))
tr_segs=$(tr_cnrs:.cnr=.cns)

# Pooled reference (normal samples)
tr_ref_cnns=$(wildcard targeted/TR_*_N.*.cnn)


# ------------------------------------------------------------------------------
# Exome samples ("EX")

ex_tcnn=$(wildcard exome/EX*.targetcoverage.cnn)
ex_cnrs=$(patsubst exome/%.targetcoverage.cnn,build/%.cnr,$(ex_tcnn))
ex_segs=$(ex_cnrs:.cnr=.cns)

# Reference: normal male samples
ex_ref_cnns=$(wildcard exome/EX_*_N.*.cnn)


# ------------------------------------------------------------------------------
#  Action!

.PHONY: picard
picard: $(picard_targets) p2-5_5-metrics.txt p2-9_2-metrics.txt p2-20-metrics.csv

.PHONY: tr
tr: heatmap-tr.pdf tr-metrics.csv

.PHONY: ex
ex: $(ex_segs) ex-metrics.csv heatmap-exome.pdf


all: picard tr ex


.PHONY: clean
clean:
	# Picard
	rm -rf build/p*targetcoverage.cnn $(picard_targets) p*-scatter.pdf p*-diagram.pdf

.PHONY: allclean
allclean:
	# Picard
	rm -rf build/p*targetcoverage.cnn $(picard_targets) p*-scatter.pdf p*-diagram.pdf
	# Targeted
	rm -vf build/TR_* heatmap-tr.pdf
	# Exome
	rm -vf build/EX_* heatmap-exome.pdf

# ------------------------------------------------------------------------------
# Standard workflow

# == Build pooled references from normal samples

build/reference-picard.cnn: $(picard_normals)
	$(cnvkit) import-picard $^ -d build/
	$(cnvkit) reference build/p2-*_5.*targetcoverage.cnn -o $@

# * Non-picard references: uncomment to rebuild, but need hg19 FASTA genome

# reference-tr.cnn: $(tr_ref_cnns)
#         $(cnvkit) reference $^ -f /db/ucsc.hg19.fasta -y -o $@

# reference-exome.cnn: $(ex_ref_cnns)
#         $(cnvkit) reference $^ -f /db/hg19-nochr.fa -y -o $@


# == Import Picard files (including GC values, used in the reference)

$(picard_tumor_cnns): build/%.cnn: picard/%.csv
	$(cnvkit) import-picard $^ -d build/


# == Build components

$(picard_cnrs): %.cnr: %.targetcoverage.cnn %.antitargetcoverage.cnn build/reference-picard.cnn
	$(cnvkit) fix $^ -o $@

$(tr_cnrs): build/%.cnr: targeted/%.targetcoverage.cnn targeted/%.antitargetcoverage.cnn reference-tr.cnn
	$(cnvkit) fix $^ -o $@

$(ex_cnrs): build/%.cnr: exome/%.targetcoverage.cnn exome/%.antitargetcoverage.cnn reference-exome.cnn
	$(cnvkit) fix $^ -o $@

$(picard_segs) $(tr_segs) $(ex_segs): %.cns: %.cnr
	$(cnvkit) segment $< -o $@


# == Results

heatmap-picard.pdf: $(picard_segs)
	$(cnvkit) heatmap -d $^ -o $@

heatmap-tr.pdf: $(tr_segs)
	$(cnvkit) heatmap -d $^ -o $@

heatmap-exome.pdf: $(ex_segs)
	$(cnvkit) heatmap -d -o $@ $^


tr-metrics.csv: $(tr_segs)
	$(cnvkit) metrics $(tr_cnrs) -s $^ -o $@

ex-metrics.csv: $(ex_segs)
	$(cnvkit) metrics $(ex_cnrs) -s $^ -o $@


# ------------------------------------------------------------------------------
# Demo other features with Picard outputs

# == Plot coverages with CBS calls ==

all-scatters.pdf: p2-5_5-scatter.pdf p2-9_2-scatter.pdf p2-9_2-chr1-scatter.pdf p2-9_2-chr9p-scatter.pdf p2-9_2-SMARCA2-PTPRD-scatter.pdf
	which pdfunite && pdfunite $^ $@ || touch $@

p2-5_5-scatter.pdf p2-9_2-scatter.pdf: %-scatter.pdf: build/%.cns build/%.cnr
	$(cnvkit) scatter -s $^ -t -o $@

p2-9_2-chr1-scatter.pdf: build/p2-9_2.cns build/p2-9_2.cnr
	$(cnvkit) scatter -s $^ -c chr1 -t -o $@

p2-9_2-chr9p-scatter.pdf: build/p2-9_2.cns build/p2-9_2.cnr
	$(cnvkit) scatter -s $^ -r 'chr9:150000-45000000' -o $@

p2-9_2-SMARCA2-PTPRD-scatter.pdf: build/p2-9_2.cns build/p2-9_2.cnr
	$(cnvkit) scatter -s $^ -g SMARCA2,PTPRD -w 4e6 -o $@

# == Draw diagrams

all-diagrams.pdf: $(foreach sfx,- -cbs- -both-,$(addsuffix $(sfx)diagram.pdf,p2-5_5 p2-9_2 p2-20_1 p2-20_2))
	which pdfunite && pdfunite $^ $@ || touch $@

p2-5_5-diagram.pdf p2-9_2-diagram.pdf p2-20_1-diagram.pdf p2-20_2-diagram.pdf: %-diagram.pdf: build/%.cnr
	$(cnvkit) diagram -y $< -o $@

p2-5_5-cbs-diagram.pdf p2-9_2-cbs-diagram.pdf p2-20_1-cbs-diagram.pdf p2-20_2-cbs-diagram.pdf: %-cbs-diagram.pdf: build/%.cns
	$(cnvkit) diagram -y --segment=$< -o $@

p2-5_5-both-diagram.pdf p2-9_2-both-diagram.pdf p2-20_1-both-diagram.pdf p2-20_2-both-diagram.pdf: %-both-diagram.pdf: build/%.cns build/%.cnr
	$(cnvkit) diagram -y --segment=$^ -t 0.5 -o $@

# == Text reporting ==

p2-9_2-breaks.txt: build/p2-9_2.cnr build/p2-9_2.cns
	$(cnvkit) breaks $^ -o $@

p2-9_2-gainloss.txt: build/p2-9_2.cns build/p2-9_2.cnr
	$(cnvkit) gainloss -y -m 2 -s $^ -o $@

gender-picard.txt: build/p2-5_5.cnr  build/p2-5_5.cns build/p2-9_2.cnr build/p2-9_2.cns build/p2-20_5.cnr build/p2-20_5.cns
	$(cnvkit) gender -y $^ -o $@

p2-5_5-metrics.txt p2-9_2-metrics.txt: %-metrics.txt: build/%.cnr build/%.cns
	$(cnvkit) metrics  $< -s $(lastword $^) -o $@

p2-20-metrics.csv: $(addprefix build/p2-20_,$(addsuffix .cns,1 2 3 4 5))
	$(cnvkit) metrics $(^:.cns=.cnr) -s $^ -o $@


# == Export to other formats

p2-9_2.cdt: build/p2-9_2.cnr
	$(cnvkit) export cdt $^ -o $@

p2-all.cdt: $(picard_cnrs)
	$(cnvkit) export cdt $^ -o $@

p2-all-jtv.txt: $(picard_cnrs)
	$(cnvkit) export jtv $^ -o $@

p2-all-cbs.seg: $(picard_segs)
	$(cnvkit) export seg $^ -o $@
