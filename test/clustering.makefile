# Clustered reference method
example_cnn_dir=~/code/cnvkit-examples/targeted

cnr_of_interest=TR_101.clustered.cnr TR_77.clustered.cnr TR_55.clustered.cnr TR_42.clustered.cnr

all: $(cnr_of_interest)

clean:
	rm -vf $(cnr_of_interest) ref-clustered.cnn

$(cnr_of_interest): \
	TR_%.clustered.cnr: \
	$(example_cnn_dir)/TR_%_T.targetcoverage.cnn \
	$(example_cnn_dir)/TR_%_T.antitargetcoverage.cnn \
	ref-clustered.cnn
	cnvkit.py fix $^ --cluster -o $@

ref-clustered.cnn: $(wildcard $(example_cnn_dir)/TR_*targetcoverage.cnn)
	cnvkit.py reference $^ --cluster -o $@
