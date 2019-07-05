# bin-level testing
# Pick up from clustering .cnr results

cnr_of_interest=TR_101.clustered.cnr TR_77.clustered.cnr TR_55.clustered.cnr TR_42.clustered.cnr
cns_of_interest=$(cnr_of_interest:.cnr=.cns)
# TR_101.clustered.cnr TR_77.clustered.cnr TR_55.clustered.cnr TR_42.clustered.cnr
out_cns=$(cnr_of_interest:.cnr=.bintest.cns)

all: $(out_cns)

$(out_cns): %.bintest.cns: %.cns %.cnr
	cnvkit.py bintest -s $^ -t -o $@

$(cns_of_interest): %.cns: %.cnr
	cnvkit.py segment $< -o $@
