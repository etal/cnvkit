version 1.0

# TODO - split out separate workflows for reference-only, all-in-one, reference-given, reference-free, etc.
workflow CNVkitBatch {
	call Main

	meta {
		author: "Eric Talevich"
		email: "etalevich@dnanexus.com"
		description: "Call copy number variants from DNA sequencing, and/or build a reference profile"
	}
}

task Main {
	Array[File] case_bams
	Array[File] normal_bams
	Array[File]? snv_vcfs
	File? cn_reference
	File? baits=None
	File? fasta=None
	File? annotation=None,
	String seq_method = 'hybrid'
	String segment_method = 'cbs'
	Boolean haploid_x_reference = false
	Boolean drop_low_coverage = false
	Array[File] exclude_access
	Int? antitarget_avg_size
	Int? target_avg_size
	Float? purity
	Int? ploidy
	Boolean do_parallel = true

	command {
		cnvkit.py version
		cnvkit.py access ${fasta} --output access.bed
		cnvkit.py target
		cnvkit.py antitarget ${baits} --access access.bed --output antitargets.bed
		cnvkit.py autobin
	}
	output {
		File cn_reference
		Array[File] copy_ratios
		Array[File] copy_segments
		Array[File] call_segments
		Array[File] genemetrics
		Array[File] cnv_beds
		Array[File] cnv_vcfs
		File seg
		File metrics
		File sexes
		File heatmap_pdf
		Array[File] scatters_png
	}
	runtime {
		docker: "etal/cnvkit:latest"
		cpu: 1
		memory: "512 MB"
	}
}

task RunCoverage {
	File bam
	File targets
	File antitargets
	Boolean do_parallel

	command {
		cnvkit.py coverage
	}
	output {
		Array[File] coverages
	}
	runtime {
		docker: "etal/cnvkit:latest"
		cpu: 8
		memory: "2 GB"
	}
}

task RunReference {
	Array[File] coverages
	File fasta
	File targets
	File antitargets
	Boolean haploid_x_reference

	command {
		cnvkit.py reference --fasta --output
	}
	output {
		File cn_reference
	}
	runtime {
		docker: "etal/cnvkit:latest"
		cpu: 8
		memory: "2 GB"
	}
}

task RunSample {
	File sample_bam
	String seq_method
	String segment_method
	File cn_reference
	File vcf
	Float purity
	Int ploidy
    Boolean drop_low_coverage
	Boolean haploid_x_reference
	Boolean do_parallel

	command {
		cnvkit.py batch --method --reference --segment-method --scatter
		ls -Altr
		cnvkit.py export vcf
		cnvkit.py export bed
		cnvkit.py genemetrics --segment --min-probes 3 --threshold 0.2 --output
	}
	output {
		File copy_ratios
		File copy_segments
		File call_segments
		File genemetrics
		File sample_vcf
		File sample_bed
		File scatter_pdf
	}
	runtime {
		docker: "etal/cnvkit:latest"
		cpu: 8
		memory: "2 GB"
	}
}

task AggregateOutputs {
	Array[File] copy_ratios
	Array[File] copy_segments
	Boolean haploid_x_reference

	command {
		cnvkit.py metrics -s --output
		cnvkit.py export seg --output
		cnvkit.py heatmap -d --output
		cnvkit.py sex --output
	}
	output {
	}
	runtime {
		docker: "etal/cnvkit:latest"
		cpu: 8
		memory: "2 GB"
	}
}