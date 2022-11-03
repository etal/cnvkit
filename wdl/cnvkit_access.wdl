# usage: cnvkit.py access [-h] [-s MIN_GAP_SIZE] [-x EXCLUDE] [-o FILENAME]
# 												fa_fname
# 
# positional arguments:
# 	fa_fname              Genome FASTA file name
# 
# options:
# 	-h, --help            show this help message and exit
# 	-s MIN_GAP_SIZE, --min-gap-size MIN_GAP_SIZE
# 												Minimum gap size between accessible sequence regions.
# 												Regions separated by less than this distance will be
# 												joined together. [Default: 5000]
# 	-x EXCLUDE, --exclude EXCLUDE
# 												Additional regions to exclude, in BED format. Can be
# 												used multiple times.
# 	-o FILENAME, --output FILENAME
# 												Output file name

task Access {
  input {
    File fa_fname
    Int? min_gap_size=5000
    Array[File]? exclude_bed
    String? output_filename
  }
  parameter_meta {
    fa_fname: "Genome FASTA file name"
  	min_gap_size: "Minimum gap size between accessible sequence regions. Regions separated by less than this distance will be joined together. [Default: 5000]"
  	exclude_bed: "Additional regions to exclude, in BED format. Can be used multiple times."
  	output_fname: "Output file name"
  }

  output {
    File output_bed
  }

  command <<<
    cnvkit.py access -o $output_bed
# usage: cnvkit.py access [-h] [-s MIN_GAP_SIZE] [-x EXCLUDE] [-o FILENAME]
# 												fa_fname
  >>>

  runtime {
    docker: "etal/cnvkit:latest"
  }
}
