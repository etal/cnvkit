version 1.1

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
    Int min_gap_size = 5000
    Array[File]? exclude_bed
    String output_filename = "out.access.bed"
  }
  parameter_meta {
    fa_fname: "Genome FASTA file name"
    min_gap_size: "Minimum gap size between accessible sequence regions. Regions separated by less than this distance will be joined together. [Default: 5000]"
    exclude_bed: "Additional regions to exclude, in BED format. Can be used multiple times."
    output_filename: "Output file name"
  }

  command <<<
    cnvkit.py access ~{fa_fname} \
      -s ~{min_gap_size} \
      ~{sep(" ", prefix("-x ",  exclude_bed))} \
      -o ~{output_filename}
  >>>

  output {
    File output_bed = "$output_filename"
  }

  runtime {
    docker: "etal/cnvkit:latest"
  }
}
