version 1.1
# usage: cnvkit.py antitarget [-h] [-g FILENAME] [-a AVG_SIZE] [-m MIN_SIZE]
#                             [-o FILENAME]
#                             targets

# positional arguments:
#   targets               BED or interval file listing the targeted regions.

# options:
#   -h, --help            show this help message and exit
#   -g FILENAME, --access FILENAME
#                         Regions of accessible sequence on chromosomes (.bed),
#                         as output by genome2access.py.
#   -a AVG_SIZE, --avg-size AVG_SIZE
#                         Average size of antitarget bins (results are
#                         approximate). [Default: 150000]
#   -m MIN_SIZE, --min-size MIN_SIZE
#                         Minimum size of antitarget bins (smaller regions are
#                         dropped). [Default: 1/16 avg size, calculated]
#   -o FILENAME, --output FILENAME
#                         Output file name.

task Antitarget {
  input {
    File targets_bed
    File access_bed
    Int avg_size=150000
    Int? min_size
    String output_filename = "out.antitarget.bed"
  }
  parameter_meta {
    targets_bed: "BED or interval file listing the targeted regions."
    access_bed: "Regions of accessible sequence on chromosomes (.bed), as output by the 'access' command."
    avg_size: "Average size of antitarget bins (results are approximate). [Default: 150000]"
    min_size: "Minimum size of antitarget bins (smaller regions aredropped). [Default: 1/16 avg size, calculated]"
    output_filename: "Output file name."
  }

  command <<<
    cnvkit.py antitarget ~{targets_bed} \
      ~{"-g " + access_bed} \
      ${~{avg_size}:+-a ~{avg_size}} \
      ${~{min_size}:+-m ~{min_size}} \
      -o ~{output_filename}
  >>>

  output {
    File output_bed = "${output_filename}"
  }

  runtime {
    docker: "etal/cnvkit:latest"
  }
}

