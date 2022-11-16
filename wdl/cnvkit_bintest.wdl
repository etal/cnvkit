version 1.1

task X {
  input {
    File cnarray
    File? segments
    Float alpha = 0.005
    Boolean test_only_targeted = false
    String output_filename = "out.bintest.cns"
  }
  parameter_meta {
    cnarray: "Bin-level log2 ratios (.cnr file), as produced by 'fix'."
    segments: "Segmentation calls (.cns), the output of the 'segment' command)."
    alpha: "Significance threhold. [Default: 0.005]"
    test_only_targeted: "Test target bins only; ignore off-target bins."
    output_filename: "Output filename."
  }

  # usage: cnvkit.py bintest [-h] [-s FILENAME] [-a ALPHA] [-t] [-o OUTPUT]
  #                          cnarray
  command <<<
    cnvkit.py bintest ~{cnarray} \
      ~{"-s " + segments} \
      -a ~{alpha} \
      ~{if test_only_targeted then "-t" else ""} \
      -o ~{output_filename}
  >>>

  output {
    File output_cns = "$output_filename"
  }

  runtime {
    docker: "etal/cnvkit:latest"
  }
}

