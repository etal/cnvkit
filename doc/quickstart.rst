Quick start
===========

Installation
------------

Download the source code from GitHub:

http://github.com/etal/cnvkit

And read the README file.


Create targets from baits
-------------------------

If your baits are not already annotated with short gene names, download the
RefSeq gene annotations (refFlat.txt) from `UCSC Genome Bioinformatics
<hgdownload.soe.ucsc.edu/downloads.html>`_.
We'll apply these gene names to the bait intervals.

Assume your baits are in a BED file, ``baits.bed``. Run the ``target`` command
like this::

    cnvkit.py target baits.bed --refflat refFlat.txt --split -o my_targets.bed


Create antitargets from targets
-------------------------------

If your reference genome is the UCSC human genome hg19, a BED file of the
sequencing-accessible regions is included in the CNVkit distribution as
``data/access-10kb.hg19.bed``. Use it with the target file created above to
produce your antitarget file::

    cnvkit.py antitarget my_targets.bed -r data/access-10kb.hg19.bed -o my_antitargets.bed

If you're not using hg19, consider building the "access" file yourself for your
reference genome of interest (say, ``mm10``) using the bundled script
``genome2access.py``::

    genome2access.py mm10.fasta -s 10000 -o access-10kb.mm10.bed
    cnvkit.py antitarget my_targets.bed -r access-10kb.mm10.bed -o my_antitargets.bed


Build a reference from normal samples and infer tumor copy ratios
-----------------------------------------------------------------

We will use the target and antitarget files created above, along with the
reference genome you downloaded, to create a pooled reference of per-bin copy
number estimates from several normal samples. This reference is then used for
processing all tumor samples that were sequenced with the same platform and
library prep. All of these steps are automated with the ``batch`` command
(assuming normal samples share the suffix "Normal.bam" and tumor samples
"Tumor.bam")::

    cnvkit.py batch *Tumor.bam -n *Normal.bam -t my_targets.bed -a my_antitargets.bed -f hg19.fasta 

If you have no normal samples to use for the reference, you can create a "flat" 
reference which assumes equal coverage in all bins by using the ``-n`` flag
without specifying any additional BAM files::

    cnvkit.py batch *Tumor.bam -t my_targets.bed -a my_antitargets.bed -f hg19.fasta

In either case, you should run this command with the reference genome sequence
file (in FASTA format) to store GC and RepeatMasker information for bias
corrections. Also see the built-in help message for additional options::

    cnvkit.py batch -h


Process more tumor samples
--------------------------

Finally, you can reuse the reference file you've previously constructed to
extract copy number information from additional tumor sample BAM files, without
repeating the steps above. Assuming the new tumor samples share the suffix
"Tumor.bam" (and let's also use all available CPUs and output to another
directory)::

    cnvkit.py batch *Tumor.bam -r my_reference.cnn -p 0 -d output_dir/

