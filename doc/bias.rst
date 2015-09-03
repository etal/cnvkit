Bias corrections
================

The sequencing coverage depth obtained at different genomic regions is variable,
particularly for targeted capture. Much of this variability is due to known
biochemical effects related to library prep, target capture and sequencing.
Normalizing the read depth in each on-- and off-target bin to the expected read
depths derived from a reference of normal samples removes much of
the biases attributable to GC content, target density. However, these biases
also vary between samples, and must still be corrected even when a normal
reference is available.

To correct each of these known effects, CNVkit calculates the relationship
between observed bin-level read depths and the values of some known biasing
factor, such as GC content. This relationship is fitted using a simple rolling
median, then subtracted from the original read depths in a sample to yield
corrected estimates.

In the case of many similarly sized target regions, there is the potential for
the bias value to be identical for many targets, including some spatially near
each other.
To ensure that the calculated biases are independent of genomic position, the
probes are randomly shuffled before being sorted by bias value.

The GC content and repeat-masked fraction of each bin are calculated during
generation of the :ref:`reference` from the user-supplied genome. The bias
corrections are then performed in the :ref:`reference` and :ref:`fix` commands.


GC content
----------

Genomic regions with extreme GC content, the fraction of sequence composed of
guanine or cytosine bases, are less amenable to hybridization, amplification and
sequencing, and will generally appear to have lower coverage than regions of
average GC content.

To correct this bias in each sample, CNVkit calculates the association between
each bin's GC content (stored in the reference) and observed read depth, fits a
trendline through the bin read depths ordered by GC value, and subtracts this
trend from the original read depths.


Sequence repeats
----------------

Repetitive elements in the genome can be masked out with `RepeatMasker
<http://www.repeatmasker.org/>`_ -- and the genome sequences provided by the
`UCSC Genome Bioinformatics Site <http://genome.ucsc.edu/>`_ have this masking
applied already. The fraction of each genomic bin masked out for repetitiveness
indicates both low mappability and the susceptibility to Cot-1 blocking, both of
which can reduce the bin's observed coverage.

CNVkit removes the association between repeat-masked fraction and bin read
depths for each sample similarly to the GC correction.


Targeting density
-----------------

In hybridization capture, two biases occur near the edge of each baited region:

- Within the baited region, read depth is lower at the "shoulders" where
  sequence fragments are not completely captured.
- Just outside the baited region, in the "flanks", read depth is elevated
  to nearly that of the adjacent baited sites due to the same effect.
  If two targets are very close together, the sequence fragments captured for
  one target can increase the read depth in the adjacent target.

CNVkit roughly estimates the potential for these two biases based on the size
and position of each baited region and its immediate neighbors.
The biases are modeled together as a linear decrease in read depth from inside
the target region to the same distance outside.
These biases occur within a distance of the interval edges equal to the sequence
fragment size (also called the insert size for paired-end sequencing reads).
Density biases are calculated from the start and end positions of a bin and its
neighbors within a fixed window around the target's genomic coordinates equal to
the sequence fragment size.

*Shoulder effect:* Letting *i* be the average insert size and *t* be the target
interval size, the negative bias at interval shoulders is calculated as
:math:`i/4t` at each side of the interval, or :math:`i/2t` for the whole interval.
When the interval is smaller than the sequence fragment size, the portion of the
fragment extending beyond the opposite edge of the interval should not be
counted in this calculation.
Thus, if :math:`t < i`, the negative bias value must be increased (absolute
value reduced) by :math:`\frac{(i-t)^2}{2it}`.

*Flank effect:* Additionally letting *g* be the size of the gap between
consecutive intervals, the positive bias that occurs when the gap is smaller
than the insert size (:math:`g<i`) is :math:`\frac{(i-g)^2}{4it}`.
If the target interval and gap together are smaller than the insert size, the
reads flanking the neighboring interval may extend beyond the target, and this
flanking portion beyond the target should not be counted.
Thus, if :math:`t+g < i`, the positive value must be reduced by
:math:`\frac{(i-g-t)^2}{4it}`.
If a target has no close neighbors (:math:`g>i`, the common case), the "flank"
bias value is 0.

These values are combined into a single value by subtracting the estimated
shoulder biases from the flank biases.
The result is a negative number between -1 and 0, or 0 for a target with
immediately adjacent targets on both sides.  Thus, subdividing a large targeted
interval into a consecutive series of smaller targets does not change the net
"density" calculation value.

The association between targeting density and bin read depths is then fitted and
subtracted, as with GC and RepeatMasker.

CNVkit applies the density bias correction to only the on-target bins; the
negative "shoulder" bias is not expected to occur in off-target regions because
those regions are not specifically captured by baits, and the positive "flank"
bias from neighboring targets is avoided by allocating off-target bins around
existing targets with a margin of twice the expected insert size.
