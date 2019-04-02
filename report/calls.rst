Variant calls obtained from Stacks, using the populations_ command.
Used parameters:

* max_locus_mm (Number of mismatches allowed between sample loci): {{ snakemake.wildcards.max_locus_mm }}
* max_individual_mm (Number of mismatches allowed between stacks): {{ snakemake.wildcards.max_individual_mm }}
* min_reads (Minimum depth of coverage required to create a stack): {{ snakemake.wildcards.min_reads }}


.. _populations: http://catchenlab.life.illinois.edu/stacks/comp/populations.php
