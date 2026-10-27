"""busco2summary - detailed multi-sample BUSCO post-processing.

Summaries, S/D/F/M matrices, per-gene statistics (with fragmented BUSCOs as a
first-class signal), occupancy-based gene-ID lists and sequence extraction for
BUSCO v3-v6 (genome / transcriptome / protein modes) and compleasm outputs.
"""

__version__ = "1.3.0"

#: Status codes used everywhere in the package.
STATUSES = ("S", "D", "F", "M")
