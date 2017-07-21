# Binding factor(s) to find and mark CpG islands.

# based on Wachter et al., 2014 (Adrian Bird lab). https://elifesciences.org/content/3/e03397
# TODO have references as meta-data for models (as per Bio-models).

#   Mammal genomes CpG island > 500bp, GC% > 60%, CpG density > 1/10 nucleotides.

# CpG islands almost synonymous with bivalent promoters.

#  GC% alone not enough to recruit H3K4me3 and H3K27me3.
#  CpG alone (in otherwise AT rich region) also not sufficient - seemingly due to methylation in that context.
# CpGs in GC-rich (>60%) avoid methylation 
# other papers (which?) suggest the lack of methylation is actively maintained by binding of TFs.

# methylated CpGs cannot recruit the histone marks H3K4me3 and H3K27me3
# AT rich regions methylate CpGs
# GC rich regions do not methylate CpGs

# what binding factor mechanisms would work here? 

# %GC e.g. match to SSSSSSSS....SSS with mismatch.rate < 0.4
# CpG  e.g. match to CG
# But what about density of CpG?  Needs to be two-stage process?  Region matching with % overlap - not easy with hits() objects








