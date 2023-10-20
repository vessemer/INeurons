import addict


PROTO = addict.Dict({})
PROTO.SEED = 42

PROTO.EXPERIMENT_ID = 'SUBSET_1'
PROTO.DESCRIPTION = 'WO OPEN Dataset, WO BANOBO'

# Datasets in use:
PROTO.DATASETS = addict.Dict({
    'E-MATB': True, # to use E-MATB ineuron dataset;
    'MENDELEY': True, # to use Mendeley Data scRNAseq;
    'OPEN': False, # to use Open Data scRNAseq (17-18 th gestation week)
})

# Subset selection:
PROTO.SUBSET = addict.Dict({
    'APE_HUMAN_OVERLAP': True, # to use only overlap of gene sets;
    'BANOBO': False, # whether to include banobo or drop it;
    'TIMEPOINTS': [14, 28, 35, 119, 126], # timepoints to consider.
})

# Pseudo bulk generation
PSEUDO_BULKS = [ 'within_group_replace', 'wo_replace', 'decoupler' ]
PROTO.PSEUDO_BULK = addict.Dict({})
PROTO.PSEUDO_BULK.TYPE = PSEUDO_BULKS[-1]
PROTO.PSEUDO_BULK.SIZE = 75
PROTO.PSEUDO_BULK.REPS = 2


# GLM contrast
PROTO.CONTRAST = ''


# Quality Control
PROTO.QC = addict.Dict({})
PROTO.QC.MIN_GENES = 500
PROTO.QC.MAX_GENES = 5000

PROTO.QC.MIN_COUNTS= 750
PROTO.QC.MAX_COUNTS= 10000

PROTO.QC.MT_PERCENTAGE = 5
PROTO.QC.RB_PERCENTAGE = 35
PROTO.QC.MIN_CELLS_PCT = .5
