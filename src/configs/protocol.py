import addict


PROTO = addict.Dict({})
PROTO.SEED = 42

PROTO.EXPERIMENT_ID = 'SUBSET_1'
PROTO.DESCRIPTION = 'WO OPEN Dataset, WO BANOBO'

INTEGRATES = [ 'harmony', 'scanorama' ]
PROTO.INTEGRATE = addict.Dict({
    'TYPE': INTEGRATES[0],
    'FIELD': 'batch' })

# Datasets in use:
PROTO.DATASETS = addict.Dict({
    'E-MATB': True, # to use E-MATB ineuron dataset;
    'MENDELEY': False, # to use Mendeley Data scRNAseq;
    'OPEN': False, # to use Open Data scRNAseq (17-18 th gestation week)
    'BATCH_FROM_TIMEPOINT': False, # to infer E-MATB dataset batch info from timepoints or from the cell 
})

# Subset selection:
HVG_ALGOS = [ 'triku', 'hvgsc' ]
PROTO.SUBSET = addict.Dict({
    'APE_HUMAN_OVERLAP': True, # to use only overlap of gene sets;
    'BANOBO': False, # whether to include banobo or drop it;
    'TIMEPOINTS': [14, 28, 35, 119, 126], # timepoints to consider;
    'MITO_GENES': False, # whether to include mito-genes;
    'RIBO_GENES': False, # whether to include ribo-genes;
    'USE_HIGHLY_VARIABLE': True, # to subset HV genes;
    'NB_HIGHLY_VARIABLE': 5000, # None;
    'HVG_ALGOS': HVG_ALGOS[0], # Whether to use triku, or sc HVG;
    'BATCH_KEY': '10X_date', # or None
})

# Pseudo bulk generation
PSEUDO_BULKS = [ 'within_group_replace', 'wo_replace', 'decoupler' ]
PROTO.PSEUDO_BULK = addict.Dict({})
PROTO.PSEUDO_BULK.TYPE = PSEUDO_BULKS[-1]
PROTO.PSEUDO_BULK.SIZE = 75
PROTO.PSEUDO_BULK.REPS = 2


# GLM contrast
PSEUDOBULK_MODES = [ 'sum', 'median' ]
PROTO.GLM = addict.Dict({
    'CONTRAST': '',
    'PSEUDOBULK_MODE': PSEUDOBULK_MODES[0],
    'MIN_CELLS': 30,
    'DEG_THRESHOLD': { 'FDR': .01, 'LFC': 1.5 },
    'FILTER': False,
})

COLLECTIONS = [ 
    'go_biological_process', 'hallmark', 'kegg_pathways',
    'chemical_and_genetic_perturbations',
]
PROTO.ENRICHMENT = addict.Dict({
    'COLLECTION': COLLECTIONS[2],
})

# Quality Control
PROTO.QC = addict.Dict({})
PROTO.QC.MIN_GENES = 500
PROTO.QC.MAX_GENES = 5000

PROTO.QC.MIN_COUNTS= 750
PROTO.QC.MAX_COUNTS= 10000

PROTO.QC.MT_PERCENTAGE = 5
PROTO.QC.RB_PERCENTAGE = 35
PROTO.QC.MIN_CELLS_PCT = .5
