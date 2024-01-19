import addict


PROTO = addict.Dict({})
PROTO.SEED = 42

PROTO.EXPERIMENT_ID = 'EMATB_MENDELEY_1'
PROTO.DESCRIPTION = 'WO OPEN Dataset, WO BANOBO, W d2 & d5'

INTEGRATES = [ 'harmony', 'scanorama' ]
PROTO.INTEGRATE = addict.Dict({
    'TYPE': INTEGRATES[0],
    'FIELD': 'batch' })

# Datasets in use:
PROTO.DATASETS = addict.Dict({
    'E-MATB': True, # to use E-MATB ineuron dataset;
    'MENDELEY': True, # to use Mendeley Data scRNAseq;
    'OPEN': False, # to use Open Data scRNAseq (17-18 th gestation week)
    'BATCH_FROM_TIMEPOINT': False, # to infer E-MATB dataset batch info from timepoints or from the cell 
})

# Subset selection:
HVG_ALGOS = [ 'triku', 'hvgsc' ]
EMATB_CLUSTERS = [ 'C{}'.format(i+1) for i in range(8) ]
PROTO.SUBSET = addict.Dict({
    'EMATB_CLUSTERS': EMATB_CLUSTERS[3:-1],
    'APE_HUMAN_OVERLAP': True, # to use only overlap of gene sets;
    'BANOBO': False, # whether to include banobo or drop it;
    'TIMEPOINTS': [ 5, 14, 28, 35 ], # timepoints to consider;
    'MITO_GENES': False, # whether to include mito-genes;
    'RIBO_GENES': False, # whether to include ribo-genes;
    'USE_HIGHLY_VARIABLE': True, # to subset HV genes;
    'NB_HIGHLY_VARIABLE': 7000, # None;
    'HVG_ALGOS': HVG_ALGOS[0], # Whether to use triku, or sc HVG;
    'BATCH_KEY': None, # 'batch', # '10X_date', # or None
})

# Pseudo bulk generation
PSEUDO_BULKS = [ 'within_group_replace', 'wo_replace', 'decoupler' ]
PROTO.PSEUDO_BULK = addict.Dict({
    'TYPE': PSEUDO_BULKS[-1],
    'SIZE': 75,
    'REPS': 2,
})

# GLM contrast
PSEUDOBULK_MODES = [ 'sum', 'median' ]
PROTO.GLM = addict.Dict({
    'CONTRAST': 'isHumanTrue.DayX-isHumanFalse.DayX',
    'PSEUDOBULK_MODE': PSEUDOBULK_MODES[0],
    'MIN_CELLS': 30,
    'DEG_THRESHOLD': { 'FDR': .01, 'LFC': 1.5 },
    'FILTER': False,
    'INCLUDE_TIME': True,
})

# GSEA & ORA
REACTOMES = [ 'C2', 'MSigDB' ]
COLLECTIONS = [ 
    'go_biological_process', 'hallmark', 'kegg_pathways',
    'chemical_and_genetic_perturbations', 'biocarta_pathways',
]
PROTO.ENRICHMENT = addict.Dict({
    'GENE_SETS': [
        'KEGG_AXON_GUIDANCE',
        'WP_INTERACTOME_OF_POLYCOMB_REPRESSIVE_COMPLEX_2_PRC2',
        'BENPORATH_PRC2_TARGETS',
        'BIOCARTA_PRC2_PATHWAY',
    ],
    'REACTOME': REACTOMES[1],
    'COLLECTION': COLLECTIONS[0],
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

