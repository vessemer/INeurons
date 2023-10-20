from pathlib import Path
import addict
import os

from .protocol import PROTO


# Init Protocol constants
PROTO = PROTO.copy()


# Initialise paths variables
PATHS = addict.Dict({})
PATHS.ROOT = Path(os.environ['HOME'])/'cdata'/'HumanTechnopole'
PATHS.DATA = PATHS.ROOT/'datasets'

PATHS.EMATB = PATHS.DATA/'E-MATB'
PATHS.OPEN = PATHS.DATA/'OPEN'
PATHS.MENDELEY = PATHS.DATA/'MENDELEY'/'multi-lines'

PATHS.CSV = PATHS.DATA/'INTERMEDIATES'
PATHS.LOGS = PATHS.DATA/'LOGS'/PROTO.EXPERIMENT_ID
PATHS.APPENDIX = PATHS.LOGS/'appendix'


# Map E-MATB time signature to days
TIMEPOINT_MAP = {
    'd35': 35,
    'w5': 35,
    'w4': 28,
    'w2': 14,
    'd5': 5,
    'd14': 14,
    '.w4': 28,
    'd35b2': 35,
    '.d5': 5,
    'd35b1': 35,
    'd35b3': 35,
    'd35b4': 35,
    'w17': 119,
    'w18': 126,
    'h6/12': 0,
    'd2': 2,
    'd1': 1,
    'd5': 5,
    'h0': 0,
}


# Map E-MATB time signature to different batches
BATCH_MAP = {
    'd5': 0,
    'd14': 0,
    'd35': 0,
    'h6/12': 0,
    'd2': 0,
    'd1': 0,
    'd5': 0,
    'h0': 0,
    'w5': 1,
    'w4': 1,
    'w2': 1,
    'w17': 1,
    'w18': 1,
    '.w4': 2,
    '.d5': 2,
    'd35b1': 3,
    'd35b2': 4,    
    'd35b3': 5,
    'd35b4': 6,
}



APES = {
    'SandraA', 'JoC', 'ciPS01', 'BmRNA' }
APES = { el.lower() for el in APES }

BONOBO = { 'BmRNA' }
BONOBO = { el.lower() for el in BONOBO }

HUMANS = {
    '409B2', '409B2_CLONE', 'SC102A1', 'HmRNA', 
    'H9', 'S_371', 'S_372', 'S_370', 'S_368'}
HUMANS = { el.lower() for el in HUMANS }


ESC = { 'H9', 'S_371', 'S_372', 'S_370', 'S_368' }
ESC = { el.lower() for el in ESC }


# SandraA, JoC, ciPS01 Shimp 
# BmRNA Bonobo
# 409B2, SC102A, HmRNA IPC Human
# H9 ESC Human