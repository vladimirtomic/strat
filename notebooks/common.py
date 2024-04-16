from csv import QUOTE_NONE
from collections import Counter
from datetime import datetime
import matplotlib.pyplot as plt
import matplotlib.ticker as plticker
import numpy as np
import pandas as pd
from PIL import Image, ImageDraw
import seaborn as sns
from string2string.alignment import NeedlemanWunsch


DIRECTIONS = ['fwd', 'rev']

COMPLEMENT = {
    'A': 'T',
    'C': 'G',
    'G': 'C',
    'T': 'A'
}

COLUMNS = [
    'direction',
    'id',
    'prefix_flank',
    'ins',
    'suffix_flank',
    'prefix_flank_q',
    'ins_q',
    'suffix_flank_q',
]

MATCH_WEIGHT = 10
MISMATCH_WEIGHT = -8
GAP_WEIGHT = -9
NW = NeedlemanWunsch(
    match_weight=MATCH_WEIGHT,
    mismatch_weight=MISMATCH_WEIGHT,
    gap_weight=GAP_WEIGHT,
    gap_char=''
)

COLORS = {
    'A': '#3DA853',  # green
    'C': '#4285F4',  # blue
    'G': '#F8BC07',  # yellow
    'T': '#EA4334',  # red
    ' ': 'white'
}


def rev_comp(seq, comps):
    return ''.join(comps.get(n, n) for n in reversed(seq))


def load(input_path, columns):
    df = pd.read_csv(input_path, sep='\t', header=None, dtype=str, quoting=QUOTE_NONE)
    df.columns = columns

    return df


def lengths(df, columns_seq, columns_len):
    for s, l in zip(columns_seq, columns_len):
        df[l] = df[s].str.len()
    return df


def freqs(seq, targets):
    counts = Counter(seq)
    res = {}
    for target in sorted(targets):
        res[target] = 0 if target not in counts else counts[target]
    for count in counts:
        if count not in targets:
            raise Exception(f'Detected {count} in input seq {seq}!')
    return res
