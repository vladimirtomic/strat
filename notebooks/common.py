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

COLUMNS_PREPARED = [
    'direction',
    'id',
    'prefix_flank',
    'ins',
    'suffix_flank',
    'prefix_flank_q',
    'ins_q',
    'suffix_flank_q',
]

COLUMNS_PROCESSED = [
    'id',
    'direction',
    'len_ins',
    'len_ins_aln',
    'len_ins_ext',
    'len_ins_ext_aln',
    'ins',
    'ins_aln',
    'ins_ext',
    'ins_ext_aln',
]

COLUMNS_SEQ = [
    'ins',
    'ins_aln',
    'ins_ext',
    'ins_ext_aln',
]

COLUMNS_LEN = [
    'len_ins',
    'len_ins_aln',
    'len_ins_ext',
    'len_ins_ext_aln',
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
    'I': '#FFFFFF',  # white
    ' ': '#FFFFFF',  # white
}


def rev_comp(seq, comps=COMPLEMENT):
    return ''.join(comps.get(n, n) for n in reversed(seq))


def load_tsv(input_path, columns=None):
    if columns:
        header = None
    else:
        header = 0
    df = pd.read_csv(input_path, sep='\t', header=header, quoting=QUOTE_NONE)
    if columns:
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


def prepare_for_plotting(df, col_seq, col_len, col_cnt, col_cov):
    results = []
    for i in range(df[col_len].max()):
        cond = df[col_len] >= i + 1
        row = dict(df[cond][col_seq].str[i].value_counts())
        
        cond = df[col_len] == i + 1
        row[col_cnt] = sum(cond)
        
        cond = df[col_len] >= i + 1
        row[col_cov] = sum(cond)
        
        results.append(row)
    
    result_df = pd.DataFrame(results).fillna(0).astype(int)
    return result_df


def plot(df, col_seq, width, output_path):
    output_image = f'{output_path}{col_seq}.png'
    col_len = 'len_' + col_seq
    col_cnt = 'cnt_' + col_seq
    col_cov = 'cov_' + col_seq
    
    cond = df['direction'] == 'fwd'
    cond &= df[col_len] <= width
    try:
        df_prep_fwd = prepare_for_plotting(df[cond][[col_seq, col_len]], col_seq, col_len, col_cnt, col_cov)
    except Exception:
        df_prep_fwd = None

    cond = df['direction'] == 'rev'
    cond &= df[col_len] <= width
    try:
        df_prep_rev = prepare_for_plotting(df[cond][[col_seq, col_len]], col_seq, col_len, col_cnt, col_cov)
    except Exception:
        df_prep_rev = None

    colors = COLORS
    color_set = colors.keys()
    width = width * 10 + 12
    height = 1002
    half = 500
    image = Image.new('RGB', (width, height), 'grey')
    draw = ImageDraw.Draw(image)
    reach_max = max(
        df_prep_fwd[col_cov].max() if df_prep_fwd is not None else 0,
        df_prep_rev[col_cov].max() if df_prep_rev is not None else 0
    )
    for i in range(width):
        x = i + 1
        N = 'CAG'[i%3]
        colors_ordered = sorted(color_set - set(N)) + [N]

        if df_prep_fwd is not None and i in df_prep_fwd.index:
            row = df_prep_fwd.iloc[i]
            # cnt = row[col_cnt]
            cov = row[col_cov]
            reach = half * cov / reach_max
            bottom = 500
            for j, n in enumerate(colors_ordered):
                if not n in row:
                    continue
                cnt = row[n]
                freq = int((reach * row[n] / cov).round())
                if j == len(colors_ordered) - 1:
                    color = 'black'
                else:
                    color = colors[n]
                draw.line([(10*x, bottom-freq), (10*x, bottom)], width=8, fill=color)
                bottom -= freq
    
        if df_prep_rev is not None and i in df_prep_rev.index:
            row = df_prep_rev.iloc[i]
            # cnt = row[col_cnt]
            cov = row[col_cov]
            reach = half * cov / reach_max
            bottom = 502
            for j, n in enumerate(colors_ordered):
                if not n in row:
                    continue
                cnt = row[n]
                freq = int((reach * row[n] / cov).round())
                if j == len(colors_ordered) - 1:
                    color = 'black'
                else:
                    color = colors[n]
                draw.line([(10*x, bottom), (10*x, bottom+freq)], width=8, fill=color)
                bottom += freq
    
        if i % 3 == 0:
            draw.line([(10*i+5, 0), (10*i+5, height)], width=2, fill='#AAAAAA')
    
        if i % 30 == 0:
            draw.line([(10*i+5, 0), (10*i+5, height)], width=2, fill='white')

        if i % 300 == 0:
            draw.line([(10*i+5, 0), (10*i+5, height)], width=2, fill='black')

    image.save(output_image)


def plot_range(input_path, col_seq, start, stop, output_path):
    df = load_tsv(input_path)
    # df[col_seq] = df.apply(lambda x: rev_comp(x[col_seq]), axis=1)
    col_len = 'len_' + col_seq
    df[col_len] = df[col_len].astype(int)
    cond = df[col_len] >= start
    cond &= df[col_len] < stop
    width = stop - 1
    width = min(width, df[cond][col_len].max())
    if len(df[cond]):
        plot(df[cond], col_seq, width, output_path)
