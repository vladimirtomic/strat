import argparse
from csv import QUOTE_NONE
from datetime import datetime
import matplotlib.pyplot as plt
import matplotlib.ticker as plticker
import pandas as pd
import seaborn as sns


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

COLUMNS_SEQ = [
    'prefix_flank',
    'ins',
    'suffix_flank',
]

COLUMNS_LEN = ['len_' + s for s in COLUMNS_SEQ]

COLUMNS_SEQ_EXT = ['ins_ext']
COLUMNS_LEN_EXT = ['len_ins_ext']

def parse_arguments():
    parser = argparse.ArgumentParser(description='STRAT Process - Short Tandem Repeat Analysis Tool - process on-target inserts collected by STRAT Prepare')
    
    # Define string parameters
    parser.add_argument('--motif', type=str, required=True, help='Repeat region motif (ex. "CAG")')
    parser.add_argument('--threshold', type=int, required=True, help='Minimum number of inserts of each size (ex. "100")')
    parser.add_argument('--input_path', type=str, required=True, help='Path to TSV file containing output of STRAT Prepare (ex. "/data/fastqs/")')
    parser.add_argument('--output_path', type=str, required=True, help='Path to directory to store output files (ex. "/data/outputs/")')
    
    return parser


def rev_comp(seq, comps):
    return ''.join(comps.get(n, n) for n in reversed(seq))


def load(input_path, columns):
    df = pd.read_csv(input_path, sep='\t', header=None, dtype=str, quoting=QUOTE_NONE)
    df.columns = columns

    return df


def extend_ins(row):
    prefix = row['prefix_flank']
    ins = row['ins']
    suffix = row['suffix_flank']
    motif = MOTIFS[row['direction']]
    target = 2 * motif
    window = len(target)
    for s in range(window - 1):
        # i in (window-1)..1
        i = window - 1 - s
        if prefix[-i:] + ins[:window - i] == target:
            ins = prefix[-i:] + ins
            break
    
    for s in range(window - 1):
        # i in 1..(window-1)
        i = s + 1
        if ins[-i:] + suffix[:window - i] == target:
            ins = ins + suffix[:window - i]
            break

    return ins


def lengths(df, columns_seq, columns_len):
    for s, l in zip(columns_seq, columns_len):
        df[l] = df[s].str.len()
        df[l + '_adj'] = (df[l] / 3).round().astype(int)
    return df


def get_abundant_lengths(df, column_seq, column_len, threshold):
    dfg = df.groupby(column_len)[column_seq].count().reset_index()

    cond = dfg[column_seq] > threshold
    dfg = dfg[cond]
    return dfg


def consensus_string(strings):
    if not strings or not all(len(strings[0]) == len(s) for s in strings):
        raise ValueError("Input strings must be non-empty and of equal length")

    consensus = ''
    for i in range(len(strings[0])):
        # Create a dictionary to count occurrences of each character at the current position
        char_count = {}
        for s in strings:
            char = s[i]
            char_count[char] = char_count.get(char, 0) + 1

        # Find the most frequent character at the current position
        most_frequent_char = max(char_count, key=char_count.get)

        # Append the most frequent character to the consensus string
        consensus += most_frequent_char

    return consensus


def get_consensus_strings(df, abundant_lengths, column_len, column_seq, directions):
    consensi = []
    for l in sorted(abundant_lengths):
        for direction in directions:
            cond = df[column_len] == l
            cond &= df['direction'] == direction
            strings = list(df[cond][column_seq])
            if strings:
                consensi.append([direction, l, len(strings), consensus_string(strings)])
    
    df_consensus = pd.DataFrame(consensi, columns=['direction', column_len, 'count', column_seq])
    return df_consensus


def plot_histogram(df, x, hue, output_histogram):
    fig, ax = plt.subplots(figsize=(16, 10))
    gfg = sns.histplot(df, x=x, discrete=True, hue=hue)
    # gfg.set_xlim(0, 1000)
    gfg.set_yscale("log")
    loc = plticker.MultipleLocator(base=5)
    gfg.xaxis.set_major_locator(loc)
    gfg.set_xticklabels(gfg.get_xticklabels(), rotation=90)
    fig.savefig(output_histogram)


def main():
    # Parse the command-line arguments
    args = parse_arguments().parse_args()

    # Access the values of the parameters
    motif = args.motif
    threshold = args.threshold
    input_path = args.input_path
    output_path = args.output_path

    # Generate reverse complement motif
    global MOTIFS
    MOTIFS = {
        'fwd': motif,
        'rev': rev_comp(motif, COMPLEMENT)
    }

    # Log provided parameters
    print(f'motif: {motif}')
    print(f'threshold: {threshold}')
    print(f'input_path: {input_path}')
    print(f'output_path: {output_path}')

    print(f'{datetime.now()} - STRAT Process - Start')

    # Load on-target inserts
    df = load(input_path, COLUMNS)
    print(f'{datetime.now()} - STRAT Process - Loaded {len(df)} rows')
    
    # Extend on-target inserts
    df['ins_ext'] = df.apply(extend_ins, axis=1)
    extended = sum(df['ins'] != df['ins_ext'])
    print(f'{datetime.now()} - STRAT Process - Extended {extended} inserts')
    
    # Calculate lengths of inserts
    df = lengths(df, COLUMNS_SEQ + COLUMNS_SEQ_EXT, COLUMNS_LEN + COLUMNS_LEN_EXT)

    # Keep abundant insert lengths
    dfg = get_abundant_lengths(df, 'ins_ext', 'len_ins_ext', threshold)
    abundant_lengths = set(dfg['len_ins_ext'])

    cond = df['len_ins_ext'].isin(abundant_lengths)
    df = df[cond]
    print(f'{datetime.now()} - STRAT Process - Kept {len(df)} abundant insert length rows')

    # Generate consensus strings per extended insert size
    dfc = get_consensus_strings(df, abundant_lengths, 'len_ins_ext', 'ins_ext', DIRECTIONS)
    output_consensus = f'{output_path}consensus.ontarget.ext.tsv'
    dfc.to_csv(output_consensus, index=False, sep='\t')
    print(f'{datetime.now()} - STRAT Process - Written {len(dfc)} consensus inserts to {output_consensus}')

    # Plot histogram
    output_histogram = f'{output_path}inserts.ontarget.ext.png'
    plot_histogram(df.sort_values('direction'), 'len_ins_ext_adj', 'direction', output_histogram)
    print(f'{datetime.now()} - STRAT Process - Plotted insert length histogram to {output_histogram}')

    print(f'{datetime.now()} - STRAT Process - End')


if __name__ == "__main__":
    main()
