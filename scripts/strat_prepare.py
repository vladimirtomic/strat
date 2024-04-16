import argparse
from datetime import datetime
import gzip
from multiprocessing.dummy import Pool
import regex
from os import listdir
from os import sched_getaffinity
from os.path import isfile, join


DIRECTIONS = ['fwd', 'rev']

COMPLEMENT = {
    'A': 'T',
    'C': 'G',
    'G': 'C',
    'T': 'A'
}


def parse_arguments():
    parser = argparse.ArgumentParser(description='STRAT Prepare - Short Tandem Repeat Analysis Tool - prepare on-target reads from input FASTQ(.GZ) files')
    
    # Define string parameters
    parser.add_argument('--prefix', type=str, required=True, help='Prefix flank sequence of the repeat region (ex. "AGAAAGAAATGGTCTGTGATCCCCC")')
    parser.add_argument('--suffix', type=str, required=True, help='Suffix flank sequence of the repeat region (ex. "CATTCCCGGCTACAAGGACCCTTCG")')
    parser.add_argument('--motif', type=str, required=True, help='Repeat region motif (ex. "CAG")')
    parser.add_argument('--tolerance', type=str, required=True, help='Regex error tolerance configuration (ex. "{e<=5}")')
    parser.add_argument('--cores', type=str, required=True, help='Number of CPU cores to utilize (ex. 2)')
    parser.add_argument('--input_path', type=str, required=True, help='Path to directory containing input FASTQ(.GZ) files (ex. "/data/fastqs/")')
    parser.add_argument('--output_path', type=str, required=True, help='Path to directory to store output files (ex. "/data/outputs/")')
    
    return parser


def rev_comp(seq, comps):
    return ''.join(comps.get(n, n) for n in reversed(seq))


def chop(string, start, end):
    prefix_flank = string[start-10:start]
    suffix_flank = string[end:end+10]
    ins = string[start:end]
    return prefix_flank, ins, suffix_flank


class ReadProcessor:
    def __init__(
        self,
        prefix,
        suffix,
        tolerance
    ):
        self.prefix = {
            'fwd': f'({prefix})' + tolerance,
            'rev': f'({rev_comp(suffix, COMPLEMENT)})' + tolerance
        }
        self.suffix = {
            'fwd': f'({suffix})' + tolerance,
            'rev': f'({rev_comp(prefix, COMPLEMENT)})' + tolerance
        }

    def process_read(self, id, seq, opt, qual):
        row = ['offtarget', id, seq, qual]
        ress = []
        for k in DIRECTIONS:
            prefixes = regex.findall(self.prefix[k], seq)
            suffixes = regex.findall(self.suffix[k], seq)
            if len(prefixes) == 1 and len(suffixes) == 1:
                start = seq.index(prefixes[0]) + len(prefixes[0])
                end = seq.index(suffixes[0])
                if end > start:
                    direction = k
                    prefix_flank, ins, suffix_flank = chop(seq, start, end)
                    prefix_flank_q, ins_q, suffix_flank_q = chop(qual, start, end)

                    ress.append([
                        'ontarget',
                        direction,
                        id,
                        prefix_flank, ins, suffix_flank,
                        prefix_flank_q, ins_q, suffix_flank_q
                    ])

        if not ress:
            ress.append(row)

        return ress


def process_fastq(fastq_path, output_path, read_processor):
    print(f'{datetime.now()} - {fastq_path}')
    fastq_name = fastq_path.split('/')[-1]
    gzipped = fastq_name.endswith('.gz')
    ontarget_output_path = f'{output_path}{fastq_name}.ontarget.tsv'
    # offtarget_output_path = f'{output_path}{fastq_name}.offtarget.tsv'
    openner = gzip.open if gzipped else open

    with openner(fastq_path, 'rt') as f, open(ontarget_output_path, 'wt') as o:
        for i, line in enumerate(f):
            line = line.strip()
            if i%4 == 0:
                if line.startswith('@'):
                    id = line.split(' ')[0]
                else:
                    print(f'Error in {fastq_path} line {i} - not an ID line')
                    raise
            elif i%4 == 1:
                seq = line
            elif i%4 == 2:
                if line.startswith('+'):
                    opt = line
                else:
                    print(f'Error in {fastq_path} line {i} - not a + line')
                    raise
            elif i%4 == 3:
                qual = line
                ress = read_processor.process_read(id, seq, opt, qual)
                for res in ress:
                    if res[0] == 'ontarget':
                        o.write('\t'.join(res[1:]) + '\n')
                    else:
                        # g.write('\t'.join(res[1:]) + '\n')
                        pass


def main():
    # Parse the command-line arguments
    args = parse_arguments().parse_args()

    # Access the values of the parameters
    prefix = args.prefix
    suffix = args.suffix
    motif = args.motif
    tolerance = args.tolerance
    cores = args.cores
    input_path = args.input_path
    output_path = args.output_path

    # Log provided parameters
    print(f'prefix: {prefix}')
    print(f'suffix: {suffix}')
    print(f'motif: {motif}')
    print(f'tolerance: {tolerance}')
    print(f'cores: {cores}')
    print(f'input_path: {input_path}')
    print(f'output_path: {output_path}')

    fastq_paths = sorted(join(input_path, f) for f in listdir(input_path) if 'fastq' in f and isfile(join(input_path, f)))
    # fastq_paths = fastq_paths[:1]
    print([len(fastq_paths), fastq_paths[0], fastq_paths[-1]])


    print(f'{datetime.now()} - STRAT Prepare - Start')

    inputs = [(fastq_path, output_path, ReadProcessor(prefix, suffix, tolerance)) for fastq_path in fastq_paths]

    with Pool(2) as p:
        p.starmap(process_fastq, inputs)

    print(f'{datetime.now()} - STRAT Prepare - End')

if __name__ == "__main__":
    main()
