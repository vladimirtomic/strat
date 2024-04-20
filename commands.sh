# v1.0

# STRAT prepare

docker run \
--rm \
-it \
-v /Users/vladimirtomic/Documents/repos/strat/:/opt/repos/strat/ \
-v /Users/vladimirtomic/Documents/projects/ONT/data/:/opt/data/ \
strat:latest

# 2h 15min
python /opt/repos/strat/scripts/strat_prepare.py \
--prefix 'AGAAAGAAATGGTCTGTGATCCCCC' \
--suffix 'CATTCCCGGCTACAAGGACCCTTCG' \
--motif 'CAG' \
--tolerance '{e<=5}' \
--cores 2 \
--input_path '/opt/data/pcr2persons/fastq/guppy/' \
--output_path '/opt/data/pcr2persons/output/guppy/'

cat /opt/data/pcr2persons/output/guppy/*.ontarget.tsv 1> /opt/data/workdir/pcr2persons.guppy.ontarget.tsv

# 2h 15min
python /opt/repos/strat/scripts/strat_prepare.py \
--prefix 'AGAAAGAAATGGTCTGTGATCCCCC' \
--suffix 'CATTCCCGGCTACAAGGACCCTTCG' \
--motif 'CAG' \
--tolerance '{e<=5}' \
--cores 2 \
--input_path '/opt/data/pcr2persons/fastq/dorado/' \
--output_path '/opt/data/pcr2persons/output/dorado/'

cat /opt/data/pcr2persons/output/dorado/*.ontarget.tsv 1> /opt/data/workdir/pcr2persons.dorado.ontarget.tsv

# 6h 10min
python /opt/repos/strat/scripts/strat_prepare.py \
--prefix 'AGAAAGAAATGGTCTGTGATCCCCC' \
--suffix 'CATTCCCGGCTACAAGGACCCTTCG' \
--motif 'CAG' \
--tolerance '{e<=5}' \
--cores 2 \
--input_path '/opt/data/jovan/fastq/guppy/' \
--output_path '/opt/data/jovan/output/guppy/'

cat /opt/data/jovan/output/guppy/*.ontarget.tsv 1> /opt/data/workdir/jovan.guppy.ontarget.tsv

# 3h 15min
python /opt/repos/strat/scripts/strat_prepare.py \
--prefix 'AGAAAGAAATGGTCTGTGATCCCCC' \
--suffix 'CATTCCCGGCTACAAGGACCCTTCG' \
--motif 'CAG' \
--tolerance '{e<=5}' \
--cores 2 \
--input_path '/opt/data/jovan/fastq/dorado/' \
--output_path '/opt/data/jovan/output/dorado/'

cat /opt/data/jovan/output/dorado/*.ontarget.tsv 1> /opt/data/workdir/jovan.dorado.ontarget.tsv

# 2h 43min
python /opt/repos/strat/scripts/strat_prepare.py \
--prefix 'AGAAAGAAATGGTCTGTGATCCCCC' \
--suffix 'CATTCCCGGCTACAAGGACCCTTCG' \
--motif 'CAG' \
--tolerance '{e<=5}' \
--cores 2 \
--input_path '/opt/data/dm108/fastq/guppy/' \
--output_path '/opt/data/dm108/output/guppy/'

cat /opt/data/dm108/output/guppy/*.ontarget.tsv 1> /opt/data/workdir/dm108.guppy.ontarget.tsv

# 1h 46min
python /opt/repos/strat/scripts/strat_prepare.py \
--prefix 'AGAAAGAAATGGTCTGTGATCCCCC' \
--suffix 'CATTCCCGGCTACAAGGACCCTTCG' \
--motif 'CAG' \
--tolerance '{e<=5}' \
--cores 2 \
--input_path '/opt/data/bc3_1/fastq/guppy/' \
--output_path '/opt/data/bc3_1/output/guppy/'

cat /opt/data/bc3_1/output/guppy/*.ontarget.tsv 1> /opt/data/workdir/bc3_1.guppy.ontarget.tsv

# 1h 17min
python /opt/repos/strat/scripts/strat_prepare.py \
--prefix 'AGAAAGAAATGGTCTGTGATCCCCC' \
--suffix 'CATTCCCGGCTACAAGGACCCTTCG' \
--motif 'CAG' \
--tolerance '{e<=5}' \
--cores 2 \
--input_path '/opt/data/bc3_1/fastq/dorado/' \
--output_path '/opt/data/bc3_1/output/dorado/'

cat /opt/data/bc3_1/output/dorado/*.ontarget.tsv 1> /opt/data/workdir/bc3_1.dorado.ontarget.tsv

# 3h 18min
python /opt/repos/strat/scripts/strat_prepare.py \
--prefix 'AGAAAGAAATGGTCTGTGATCCCCC' \
--suffix 'CATTCCCGGCTACAAGGACCCTTCG' \
--motif 'CAG' \
--tolerance '{e<=5}' \
--cores 2 \
--input_path '/opt/data/bc3_2/fastq/guppy/' \
--output_path '/opt/data/bc3_2/output/guppy/'

cat /opt/data/bc3_2/output/guppy/*.ontarget.tsv 1> /opt/data/workdir/bc3_2.guppy.ontarget.tsv

# 0h 4min
python /opt/repos/strat/scripts/strat_prepare.py \
--prefix 'AGAAAGAAATGGTCTGTGATCCCCC' \
--suffix 'CATTCCCGGCTACAAGGACCCTTCG' \
--motif 'CAG' \
--tolerance '{e<=5}' \
--cores 2 \
--input_path '/opt/data/bc3_2/fastq/dorado/' \
--output_path '/opt/data/bc3_2/output/dorado/'

cat /opt/data/bc3_2/output/dorado/*.ontarget.tsv 1> /opt/data/workdir/bc3_2.dorado.ontarget.tsv

# 3h 15min
python /opt/repos/strat/scripts/strat_prepare.py \
--prefix 'AGAAAGAAATGGTCTGTGATCCCCC' \
--suffix 'CATTCCCGGCTACAAGGACCCTTCG' \
--motif 'CAG' \
--tolerance '{e<=5}' \
--cores 2 \
--input_path '/opt/data/bc3_3/fastq/guppy/' \
--output_path '/opt/data/bc3_3/output/guppy/'

cat /opt/data/bc3_3/output/guppy/*.ontarget.tsv 1> /opt/data/workdir/bc3_3.guppy.ontarget.tsv

# 0h 19min
python /opt/repos/strat/scripts/strat_prepare.py \
--prefix 'AGAAAGAAATGGTCTGTGATCCCCC' \
--suffix 'CATTCCCGGCTACAAGGACCCTTCG' \
--motif 'CAG' \
--tolerance '{e<=5}' \
--cores 2 \
--input_path '/opt/data/bc3_3/fastq/dorado/' \
--output_path '/opt/data/bc3_3/output/dorado/'

cat /opt/data/bc3_3/output/dorado/*.ontarget.tsv 1> /opt/data/workdir/bc3_3.dorado.ontarget.tsv

# 1h 31min
python /opt/repos/strat/scripts/strat_prepare.py \
--prefix 'AGAAAGAAATGGTCTGTGATCCCCC' \
--suffix 'CATTCCCGGCTACAAGGACCCTTCG' \
--motif 'CAG' \
--tolerance '{e<=5}' \
--cores 2 \
--input_path '/opt/data/bc6_1/fastq/guppy/' \
--output_path '/opt/data/bc6_1/output/guppy/'

cat /opt/data/bc6_1/output/guppy/*.ontarget.tsv 1> /opt/data/workdir/bc6_1.guppy.ontarget.tsv

# 1h 32min
python /opt/repos/strat/scripts/strat_prepare.py \
--prefix 'AGAAAGAAATGGTCTGTGATCCCCC' \
--suffix 'CATTCCCGGCTACAAGGACCCTTCG' \
--motif 'CAG' \
--tolerance '{e<=5}' \
--cores 2 \
--input_path '/opt/data/bc6_2/fastq/guppy/' \
--output_path '/opt/data/bc6_2/output/guppy/'

cat /opt/data/bc6_2/output/guppy/*.ontarget.tsv 1> /opt/data/workdir/bc6_2.guppy.ontarget.tsv

# 1h 31min
python /opt/repos/strat/scripts/strat_prepare.py \
--prefix 'AGAAAGAAATGGTCTGTGATCCCCC' \
--suffix 'CATTCCCGGCTACAAGGACCCTTCG' \
--motif 'CAG' \
--tolerance '{e<=5}' \
--cores 2 \
--input_path '/opt/data/bc6_3/fastq/guppy/' \
--output_path '/opt/data/bc6_3/output/guppy/'

cat /opt/data/bc6_3/output/guppy/*.ontarget.tsv 1> /opt/data/workdir/bc6_3.guppy.ontarget.tsv

# 1h 27min
python /opt/repos/strat/scripts/strat_prepare.py \
--prefix 'AGAAAGAAATGGTCTGTGATCCCCC' \
--suffix 'CATTCCCGGCTACAAGGACCCTTCG' \
--motif 'CAG' \
--tolerance '{e<=5}' \
--cores 2 \
--input_path '/opt/data/bc6_4/fastq/guppy/' \
--output_path '/opt/data/bc6_4/output/guppy/'

cat /opt/data/bc6_4/output/guppy/*.ontarget.tsv 1> /opt/data/workdir/bc6_4.guppy.ontarget.tsv

# 0h 54min
python /opt/repos/strat/scripts/strat_prepare.py \
--prefix 'AGAAAGAAATGGTCTGTGATCCCCC' \
--suffix 'CATTCCCGGCTACAAGGACCCTTCG' \
--motif 'CAG' \
--tolerance '{e<=5}' \
--cores 2 \
--input_path '/opt/data/bc6_5/fastq/guppy/' \
--output_path '/opt/data/bc6_5/output/guppy/'

cat /opt/data/bc6_5/output/guppy/*.ontarget.tsv 1> /opt/data/workdir/bc6_5.guppy.ontarget.tsv

# 0h 57min
python /opt/repos/strat/scripts/strat_prepare.py \
--prefix 'AGAAAGAAATGGTCTGTGATCCCCC' \
--suffix 'CATTCCCGGCTACAAGGACCCTTCG' \
--motif 'CAG' \
--tolerance '{e<=5}' \
--cores 2 \
--input_path '/opt/data/bc6_6/fastq/guppy/' \
--output_path '/opt/data/bc6_6/output/guppy/'

cat /opt/data/bc6_6/output/guppy/*.ontarget.tsv 1> /opt/data/workdir/bc6_6.guppy.ontarget.tsv

# STRAT process

python /opt/repos/strat/scripts/strat_process.py \
--motif 'CAG' \
--threshold '100' \
--input_path '/opt/data/pcr2persons/output/guppy/guppy.ontarget.tsv' \
--output_path '/opt/data/pcr2persons/output/guppy/'

python /opt/repos/strat/scripts/strat_process.py \
--motif 'CAG' \
--threshold '100' \
--input_path '/opt/data/pcr2persons/output/dorado/dorado.ontarget.tsv' \
--output_path '/opt/data/pcr2persons/output/dorado/'

python /opt/repos/strat/scripts/strat_process.py \
--motif 'CAG' \
--threshold '100' \
--input_path '/opt/data/jovan/output/guppy/guppy.ontarget.tsv' \
--output_path '/opt/data/jovan/output/guppy/'

python /opt/repos/strat/scripts/strat_process.py \
--motif 'CAG' \
--threshold '100' \
--input_path '/opt/data/jovan/output/dorado/dorado.ontarget.tsv' \
--output_path '/opt/data/jovan/output/dorado/'

python /opt/repos/strat/scripts/strat_process.py \
--motif 'CAG' \
--threshold '10' \
--input_path '/opt/data/dm108/output/guppy/guppy.ontarget.tsv' \
--output_path '/opt/data/dm108/output/guppy/'

python /opt/repos/strat/scripts/strat_process.py \
--motif 'CAG' \
--threshold '100' \
--input_path '/opt/data/bc3_1/output/guppy/guppy.ontarget.tsv' \
--output_path '/opt/data/bc3_1/output/guppy/'

python /opt/repos/strat/scripts/strat_process.py \
--motif 'CAG' \
--threshold '10' \
--input_path '/opt/data/bc3_1/output/dorado/dorado.ontarget.tsv' \
--output_path '/opt/data/bc3_1/output/dorado/'

python /opt/repos/strat/scripts/strat_process.py \
--motif 'CAG' \
--threshold '100' \
--input_path '/opt/data/bc3_2/output/guppy/guppy.ontarget.tsv' \
--output_path '/opt/data/bc3_2/output/guppy/'

python /opt/repos/strat/scripts/strat_process.py \
--motif 'CAG' \
--threshold '10' \
--input_path '/opt/data/bc3_2/output/dorado/dorado.ontarget.tsv' \
--output_path '/opt/data/bc3_2/output/dorado/'

python /opt/repos/strat/scripts/strat_process.py \
--motif 'CAG' \
--threshold '100' \
--input_path '/opt/data/bc3_3/output/guppy/guppy.ontarget.tsv' \
--output_path '/opt/data/bc3_3/output/guppy/'

python /opt/repos/strat/scripts/strat_process.py \
--motif 'CAG' \
--threshold '10' \
--input_path '/opt/data/bc3_3/output/dorado/dorado.ontarget.tsv' \
--output_path '/opt/data/bc3_3/output/dorado/'

# STRAT process new

python /opt/repos/strat/scripts/strat_process.py \
--motif 'CAG' \
--threshold '1' \
--input_path '/opt/data/workdir/pcr2persons.guppy.ontarget.tsv' \
--output_path '/opt/data/workdir/' \
1> /opt/data/workdir/pcr2persons.guppy.std.out \
2> /opt/data/workdir/pcr2persons.guppy.std.err

python /opt/repos/strat/scripts/strat_process.py \
--motif 'CAG' \
--threshold '1' \
--input_path '/opt/data/workdir/pcr2persons.dorado.ontarget.tsv' \
--output_path '/opt/data/workdir/' \
1> /opt/data/workdir/pcr2persons.dorado.std.out \
2> /opt/data/workdir/pcr2persons.dorado.std.err

python /opt/repos/strat/scripts/strat_process.py \
--motif 'CAG' \
--threshold '1' \
--input_path '/opt/data/workdir/jovan.guppy.ontarget.tsv' \
--output_path '/opt/data/workdir/' \
1> /opt/data/workdir/jovan.guppy.std.out \
2> /opt/data/workdir/jovan.guppy.std.err

python /opt/repos/strat/scripts/strat_process.py \
--motif 'CAG' \
--threshold '1' \
--input_path '/opt/data/workdir/jovan.dorado.ontarget.tsv' \
--output_path '/opt/data/workdir/' \
1> /opt/data/workdir/jovan.dorado.std.out \
2> /opt/data/workdir/jovan.dorado.std.err

python /opt/repos/strat/scripts/strat_process.py \
--motif 'CAG' \
--threshold '1' \
--input_path '/opt/data/workdir/dm108.guppy.ontarget.tsv' \
--output_path '/opt/data/workdir/' \
1> /opt/data/workdir/dm108.guppy.std.out \
2> /opt/data/workdir/dm108.guppy.std.err

python /opt/repos/strat/scripts/strat_process.py \
--motif 'CAG' \
--threshold '1' \
--input_path '/opt/data/workdir/bc3_1.guppy.ontarget.tsv' \
--output_path '/opt/data/workdir/' \
1> /opt/data/workdir/bc3_1.guppy.std.out \
2> /opt/data/workdir/bc3_1.guppy.std.err

python /opt/repos/strat/scripts/strat_process.py \
--motif 'CAG' \
--threshold '1' \
--input_path '/opt/data/workdir/bc3_1.dorado.ontarget.tsv' \
--output_path '/opt/data/workdir/' \
1> /opt/data/workdir/bc3_1.dorado.std.out \
2> /opt/data/workdir/bc3_1.dorado.std.err

python /opt/repos/strat/scripts/strat_process.py \
--motif 'CAG' \
--threshold '1' \
--input_path '/opt/data/workdir/bc3_2.guppy.ontarget.tsv' \
--output_path '/opt/data/workdir/' \
1> /opt/data/workdir/bc3_2.guppy.std.out \
2> /opt/data/workdir/bc3_2.guppy.std.err

python /opt/repos/strat/scripts/strat_process.py \
--motif 'CAG' \
--threshold '1' \
--input_path '/opt/data/workdir/bc3_2.dorado.ontarget.tsv' \
--output_path '/opt/data/workdir/' \
1> /opt/data/workdir/bc3_2.dorado.std.out \
2> /opt/data/workdir/bc3_2.dorado.std.err

python /opt/repos/strat/scripts/strat_process.py \
--motif 'CAG' \
--threshold '1' \
--input_path '/opt/data/workdir/bc3_3.guppy.ontarget.tsv' \
--output_path '/opt/data/workdir/' \
1> /opt/data/workdir/bc3_3.guppy.std.out \
2> /opt/data/workdir/bc3_3.guppy.std.err

python /opt/repos/strat/scripts/strat_process.py \
--motif 'CAG' \
--threshold '1' \
--input_path '/opt/data/workdir/bc3_3.dorado.ontarget.tsv' \
--output_path '/opt/data/workdir/' \
1> /opt/data/workdir/bc3_3.dorado.std.out \
2> /opt/data/workdir/bc3_3.dorado.std.err

python /opt/repos/strat/scripts/strat_process.py \
--motif 'CAG' \
--threshold '1' \
--input_path '/opt/data/workdir/bc6_1.guppy.ontarget.tsv' \
--output_path '/opt/data/workdir/' \
1> /opt/data/workdir/bc6_1.guppy.std.out \
2> /opt/data/workdir/bc6_1.guppy.std.err

python /opt/repos/strat/scripts/strat_process.py \
--motif 'CAG' \
--threshold '1' \
--input_path '/opt/data/workdir/bc6_2.guppy.ontarget.tsv' \
--output_path '/opt/data/workdir/' \
1> /opt/data/workdir/bc6_2.guppy.std.out \
2> /opt/data/workdir/bc6_2.guppy.std.err

python /opt/repos/strat/scripts/strat_process.py \
--motif 'CAG' \
--threshold '1' \
--input_path '/opt/data/workdir/bc6_3.guppy.ontarget.tsv' \
--output_path '/opt/data/workdir/' \
1> /opt/data/workdir/bc6_3.guppy.std.out \
2> /opt/data/workdir/bc6_3.guppy.std.err

python /opt/repos/strat/scripts/strat_process.py \
--motif 'CAG' \
--threshold '1' \
--input_path '/opt/data/workdir/bc6_4.guppy.ontarget.tsv' \
--output_path '/opt/data/workdir/' \
1> /opt/data/workdir/bc6_4.guppy.std.out \
2> /opt/data/workdir/bc6_4.guppy.std.err

python /opt/repos/strat/scripts/strat_process.py \
--motif 'CAG' \
--threshold '1' \
--input_path '/opt/data/workdir/bc6_5.guppy.ontarget.tsv' \
--output_path '/opt/data/workdir/' \
1> /opt/data/workdir/bc6_5.guppy.std.out \
2> /opt/data/workdir/bc6_5.guppy.std.err

python /opt/repos/strat/scripts/strat_process.py \
--motif 'CAG' \
--threshold '1' \
--input_path '/opt/data/workdir/bc6_6.guppy.ontarget.tsv' \
--output_path '/opt/data/workdir/' \
1> /opt/data/workdir/bc6_6.guppy.std.out \
2> /opt/data/workdir/bc6_6.guppy.std.err
