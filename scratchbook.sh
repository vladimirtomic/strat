# ----------
# 2024-07-13
# ----------

docker run \
--rm \
-it \
-v /Users/vladimirtomic/Documents/repos/strat/:/opt/repos/strat/ \
-v /Users/vladimirtomic/Documents/projects/ONT/data/:/opt/data/ \
-v /Users/vladimirtomic/Documents/projects/ONT/data_out/:/opt/data_out/ \
strat:latest

# ----------
# 2024-06-15
# ----------

docker run \
--rm \
-it \
-p 8888:8888 \
-v /Users/vladimirtomic/Documents/repos/strat/:/opt/repos/strat/ \
-v /Users/vladimirtomic/Documents/projects/ONT/data/:/opt/data/ \
-v /Users/vladimirtomic/Documents/projects/ONT/data_out/:/opt/data_out/ \
strat:latest

# Run notebooks inside the docker container
jupyter notebook --ip 0.0.0.0 --port 8888 --no-browser --allow-root

# ----------
# 2024-05-18
# ----------

docker run \
--rm \
-it \
-v /Users/vladimirtomic/Documents/repos/strat/:/opt/repos/strat/ \
-v /Users/vladimirtomic/Documents/projects/ONT/data/:/opt/data/ \
strat:latest

# 21min
python /opt/repos/strat/scripts/strat_prepare.py \
--prefix 'TAAGATAATATATTTTTAAAAAATG' \
--suffix 'TAAAGCCAGGTTTTCTAACATGAAG' \
--motif 'CAG' \
--tolerance '{e<=5}' \
--cores 2 \
--input_path '/opt/data/sca8_12/fastq/guppy/' \
--output_path '/opt/data/sca8_12/output/guppy/'

cat /opt/data/sca8_12/output/guppy/*.ontarget.tsv 1> /opt/data/workdir/sca8_12.guppy.ontarget.tsv

python /opt/repos/strat/scripts/strat_process.py \
--motif 'CAG' \
--threshold '1' \
--input_path '/opt/data/workdir/sca8_12.guppy.ontarget.tsv' \
--output_path '/opt/data/workdir/' \
1> /opt/data/workdir/sca8_12.guppy.std.out \
2> /opt/data/workdir/sca8_12.guppy.std.err

# 18min
python /opt/repos/strat/scripts/strat_prepare.py \
--prefix 'TAAGATAATATATTTTTAAAAAATG' \
--suffix 'TAAAGCCAGGTTTTCTAACATGAAG' \
--motif 'CAG' \
--tolerance '{e<=5}' \
--cores 2 \
--input_path '/opt/data/sca8_15/fastq/guppy/' \
--output_path '/opt/data/sca8_15/output/guppy/'

cat /opt/data/sca8_15/output/guppy/*.ontarget.tsv 1> /opt/data/workdir/sca8_15.guppy.ontarget.tsv

python /opt/repos/strat/scripts/strat_process.py \
--motif 'CAG' \
--threshold '1' \
--input_path '/opt/data/workdir/sca8_15.guppy.ontarget.tsv' \
--output_path '/opt/data/workdir/' \
1> /opt/data/workdir/sca8_15.guppy.std.out \
2> /opt/data/workdir/sca8_15.guppy.std.err

# 25min
python /opt/repos/strat/scripts/strat_prepare.py \
--prefix 'TAAGATAATATATTTTTAAAAAATG' \
--suffix 'TAAAGCCAGGTTTTCTAACATGAAG' \
--motif 'CAG' \
--tolerance '{e<=5}' \
--cores 2 \
--input_path '/opt/data/sca8_16/fastq/guppy/' \
--output_path '/opt/data/sca8_16/output/guppy/'

cat /opt/data/sca8_16/output/guppy/*.ontarget.tsv 1> /opt/data/workdir/sca8_16.guppy.ontarget.tsv

python /opt/repos/strat/scripts/strat_process.py \
--motif 'CAG' \
--threshold '1' \
--input_path '/opt/data/workdir/sca8_16.guppy.ontarget.tsv' \
--output_path '/opt/data/workdir/' \
1> /opt/data/workdir/sca8_16.guppy.std.out \
2> /opt/data/workdir/sca8_16.guppy.std.err

# 19min
python /opt/repos/strat/scripts/strat_prepare.py \
--prefix 'TAAGATAATATATTTTTAAAAAATG' \
--suffix 'TAAAGCCAGGTTTTCTAACATGAAG' \
--motif 'CAG' \
--tolerance '{e<=5}' \
--cores 2 \
--input_path '/opt/data/sca8_19/fastq/guppy/' \
--output_path '/opt/data/sca8_19/output/guppy/'

cat /opt/data/sca8_19/output/guppy/*.ontarget.tsv 1> /opt/data/workdir/sca8_19.guppy.ontarget.tsv

python /opt/repos/strat/scripts/strat_process.py \
--motif 'CAG' \
--threshold '1' \
--input_path '/opt/data/workdir/sca8_19.guppy.ontarget.tsv' \
--output_path '/opt/data/workdir/' \
1> /opt/data/workdir/sca8_19.guppy.std.out \
2> /opt/data/workdir/sca8_19.guppy.std.err

# 13min
python /opt/repos/strat/scripts/strat_prepare.py \
--prefix 'TAAGATAATATATTTTTAAAAAATG' \
--suffix 'TAAAGCCAGGTTTTCTAACATGAAG' \
--motif 'CAG' \
--tolerance '{e<=5}' \
--cores 2 \
--input_path '/opt/data/sca8_20/fastq/guppy/' \
--output_path '/opt/data/sca8_20/output/guppy/'

cat /opt/data/sca8_20/output/guppy/*.ontarget.tsv 1> /opt/data/workdir/sca8_20.guppy.ontarget.tsv

python /opt/repos/strat/scripts/strat_process.py \
--motif 'CAG' \
--threshold '1' \
--input_path '/opt/data/workdir/sca8_20.guppy.ontarget.tsv' \
--output_path '/opt/data/workdir/' \
1> /opt/data/workdir/sca8_20.guppy.std.out \
2> /opt/data/workdir/sca8_20.guppy.std.err

# ----------
# 2024-05-16
# ----------

docker run \
--rm \
-it \
-v /Users/vladimirtomic/Documents/repos/strat/:/opt/repos/strat/ \
-v /Users/vladimirtomic/Documents/projects/ONT/data/:/opt/data/ \
strat:latest

# 50min
python /opt/repos/strat/scripts/strat_prepare.py \
--prefix 'TAAGATAATATATTTTTAAAAAATG' \
--suffix 'TAAAGCCAGGTTTTCTAACATGAAG' \
--motif 'CAG' \
--tolerance '{e<=5}' \
--cores 5 \
--input_path '/opt/data/sca8_11/fastq/guppy/' \
--output_path '/opt/data/sca8_11/output/guppy/'

cat /opt/data/sca8_11/output/guppy/*.ontarget.tsv 1> /opt/data/workdir/sca8_11.guppy.ontarget.tsv

python /opt/repos/strat/scripts/strat_process.py \
--motif 'CAG' \
--threshold '1' \
--input_path '/opt/data/workdir/sca8_11.guppy.ontarget.tsv' \
--output_path '/opt/data/workdir/' \
1> /opt/data/workdir/sca8_11.guppy.std.out \
2> /opt/data/workdir/sca8_11.guppy.std.err

# ----------
# 2024-05-13
# ----------

docker run \
--rm \
-it \
-v /Users/vladimirtomic/Documents/repos/strat/:/opt/repos/strat/ \
-v /Users/vladimirtomic/Documents/projects/ONT/data/:/opt/data/ \
strat:latest

# 2h 3min
python /opt/repos/strat/scripts/strat_prepare.py \
--prefix 'AGAAAGAAATGGTCTGTGATCCCCC' \
--suffix 'CATTCCCGGCTACAAGGACCCTTCG' \
--motif 'CAG' \
--tolerance '{e<=5}' \
--cores 1 \
--input_path '/opt/data/bc6_1/fastq/dorado/' \
--output_path '/opt/data/bc6_1/output/dorado/'

cat /opt/data/bc6_1/output/dorado/*.ontarget.tsv 1> /opt/data/workdir/bc6_7.dorado.ontarget.tsv

python /opt/repos/strat/scripts/strat_process.py \
--motif 'CAG' \
--threshold '1' \
--input_path '/opt/data/workdir/bc6_7.dorado.ontarget.tsv' \
--output_path '/opt/data/workdir/' \
1> /opt/data/workdir/bc6_7.dorado.std.out \
2> /opt/data/workdir/bc6_7.dorado.std.err

# 4h 34min
python /opt/repos/strat/scripts/strat_prepare.py \
--prefix 'AGAAAGAAATGGTCTGTGATCCCCC' \
--suffix 'CATTCCCGGCTACAAGGACCCTTCG' \
--motif 'CAG' \
--tolerance '{e<=5}' \
--cores 1 \
--input_path '/opt/data/bc6_2/fastq/dorado/' \
--output_path '/opt/data/bc6_2/output/dorado/'

cat /opt/data/bc6_2/output/dorado/*.ontarget.tsv 1> /opt/data/workdir/bc6_8.dorado.ontarget.tsv

python /opt/repos/strat/scripts/strat_process.py \
--motif 'CAG' \
--threshold '1' \
--input_path '/opt/data/workdir/bc6_8.dorado.ontarget.tsv' \
--output_path '/opt/data/workdir/' \
1> /opt/data/workdir/bc6_8.dorado.std.out \
2> /opt/data/workdir/bc6_8.dorado.std.err

# 4h 32min
python /opt/repos/strat/scripts/strat_prepare.py \
--prefix 'AGAAAGAAATGGTCTGTGATCCCCC' \
--suffix 'CATTCCCGGCTACAAGGACCCTTCG' \
--motif 'CAG' \
--tolerance '{e<=5}' \
--cores 1 \
--input_path '/opt/data/bc6_3/fastq/dorado/' \
--output_path '/opt/data/bc6_3/output/dorado/'

cat /opt/data/bc6_3/output/dorado/*.ontarget.tsv 1> /opt/data/workdir/bc6_9.dorado.ontarget.tsv

python /opt/repos/strat/scripts/strat_process.py \
--motif 'CAG' \
--threshold '1' \
--input_path '/opt/data/workdir/bc6_9.dorado.ontarget.tsv' \
--output_path '/opt/data/workdir/' \
1> /opt/data/workdir/bc6_9.dorado.std.out \
2> /opt/data/workdir/bc6_9.dorado.std.err

# 1h 54min
python /opt/repos/strat/scripts/strat_prepare.py \
--prefix 'AGAAAGAAATGGTCTGTGATCCCCC' \
--suffix 'CATTCCCGGCTACAAGGACCCTTCG' \
--motif 'CAG' \
--tolerance '{e<=5}' \
--cores 1 \
--input_path '/opt/data/bc6_4/fastq/dorado/' \
--output_path '/opt/data/bc6_4/output/dorado/'

cat /opt/data/bc6_4/output/dorado/*.ontarget.tsv 1> /opt/data/workdir/bc6_10.dorado.ontarget.tsv

python /opt/repos/strat/scripts/strat_process.py \
--motif 'CAG' \
--threshold '1' \
--input_path '/opt/data/workdir/bc6_10.dorado.ontarget.tsv' \
--output_path '/opt/data/workdir/' \
1> /opt/data/workdir/bc6_10.dorado.std.out \
2> /opt/data/workdir/bc6_10.dorado.std.err

# 1h 57min
python /opt/repos/strat/scripts/strat_prepare.py \
--prefix 'AGAAAGAAATGGTCTGTGATCCCCC' \
--suffix 'CATTCCCGGCTACAAGGACCCTTCG' \
--motif 'CAG' \
--tolerance '{e<=5}' \
--cores 1 \
--input_path '/opt/data/bc6_5/fastq/dorado/' \
--output_path '/opt/data/bc6_5/output/dorado/'

cat /opt/data/bc6_5/output/dorado/*.ontarget.tsv 1> /opt/data/workdir/bc6_11.dorado.ontarget.tsv

python /opt/repos/strat/scripts/strat_process.py \
--motif 'CAG' \
--threshold '1' \
--input_path '/opt/data/workdir/bc6_11.dorado.ontarget.tsv' \
--output_path '/opt/data/workdir/' \
1> /opt/data/workdir/bc6_11.dorado.std.out \
2> /opt/data/workdir/bc6_11.dorado.std.err

# 4h 36min
python /opt/repos/strat/scripts/strat_prepare.py \
--prefix 'AGAAAGAAATGGTCTGTGATCCCCC' \
--suffix 'CATTCCCGGCTACAAGGACCCTTCG' \
--motif 'CAG' \
--tolerance '{e<=5}' \
--cores 1 \
--input_path '/opt/data/bc6_6/fastq/dorado/' \
--output_path '/opt/data/bc6_6/output/dorado/'

cat /opt/data/bc6_6/output/dorado/*.ontarget.tsv 1> /opt/data/workdir/bc6_12.dorado.ontarget.tsv

python /opt/repos/strat/scripts/strat_process.py \
--motif 'CAG' \
--threshold '1' \
--input_path '/opt/data/workdir/bc6_12.dorado.ontarget.tsv' \
--output_path '/opt/data/workdir/' \
1> /opt/data/workdir/bc6_12.dorado.std.out \
2> /opt/data/workdir/bc6_12.dorado.std.err

# ----------
# 2024-04-15
# ----------

docker run \
--rm \
-it \
-v /Users/vladimirtomic/Documents/repos/strat/:/opt/repos/strat/ \
-v /Users/vladimirtomic/Documents/projects/ONT/data/:/opt/data/ \
strat:latest

# 1h 7min
python /opt/repos/strat/scripts/strat_prepare.py \
--prefix 'AGAAAGAAATGGTCTGTGATCCCCC' \
--suffix 'CATTCCCGGCTACAAGGACCCTTCG' \
--motif 'CAG' \
--tolerance '{e<=5}' \
--cores 2 \
--input_path '/opt/data/bc6_2/fastq/guppy/' \
--output_path '/opt/data/bc6_2/output/guppy/'

cat /opt/data/bc6_2/output/guppy/*.ontarget.tsv 1> /opt/data/workdir/bc6_8.guppy.ontarget.tsv

python /opt/repos/strat/scripts/strat_process.py \
--motif 'CAG' \
--threshold '1' \
--input_path '/opt/data/workdir/bc6_8.guppy.ontarget.tsv' \
--output_path '/opt/data/workdir/' \
1> /opt/data/workdir/bc6_8.guppy.std.out \
2> /opt/data/workdir/bc6_8.guppy.std.err

# 1h 5min
python /opt/repos/strat/scripts/strat_prepare.py \
--prefix 'AGAAAGAAATGGTCTGTGATCCCCC' \
--suffix 'CATTCCCGGCTACAAGGACCCTTCG' \
--motif 'CAG' \
--tolerance '{e<=5}' \
--cores 2 \
--input_path '/opt/data/bc6_3/fastq/guppy/' \
--output_path '/opt/data/bc6_3/output/guppy/'

cat /opt/data/bc6_3/output/guppy/*.ontarget.tsv 1> /opt/data/workdir/bc6_9.guppy.ontarget.tsv

python /opt/repos/strat/scripts/strat_process.py \
--motif 'CAG' \
--threshold '1' \
--input_path '/opt/data/workdir/bc6_9.guppy.ontarget.tsv' \
--output_path '/opt/data/workdir/' \
1> /opt/data/workdir/bc6_9.guppy.std.out \
2> /opt/data/workdir/bc6_9.guppy.std.err

# 1h 2min
python /opt/repos/strat/scripts/strat_prepare.py \
--prefix 'AGAAAGAAATGGTCTGTGATCCCCC' \
--suffix 'CATTCCCGGCTACAAGGACCCTTCG' \
--motif 'CAG' \
--tolerance '{e<=5}' \
--cores 2 \
--input_path '/opt/data/bc6_4/fastq/guppy/' \
--output_path '/opt/data/bc6_4/output/guppy/'

cat /opt/data/bc6_4/output/guppy/*.ontarget.tsv 1> /opt/data/workdir/bc6_10.guppy.ontarget.tsv

python /opt/repos/strat/scripts/strat_process.py \
--motif 'CAG' \
--threshold '1' \
--input_path '/opt/data/workdir/bc6_10.guppy.ontarget.tsv' \
--output_path '/opt/data/workdir/' \
1> /opt/data/workdir/bc6_10.guppy.std.out \
2> /opt/data/workdir/bc6_10.guppy.std.err

# 0h 48min
python /opt/repos/strat/scripts/strat_prepare.py \
--prefix 'AGAAAGAAATGGTCTGTGATCCCCC' \
--suffix 'CATTCCCGGCTACAAGGACCCTTCG' \
--motif 'CAG' \
--tolerance '{e<=5}' \
--cores 2 \
--input_path '/opt/data/bc6_5/fastq/guppy/' \
--output_path '/opt/data/bc6_5/output/guppy/'

cat /opt/data/bc6_5/output/guppy/*.ontarget.tsv 1> /opt/data/workdir/bc6_11.guppy.ontarget.tsv

python /opt/repos/strat/scripts/strat_process.py \
--motif 'CAG' \
--threshold '1' \
--input_path '/opt/data/workdir/bc6_11.guppy.ontarget.tsv' \
--output_path '/opt/data/workdir/' \
1> /opt/data/workdir/bc6_11.guppy.std.out \
2> /opt/data/workdir/bc6_11.guppy.std.err

# 0h 52min
python /opt/repos/strat/scripts/strat_prepare.py \
--prefix 'AGAAAGAAATGGTCTGTGATCCCCC' \
--suffix 'CATTCCCGGCTACAAGGACCCTTCG' \
--motif 'CAG' \
--tolerance '{e<=5}' \
--cores 2 \
--input_path '/opt/data/bc6_6/fastq/guppy/' \
--output_path '/opt/data/bc6_6/output/guppy/'

cat /opt/data/bc6_6/output/guppy/*.ontarget.tsv 1> /opt/data/workdir/bc6_12.guppy.ontarget.tsv

python /opt/repos/strat/scripts/strat_process.py \
--motif 'CAG' \
--threshold '1' \
--input_path '/opt/data/workdir/bc6_12.guppy.ontarget.tsv' \
--output_path '/opt/data/workdir/' \
1> /opt/data/workdir/bc6_12.guppy.std.out \
2> /opt/data/workdir/bc6_12.guppy.std.err

# ----------
# 2024-04-03
# ----------

docker run \
--rm \
-it \
-v /Users/vladimirtomic/Documents/repos/strat/:/opt/repos/strat/ \
-v /Users/vladimirtomic/Documents/projects/ONT/data/:/opt/data/ \
strat:latest

# 0h 43min
python /opt/repos/strat/scripts/strat_prepare.py \
--prefix 'AGAAAGAAATGGTCTGTGATCCCCC' \
--suffix 'CATTCCCGGCTACAAGGACCCTTCG' \
--motif 'CAG' \
--tolerance '{e<=5}' \
--cores 2 \
--input_path '/opt/data/bc6_1/fastq/guppy/' \
--output_path '/opt/data/bc6_1/output/guppy/'

cat /opt/data/bc6_1/output/guppy/*.ontarget.tsv 1> /opt/data/workdir/bc6_7.guppy.ontarget.tsv

python /opt/repos/strat/scripts/strat_process.py \
--motif 'CAG' \
--threshold '1' \
--input_path '/opt/data/workdir/bc6_7.guppy.ontarget.tsv' \
--output_path '/opt/data/workdir/' \
1> /opt/data/workdir/bc6_7.guppy.std.out \
2> /opt/data/workdir/bc6_7.guppy.std.err

# ----------
# 2024-03-18
# ----------

docker run \
--rm \
-it \
-v /Users/vladimirtomic/Documents/repos/strat/:/opt/repos/strat/ \
-v /Users/vladimirtomic/Documents/projects/ONT/data/:/opt/data/ \
strat:latest

# 1h 31min
python /opt/repos/strat/scripts/strat_prepare.py \
--prefix 'AGAAAGAAATGGTCTGTGATCCCCC' \
--suffix 'CATTCCCGGCTACAAGGACCCTTCG' \
--motif 'CAG' \
--tolerance '{e<=5}' \
--cores 2 \
--input_path '/opt/data/bc6_1/fastq/guppy/' \
--output_path '/opt/data/bc6_1/output/guppy/'

# 1h 32min
python /opt/repos/strat/scripts/strat_prepare.py \
--prefix 'AGAAAGAAATGGTCTGTGATCCCCC' \
--suffix 'CATTCCCGGCTACAAGGACCCTTCG' \
--motif 'CAG' \
--tolerance '{e<=5}' \
--cores 2 \
--input_path '/opt/data/bc6_2/fastq/guppy/' \
--output_path '/opt/data/bc6_2/output/guppy/'

# 1h 31min
python /opt/repos/strat/scripts/strat_prepare.py \
--prefix 'AGAAAGAAATGGTCTGTGATCCCCC' \
--suffix 'CATTCCCGGCTACAAGGACCCTTCG' \
--motif 'CAG' \
--tolerance '{e<=5}' \
--cores 2 \
--input_path '/opt/data/bc6_3/fastq/guppy/' \
--output_path '/opt/data/bc6_3/output/guppy/'

# 1h 27min
python /opt/repos/strat/scripts/strat_prepare.py \
--prefix 'AGAAAGAAATGGTCTGTGATCCCCC' \
--suffix 'CATTCCCGGCTACAAGGACCCTTCG' \
--motif 'CAG' \
--tolerance '{e<=5}' \
--cores 2 \
--input_path '/opt/data/bc6_4/fastq/guppy/' \
--output_path '/opt/data/bc6_4/output/guppy/'

# 0h 54min
python /opt/repos/strat/scripts/strat_prepare.py \
--prefix 'AGAAAGAAATGGTCTGTGATCCCCC' \
--suffix 'CATTCCCGGCTACAAGGACCCTTCG' \
--motif 'CAG' \
--tolerance '{e<=5}' \
--cores 2 \
--input_path '/opt/data/bc6_5/fastq/guppy/' \
--output_path '/opt/data/bc6_5/output/guppy/'

# 0h 57min
python /opt/repos/strat/scripts/strat_prepare.py \
--prefix 'AGAAAGAAATGGTCTGTGATCCCCC' \
--suffix 'CATTCCCGGCTACAAGGACCCTTCG' \
--motif 'CAG' \
--tolerance '{e<=5}' \
--cores 2 \
--input_path '/opt/data/bc6_6/fastq/guppy/' \
--output_path '/opt/data/bc6_6/output/guppy/'

cat /opt/data/bc6_1/output/guppy/*.ontarget.tsv 1> /opt/data/workdir/bc6_1.guppy.ontarget.tsv
cat /opt/data/bc6_2/output/guppy/*.ontarget.tsv 1> /opt/data/workdir/bc6_2.guppy.ontarget.tsv
cat /opt/data/bc6_3/output/guppy/*.ontarget.tsv 1> /opt/data/workdir/bc6_3.guppy.ontarget.tsv
cat /opt/data/bc6_4/output/guppy/*.ontarget.tsv 1> /opt/data/workdir/bc6_4.guppy.ontarget.tsv
cat /opt/data/bc6_5/output/guppy/*.ontarget.tsv 1> /opt/data/workdir/bc6_5.guppy.ontarget.tsv
cat /opt/data/bc6_6/output/guppy/*.ontarget.tsv 1> /opt/data/workdir/bc6_6.guppy.ontarget.tsv

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

# ----------
# 2024-02-25
# ----------

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

# ----------
# 2024-02-18
# ----------

# Build bwa-mem2 docker image
docker build -t bwa-mem2:latest .

docker run \
--rm \
-it \
-v /Users/vladimirtomic/Documents/projects/ONT/data/:/opt/data/ \
bwa-mem2:latest

/opt/bwa-mem2-2.2.1_x64-linux/bwa-mem2 mem -k 1 -P -L 100,100 -T 5 chrI_10.fa seq_10.fastq

# ----------
# 2024-02-14
# ----------

docker run \
--rm \
-it \
-v /Users/vladimirtomic/Documents/repos/strat/:/opt/repos/strat/ \
-v /Users/vladimirtomic/Documents/projects/ONT/data/:/opt/data/ \
strat:latest

# 1h 17min
python /opt/repos/strat/scripts/strat_prepare.py \
--prefix 'AGAAAGAAATGGTCTGTGATCCCCC' \
--suffix 'CATTCCCGGCTACAAGGACCCTTCG' \
--motif 'CAG' \
--tolerance '{e<=5}' \
--cores 2 \
--input_path '/opt/data/bc3_1/fastq/dorado/' \
--output_path '/opt/data/bc3_1/output/dorado/'

# 0h 4min
python /opt/repos/strat/scripts/strat_prepare.py \
--prefix 'AGAAAGAAATGGTCTGTGATCCCCC' \
--suffix 'CATTCCCGGCTACAAGGACCCTTCG' \
--motif 'CAG' \
--tolerance '{e<=5}' \
--cores 2 \
--input_path '/opt/data/bc3_2/fastq/dorado/' \
--output_path '/opt/data/bc3_2/output/dorado/'

# 0h 19min
python /opt/repos/strat/scripts/strat_prepare.py \
--prefix 'AGAAAGAAATGGTCTGTGATCCCCC' \
--suffix 'CATTCCCGGCTACAAGGACCCTTCG' \
--motif 'CAG' \
--tolerance '{e<=5}' \
--cores 2 \
--input_path '/opt/data/bc3_3/fastq/dorado/' \
--output_path '/opt/data/bc3_3/output/dorado/'

python /opt/repos/strat/scripts/strat_process.py \
--motif 'CAG' \
--threshold '10' \
--input_path '/opt/data/bc3_1/output/dorado/dorado.ontarget.tsv' \
--output_path '/opt/data/bc3_1/output/dorado/'

python /opt/repos/strat/scripts/strat_process.py \
--motif 'CAG' \
--threshold '10' \
--input_path '/opt/data/bc3_2/output/dorado/dorado.ontarget.tsv' \
--output_path '/opt/data/bc3_2/output/dorado/'

python /opt/repos/strat/scripts/strat_process.py \
--motif 'CAG' \
--threshold '10' \
--input_path '/opt/data/bc3_3/output/dorado/dorado.ontarget.tsv' \
--output_path '/opt/data/bc3_3/output/dorado/'

# ----------
# 2024-02-04
# ----------

docker run \
--rm \
-it \
-p 8888:8888 \
-v /Users/vladimirtomic/Documents/repos/strat/:/opt/repos/strat/ \
-v /Users/vladimirtomic/Documents/projects/ONT/data/:/opt/data/ \
strat:latest

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
--threshold '100' \
--input_path '/opt/data/bc3_2/output/guppy/guppy.ontarget.tsv' \
--output_path '/opt/data/bc3_2/output/guppy/'

python /opt/repos/strat/scripts/strat_process.py \
--motif 'CAG' \
--threshold '100' \
--input_path '/opt/data/bc3_3/output/guppy/guppy.ontarget.tsv' \
--output_path '/opt/data/bc3_3/output/guppy/'

# Run notebooks inside the docker container
jupyter notebook --ip 0.0.0.0 --port 8888 --no-browser --allow-root

docker run \
--rm \
-it \
-v /Users/vladimirtomic/Documents/repos/strat/:/opt/repos/strat/ \
-v /Users/vladimirtomic/Documents/projects/ONT/data/:/opt/data/ \
strat:latest

# 1h 46min
python /opt/repos/strat/scripts/strat_prepare.py \
--prefix 'AGAAAGAAATGGTCTGTGATCCCCC' \
--suffix 'CATTCCCGGCTACAAGGACCCTTCG' \
--motif 'CAG' \
--tolerance '{e<=5}' \
--cores 2 \
--input_path '/opt/data/bc3_1/fastq/guppy/' \
--output_path '/opt/data/bc3_1/output/guppy/'

# 3h 18min
python /opt/repos/strat/scripts/strat_prepare.py \
--prefix 'AGAAAGAAATGGTCTGTGATCCCCC' \
--suffix 'CATTCCCGGCTACAAGGACCCTTCG' \
--motif 'CAG' \
--tolerance '{e<=5}' \
--cores 2 \
--input_path '/opt/data/bc3_2/fastq/guppy/' \
--output_path '/opt/data/bc3_2/output/guppy/'

# 3h 15min
python /opt/repos/strat/scripts/strat_prepare.py \
--prefix 'AGAAAGAAATGGTCTGTGATCCCCC' \
--suffix 'CATTCCCGGCTACAAGGACCCTTCG' \
--motif 'CAG' \
--tolerance '{e<=5}' \
--cores 2 \
--input_path '/opt/data/bc3_3/fastq/guppy/' \
--output_path '/opt/data/bc3_3/output/guppy/'

# ----------
# 2024-01-28
# ----------

cat *.ontarget.tsv 1> guppy.ontarget.tsv
cat *.ontarget.tsv 1> dorado.ontarget.tsv

# hifi
pcr2persons guppy 5, 16, (36, 37), (80, 81)
pcr2persons dorado 5, 16, (36, 37), (79, 80)
jovan guppy 5, 11, 21, 33, 37, 68
jovan dorado 5, 11, 21, 33, 37, 68
dm108 guppy 5, 14, 20?, 86?, 295-420

# ontarget imp
pcr2persons guppy 5, 16, 36, 81
pcr2persons dorado 5, 16, 36, 81
jovan guppy 5, 11, 21, 33, 37, 68
jovan dorado 5, 11, 21, 33, 37, 68
dm108 guppy 5, 14, 19/20?, 85, 350, 380 240-470

# ----------
# 2024-01-27
# ----------

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

# 2h 15min
python /opt/repos/strat/scripts/strat_prepare.py \
--prefix 'AGAAAGAAATGGTCTGTGATCCCCC' \
--suffix 'CATTCCCGGCTACAAGGACCCTTCG' \
--motif 'CAG' \
--tolerance '{e<=5}' \
--cores 2 \
--input_path '/opt/data/pcr2persons/fastq/dorado/' \
--output_path '/opt/data/pcr2persons/output/dorado/'

# 6h 10min
python /opt/repos/strat/scripts/strat_prepare.py \
--prefix 'AGAAAGAAATGGTCTGTGATCCCCC' \
--suffix 'CATTCCCGGCTACAAGGACCCTTCG' \
--motif 'CAG' \
--tolerance '{e<=5}' \
--cores 2 \
--input_path '/opt/data/jovan/fastq/guppy/' \
--output_path '/opt/data/jovan/output/guppy/'

# 3h 15min
python /opt/repos/strat/scripts/strat_prepare.py \
--prefix 'AGAAAGAAATGGTCTGTGATCCCCC' \
--suffix 'CATTCCCGGCTACAAGGACCCTTCG' \
--motif 'CAG' \
--tolerance '{e<=5}' \
--cores 2 \
--input_path '/opt/data/jovan/fastq/dorado/' \
--output_path '/opt/data/jovan/output/dorado/'

# 2h 43min
python /opt/repos/strat/scripts/strat_prepare.py \
--prefix 'AGAAAGAAATGGTCTGTGATCCCCC' \
--suffix 'CATTCCCGGCTACAAGGACCCTTCG' \
--motif 'CAG' \
--tolerance '{e<=5}' \
--cores 2 \
--input_path '/opt/data/dm108/fastq/guppy/' \
--output_path '/opt/data/dm108/output/guppy/'

# ----------
# 2024-01-21
# ----------

# Build STRAT docker image
docker build -t strat:latest .

# Start interactive docker session with STRAT docker image
docker run \
--rm \
-it \
-p 8888:8888 \
-v /Users/vladimirtomic/Documents/repos/strat/:/opt/repos/strat/ \
-v /Users/vladimirtomic/Documents/projects/ONT/data/:/opt/data/ \
strat:latest

# Run notebooks inside the docker container
jupyter notebook --ip 0.0.0.0 --port 8888 --no-browser --allow-root

# ----------
# 2024-01-20
# ----------

# pyenv install 3.8.18
# pyenv virtualenv 3.8.18 3.8.18-strat
pip install --upgrade pip
pip install pyalign==0.4.4
pip install regex==2023.12.25
pip install jupyter==1.0.0
pip install pandas==2.0.3
pip install seaborn==0.13.1

# ----------
# 2024-01-14
# ----------

cat * 1> guppy.hifi.tsv

# ----------
# 2024-01-08
# ----------

docker build Dockerfile -t squigulator:0.2.2

docker run \
-v "/Users/vladimirtomic/Documents/projects/ONT/data/dmpk":"/data/dmpk" \
--rm \
squigulator:0.2.2 \
/data/dmpk/dmpk_fwd_050_0005.fa \
-x dna-r10-min \
-o /data/dmpk/dmpk_fwd_050_0005_reads.blow5 \
--sample-rate 5000 \
--full-contigs

# ----------
# 2023-12-11
# ----------

unzip '*.zip' -d Plazmid
