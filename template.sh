# ----------
# 2024-05-18
# ----------

docker run \
--rm \
-it \
-v /Users/vladimirtomic/Documents/repos/strat/:/opt/repos/strat/ \
-v /Users/vladimirtomic/Documents/projects/ONT/data/:/opt/data/ \
strat:latest

export DATASET=bc7_1_18
export DATASET=bc7_1_19
export DATASET=bc7_1_20
export DATASET=bc7_1_21
export DATASET=bc7_1_22
export DATASET=bc7_1_23
export DATASET=bc7_1_24
export DATASET=bc7_2_18
export DATASET=bc7_2_19
export DATASET=bc7_2_20
export DATASET=bc7_2_21
export DATASET=bc7_2_22
export DATASET=bc7_2_23
export DATASET=bc7_2_24
export CALLER=guppy

# Xmin
python /opt/repos/strat/scripts/strat_prepare.py \
--prefix 'AGAAAGAAATGGTCTGTGATCCCCC' \
--suffix 'CATTCCCGGCTACAAGGACCCTTCG' \
--motif 'CAG' \
--tolerance '{e<=5}' \
--cores 2 \
--input_path "/opt/data/$DATASET/fastq/$CALLER/" \
--output_path "/opt/data/$DATASET/output/$CALLER/" \
1> /opt/data/workdir/$DATASET.$CALLER.strat_prepare.std.out \
2> /opt/data/workdir/$DATASET.$CALLER.strat_prepare.std.err

cat /opt/data/$DATASET/output/$CALLER/*.ontarget.tsv 1> /opt/data/workdir/$DATASET.$CALLER.ontarget.tsv

python /opt/repos/strat/scripts/strat_process.py \
--motif 'CAG' \
--threshold '1' \
--input_path "/opt/data/workdir/$DATASET.$CALLER.ontarget.tsv" \
--output_path "/opt/data/workdir/" \
1> /opt/data/workdir/$DATASET.$CALLER.strat_process.std.out \
2> /opt/data/workdir/$DATASET.$CALLER.strat_process.std.err
