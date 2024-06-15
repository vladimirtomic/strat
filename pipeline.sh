python /opt/repos/strat/scripts/strat_prepare.py \
--prefix "$PREFIX" \
--suffix "$SUFFIX" \
--motif "$MOTIF" \
--tolerance "$TOLERANCE" \
--cores "$CORES" \
--input_path "/opt/data/$DATASET/fastq/$CALLER/" \
--output_path "/opt/data/$DATASET/output/$CALLER/" \
1> /opt/data/workdir/$DATASET.$CALLER.strat_prepare.std.out \
2> /opt/data/workdir/$DATASET.$CALLER.strat_prepare.std.err

cat /opt/data/$DATASET/output/$CALLER/*.ontarget.tsv 1> /opt/data/workdir/$DATASET.$CALLER.ontarget.tsv

python /opt/repos/strat/scripts/strat_process.py \
--motif "$MOTIF" \
--threshold '1' \
--input_path "/opt/data/workdir/$DATASET.$CALLER.ontarget.tsv" \
--output_path "/opt/data/workdir/" \
1> /opt/data/workdir/$DATASET.$CALLER.strat_process.std.out \
2> /opt/data/workdir/$DATASET.$CALLER.strat_process.std.err
