# python /opt/repos/strat/scripts/strat_prepare.py \
# --prefix "$PREFIX" \
# --suffix "$SUFFIX" \
# --motif "$MOTIF_PRIM" \
# --tolerance "$TOLERANCE" \
# --cores "$CORES" \
# --input_path "/opt/data/$DATASET/fastq/$CALLER/" \
# --output_path "/opt/data/$DATASET/output/$CALLER/" \
# 1> /opt/data_out/workdir/$DATASET.$CALLER.strat_prepare.std.out \
# 2> /opt/data_out/workdir/$DATASET.$CALLER.strat_prepare.std.err

# cat /opt/data/$DATASET/output/$CALLER/*.ontarget.tsv 1> /opt/data_out/workdir/$DATASET.$CALLER.ontarget.tsv

# python /opt/repos/strat/scripts/strat_process.py \
# --motif_prim "$MOTIF_PRIM" \
# --motif_scnd "$MOTIF_SCND" \
# --threshold '1' \
# --input_path "/opt/data_out/workdir/$DATASET.$CALLER.ontarget.tsv" \
# --output_path "/opt/data_out/workdir/" \
# 1> /opt/data_out/workdir/$DATASET.$CALLER.strat_process.std.out \
# 2> /opt/data_out/workdir/$DATASET.$CALLER.strat_process.std.err

python /opt/repos/strat/scripts/strat_process_kmers.py \
--motif_prim "$MOTIF_PRIM" \
--motif_scnd "$MOTIF_SCND" \
--threshold '1' \
--input_path "/opt/data_out/kmers/$DATASET.$CALLER.kmers.tsv" \
--output_path "/opt/data_out/kmers_processed/" \
1> /opt/data_out/kmers_processed/$DATASET.$CALLER.strat_process_kmers.std.out \
2> /opt/data_out/kmers_processed/$DATASET.$CALLER.strat_process_kmers.std.err
