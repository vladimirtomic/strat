# python /opt/repos/strat/scripts/strat_prepare.py \
# --prefix "$PREFIX" \
# --suffix "$SUFFIX" \
# --motif "$MOTIF" \
# --tolerance "$TOLERANCE" \
# --cores "$CORES" \
# --input_path "/opt/data/$DATASET/fastq/$CALLER/" \
# --output_path "/opt/data/$DATASET/output/$CALLER/" \
# 1> /opt/data/workdir/$DATASET.$CALLER.strat_prepare.std.out \
# 2> /opt/data/workdir/$DATASET.$CALLER.strat_prepare.std.err

# cat /opt/data/$DATASET/output/$CALLER/*.ontarget.tsv 1> /opt/data/workdir/$DATASET.$CALLER.ontarget.tsv

python /opt/repos/strat/scripts/strat_process.py \
--motif_prim "$MOTIF_PRIM" \
--motif_scnd "$MOTIF_SCND" \
--threshold '1' \
--input_path "/opt/data_out/L25E5M0_workdir/$DATASET.$CALLER.ontarget.tsv" \
--output_path "/opt/data_out/L25E5M0_workdir_new/" \
1> /opt/data_out/L25E5M0_workdir_new/$DATASET.$CALLER.strat_process.std.out \
2> /opt/data_out/L25E5M0_workdir_new/$DATASET.$CALLER.strat_process.std.err
