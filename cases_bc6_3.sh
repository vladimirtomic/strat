export CALLER=guppy
export TOLERANCE="{e<=5}"
export CORES=2

export MOTIF_PRIM=CAG
export MOTIF_SCND=CAG
export PREFIX=AGAAAGAAATGGTCTGTGATCCCCC
export SUFFIX=CATTCCCGGCTACAAGGACCCTTCG

export DATASET=bc6_3_05
bash pipeline.sh

export DATASET=bc6_3_06
bash pipeline.sh

export DATASET=bc6_3_07
bash pipeline.sh

export DATASET=bc6_3_08
bash pipeline.sh

export DATASET=bc6_3_09
bash pipeline.sh

export DATASET=bc6_3_10
bash pipeline.sh
