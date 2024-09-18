export CALLER=guppy
export TOLERANCE="{e<=5}"
export CORES=2

export MOTIF_PRIM=CAG
export MOTIF_SCND=CAG
export PREFIX=AGAAAGAAATGGTCTGTGATCCCCC
export SUFFIX=CATTCCCGGCTACAAGGACCCTTCG

export DATASET=bc6_4_05
bash pipeline.sh

export DATASET=bc6_4_06
bash pipeline.sh

export DATASET=bc6_4_07
bash pipeline.sh

export DATASET=bc6_4_08
bash pipeline.sh

export DATASET=bc6_4_09
bash pipeline.sh

export DATASET=bc6_4_10
bash pipeline.sh
