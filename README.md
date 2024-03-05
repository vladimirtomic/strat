# STRAT

Short tandem repeats analysis tool (STRAT).

<style>
    .color-square {
        display: inline-block;
        width: 40px;
        height: 20px;
        margin-right: 5px;
        border: 1px solid black;
        text-align: center;
        line-height: 20px;
    }

    .green {
        background-color: #3DA853;
    }

    .blue {
        background-color: #4285F4;
    }

    .yellow {
        background-color: #F8BC07;
    }

    .red {
        background-color: #EA4334;
    }

    .white {
        background-color: white;
    }
</style>

<div class="color-square green">A</div>
<div class="color-square blue">C</div>
<div class="color-square yellow">G</div>
<div class="color-square red">T</div>
<div class="color-square white">DEL</div>

## 1. Setup

```
cd /repos
git clone https://github.com/vladimirtomic/strat.git
```

### 1.1 Build STRAT docker image

```
cd /repos/strat/dockerfiles/strat
docker build -t strat:latest .
```

### 1.2 Run STRAT docker container for notebooks

```
docker run \
--rm \
-it \
-p 8888:8888 \
-v /repos/strat/:/opt/repos/strat/ \
-v /data/:/opt/data/ \
strat:latest
```

```
jupyter notebook --ip 0.0.0.0 --port 8888 --no-browser --allow-root
```

### 1.3 Run STRAT docker container to execute scripts

```
docker run \
--rm \
-it \
-v /repos/strat/:/opt/repos/strat/ \
-v /data/:/opt/data/ \
strat:latest
```

#### 1.3.1 Execute STRAT prepare script

```
python /opt/repos/strat/scripts/strat_prepare.py \
--prefix 'AGAAAGAAATGGTCTGTGATCCCCC' \
--suffix 'CATTCCCGGCTACAAGGACCCTTCG' \
--motif 'CAG' \
--tolerance '{e<=5}' \
--cores 2 \
--input_path '/opt/data/bc3_1/fastq/guppy/' \
--output_path '/opt/data/bc3_1/output/guppy/'
```

#### 1.3.2 Merge results of STRAT prepare

```
cat /opt/data/bc3_1/output/guppy/*.ontarget.tsv 1> /opt/data/workdir/bc3_1.guppy.ontarget.tsv
```

#### 1.3.3 Execute STRAT process script

```
python /opt/repos/strat/scripts/strat_process.py \
--motif 'CAG' \
--threshold '1' \
--input_path '/opt/data/workdir/bc3_1.guppy.ontarget.tsv' \
--output_path '/opt/data/workdir/' \
1> /opt/data/workdir/bc3_1.guppy.std.out \
2> /opt/data/workdir/bc3_1.guppy.std.err
```
