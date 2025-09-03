#!/bin/bash

#SBATCH --job-name=wgd.viz
#SBATCH --cpus-per-task=1
#SBATCH --mem=10G
#SBATCH --time=0-01:00:00
#SBATCH --output=slurm.%x-%j.out

# input parameter
KS=${1:-wgd.ksd/Clusia_multiflora.paranome.tsv.ks.tsv}
SYNDIR=${2:-wgd.syn}
OUTDIR=${3:-wgd.viz}

# load envs and add haphic to PATH
module load conda
conda activate wgd
export NUMEXPR_MAX_THREADS=1

#### The visualization of KS age distribution and ELMM analysis (fit mixture model), mem=1G, time=1h
#/usr/bin/time wgd viz -d "$KS"

#### Substitution rate correction and integrate GMM anchor points (KS: orthologous ks dist, global_MRBH.tsv.ks.tsv, cpu=1, mem=10G, time=1h; Cmu only on login node: time=3h)
#wgd viz -d "$KS" -fa Clusia_multiflora -epk wgd.ksd/Clusia_multiflora.paranome.tsv.ks.tsv -sp speciestree.cmu.nw --reweight --plotapgmm -ap wgd_syn/iadhore-out/anchorpoints.txt -o wgd.viz.rate
#wgd viz -d wgd.ksd/global_MRBH.homoeologs.tsv.ks.tsv -fa Clusia_multiflora_H1 -epk wgd.ksd/Clusia_multiflora_H1.paranome.tsv.ks.tsv -sp speciestree.homoeologs.nw --reweight --xlim 0 3 --plotapgmm -ap wgd.syn.H1/iadhore-out/anchorpoints.txt -o wgd.viz.rate.H1
#wgd viz -d wgd.ksd/global_MRBH.homoeologs.tsv.ks.tsv -fa Clusia_multiflora_H2 -epk wgd.ksd/Clusia_multiflora_H2.paranome.tsv.ks.tsv -sp speciestree.homoeologs.nw --reweight --xlim 0 3 --plotapgmm -ap wgd.syn.H2/iadhore-out/anchorpoints.txt -o wgd.viz.rate.H2

#### Substitution rate correction and integrate GMM & ELMM models
#/usr/bin/time wgd viz -d wgd.ksd/global_MRBH.tsv.ks.tsv -fa Clusia_multiflora -epk wgd.ksd/Clusia_multiflora.paranome.tsv.ks.tsv -sp speciestree.cmu.nw --reweight --plotapgmm -ap wgd_syn/iadhore-out/anchorpoints.txt --plotelmm -o wgd.viz.mixed

#### Substitution rate correction (classic plots)
#/usr/bin/time wgd viz -d wgd.ksd/global_MRBH.cmu.tsv.ks.tsv -epk wgd.ksd/Clusia_multiflora.paranome.tsv.ks.tsv -sp speciestree.cmu.nw --reweight --toparrow --plotkde --classic --xlim 0 3 --adjustortho --okalpha 0.6 --plotapgmm -ap wgd.syn/iadhore-out/anchorpoints.txt --spair "Clusia_multiflora;Vitis_vinifera" --spair "Clusia_multiflora;Hypericum_perforatum" --spair "Clusia_multiflora;Garcinia_oblongifolia" --spair "Clusia_multiflora;Clusia_multiflora" -o wgd.viz.rate.classic
/usr/bin/time wgd viz -d wgd.ksd/global_MRBH.homoeologs.tsv.ks.tsv -epk wgd.ksd/Clusia_multiflora_H1.paranome.tsv.ks.tsv -sp speciestree.homoeologs.nw --reweight --toparrow --plotkde --classic --xlim 0 2 --adjustortho --okalpha 0.6 --plotapgmm -ap wgd.syn.H1/iadhore-out/anchorpoints.txt --spair "Clusia_multiflora_H1;Vitis_vinifera" --spair "Clusia_multiflora_H1;Hypericum_perforatum" --spair "Clusia_multiflora_H1;Garcinia_oblongifolia" --spair "Clusia_multiflora_H1;Clusia_multiflora_H2" --spair "Clusia_multiflora_H1;Clusia_multiflora_H1" -o wgd.viz.rate.classic.H1
