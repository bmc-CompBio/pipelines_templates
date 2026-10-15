# conda env create --file sra_kallisto_paired.yaml

conda activate sra_kallisto_paired
# conda install cmake

# git clone https://github.com/pachterlab/kallisto
# cd kallisto
# mkdir -p build
# cd build
# cmake .. -DCMAKE_INSTALL_PREFIX=~/miniconda3/envs/sra_kallisto_dev/ -DUSE_HDF5=ON
# make
# make install

snakemake -np
snakemake --dag | dot -Tsvg > dag.svg

## on local machine
snakemake --cores 16

## on the cluster
snakemake --profile slurm/ 
