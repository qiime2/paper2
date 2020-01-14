# paper2

## development setup

    conda env create -n protocols --file qiime2-$RELEASE-py36-$ARCH-conda.yml
    conda activate protocols
    conda install -c conda-forge songbird redbiom -y
    conda install -c bioconda bowtie2 -y
    pip install https://github.com/knights-lab/SHOGUN/archive/master.zip
    pip install https://github.com/qiime2/q2-shogun/archive/master.zip
    conda install -c conda-forge sphinx -y
