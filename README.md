# paper2

## development setup

    conda env create -n protocols --file qiime2-$RELEASE-py36-$ARCH-conda.yml
    conda activate protocols
    conda install -c conda-forge songbird redbiom -y
    conda install -c bioconda bowtie2 -y
    pip install https://github.com/knights-lab/SHOGUN/archive/master.zip
    pip install https://github.com/qiime2/q2-shogun/archive/master.zip
    conda install -c conda-forge sphinx awscli -y
    conda install cytoolz

## building

    export AWS_ACCESS_KEY_ID='MY_KEY_ID'
    export AWS_SECRET_ACCESS_KEY='MY_SECRET_ACCESS_KEY'

    make html

    # once you have spot-checked the resulting HTML deploy the build:
    make deploy

    # and commit/push the latest `docs/` dir changes to this repo
    git add docs/*
    git commit -m "Some clever message here"
    git push origin master
