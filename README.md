# paper2

## development setup

The following is not necessary to run if you are following along with the
tutorial. Only run these commands if you are publishing a new version or
release of this document!

    conda env create -n protocols --file qiime2-$RELEASE-py36-$ARCH-conda.yml
    conda activate protocols
    conda install -c conda-forge songbird redbiom -y
    conda install -c bioconda bowtie2 -y
    pip install https://github.com/knights-lab/SHOGUN/archive/v1.0.7.zip
    pip install https://github.com/qiime2/q2-shogun/archive/master.zip
    conda install cytoolz -y
    conda install -c conda-forge sphinx awscli -y

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
