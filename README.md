# P6
Peter's Parse and Processing of PreGen Particulars via Pandas

- **adapted from https://github.com/VarenyaJ/P5/tree/exp/prioritize-conda**

## Install Conda
```bash
mkdir -p $HOME/miniconda3
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-<YOUR_SYSTEM>.sh -O $HOME/miniconda3/miniconda.sh
bash $HOME/miniconda3/miniconda.sh -b -u -p $HOME/miniconda3
#rm $HOME/miniconda3/miniconda.sh
source $HOME/.bashrc
source $HOME/miniconda3/bin/activate
conda init --all
conda --version
conda info
```


## Setup Conda (and optionally install Mamba)
```bash
# 1. Activate Conda for Shell
source $HOME/.bashrc && source $HOME/miniconda3/bin/activate && source $HOME/.bashrc && conda init --all && conda --version && conda info && conda list envs && which conda && conda --version

# 2. Setup Conda-Forge
conda update -n base -c defaults conda && conda install -n base -c conda-forge mamba conda-lock && conda list --show-channel-urls

# 3. Initialize conda
eval "$(conda shell.bash hook)" || echo 'no conda :('

# 4. OPTIONAL: Initialize mamba
eval "$(mamba shell hook --shell bash)" || echo 'no mamba :('
```

## Setup Project
```bash
cd $HOME
git clone https://github.com/VarenyaJ/P6.git
cd $HOME/P6/
git checkout main

# 1. Clear caches
conda clean --all -y
pip cache purge

# 2. Remove old env
conda deactivate
conda env remove -n p6 -y

# 3. Create new env
conda env create -f requirements/environment.yml -n p6 -y || mamba env create -f requirements/environment.yml -n p6 -y

# 4. TODO: Validate packages work with a little test
conda activate p6
python - <<EOF
import phenopackets, hpo-toolkit, pandas
print("OK:", phenopackets.__version__, hpo-toolkit.__version__, pandas.__version__)
EOF

# 4.5 TODO: Install and Verify Package
pip install -e .
python -c "import p6; print(p6.__version__)"
sample --help

pytest --maxfail=1 -q
```

# TODO:

### 5. Install lock tool & generate lock
```
conda install -n p6 -c conda-forge conda-lock -y || mamba install -n p6 -c conda-forge conda-lock -y
conda-lock lock -f requirements/environment.yml \
 -p linux-64 -p osx-64 -p osx-arm64
```

### Create lock
```bash
conda env create --yes -f requirements/environment.yml || mamba env create --yes -f requirements/environment.yml
conda-lock lock -f requirements/environment.yml -p linux-64 -p osx-64 -p win-64 --name p6
```

### Commit the generated lock files and update them via:
```bash
conda env update --yes -f requirements/environment.yml || mamba env update --yes -f requirements/environment.yml
conda-lock lock --update-lock-file
```
