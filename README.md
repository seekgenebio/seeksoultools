# seeksoultools
## Installation
```
git clone git@github.com:seekgenebio/seeksoultools.git
cd seeksoultools
conda env create -n seeksoultools -f conda_dependencies.yml
conda activate seeksoultools
pip install .
cd src/search-pattern
maturin build
pip install target/wheels/search_pattern-*.whl
cd ../simpleqc
maturin build
pip install target/wheels/search_pattern-*.whl
```
