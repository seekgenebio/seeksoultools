# seeksoultools

SeekSoulTools is a comprehensive software suite developed by SeekGene for processing single-cell transcriptome data, supporting various single-cell sequencing data analysis workflows.

## Features

The tools includes the following main modules:

- **rna module**: Identifies cell barcodes, performs alignment and quantification, generates cell expression matrices for downstream analysis, and supports cell clustering and differential analysis. Compatible with SeekOne® series kits and various custom designs.
- **fast module**: Specifically designed for SeekOne® DD single-cell full-sequence transcriptome kits and FFPE samples, enabling barcode extraction, paired-end read alignment, quantification, and full-sequence specific metrics.
- **vdj module**: Tailored for SeekOne® DD single-cell immune analysis kits, facilitating the assembly, filtering, and annotation of immune receptors.
- **multivdj module**: Integrated analysis pipeline for combined RNA-seq and V(D)J sequencing data, supporting simultaneous processing of RNA, TCR, and BCR data with comprehensive reporting.
- **utils module**: Contains additional utility tools to assist in data processing and analysis.

## System Requirements

- Linux operating system
- Conda package manager
- Minimum 8GB RAM (16GB recommended)
- Sufficient disk space for data processing (at least 100GB recommended)

## Installation

1. Clone the repository:
```bash
git clone git@github.com:seekgenebio/seeksoultools.git
cd seeksoultools
```

2. Create and activate conda environment:
```bash
conda env create -n seeksoultools -f conda_dependencies.yml
conda activate seeksoultools
```

3. Install the main package:
```bash
pip install .
```

4. Install additional dependencies:
```bash
cd src/search-pattern
maturin build
pip install target/wheels/search_pattern-*.whl
cd ../simpleqc
maturin build
pip install target/wheels/search_pattern-*.whl
```

## Usage

Each module can be accessed through the command-line interface. For example:

```bash
# RNA-seq analysis
seeksoultools rna --help

# Full-length-seq analysis
seeksoultools fast --help

# V(D)J analysis
seeksoultools vdj --help

# Combined RNA and V(D)J analysis
seeksoultools multivdj --help
```

For detailed usage instructions and examples, please refer to the official documentation.


## Support

For technical support or questions, please contact SeekGene support team or open an issue on GitHub.
