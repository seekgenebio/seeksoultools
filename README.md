# SeekSoulTools

SeekSoulTools is a professional bioinformatics toolkit developed by SeekGene for comprehensive single-cell transcriptome analysis. It provides robust solutions for RNA-seq, V(D)J sequencing, and full-length transcriptome analysis, featuring high accuracy, efficiency, and user-friendly interfaces. The toolkit is specifically optimized for SeekOne® series kits while maintaining compatibility with various custom designs.

## Features

The tools includes the following main modules:

- **rna module**: Identifies cell barcodes, performs alignment and quantification, generates cell expression matrices for downstream analysis, and supports cell clustering and differential analysis. Compatible with SeekOne® series kits and various custom designs.
- **fast module**: Specifically designed for SeekOne® DD single-cell full-sequence transcriptome kits and FFPE samples, enabling barcode extraction, paired-end read alignment, quantification, and full-sequence specific metrics.
- **vdj module**: Tailored for SeekOne® DD single-cell immune analysis kits, facilitating the assembly, filtering, and annotation of immune receptors.
- **multivdj module**: Integrated analysis pipeline for combined RNA-seq and V(D)J sequencing data, supporting simultaneous processing of RNA, TCR, and BCR data with comprehensive reporting.
- **utils module**: Contains additional utility tools to assist in data processing and analysis.

## System Requirements

- Linux operating system
- [Conda](https://docs.conda.io/en/latest/) package manager
- Minimum 8GB RAM (16GB recommended)
- Sufficient disk space for data processing (at least 100GB recommended)

## Installation

1. Clone the repository:
```bash
git clone https://github.com/seekgenebio/seeksoultools.git
cd seeksoultools
```

2. Create and activate conda environment:

For users in China:
```bash
conda env create -n seeksoultools -f conda_dependencies.zh.yml
conda activate seeksoultools
```

For international users:
```bash
conda env create -n seeksoultools -f conda_dependencies.yml
conda activate seeksoultools
```

Note: If you encounter slow download speeds or installation issues with pip packages, you can try using alternative PyPI mirrors in your region:
```bash
# Use PyPI mirrors
pip config set global.index-url https://pypi.org/simple
```

3. Install the  package:
```bash
pip install . src/simpleqc/target/wheels/simpleqc-0.1.0-py3-none-manylinux_2_17_x86_64.manylinux2014_x86_64.whl src/search-pattern/target/wheels/search_pattern-0.1.0-py3-none-manylinux_2_5_x86_64.manylinux1_x86_64.whl
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

