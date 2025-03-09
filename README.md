# Bioinfo-Toolkit

A comprehensive collection of Python scripts for bioinformatics analysis, focusing on processing and manipulating common file formats such as FASTA and FASTQ.

## Features

* Process and manipulate FASTA files
* Handle paired-end FASTQ reads
* Convert between bioinformatics file formats
* Generate sample datasets for testing
* Benchmark script performance
* Memory tracking utilities

## Installation

### From Source

1. Clone this repository:
   ```bash
   git clone https://github.com/yourusername/bioinfo-toolkit.git
   cd bioinfo-toolkit
   ```

2. Install the package in development mode:
   ```bash
   pip install -e .
   ```

## Usage

You can run any script in the toolkit using the main entry point:

```bash
bioinfo-toolkit <script_name> [arguments]
```

To see a list of all available scripts:

```bash
bioinfo-toolkit
```

## Available Scripts

### FASTA Processing Tools

* **fasta_length_filter**: Filter FASTA sequences by length
* **fasta_n_filter**: Filter FASTA sequences by N content
* **fasta_name_filter**: Filter FASTA sequences by name/header
* **fasta_length_histo**: Generate length histogram of FASTA sequences
* **sort_fasta_by_length**: Sort FASTA sequences by length
* **gener_sample_fasta**: Generate sample FASTA datasets for testing

### FASTQ Processing Tools

* **get_pairs**: Process paired-end reads and separate them into pairs and singletons
* **get_pairs_v3**: Advanced implementation of get_pairs with improved memory usage
* **split_paired_fastq**: Split interleaved FASTQ files into separate files
* **fastq_umi_merge**: Merge FASTQ files with UMIs (Unique Molecular Identifiers)
* **fastq_name_filter**: Filter FASTQ sequences by name/header
* **gener_sample_fastq**: Generate sample FASTQ datasets for testing

### File Format Conversion Tools

* **fastq_to_fasta**: Convert FASTQ files to FASTA format

### Utility Tools

* **script_template**: Template for creating new scripts for the toolkit

### Benchmarking and Testing Tools

* **test_get_pairs**: Test and compare different implementations of get_pairs
* **benchmark_get_pairs**: Benchmark the performance of get_pairs
* **create_test_datasets**: Create standardized test datasets for benchmarking

## Development

### Setting Up Development Environment

1. Install development dependencies:
   ```bash
   pip install -e ".[dev]"
   ```

2. Run tests:
   ```bash
   pytest
   ```

3. Check code style:
   ```bash
   flake8 src tests
   ```

4. Format code:
   ```bash
   black src tests
   ```

## Repository Structure

```
bioinfo-toolkit/
├── bioinfo-toolkit        # Command-line entry point
├── pyproject.toml         # Project metadata and build configuration
├── src/                   # Source code
│   └── bioinfotoolkit/    # Main package
│       ├── cli.py         # Command-line interface
│       ├── scripts/       # Categorized scripts
│       │   ├── fasta/     # FASTA processing scripts
│       │   ├── fastq/     # FASTQ processing scripts
│       │   ├── conversion/# Format conversion scripts
│       │   ├── tools/     # Utility tools
│       │   └── benchmark/ # Benchmarking scripts
│       └── utils/         # Utility modules
├── tests/                 # Test directory
│   ├── unit/              # Unit tests
│   ├── integration/       # Integration tests
│   ├── functional/        # Functional tests
│   └── performance/       # Performance tests
└── data/                  # Data directory
    └── test_data/         # Test datasets
```

## License

This project is licensed under the [MIT License](LICENSE).

## Contributors

* Pierre Pericard - Main developer and maintainer
