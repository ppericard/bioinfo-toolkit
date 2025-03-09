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

1. Clone this repository:
   ```bash
   git clone https://github.com/yourusername/bioinfo-toolkit.git
   cd bioinfo-toolkit
   ```

2. Install the dependencies:
   ```bash
   pip install -r requirements.txt
   ```

## Usage

You can run any script in the toolkit using the main entry point:

```bash
python bioinfo-toolkit.py <script_name> [arguments]
```

To see a list of all available scripts:

```bash
python bioinfo-toolkit.py
```

## Available Scripts

### FASTA Processing Tools

* **fasta_length_filter**: Filter FASTA sequences by length
* **fasta_n_filter**: Filter FASTA sequences by N content
* **fasta_name_filter**: Filter FASTA sequences by name/header
* **fasta_length_histo**: Generate length histogram of FASTA sequences
* **sort_fasta_by_length**: Sort FASTA sequences by length
* **gener_sample_fasta**: Generate sample FASTA datasets for testing

#### Examples

```bash
# Filter sequences by length
python bioinfo-toolkit.py fasta_length_filter -i input.fasta -o output.fasta -m 300 -M 1000

# Generate a histogram of sequence lengths
python bioinfo-toolkit.py fasta_length_histo -i input.fasta -o length_histogram.png
```

### FASTQ Processing Tools

* **get_pairs**: Process paired-end reads and separate them into pairs and singletons
* **get_pairs_v3**: Advanced implementation of get_pairs with improved memory usage
* **split_paired_fastq**: Split interleaved FASTQ files into separate files
* **fastq_umi_merge**: Merge FASTQ files with UMIs (Unique Molecular Identifiers)
* **fastq_name_filter**: Filter FASTQ sequences by name/header
* **gener_sample_fastq**: Generate sample FASTQ datasets for testing

#### Examples

```bash
# Process paired-end reads
python bioinfo-toolkit.py get_pairs -l reads_1.fastq -r reads_2.fastq -o output_dir

# Filter FASTQ reads by name
python bioinfo-toolkit.py fastq_name_filter -i input.fastq -o filtered.fastq -l names_list.txt
```

### File Format Conversion Tools

* **fastq_to_fasta**: Convert FASTQ files to FASTA format

#### Examples

```bash
# Convert FASTQ to FASTA
python bioinfo-toolkit.py fastq_to_fasta -i input.fastq -o output.fasta
```

### Utility Tools

* **script_template**: Template for creating new scripts for the toolkit

#### Examples

```bash
# Create a new script from template
python bioinfo-toolkit.py script_template -n my_new_script -c fasta
```

### Benchmarking and Testing Tools

* **test_get_pairs**: Test and compare different implementations of get_pairs
* **benchmark_get_pairs**: Benchmark the performance of get_pairs
* **create_test_datasets**: Create standardized test datasets for benchmarking

#### Examples

```bash
# Create test datasets for benchmarking
python bioinfo-toolkit.py create_test_datasets -c

# Benchmark get_pairs implementations
python bioinfo-toolkit.py test_get_pairs -v
```

## Utility Modules

The toolkit includes several utility modules that provide common functionality:

* **bioinfo_logger**: Standardized logging for bioinformatics scripts
* **fastq_utils**: Utility functions for FASTQ file handling
* **memory_tracker**: Track and analyze memory usage of Python scripts

## Testing

The repository includes a comprehensive test suite to ensure code quality and functionality. The tests are organized into categories:

- **Unit Tests**: Test individual functions and classes
- **Integration Tests**: Test interactions between components
- **Functional Tests**: Test complete workflows
- **Performance Tests**: Benchmark performance of different implementations

### Running Tests

You can run the tests using the included `run_tests.py` script:

```bash
# Run all tests
python run_tests.py

# Run specific test categories
python run_tests.py --unit
python run_tests.py --integration
python run_tests.py --functional
python run_tests.py --performance

# Generate code coverage report
python run_tests.py --coverage
```

Alternatively, you can use pytest directly:

```bash
# Install development dependencies
pip install -r requirements-dev.txt

# Run all tests
pytest

# Run specific test categories
pytest tests/unit
pytest tests/integration
pytest tests/functional
pytest tests/performance
```

## Repository Structure

```
bioinfo-toolkit/
├── bioinfo-toolkit.py      # Main entry point
├── run_tests.py            # Test runner script
├── src/                    # Source code
│   ├── scripts/            # Categorized scripts
│   │   ├── fasta/          # FASTA processing scripts
│   │   │   ├── fasta_length_filter.py
│   │   │   ├── fasta_n_filter.py
│   │   │   ├── fasta_name_filter.py
│   │   │   ├── fasta_length_histo.py
│   │   │   ├── sort_fasta_by_length.py
│   │   │   └── gener_sample_fasta.py
│   │   ├── fastq/          # FASTQ processing scripts
│   │   │   ├── get_pairs.py
│   │   │   ├── get_pairs_v3.py
│   │   │   ├── split_paired_fastq.py
│   │   │   ├── fastq_umi_merge.py
│   │   │   ├── fastq_name_filter.py
│   │   │   └── gener_sample_fastq.py
│   │   ├── conversion/     # Format conversion scripts
│   │   │   └── fastq_to_fasta.py
│   │   ├── tools/          # Utility tools
│   │   │   └── script_template.py
│   │   └── benchmark/      # Benchmarking scripts
│   │       ├── test_get_pairs.py
│   │       ├── benchmark_get_pairs.py
│   │       └── create_test_datasets.py
│   ├── archive/            # Archived deprecated scripts
│   │   └── tools/          # Archived tools
│   │       └── atomicblastplus.py
│   └── utils/              # Utility modules
│       ├── bioinfo_logger.py
│       ├── fastq_utils.py
│       └── memory_tracker.py
├── data/                   # Data directory
│   ├── test_data/          # Test datasets
│   └── benchmarks/         # Benchmark results
├── tests/                  # Test directory
│   ├── unit/               # Unit tests
│   ├── integration/        # Integration tests
│   ├── functional/         # Functional tests
│   └── performance/        # Performance tests
├── requirements.txt        # Dependencies
└── requirements-dev.txt    # Development dependencies
```

## Dependencies

### Runtime Dependencies
* psutil >= 5.9.0
* pandas >= 1.3.0
* matplotlib >= 3.5.0

### Development Dependencies
* pytest >= 7.0.0
* pytest-cov >= 4.1.0
* pytest-mock >= 3.10.0
* flake8 >= 6.0.0
* black >= 23.0.0
* memory_profiler >= 0.61.0

## License

This project is licensed under the [MIT License](LICENSE).

## Contributors

* Pierre Pericard - Main developer and maintainer
