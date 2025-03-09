# Bioinfo-Toolkit

A collection of Python scripts for bioinformatics analysis, focusing on processing and manipulating common file formats such as FASTA and FASTQ.

## Features

* Process paired-end FASTQ reads
* Filter FASTA sequences by length
* Convert between FASTQ and FASTA formats
* Submit parallel BLAST+ jobs
* Generate sample datasets
* And more

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

### Examples

#### Processing Paired-End Reads

The `get_pairs` script separates paired and unpaired reads from two FASTQ files:

```bash
python bioinfo-toolkit.py get_pairs -l reads_1.fastq -r reads_2.fastq -o output_dir
```

#### Converting FASTQ to FASTA

Convert FASTQ files to FASTA format:

```bash
python bioinfo-toolkit.py fastq_to_fasta -i input.fastq -o output.fasta
```

#### Filtering FASTA by Length

Filter FASTA sequences by length:

```bash
python bioinfo-toolkit.py fasta_length_filter -i input.fasta -o output.fasta -m 300 -M 1000
```

## Testing

The repository includes a comprehensive test suite to ensure code quality and functionality. The tests are organized into the following categories:

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

### Continuous Integration

This repository uses GitHub Actions for continuous integration. Tests are automatically run on both Linux and Windows environments with different Python versions to ensure cross-platform compatibility.

## Performance Benchmarking

The toolkit includes benchmarking tools that allow you to compare the performance of different implementations:

1. Create test datasets:
   ```bash
   python bioinfo-toolkit.py create_test_datasets -c
   ```

2. Run the benchmarks:
   ```bash
   python bioinfo-toolkit.py test_get_pairs
   ```

3. View detailed results including memory usage:
   ```bash
   python bioinfo-toolkit.py test_get_pairs -v
   ```

This will generate benchmark results comparing execution time and memory usage across different dataset sizes.

## Repository Structure

```
bioinfo-toolkit/
├── bioinfo-toolkit.py     # Main entry point
├── run_tests.py           # Test runner script
├── src/                   # Source code
│   ├── scripts/           # Categorized scripts
│   │   ├── fasta/         # FASTA processing scripts
│   │   ├── fastq/         # FASTQ processing scripts
│   │   ├── conversion/    # Format conversion scripts
│   │   ├── tools/         # Utility tools
│   │   └── benchmark/     # Benchmarking scripts
│   └── utils/             # Utility modules
├── data/                  # Data directory
│   ├── test_data/         # Test datasets
│   └── benchmarks/        # Benchmark results
├── tests/                 # Test directory
│   ├── unit/              # Unit tests
│   ├── integration/       # Integration tests
│   ├── functional/        # Functional tests
│   └── performance/       # Performance tests
├── requirements.txt       # Dependencies
└── requirements-dev.txt   # Development dependencies
```

## Contributors

* Pierre Pericard (pierre.pericard@ed.univ-lille1.fr)
