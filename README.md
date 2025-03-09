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

## Main Scripts

* **get_pairs**: Separates paired reads and singletons from two paired FASTQ files (left and right)
* **fastq_to_fasta**: Converts FASTQ files to FASTA format
* **fasta_length_filter**: Filters FASTA sequences by length
* **atomicblastplus**: Submits a massively parallel Blast+ job-array to a computer cluster

## Usage Examples

### Processing Paired-End Reads

The `get_pairs.py` script separates paired and unpaired reads from two FASTQ files:

```bash
python bin/get_pairs.py -l reads_1.fastq -r reads_2.fastq -o output_dir
```

### Converting FASTQ to FASTA

Convert FASTQ files to FASTA format:

```bash
python bin/fastq_to_fasta.py -i input.fastq -o output.fasta
```

### Filtering FASTA by Length

Filter FASTA sequences by length:

```bash
python bin/fasta_length_filter.py -i input.fasta -o output.fasta -m 300 -M 1000
```

## Performance Benchmarking

The toolkit includes benchmarking tools that allow you to compare the performance of different implementations:

1. Create test datasets:
   ```bash
   python bin/test_get_pairs.py -c
   ```

2. Run the benchmarks:
   ```bash
   python bin/test_get_pairs.py
   ```

3. View detailed results including memory usage:
   ```bash
   python bin/test_get_pairs.py -v
   ```

This will generate benchmark results comparing execution time and memory usage across different dataset sizes.

## File Structure

* `bin/`: Contains the main scripts
* `test_data/`: Generated test datasets
* `benchmark_results/`: Benchmark output files

## Contributors

* Pierre Pericard (pierre.pericard@ed.univ-lille1.fr)
