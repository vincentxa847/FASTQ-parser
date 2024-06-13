# FASTQ parser

## The script opens a series of FASTQ files and outputs the corresponding FASTA files. This project demonstrates the following programming fundamentals:
- String comparisons
- Flowcontrol
- Simple and Nested Looping
- File I/O
- Complex objects besides lists
- Functions
- Exceptions
- Object-based Programming

## Running
The script can be run as follows, allowing user to apply a quality filter to remove low-quality sequences :
```
 python <program>.py file1.fastq file2.fastq ...--quality_filter=45
```
## Requirement
Codes is tested to work under
Python 3.11.1
