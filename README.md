# Tool for Transforming SMAP to VCF Format

### Overview 
The SMAP to VCF converter tool is a stand alone python script that converts insertions and deletions in a SMAP file to VCF format.
The QUAL score is the SMAP confidence score for the given call. This output VCF file can be used for further downstream analysis using any tools that take VCF files as input.

###Usage

smap_to_vcf.py [-h] [-s SMAPPATH] [-n SAMPLE]

optional arguments:
*  -h, --help      show this help message and exit
*  -s SMAPPATH     Path to smap file to convert (required)
*  -n SAMPLE       Sample ID name for genotype data (optional, default "Sample1")

Note:  `python smap_to_vcf.py -h` to see usage on command line

### Requirements
This tool was designed to run with Python 2.7.  

### License
We offer this tool for open source use under the [MIT Software License](https://opensource.org/licenses/MIT). 
