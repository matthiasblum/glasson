# glasson
Glasson (from the French *glaçon* for *ice cube*) is a binary, compressed file format for storing Hi-C contact maps.

*Seriously? Another format for Hi-C?*

Yes, unfortunately.

## Background

There are many formats for Hi-C data, but this one is mine. It is my best friend. It is my li&hellip;

There are many formats for Hi-C data. 
In fact, every PhD student, researcher, or engineer who analysed Hi-C data at one point has probably developed their own format. 
However, for sake of shortness, we are going to restrict the list of formats to major ones:

- Dense: *n &times; n* matrices stored in text files, sometimes gzip-compressed.
- [*hiclib*](https://bitbucket.org/mirnylab/hiclib) format: HDF5-based, used in [Imakaev et al., Nat Methods. 2012](http://www.nature.com/nmeth/journal/v9/n10/full/nmeth.2148.html).
- COO sparse matrices stored in text files, used in Hi-C Pro ([Servant et al., Genome Biol. 2015](http://genomebiology.biomedcentral.com/articles/10.1186/s13059-015-0831-x)).
- *hic* format, created for Juicebox ([Durand et al., Cell Syst. 2016](http://www.sciencedirect.com/science/article/pii/S240547121500054X)), described in the Juicer paper ([Durand et al., Cell Syst. 2016](http://www.sciencedirect.com/science/article/pii/S2405471216302198)).
- [*cool*](https://github.com/mirnylab/cooler/) format, HDF5-based, developed by the [Mirny lab](http://mirnylab.mit.edu/).

Text files are convenient because they are human readable, and can be easily loaded in any programming language, but they do not allow a fast access to a subset of a contact map.

HDF5 files can contains multiple contact maps, allow fast queries, and there are interfaces to the HDF5 format in almost every programming language.
However, remote access is not supported: one has to download the entire file to extract any information.

Finally, the `hic` format does allow remote access, and can embed contact maps at multiple resolutions. Nevertheless, it was not described yet when I started to work on Hi-C data.

The `glasson` format is a generic binary file format aiming to:

- store intra-chromosomal contact maps, at multiple resolutions if desired
- allow fast queries on intervals
- enable remote access
- use as little disk space as possible

## Installation

This implementation requires Python 2.7/3.4+ to create/read `glasson` files. No additional library is required. It is recommended to use a virtual environment.

    git clone https://github.com/matthiasblum/glasson.git
    cd glasson
    python setup.py install

### Store

TODO

### Extract

TODO

## Sweet! Should I use keep my Hi-C data in glasson files?

For Crick's sake, no! Although you are more than welcome to visualize your data on the QC Genomics genome browser, `glasson` is not meant to become the standard format for storing Hi-C data.

## Format specification
 
**Header**

| Field     | Description | Type |
|-----------|-------------|------|
|  signature  | 3 bytes string "GLA"  | char\[3\] |
|  idx_off  | file offset of the start of the index | unsigned long long |

**Body**

*For one row of one contact map of one chromosome:*

| Field | Description | Type |
|-------|-------------|------|
| n_items | number of non-zero items on the row | unsigned int |
| l_zstr | length of the zip-compressed string | unsigned int |
| zstr | zip-compressed string | char\[l_zstr\] |

`zstr` contains:

- an array `unsigned int[n_items]` representing to column indices in the matrix.
- an array `float[n_items]` representing the associated values.

**Index**

The index contains a list of chromosomes, and a list of contact maps for each chromosome in the following order, until the end of the file:

- chrom 1
    - map 1
    - map 2
    - &hellip;
- chrom 2
    - map 1
    - &hellip;
- &hellip;

*For one chromosome:*

| Field | Description | Type |
|-------|-------------|------|
| l_chrom   | length of chromosome | unsigned int     |
| l_chrom_name | length of chromosome name | unsigned int     |
| n_maps | number of contact maps | unsigned int |
| chrom_name | chromosome name | char\[l_chrom_name\]     |

*For one contact map:*

| Field | Description | Type |
|-------|-------------|------|
| res | resolution in bp | unsigned int     |
| n_rows | number of non-empty rows | unsigned int     |
| rows | array of non-empty row indices in the matrix | unsigned int\[n_rows\] |
| row_off | array of file offsets for non-empty rows | unsigned long long\[n_rows\] |

*Other chromosomes, and their contact maps, until the end of the file.*

## License

Public domain.

