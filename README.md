GERMLINE
--------
http://www1.cs.columbia.edu/~gusev/germline/

# About
GERMLINE is a program for discovering long shared segments of Identity by Descent (IBD) between pairs of individuals in a large population. It takes as input genotype or haplotype marker data for individuals (as well as an optional known pedigree) and generates a list of all pairwise segmental sharing.

GERMLINE uses a novel hashing & extension algorithm which allows for segment identification in haplotype data in time proportional to the number of individuals. With genotype data, GERMLINE implements several pedigree-based phasing techniques to impute data for related individuals, and then iteratively uses identified IBD segments to infer additional missing information. GERMLINE can identify shared segments of any specified length, as well as allow for any number of mismatching markers.

The program has been developed in Itsik Pe'er's Lab of Computational Genetics at Columbia University. It is built in C++ and tested in the Red Hat Linux environment; the source is distributed here in a tar.gz package under the GPL license. 

# Usage
Compile germline with make.
The executable is run as `germline <options>` which prompts the user for input/output file information and runs the algorithm.

# Input
GERMLINE accepts as input the following formats:

    * Plink / ped+map
    * PHASE / HapMap

GERMLINE also accepts an optional genetic map, formatted according to the HapMap standard described here ( http://ftp.hapmap.org/recombination/ ).

# Output
Upon completion, GERMLINE generates a .match file in the specified location. The first five lines contain meta-data detailing the run settings and executions time. The following rows detail the identified pairwise shared segments, one per row, with each row containing the following fields:

* Family ID 1
* Individual ID 1
* Family ID 2
* Individual ID 2
* Chromosome
* Segment start (bp)
* Segment end (bp)
* Segment start (SNP)
* Segment end (SNP)
* Total SNPs in segment
* Genetic length of segment
* Units for genetic length (cM or MB)
* Mismatching SNPs in segment
* 1 if Individual 1 is homozygous in match; 0 otherwise
* 1 if Individual 2 is homozygous in match; 0 otherwise


## Embark fork of Germline

Embark (Adam G.) set out to fix Germline's -haploid and -homoz-only -w_extend issues in May 2018.

### Bugs in Columbia's released version (1.5.1)

* When running germline with `-haploid` and `-w_extend` flags and with `bits` (segment size) anything greater than 1, match tracts extend beyond the correct start and end coordinates. To demo this, run the tests  after setting germline PATH to version 1.5.1 (this can be downloaded from Columbia's website).

* (STILL TO FIX) When running germline with `-homoz-only` and `-w_extend` flags and with `bits` (segment size) anything greater than 1, match tracts do not line up with the correct start and end coordinates.

* (STILL TO FIX) When running germline with `-haploid -w_extend` and `bits` anything greater than 1, the "snp_count" column is incorrect

* (STILL TO FIX) When running germline with `-haploid -w_extend bits=1`, with a tract ending at the end of the chromosome, germline does not report that the tract extends all the way to the end of the chromosome, instead reporting a tract one SNP short. To demo this,  add 1 to the BITS_TO_TEST variable in test/testing_constants_and_functions.sh, then run the tests.

### Feature wishlist

* Modify Germline to take a list of family_ids_we_want_outputs_for, so we can skip grepping new-individual match tracts out of output files
* Modify Germline so output files use consistent delimiters (e.g., all tabs rather than a mix of tabs and spaces)

### Steps to update and release a new version

* Change the source code on a branch and test test test (`cd test && bash run_tests.sh`)
* Sync the code to an EC2 and run the tests there to confirm they work on Embark's ubuntu EC2s
* Document the updated tests
* Increment the VERSION constant in [GERMLINE_0001.cpp](GERMLINE_0001.cpp) following semantic versioning
* Update the change log
* Run the `Build GERMLINE` GitHub Action with `make_release: true` to build and create a release.

### Building GERMLINE

To build GERMLINE, run `make` in the germline directory.
This will create the `bin/germline` and `bin/parse_bmatch` executables.

GERMLINE can also be built and tested using Docker / docker compose. To do this, run `docker compose up` from the repository root.

### Change log

#### v1.6.0-embark

* Support for using only subsets of samples using `-old_samples sample_list` and `-new_samples sample_list`: Samples for the old list will be compared to the new list
* Only loads the used samples in-memory. In many cases will use orders of magnitude less memory allowing for parallel execution on smaller-memory machines
* Comparison time reduced to whatever is needed from both lists. In many cases a order of magnitude faster

#### v1.5.2

* replicates stream-reading fix that David Riccardi made in 2016 (but didn't make into version control)
* `germline -haploid -w_extend -bits 41` (or any bits value > 1) now produces correct match tracts
* exit codes follow shell-script convention: 1 for error, 0 for success
* `germline -version` now reports the version number
* added tests

