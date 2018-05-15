## Embark fork of Germline

Embark (Adam G.) set out to fix Germline's -haploid and -homoz-only -w_extend issues in May 2018.

### Bugs in Columbia's released version (1.5.1)

* When running germline with `-haploid` and `-w_extend` flags and with `bits` (segment size) anything greater than 1, match tracts extend beyond the correct start and end coordinates. To demo this, run 'bash demo_haploid_w_extend.sh' after setting germline PATH to version 1.5.1 (this can be downloaded from Columbia's website).

* When running germline with `-homoz-only` and `-w_extend` flags and with `bits` (segment size) anything greater than 1, match tracts do not line up with the correct start and end coordinates.

* When running germline with `-haploid -w_extend` and `bits` anything greater than 1, the "snp_count" column is incorrect

* When running germline with `-haploid -w_extend bits=1`, with a tract ending at the end of the chromosome, germline does not report that the tract extends all the way to the end of the chromosome, instead reporting a tract one SNP short.

### Steps to update and release a new version

* Change the source code on a branch and test test test (`bash run_tests.sh in the /test directory`)
* Sync the code to an EC2 and run the tests there to confirm they work on Embark's ubuntu EC2s
* Document the updated tests
* Increment the version number in Germline_0001.cpp following semantic versioning
* Update the change log
* Sync the updated version to S3, and update the ancestry assignment EC2 provisioning script to use the new version
* Update the AMI for our EC2s using the updated provisioning script, so new EC2s come with the latest Germline version

### Change log

#### Version 1.5.2

* `germline -haploid -w_extend -bits 41` (or any bits value > 1) now produces correct match tracts
* exit codes follow shell-script convention: 1 for error, 0 for success
* `germline -version` now reports the version number
* added tests
