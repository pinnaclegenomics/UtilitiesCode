# Utility Scripts

## Table of content

- [Introduction to the repository](#intro)
- [Version](#version)
- [Systems Requirements](#req)
- [Scripts provided by this repository](#scripts)



<a name="intro"/>

## Introduction to the repository

This repository contains a collection of scripts to convert/filter
mutation files. More specifically, these scripts are meant to be
a file format/filter bridge between the output of popular mutation
callers (e.g. Strelka, Manta, Caveman, Pindel and Ascat) and the
mutational signatures and annotation packages ```siganture.tools.lib```
and ```annotate.tools.lib```.

<a name="version"/>

## Version

1.0.1

- Added a function for converting Canvas vcf to copy number text files,
available via the ```prepareData``` script.

1.0

- Implementation of functions for filtering mutations from Caveman, Pindel,
Brass, Ascat, Strelka and Manta.
- The ```prepareData``` script can filter multiple mutations types at the
same time and organise them according to sample name.

<a name="req"/>

## Systems Requirements

Each individual script has its own requirements. In general, most
of the scripts require R 3.5.1 or higher, and the installation of
R packages. The exact requirements are given in the README.md files
in each of the script folders.

<a name="req"/>

## Scripts

Scripts available in this repository:

- **```mantaVcfToBedpe```**: this script converts a VCF file from Manta
into a BEDPE format (i.e. two breakpoints of the same rearrangement on
the same line). Moreover, a column indicating the strucutural variant
type is added (```svclass```), and the mutations can be filtered by PASS
and by the minimum value of the PR parameter. A flag for germline or
somatic variants must be specified. 
- **```strelkaVcfFilter```**: This script filters a VCF file from Strelka by PASS
and by the minimum value of the SomaticEVS parameter. 
- **```cavemanVcfFilter```**: This script filters a VCF file from Caveman by PASS
and by the ASMD and CLPM parameters. 
- **```pindelVcfFilter```**: This script filters a VCF file from Pindel by PASS
and by the QUAL and REP parameters. 
- **```brassBedpeFilter```**: This script filters a Brass BEDPE file to remove
possible artefacts. 
- **```ascatToTsv```**: This script converts an ascat comma separated values (CSV) file
into a tab separated values (TSV) file and adds headers.

## Installation

The ```utility.scripts``` R package needs to be installed using ```devtools:install()```
from the main package directory. After
the package is installed, the command line scripts can be used.

The script ```INSTALL.sh``` in the repository main directory can be used to quickly install
the scripts so that they are available in your command line. Note that ```INSTALL.sh``` will not
install the ```utility.scripts``` R package for you.

The ```INSTALL.sh```
will check whether the directory ```~/bin``` exists and is present in the ```PATH``` variable.
If not, ```~/bin``` will be created and added to ```PATH``` in either the ```.bashrc``` or the 
```.profile``` file. The above scripts will then be copied into ```~/bin```, and should be
available at command line at the next login or after typing ```PATH=~/bin:$PATH```.

