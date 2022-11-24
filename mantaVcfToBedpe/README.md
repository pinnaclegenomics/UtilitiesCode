# Utility Scripts

## mantaVcfToBedpe

This script converts a VCF file from Manta
into a BEDPE format (i.e. two breakpoints of the same rearrangement on
the same line). Moreover, a column indicating the strucutural variant
type is added (```svclass```), and the mutations can be filtered by PASS
and by the minimum value of the PR parameter. A flag for germline or
somatic variants must be specified. 

## Systems Requirements

This script requires R 3.5.1 or higher, and the installation of the
following R packages:

- ```getopt```, for supporting the command line options
- ```StructuralVariantAnnotation```, for the actual conversion from VCF
to BEDPE.
- ```VariantAnnotation```, for reading VCF files.

A version of the ```StructuralVariantAnnotation``` that works in R 3.5.1
can be installed as follows:

```
devtools::install_github("PapenfussLab/StructuralVariantAnnotation@cb584be9474318e7f91c2cd8fca99f1212cab021")
```
## Usage

For the help on how to use this file type:

```
./mantaVcfToBedpe.R -h
```

#### Install in your PATH

You can move this script to a directory in your PATH, so that it can be used from
the command line like any other program.

For example, you can create a ```bin``` directory in your home directory:

```
mkdir ~/bin
```

And add it to your path:

```
export PATH=~/bin:$PATH
```

You can add the above line to your ```~/.bashrc``` file so you don't have to retype it every time.
Then, you can simply copy the files and even remove the ```.R``` extension:

```
cp mantaVcfToBedpe.R ~/bin/mantaVcfToBedpe
```

You should now have the scripts working anywhere in your command line.
