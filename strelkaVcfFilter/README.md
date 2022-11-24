# Utility Scripts

## strelkaVcfFilter

This script filters a VCF file from Strelka by PASS
and by the minimum value of the SomaticEVS parameter. 

## Systems Requirements

This script requires R 3.5.1 or higher, and the installation of the
following R packages:

- ```getopt```, for supporting the command line options
- ```VariantAnnotation```, for reading VCF files.

## Usage

For the help on how to use this file type:

```
./strelkaVcfFilter.R -h
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
cp strelkaVcfFilter.R ~/bin/strelkaVcfFilter
```

You should now have the scripts working anywhere in your command line.
