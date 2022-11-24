# Utility Scripts

## brassBedpeFilter

This script filters a Brass BEDPE file to remove
possible artefacts. Filters applied:

1. Remove if no assembly score is given
2. Remove if the distance between breakpoints is <1000bp
3. Remove if it is an inversion shorter than 5000bp and readpair count is less than 6

## Systems Requirements

This script requires R 3.5.1 or higher, and the installation of the
following R packages:

- ```getopt```, for supporting the command line options
- ```R.utils```, for the gunzip function.

## Usage

For the help on how to use this file type:

```
./brassBedpeFilter.R -h
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
cp brassBedpeFilter.R ~/bin/brassBedpeFilter
```

You should now have the scripts working anywhere in your command line.
