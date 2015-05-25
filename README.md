# CPPPOA

CPPPOA is tool for generating consensus from multiple sequences. It is C++ implementation of Partial Order Alignment algorithm for Multiple sequence alignment based on simple python implementation (https://github.com/ljdursi/poapy) and explanation of this implementation given in blog post (http://simpsonlab.github.io/2015/05/01/understanding-poa/).

Algorithm is described in following papers:

[Multiple sequence alignment using partial order graphs](http://bioinformatics.oxfordjournals.org/content/18/3/452.short) (2002) by Lee, Grasso, and Sharlow
 and
[Generating consensus sequences from partial order multiple sequence alignment graphs](http://bioinformatics.oxfordjournals.org/content/19/8/999.short) (2003) by Lee

It was created because of neccessity for having fast POA implementation which provides programming interface for generating consensus sequence (so that programmer doesn't have to store intermediate outputs of algorithm in a file).

The tool should be compatible with most UNIX flavors and has been successfully test on the following operating systems:

- Mac OS X 10.10.3
- Ubuntu 14.04 LTS

## Requirements
- g++ (4.8.2. or higher)
- [GNU Make][2]
- [Doxygen][1] (optional)

## Installation

To install the CPPPOA run the following commands from the folder where you want to install the tool:

	git clone https://github.com/mculinovic/cpppoa.git
	cd cpppoa/
	make

Running the `make` command without arguments will build the release version of the tool and the example binary file `release/cpppoa`.

To build the debug version of the tool use:

	make debug

To build both the debug and release versions use:

	make all

To delete all files generated during the build process, for both debug and release, use:

	make clean

## Documentation

The documentation for this tool was written to work with the doxygen documentation generator. To successfully generate the documentation, the doxygen executable must be in your `PATH` variable.

To create the documentation in HTML and LaTeX format run the following command from the root of the tool:

	make docs

HTML documentation is placed in `docs/html`, while the LaTeX documentation is placed in `docs/latex`. To view the HTML documentation open `docs/html/index.html` in any web browser. The PDF documentation is obtainable by compiling the generated LaTeX code with the provided makefile.

Use the following commands from the root of the project to create the PDF version of the documentation:

	cd docs/latex/
	make
	open refman.pdf


## Usage

To use cpppoa as programming interface just do following - include poa.hpp and compile it together with cpppoa:


```
...
#include "poa.hpp"
...
```
```
...
string seq("PESLLYGRFTIESDVW");
string seq2("PEAALYGRFTIKSDVW");
vector<string> sequences;
sequences.emplace_back(seq);
sequences.emplace_back(seq2);

string consensus = poa_consensus(sequences);
...
```

## Contributors

- [Marko Culinovic](marko.culinovic@gmail.com)

[1]: http://www.stack.nl/~dimitri/doxygen/ "Doxygen"
[2]: http://www.gnu.org/software/make/ "GNU Make"
