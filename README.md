# Summary

This package is a library and some associated scripts and executables
for performing speech recognition using conditional random fields and
segmental conditional random fields. It has been developed primarily
with research in mind, as opposed to production speech recognition.

# Installation

## Dependencies

- OpenFST (dynamically linked)
- Quicknet3 (statically linked)
- cblas (dynamically linked)
- Kaldi installation (demo)
- TIMIT dataset (demo)
- autoconf 2.68 (development build)
- automake 1.15 (development build)
- gnulib/gnulib-tool (development build)

## For End Users

A source distribution packaged with autotools is available for download
[here](https://cse.osu.edu/~stiff/software/asr-craft-1.1.tar.gz). Download it to
a location of your choice and deflate (e.g. with `tar xvfz`).

From the top-level directory (containing this README and the script
`configure`), the basic installation is as simple as:

~~~~
./configure
make
make install
~~~~

Note, however, the following caveats:

Compilation depends on three external libraries for linking:
[OpenFST](http://www.openfst.org/twiki/bin/view/FST/WebHome), 
[Quicknet3](http://www1.icsi.berkeley.edu/Speech/qn.html), and
[cblas](http://www.netlib.org/blas/). If these are not installed in
standard locations (e.g. /usr/lib/) on your target system, the configure
script provides flags as a convenience for specifying their location,
`--with-libquicknet3-prefix=DIR` and `--with-libfst-prefix=DIR` (see
`./configure --help` for details).

This package also comes with an optionally-installed demonstration
of its Segmental CRF capabilities, for the TIMIT phone-recognition
task. This demonstration depends on the existence of a functional
installation of the [Kaldi toolkit](http://kaldi-asr.org/) and the TIMIT
dataset. If the location of the Kaldi TIMIT example recipe directory is
provided at configure time via the `--with-kaldi-timit-s5=DIR` argument,
the demonstration scripts will be installed. *NOTE: this will silently
overwrite the run.sh installed in that directory. Any important work in
that file should be backed up before installation of this package.*
After installation, you will have to modify the run.sh script to point
to the correct location of your TIMIT data directory. With the demo
scripts installed, you can then run the script `scrf-timit-demo.sh` (which
should have been installed to ${prefix}/bin). `scrf-timit-demo.sh` takes
a single argument which is the name of a directory of your choosing,
which will be used to store files generated over the course of the
experiment.

## For Developers

This package has been built using autoconf version 2.68 (or higher) and
automake version 1.15. You may be able to build it with older software, but
it may take some coaxing. Some m4 macros were undefined when attempting to
build using automake version 1.11, for example.

Dependencies:
ASR-CRaFT statically compiles dependencies into a library and executables.
It requires the QuickNet3 library (see above), CBLAS (for
example, via ATLAS), and OpenFST (see above) to build. If these libraries
and headers are installed in non-standard locations, locations can be
provided to the configure script as described above for end-users.

The configure script depends on [gnulib](https://www.gnu.org/software/gnulib/)
for a couple of m4 macros to facilitate finding library dependencies and
setting appropriate rpath variables for the compiled executables. Thus,
for development you will have to install gnulib-tool and invoke it as
described in ${topsrcdir}/m4/gnulib-cache.m4. [Note that the gnulib source
code contained in this package is distributed under the GPL in accordance
with the [Autoconf exception](https://www.gnu.org/licenses/autoconf-exception-3.0.en.html)].

### To build

Ideally, you should be able to:

1. Clone the source
2. Run `autoreconf` in the top-level directory
   (you may need to run `automake --add-missing` first).
3. Run `./configure`
   (see above about non-standard library and header locations; also, the
   .gitignore file is set up to allow a parallel build from a ./build/
   directory. If you want to keep your source tree free of object files,
   you may want to create the build/ directory, and from build/, run
   `../configure [options]`).
4. Run `make`.
5. Run `make install`.
6. Optionally, run `make dist` to package a new source distribution.