GenSolFlux
----------

GenSolFlux is a simple program to generate N-tuples of boosted dark matter momenta coming from the Sun. Solar positions are generated using the SolTrack package over a user-specified time range. For each solar position, the corresponding dark matter momentum coming from that position to a detector located at a user-specified latitude and longitude is calculated in the detector reference frame. The output is a GENIE-readable ROOT file containing the N-tuples of dark matter momenta.

Prerequisites
-------------

GENIE v3 with --enable-boosted-dark-matter configure option
SolTrack

Installation
------------

Installation proceeds as follows from the source directory

$ ./confiugre --with-soltrack-inc=<path-to-soltrack-include> --with-soltrack-lib=<path-to-soltrack-library>

If the SolTrack directories are not specified, they will be searched for.

$ make
$ make install

Author
------

Joshua Berger (josh.berger at pitt.edu)
Based in part on files from the GENIE software suite

License
-------

GPL v3
