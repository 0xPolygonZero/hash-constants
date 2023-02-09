# Generation of MDS matrices and precomputed constants for the Poseidon hash function

This repository consists of scripts for finding secure MDS matrices
whose entries admit fast computation, as well as for generating the
precomputed constants for such matrices.

**Caveat emptor:** None of these files were written with the
expectation of release or even re-use; as such they well below the
usual standard one might expect of a public repository, and we have no
intention to improve them.

Nevertheless, please feel free to raise any issues encountered --
especially bugs! -- and we will do our best to address them.

This repository contains five files:

- `mds_search.sage`: Functions for discovering MDS matrices; examples
  for the Crandall and Goldilocks fields are at end of file.
  
- `mds_fft_search.sage`: Partly copypasta from `mds_search.sage`
  adapted to find MDS matrices with nice FFT values.

- `calc_round_numbers.py`: Used to compute the number of full and
  partial rounds needed to obtain sufficient security. Will calculate
  for Crandall and Goldilocks fields by default; specify a prime on
  the command line to use that instead.

- `poseidon_constants.sage`: Used to produce precomputed constants for
  the 'fast' Poseidon implementation in Rust. Produces values for
  Crandall field by default; specify 'goldilocks' on the command line
  to get those instead. Uses width 8 by default, specify 'crandall 12'
  or 'goldilocks 12' on the command line for width 12. Edit 'main'
  function directly to use a different MDS matrix.

- `mds_security.sage`: Checks the security of a given MDS matrix wrt
  to several tests; primarily to avoid 'infinitely long subspace trails'.

The third and fourth files are originally from
[Hadeshash](https://extgit.iaik.tugraz.at/krypto/hadeshash) reference
implementation of Poseidon. They have been _substantially_ cleaned up,
as well as modified for our particular case.

The last file contains functions originally from
[linear-layer-tool](https://extgit.iaik.tugraz.at/krypto/linear-layer-tool)
by the authors of Poseidon and its reference implementation. It has
also been cleaned up and adapted to our case.
