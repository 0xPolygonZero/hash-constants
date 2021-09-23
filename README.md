# Generation of MDS matrices and precomputed constants for the Poseidon hash function

This repository contains three files:

- `mds_search.sage`: Functions for discovering MDS matrices; examples
  for the Crandall and Goldilocks fields are at end of file.

- `calc_round_numbers.py`: Used to compute the number of full and
  partial rounds needed to obtain sufficient security. Will calculate
  for Crandall and Goldilocks fields by default; specify a prime on
  the command line to use that instead.

- `poseidonperm_x3_64_24_optimized.sage`: Used to produce precomputed
  constants for the 'fast' Poseidon implementation in Rust. Produces
  values for Crandall field by default; specify 'goldilocks' on the
  command line to get those instead. Uses width 8 by default, specify
  'crandall 12' or 'goldilocks 12' on the command line for
  width 12. Edit 'main' function directly to use a different MDS matrix.

The latter two files are originally from
[Hadeshash](https://extgit.iaik.tugraz.at/krypto/hadeshash) reference
implementation of Poseidon. They have been _substantially_ cleaned up,
as well as modified for our particular case.
