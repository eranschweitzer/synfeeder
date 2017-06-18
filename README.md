#Distribution Feeder Generation

This is a collection of code to automatically generate distribution feeders

##File Structure
- `src` The main functions used to to generate the feeders. Of these the key ones are
  - `single_feeder_gen` produces a single radial feeder
  - `feeder_main` is a script to produce many feeders using an input csv file
- `data` contains various distributions, limiting functions, libraries and similar data needed by the algorithm
- `CIM` contains some preliminary code to convert csv output of feeders to CIM documents
- `docs` contatins documentation

## Installation
Simply add the `src` folder to your Matlab path.

## Documentation
There is (as of yet) no user manual. Instead the FEN-report in `docs` as well as [reference \[1\]][1] 

> E. Schweitzer, A. Scaglione, A. Monti and G. A. Pagani, 
 "Automated Generation Algorithm for Synthetic Medium Voltage Radial Distribution Systems," 
  IEEE Journal on Emerging and Selected Topics in Circuits and Systems, vol. 7, no. 2, pp. 271-284, June 2017.

[1]: https://dx.doi.org/10.1109/JETCAS.2017.2682934
