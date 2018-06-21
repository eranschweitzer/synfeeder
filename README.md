# Distribution Feeder Generation

This code repository automatically generates distribution feeders.
The data comes from a medium voltage network in the Netherlands.
As a result, at present all models are assumed balanced three phase, and underground cables as conductors. 

Please cite [\[1\]][1] when using this tool.

## File Structure
- `src` The main functions used to to generate the feeders. Of these the key ones are
  - `single_feeder_gen` produces a single radial feeder
  - `feeder_main` is a script to produce many feeders using an input csv file
- `data` contains various distributions, limiting functions, libraries and similar data needed by the algorithm
- `CIM` contains some preliminary code to convert csv output of feeders to CIM documents
- `docs` contatins documentation

## Installation
Simply add the `src` folder to your Matlab path.

## Usage
The basic functionality creates two Matlab structures, `n` and `e`, with information about the nodes and buses of the system.
```matlab
  >> [n,e] = single_feeder_gen(N, Stotal, Pinj_total);
```
If no inputs are passed a KDE based on the data is used to sample inputs. 
The function `inputs_sample()` is provided as a convinience to generate samples.
```matlab
  >> [N, Stotal, Pinj_total] = inputs_sample(n, use_pinj);
```
Here `n` is the number of desired samples.
Boolean `use_pinj` determines whether injections should be considered, if it is `false` then `Pinj_total = 0` always.

*Note:* Since `Pinj_total` is a fixed net injection (rather than a combination of load and generation) it is not particularly flexible. 
As such we do not necessarily advise using it. 
In the case where `single_feeder_gen()` is called with no arguments, `Pinj_total` is set to zero.

### Matpower Format
Structures `n` and `e` can be converted to a [MATPOWER][2] case, `mpc`, using the `matpower_fmt()` function:
```matlab
  >> mpc = matpower_fmt(n,e,freq);
```
The third argument, `freq`, is the nominal system frequency, which is needed to convert the line capacitance from &mu;F to per-unit susceptance.
If it is not provided 50 Hz is used, since this was the default frequency in the system the data came from.
The feeder powerflow can then be easily solved (assuming MATPOWER is installed):
```matlab
  >> r = runpf(mpc);
```

The created feeders are radial, which occasionally might cause convergence problems with the normal Newton-Raphson algorithm (so far solutions appear to be stable).
Matpower comes with a few radial algorithms specifically for these situations.
To use the current summation method, for example use:
```matlab
  >> mpopt = mpoptions;
  >> mpopt.pf.alg = 'ISUM';
  >> mpopt.pf.radial.max_it = 500;
```
From experience we recommend increasing the iteration number to greater than the default 20.

While the algorithm produces radial systems, parallel elements can be produced.
To check if any parallel elements exists try `any(e.num_parallel > 1)`.
The utility function `parallel_branch_join()` is provided to combine parallel branches to enable use of the radial algorithms.

A simple set of commands to create and solve a feeder using a radial powerflow is:
```matlab
  >> [n,e] = single_feeder_gen();
  >> mpc = matpower_fmt(n,e);
  >> mpc2 = parallel_branch_join(mpc);
  >> r = runpf(mpc2, mpopt);
```
## Documentation
There is (as of yet) no user manual. Instead the [FEN-report][3] in `docs` as well as reference [\[1\]][1] provide a fairly detailed overivew of the work.

## Known Issues
- The distribution transformer in real systems generally has taps and varies these taps to get the load on the feeder to a desired point.
Currently the tap entry in the `branch` matrix is determined based on the ratio between the system voltages (HV=110 kV, MV=10 kV) and the nominal voltages of the transformer.
The voltage behavior over the transfomer is the least predictable in the output.
Eventually, there should probably be a step to come up with good settings but currently there are two quick ways to play with the set point: 
  1. **Tap Setting** Reducing the tap setting from 1 to 0.98 or so will help raise the voltage on the low voltage side. Conversely, increasing the tap from 0.92 to 0.94 will lower the low voltage terminal.
  2. **Reactive Support** Adding a shunt reactance to support the voltage at the low-voltage bus can help raise or lower the voltage.
For example, take a look at `e.qdownstream(1)` which is roughly the reactive power in the transformer in MVAr, and add a fraction of this to `mpc.bus(2,BS)` if trying to rais the low voltage terminal.

## To Do
- Options argument to change some of the defaults in the `single_feeder_gen` function.
- Better interface to the various distributions so that new ones, possibly even non-parametric, can be used in the future.
- Handle tap and reactive support options automatically, or as an option.
- Add a module that connects multiple feeders with normally-open branches (partially started, see [FEN-report][3]).
- Time-series modeling.
- Translation to more modeling languages such as GridLAB-D (started already) and OpenDSS.
- Include modeling of single phases rather then only balanced 3-phase. 
- ...

## Publications
1. E. Schweitzer, A. Scaglione, A. Monti and G. A. Pagani, 
 "Automated Generation Algorithm for Synthetic Medium Voltage Radial Distribution Systems," 
  IEEE Journal on Emerging and Selected Topics in Circuits and Systems, vol. 7, no. 2, pp. 271-284, June 2017.

[1]: https://dx.doi.org/10.1109/JETCAS.2017.2682934
[2]: http://www.pserc.cornell.edu/matpower/
[3]: ./docs/FEN_report.pdf
