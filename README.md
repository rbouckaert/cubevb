# CubeVB

This is a package for [BEAST 2](http://beast2.org) for variational Bayes on a cube.

It also contains a Matrix based summary tree topology provider that works with TreeAnnotator in BEAST v2.7.6.

## Build

Get code for beast2, BeastFX and cubevb. Then run

```
ant install
```

to install the package.

## Paper

Remco R. Bouckaert.
Variational Bayesian Phylogenies through Matrix Representation of Tree Space
2023 (awaiting preprint ling).

## Data

Data used in the paper:


* [tree.tgz](https://github.com/rbouckaert/cubevb/releases/download/v1.0.0/trees.tgz) : trees used in Section 3 of the paper
* [wcss.tgz](https://github.com/rbouckaert/cubevb/releases/download/v1.0.0/wcss.tgz) : ground truth + BEAST XML + summary log files for well calibrated simulation study of Section 3.2 of the paper
* [performance.tgz](https://github.com/rbouckaert/cubevb/releases/download/v1.0.0/performance.tgz) : ground truth + BEAST XML + summary log for performance measurement used in section 3.3 of the paper
