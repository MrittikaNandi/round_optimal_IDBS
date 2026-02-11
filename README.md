This repository contains the scripts accompanying the article

**Round-Optimal Identity-Based Blind Signature from Module Lattice
Assumptions**

# Description of files
Short description of the content:
* MLWE_security.py, model_BKZ.py, MSIS_security.py, proba.py are MSIS/MLWE security estimator scipts borrowed from  [BDK+18]
* signature_params.py is the main file that outputs signature size and the concrete bit securities. scripts are borowed from [AKSY22] with necessary changes

The inline comments have the detailed description relating to the corresponding parameters and where they come from. 

The estimates for the Zero-knowledge proof are based on

* [LNP22] Lattice-Based Zero-Knowledge Proofs and Applications: Shorter, Simpler, and More General _Vadim Lyubashevsky and Ngoc Khanh Nguyen and Maxime Plancon_ https://eprint.iacr.org/2022/284.pdf

The estimates on the entropy of a discrete Gaussian variable is taken from

* [ESLR22] Efficient Hybrid Exact/Relaxed Lattice Proofs and Applications to Rounding and VRFs _Muhammed F. Esgin and Ron Steinfeld and Dongxi Liu and Sushmita Ruj_ https://eprint.iacr.org/2022/141.pdf


The security estimates for module SIS, moduel LWE, and NTRU are borrowed from
* [BDK+18] CRYSTALS-Dilithium – Algorithm Specifications and Supporting Documentation (https://pq-crystals.org/dilithium/resources.shtml) and the accompanying scripts  
(https://github.com/lducas/leaky-LWE-Estimator/tree/NIST-round3/NIST-round3)

The parameters with necessary modifications follow from:

* [AKSY22] Practical, Round-Optimal Lattice-Based Blind Signatures _Shweta Agrawal and Elena Kirshanova and Damien Stehlé and Anshu Yadav_
(https://eprint.iacr.org/2021/1565) and the accompanying scripts (https://gitlab.com/ElenaKirshanova/onemoresis_estimates/) 

# Requirements

* Install the module 'sympy' if not already present
```
pip install sympy
```
* [Python 3] (https://www.python.org/downloads/)

# Experiments
To compute signature sizes and security estimates run
```
python3 signature_params.py
```
