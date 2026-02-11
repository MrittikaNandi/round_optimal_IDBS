This repository contains the scripts accompanying the article

**Round-Optimal Identity-Based Blind Signature from Lattice
Assumptions**


The estimates for the Zero-knowledge proof are based on

* [LNP] Lattice-Based Zero-Knowledge Proofs and Applications: Shorter, Simpler, and More General _Vadim Lyubashevsky and Ngoc Khanh Nguyen and Maxime Plancon_ https://eprint.iacr.org/2022/284.pdf

The estimates on the entropy of a discrete Gaussian variable is taken from

* Efficient Hybrid Exact/Relaxed Lattice Proofs and Applications to Rounding and VRFs _Muhammed F. Esgin and Ron Steinfeld and Dongxi Liu and Sushmita Ruj_ https://eprint.iacr.org/2022/141.pdf


Our security estimates for module SIS, moduel LWE, and NTRU are borrowed from
* [Dil] CRYSTALS-Dilithium – Algorithm Specifications and Supporting Documentation (https://pq-crystals.org/dilithium/resources.shtml) and the accompanying scripts  
(https://github.com/lducas/leaky-LWE-Estimator/tree/NIST-round3/NIST-round3)

Our parameters with necessary modifications follow from:

* [AKSY 22] Practical, Round-Optimal Lattice-Based Blind Signatures _Shweta Agrawal and Elena Kirshanova and Damien Stehlé and Anshu Yadav_
(https://eprint.iacr.org/2021/1565) and the accompanying scripts (https://gitlab.com/ElenaKirshanova/onemoresis_estimates/) 

# Requirements

* [Python 3] (https://www.python.org/downloads/)

# Description of files
Short description of the content:
* MLWE_security.py, model_BKZ.py, MSIS_security.py, proba.py are MSIS/MLWE security estimator scipts borrowed from  [Dil]
* signature_params.py is the main file that outputs signature size and the concrete bit securities. Please read the inline comments to understand where the figures come from. scripts are borowed from [AKSY 22] with necessary changes


# Experiments
To compute signature sizes and security estimates run
```
python3 signature_params.py
```
