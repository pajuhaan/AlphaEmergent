# Emergent Fine-Structure Constant Calculation: AlphaEmergent


This repository contains the code to compute the fine-structure constant $\alpha$ as an emergent and parameter-free invariant within Relator theory. The framework solves a closed root equation that couples the Coulombic (scalar) and inductive (vector) channels in a unified C ⊕ R3 geometry. This computation does not rely on any measured constants like $e$, $c$, or $\hbar$, and produces a result matching the CODATA 2022 $\alpha$ value with an extraordinary precision.

The method relies on a geometric and gauge-invariant construction that predicts $\alpha$ without any empirical inputs. It is based on a unique approach to coupling quantum electrodynamics and relativistic mechanics using Relator geometry.

For more details, see the full paper: [Emergent Alpha](https://zenodo.org/records/16944533).

## Features

* **Precision**: Sub-ppb accuracy matching the CODATA 2022 value of $\alpha$.
* **No Fitted Constants**: The value of $\alpha$ is emergent from geometry, with no empirical tuning.
* **Pure Geometry**: Derived from Relator theory and geometry, offering a novel approach to the fine-structure constant.
* **Self-Contained**: No need for external constants such as $e$, $c$, or $\hbar$.

## How It Works

This framework implements the closed root equation for $\alpha$, coupling the Coulombic and inductive (Λ) sectors. The main steps include:

1. **Coulombic Base Calculation**: Using the scalar channel to compute the baseline for $\alpha$.
2. **Inductive Channel Correction**: Incorporating the effects of the vector channel via a logarithmic correction.
3. **Geometric Closure**: Ensuring the correctness of the derived $\alpha$ by locking the scalar and vector sectors using geometric relations.
4. **Convergence and Precision**: Iterating to convergence, with built-in checks for numerical stability and precision.

## Installation

To run this code, you will need Python 3.x and the following libraries:

* `mpmath` (for arbitrary precision arithmetic)
* `numpy` (for mathematical operations)

Install dependencies using pip:

```
pip install mpmath numpy
```

## Usage

After installing the necessary libraries, you can run the code using:

```bash
python AlphaEmergent.py
```

The script will compute the emergent $\alpha$ and output the results to the console.

## Results

This code computes the emergent value of the fine-structure constant and outputs it alongside the associated uncertainties. The precision achieved by this approach is comparable to the latest CODATA values.

### Example Output:

```
alpha_emergent     = 0.007297352564332633809798
alpha_em^-1        = 137.0359991769773
Λ_eff (final)      = 0.6916840202847290215451
K (spectral)       = 0.002231538916531970186409
P^(IR)_χ(ℓ0)       = 0.08577919258455560110975
∆Λ_OUT (η0)        = -0.01396715806205758007151
∆Λ^(UV→IR)         = 0.05448534958655209065325
C_log(α_em)        = 0.333333333333333333   (∆ vs 1/3 = -2.71498346735e-30)
[Context] α_ref    = 0.007297352564311  →  ∆α(ppb) = 0.00296461074167
```


## Credits

This work is based on the paper *Emergent Fine-Structure Constant from Relator Geometry* by Mehrdad Pajuhaan.
