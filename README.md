# singular-optics-phase-retrieval
MATLAB implementation of three-plane phase retrieval for optical vortex fields using Laguerre-Gaussian beams, partial Fourier transforms, and Fresnel propagation. Includes Gerchberg-Saxton reconstruction, twist parameter estimation, and support for noisy intensity measurements.

# Three-Plane Phase Retrieval using Partial Fourier Transform and Fresnel Propagation

MATLAB implementation accompanying the work:

**P. A. A. Yasir and J. S. Ivan,**
*Estimation of phases with dislocations in paraxial wave fields from intensity measurements*,
**Physical Review A** **97**, 023817 (2018).

https://doi.org/10.1103/PhysRevA.97.043838

This repository implements a three-plane phase retrieval algorithm for paraxial optical fields carrying orbital angular momentum (OAM). The method combines a **partial Fourier transform**, **Fresnel propagation**, and a **Gerchberg-Saxton (GS)** iterative reconstruction algorithm to recover the phase of a complex wave field from intensity-only measurements.

---

## Physical model

The optical system follows the sequence

```text
Laguerre-Gaussian superposition
            │
            ▼
  Partial Fourier Transform
            │
            ▼
    Rotated HG-like field
            │
            ▼
    Fresnel propagation
            │
            ▼
   Intensity measurements at
       three transverse planes
            │
            ▼
  Gerchberg-Saxton reconstruction
            │
            ▼
   Twist parameter estimation
```

The partial Fourier transform naturally creates or annihilates phase dislocations and allows the orientation of optical vortices to be distinguished. The reconstruction algorithm uses intensity information from three transverse planes to recover the original complex field.

---

## Features

* Generation of Laguerre-Gaussian (LG) beams.
* Random superpositions of LG modes.
* Partial Fourier transform implemented through an explicit DFT matrix.
* Fresnel propagation using the transfer-function method.
* Three-plane Gerchberg-Saxton phase retrieval.
* Third-order finite-difference derivatives.
* Twist parameter (orbital angular momentum-related quantity) evaluation.
* Addition of Gaussian noise and circular aperture masking.
* Reproducible simulation pipeline suitable for research and teaching.

---

## Repository structure

| File                            | Description                                                |
| ------------------------------- | ---------------------------------------------------------- |
| `Demo_3plane_phase_retrieval.m` | Main demonstration script.                                 |
| `LG_beam.m`                     | Generates normalized Laguerre-Gaussian modes.              |
| `random_lg_superposition.m`     | Creates random superpositions of LG modes.                 |
| `lg_mode_table.m`               | Generates the mapping between `(j,m)` and `(l,p)` indices. |
| `laguerre.m`                    | Associated Laguerre polynomial evaluation.                 |
| `laguerre_matrix.m`             | Matrix implementation of associated Laguerre polynomials.  |
| `PFTransform.m`                 | Partial Fourier transform.                                 |
| `IPFTransform.m`                | Inverse partial Fourier transform.                         |
| `fresnel_propagation.m`         | Fresnel propagation using the transfer-function approach.  |
| `gerchberg_saxton_3plane.m`     | Three-plane Gerchberg-Saxton reconstruction algorithm.     |
| `replace_amplitude.m`           | Replaces field amplitude while preserving phase.           |
| `third_order_gradient.m`        | Third-order finite-difference derivatives.                 |
| `twist_parameter.m`             | Computes the twist parameter of a complex field.           |
| `make_grid.m`                   | Generates the computational grid.                          |
| `circ_aper.m`                   | Generates a circular aperture mask.                        |
| `normalize_intensity.m`         | Intensity normalization utility.                           |
| `add_gaussian_noise.m`          | Adds Gaussian noise with a specified SNR.                  |

---

## Requirements

* MATLAB R2018b or later (earlier versions may also work).
* No external toolboxes are required beyond standard MATLAB functionality.

---

## Running the code

Clone the repository and run

```matlab
Demo_3plane_phase_retrieval
```

The script:

1. Generates a random superposition of LG modes.
2. Applies the partial Fourier transform.
3. Propagates the field using Fresnel diffraction.
4. Simulates intensity measurements at three planes.
5. Adds optional Gaussian noise and aperture truncation.
6. Performs three-plane Gerchberg-Saxton reconstruction.
7. Computes and compares the twist parameters of the original and reconstructed fields.
8. Displays intensity, phase, and convergence plots.

---

## Typical output

The demonstration script produces:

* Input and reconstructed intensity distributions.
* Input and reconstructed phase distributions.
* Iterative convergence of the GS algorithm.
* Twist parameter comparison between the original and reconstructed fields.
* Final intensity and amplitude correlation statistics.

---

## Notes

* The implementation follows the conventions and notation used in the accompanying publication.
* The inverse partial Fourier transform is implemented as the Hermitian adjoint of the forward transform matrix.
* The code is intended primarily for numerical investigations of singular optics and phase retrieval problems.

---

## License

This repository is released under the MIT License. See the `LICENSE` file for details.
