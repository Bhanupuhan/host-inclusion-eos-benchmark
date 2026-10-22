# Endpoint and incremental calculations

Moulas et al. (2023), https://doi.org/10.2475/001c.68195, provide multiple approaches. Section 4.2 examines incremental calculations, but Eq. (3) is an endpoint relation. Section 5 explicitly discusses a single-step calculation when the shear modulus is constant. Therefore an endpoint solution is not a distinction unique to this implementation.

Using the repository EoS, constant host G = 531 GPa, trapping at 5.5 GPa / 1473.15 K and final conditions 0 GPa / 298 K, `python/compare_moulas_endpoint.py` gives:

| Calculation | Inclusion pressure, GPa |
|---|---:|
| Moulas et al. Eq. (3), infinite spherical host | 0.840733784 |
| Log-volume mismatch relation, infinite spherical host | 0.840759621 |
| Author-reported 3D total-form UMAT, A/B/C | 0.8409 |
| Recorded EosFit reference | 0.842 |

The Eq. (3) comparison evaluates the published equation with the repository parameters; it is not a reproduction of the paper's examples or a run of its authors' software. It assumes an infinite spherical host and is not mathematically identical to the finite 3D small-strain model. The approximately 0.020% discrepancy does not establish which approximation contributes how much.

The total-form UMAT evaluates stress from current strain and temperature and still requires nonlinear equilibrium iterations. Its A/B/C outputs agree at the supplied four-decimal precision; they do not prove arbitrary-path agreement.

The `*Elastic` prestress model is a separate linear relaxation problem, beginning with inclusion pressure 0.997681 GPa and an unstressed host. Its 0.8468 GPa result should first be compared with its matching Lame solution (0.8466 GPa).

The legacy UMAT gives 0.8577 versus 0.7803 GPa across two step setups. The reports alone do not isolate initialisation, strain history or integration error. This discrepancy is unrelated to any test of Moulas's code. No claim of a defect in that work follows from it.
