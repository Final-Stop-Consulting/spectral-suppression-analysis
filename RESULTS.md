# Spectral Suppression — Analytical Results (ℓ=2 coupling)

Results of the fundamental ℓ=2 Lamb tidal-coupling calculation. Every value below is reproduced by `run_analysis.py` and cross-checked by `verify.py`.

## Setup and assumptions
- Free, homogeneous, **isotropic** elastic sphere; fundamental ℓ=2 spheroidal ("quadrupolar") Lamb mode.
- Displacement from two scalar potentials (compressional j_ℓ(hr), shear-poloidal j_ℓ(kr)), with
  `U(r)=A h j_ℓ'(hr)+B ℓ(ℓ+1) j_ℓ(kr)/r`, `V(r)=A j_ℓ(hr)/r+B[j_ℓ(kr)/r+k j_ℓ'(kr)]`, h=ω/v_l, k=ω/v_t.
- Traction-free surface (σ_rr=σ_rθ=0 at R) → 2×2 secular determinant.
- Reference particle: fused silica, R=100 nm, v_l=5970, v_t=3760 m/s, ρ=2200 kg/m³ (Poisson ν=0.171).
- Long-wavelength (uniform-field) coupling assumed and justified: λ_grav≈2 cm ≫ R.

## 1. Eigenvalue (the Lamb frequency)
- Reduced eigenfrequency **η ≡ ω₀R/v_t = 2.6282**.
- **ω₀ = 9.88×10¹⁰ rad/s, f₀ = 15.73 GHz**.
- Eigenvector **A/B = −2.0559**, fixed by both boundary conditions independently (agree to 4×10⁻⁸ → the root is genuine, not a determinant artifact).

## 2. Geometric factor α
Induced symmetric-trace-free mass quadrupole δQ_ij = ρ∫(u_i x_j + x_i u_j − ⅔δ_ij x·u)d³x, mode normalized to effective mass = M, so |⟨Q^body⟩| = α M R x_zpf with x_zpf=√(ħ/2Mω₀).

- **α = 0.837** (Frobenius / basis-independent tensor norm)
- α = 0.683 (zz-component convention)
- → α is order unity (~0.7–0.8), computed from the actual eigenfunction.
- Quadrupole verified **trace-free analytically**: Q_xx=Q_yy=−Q_zz/2 (b₁=−a₁, b₂=a₂ in the angular integrals).

## 3. Matrix-element ratio (coupling strength, internal vs COM)
`R_mat = |Q^body|²/|Q^COM|² = (3/2) α² M R² ω_COM²/(ħ ω₀)` (one self-consistent Frobenius tensor convention for both channels; COM is the two-phonon ⟨2|z²|0⟩ matrix element).

- ω_COM = 10⁵ rad/s: **R_mat ≈ 93 (≈10²)**
- ω_COM = 10⁶ rad/s: R_mat ≈ 9.3×10³
→ internal-mode coupling is ~2 orders of magnitude stronger *per unit noise PSD*.

## 4. Spectral ratio (noise PSD, n=5)
`S_EE(ω₀)/S_EE(ω_COM) = (ω₀/ω_COM)^n`
- n=1 (correspondence-principle floor): 9.9×10⁵ ≈ **10⁶**
- n=5 (graviton bath / vacuum FDT / dimensional): 9.4×10²⁹ ≈ **10³⁰**

## 5. Combined advantage (ratio of decoherence/heating rates)
`Γ_internal/Γ_COM = R_mat × (ω₀/ω_COM)^n`  (∝ ω_COM^−3, so larger at lower trap frequency)
- **n=5, ω_COM=10⁵ rad/s: 8.75×10³¹ ≈ 10³²**
- n=5, ω_COM=10⁶ rad/s: 8.75×10²⁸ ≈ 10²⁹
- Across typical trap frequencies (ω_COM ≈ 3×10⁴–10⁵), the combined advantage spans **~10³²–10³³**.

## 6. Rate scalings (parameter-free exponents)
- |Q^COM|² ∝ ω_COM⁻² and S_EE ∝ ω⁵ ⟹ **Γ_COM ∝ ω³** (recovers Toroš et al.).
- |Q^body|² ∝ ω₀⁻¹ and S_EE ∝ ω⁵ ⟹ **Γ_internal ∝ ω⁴** (steeper; matrix-element and spectral advantages compound in the same direction).

## 7. Material spread (ℓ=2 fundamental, isotropic-equivalent constants — approximate for anisotropic crystals)
| Material | ν | η | f₀ @ R=100 nm | f₀ @ matched mass | α |
|---|---|---|---|---|---|
| Silica | 0.17 | 2.628 | 15.7 GHz | 15.7 GHz | 0.84 |
| Silicon | ~0.04* | 2.602 | 24.2 GHz | 24.7 GHz | 0.84 |
| Sapphire | 0.29 | 2.645 | 25.4 GHz | 31.0 GHz | 0.84 |
| Diamond | 0.10 | 2.615 | 50.0 GHz | 58.4 GHz | 0.84 |

*Isotropic-equivalent speeds understate Si's effective Poisson ratio; treat as order-of-magnitude. The ℓ=2 fundamentals span ~16–60 GHz, all inside the band where (for any n≥1) the tidal coupling has its maximum spectral weight.

## Verification checklist (all passed)
1. Secular root η=2.6282 stable to 200 bisections; consistent with the standard Lamb ℓ=2 spheroidal spectrum (η≈2.6 for low Poisson ratio).
2. Eigenvector satisfies σ_rr=0 and σ_rθ=0 **independently** to 4×10⁻⁸.
3. Displacement potentials are exact Navier solutions by construction; only the boundary condition is solved numerically.
4. Induced quadrupole **trace-free analytically**; axisymmetric Q_xx=Q_yy=−Q_zz/2.
5. α reproduced to 8 figures by two independent code paths.
6. Γ_COM ∝ ω³ recovers the published Toroš graviton-bath scaling (independent anchor).
7. Dimensional consistency: α, R_mat, and the combined advantage are all dimensionless.

## Key values
- η = 2.628, f₀ = 15.7 GHz (silica, 100 nm)
- α ≈ 0.84 (Frobenius)
- R_mat ≈ 90 (≈10²) at ω_COM=10⁵ rad/s
- Spectral suppression 10³⁰ at n=5
- Combined internal-mode advantage ≈ 10³² (10³²–10³³ across trap frequencies)
- Γ_internal ∝ ω⁴ vs Γ_COM ∝ ω³
