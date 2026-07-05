# Results

Every number below is reproduced by `run_analysis.py` and asserted by `verify.py` (23/23 checks pass). Reference apparatus: 100 nm fused-silica sphere, T = 100 mK, integration time τ = 10⁵ s, unless noted.

## 1. Elastic eigenmodes (ℓ=2 fundamental, R = 100 nm)
| Material | ν | η = ω₀R/v_t | f₀ | α |
|---|---|---|---|---|
| Silica | 0.17 | 2.628 | 15.7 GHz | 0.837 |
| Silicon | 0.04 | 2.602 | 24.2 GHz | 0.820 |
| Sapphire | 0.29 | 2.645 | 25.4 GHz | 0.846 |
| Diamond | 0.10 | 2.615 | 50.0 GHz | 0.829 |

Isotropic-equivalent elastic constants (crystals are anisotropic → frequencies good to ~10%). Eigenvector A/B = −2.056 for silica; boundary conditions vanish to machine precision; induced quadrupole verified trace-free (Q_xx = Q_yy = −Q_zz/2).

## 2. Tidal coupling and the channel ratio
- Geometric factor **α = 0.837** (Frobenius), 0.683 (zz) — order unity, from the eigenfunction.
- COM two-phonon Frobenius factor **D = 10** exactly (= ⅔⟨r⁴⟩ = ⅔·15), verified by direct Fock-space summation.
- **Channel ratio |Q_body|²/|Q_COM|²** at ω_COM = 10⁵ rad/s is O(10²), convention-dependent by mode counting:
  - single-mode (m=0 internal vs 1D COM, D=4/3): **93**
  - full-channel (5 ℓ=2 m-modes vs 3D COM, D=10): **62**
  - naive (D=1): 124
  Advantage is robustly ~10²; quote the convention, not a false-precision "90."
- Spectral ratio (ω₀/ω_COM)ⁿ: 10⁶ (n=1), **10³⁰ (n=5)**.
- Combined advantage (single, n=5) ≈ **9×10³¹ ≈ 10³²**.
- Rate scalings: Γ_COM ∝ ω³ (recovers Toroš), Γ_internal ∝ ω⁴.

## 3. First finite bound on S_EE at 16 GHz (headline result)
`Γ_min = √(γ n̄_th/τ)`, γ = ω₀/Q (thermal-limited, no clamping loss); `S_EE^min = Γ_min ħ²/|Q_body|²`.

| Q | S_EE^min (s⁻³) |
|---|---|
| 10⁴ | 7.4×10⁷ |
| **10⁵** | **2.33×10⁷** |
| 10⁶ | 7.4×10⁶ |

- Vacuum reference `S_EE^vac = t_Pl² ω₀⁵ = 2.74×10⁻³² s⁻³`.
- **Gap: bound/vacuum = 8.5×10³⁸ ≈ 10³⁹.** The bound therefore does *not* test the linearized-gravity vacuum; it is the first finite constraint in the band and bounds enhanced-noise states (squeezed / thermal-graviton / large-amplitude).

## 4. No-go (reach never closes the vacuum gap)
Reach ∝ R⁻⁶, band ω₀ ∝ 1/R ⇒ reach/vacuum ratio ∝ ~1/R:
| R | f₀ | gap/vacuum |
|---|---|---|
| 100 nm | 16 GHz | 8.5×10³⁸ |
| 1 mm | 1.6 MHz | 1.4×10³⁶ |
| 100 mm | 16 kHz | 1.4×10³⁴ |
| 300 mm | 5.2 kHz | 4.5×10³³ |

Even rebuilding a bar detector (R ≈ 30 cm, f₀ ≈ 5 kHz) leaves the gap at ~10³³: **no levitated or bar-scale mesoscopic experiment closes the near-field tidal vacuum gap.** State this — it pre-empts the fatal objection.

## 5. Material discriminant (fixed R, normalized to silica)
| channel | silica | silicon | sapphire | diamond |
|---|---|---|---|---|
| tidal n=1 (floor) | 1.0 | 1.06 | 1.81 | 1.60 |
| tidal n=3 | 1.0 | 2.5 | 4.7 | 16.1 |
| tidal n=5 | 1.0 | 5.9 | 12.4 | **163** |
| thermal loss | 1.0 | 2.6×10⁻² | 1.5×10⁻² | **2.3×10⁻⁷** |
| gas (sudden) | 1.0 | 0.61 | 0.34 | 0.20 |

- Γ_tidal ∝ ρ ω₀^{n−1}. At n=1 this is the density ratio (tracks ρ, not f₀ — sapphire highest); for n ≥ 1 it is flat-or-rising with ω₀.
- Intrinsic backgrounds fall steeply. log-space Pearson(tidal n=5, thermal loss) = **−0.970**; (tidal, gas) = **−0.976** — near-orthogonal.
- **The sign of the material trend is the discriminant and rests on n ≥ 1 alone.** Diamond vs silica: tidal says diamond ~163× hotter, thermal loss says diamond ~4×10⁶× colder — a GHz mode that heats faster in diamond than silica cannot be thermal loss.

## 6. Cutoff threshold (silica = diamond sign inversion)
For S_EE ∝ ωⁿ e^{−ω/ω_c}: `ω_c*(n) = (ω_di − ω_si)/ln[(ρ_di/ρ_si)(ω_di/ω_si)^{n−1}]`.
| n | ω_c* (rad/s) |
|---|---|
| 1 | 4.6×10¹¹ |
| 3 | 7.7×10¹⁰ |
| 5 | **4.22×10¹⁰** |
| 6 | 3.4×10¹⁰ |

Earth c/r_s (r_s = 8.87 mm) = 3.38×10¹⁰ rad/s, just below ω_c*(n=5). The measured silica–diamond sign thus bounds a phenomenological cutoff scale against ~4×10¹⁰ rad/s. (Treat the cutoff phenomenologically; a corpuscular interpretation is one branch, not a premise — 𝒞_Earth ~ 10⁻⁹.)

## 7. Readout ε-threshold
Mayor et al. 2025 cooperativity C ≈ 1400 ⇒ **ε ≥ 1/√C ≈ 2.7%**. The mode stays single-phonon-resolved for ℓ=2-to-transducer overlap above ~3%. Caveats: (i) Mayor's C ≈ 1400 is the *classical* cooperativity at a 3 K operating point; the same device's quantum cooperativity is C/n̄_th ≈ 175 ⇒ threshold ≈ 8%. At 100 mK, n̄_th < 1 and the two coincide, restoring ~3%. (ii) Mayor's C is an optical readout of a *breathing* (volume-changing) mode; the ℓ=2 is volume-preserving, so a strain-coupled or polarizability-dispersive readout is required, and its overlap is the one open parameter.

## Key values
- η = 2.628, f₀ = 15.7 GHz, α = 0.84 (silica, 100 nm)
- channel ratio ~10² (62 full / 93 single-mode); combined advantage ~10³² at n=5
- **first finite bound S_EE^min ≈ 2.3×10⁷ s⁻³ (Q=10⁵), 10³⁹ above vacuum, gap never closes**
- material discriminant: tidal rises (diamond 163× at n=5), backgrounds fall (Pearson −0.97); sign rests on n ≥ 1
- cutoff threshold ω_c*(n=5) = 4.2×10¹⁰ rad/s; ε-threshold 2.7% (classical C ≈ 1400; quantum C_q ≈ 175 ⇒ ~8%; coincide at 100 mK)
