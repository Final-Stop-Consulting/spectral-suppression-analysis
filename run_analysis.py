#!/usr/bin/env python3
"""
run_analysis.py
===============
Reproduces every number and table in the accompanying manuscript,
"Spectral Suppression of Center-of-Mass Coupling to Gravitational Tidal
Fluctuations" (2026).

Usage
-----
    python run_analysis.py
"""
import numpy as np

import lamb_mode as lm
import coupling as cp

R = 100e-9              # reference particle radius (m)
TWO_PI = 2 * np.pi


def banner(s):
    print("\n" + s + "\n" + "-" * len(s))


def main():
    banner("1. Fundamental l=2 Lamb mode (fused silica, R = 100 nm)")
    g = cp.geometric_factor(lm.SILICA, R)
    f0 = g["omega0"] / TWO_PI
    print(f"  Poisson ratio nu     = {lm.SILICA.poisson:.4f}")
    print(f"  reduced freq  eta    = {g['eta']:.4f}      (omega0 R / v_t)")
    print(f"  omega0               = {g['omega0']:.4e} rad/s")
    print(f"  f0                   = {f0/1e9:.3f} GHz")
    print(f"  eigenvector  A/B     = {g['A']:.4f}")
    print(f"  particle mass M      = {g['M']:.4e} kg")

    banner("2. Geometric factor (full angular calculation)")
    print(f"  alpha (Frobenius)    = {g['alpha_frob']:.4f}")
    print(f"  alpha (zz-component) = {g['alpha_zz']:.4f}")
    print("  => order unity, computed from the eigenfunction (not assumed).")

    banner("3. Matrix-element ratio  |Q_body|^2 / |Q_COM|^2")
    for wcom in (1e5, 1e6):
        rm = cp.matrix_ratio(lm.SILICA, R, wcom,
                             alpha=g["alpha_frob"], w0=g["omega0"], M=g["M"])
        print(f"  omega_COM = {wcom:.0e} rad/s : R_mat = {rm:.3e}  (~1e{np.log10(rm):.1f})")

    banner("4. Spectral ratio  S_EE(omega0)/S_EE(omega_COM) = (omega0/omega_COM)^n")
    for n in (1, 5):
        sr = cp.spectral_ratio(g["omega0"], 1e5, n)
        print(f"  n = {n} (omega_COM=1e5): {sr:.3e}  (~1e{np.log10(sr):.1f})")

    banner("5. Combined advantage  Gamma_internal / Gamma_COM = R_mat x (omega0/omega_COM)^n")
    for wcom in (3e4, 1e5, 1e6):
        comb = cp.combined_advantage(lm.SILICA, R, wcom, 5)
        print(f"  n = 5, omega_COM = {wcom:.0e} rad/s : {comb:.3e}  (~1e{np.log10(comb):.1f})")
    print("  Combined advantage ~ omega_COM^-3; ~1e32 at 1e5 rad/s, up to ~1e33 at low trap freq.")

    banner("6. Rate scalings (parameter-free exponents)")
    print("  |Q_COM|^2  ~ omega^-2,  S_EE ~ omega^5  =>  Gamma_COM      ~ omega^3   (recovers Toros et al.)")
    print("  |Q_body|^2 ~ omega^-1,  S_EE ~ omega^5  =>  Gamma_internal ~ omega^4")

    banner("7. Material spread (l=2 fundamental)")
    M_silica = lm.SILICA.rho * (4 / 3) * np.pi * R**3
    print(f"  {'Material':9}{'nu':>7}{'eta':>8}{'f0@100nm':>11}{'f0@matchedM':>13}{'alpha':>8}")
    for mat in lm.MATERIALS:
        gm = cp.geometric_factor(mat, R)
        Rm = (3 * M_silica / (4 * np.pi * mat.rho))**(1 / 3)
        w0m = cp.omega0(mat, Rm)
        print(f"  {mat.name:9}{mat.poisson:7.2f}{gm['eta']:8.3f}"
              f"{gm['omega0']/TWO_PI/1e9:9.1f}G{w0m/TWO_PI/1e9:11.1f}G{gm['alpha_frob']:8.3f}")
    print("  (Crystals are anisotropic; isotropic-equivalent constants -> order-of-magnitude.)")


if __name__ == "__main__":
    main()
