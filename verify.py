#!/usr/bin/env python3
"""
verify.py
=========
Independent self-checks of the spectral-suppression calculation. Each test
either confirms a numerical result against an independent computation or checks
an analytic identity. Run:

    python verify.py        # prints PASS/FAIL, exits nonzero on any failure
"""
import sys
import numpy as np
from scipy.integrate import quad

import lamb_mode as lm
import coupling as cp

R = 100e-9
C_LIGHT = 2.99792458e8
results = []


def check(name, condition, detail=""):
    results.append((name, bool(condition), detail))
    print(f"  [{'PASS' if condition else 'FAIL'}] {name}" + (f"  ({detail})" if detail else ""))


print("Verification of the l=2 Lamb tidal-coupling calculation\n")

# 1. Fundamental eigenvalue reproduces the expected Lamb value.
g = cp.geometric_factor(lm.SILICA, R)
check("eta = 2.628 (silica fundamental l=2)",
      abs(g["eta"] - 2.6282) < 1e-3, f"eta={g['eta']:.4f}")

# 2. Analytic Bessel-derivative formula agrees with finite differences.
h = g["omega0"] / lm.SILICA.v_l
k = g["omega0"] / lm.SILICA.v_t
r0 = 0.37 * R
Up_a, Vp_a = lm.dUV(r0, g["A"], g["B"], h, k)
dr = r0 * 1e-6
Up_fd = (lm.UV(r0 + dr, g["A"], g["B"], h, k)[0] - lm.UV(r0 - dr, g["A"], g["B"], h, k)[0]) / (2 * dr)
Vp_fd = (lm.UV(r0 + dr, g["A"], g["B"], h, k)[1] - lm.UV(r0 - dr, g["A"], g["B"], h, k)[1]) / (2 * dr)
check("analytic U',V' match finite difference",
      abs(Up_a - Up_fd) / abs(Up_fd) < 1e-5 and abs(Vp_a - Vp_fd) / abs(Vp_fd) < 1e-5)

# 3. Both traction-free boundary conditions vanish at the eigenfrequency.
s_rr, s_rt = lm.tractions(R, g["A"], g["B"], h, k, lm.SILICA)
scale = lm.SILICA.lame_mu * abs(lm.UV(R, g["A"], g["B"], h, k)[0]) / R
check("sigma_rr(R) = 0 and sigma_r_theta(R) = 0",
      abs(s_rr) / scale < 1e-6 and abs(s_rt) / scale < 1e-6,
      f"|srr|/s={abs(s_rr)/scale:.1e}, |srt|/s={abs(s_rt)/scale:.1e}")

# 4. Eigenvector from sigma_rr=0 also satisfies sigma_r_theta=0 (genuine root).
A_rt = -lm.tractions(R, 0, 1, h, k, lm.SILICA)[1] / lm.tractions(R, 1, 0, h, k, lm.SILICA)[1]
check("eigenvector consistent across both BCs",
      abs(A_rt - g["A"]) < 1e-5, f"dA={abs(A_rt-g['A']):.1e}")

# 5. Induced quadrupole is trace-free with Q_xx=Q_yy=-Q_zz/2.
#    Equivalent to the angular identities  b1=-a1  and  b2=a2.
def Y(t):  return np.sqrt(5/(16*np.pi))*(3*np.cos(t)**2-1)
def dY(t): return np.sqrt(5/(16*np.pi))*(-6*np.cos(t)*np.sin(t))
a1 = quad(lambda t: Y(t)*np.cos(t)**2*np.sin(t), 0, np.pi)[0]
a2 = quad(lambda t: dY(t)*np.cos(t)*np.sin(t)**2, 0, np.pi)[0]
b1 = quad(lambda t: Y(t)*np.sin(t)**3, 0, np.pi)[0]
b2 = quad(lambda t: dY(t)*np.sin(t)**2*np.cos(t), 0, np.pi)[0]
check("quadrupole trace-free (b1=-a1, b2=a2)",
      abs(b1 + a1) < 1e-9 and abs(b2 - a2) < 1e-9)

# 6. alpha reproduces the reported value and is order unity.
check("alpha_frob = 0.84 (order unity)",
      abs(g["alpha_frob"] - 0.8367) < 1e-3 and 0.5 < g["alpha_frob"] < 1.5,
      f"alpha={g['alpha_frob']:.4f}")

# 7. Matrix ratio scales as omega_COM^2 (two-phonon zero-point amplitude).
r1 = cp.matrix_ratio(lm.SILICA, R, 1e5, alpha=g["alpha_frob"], w0=g["omega0"], M=g["M"])
r2 = cp.matrix_ratio(lm.SILICA, R, 1e6, alpha=g["alpha_frob"], w0=g["omega0"], M=g["M"])
check("R_mat ~ omega_COM^2", abs(r2 / r1 - 100) < 1e-6)

# 8. Combined advantage scales as omega_COM^-3 (=> Gamma_COM ~ omega^3, Toros).
c1 = cp.combined_advantage(lm.SILICA, R, 1e5, 5)
c2 = cp.combined_advantage(lm.SILICA, R, 1e6, 5)
check("combined ~ omega_COM^-3 (Toros Gamma_COM ~ omega^3)",
      abs(c1 / c2 - 1000) < 1e-3, f"ratio={c1/c2:.1f}")

# 9. Spectral suppression at n=5 is ~1e30.
sr = cp.spectral_ratio(g["omega0"], 1e5, 5)
check("spectral suppression ~ 1e30 at n=5", 5e29 < sr < 2e30, f"{sr:.2e}")

# 10. Long-wavelength (uniform-field) coupling is valid: lambda_grav >> R.
lam_grav = C_LIGHT / (g["omega0"] / (2*np.pi))
check("long-wavelength expansion valid (lambda_grav/R > 1e4)",
      lam_grav / R > 1e4, f"lambda/R={lam_grav/R:.1e}")

n_fail = sum(1 for _, ok, _ in results if not ok)
print(f"\n{len(results)-n_fail}/{len(results)} checks passed.")
sys.exit(1 if n_fail else 0)
