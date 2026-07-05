"""
sensitivity.py
==============
Observable-level results: the first finite bound on the internal-band tidal
noise, the R^-6 no-go against the vacuum prediction, the intrinsic backgrounds,
the four-material discriminant, and the cutoff / readout thresholds.
"""
from __future__ import annotations
import numpy as np
import lamb_mode as lm
import coupling as cp

HBAR = 1.054571817e-34
KB   = 1.380649e-23
C    = 2.99792458e8
G    = 6.67430e-11
T_PL = np.sqrt(HBAR * G / C**5)          # Planck time
R_S_EARTH = 8.87e-3                       # Earth Schwarzschild radius (m)


def thermal_occupation(w, T):
    """Bose occupation n_th of a mode at angular frequency w, temperature T."""
    return 1.0 / np.expm1(HBAR * w / (KB * T))


def bound(mat=lm.SILICA, R=100e-9, T=0.1, tau=1e5, Q=1e5):
    """First finite bound on S_EE at the internal-mode frequency.
    Gamma_min = sqrt(gamma n_th / tau), gamma = omega0/Q (thermal-limited,
    no clamping loss); S_EE^min = Gamma_min hbar^2 / |Q_body|^2. Units s^-3."""
    g = cp.geometric_factor(mat, R)
    w0 = g["omega0"]
    Qbody2 = g["alpha_frob"]**2 * g["M"] * R**2 * HBAR / (2 * w0)
    gamma = w0 / Q
    Gamma_min = np.sqrt(gamma * thermal_occupation(w0, T) / tau)
    return Gamma_min * HBAR**2 / Qbody2


def vacuum_reference(w):
    """Linearized-gravity vacuum tidal spectral density S_EE^vac = t_Pl^2 w^5 (s^-3)."""
    return T_PL**2 * w**5


def gap_to_vacuum(mat=lm.SILICA, R=100e-9, **kw):
    """Ratio of the achievable bound to the vacuum prediction."""
    w0 = cp.omega0(mat, R)
    return bound(mat, R, **kw) / vacuum_reference(w0)


# --- intrinsic backgrounds and tidal heating (fixed R, per-material trends) ---
def tidal_rate(mat, R, n):
    """Gamma_tidal ~ rho omega0^(n-1) (fixed R). Relative scaling only."""
    return mat.rho * cp.omega0(mat, R)**(n - 1)


def thermal_loss(mat, R, T=0.1, Q=1e5):
    """Intrinsic thermal-loss heating ~ (omega0/Q) n_th (relative)."""
    w0 = cp.omega0(mat, R)
    return (w0 / Q) * thermal_occupation(w0, T)


def gas_rate(mat, R):
    """Sudden-impulse residual-gas heating of the l=2 mode ~ 1/(rho omega0) (relative)."""
    return 1.0 / (mat.rho * cp.omega0(mat, R))


def cutoff_threshold(matA, matB, R, n):
    """Sign-inversion cutoff scale omega_c* where Gamma_tidal(A)=Gamma_tidal(B)
    for S_EE ~ omega^n exp(-omega/omega_c). Solves
    omega_c* = (wB-wA)/ln[(rhoB/rhoA)(wB/wA)^(n-1)]."""
    wA, wB = cp.omega0(matA, R), cp.omega0(matB, R)
    return (wB - wA) / np.log((matB.rho / matA.rho) * (wB / wA)**(n - 1))


def epsilon_threshold(C):
    """Readout overlap threshold below which the mode is no longer single-phonon
    resolved: epsilon >= 1/sqrt(C) for cooperativity C.

    C is the demonstrated *classical* optomechanical cooperativity (Mayor et al.
    2025: C ~= 1400 at a 3 K operating point). The same device's quantum
    cooperativity there is C/n_th ~= 175 (threshold ~= 8%); at 100 mK, where
    n_th < 1, classical and quantum cooperativity coincide."""
    return 1.0 / np.sqrt(C)
