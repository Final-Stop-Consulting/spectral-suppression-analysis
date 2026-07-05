"""
coupling.py
===========
Gravitational tidal coupling of the fundamental l=2 Lamb mode: the induced mass
quadrupole, the dimensionless geometric factor alpha, the internal-vs-COM
matrix-element (channel) ratio, the spectral-density ratio, and the combined
advantage.
"""
from __future__ import annotations
import numpy as np
from functools import lru_cache
from scipy.integrate import quad

from lamb_mode import Material, fundamental_eta, eigenvector, UV

HBAR = 1.054571817e-34          # J s

# Exact angular integrals for the m=0, l=2 mode (closed form):
_A1 = 2 * np.sqrt(5) / (15 * np.sqrt(np.pi))
_A2 = -2 * np.sqrt(5) / (5 * np.sqrt(np.pi))


@lru_cache(maxsize=None)
def omega0(mat: Material, R, l=2) -> float:
    """Angular frequency of the fundamental l mode."""
    return fundamental_eta(mat, R, l) * mat.v_t / R


@lru_cache(maxsize=None)
def geometric_factor(mat: Material, R, l=2):
    """Dimensionless geometric factor alpha (Frobenius convention) for the m=0
    mode, with the eigenfrequency and mass. Induced STF quadrupole, mode
    normalized to effective mass M: |Q_body| = alpha * M * R * x_zpf,
    x_zpf = sqrt(hbar/(2 M omega0)). Returns dict(eta, omega0, alpha_frob,
    alpha_zz, A, B, M)."""
    eta = fundamental_eta(mat, R, l)
    w0 = eta * mat.v_t / R
    h, k = w0 / mat.v_l, w0 / mat.v_t
    A, B = eigenvector(eta, mat, R, l)
    U = lambda r: UV(r, A, B, h, k, l)[0]
    V = lambda r: UV(r, A, B, h, k, l)[1]
    IU = quad(lambda r: U(r) * r**3, 0, R, limit=400)[0]
    IV = quad(lambda r: V(r) * r**3, 0, R, limit=400)[0]
    Nkin = mat.rho * quad(lambda r: (U(r)**2 + l*(l+1)*V(r)**2) * r**2, 0, R, limit=400)[0]
    Qzz = mat.rho * 2 * np.pi * (2 * _A1 * IU - 2 * _A2 * IV)
    Qfrob = np.sqrt(1.5) * abs(Qzz)
    M = mat.rho * (4 / 3) * np.pi * R**3
    alpha_frob = Qfrob / (R * np.sqrt(Nkin * M))
    alpha_zz = abs(Qzz) / (R * np.sqrt(Nkin * M))
    return dict(eta=eta, omega0=w0, alpha_frob=alpha_frob, alpha_zz=alpha_zz, A=A, B=B, M=M)


def com_two_phonon_factor() -> float:
    """The 3D isotropic center-of-mass two-phonon Frobenius factor
    D = sum_{ij, N=2} |<f|Q_ij/(M x_zpf^2)|0>|^2 = (2/3) <0|r^4|0> = (2/3)(15) = 10.
    (Verified by direct Fock-space summation.)"""
    return 10.0


# Coefficient c in  R_mat = c * alpha^2 M R^2 omega_com^2 /(hbar omega0),
# by mode-counting convention:
_RMAT_COEFF = {"single": 1.5,   # one m=0 internal mode vs 1D COM (D=4/3): ~93
               "full":   1.0,   # five l=2 m-modes vs 3D COM (D=10): ~62
               "naive":  2.0}   # single internal vs COM with D=1: ~124


def channel_ratio(mat: Material, R, omega_com, convention="single", l=2,
                  alpha=None, w0=None, M=None):
    """Matrix-element-only ratio |Q_body|^2/|Q_COM|^2 at equal noise PSD.
    convention: 'single' (appendix number, ~93), 'full' (total-to-total, ~62),
    or 'naive' (D=1, ~124). The advantage is O(10^2) and robust; the precise
    value is a mode-counting convention."""
    if alpha is None or w0 is None or M is None:
        g = geometric_factor(mat, R, l)
        alpha, w0, M = g["alpha_frob"], g["omega0"], g["M"]
    return _RMAT_COEFF[convention] * alpha**2 * M * R**2 * omega_com**2 / (HBAR * w0)


# Back-compat alias (single-mode convention = the appendix number).
def matrix_ratio(mat, R, omega_com, l=2, alpha=None, w0=None, M=None):
    return channel_ratio(mat, R, omega_com, "single", l, alpha, w0, M)


def spectral_ratio(w0, omega_com, n):
    """S_EE(omega0)/S_EE(omega_com) = (omega0/omega_com)^n for S_EE ~ omega^n."""
    return (w0 / omega_com)**n


def combined_advantage(mat: Material, R, omega_com, n, convention="single", l=2):
    """Ratio of decoherence/heating rates = channel_ratio * (omega0/omega_com)^n."""
    g = geometric_factor(mat, R, l)
    r = channel_ratio(mat, R, omega_com, convention, l,
                      alpha=g["alpha_frob"], w0=g["omega0"], M=g["M"])
    return r * spectral_ratio(g["omega0"], omega_com, n)
