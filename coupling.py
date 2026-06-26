"""
coupling.py
===========
Gravitational tidal coupling of the fundamental l=2 Lamb mode, and the
spectral-suppression / matrix-element comparison against center-of-mass motion.

Defines the dimensionless geometric factor alpha, the matrix-element ratio,
the spectral-density ratio, and the combined internal-vs-COM advantage that is
the central result of the paper.
"""
from __future__ import annotations
import numpy as np
from scipy.integrate import quad

from lamb_mode import Material, fundamental_eta, eigenvector, UV

HBAR = 1.054571817e-34          # J s

# Exact angular integrals for the m=0, l=2 mode (sympy-evaluated, closed form):
#   a1 = int Y20 cos^2(theta) sin(theta) dtheta = 2*sqrt(5)/(15*sqrt(pi))
#   a2 = int (dY20/dtheta) cos(theta) sin^2(theta) dtheta = -2*sqrt(5)/(5*sqrt(pi))
# The induced quadrupole is trace-free with Q_xx = Q_yy = -Q_zz/2 (proved in the
# test suite), so the Frobenius norm is sqrt(3/2)|Q_zz|.
_A1 = 2 * np.sqrt(5) / (15 * np.sqrt(np.pi))
_A2 = -2 * np.sqrt(5) / (5 * np.sqrt(np.pi))


def omega0(mat: Material, R, l=2) -> float:
    """Angular frequency of the fundamental l mode."""
    return fundamental_eta(mat, R, l) * mat.v_t / R


def geometric_factor(mat: Material, R, l=2):
    """Dimensionless geometric factor alpha (Frobenius convention) and the
    eigenfrequency.  The induced STF mass quadrupole of the mode, normalized so
    that the modal effective mass equals the particle mass M, is

        |Q_body| = alpha * M * R * x_zpf,    x_zpf = sqrt(hbar / (2 M omega0)).

    Returns dict(eta, omega0, alpha_frob, alpha_zz, A, B, M).
    """
    eta = fundamental_eta(mat, R, l)
    w0 = eta * mat.v_t / R
    h, k = w0 / mat.v_l, w0 / mat.v_t
    A, B = eigenvector(eta, mat, R, l)

    U = lambda r: UV(r, A, B, h, k, l)[0]
    V = lambda r: UV(r, A, B, h, k, l)[1]

    IU = quad(lambda r: U(r) * r**3, 0, R, limit=400)[0]
    IV = quad(lambda r: V(r) * r**3, 0, R, limit=400)[0]
    Nkin = mat.rho * quad(lambda r: (U(r)**2 + l * (l + 1) * V(r)**2) * r**2,
                          0, R, limit=400)[0]          # int rho |u|^2 d^3x

    Qzz = mat.rho * 2 * np.pi * (2 * _A1 * IU - 2 * _A2 * IV)   # induced Q_zz (raw)
    Qfrob = np.sqrt(1.5) * abs(Qzz)

    M = mat.rho * (4 / 3) * np.pi * R**3
    alpha_frob = Qfrob / (R * np.sqrt(Nkin * M))
    alpha_zz = abs(Qzz) / (R * np.sqrt(Nkin * M))
    return dict(eta=eta, omega0=w0, alpha_frob=alpha_frob, alpha_zz=alpha_zz,
                A=A, B=B, M=M)


def matrix_ratio(mat: Material, R, omega_com, l=2, alpha=None, w0=None, M=None):
    """Matrix-element-only ratio |Q_body|^2 / |Q_COM|^2 (Frobenius convention):

        R_mat = (3/2) alpha^2 M R^2 omega_com^2 / (hbar omega0).

    The COM quadrupole is the two-phonon <2|z^2|0> matrix element."""
    if alpha is None or w0 is None or M is None:
        g = geometric_factor(mat, R, l)
        alpha, w0, M = g["alpha_frob"], g["omega0"], g["M"]
    return 1.5 * alpha**2 * M * R**2 * omega_com**2 / (HBAR * w0)


def spectral_ratio(w0, omega_com, n):
    """S_EE(omega0) / S_EE(omega_com) = (omega0 / omega_com)^n for S_EE ~ omega^n."""
    return (w0 / omega_com)**n


def combined_advantage(mat: Material, R, omega_com, n, l=2):
    """Ratio of decoherence/heating rates, Gamma_internal / Gamma_COM
       = R_mat * (omega0/omega_com)^n.  Scales as omega_com^(-3) at n=5."""
    g = geometric_factor(mat, R, l)
    rmat = matrix_ratio(mat, R, omega_com, l,
                        alpha=g["alpha_frob"], w0=g["omega0"], M=g["M"])
    return rmat * spectral_ratio(g["omega0"], omega_com, n)
