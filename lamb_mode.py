"""
lamb_mode.py
============
Spheroidal ("Lamb") elastic eigenmodes of a free, homogeneous, isotropic sphere,
and the mass-quadrupole coupling of the fundamental l=2 mode to a gravitational
tidal field.

This module underlies the calculation in the accompanying manuscript,
"Spectral Suppression of Center-of-Mass Coupling to Gravitational Tidal
Fluctuations" (2026).

Physics
-------
The spheroidal displacement field of degree l is built from two scalar
potentials (Lamb 1882):

    phi = A j_l(h r) Y_lm   (compressional,  h = omega / v_l)
    psi = B j_l(k r) Y_lm   (shear-poloidal, k = omega / v_t)

giving radial functions

    U(r) = A h j_l'(h r) + B l(l+1) j_l(k r) / r
    V(r) = A j_l(h r)/r  + B [ j_l(k r)/r + k j_l'(k r) ]

with displacement  u = U(r) Y_lm r_hat + V(r) grad_1 Y_lm.

The traction-free surface conditions sigma_rr(R) = sigma_r_theta(R) = 0 give a
2x2 homogeneous system in (A, B); the vanishing of its determinant is the
secular equation that quantizes omega. Derivatives of the spherical Bessel
functions are evaluated analytically (the spherical Bessel ODE supplies j'').

References
----------
H. Lamb, Proc. London Math. Soc. s1-13, 189 (1882).
A. C. Eringen and E. S. Suhubi, *Elastodynamics* Vol. 2 (Academic, 1975).
L. Saviot and D. B. Murray, Phys. Rev. B 79, 214101 (2009).
"""
from __future__ import annotations
from dataclasses import dataclass
import numpy as np
from scipy.special import spherical_jn


# --------------------------------------------------------------------------- #
#  Spherical Bessel function and its first two derivatives (analytic)
# --------------------------------------------------------------------------- #
def j(n: int, x):
    """Spherical Bessel function j_n(x)."""
    return spherical_jn(n, x)


def jp(n: int, x):
    """First derivative j_n'(x)."""
    return spherical_jn(n, x, derivative=True)


def jpp(n: int, x):
    """Second derivative j_n''(x), from the spherical Bessel ODE
        x^2 j'' + 2x j' + (x^2 - n(n+1)) j = 0.
    """
    return -(2.0 / x) * jp(n, x) + (n * (n + 1) / x**2 - 1.0) * j(n, x)


# --------------------------------------------------------------------------- #
#  Material
# --------------------------------------------------------------------------- #
@dataclass(frozen=True)
class Material:
    """Isotropic elastic material.  Speeds in m/s, density in kg/m^3."""
    name: str
    v_l: float          # longitudinal (compressional) sound speed
    v_t: float          # transverse (shear) sound speed
    rho: float          # mass density

    @property
    def poisson(self) -> float:
        vl, vt = self.v_l, self.v_t
        return (vl**2 - 2 * vt**2) / (2 * (vl**2 - vt**2))

    @property
    def lame_lambda(self) -> float:
        return self.rho * (self.v_l**2 - 2 * self.v_t**2)

    @property
    def lame_mu(self) -> float:
        return self.rho * self.v_t**2


# Isotropic-equivalent bulk constants (approximate for anisotropic crystals).
SILICA   = Material("Silica",   5970.0,  3760.0, 2200.0)
SILICON  = Material("Silicon",  8433.0,  5843.0, 2330.0)
SAPPHIRE = Material("Sapphire", 11100.0, 6040.0, 3980.0)
DIAMOND  = Material("Diamond",  18000.0, 12000.0, 3515.0)
MATERIALS = [SILICA, SILICON, SAPPHIRE, DIAMOND]


# --------------------------------------------------------------------------- #
#  Eigenfunction radial functions and surface tractions
# --------------------------------------------------------------------------- #
def UV(r, A, B, h, k, l=2):
    """Radial displacement functions U(r), V(r)."""
    U = A * h * jp(l, h * r) + B * l * (l + 1) * j(l, k * r) / r
    V = A * j(l, h * r) / r + B * (j(l, k * r) / r + k * jp(l, k * r))
    return U, V


def dUV(r, A, B, h, k, l=2):
    """Analytic radial derivatives U'(r), V'(r)."""
    Up = (A * h**2 * jpp(l, h * r)
          + B * l * (l + 1) * (k * jp(l, k * r) / r - j(l, k * r) / r**2))
    Vp = (A * (h * jp(l, h * r) / r - j(l, h * r) / r**2)
          + B * (k * jp(l, k * r) / r - j(l, k * r) / r**2 + k**2 * jpp(l, k * r)))
    return Up, Vp


def tractions(r, A, B, h, k, mat: Material, l=2):
    """Surface tractions (sigma_rr, sigma_r_theta) with the angular factor
    divided out:  sigma_rr / Y_lm  and  sigma_r_theta / (d_theta Y_lm)."""
    lam, mu = mat.lame_lambda, mat.lame_mu
    U, V = UV(r, A, B, h, k, l)
    Up, Vp = dUV(r, A, B, h, k, l)
    s_rr = lam * (Up + 2 * U / r - l * (l + 1) * V / r) + 2 * mu * Up
    s_rt = mu * (U / r + Vp - V / r)
    return s_rr, s_rt


def secular_determinant(xi, mat: Material, R, l=2):
    """Determinant of the 2x2 traction-free system as a function of the
    reduced frequency xi = omega R / v_t."""
    omega = xi * mat.v_t / R
    h, k = omega / mat.v_l, omega / mat.v_t
    s_rr_A, s_rt_A = tractions(R, 1.0, 0.0, h, k, mat, l)
    s_rr_B, s_rt_B = tractions(R, 0.0, 1.0, h, k, mat, l)
    return s_rr_A * s_rt_B - s_rr_B * s_rt_A


def fundamental_eta(mat: Material, R, l=2, xi_lo=0.05, xi_hi=4.0, n_scan=4000):
    """Lowest root xi = eta of the secular equation (the fundamental l mode),
    by sign-change scan + bisection."""
    xs = np.linspace(xi_lo, xi_hi, n_scan)
    ds = np.array([secular_determinant(x, mat, R, l) for x in xs])
    for i in range(len(xs) - 1):
        if ds[i] * ds[i + 1] < 0:
            a, b = xs[i], xs[i + 1]
            for _ in range(200):
                m = 0.5 * (a + b)
                if (secular_determinant(a, mat, R, l)
                        * secular_determinant(m, mat, R, l)) <= 0:
                    b = m
                else:
                    a = m
            return 0.5 * (a + b)
    raise RuntimeError("No fundamental root found in scan window.")


def eigenvector(eta, mat: Material, R, l=2):
    """Amplitude ratio (A, B) at the eigenfrequency, from sigma_rr = 0.
    Returns (A, B) with B = 1.  Both boundary conditions agree at a true root
    (checked in the test suite)."""
    omega = eta * mat.v_t / R
    h, k = omega / mat.v_l, omega / mat.v_t
    s_rr_A, _ = tractions(R, 1.0, 0.0, h, k, mat, l)
    s_rr_B, _ = tractions(R, 0.0, 1.0, h, k, mat, l)
    A = -s_rr_B / s_rr_A
    return A, 1.0
