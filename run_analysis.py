#!/usr/bin/env python3
"""
run_analysis.py  --  master driver
==================================
Reproduces every number and table supporting the manuscript:
elastic eigenmodes, the tidal coupling and channel ratio, the first finite
bound and its no-go, the material discriminant with intrinsic backgrounds, and
the cutoff / readout thresholds.

    python run_analysis.py
"""
import numpy as np
import lamb_mode as lm
import coupling as cp
import sensitivity as sn

R, T, TAU = 100e-9, 0.1, 1e5
MATS = [lm.SILICA, lm.SILICON, lm.SAPPHIRE, lm.DIAMOND]
TWO_PI = 2 * np.pi


def h(s): print("\n" + s + "\n" + "-" * len(s))


def main():
    h("1. Elastic eigenmodes (l=2 fundamental, R = 100 nm)")
    print(f"  {'Material':9}{'nu':>6}{'eta':>8}{'f0 (GHz)':>10}{'alpha':>8}")
    G = {}
    for m in MATS:
        g = cp.geometric_factor(m, R); G[m.name] = g
        print(f"  {m.name:9}{m.poisson:6.2f}{g['eta']:8.3f}{g['omega0']/TWO_PI/1e9:10.2f}{g['alpha_frob']:8.3f}")

    g = G["Silica"]; w0 = g["omega0"]
    h("2. Tidal coupling (silica)")
    print(f"  geometric factor alpha = {g['alpha_frob']:.3f}   (f0 = {w0/TWO_PI/1e9:.2f} GHz)")
    print(f"  COM two-phonon factor  D = {cp.com_two_phonon_factor():.0f}  (= 2/3 * <r^4>)")
    print("  channel ratio |Q_body|^2/|Q_COM|^2 at omega_COM = 1e5 rad/s:")
    for conv in ("single", "full", "naive"):
        print(f"     {conv:7}: {cp.channel_ratio(lm.SILICA, R, 1e5, conv):6.1f}")
    print("  spectral ratio (omega0/omega_COM)^n:  n=1 -> %.2e   n=5 -> %.2e"
          % (cp.spectral_ratio(w0, 1e5, 1), cp.spectral_ratio(w0, 1e5, 5)))
    print("  combined advantage (single, n=5) = %.2e" % cp.combined_advantage(lm.SILICA, R, 1e5, 5))
    print("  rate scalings: Gamma_COM ~ omega^3 (recovers Toros); Gamma_internal ~ omega^4")

    h("3. First finite bound on S_EE at 16 GHz (silica)")
    for Q in (1e4, 1e5, 1e6):
        print(f"  Q = {Q:.0e}:  S_EE^min = {sn.bound(Q=Q):.2e} s^-3")
    print(f"  vacuum reference S_EE^vac = t_Pl^2 omega0^5 = {sn.vacuum_reference(w0):.2e} s^-3")
    print(f"  gap: bound/vacuum = {sn.gap_to_vacuum():.2e}  (~1e39; bounds enhanced states only)")

    h("4. No-go: reach/vacuum ratio ~ 1/R (never closes)")
    print(f"  {'R':>10}{'f0 (Hz)':>13}{'S_EE^min':>12}{'gap/vac':>12}")
    for Rr in (1e-7, 1e-3, 0.1, 0.3):
        w = cp.omega0(lm.SILICA, Rr)
        print(f"  {Rr*1e3:8.3f}mm{w/TWO_PI:13.2e}{sn.bound(R=Rr):12.2e}{sn.gap_to_vacuum(R=Rr):12.2e}")
    print("  Even at bar-detector scale (R = 30 cm, f0 ~ 5 kHz) the gap is ~1e33: it never closes.")

    h("5. Material discriminant (fixed R; normalized to silica)")
    print(f"  {'channel':16}" + "".join(f"{m.name:>10}" for m in MATS))
    for n in (1, 3, 5):
        row = [sn.tidal_rate(m, R, n) / sn.tidal_rate(lm.SILICA, R, n) for m in MATS]
        print(f"  tidal n={n:<9}" + "".join(f"{x:10.2f}" for x in row))
    row = [sn.thermal_loss(m, R) / sn.thermal_loss(lm.SILICA, R) for m in MATS]
    print(f"  thermal loss    " + "".join(f"{x:10.2e}" for x in row))
    row = [sn.gas_rate(m, R) / sn.gas_rate(lm.SILICA, R) for m in MATS]
    print(f"  gas             " + "".join(f"{x:10.2f}" for x in row))
    def pear(a, b):
        a = np.log(a) - np.mean(np.log(a)); b = np.log(b) - np.mean(np.log(b))
        return (a @ b) / np.sqrt((a @ a) * (b @ b))
    tid5 = [sn.tidal_rate(m, R, 5) for m in MATS]
    th = [sn.thermal_loss(m, R) for m in MATS]; gs = [sn.gas_rate(m, R) for m in MATS]
    print(f"  log-Pearson(tidal n=5, thermal loss) = {pear(tid5, th):+.3f}")
    print(f"  log-Pearson(tidal n=5, gas)          = {pear(tid5, gs):+.3f}")
    print("  Tidal is flat-or-rising (n>=1); backgrounds fall. Sign discriminant rests on n>=1.")

    h("6. Cutoff threshold (silica = diamond sign inversion)")
    for n in (1, 3, 5, 6):
        print(f"  n = {n}:  omega_c* = {sn.cutoff_threshold(lm.SILICA, lm.DIAMOND, R, n):.2e} rad/s")
    print(f"  Earth c/r_s (r_s = 8.87 mm) = {sn.C/sn.R_S_EARTH:.2e} rad/s  (just below omega_c*(n=5))")

    h("7. Readout epsilon-threshold")
    print(f"  cooperativity C = 1400 (Mayor et al. 2025) -> epsilon >= 1/sqrt(C) = {sn.epsilon_threshold(1400)*100:.1f}%")


if __name__ == "__main__":
    main()
