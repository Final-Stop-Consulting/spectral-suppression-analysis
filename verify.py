#!/usr/bin/env python3
"""
verify.py  --  master verification suite
========================================
Independent checks of the whole calculation: elastic eigenmodes, the tidal
coupling and channel ratio, the first finite bound and its no-go, the material
discriminant with background orthogonality, and the cutoff / readout thresholds.

    python verify.py     # prints PASS/FAIL; exits nonzero on any failure
"""
import sys
import numpy as np
from scipy.integrate import quad
import lamb_mode as lm
import coupling as cp
import sensitivity as sn

R, T, TAU = 100e-9, 0.1, 1e5
res = []
def ck(name, cond, detail=""):
    res.append(bool(cond))
    print(f"  [{'PASS' if cond else 'FAIL'}] {name}" + (f"  ({detail})" if detail else ""))

print("Eigenmodes and coupling")
g = cp.geometric_factor(lm.SILICA, R); w0 = g["omega0"]
h, k = w0/lm.SILICA.v_l, w0/lm.SILICA.v_t
ck("eta = 2.628 (silica l=2 fundamental)", abs(g["eta"]-2.6282) < 1e-3, f"eta={g['eta']:.4f}")
# analytic vs finite-difference derivatives
r0 = 0.37*R; dr = r0*1e-6
Up_a, Vp_a = lm.dUV(r0, g["A"], g["B"], h, k)
Up_fd = (lm.UV(r0+dr,g["A"],g["B"],h,k)[0]-lm.UV(r0-dr,g["A"],g["B"],h,k)[0])/(2*dr)
Vp_fd = (lm.UV(r0+dr,g["A"],g["B"],h,k)[1]-lm.UV(r0-dr,g["A"],g["B"],h,k)[1])/(2*dr)
ck("analytic U',V' match finite difference", abs(Up_a-Up_fd)/abs(Up_fd)<1e-5 and abs(Vp_a-Vp_fd)/abs(Vp_fd)<1e-5)
s_rr, s_rt = lm.tractions(R, g["A"], g["B"], h, k, lm.SILICA)
scale = lm.SILICA.lame_mu*abs(lm.UV(R,g["A"],g["B"],h,k)[0])/R
ck("both traction-free BCs vanish at eigenfrequency", abs(s_rr)/scale<1e-6 and abs(s_rt)/scale<1e-6)
A_rt = -lm.tractions(R,0,1,h,k,lm.SILICA)[1]/lm.tractions(R,1,0,h,k,lm.SILICA)[1]
ck("eigenvector consistent across both BCs", abs(A_rt-g["A"])<1e-5)
def Y(t): return np.sqrt(5/(16*np.pi))*(3*np.cos(t)**2-1)
def dY(t): return np.sqrt(5/(16*np.pi))*(-6*np.cos(t)*np.sin(t))
a1=quad(lambda t:Y(t)*np.cos(t)**2*np.sin(t),0,np.pi)[0]; a2=quad(lambda t:dY(t)*np.cos(t)*np.sin(t)**2,0,np.pi)[0]
b1=quad(lambda t:Y(t)*np.sin(t)**3,0,np.pi)[0]; b2=quad(lambda t:dY(t)*np.sin(t)**2*np.cos(t),0,np.pi)[0]
ck("induced quadrupole trace-free (b1=-a1, b2=a2)", abs(b1+a1)<1e-9 and abs(b2-a2)<1e-9)
ck("alpha_frob = 0.84 (order unity)", abs(g["alpha_frob"]-0.8367)<1e-3, f"alpha={g['alpha_frob']:.4f}")
ck("long-wavelength valid (lambda_grav/R > 1e4)", (sn.C/(w0/2/np.pi))/R>1e4)
ck("spectral suppression ~1e30 at n=5", 5e29 < cp.spectral_ratio(w0,1e5,5) < 2e30)

print("\nChannel ratio (reconciled)")
ck("COM two-phonon factor D = 10 (=2/3*<r^4>)", abs(cp.com_two_phonon_factor()-(2/3)*15)<1e-9)
ck("channel ratio single ~93", 85<cp.channel_ratio(lm.SILICA,R,1e5,"single")<100, f"{cp.channel_ratio(lm.SILICA,R,1e5,'single'):.1f}")
ck("channel ratio full ~62", 55<cp.channel_ratio(lm.SILICA,R,1e5,"full")<70, f"{cp.channel_ratio(lm.SILICA,R,1e5,'full'):.1f}")
ck("combined ~ omega_COM^-3 (Toros Gamma_COM~omega^3)", abs(cp.combined_advantage(lm.SILICA,R,1e5,5)/cp.combined_advantage(lm.SILICA,R,1e6,5)-1000)<1e-3)

print("\nThe bound, no-go, and vacuum reference")
ck("bound S_EE^min(Q=1e5) ~ 2.3e7 s^-3", 2.0e7<sn.bound(Q=1e5)<2.7e7, f"{sn.bound(Q=1e5):.2e}")
ck("bound ~ Q^-1/2 : S(1e4)/S(1e5) ~ sqrt(10)", abs(sn.bound(Q=1e4)/sn.bound(Q=1e5)-np.sqrt(10))<0.05)
ck("vacuum S_EE^vac(16GHz) ~ 2.7e-32 s^-3", 2.5e-32<sn.vacuum_reference(w0)<3.0e-32, f"{sn.vacuum_reference(w0):.2e}")
ck("gap bound/vacuum ~ 1e39", 5e38<sn.gap_to_vacuum()<2e39, f"{sn.gap_to_vacuum():.2e}")
ck("no-go: gap still >>1 at R=30 cm (bar scale)", sn.gap_to_vacuum(R=0.3)>1e30, f"{sn.gap_to_vacuum(R=0.3):.1e}")

print("\nMaterial discriminant")
tid5=[sn.tidal_rate(m,R,5) for m in lm.MATERIALS]; th=[sn.thermal_loss(m,R) for m in lm.MATERIALS]
ck("tidal n=5 diamond/silica > 100 (rises)", tid5[3]/tid5[0]>100, f"{tid5[3]/tid5[0]:.0f}")
ck("thermal loss diamond/silica < 1e-6 (falls)", th[3]/th[0]<1e-6, f"{th[3]/th[0]:.1e}")
def pear(x,y):
    x=np.log(x)-np.mean(np.log(x)); y=np.log(y)-np.mean(np.log(y)); return (x@y)/np.sqrt((x@x)*(y@y))
ck("log-Pearson(tidal n=5, thermal) < -0.9 (orthogonal)", pear(tid5,th)<-0.9, f"{pear(tid5,th):.3f}")

print("\nCutoff and readout thresholds")
wc5=sn.cutoff_threshold(lm.SILICA,lm.DIAMOND,R,5)
ck("cutoff omega_c*(n=5) ~ 4.2e10 rad/s", 4.0e10<wc5<4.5e10, f"{wc5:.2e}")
ck("Earth c/r_s just below omega_c*(n=5)", 3e10<sn.C/sn.R_S_EARTH<wc5, f"{sn.C/sn.R_S_EARTH:.2e}")
ck("epsilon_threshold = 1/sqrt(1400) ~ 2.7%", abs(sn.epsilon_threshold(1400)-0.0267)<0.001)

nf = res.count(False)
print(f"\n{len(res)-nf}/{len(res)} checks passed.")
sys.exit(1 if nf else 0)
