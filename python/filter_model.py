"""Independent re-simulation of the ADS filter designs with ABCD matrices.

Component models match the ADS ones used in the project:
* INDQ / CAPQ in "proportional to freq" mode: Q(f) = Q0 * f / f0, which
  gives a frequency-independent loss resistance/conductance
      R_L = 2*pi*f0*L / Q0      (series with the inductor)
      G_C = 2*pi*f0*C / Q0      (in parallel with the capacitor)
* Ideal parts are obtained with Q0 = inf.

Values are in nH / pF; all functions are vectorised over frequency (Hz).
"""

import numpy as np

Z0 = 50.0
QL, QC, F_Q = 100.0, 1000.0, 100e6           # project settings: QL=100, QC=1000 at 100 MHz


# ----------------------------------------------------------------- elements
def z_ind(f, L_nH, q=QL):
    L = L_nH * 1e-9
    return 1j * 2 * np.pi * f * L + (2 * np.pi * F_Q * L / q if np.isfinite(q) else 0.0)


def y_cap(f, C_pF, q=QC):
    C = C_pF * 1e-12
    return 1j * 2 * np.pi * f * C + (2 * np.pi * F_Q * C / q if np.isfinite(q) else 0.0)


def series(z):
    one, zero = np.ones_like(z), np.zeros_like(z)
    return np.moveaxis(np.array([[one, z], [zero, one]]), -1, 0)


def shunt(y):
    one, zero = np.ones_like(y), np.zeros_like(y)
    return np.moveaxis(np.array([[one, zero], [y, one]]), -1, 0)


def cascade(*stages):
    out = stages[0]
    for s in stages[1:]:
        out = out @ s
    return out


def s_params(abcd, z0=Z0):
    """S11 and S21 of a reciprocal two-port between equal terminations."""
    a, b, c, d = abcd[:, 0, 0], abcd[:, 0, 1], abcd[:, 1, 0], abcd[:, 1, 1]
    den = a + b / z0 + c * z0 + d
    return (a + b / z0 - c * z0 - d) / den, 2.0 / den


def tank(f, L, C, ql, qc):
    """Series arm: inductor in parallel with a capacitor (impedance)."""
    return 1.0 / (1.0 / z_ind(f, L, ql) + y_cap(f, C, qc))


# ------------------------------------------------------------------ designs
# 1) 7th-order elliptic low-pass from the ADS Filter DesignGuide
#    (Ap = 0.11 dB, As = 40 dB, Fp = 161 MHz, Fs = 188 MHz, 50 ohm, shunt-first)
LP_PROTOTYPE = dict(shunt_C=[20.009399, 25.637365, 20.884861, 12.042828],
                    arm_L=[57.670514, 29.72434, 33.868581],
                    arm_C=[4.789392, 26.560229, 18.853142])

# 2) Low-pass after adding shunt transmission zeros (L4, L5) and optimisation
#    with finite-Q parts. Shunt branches: C4+L4 (input), C5, C6, C7+L5 (output).
LP_OPTIMISED = dict(L=[27.4141, 9.95827, 25.1808, 24.5999, 16.0994],
                    C=[20.2239, 3.68688, 27.6634, 14.2437, 3.60615, 28.5257, 12.1215])

# 3) High-pass after optimisation. Series arms L1||C1, L2||C2, L3||C3;
#    shunt inductors L4+L8 (input), L5, L6, L7+L9 (output).
HP_OPTIMISED = dict(L=[182.657, 47.6741, 60.3074, 50.201, 62.4427, 64.8717, 82.1298, 44.1952, 46.84],
                    C=[24.3573, 45.4245, 42.8266])


def lowpass_prototype(f, d=LP_PROTOTYPE, ql=np.inf, qc=np.inf):
    sc, al, ac = d["shunt_C"], d["arm_L"], d["arm_C"]
    return cascade(shunt(y_cap(f, sc[0], qc)), series(tank(f, al[0], ac[0], ql, qc)),
                   shunt(y_cap(f, sc[1], qc)), series(tank(f, al[1], ac[1], ql, qc)),
                   shunt(y_cap(f, sc[2], qc)), series(tank(f, al[2], ac[2], ql, qc)),
                   shunt(y_cap(f, sc[3], qc)))


def lowpass_optimised(f, d=LP_OPTIMISED, ql=QL, qc=QC):
    L, C = d["L"], d["C"]

    def shunt_lc(Lx, Cx):                      # series L-C to ground -> transmission zero
        return shunt(1.0 / (z_ind(f, Lx, ql) + 1.0 / y_cap(f, Cx, qc)))

    return cascade(shunt_lc(L[3], C[3]), series(tank(f, L[0], C[0], ql, qc)),
                   shunt(y_cap(f, C[4], qc)), series(tank(f, L[1], C[1], ql, qc)),
                   shunt(y_cap(f, C[5], qc)), series(tank(f, L[2], C[2], ql, qc)),
                   shunt_lc(L[4], C[6]))


def highpass_optimised(f, d=HP_OPTIMISED, ql=QL, qc=QC):
    L, C = d["L"], d["C"]
    return cascade(shunt(1.0 / (z_ind(f, L[3], ql) + z_ind(f, L[7], ql))),
                   series(tank(f, L[0], C[0], ql, qc)),
                   shunt(1.0 / z_ind(f, L[4], ql)),
                   series(tank(f, L[1], C[1], ql, qc)),
                   shunt(1.0 / z_ind(f, L[5], ql)),
                   series(tank(f, L[2], C[2], ql, qc)),
                   shunt(1.0 / (z_ind(f, L[6], ql) + z_ind(f, L[8], ql))))


def bandpass(f, lp=LP_OPTIMISED, hp=HP_OPTIMISED, ql=QL, qc=QC):
    """Band-pass = low-pass section followed by the high-pass section (as in ADS)."""
    return lowpass_optimised(f, lp, ql, qc) @ highpass_optimised(f, hp, ql, qc)


def db(x):
    return 20.0 * np.log10(np.abs(x))


def group_delay(f, s21):
    phase = np.unwrap(np.angle(s21))
    return -np.gradient(phase, 2 * np.pi * f)
