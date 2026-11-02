"""Reproduce, check and extend the ADS filter results.

    python analyze.py

1. Validates the Python model against the marker values read in ADS.
2. Plots the low-pass, high-pass and band-pass responses (finite-Q parts).
3. Computes the band-pass group delay.
4. Monte-Carlo tolerance analysis: how much rejection margin survives
   realistic component tolerances (beyond the original assignment).
"""

import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

import filter_model as fm  # noqa: E402

DOCS = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "docs")
PASS = (119e6, 161e6)
F = np.linspace(1e6, 300e6, 3000)

# Marker values read from the ADS data displays (finite-Q, before real parts)
ADS_MARKERS = {
    "Low-pass (optimised)": (fm.lowpass_optimised, [
        (119e6, "S21", -0.249), (119e6, "S11", -16.509), (161e6, "S21", -1.032),
        (161e6, "S11", -16.348), (188e6, "S21", -40.000), (230e6, "S21", -40.267)]),
    "High-pass (optimised)": (fm.highpass_optimised, [
        (67e6, "S21", -40.135), (109e6, "S21", -40.005), (119e6, "S21", -1.910),
        (119e6, "S11", -16.291), (161e6, "S21", -0.303), (161e6, "S11", -16.295)]),
    "Low-pass prototype (ideal)": (fm.lowpass_prototype, [
        (119e6, "S21", -0.098), (119e6, "S11", -16.516), (161e6, "S21", -0.110),
        (161e6, "S11", -16.019), (188e6, "S21", -40.765), (200e6, "S21", -63.813)]),
}


def validate():
    rows = ["| Circuit | Frequency | Quantity | ADS (dB) | Python (dB) | Difference (dB) |",
            "|---|---:|---|---:|---:|---:|"]
    worst = 0.0
    for name, (circuit, markers) in ADS_MARKERS.items():
        for f, q, ads in markers:
            s11, s21 = fm.s_params(circuit(np.array([f])))
            py = fm.db(s21 if q == "S21" else s11)[0]
            worst = max(worst, abs(py - ads))
            rows.append(f"| {name} | {f / 1e6:.0f} MHz | {q} | {ads:.2f} | {py:.2f} | {py - ads:+.2f} |")
    print("\n".join(rows))
    print(f"\nworst-case difference: {worst:.2f} dB")
    return rows


# Stage 2 of the original flow ("add Q"): values as typed into the VAR block.
# The tank capacitors were given the prototype's C1..C3 (which are mostly
# shunt values) and the shunt capacitors C4..C7 the remaining ones.
STAGE2_AS_BUILT = dict(shunt_C=[26.6, 20.9, 18.9, 12.0], arm_L=[57.7, 25.7, 33.9], arm_C=[20.0, 4.8, 25.6])


def plot_stage2():
    fig, ax = plt.subplots(figsize=(8, 4.2))
    for d, label, ls in ((fm.LP_PROTOTYPE, "prototype, ideal parts", "-"),
                         (fm.LP_PROTOTYPE, "prototype values + finite Q (correct mapping)", "--"),
                         (STAGE2_AS_BUILT, "stage 2 as built (capacitor values swapped)", ":")):
        q = (np.inf, np.inf) if "ideal" in label else (fm.QL, fm.QC)
        _, s21 = fm.s_params(fm.lowpass_prototype(F, d, *q))
        ax.plot(F / 1e6, fm.db(s21), ls, label=label, lw=1.6)
    ax.set_ylim(-80, 3)
    style(ax, "Low-pass |S21|: what really broke in stage 2")
    ax.legend(loc="lower left", fontsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(DOCS, "python_stage2_diagnosis.png"), dpi=130)
    plt.close(fig)


def style(ax, title):
    ax.axvspan(*np.array(PASS) / 1e6, color="#2e86de", alpha=0.08, label="passband 119–161 MHz")
    ax.set_xlabel("Frequency (MHz)")
    ax.set_ylabel("dB")
    ax.set_title(title)
    ax.grid(alpha=0.3)
    ax.set_xlim(0, 300)


def plot_sections():
    fig, axs = plt.subplots(1, 2, figsize=(12, 4.2))
    for ax, circuit, title in ((axs[0], fm.lowpass_optimised, "Low-pass section (elliptic, N = 7)"),
                               (axs[1], fm.highpass_optimised, "High-pass section")):
        s11, s21 = fm.s_params(circuit(F))
        ax.plot(F / 1e6, fm.db(s21), label="|S21|")
        ax.plot(F / 1e6, fm.db(s11), label="|S11|")
        ax.axhline(-40, color="k", ls=":", lw=0.8)
        ax.set_ylim(-80, 3)
        style(ax, title)
        ax.legend(loc="lower left")
    fig.tight_layout()
    fig.savefig(os.path.join(DOCS, "python_sections.png"), dpi=130)
    plt.close(fig)


def plot_bandpass():
    s11, s21 = fm.s_params(fm.bandpass(F))
    gd = fm.group_delay(F, s21)
    fig, axs = plt.subplots(2, 1, figsize=(8, 7), gridspec_kw=dict(height_ratios=[2, 1]), sharex=True)
    axs[0].plot(F / 1e6, fm.db(s21), label="|S21|")
    axs[0].plot(F / 1e6, fm.db(s11), label="|S11|")
    for f_stop in (67e6, 109e6, 188e6, 230e6):
        axs[0].plot(f_stop / 1e6, fm.db(fm.s_params(fm.bandpass(np.array([f_stop])))[1])[0], "kv")
    axs[0].axhline(-40, color="k", ls=":", lw=0.8, label="−40 dB rejection target")
    axs[0].set_ylim(-90, 3)
    style(axs[0], "Band-pass filter (LP + HP cascade, QL = 100, QC = 1000)")
    axs[0].legend(loc="lower right", fontsize=8)
    band = (F > 112e6) & (F < 168e6)
    axs[1].plot(F[band] / 1e6, gd[band] * 1e9, color="C2")
    axs[1].axvspan(*np.array(PASS) / 1e6, color="#2e86de", alpha=0.08)
    axs[1].set_ylabel("group delay (ns)")
    axs[1].set_xlabel("Frequency (MHz)")
    axs[1].grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(os.path.join(DOCS, "python_bandpass.png"), dpi=130)
    plt.close(fig)

    inband = (F >= PASS[0]) & (F <= PASS[1])
    il = -fm.db(s21[inband])
    print(f"band-pass (finite Q): IL {il.min():.2f}..{il.max():.2f} dB, "
          f"worst return loss {-fm.db(s11[inband]).max():.1f} dB, "
          f"group delay {gd[inband].min() * 1e9:.1f}..{gd[inband].max() * 1e9:.1f} ns")


def perturb(d, tol, rng):
    return {k: [v * (1 + rng.uniform(-tol, tol)) for v in vals] for k, vals in d.items()}


def monte_carlo(n=2000, seed=1):
    """Uniformly distributed tolerances on every L and C."""
    rng = np.random.default_rng(seed)
    f = np.array([67e6, 109e6, 119e6, 140e6, 161e6, 188e6, 230e6])
    tolerances = [0.01, 0.02, 0.05]
    results = {}
    for tol in tolerances:
        rej, il = [], []
        for _ in range(n):
            _, s21 = fm.s_params(fm.bandpass(f, perturb(fm.LP_OPTIMISED, tol, rng),
                                             perturb(fm.HP_OPTIMISED, tol, rng)))
            s = -fm.db(s21)
            rej.append(min(s[1], s[5]))              # worst rejection at the near stopband edges
            il.append(max(s[2], s[3], s[4]))         # worst insertion loss in the passband
        results[tol] = (np.array(rej), np.array(il))

    fig, axs = plt.subplots(1, 2, figsize=(12, 4))
    for tol in tolerances:
        rej, il = results[tol]
        axs[0].hist(rej, bins=50, alpha=0.55, label=f"±{tol:.0%}")
        axs[1].hist(il, bins=50, alpha=0.55, label=f"±{tol:.0%}")
    axs[0].axvline(40, color="k", ls="--", lw=1)
    axs[0].set_xlabel("worst rejection at 109 / 188 MHz (dB)")
    axs[1].set_xlabel("worst passband insertion loss (dB)")
    for ax in axs:
        ax.set_ylabel("count")
        ax.grid(alpha=0.3)
        ax.legend(title="L, C tolerance")
    fig.suptitle(f"Monte-Carlo tolerance analysis ({n} samples per tolerance)")
    fig.tight_layout()
    fig.savefig(os.path.join(DOCS, "python_monte_carlo.png"), dpi=130)
    plt.close(fig)

    rows = ["| Tolerance | P(rejection ≥ 40 dB) | P(rejection ≥ 35 dB) | Median insertion loss (dB) | 95th percentile insertion loss (dB) |",
            "|---|---:|---:|---:|---:|"]
    for tol in tolerances:
        rej, il = results[tol]
        rows.append(f"| ±{tol:.0%} | {np.mean(rej >= 40):.0%} | {np.mean(rej >= 35):.0%} | "
                    f"{np.median(il):.2f} | {np.percentile(il, 95):.2f} |")
    print("\n".join(rows))
    return rows


if __name__ == "__main__":
    sys.stdout.reconfigure(encoding="utf-8")
    os.makedirs(DOCS, exist_ok=True)
    validate()
    plot_sections()
    plot_stage2()
    plot_bandpass()
    monte_carlo()
