# /// script
# requires-python = ">=3.11"
# dependencies = [
#   "numpy==2.5.4", "matplotlib==3.11.2", "scipy==1.18.1",
#   "pymittagleffler==0.2.1",
# ]
# ///
import argparse
import io
import json
import math
from decimal import Decimal, localcontext
from importlib.metadata import version
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from pymittagleffler import mittag_leffler
from scipy.special import erfcx

parser = argparse.ArgumentParser()
parser.add_argument("transcript", type=Path)
parser.add_argument("--figures", type=Path, required=True)
parser.add_argument("--data", type=Path, required=True)
args = parser.parse_args()
args.figures.mkdir(parents=True, exist_ok=True)
args.data.mkdir(parents=True, exist_ok=True)
text = "\n".join(line.rstrip() for line in args.transcript.read_text().splitlines())

sample_text = text.split("ML_SAMPLES_BEGIN\n", 1)[1].split("ML_SAMPLES_END", 1)[0]
samples = np.genfromtxt(io.StringIO(sample_text), delimiter=",", names=True)
assert samples.size == 12 * 302, f"Incomplete samples: {samples.size}"
assert all(np.isfinite(samples[name]).all() for name in samples.dtype.names)
bench = {}
for line in text.splitlines():
    fields = line.split()
    if len(fields) == 8 and fields[2] in ("quadrature", "exp", "expm1"):
        key = (float(fields[0]), float(fields[1]))
        bench.setdefault(key, []).append([float(fields[3]), float(fields[4])])
assert len(bench) == 11 and all(len(rows) == 5 for rows in bench.values())
assert all(np.isfinite(rows).all() and (np.asarray(rows) > 0).all() for rows in bench.values())

colors = ["#2266A1", "#D17A22", "#25856F"]
plt.rcParams.update({
    "font.size": 10, "axes.titlesize": 12, "axes.labelsize": 11,
    "axes.spines.top": False, "axes.spines.right": False,
    "axes.grid": True, "grid.alpha": 0.2, "grid.linewidth": 0.6,
    "figure.facecolor": "white", "savefig.facecolor": "white",
})

def save(fig, name):
    fig.savefig(args.figures / f"{name}.png", dpi=170, bbox_inches="tight")
    plt.close(fig)

reference = np.empty(samples.size)
imag_max = 0.0
summary = {"versions": {name: version(name) for name in ("numpy", "matplotlib", "scipy", "pymittagleffler")}}
metadata = text.split("ML_METADATA_BEGIN\n", 1)[1].split("ML_METADATA_END", 1)[0].strip().splitlines()
cpu = next(line.partition(":")[2].strip() for line in metadata if line.startswith("Model name:"))
cpu = cpu.replace("(R)", "").replace("CPU @", "@").strip()
compiler = next(line for line in metadata if "g++" in line or "clang" in line)
affinity = next(line.split("=", 1)[1] for line in metadata if line.startswith("CPU_AFFINITY="))
summary["benchmark_metadata"] = {"cpu": cpu, "compiler": compiler, "cpu_affinity": affinity}
errors = []
for alpha, beta in dict.fromkeys(zip(samples["alpha"], samples["beta"])):
    mask = (samples["alpha"] == alpha) & (samples["beta"] == beta)
    x = samples["x"][mask]
    ref = np.empty(x.shape, dtype=complex)
    if alpha == 1:
        with localcontext() as context:
            context.prec = 80
            for index, value in enumerate(x):
                dx = Decimal.from_float(float(value))
                ref[index] = float((-dx).exp() if beta == 1 else
                                   (1 - (-dx).exp()) / dx if value else Decimal(1))
        reference_name = "80-digit Decimal elementary identity"
    else:
        ref[x > 0] = mittag_leffler(-x[x > 0], float(alpha), float(beta))
        reference_name = "pymittagleffler 0.2.1; reciprocal gamma at zero"
    ref[x == 0] = 1 / math.gamma(beta)
    assert np.isfinite(ref).all()
    imag_max = max(imag_max, float(np.max(np.abs(ref.imag))))
    reference[mask] = ref.real
    error = np.abs(samples["value"][mask] - ref.real)
    relative = np.divide(error, np.abs(ref.real), out=np.full(error.shape, np.nan), where=ref.real != 0)
    rows = {"alpha": float(alpha), "beta": float(beta), "reference": reference_name, "max_abs": float(error.max()),
            "max_rel_nonzero": float(np.nanmax(relative)), "x_max_abs": float(x[error.argmax()])}
    core = (x >= 1e-6) & (x <= 100)
    rows["core_max_abs"] = float(error[core].max())
    rows["core_max_rel"] = float(np.nanmax(relative[core]))
    errors.append(rows)
summary["reference_max_imag"] = imag_max
summary["errors"] = errors
mask = (samples["alpha"] == 0.5) & (samples["beta"] == 1)
exact = erfcx(samples["x"][mask])
exact_error = np.abs(samples["value"][mask] - exact)
summary["erfcx_check"] = {"max_abs": float(exact_error.max()), "max_rel": float(np.max(exact_error / exact))}
np.savetxt(args.data / "accuracy.csv", np.column_stack([samples[name] for name in samples.dtype.names] + [reference]),
           delimiter=",", header="alpha,beta,x,value,reference", comments="", fmt="%.17g")

fig, axes = plt.subplots(1, 3, figsize=(11.8, 3.8), layout="constrained")
for col, (title, beta_func) in enumerate([(r"$\beta=\alpha$", lambda a: a),
                                        (r"$\beta=1$", lambda a: 1.0),
                                        (r"$\beta=\alpha+1$", lambda a: a+1)]):
    ax = axes[col]
    for color, alpha in zip(colors, [0.55, 0.75, 0.95]):
        mask = (samples["alpha"] == alpha) & np.isclose(samples["beta"], beta_func(alpha)) & (samples["x"] > 0)
        ax.loglog(samples["x"][mask], samples["value"][mask], color=color, label=fr"$\alpha={alpha}$", lw=1.6)
    ax.set(title=title, xlabel=r"$X$ (argument is $-X$)", xlim=(1e-8, 1e4))
    ax.legend(frameon=False, loc="lower left")
axes[0].set_ylabel(r"$E_{\alpha,\beta}(-X)$")
fig.suptitle("Values computed by MLQuadEval", fontsize=14)
save(fig, "values")

fig, axes = plt.subplots(2, 3, figsize=(11.8, 6.6), layout="constrained")
for col, (title, beta_func) in enumerate([(r"$\beta=\alpha$", lambda a: a),
                                        (r"$\beta=1$", lambda a: 1.0),
                                        (r"$\beta=\alpha+1$", lambda a: a+1)]):
    for color, alpha in zip(colors, [0.55, 0.75, 0.95]):
        mask = (samples["alpha"] == alpha) & np.isclose(samples["beta"], beta_func(alpha)) & (samples["x"] > 0)
        x = samples["x"][mask]
        absolute = np.abs(samples["value"][mask] - reference[mask])
        relative = absolute / np.abs(reference[mask])
        for row, error in enumerate([absolute, relative]):
            axes[row, col].loglog(x, np.ma.masked_where(error == 0, error), color=color, lw=1.3, label=fr"$\alpha={alpha}$")
            axes[row, col].axvspan(1e-6, 100, color="#DDE5EC", alpha=0.35, zorder=0)
            axes[row, col].set_xlim(1e-8, 1e4)
    axes[0, col].set_title(title)
    axes[0, col].legend(frameon=False, fontsize=9)
    axes[1, col].set_xlabel(r"$X$ (argument is $-X$)")
axes[0, 0].set_ylabel("Absolute difference")
axes[1, 0].set_ylabel("Relative difference")
fig.suptitle("Comparison with pymittagleffler 0.2.1\nShading marks the throughput input range; exact zeros are omitted", fontsize=13)
save(fig, "accuracy")

fig, ax = plt.subplots(figsize=(9.2, 6.0), layout="constrained")
keys = list(bench)
medians = np.array([np.median(bench[key], axis=0)[0] / 1e6 for key in keys])
mins = np.array([np.min(bench[key], axis=0)[0] / 1e6 for key in keys])
maxs = np.array([np.max(bench[key], axis=0)[0] / 1e6 for key in keys])
ax.barh(range(len(keys)), medians, color=[colors[0]]*9 + [colors[2]]*2,
        xerr=np.vstack((medians-mins, maxs-medians)), capsize=3, height=0.65)
for index, value in enumerate(medians):
    ax.text(maxs[index]+1.6, index, f"{value:.1f}", va="center", fontsize=10)
labels = [fr"$({a:.2f}, {b:.2f})$" + ("  exp" if (a,b)==(1,1) else "  expm1" if (a,b)==(1,2) else "") for a,b in keys]
ax.set_yticks(range(len(keys)), labels)
ax.invert_yaxis()
ax.set(xlabel="Million evaluations / second", ylabel=r"$(\alpha,\beta)$", xlim=(0, max(maxs)*1.12))
ax.set_title("Single-thread evaluation throughput\nMedian of five runs; whiskers show min/max", loc="left", pad=12)
ax.grid(axis="y", visible=False)
fig.supxlabel(f"{cpu} · {compiler.split()[0]} {compiler.split()[-1]} · CPU {affinity}\n"
              "-O3 -march=native · 4,096 inputs in [1e-6, 100]\n"
              "Initialization excluded · checksum accumulation included", fontsize=9)
save(fig, "throughput")
summary["throughput"] = [{"alpha": a, "beta": b, "median_eval_s": float(np.median(bench[(a,b)],axis=0)[0]),
                           "min_eval_s": float(np.min(bench[(a,b)],axis=0)[0]), "max_eval_s": float(np.max(bench[(a,b)],axis=0)[0]),
                           "median_ns_eval": float(np.median(bench[(a,b)],axis=0)[1])} for a,b in keys]

fig, ax = plt.subplots(figsize=(7.5, 4.2), layout="constrained")
u = np.linspace(-6, 6, 601)
s = (1 + 1j*u)**2
nodes = (1 + 1j*0.2*np.arange(-30, 31))**2
ax.plot(s.real, s.imag, color=colors[0], lw=1.5)
ax.scatter(nodes.real, nodes.imag, s=15, color=colors[0], zorder=3, label="61 trapezoidal nodes")
ax.plot([-38, 0], [0, 0], ls="--", color="#8C5368", lw=1.5, label="Principal-branch cut")
ax.scatter([1], [0], s=42, color=colors[1], zorder=4)
ax.annotate(r"$s(0)=1$", (1, 0), xytext=(5, 3), arrowprops={"arrowstyle":"->", "color":"#555555"})
ax.annotate("Conjugate pairs share\none real evaluation", (-18, 8.7), xytext=(-34, 2.2), arrowprops={"arrowstyle":"->", "color":"#555555"})
ax.set(xlabel=r"$\mathrm{Re}\,s$", ylabel=r"$\mathrm{Im}\,s$", xlim=(-38, 10), ylim=(-14, 14),
       title=r"Fixed contour: $s(u)=(1+iu)^2$, $u_k=0.2k$, $-30\leq k\leq30$")
ax.legend(frameon=False, loc="lower right")
save(fig, "contour")

(args.data / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
print(json.dumps(summary, indent=2))
