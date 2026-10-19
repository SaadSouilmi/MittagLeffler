# /// script
# requires-python = ">=3.12"
# dependencies = ["numpy==2.5.4", "scipy==1.18.1", "pymittagleffler==0.2.1"]
# ///
import json
import math
import sys
from decimal import Decimal, localcontext
from importlib.metadata import version
from pathlib import Path

import numpy as np
from pymittagleffler import mittag_leffler
from scipy.special import erfcx

results = Path(sys.argv[1])
old = np.genfromtxt(results / "baseline-samples.csv", delimiter=",", names=True)
new = np.genfromtxt(results / "real-samples.csv", delimiter=",", names=True)
assert old.size == new.size == 72762
assert all(np.array_equal(old[key], new[key]) for key in ("alpha", "beta", "x"))
assert np.isfinite(old["value"]).all() and np.isfinite(new["value"]).all()
reference = np.empty(new.size)
imaginary = 0.0
for alpha, beta in dict.fromkeys(zip(new["alpha"], new["beta"])):
    mask = (new["alpha"] == alpha) & (new["beta"] == beta)
    x = new["x"][mask]
    if alpha == 1:
        with localcontext() as context:
            context.prec = 80
            ref = []
            for value in x:
                dx = Decimal.from_float(float(value))
                ref.append(float((-dx).exp() if beta == 1 else
                                 (1 - (-dx).exp()) / dx if value else Decimal(1)))
            ref = np.asarray(ref)
    elif alpha == 0.5 and beta == 1:
        ref = erfcx(x)
    else:
        ref = np.empty(x.size)
        positive = x > 0
        values = mittag_leffler(-x[positive], float(alpha), float(beta))
        assert np.isfinite(values).all()
        imaginary = max(imaginary, float(np.abs(values.imag).max()))
        ref[positive] = values.real
        ref[~positive] = 1 / math.gamma(beta)
    reference[mask] = ref


def location(index):
    return {key: float(new[key][index]) for key in ("alpha", "beta", "x")}


def errors(actual, expected, mask):
    indices = np.flatnonzero(mask)
    absolute = np.abs(actual[mask] - expected[mask])
    relative = np.divide(absolute, np.abs(expected[mask]), out=np.zeros(absolute.shape),
                         where=expected[mask] != 0)
    return {"max_abs": float(absolute.max()), "max_rel_nonzero": float(relative.max()),
            "max_abs_at": location(indices[absolute.argmax()]),
            "max_rel_at": location(indices[relative.argmax()])}


summary = {
    "baseline_commit": (results / "baseline-commit.txt").read_text().strip(),
    "versions": {name: version(name) for name in ("numpy", "scipy", "pymittagleffler")},
    "sample_count": int(new.size), "reference_max_imag": imaginary,
    "cpu_affinity": (results / "cpu-affinity.txt").read_text().strip(),
    "accuracy": {}, "evaluation": [], "initialization": [],
}
for name, mask in (("original", np.arange(new.size) < 3624),
                   ("expanded", np.arange(new.size) >= 3624)):
    summary["accuracy"][name] = {
        "samples": int(mask.sum()),
        "baseline_reference": errors(old["value"], reference, mask),
        "real_reference": errors(new["value"], reference, mask),
        "real_baseline": errors(new["value"], old["value"], mask),
    }
summary["original_pairs"] = []
for alpha, beta in dict.fromkeys(zip(new["alpha"][:3624], new["beta"][:3624])):
    mask = (np.arange(new.size) < 3624) & (new["alpha"] == alpha) & (new["beta"] == beta)
    summary["original_pairs"].append({
        "alpha": float(alpha), "beta": float(beta),
        "baseline": errors(old["value"], reference, mask),
        "real": errors(new["value"], reference, mask),
        "difference": errors(new["value"], old["value"], mask),
    })

for kind in ("evaluation", "initialization"):
    tables = {}
    for variant in ("baseline", "real"):
        table = {}
        for run in range(1, 6):
            if kind == "evaluation":
                for line in (results / f"{variant}-throughput-{run}.txt").read_text().splitlines():
                    fields = line.split()
                    if len(fields) == 8 and fields[2] in ("quadrature", "exp", "expm1"):
                        key = (float(fields[0]), float(fields[1]))
                        table.setdefault(key, []).append(float(fields[4]))
            else:
                rows = np.genfromtxt(results / f"{variant}-initialization-{run}.csv", delimiter=",", names=True)
                for row in rows:
                    key = (float(row["alpha"]), float(row["beta"]))
                    table.setdefault(key, []).append(float(row["ns_init"]))
        assert len(table) == (11 if kind == "evaluation" else 9)
        assert all(len(times) == 5 and min(times) > 0 for times in table.values())
        tables[variant] = table
    assert tables["baseline"].keys() == tables["real"].keys()
    for key in tables["baseline"]:
        entry = {"alpha": key[0], "beta": key[1]}
        for variant in ("baseline", "real"):
            times = tables[variant][key]
            entry[variant] = {"median_ns": float(np.median(times)), "min_ns": min(times), "max_ns": max(times)}
        entry["speedup"] = entry["baseline"]["median_ns"] / entry["real"]["median_ns"]
        summary[kind].append(entry)

delta = np.abs(new["value"] - old["value"])
error = np.abs(new["value"] - reference)
old_error = np.abs(old["value"] - reference)
summary["checks"] = {
    "baseline_agreement": bool(np.all(delta <= 32 * np.finfo(float).eps * np.maximum(1, np.abs(old["value"])))),
    "reference_absolute_tolerance_3e_11": bool(np.all(error <= 3e-11 + 2e-13 * np.abs(reference))),
    "reference_regression_tolerance_1e_14": bool(np.all(error <= old_error + 1e-14)),
    "shortcuts_unchanged": bool(np.array_equal(new["value"][new["alpha"] == 1], old["value"][new["alpha"] == 1])),
}
np.savetxt(results / "reference.csv", np.column_stack((new["alpha"], new["beta"], new["x"], reference)),
           delimiter=",", header="alpha,beta,x,reference", comments="", fmt="%.17g")
(results / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
print(json.dumps(summary, indent=2))
assert all(summary["checks"].values()), summary["checks"]
