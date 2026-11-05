"""Reproducible M3 controls. Run: python -m prat_growth.attribution."""
import csv
import json
import time
from pathlib import Path
from dataclasses import asdict
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from .parameters import PhysicalParameters, NumericalParameters
from .analytic import general_growth_velocity
from .dynamic import DynamicALEModel, snapshot_diagnostics
from .frozen import solve_frozen
from .verification import local_audit
from .provenance import identity, cache_matches, source_hash
from .campaign import save_result, load_result
from .trajectory_checks import startup, independent_time_checks

TIMES = (300., 3600., 10800., 21600.)


def table(path, rows):
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def l95(result, index):
    d = snapshot_diagnostics(result, index)
    masses = np.maximum(d["reaction"], 0)*result.height[index]*result.widths
    cum = np.r_[0., np.cumsum(masses[::-1])]
    return float(np.interp(.95*cum[-1], cum, result.height[index]*(1-result.faces[::-1])))


def state_vector(r, i):
    return np.r_[r.epsilon[:, i], r.concentration[:, i], r.height[i], r.top_solid[i],
                 r.internal_solid[i], r.bottom_salt[i], r.bottom_solution[i], r.evaporated[i]]


def execute(out, name, p, n, closure, initial="equilibrium"):
    meta, data = out/"runs"/(name+".json"), out/"runs"/(name+".npz")
    if cache_matches(meta, data, identity(p, n, closure, initial)):
        r = load_result(out, name)
    else:
        start = time.perf_counter()
        r = DynamicALEModel(p, n, closure).run(initial)
        if r.status != 0 or abs(r.time[-1]-n.t_end) > 1e-8:
            raise RuntimeError(f"{name}: {r.message}")
        save_result(out, name, r)
        print(f"{name}: {time.perf_counter()-start:.2f} s", flush=True)
    return r


def run(out=Path("outputs/m3")):
    out.mkdir(parents=True, exist_ok=True)
    p = PhysicalParameters()
    base = NumericalParameters(n_cells=640, grid_power=3., max_step=60., rtol=1e-8,
                               atol_field=1e-12, atol_inventory=1e-12, atol_height=1e-14)
    runs, convergence = {}, []
    for count in (160, 240, 320, 480, 640):
        name = f"grid_{count}"
        runs[name] = execute(out, name, p, base.with_updates(n_cells=count), "B2")
    runs["step_30"] = execute(out, "step_30", p, base.with_updates(max_step=30), "B2")
    runs["step_120"] = execute(out, "step_120", p, base.with_updates(max_step=120), "B2")
    runs["legacy_grid_corrected"] = execute(out, "legacy_grid_corrected", p,
        base.with_updates(n_cells=160, grid_power=2.5, rtol=1e-7, atol_field=1e-10,
                          atol_inventory=1e-10, atol_height=1e-12), "B2")
    runs["T2"] = runs["grid_640"]
    for label, closure, initial in (("T0", "B1", "equilibrium"), ("T1", "B1", "quasi_steady"),
                                     ("T3", "B2", "quasi_steady")):
        runs[label] = execute(out, label, p, base, closure, initial)
    audit_rows, transient, qrows, contributions, profiles = [], [], [], [], []
    for name, r in runs.items():
        d = snapshot_diagnostics(r, -1)
        q = solve_frozen(p, r.numerical, r.height[-1], r.epsilon[:, -1], "Q6" if r.closure == "B2" else "Q3")
        convergence.append(dict(case=name, n_cells=r.numerical.n_cells, power=r.numerical.grid_power,
            max_step=r.numerical.max_step, velocity=d["velocity"], q_velocity=q.velocity,
            h=r.height[-1], epsilon_min=d["epsilon_min"], L95_interpolated=l95(r, -1),
            cells_delta=d["cells_across_delta"], internal_fraction=r.internal_solid[-1]/(r.internal_solid[-1]+r.top_solid[-1]),
            salt_residual=d["salt_residual_normalized"], bulk_residual=d["bulk_residual_normalized"]))
        m = DynamicALEModel(p, r.numerical, r.closure)
        # Audit every saved output, including startup, using independently rebuilt faces.
        for i, t in enumerate(r.time):
            audit_rows.append(dict(case=name, time_s=t, **local_audit(m, state_vector(r, i))))
    for label in ("T0", "T1", "T2", "T3"):
        r = runs[label]
        for i, t in enumerate(r.time):
            d = snapshot_diagnostics(r, i)
            transient.append(dict(case=label, time_s=t, velocity=d["velocity"],
                supersaturation=d["c_top"]-p.C_eq, h=r.height[i], epsilon_min=d["epsilon_min"]))
    # Repeat attribution on both finest grids to bound error in each contrast.
    for grid in (480, 640):
        r = runs[f"grid_{grid}"]
        for t in TIMES:
            i = int(np.argmin(abs(r.time-t)))
            h, eps = r.height[i], r.epsilon[:, i]
            A = general_growth_velocity(p, h)
            qs = {f"Q{k}": solve_frozen(p, r.numerical, h, eps, f"Q{k}") for k in range(7)}
            qs["Q0R"] = solve_frozen(p, r.numerical, h, eps, "Q0", linear_reaction=False)
            qs["Q0B"] = solve_frozen(p, r.numerical, h, eps, "Q0", linear_boundary=False)
            for label, q in qs.items():
                qrows.append(dict(n_cells=grid, time_s=t, case=label, h=h, velocity=q.velocity,
                    supersaturation=q.c_top-p.C_eq, integrated_reaction=np.sum(q.reaction*h*r.widths),
                    residual_inf=q.residual_inf, interface_residual=q.interface_residual))
                if grid == 640:
                    for j, x in enumerate(r.centers):
                        profiles.append(dict(time_s=t, case=label, distance_front_um=h*(1-x)*1e6,
                            supersaturation=q.concentration[j]-p.C_eq, reaction=q.reaction[j], epsilon=q.epsilon[j]))
            v = {k: q.velocity for k, q in qs.items()}
            dynamic_v = snapshot_diagnostics(r, i)["velocity"]
            parts = {"Q0_discretization": v["Q0"]-A,
                "reaction_nonlinearity": v["Q0R"]-v["Q0"],
                "boundary_nonlinearity": v["Q0B"]-v["Q0"],
                "nonlinear_interaction": v["Q1"]-v["Q0R"]-v["Q0B"]+v["Q0"],
                "advection_ALE": v["Q2"]-v["Q1"], "porosity": v["Q3"]-v["Q2"],
                "diffusivity": v["Q4"]-v["Q3"], "surface_area": v["Q5"]-v["Q3"],
                "D_av_interaction": v["Q6"]-v["Q4"]-v["Q5"]+v["Q3"],
                "transient_T2_minus_Q6": dynamic_v-v["Q6"]}
            for label, value in parts.items():
                contributions.append(dict(n_cells=grid, time_s=t, mechanism=label,
                    delta_velocity=value, percentage_points=100*value/A))
            assert abs(sum(parts.values())-(dynamic_v-A)) < 1e-20
            if grid == 640:
                d = snapshot_diagnostics(r, i)
                for j, x in enumerate(r.centers):
                    profiles.append(dict(time_s=t, case="T2", distance_front_um=h*(1-x)*1e6,
                        supersaturation=r.concentration[j, i]-p.C_eq, reaction=d["reaction"][j], epsilon=eps[j]))
    # Matching frozen snapshots for all four transients, not B1 against B2 structure.
    matched = []
    for label in ("T0", "T1", "T2", "T3"):
        r = runs[label]
        for t in TIMES:
            i = int(np.argmin(abs(r.time-t)))
            q = solve_frozen(p, base, r.height[i], r.epsilon[:, i], "Q3" if r.closure == "B1" else "Q6")
            v = snapshot_diagnostics(r, i)["velocity"]
            matched.append(dict(case=label, time_s=t, velocity=v, frozen_velocity=q.velocity,
                                relative_transient=100*(v/q.velocity-1)))
    for name, rows in (("convergence", convergence), ("local_audits", audit_rows), ("transients", transient),
                       ("frozen", qrows), ("contributions", contributions), ("profiles", profiles), ("matched", matched)):
        table(out/(name+".csv"), rows)
    early, early_states, early_warnings = startup(p, base)
    table(out/"startup.csv", early)
    np.savez_compressed(out/"startup_states.npz", **early_states)
    time_checks, time_warnings = independent_time_checks(runs["T2"], state_vector)
    table(out/"trajectory_residuals.csv", time_checks)
    (out/"manifest.json").write_text(json.dumps(dict(code_sha256=source_hash(), physical=asdict(p),
        numerical=asdict(base), times=TIMES, T2_data_alias="grid_640", numpy=np.__version__,
        warnings=[dict(case=name, warnings=r.solver_warnings) for name, r in runs.items()]+early_warnings+time_warnings,
        max_residuals={key: max(row[key] for row in audit_rows) for key in audit_rows[0] if key not in ("case", "time_s")}), indent=2))
    figures(out, transient, contributions, profiles, convergence)
    fig, axes = plt.subplots(1, 2, figsize=(10, 4), layout="constrained")
    for label in ("T0", "T1", "T2", "T3"):
        rows = [r for r in early if r["case"] == label and r["time_s"] > 0]
        for ax, key, factor in ((axes[0], "velocity", 86400e3), (axes[1], "supersaturation", 1e6)):
            ax.semilogx([r["time_s"] for r in rows], [r[key]*factor for r in rows], "--" if label in ("T1", "T3") else "-", label=label)
            ax.set_xlabel("Time (s)")
    axes[0].set_ylabel("Growth velocity (mm/day)")
    axes[1].set_ylabel(r"Surface excess mass fraction ($10^{-6}$)")
    axes[0].legend()
    savefig(fig, out/"startup")
    print("M3 study complete", flush=True)


def figures(out, transient, contributions, profiles, convergence):
    plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False})
    colors = {f"Q{i}": plt.get_cmap("tab10")(i) for i in range(7)}
    colors["T2"] = "black"
    fig, axes = plt.subplots(1, 2, figsize=(10, 4), layout="constrained")
    for label in ("T0", "T1", "T2", "T3"):
        rows = [r for r in transient if r["case"] == label and r["time_s"] >= 300]
        style = "--" if label in ("T1", "T3") else "-"
        axes[0].plot([r["time_s"]/3600 for r in rows], [r["velocity"]*86400e3 for r in rows], style, label=label)
        axes[1].plot([r["time_s"]/3600 for r in rows], [r["supersaturation"]*1e6 for r in rows], style, label=label)
    axes[0].axhline(general_growth_velocity(PhysicalParameters(), .001)*86400e3, color="k", ls=":", label="Prat")
    for ax in axes:
        ax.set_xlabel("Time (h)")
        ax.legend()
    axes[0].set_ylabel("Growth velocity (mm/day)")
    axes[1].set_ylabel(r"Surface excess mass fraction ($10^{-6}$)")
    savefig(fig, out/"velocity_supersaturation")
    fig, axes = plt.subplots(2, 2, figsize=(11, 7), layout="constrained")
    for t, ax in zip(TIMES, axes.flat):
        for label in ("Q0", "Q2", "Q3", "Q4", "Q5", "Q6", "T2"):
            rows = [r for r in profiles if r["case"] == label and r["time_s"] == t and r["distance_front_um"] < 3]
            ax.plot([r["distance_front_um"] for r in rows], [r["supersaturation"]*1e6 for r in rows], "--" if label == "T2" else "-", label=label, color=colors[label])
        ax.set(title="t = 5 min" if t == 300 else f"t = {t/3600:g} h", xlabel="Distance below front (micrometres)", ylabel=r"Excess mass fraction ($10^{-6}$)")
    axes[0, 0].legend(ncol=2)
    savefig(fig, out/"concentration_profiles")
    fig, axes = plt.subplots(2, 2, figsize=(11, 7), layout="constrained")
    for t, ax in zip(TIMES, axes.flat):
        for label in ("Q0", "Q3", "Q4", "Q5", "Q6", "T2"):
            rows = [r for r in profiles if r["case"] == label and r["time_s"] == t and r["distance_front_um"] < 3]
            ax.plot([r["distance_front_um"] for r in rows], [r["reaction"] for r in rows], "--" if label == "T2" else "-", label=label, color=colors[label])
        ax.set(title="t = 5 min" if t == 300 else f"t = {t/3600:g} h", xlabel="Distance below front (micrometres)", ylabel=r"Precipitation ($kg\,m^{-3}\,s^{-1}$)")
    axes[0, 0].legend(ncol=2)
    savefig(fig, out/"precipitation_profiles")
    fig, ax = plt.subplots(figsize=(10, 5), layout="constrained")
    rows = [r for r in contributions if r["n_cells"] == 640 and r["time_s"] == 21600]
    ax.barh([r["mechanism"].replace("_", " ") for r in rows], [r["percentage_points"] for r in rows])
    ax.axvline(0, color="k", lw=.7)
    ax.set_xlabel("Contribution to (V - V Prat) / V Prat (percentage points)")
    savefig(fig, out/"attribution_6h")
    fig, axes = plt.subplots(1, 2, figsize=(10, 4), layout="constrained")
    rows = [r for r in convergence if r["case"].startswith("grid_")]
    for key, label in (("velocity", "Transient"), ("q_velocity", "Frozen Q6")):
        axes[0].plot([r["n_cells"] for r in rows], [r[key]*86400e3 for r in rows], "o-", label=label)
    rows = [r for r in convergence if r["case"] in ("step_30", "step_120", "grid_640")]
    rows.sort(key=lambda r: r["max_step"])
    reference = next(r["velocity"] for r in rows if r["max_step"] == 60)
    axes[1].plot([r["max_step"] for r in rows], [(r["velocity"]/reference-1)*1e12 for r in rows], "o-")
    axes[0].set_xlabel("Number of cells")
    axes[1].set_xlabel("Maximum time step (s)")
    for ax in axes:
        ax.set_ylabel("Growth velocity (mm/day)")
        ax.ticklabel_format(useOffset=False, axis="y")
    axes[1].set_ylabel(r"Relative velocity change ($10^{-12}$)")
    axes[1].set_title("Numerically indistinguishable")
    axes[1].set_ylim(-6, 1)
    axes[0].legend()
    savefig(fig, out/"convergence")


def savefig(fig, path):
    fig.savefig(path.with_suffix(".png"), dpi=180)
    fig.savefig(path.with_suffix(".svg"))
    plt.close(fig)


if __name__ == "__main__":
    run()
