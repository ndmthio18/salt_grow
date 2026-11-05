"""Fixed-laboratory-coordinate material diagnostics for the ALE solution."""
from __future__ import annotations

from dataclasses import dataclass
import warnings

import numpy as np
from scipy.integrate import cumulative_trapezoid, solve_ivp


@dataclass
class WindowResult:
    time: np.ndarray
    states: np.ndarray
    accepted_time: np.ndarray
    nfev: int
    njev: int
    nlu: int
    warnings: list[str]


def state_from_result(result, index: int) -> np.ndarray:
    return np.r_[
        result.epsilon[:, index], result.concentration[:, index],
        result.height[index], result.top_solid[index], result.internal_solid[index],
        result.bottom_salt[index], result.bottom_solution[index], result.evaporated[index],
    ]


def restart_window(model, start_time: float, y0: np.ndarray, duration: float,
                   output_interval: float, max_step: float = 1.0,
                   rtol: float = 1e-9) -> WindowResult:
    """Restart without changing inventories; dense samples are not solver steps."""
    n = model.n.n_cells
    atol = np.r_[np.full(2*n, 1e-12), 1e-14, np.full(model.N_INVENTORIES, 1e-12)]
    with warnings.catch_warnings(record=True) as captured:
        warnings.simplefilter("always", RuntimeWarning)
        sol = solve_ivp(
            model.rhs, (start_time, start_time + duration), y0, method="BDF",
            dense_output=True, rtol=rtol, atol=atol, max_step=max_step,
        )
    if not sol.success or not np.all(np.isfinite(sol.y)):
        raise RuntimeError(sol.message)
    count = int(np.ceil(duration/output_interval))
    sample_time = start_time + np.linspace(0.0, duration, count + 1)
    return WindowResult(
        time=sample_time, states=sol.sol(sample_time), accepted_time=sol.t,
        nfev=sol.nfev, njev=sol.njev, nlu=sol.nlu,
        warnings=sorted(set(str(item.message) for item in captured)),
    )


def reconstruct(z0: float, height: float, faces: np.ndarray, centers: np.ndarray,
                cell_values: np.ndarray, top_value: float, method: str = "linear") -> float:
    """Reconstruct a bounded cell-average field at fixed physical coordinate z0."""
    if z0 < -1e-14 or z0 > height*(1 + 1e-12):
        raise ValueError(f"material coordinate z0={z0} outside [0,h={height}]")
    if z0 >= height*(1 - 1e-13):
        return float(top_value)
    if method == "piecewise_constant":
        cell = int(np.clip(np.searchsorted(height*faces, z0, side="right") - 1,
                           0, len(cell_values)-1))
        return float(cell_values[cell])
    if method != "linear":
        raise ValueError(method)
    zc = height*centers
    value = float(np.interp(z0, np.r_[zc, height], np.r_[cell_values, top_value]))
    local_min = min(float(cell_values[-1]), float(top_value)) if z0 > zc[-1] else float(np.min(cell_values))
    local_max = max(float(cell_values[-1]), float(top_value)) if z0 > zc[-1] else float(np.max(cell_values))
    return float(np.clip(value, local_min, local_max))


def evaluate_marker(model, state: np.ndarray, z0: float, method: str = "linear") -> dict[str, float]:
    p = model.p
    eps, concentration, height, _ = model.unpack(state)
    av_cell, d_cell = model.coefficients(eps)
    top = model.flux_state(eps, concentration, height)
    c_top = float(top["c_top"])
    eps_z = reconstruct(z0, height, model.faces, model.centers, eps, p.epsilon_birth, method)
    c_z = reconstruct(z0, height, model.faces, model.centers, concentration, c_top, method)
    if model.closure == "B1":
        av_z, d_z = p.a_v_ref, p.D_eff_ref
    else:
        av_z, d_z = 36.0*(1.0-eps_z)/p.d_sub, np.sqrt(eps_z)*p.D_molecular
    reaction = p.rho_solution*av_z*p.k_reaction*(c_z-p.C_eq)
    return {
        "h_m": height, "x": z0/height, "d_m": height-z0,
        "epsilon": eps_z, "supersaturation": c_z-p.C_eq,
        "reaction_kg_m3_s": reaction, "a_v_m_inv": av_z,
        "D_eff_m2_s": d_z, "velocity_m_s": float(top["velocity"]),
    }


def marker_history(model, window: WindowResult, z0: float, birth_time: float | None,
                   method: str = "linear", initial_epsilon: float | None = None) -> list[dict]:
    records = [evaluate_marker(model, window.states[:, i], z0, method)
               for i in range(window.time.size)]
    reaction = np.array([row["reaction_kg_m3_s"] for row in records])
    signed = cumulative_trapezoid(reaction, window.time, initial=0.0)
    positive = cumulative_trapezoid(np.maximum(reaction, 0.0), window.time, initial=0.0)
    dissolution = cumulative_trapezoid(np.maximum(-reaction, 0.0), window.time, initial=0.0)
    eps0 = model.p.epsilon_birth if initial_epsilon is None else initial_epsilon
    eps_reaction = eps0-signed/model.p.rho_crystal
    for i, row in enumerate(records):
        row.update({
            "time_s": float(window.time[i]), "z0_m": z0,
            "birth_time_s": "" if birth_time is None else birth_time,
            "material_age_s": float(window.time[i] if birth_time is None else window.time[i]-birth_time),
            "cumulative_signed_reaction_kg_m3": float(signed[i]),
            "cumulative_positive_precipitation_kg_m3": float(positive[i]),
            "cumulative_dissolution_kg_m3": float(dissolution[i]),
            "epsilon_reaction": float(eps_reaction[i]),
            "material_balance_error": float(row["epsilon"]-eps_reaction[i]),
            "reconstruction": method,
        })
    return records


def exact_region_average(epsilon: np.ndarray, physical_faces: np.ndarray,
                         lower: float, upper: float) -> tuple[float, float]:
    overlap = np.maximum(0.0, np.minimum(physical_faces[1:], upper)
                         - np.maximum(physical_faces[:-1], lower))
    length = float(np.sum(overlap))
    return (float(np.sum(epsilon*overlap)/length), length) if length > 0 else (np.nan, 0.0)


def synthetic_profile(h0: float, velocity: float, eps_birth: float, amplitude: float,
                      length: float, birth_time: float, age: np.ndarray) -> dict[str, np.ndarray | float]:
    z0 = h0+velocity*birth_time
    distance = velocity*age
    epsilon = eps_birth-amplitude*(1.0-np.exp(-distance/length))
    reaction = 2170.0*amplitude*velocity/length*np.exp(-distance/length)
    integrated = cumulative_trapezoid(reaction, age, initial=0.0)
    return {"z0": z0, "distance": distance, "epsilon": epsilon,
            "reaction": reaction, "epsilon_reaction": eps_birth-integrated/2170.0}
