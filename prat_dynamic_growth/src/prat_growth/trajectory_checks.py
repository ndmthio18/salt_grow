"""Resolved startup and independent derivatives of integrated trajectories."""
import warnings
import numpy as np
from scipy.integrate import solve_ivp
from .dynamic import DynamicALEModel
from .verification import local_audit


def short_trajectory(model, y, end, times, max_step):
    n = model.n.n_cells
    atol = np.r_[np.full(2*n, 1e-13), 1e-15, np.full(5, 1e-13)]
    with warnings.catch_warnings(record=True) as captured:
        warnings.simplefilter("always", RuntimeWarning)
        sol = solve_ivp(model.rhs, (0, end), y, method="BDF", t_eval=times,
                        dense_output=True, rtol=1e-10, atol=atol, max_step=max_step)
    if not sol.success or not np.all(np.isfinite(sol.y)):
        raise RuntimeError(sol.message)
    return sol, sorted(set(str(w.message) for w in captured))


def startup(p, n):
    rows, states, warning_rows = [], {}, []
    times = np.r_[0., np.geomspace(1e-8, 1., 121)]
    for name, closure, initial in (("T0", "B1", "equilibrium"), ("T1", "B1", "quasi_steady"),
                                   ("T2", "B2", "equilibrium"), ("T3", "B2", "quasi_steady")):
        m = DynamicALEModel(p, n, closure)
        sol, messages = short_trajectory(m, m.initial_state(initial), 1., times, .05)
        warning_rows.append(dict(case="startup_"+name, warnings=messages))
        states[name] = sol.y
        for i, t in enumerate(times):
            eps, c, h, _ = m.unpack(sol.y[:, i])
            _, D = m.coefficients(eps)
            ct, v = m.top_state(c[-1], eps[-1], D[-1], h)
            rows.append(dict(case=name, time_s=t, velocity=v, supersaturation=ct-p.C_eq))
    return rows, dict(time=times, **states), warning_rows


def independent_time_checks(result, state_vector):
    """Short forward windows, central differences of dense trajectories.

    dy is NOT supplied by model.rhs. This complements the algebraic audit;
    varying stencil spacing reveals differentiation/cancellation error.
    """
    m = DynamicALEModel(result.physical, result.numerical, result.closure)
    rows, warning_rows = [], []
    for t in (300., 3600., 10800., 21600.):
        i = int(np.argmin(abs(result.time-t)))
        sol, messages = short_trajectory(m, state_vector(result, i), .5, [.25], .025)
        warning_rows.append(dict(case=f"time_check_{t:g}", warnings=messages))
        for dt in (.001, .01, .05):
            center = .25
            # Four-point, fourth-order central differentiation.
            dy = (sol.sol(center-2*dt)-8*sol.sol(center-dt)
                  +8*sol.sol(center+dt)-sol.sol(center+2*dt))/(12*dt)
            rows.append(dict(snapshot_time_s=t, offset_s=center, stencil_step_s=dt,
                             **local_audit(m, sol.sol(center), dy)))
    return rows, warning_rows
