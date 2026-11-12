"""Stand-alone check of the substrate pore network: no vapour solve, no crystallisation.

Regime: evaporation-wicking.  A column of pore layers is fed by the reservoir at its bottom
(fixed C0, pressure datum, saturated) and loses water at its top at a prescribed evaporation
flux.  Salt stays in the network, so the surface pore enriches until it reaches C_sat.  Two
things are checked against theory, both from Veran-Tissoires & Prat 2014 (JFM 749) section 4:

1. STEADY PROFILE.  Once the enriched layer has filled the column the profile is stationary and
   carries zero salt flux: advection up balances back-diffusion.  Continuum:

       u C = D* dC/dz   ->   C(z) = C_bottom exp(u z / D*),      Pe = u L / D*

   The network's own discretisation (upwind advection) gives C_above/C_below = 1 + Q/d per
   throat, which is the first-order form of it; refining dz must close the gap.

2. TIME OF FIRST PRECIPITATION, their equations (4.17) and (4.18):

       tau_c = Pe^-2 [ C_cris/C0 - ln(C_cris/C0) - 1 ],     tau = t D* / L^2

   valid while the enriched layer has not reached the bottom.  This is the quantity their
   figure 11 compares against a 1-D numerical solution, so it tests storage, advection and
   diffusion together rather than any single closure.

Mapping between the network and the continuum (exact, no fitted factor):

       Q  = J A_col / rho_l          liquid flow through one column [m^3/s]
       u  = Q / (eps A_col)          interstitial velocity, eps = pore volume / bulk volume
       D* = d dz / (eps A_col)       column diffusivity from the throat conductance d
       Pe = u L / D* = (Q/d) (L/dz)  -- eps cancels, so Pe is a pure network quantity

Usage::

    python -m saltv2.cases.pnm_column                      # default sweep (from the repo root)
    python -m saltv2.cases.pnm_column --layers 12 24 48

Nothing is written unless ``--csv`` is given (into the ignored simulation_outputs/).
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np

from ..materials import brine
from ..materials.brine import C_SAT
from ..pnm import PoreLayer, boundary_layer_evaporation, network_step
from .dumux1d import DELTA_BL, RH_ENV

L_COLUMN = 22.0e-3          # substrate depth of the 2014 bead column [m, M]
C0 = 0.25                   # reservoir mass fraction [-, M]
PORE_RADIUS = 270e-6        # the dgf INSCRIBED radius (R270um): the pore cube is 2 r = 540 um [L]


def column(n_layers, L=L_COLUMN, pore_radius=PORE_RADIUS):
    """One column of ``n_layers`` pores spanning depth ``L`` (the reservoir throat closes the
    remaining dz under the bottom layer)."""
    dz = L / n_layers
    if 2 * pore_radius >= dz:
        raise ValueError(f"pore 2r = {2*pore_radius*1e6:.0f} um does not fit in dz = {dz*1e6:.0f} um")
    return PoreLayer(1, n_side=1, pore_radius=pore_radius, n_layers=n_layers, dz=dz, Lt_res=dz), dz


def continuum_parameters(layer, dz, j_evap):
    """(Q, d, eps, D_star, Pe, L) for the column, in the mapping of the module docstring."""
    pitch = layer.L_box                       # one pore per layer: the column's cross-section
    A_col = pitch * pitch
    Q = j_evap * A_col / brine.density(C0)    # [m^3/s] the nominal draw, at the reservoir's
                                              # brine density (it was the constant RHO_L until
                                              # 2026-09-17; 0.93 % at C0 = 0.25)
    _, d = layer.diffusive()[0], layer.diffusive()[0]
    d_lat = layer.diffusive()[0]
    d_throat = float(d_lat[0])                # every vertical throat is identical here
    eps = float(layer.V.sum() / (A_col * layer.n_layers * dz))
    D_star = d_throat * dz / (eps * A_col)
    L = layer.n_layers * dz
    Pe = (Q / d_throat) * (L / dz)
    return Q, d_throat, eps, D_star, Pe, L


def tau_c_analytic(Pe, ratio=C_SAT / C0):
    """Veran-Tissoires & Prat 2014, equation (4.18)."""
    return (ratio - np.log(ratio) - 1.0) / Pe ** 2


def tau_c_continuum(Pe, n_cells=2000, ratio=C_SAT / C0, dtau=None):
    """tau_c from the 1-D continuum problem itself (their equations 4.4-4.5), not from the
    closed form: implicit finite volume in z' = z/L, tau = t D*/L^2, salt flux zero at the top,
    C' = 1 at the bottom.  This is the curve their figure 11 calls "numerical computation"."""
    from scipy.sparse import diags
    from scipy.sparse.linalg import spsolve
    h = 1.0 / n_cells
    dtau = dtau if dtau is not None else 0.2 * h / max(Pe, 1e-12) * h   # advective + diffusive
    C = np.ones(n_cells)                       # C' = C / C0, bottom cell index 0, top n_cells-1
    lower = np.zeros(n_cells - 1); upper = np.zeros(n_cells - 1); diag = np.ones(n_cells)
    # interior faces: upwind advection (flow is upwards, from cell i to i+1) + diffusion
    for i in range(n_cells - 1):
        adv, dif = Pe, 1.0 / h
        diag[i] += dtau / h * (adv + dif); upper[i] = -dtau / h * dif
        diag[i + 1] += dtau / h * dif;     lower[i] = -dtau / h * (adv + dif)
    # bottom face: reservoir at C' = 1, same upwind form; top face: zero salt flux
    diag[0] += dtau / h * (1.0 / (0.5 * h))
    rhs0 = dtau / h * (Pe + 1.0 / (0.5 * h))
    A = diags([lower, diag, upper], [-1, 0, 1], format="csr")
    tau = 0.0
    while C[-1] < ratio:
        b_ = C.copy(); b_[0] += rhs0
        C = spsolve(A, b_)
        tau += dtau
        if tau > 1e4:
            return np.nan
    return tau


def run_column(n_layers, j_evap, dt_factor=0.05, t_max=None, steady=False,
               top_bc="prescribed", delta_bl=DELTA_BL):
    """Transient enrichment of the column.  Returns (t_c, C_profile_at_t_c, params).

    ``steady``: keep going after C_sat (no precipitation here) until the profile stops moving,
    and return that profile instead.

    ``top_bc``:
      ``"prescribed"``       the flux ``j_evap`` is imposed for the whole run, whatever the
          surface pore's concentration does.  This is the case the 2014 closed form (4.18)
          assumes -- constant Pe -- so it is what the comparison in README section 13.1 needs.
      ``"boundary_layer"``   the DuMuX reference case's condition instead
          (``saltv2.pnm.boundary_layer_evaporation``, README section 13.6): the surface pore
          evaporates through a 1.6 mm boundary layer driven by its OWN a_w, over the pore
          opening (2 r)^2 rather than the whole column cross-section.  ``j_evap`` is then only
          the value the reported Pe is computed at (t = 0), not a boundary condition: as the
          pore enriches a_w falls, the draw falls with it and Pe drifts down, so (4.18) does not
          apply.  That drift is the point of the option.

    ``delta_bl`` [m] is the boundary-layer thickness of that condition.  MIND THE TWO DELTAS:
    the default 1.6 mm is the PORE-SCALE layer of QSPNM and of the DuMuX reference case, while
    the 2014 experiments are matched by delta_eff = 18 mm, which stands for the whole recessed
    tube and its free convection (README section 13.1 of runs/exp_2014, plan/14 section 12).
    They differ by 11x, so a column driven at 1.6 mm evaporates about 5x faster than the 2014
    plane flux and precipitates in minutes rather than hours.  Neither is wrong; they are
    boundary conditions of different experiments.
    """
    if top_bc not in ("prescribed", "boundary_layer"):
        raise ValueError("top_bc must be 'prescribed' or 'boundary_layer'")
    layer, dz = column(n_layers)
    Q, d_throat, eps, D_star, Pe, L = continuum_parameters(layer, dz, j_evap)
    g_lat, g_res = layer.conductances()
    d_lat, d_res = layer.diffusive()
    top = int(layer.surface[0])
    A_open = 4.0 * PORE_RADIUS ** 2                # their poreCrossSectionalArea = (2 r)^2
    # Salt inventory through the same brine density network_step uses; a constant RHO_L here would
    # disagree with it by ~1 % at C0 = 0.25 (it did until 2026-09-17, when salt_step changed).
    m = brine.salt_mass(C0, layer.V)
    M = j_evap * layer.L_box ** 2                  # water leaving the column [kg/s]
    # step: a fraction of the fastest pore turnover (advective + diffusive exchange)
    dt = dt_factor * float(np.min(layer.V / (Q + d_throat + d_res.max())))
    t = 0.0
    t_c = np.nan
    t_max = t_max if t_max is not None else 1e9
    while t < t_max:
        if top_bc == "boundary_layer":
            w_top = float(brine.mass_fraction_from_salt_mass(m[top], layer.V[top]))
            M = float(boundary_layer_evaporation(w_top, A_open, delta_bl, RH_ENV))
        # liquid and salt of the step together, pore storage included (pnm.network_step); the
        # surface pore loses the liquid mass M
        q_mass = np.zeros(layer.N); q_mass[top] = M
        m = network_step(layer, m, C0, q_mass, g_lat, g_res, d_lat, d_res, dt)[0]
        t += dt
        C = brine.mass_fraction_from_salt_mass(m, layer.V)
        if not np.isfinite(t_c) and C[top] >= C_SAT:
            t_c = t
            if not steady:
                break
        if steady and np.isfinite(t_c) and t > 50 * t_c:
            break
    return dict(t_c=t_c, C=brine.mass_fraction_from_salt_mass(m, layer.V), layer=layer, dz=dz,
                Q=Q, d=d_throat, eps=eps, D_star=D_star, Pe=Pe, L=L, dt=dt, t=t, top_bc=top_bc,
                M_top=M)          # [kg/s] the water rate actually drawn at the end of the run


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--layers", type=int, nargs="+", default=[11, 22, 44, 88],
                    help="pore layers over the 22 mm column")
    ap.add_argument("--flux", type=float, nargs="+",
                    default=[1.5854e-5, 5e-6, 1.5e-6],
                    help="evaporation flux [kg m^-2 s^-1]; the first is the 2014 Table 1 value")
    ap.add_argument("--dt-factor", type=float, default=0.05)
    ap.add_argument("--csv", type=Path, default=None)
    args = ap.parse_args(argv)

    rows = []
    print(f"column {L_COLUMN*1e3:.0f} mm, pore 2r = {2*PORE_RADIUS*1e6:.0f} um, C0 = {C0}, "
          f"C_sat = {C_SAT}\n")
    print(f"{'layers':>6} {'dz [um]':>8} {'eps':>6} {'D* [m2/s]':>10} {'Pe':>8} "
          f"{'t_c PNM [h]':>12} {'1-D same dz':>13} {'1-D resolved':>13} {'eq 4.18 [h]':>13} "
          f"{'PNM vs 1-D':>10}")
    for j_evap in args.flux:
        for n in args.layers:
            r = run_column(n, j_evap, dt_factor=args.dt_factor)
            scale = r["L"] ** 2 / r["D_star"]
            t_ana = tau_c_analytic(r["Pe"]) * scale
            t_fine = tau_c_continuum(r["Pe"]) * scale               # resolved continuum
            t_same = tau_c_continuum(r["Pe"], n_cells=n) * scale    # continuum at the SAME dz
            err = (r["t_c"] - t_same) / t_same
            rows.append((j_evap, n, r["dz"], r["eps"], r["D_star"], r["Pe"], r["t_c"], t_same,
                         t_fine, t_ana, err))
            print(f"{n:6d} {r['dz']*1e6:8.1f} {r['eps']:6.3f} {r['D_star']:10.3e} {r['Pe']:8.3f} "
                  f"{r['t_c']/3600:12.4f} {t_same/3600:13.4f} {t_fine/3600:13.4f} "
                  f"{t_ana/3600:13.4f} {err:+8.2%}")
        print()
    # The steady zero-flux profile is the ceiling of the enrichment: C_surface / C_bottom is
    # exp(Pe) in the continuum and (1 + Pe/n)^n in the network's upwind form.  At the 2014 flux
    # both are far above saturation, i.e. the column precipitates long before it gets there --
    # which is why the quantity to compare is the onset time above, not the steady profile.
    # (tests/test_pnm_column.py checks the steady balance itself, at a Pe where it stays a brine.)
    n_nat = int(round(L_COLUMN / 0.7e-3))
    layer, dz = column(n_nat)
    _, _, _, _, Pe, _ = continuum_parameters(layer, dz, args.flux[0])
    ratio_discrete = (1.0 + Pe / n_nat) ** n_nat
    print(f"steady ceiling, {n_nat} layers (dz = pitch = {dz*1e6:.0f} um), Pe = {Pe:.2f}: "
          f"C_surface / C_bottom = {ratio_discrete:.1f} (network) vs {np.exp(Pe):.1f} (continuum); "
          f"saturation is reached at {C_SAT / C0:.3f}, so neither is ever attained")

    if args.csv:
        args.csv.parent.mkdir(parents=True, exist_ok=True)
        np.savetxt(args.csv, np.array(rows), delimiter=",",
                   header="j_evap,n_layers,dz,eps,D_star,Pe,t_c_pnm,t_c_1d_same_dz,t_c_1d_resolved,t_c_eq418,rel_err")
        print(f"written {args.csv}")
    return rows


if __name__ == "__main__":
    main()
