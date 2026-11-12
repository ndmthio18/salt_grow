"""Independent 1-D conservative numerical reference for the Prat (2026) growth problem with
layer-dependent porosity, specific surface and effective diffusivity.

What is solved (Prat 2026, equation numbers of the paper; derivation checked against
../salt_grow/reference/Marc Prat 2026 ... .md sections 2-3 and 7):

    (1)  d(rho eps C)/dt + d(rho C U)/dz = d(rho eps D_eff dC/dz)/dz - a_v rho k_r (C - C_eq)
    (2)  -rho_cr d(eps)/dt = a_v rho k_r (C - C_eq)
    (3)  d(rho eps)/dt + d(rho U)/dz = -a_v rho k_r (C - C_eq)          (rho constant, eq 4)
    (6)  rho U = [rho eps + rho_cr (1 - eps)] dh/dt + j                 at z = h   (rho_v dropped)
    (7)  rho C U - rho eps D_eff dC/dz = [rho C eps + rho_cr (1 - eps)] dh/dt   at z = h
    (20) dh/dt = beta (rho / rho_cr) k_r (C_top - C_eq)                 interface kinetics
         C(0) = C_eq                                                    substrate side

NOT used: eq (29) F(Da*) or any other analytic growth law.  The interface speed follows from
(20) with C_top solved from the discretised (1), (3), (6), (7).  Equations (29), (31), (19) are
evaluated separately (``prat_analytic``) only to compare with.

Frame.  xi = h(t) - z, distance below the moving interface, on a window 0 <= xi <= L with
L >> delta_C (the material below L is at C = C_eq to exponential accuracy; eq 19).  In this frame
the solid moves at +v = dh/dt and the liquid Darcy flux relative to the frame is W = U - eps v
(upward, towards xi = 0).  Exact transformations of (1)-(3):

    d(rho eps C)/dt|xi + d/dxi [ -rho C W - rho eps D dC/dxi ] = -S,    S = a_v rho k_r (C - C_eq)
    d(eps)/dt|xi + v d(eps)/dxi = -S / rho_cr
    d(eps)/dt|xi - dW/dxi      = -S / rho
    interface (xi = 0):  rho W_i C_top + rho eps_i D_i dC/dxi|0 = rho_cr (1 - eps_i) v = P_i,
                         W_i = U_evap + (rho_cr / rho)(1 - eps_i) v,   v = beta (rho/rho_cr) k_r (C_top - C_eq)
    bottom (xi = L):     C = C_eq

The material at xi sits at the lab height z_m = h(t) - xi, which decides its layer (crystal size
d_c(n), birth porosity eps_new(n)); the specific surface and the diffusivity are Prat's own
closures evaluated with the LOCAL porosity: a_v = 5 (1 - eps)/d_c (eq 16), D_eff = sqrt(eps) D_s.
Porosity is a material field advected with the solid (upwind) and lowered by (2).

Discretisation.  Finite volumes on a geometrically stretched grid (finest at the interface),
backward Euler, upwind advection, Picard iteration on (v, eps) inside each step until the step's
salt and solid balances close to ``picard_tol``.  Each step reports both balances, so the run
states a residual instead of asserting conservation.

Everything is SI.  Parameter provenance: ``PratParams``.
"""
from __future__ import annotations

from dataclasses import dataclass, field, replace

import numpy as np
from scipy.linalg import solve_banded

DAY = 86400.0


# ----------------------------------------------------------------------------- parameters
@dataclass(frozen=True)
class PratParams:
    """Prat (2026) base values; see prat2026/params.py for provenance.  [L] literature, [B] back-calibrated."""
    C_eq: float = 0.264          # NaCl solubility, mass fraction [L, paper Materials and Methods]
    rho: float = 1200.0          # saturated brine density [kg/m^3] [B, prat2026/params.py]
    rho_cr: float = 2170.0       # halite [kg/m^3] [L, below eq 2]
    beta: float = 2.0            # geometric up-scaling factor of eq 20 [L, Appendix]
    k_r: float = 2.3e-3          # precipitation coefficient at 21 C [m/s] [L, Naillon]
    D_s: float = 1.3e-9          # NaCl diffusivity [m^2/s] [L, paper text]; the project's PNM uses 1.54e-9
    dc_over_dsub: float = 0.14   # crystal / substrate grain [L, Wijnhorst]

    def D_eff(self, eps):
        """D_s,eff = D_s / tau, tau = eps^-0.5 (paper, after eq 17)."""
        return np.sqrt(eps) * self.D_s

    def a_v(self, eps, d_c):
        """Specific surface of a cubic-crystal packing, eq 16: 5 (1 - eps)/d_c."""
        return 5.0 * (1.0 - eps) / d_c


@dataclass(frozen=True)
class Layering:
    """Wijnhorst-type stack as the project writes it (saltv2/growth/crust.py): crystals shrink by
    ``r`` and layer thickness grows by ``alpha`` per layer; h_1 = c_h d_sub.  r = 1 is one scale."""
    d_sub: float
    r: float = 0.5
    alpha: float = 2.0
    c_h: float = 1.0
    n_max: int = 40

    def layer(self, z_m):
        """Layer index n >= 1 of material at lab height z_m [m] (same convention as crust.generation)."""
        z = np.maximum(np.asarray(z_m, float), 0.0)
        h1 = self.c_h * self.d_sub
        if self.alpha == 1.0:
            n = np.floor(z / h1) + 1.0
        else:
            n = np.floor(np.log(1.0 + z * (self.alpha - 1.0) / h1) / np.log(self.alpha)) + 1.0
        return np.clip(n, 1, self.n_max).astype(int)

    def d_c(self, n, p: PratParams):
        return p.dc_over_dsub * self.d_sub * self.r ** (np.asarray(n, float) - 1.0)

    def layer_bottom(self, n):
        """Lab height at which layer n starts."""
        h1 = self.c_h * self.d_sub
        if self.alpha == 1.0:
            return h1 * (n - 1)
        return h1 * (self.alpha ** (n - 1) - 1.0) / (self.alpha - 1.0)


# ----------------------------------------------------------------------------- analytic Prat
def prat_analytic(d_c, eps, U_evap, p: PratParams = PratParams()):
    """Prat (2026) closed forms for ONE crystal size and ONE porosity: Da* (27), F (29),
    delta_C (after 19), C*_top (31), dh/dt (32), U/U_evap (33).  Da >> 1 assumed as in the paper."""
    av = p.a_v(eps, d_c)
    De = p.D_eff(eps)
    Da_star = p.beta ** 2 * (1 - eps) ** 2 * (1 - p.C_eq) * p.k_r / (eps * De * av)
    F = np.sqrt(Da_star) / (1 + np.sqrt(Da_star))
    delta_C = np.sqrt(eps * De / ((1 - p.C_eq) * av * p.k_r))
    C_top_star = U_evap * p.C_eq / (p.beta * (1 - eps) * (1 - p.C_eq) * p.k_r
                                    + np.sqrt((1 - p.C_eq) * av * p.k_r * eps * De))
    v = p.rho * p.C_eq * U_evap / (p.rho_cr * (1 - eps) * (1 - p.C_eq)) * F
    U_ratio = 1 + (eps + p.rho_cr / p.rho * (1 - eps)) * v / U_evap
    return dict(Da_star=float(Da_star), F=float(F), delta_C=float(delta_C), C_top_star=float(C_top_star),
                v=float(v), U_over_Uevap=float(U_ratio), a_v=float(av), D_eff=float(De))


def prat_profile(xi, delta_C, C_top_star):
    """Eq 19 in the Da >> 1 limit: C* = C*_top exp(-xi / delta_C)."""
    return C_top_star * np.exp(-np.asarray(xi, float) / delta_C)


# ----------------------------------------------------------------------------- solver
@dataclass
class ReferenceSolver:
    """Time-dependent 1-D reference solution in the interface frame."""
    layering: Layering
    U_evap: float                       # evaporation velocity j / rho [m/s]
    eps_new: dict | float               # birth porosity per layer {n: eps} or one value
    p: PratParams = PratParams()
    L: float | None = None              # window depth [m]; default 60 delta_C of layer 1
    dxi0: float | None = None           # first cell [m]; default delta_C(layer n_ref) / n_band
    n_band: int = 40                    # cells per delta_C of the reference layer
    n_ref_layer: int = 1                # the layer whose delta_C sets the first cell
    stretch: float = 1.03               # geometric growth of the cell size
    evolve_eps: bool = True             # eq (2); False keeps eps at its birth value (Prat's constant-eps model)
    liquid_consumption: bool = True     # eq (3)/(4): W varies through the band; False: W = W_i everywhere
    local_props: bool = True            # a_v = 5(1-eps)/d_c and D_eff = sqrt(eps) D_s with the LOCAL eps (Prat's
                                        # closures evaluated where the material is); False: both at the birth eps,
                                        # i.e. the constant-property band of Prat's analytic solution
    transport_eps_fixed: bool = False   # eps in the accumulation and diffusion terms of (1) held at the birth value
                                        # (Prat's linearised band: constant eps_s in eq 5) while (2)-(3) still shrink
                                        # the pores for the liquid balance -- the hybrid the analytic law is derived from
    picard_tol: float = 1e-10
    picard_max: int = 30
    h0: float = 0.0                     # initial height [m] (the material at z_m < h0 is given its layer's eps_new)
    eps_stop: float = 1e-3              # a band porosity at or below this is a closed pore: the run stops and
                                        # records ``closure_event`` (no clamping, no continuation)

    # state (set by setup)
    xi: np.ndarray = field(init=False, repr=False)
    dxi: np.ndarray = field(init=False, repr=False)
    xf: np.ndarray = field(init=False, repr=False)
    C: np.ndarray = field(init=False, repr=False)
    eps: np.ndarray = field(init=False, repr=False)
    h: float = field(init=False)
    t: float = field(init=False)
    v: float = field(init=False)
    C_top: float = field(init=False)

    def __post_init__(self):
        p = self.p
        lay = self.layering
        e1 = self._eps_new(1)
        ref = prat_analytic(lay.d_c(self.n_ref_layer, p), self._eps_new(self.n_ref_layer), self.U_evap, p)
        a1 = prat_analytic(lay.d_c(1, p), e1, self.U_evap, p)
        if self.L is None:
            self.L = 60.0 * a1["delta_C"]
        if self.dxi0 is None:
            self.dxi0 = ref["delta_C"] / self.n_band
        # geometric grid: dxi_k = dxi0 * stretch^k until the window is covered
        widths = [self.dxi0]
        while sum(widths) < self.L:
            widths.append(widths[-1] * self.stretch)
        widths[-1] -= sum(widths) - self.L                     # trim the last cell onto L exactly
        self.dxi = np.array(widths)
        self.xf = np.concatenate([[0.0], np.cumsum(self.dxi)])  # faces, xf[0] = 0 (interface), xf[-1] = L
        self.xi = 0.5 * (self.xf[:-1] + self.xf[1:])
        self.h = float(self.h0); self.t = 0.0
        n = lay.layer(self.z_m())
        self.eps = np.array([self._eps_new(k) for k in n], float)
        self.C = np.full(self.xi.size, p.C_eq)
        self.v = 0.0; self.C_top = p.C_eq
        self.history = []
        self.closure_event = None
        self.cum = dict(salt_in=0.0, salt_out=0.0, salt_acc=0.0, solid_in=0.0, solid_out=0.0, solid_src=0.0, solid_acc=0.0)

    # ---- helpers
    def _eps_new(self, n):
        if isinstance(self.eps_new, dict):
            return float(self.eps_new[int(n)]) if int(n) in self.eps_new else float(self.eps_new[max(self.eps_new)])
        return float(self.eps_new)

    def z_m(self):
        return self.h - self.xi

    def props(self, eps):
        """(a_v, D_eff, d_c, layer) of the material now at every cell."""
        n = self.layering.layer(self.z_m())
        d_c = self.layering.d_c(n, self.p)
        eps_birth = np.array([self._eps_new(k) for k in n], float)
        e = np.maximum(eps, 0.0) if self.local_props else eps_birth
        return self.p.a_v(e, d_c), self.p.D_eff(e), d_c, n, eps_birth

    def interface_eps(self):
        return self._eps_new(int(self.layering.layer(self.h)))

    # ---- one step
    def step(self, dt):
        p = self.p; N = self.xi.size
        rho, rho_cr, C_eq = p.rho, p.rho_cr, p.C_eq
        eps_n, C_n = self.eps.copy(), self.C.copy()
        eps_i = self.interface_eps()                           # birth porosity of the material being laid down
        v = self.v; eps_new = eps_n.copy(); C_new = C_n.copy()
        av, D, d_c, n_layer, eps_birth = self.props(eps_n)
        dist = np.diff(self.xi)                                # centre-to-centre distances
        for it in range(self.picard_max):
            # --- porosity (eq 2 + solid advection, implicit upwind, sink from the latest C):
            #     (1 + a_i) eps_i - a_i eps_{i-1} = eps_i^n - dt S_i / rho_cr,  eps_{-1} = eps_i (interface)
            S = av * rho * p.k_r * (C_new - C_eq)
            if self.evolve_eps:
                a_adv = v * dt / self.dxi
                ab_e = np.zeros((2, N)); ab_e[0] = 1.0 + a_adv; ab_e[1, :-1] = -a_adv[1:]
                rhs_e = eps_n - dt * S / rho_cr; rhs_e[0] += a_adv[0] * eps_i
                eps_new = solve_banded((1, 0), ab_e, rhs_e)
            else:
                eps_new = np.array([self._eps_new(k) for k in n_layer], float)   # Prat's constant-eps model
            if np.any(eps_new <= self.eps_stop):
                # a pore closed during this step: the step is abandoned, the event recorded, the run stops
                # ``layer`` is the layer AT THE INTERFACE when the band closes (the material whose arrival
                # triggers it); ``layer_material`` is the layer of the cell that closed, which can be the
                # previous layer's last band thickness just below the boundary
                i_c = int(np.argmax(eps_new <= self.eps_stop))
                self.closure_event = dict(t=self.t + dt, h=self.h + v * dt, layer=int(self.layering.layer(self.h + v * dt)),
                                          layer_material=int(n_layer[i_c]), xi=float(self.xi[i_c]),
                                          eps_min=float(eps_new.min()), picard=it + 1)
                return None
            av, D, d_c, n_layer, eps_birth = self.props(eps_new)
            # --- liquid flux relative to the frame (eq 3/4 integrated from the interface, discrete-consistent):
            #     W_{i+1/2} = W_{i-1/2} - v (eps_i - eps_{i-1}) + dxi_i S_i (1/rho - 1/rho_cr)
            W_i = self.U_evap + rho_cr / rho * (1.0 - eps_i) * v
            W = np.empty(N + 1); W[0] = W_i
            if self.liquid_consumption:
                W[1:] = W_i - v * (eps_new - eps_i) + np.cumsum(self.dxi * S * (1.0 / rho - 1.0 / rho_cr))
            else:
                W[:] = W_i
            # --- salt (eq 1 in the frame), backward Euler, upwind towards the interface, tridiagonal
            eps_tr = eps_birth if self.transport_eps_fixed else eps_new      # eps of the transport operator
            eps_tr_n = eps_birth if self.transport_eps_fixed else eps_n
            epsD = eps_tr * D
            epsD_f = 0.5 * (epsD[:-1] + epsD[1:])              # interior faces 1..N-1
            # interface closure (eqs 7 + 20):  a C_top + b (C_0 - C_top) = c (C_top - C_eq)
            a_ = rho * W[0]; b_ = rho * eps_i * p.D_eff(eps_i) / self.xi[0]; c_ = rho * p.beta * (1 - eps_i) * p.k_r
            den = b_ + c_ - a_
            kappa = c_ * b_ / den                               # P_i = kappa (C_0 - C_eq) + P_conv
            P_conv = c_ * a_ * C_eq / den
            acc = rho * eps_tr * self.dxi / dt
            sink = av * rho * p.k_r * self.dxi
            diag = acc + sink
            rhs = rho * eps_tr_n * C_n * self.dxi / dt + sink * C_eq
            lower = np.zeros(N); upper = np.zeros(N)
            # interior face f between cell i (above) and i+1 (below): flux in +xi = -rho W_f C_{i+1} - Df (C_{i+1} - C_i)
            Wf = rho * W[1:N]; Df = rho * epsD_f / dist
            diag[:-1] += Df; upper[:-1] += -(Wf + Df)           # cell i gains -F_f
            diag[1:] += Wf + Df; lower[1:] += -Df               # cell i+1 loses -F_f
            # interface: cell 0 loses P_i = kappa (C_0 - C_eq) + P_conv
            diag[0] += kappa; rhs[0] += kappa * C_eq - P_conv
            # bottom: Dirichlet C = C_eq at xi = L
            DL = rho * epsD[-1] / (self.L - self.xi[-1])
            diag[-1] += DL; rhs[-1] += rho * W[-1] * C_eq + DL * C_eq
            ab = np.zeros((3, N)); ab[0, 1:] = upper[:-1]; ab[1] = diag; ab[2, :-1] = lower[1:]
            C_next = solve_banded((1, 1), ab, rhs)
            C_top = (b_ * C_next[0] + c_ * C_eq) / den
            v_next = p.beta * rho / rho_cr * p.k_r * (C_top - C_eq)
            converged = (abs(v_next - v) <= self.picard_tol * max(abs(v_next), 1e-30)
                         and np.allclose(C_next, C_new, rtol=1e-11, atol=0.0))
            v, C_new = v_next, C_next
            if converged and it > 0:
                break
        # --- balances of the step (kg/m^2)
        S = av * rho * p.k_r * (C_new - C_eq)
        P_i = c_ * (C_top - C_eq)
        salt_in = (rho * W[-1] * C_eq + DL * (C_eq - C_new[-1])) * dt                     # through xi = L
        salt_acc = float(np.sum(rho * (eps_tr * C_new - eps_tr_n * C_n) * self.dxi))
        salt_out = (P_i + float(np.sum(S * self.dxi))) * dt
        res_salt = (salt_in - salt_acc - salt_out) / max(abs(salt_in), 1e-300)
        solid_in = rho_cr * (1 - eps_i) * v * dt                                            # material born at the interface
        solid_out = rho_cr * (1 - eps_new[-1]) * v * dt                                     # material leaving at L
        solid_acc = float(np.sum(rho_cr * (eps_n - eps_new) * self.dxi))                   # pores filled inside
        solid_src = float(np.sum(S * self.dxi)) * dt if self.evolve_eps else 0.0
        # discrete solid balance of the upwind scheme: acc = in - out + src (exact up to round-off)
        res_solid = (solid_in - solid_out + solid_src - solid_acc) / max(abs(solid_in), 1e-300)
        for k_, v_ in dict(salt_in=salt_in, salt_out=salt_out, salt_acc=salt_acc, solid_in=solid_in, solid_out=solid_out,
                           solid_src=solid_src, solid_acc=solid_acc).items():
            self.cum[k_] += v_
        self.C, self.eps, self.v, self.C_top = C_new, eps_new, v, C_top
        if np.any(eps_new <= self.eps_stop):
            i_c = int(np.argmax(eps_new <= self.eps_stop))
            self.closure_event = dict(t=self.t, h=self.h, layer=int(self.layering.layer(self.h)),
                                      layer_material=int(n_layer[i_c]), xi=float(self.xi[i_c]), eps_min=float(eps_new.min()))
        self.h += v * dt; self.t += dt
        S_int = float(np.sum(S * self.dxi))
        F_num = P_i / (P_i + S_int) if (P_i + S_int) > 0 else np.nan        # interface share of all precipitation
        F_prat = P_i * (1 - C_eq) / (C_eq * rho * self.U_evap)               # Prat's F: v / v_max of eq 29
        rec = dict(t=self.t, h=self.h, v=v, C_top_star=C_top - C_eq, P_i=P_i, S_int=float(np.sum(S * self.dxi)),
                   F_num=F_num, F_prat=F_prat, eps_final=float(eps_new[-1]), eps_interface=eps_i,
                   layer_top=int(self.layering.layer(self.h)), res_salt=res_salt, res_solid=res_solid,
                   picard=it + 1, W_i=W[0], W_L=W[-1])
        self.history.append(rec)
        return rec

    def run(self, t_end, dt, h_end=None, record_every=1):
        """Advance to t_end (or until h >= h_end, or until a pore closes).  ``record_every`` thins the history."""
        k = 0
        while self.t < t_end - 1e-12 and (h_end is None or self.h < h_end):
            rec = self.step(min(dt, t_end - self.t)); k += 1
            if rec is None:                                    # pore closure: the step was abandoned
                break
            if record_every > 1 and k % record_every:
                self.history.pop()
            if self.closure_event is not None:
                break
        return self.history

    def cumulative_balance(self):
        """Run-integrated salt and solid balances (kg/m^2) and their relative residuals."""
        c = self.cum
        return dict(c, res_salt=(c['salt_in'] - c['salt_acc'] - c['salt_out']) / max(abs(c['salt_in']), 1e-300),
                    res_solid=(c['solid_in'] - c['solid_out'] + c['solid_src'] - c['solid_acc']) / max(abs(c['solid_in']), 1e-300))

    def run_to_steady(self, dt, tol=1e-9, t_max=None, min_transits=3):
        """Advance until v and eps_final stop changing (relative change per band-transit time below tol)."""
        a1 = prat_analytic(self.layering.d_c(1, self.p), self._eps_new(1), self.U_evap, self.p)
        transit = self.L / max(a1["v"], 1e-30)                 # time for material to cross the window
        t_max = 20 * transit if t_max is None else t_max
        last = None
        while self.t < t_max:
            for _ in range(max(1, int(round(transit / dt)))):
                if self.step(dt) is None:
                    raise RuntimeError(f"pore closure in a homogeneous run: {self.closure_event}")
            cur = (self.v, self.eps[-1])
            if last is not None and self.t > min_transits * transit:
                if abs(cur[0] - last[0]) <= tol * abs(cur[0]) and abs(cur[1] - last[1]) <= tol * max(abs(cur[1]), 1e-30):
                    return self.history[-1]
            last = cur
        return self.history[-1]


def homogeneous_steady(d_sub, eps_new, U_evap, p=PratParams(), evolve_eps=True, liquid_consumption=True,
                       local_props=True, transport_eps_fixed=False, n_band=40, stretch=1.03, L_over_deltaC=60.0, dt_over_transit=0.5, tol=1e-9):
    """Steady travelling-wave solution of ONE material (r = 1): returns the last record plus the profile."""
    lay = Layering(d_sub=d_sub, r=1.0, alpha=1.0, c_h=1.0)
    a = prat_analytic(lay.d_c(1, p), eps_new, U_evap, p)
    s = ReferenceSolver(lay, U_evap, eps_new, p, L=L_over_deltaC * a["delta_C"], n_band=n_band, stretch=stretch,
                        evolve_eps=evolve_eps, liquid_consumption=liquid_consumption, local_props=local_props,
                        transport_eps_fixed=transport_eps_fixed)
    dt = dt_over_transit * a["delta_C"] / a["v"]
    rec = s.run_to_steady(dt, tol=tol)
    rec = dict(rec); rec.update(analytic=a, xi=s.xi.copy(), C=s.C.copy(), eps=s.eps.copy(), N=int(s.xi.size), dt=dt)
    return rec


def eps_new_for_final(d_sub, eps_final_target, U_evap, p=PratParams(), xtol=1e-7,
                      scan=(0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.99), **kw):
    """eps_final interpretation: the birth porosity whose steady final porosity equals the target.

    A birth porosity whose band closes (eps -> 0) counts as eps_final = 0, so the bracket function
    g(eps_new) = eps_final - target is continuous.  g need not be monotone in the consistent model
    (the local a_v = 5(1-eps)/d_c feedback closes the band again at very open births), so the
    bracket is found by a coarse scan upward from the target and Brent's method is run inside the
    first sign change.  Returns (None, hist, None) when no birth porosity in [target, 0.99] reaches
    the target: the eps_final interpretation is then INFEASIBLE for this material."""
    from scipy.optimize import brentq
    hist = []; cache = {}

    def g(e):
        if e in cache:
            return cache[e]
        try:
            r = homogeneous_steady(d_sub, e, U_evap, p, **kw)
            val = r["eps_final"] - eps_final_target
            hist.append((e, r["eps_final"], r["F_num"], r))
        except RuntimeError:                                    # pore closure: the material ends at eps = 0
            val = 0.0 - eps_final_target
            hist.append((e, 0.0, np.nan, None))
        cache[e] = val
        return val

    lo = eps_final_target; f_lo = g(lo)
    if f_lo >= 0:
        return lo, [h[:3] for h in hist], hist[-1][3]
    for hi in [x for x in scan if x > lo]:
        f_hi = g(hi)
        if f_hi >= 0:
            e = brentq(g, lo, hi, xtol=xtol)
            r = homogeneous_steady(d_sub, e, U_evap, p, **kw)
            return e, [h[:3] for h in hist], r
        lo, f_lo = hi, f_hi
    return None, [h[:3] for h in hist], None
