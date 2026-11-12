"""The substrate pore network: geometry, throat conductances.  Regime: evaporation-wicking.

``PoreLayer`` is either the generated lattice (the 5x5 QSPNM dgf scaled by lambda, optionally
stacked into a column) or, through ``PoreLayer.from_dgf``, a DuMuX network read as it is -- the
one the transport and the pressure solve were validated against.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np

from ..materials.brine import D_NACL, MU
from .dgf import parse_dgf_2d

GRIDS = Path(__file__).resolve().parent / "grids"
# verbatim copy of QSPNM_enhanced/.trees/isolated_clusters/QSPNM_enhanced/grids_syn/<same name>
# (copied 2026-09-03 so the model is self-contained on other machines)
DGF = str(GRIDS / "N5x5_R270um_SD10um_L3.5mmx3.5mm_n1_top.dgf")
# verbatim copy of dumux/pnm2pncmin2/appl/grids/pnm_1d_homogeneous.dgf (byte-identical to the
# QSPNM test grid of the same name): 6 cubic pores 0.6 mm apart, the bottom one the reservoir
DGF_DUMUX_1D = str(GRIDS / "pnm_1d_homogeneous.dgf")
DGF_DUMUX_5X5_VARTHROAT = str(
    GRIDS / "N5x5_R270um_SD10um_L3.5mmx3.5mm_n1_top_varthroat.dgf"
)
K_SQUARE = 0.5623            # Patzek & Silin prefactor of a square duct [-, L] (DuMuX)
G_SQUARE = 1.0 / 16.0        # shape factor of a square [-, L]


class PoreLayer:
    """Surface pores on a periodic square box, lateral throats, one reservoir
    throat per pore. lam scales r -> r/lam and multiplies the count by lam^2.

    ``n_layers > 1`` stacks that pattern downwards into a column of pore layers
    joined by vertical throats; only the BOTTOM layer keeps a reservoir throat,
    so the substrate depth between the surface and the reservoir is resolved
    instead of being lumped into one long throat.  Pores are numbered layer by
    layer from the surface, so indices 0 .. n_side^2 - 1 are the surface pores
    (``surface``) that the drivers open into the gas domain, and every array is
    over ALL pores.  Regime: evaporation-wicking -- the network stays
    liquid-saturated, the reservoir holds the bottom layer at C0 and at the
    pressure datum; saturation is not a variable (that is the drying stage)."""

    def __init__(self, lam=1, L_box=3.5e-3, seed=0, Lt_res=0.7e-3, n_side=None, cv=0.0,
                 pitch_factor=1.0, pore_radius=None, n_layers=1, dz=None, throat_shape="square"):
        """n_side=None: 5*lam pores per side in the fixed L_box (total plane area
        conserved).  n_side=k: k pores per side at the same pitch 0.7 mm/lam, so
        the periodic box shrinks to k*0.7 mm/lam (small box for large lam).
        cv > 0: pore radii lognormal with this coefficient of variation around the
        dgf median (270 um / lam), rescaled so the total opening area equals the
        cv = 0 value, and clipped to 2r <= 0.95 pitch (no overlap on the lattice).
        ``pore_radius`` [m], when supplied, assigns that same physical radius
        to every pore.  It is intended for controlled one-/few-pore studies;
        it cannot be combined with ``cv > 0``.  Throat radii then scale with
        the pores (lateral: 0.485 * min(r_a, r_b);
        reservoir: 0.485 * r_i; 0.485 = 131/270 of the dgf).  pitch_factor > 1
        widens the lattice pitch (lower opening fraction) to admit wider cv.
        ``throat_shape`` "square" (default since 2026-09-17) or "circle": the pores are QSPNM
        CUBES, so a square throat is the consistent cross-section and it is what the DuMuX
        reference case (pnm2pncmin2) uses.  Square: area A = 4 r^2, shape factor
        G = 1/16, conductance 0.5623 A^2 G / (mu L) (Patzek & Silin, DuMuX
        transmissibility1p.hh), diffusive conductance D A / L.  Circle: pi r^2, pi r^4 / (8 mu L),
        D pi r^2 / L -- what every run before that date used.
        ``n_layers`` k > 1 stacks k copies of the layer downwards, ``dz`` apart
        (default: the lateral pitch), joined by vertical throats built with the
        same template as the lateral ones.  The reservoir throat of length
        ``Lt_res`` then hangs under the BOTTOM layer only, so the depth it
        stands for is the remaining one: the column itself resolves
        (k - 1) * dz.  Gravity is not in the network (unchanged from the
        single-layer model): p is the viscous pressure with the reservoir at 0."""
        V, E = parse_dgf_2d(DGF)
        base_r = V[:25, 2]                               # 270 +- 10 um
        base_rt = E[:40, 2].mean()                       # ~131 um lateral throats
        s = (L_box / 5) / lam * pitch_factor             # pitch 0.7 mm / lam (x pitch_factor)
        n = 5 * lam if n_side is None else int(n_side)
        L_box = n * s
        rng = np.random.default_rng(seed)
        if pore_radius is not None:
            if cv > 0.0:
                raise ValueError("pore_radius and cv > 0 cannot be combined")
            if float(pore_radius) <= 0.0:
                raise ValueError("pore_radius must be positive")
            r = np.full(n * n, float(pore_radius))
        elif cv > 0.0:
            sigma = np.sqrt(np.log(1.0 + cv ** 2))
            r = np.median(base_r) / lam * rng.lognormal(-0.5 * sigma ** 2, sigma, size=n * n)
            r *= np.sqrt(np.sum(base_r ** 2) / 25.0 * n * n / lam ** 2 / np.sum(r ** 2))   # same total opening area
            r = np.minimum(r, 0.475 * s)
        else:
            r = rng.choice(base_r, size=n * n, replace=True) / lam
        self.cv = cv
        if throat_shape not in ("square", "circle"):
            raise ValueError("throat_shape must be 'square' or 'circle'")
        self.throat_shape = throat_shape
        n_layers = int(n_layers)
        if n_layers < 1:
            raise ValueError("n_layers must be at least 1")
        dz = float(s if dz is None else dz)
        if dz <= 0.0:
            raise ValueError("dz must be positive")
        uniform = pore_radius is None and cv <= 0.0     # one template radius for every throat
        # deeper layers repeat the surface pattern's statistics; the surface draw above is
        # untouched, so a one-layer network is exactly what it was before n_layers existed
        if n_layers > 1:
            if pore_radius is not None:
                r_deep = np.full(n * n * (n_layers - 1), float(pore_radius))
            elif cv > 0.0:
                sigma = np.sqrt(np.log(1.0 + cv ** 2))
                r_deep = np.median(base_r) / lam * rng.lognormal(-0.5 * sigma ** 2, sigma,
                                                                 size=n * n * (n_layers - 1))
                r_deep *= np.sqrt(np.sum(base_r ** 2) / 25.0 * n * n / lam ** 2
                                  / np.sum(r_deep.reshape(n_layers - 1, -1) ** 2, axis=1).mean())
                r_deep = np.minimum(r_deep, 0.475 * s)
            else:
                r_deep = rng.choice(base_r, size=n * n * (n_layers - 1), replace=True) / lam
            r = np.concatenate([r, r_deep])
        ix, iy = np.meshgrid(np.arange(n), np.arange(n), indexing="ij")
        x1 = (ix.ravel() + 0.5) * s
        y1 = (iy.ravel() + 0.5) * s
        self.n_layers = n_layers
        self.dz = dz
        self.n_side = n                                  # pores per side of the periodic box
        self.n_surface = n * n
        self.x = np.tile(x1, n_layers)
        self.y = np.tile(y1, n_layers)
        self.z = -np.repeat(np.arange(n_layers) * dz, n * n)     # 0 at the surface, down is negative
        self.r = r
        self.V = (2 * r) ** 3                            # cubic pores as in QSPNM
        self.A_open = (2 * r) ** 2
        self.N = n * n * n_layers
        self.surface = np.arange(n * n)                  # the pores that open into the gas domain
        self.bottom = np.arange(n * n) + n * n * (n_layers - 1)
        self.L_box = L_box
        self.lam = lam

        def throat_radius(a, b):
            return base_rt / lam if uniform else (base_rt / np.median(base_r)) * min(r[a], r[b])

        # lateral throats: 4-neighbour periodic lattice, in every layer
        idx0 = np.arange(n * n).reshape(n, n)
        t = []
        for k in range(n_layers if n >= 2 else 0):        # n = 1: a periodic roll would be a
            off = k * n * n                               # self-loop (it cancels, but says nothing)
            for dx, dy in ((1, 0), (0, 1)):
                j_idx = np.roll(idx0, -dx, axis=0)
                j_idx = np.roll(j_idx, -dy, axis=1)
                for a0, b0 in zip(idx0.ravel(), j_idx.ravel()):
                    a, b = a0 + off, b0 + off
                    t.append((a, b, throat_radius(a, b), max(s - r[a] - r[b], 1e-5)))
        # vertical throats: same template, length dz minus the two pore half-sizes
        for k in range(n_layers - 1):
            for a0 in range(n * n):
                a, b = a0 + k * n * n, a0 + (k + 1) * n * n
                t.append((a, b, throat_radius(a, b), max(dz - r[a] - r[b], 1e-5)))
        self.lat = t
        # reservoir throat: template = the dgf boundary throat (131 um, 0.7 mm); scales with r_i if
        # cv > 0.  Only the BOTTOM layer touches the reservoir; the conductance of every other pore
        # is zero, which is how solve_pressure / salt_step see "no reservoir here".
        rt_res = np.zeros(self.N)
        rt_res[self.bottom] = (base_rt / lam if uniform
                               else (base_rt / np.median(base_r)) * r[self.bottom])
        self.rt_res = rt_res
        self.Lt_res = np.full(self.N, Lt_res)

    @classmethod
    def from_dgf(cls, path, L_box=None, reservoir_label=1, throat_shape="square"):
        """Build the substrate from a DuMuX dgf instead of the generated lattice.

        Intended for ``dumux/pnm2pncmin2/appl/grids/pnm_1d_homogeneous.dgf``, the network the
        model was validated against: a single column of pores,
        the bottom one flagged ``reservoir_label`` and therefore NOT a degree of freedom but the
        Dirichlet boundary, with the throat that reaches it becoming the reservoir throat of the
        pore above.  Pore radii, pore volumes, throat radii and throat LENGTHS are the file's, not
        this class's templates -- that is the point of using it.

        Only a column is accepted (one pore per height).  ``L_box`` is the lateral size of the
        periodic voxel box the surface pore opens into; a 1-D file cannot state it, so the default
        is the column's own pore spacing, i.e. the packing is taken to be isotropic.
        """
        V_all, T_all = parse_dgf_2d(path)
        res = np.flatnonzero(V_all[:, 4] == reservoir_label)
        if res.size != 1:
            raise ValueError(f"{path}: expected one pore labelled {reservoir_label} (the "
                             f"reservoir), found {res.size}")
        res = int(res[0])
        keep = np.ones(len(V_all), bool); keep[res] = False
        idx = -np.ones(len(V_all), int); idx[keep] = np.arange(int(keep.sum()))
        height = V_all[keep, 1]                       # the dgf's second coordinate is the vertical
        self = cls.__new__(cls)
        self.reservoir_pore = res                     # its index in the FILE, for reading DuMuX output
        self.height = height                          # heights as the file states them [m]
        self.throat_shape = throat_shape
        self.cv = 0.0
        self.lam = None                               # not a scaled lattice; nothing reads it
        self.r = V_all[keep, 2]
        self.V = V_all[keep, 3]                       # the file's pore volumes, not (2 r)^3
        self.A_open = (2 * self.r) ** 2               # the QSPNM cube face, as coupling.interface wants
        self.N = int(keep.sum())
        self.z = height - height.max()                # 0 at the surface, down is negative
        levels = np.unique(self.z)
        if levels.size != self.N:
            raise ValueError(f"{path}: from_dgf builds one column, but {self.N} pores sit on "
                             f"{levels.size} distinct heights")
        self.n_layers = self.N
        self.n_side = 1
        self.n_surface = 1
        self.dz = float(np.median(np.diff(levels))) if self.N > 1 else 0.0
        self.L_box = float(self.dz if L_box is None else L_box)
        self.x = np.full(self.N, 0.5 * self.L_box)
        self.y = np.full(self.N, 0.5 * self.L_box)
        self.surface = np.flatnonzero(self.z == 0.0)
        self.bottom = np.flatnonzero(self.z == levels[0])
        self.lat, self.rt_res, self.Lt_res = [], np.zeros(self.N), np.ones(self.N)
        for a, b, rt, Lt, _ in T_all:
            a, b = int(a), int(b)
            if a == res or b == res:
                other = idx[b if a == res else a]
                self.rt_res[other] = rt
                self.Lt_res[other] = Lt
            else:
                self.lat.append((idx[a], idx[b], rt, Lt))
        return self

    @classmethod
    def from_dgf_dirichlet_boundary(cls, path, boundary_axis=1, boundary_side="min",
                                    throat_shape="square"):
        """Read a 2-D DGF whose complete lower/upper row is a Dirichlet reservoir.

        The selected boundary vertices are eliminated from the dynamic pore unknowns.  Edges
        wholly inside that fixed-pressure row are discarded, while an edge crossing from the
        row into the dynamic network becomes that dynamic pore's reservoir throat.  All other
        edges retain their file order and exact DGF geometry.

        This is the multi-pore analogue of :meth:`from_dgf`.  It is intentionally strict: each
        dynamic pore may have at most one edge to the reservoir, because the current pressure and
        transport kernels store one reservoir conductance per pore.  The original file indices
        and edge indices are retained for pore-by-pore and throat-by-throat comparison with
        DuMuX output.
        """
        if boundary_axis not in (0, 1):
            raise ValueError("boundary_axis must be 0 or 1 for a 2-D DGF")
        if boundary_side not in ("min", "max"):
            raise ValueError("boundary_side must be 'min' or 'max'")
        if throat_shape not in ("square", "circle"):
            raise ValueError("throat_shape must be 'square' or 'circle'")

        V_all, T_all = parse_dgf_2d(path)
        coordinate = V_all[:, boundary_axis]
        boundary_value = coordinate.min() if boundary_side == "min" else coordinate.max()
        scale = max(float(np.ptp(coordinate)), 1.0)
        reservoir = np.flatnonzero(np.isclose(coordinate, boundary_value, rtol=0.0,
                                              atol=1e-12 * scale))
        if reservoir.size == 0:
            raise ValueError(f"{path}: no vertices found on the {boundary_side} boundary")

        is_reservoir = np.zeros(len(V_all), dtype=bool)
        is_reservoir[reservoir] = True
        keep = ~is_reservoir
        file_index = np.flatnonzero(keep)
        file_to_dynamic = np.full(len(V_all), -1, dtype=int)
        file_to_dynamic[file_index] = np.arange(file_index.size)

        self = cls.__new__(cls)
        self.throat_shape = throat_shape
        self.cv = 0.0
        self.lam = None
        self.dgf_vertices = np.array(V_all, copy=True)
        self.dgf_edges = np.array(T_all, copy=True)
        self.file_vertex_index = file_index
        self.file_to_dynamic = file_to_dynamic
        self.reservoir_file_indices = reservoir
        self.reservoir_pores = reservoir
        self.height = V_all[keep, 1]
        self.x = V_all[keep, 0]
        self.y = V_all[keep, 1]
        self.z = self.height - self.height.max()
        self.r = V_all[keep, 2]
        self.V = V_all[keep, 3]
        self.A_open = (2.0 * self.r) ** 2
        self.N = int(file_index.size)

        levels = np.unique(self.height)
        lateral = np.unique(self.x)
        self.surface = np.flatnonzero(self.height == levels[-1])
        self.bottom = np.flatnonzero(self.height == levels[0])
        self.n_layers = int(levels.size)
        self.n_side = int(lateral.size)
        self.n_surface = int(self.surface.size)
        self.dz = float(np.median(np.diff(levels))) if levels.size > 1 else 0.0
        pitch = float(np.median(np.diff(lateral))) if lateral.size > 1 else self.dz
        self.L_box = float(np.ptp(lateral) + pitch) if lateral.size > 1 else float(pitch)

        self.lat = []
        self.lat_file_edge_index = []
        self.rt_res = np.zeros(self.N)
        self.Lt_res = np.ones(self.N)
        self.res_file_edge_index = np.full(self.N, -1, dtype=int)
        reservoir_internal_edges = []
        for edge_idx, (a_value, b_value, rt, Lt, _) in enumerate(T_all):
            a, b = int(a_value), int(b_value)
            a_res, b_res = bool(is_reservoir[a]), bool(is_reservoir[b])
            if a_res and b_res:
                reservoir_internal_edges.append(edge_idx)
                continue
            if a_res != b_res:
                dynamic_file_idx = b if a_res else a
                dynamic_idx = int(file_to_dynamic[dynamic_file_idx])
                if self.res_file_edge_index[dynamic_idx] >= 0:
                    raise ValueError(
                        f"{path}: dynamic pore {dynamic_file_idx} has more than one reservoir edge"
                    )
                self.rt_res[dynamic_idx] = rt
                self.Lt_res[dynamic_idx] = Lt
                self.res_file_edge_index[dynamic_idx] = edge_idx
                continue
            self.lat.append((int(file_to_dynamic[a]), int(file_to_dynamic[b]), rt, Lt))
            self.lat_file_edge_index.append(edge_idx)

        self.lat_file_edge_index = np.asarray(self.lat_file_edge_index, dtype=int)
        self.reservoir_internal_file_edge_indices = np.asarray(reservoir_internal_edges, dtype=int)
        return self

    def throat_area(self, rt):
        """Cross-sectional area of a throat of inscribed radius rt [m^2]."""
        rt = np.asarray(rt, float)
        return 4.0 * rt ** 2 if self.throat_shape == "square" else np.pi * rt ** 2

    def conductances(self):
        """Single-phase throat conductance [m^3 / (Pa s)]: Patzek & Silin for a square,
        Hagen-Poiseuille for a circle."""
        def g(rt, L):
            if self.throat_shape == "square":
                return K_SQUARE * self.throat_area(rt) ** 2 * G_SQUARE / (MU * L)
            return np.pi * rt ** 4 / (8 * MU * L)
        g_lat = np.array([g(rt, L) for _, _, rt, L in self.lat])
        g_res = g(self.rt_res, self.Lt_res)
        return g_lat, g_res

    def diffusive(self):
        """Throat diffusive conductance D A / L [m^3/s] (multiplied by the molar density of the
        brine where it is used, as in DuMuX's Fick law)."""
        d_lat = np.array([D_NACL * self.throat_area(rt) / L for _, _, rt, L in self.lat])
        d_res = D_NACL * self.throat_area(self.rt_res) / self.Lt_res
        return d_lat, d_res
