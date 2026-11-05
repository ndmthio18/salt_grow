from __future__ import annotations

from dataclasses import dataclass
from math import sqrt
from typing import Literal

import numpy as np
from scipy.integrate import solve_ivp

from .analytic import dimensionless_groups
from .grid import cells_across_top_layer, top_refined_grid
from .parameters import NumericalParameters, PhysicalParameters

Closure = Literal["B1", "B2"]


@dataclass
class DynamicResult:
    closure: Closure
    physical: PhysicalParameters
    numerical: NumericalParameters
    faces: np.ndarray
    centers: np.ndarray
    widths: np.ndarray
    time: np.ndarray
    epsilon: np.ndarray
    concentration: np.ndarray
    height: np.ndarray
    top_solid: np.ndarray
    internal_solid: np.ndarray
    bottom_salt: np.ndarray
    bottom_solution: np.ndarray
    evaporated: np.ndarray
    status: int
    message: str


class DynamicALEModel:
    """Conservative ALE finite-volume model per unit cross-sectional area."""

    N_INVENTORIES = 5

    def __init__(self, physical: PhysicalParameters, numerical: NumericalParameters,
                 closure: Closure):
        if closure not in ("B1", "B2"):
            raise ValueError(f"unknown closure {closure}")
        self.p = physical
        self.n = numerical
        self.closure = closure
        self.faces, self.centers, self.widths = top_refined_grid(
            numerical.n_cells, numerical.grid_power
        )

    def coefficients(self, epsilon: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        if self.closure == "B1":
            return (
                np.full_like(epsilon, self.p.a_v_ref),
                np.full_like(epsilon, self.p.D_eff_ref),
            )
        a_v = 36.0 * (1.0 - epsilon) / self.p.d_sub
        d_eff = np.sqrt(epsilon) * self.p.D_molecular
        return a_v, d_eff

    def unpack(self, y: np.ndarray) -> tuple[np.ndarray, np.ndarray, float, np.ndarray]:
        n = self.n.n_cells
        epsilon = y[:n]
        concentration = y[n:2 * n]
        height = float(y[2 * n])
        inventories = y[2 * n + 1:]
        return epsilon, concentration, height, inventories

    def top_state(self, concentration_last: float, epsilon_last: float,
                  d_eff_last: float, height: float) -> tuple[float, float]:
        p = self.p
        distance = height * (self.faces[-1] - self.centers[-1])
        conductance = p.rho_solution * epsilon_last * d_eff_last / distance
        kinetic_mass = p.beta * p.rho_solution * (1.0 - p.epsilon_birth) * p.k_reaction
        # j*C - beta*rho*(1-e_birth)*kr*(1-C)*(C-Ceq) - K*(C-C_N) = 0
        qa = kinetic_mass
        qb = p.evaporation_mass_flux - kinetic_mass * (1.0 + p.C_eq) - conductance
        qc = kinetic_mass * p.C_eq + conductance * concentration_last
        discriminant = qb * qb - 4.0 * qa * qc
        if discriminant < 0.0:
            raise RuntimeError(f"negative top-boundary discriminant: {discriminant}")
        root = sqrt(discriminant)
        candidates = np.array([(-qb - root) / (2.0 * qa), (-qb + root) / (2.0 * qa)])
        physical = candidates[(candidates >= 0.0) & (candidates <= 1.0)]
        if physical.size:
            c_top = float(physical[np.argmin(np.abs(physical - p.C_eq))])
        else:
            c_top = float(candidates[np.argmin(np.abs(candidates - p.C_eq))])
        velocity = p.beta * p.rho_solution / p.rho_crystal * p.k_reaction * (c_top - p.C_eq)
        return c_top, velocity

    def flux_state(self, epsilon: np.ndarray, concentration: np.ndarray, height: float):
        p = self.p
        n = self.n.n_cells
        volumes = height * self.widths
        a_v, d_eff = self.coefficients(epsilon)
        reaction = p.rho_solution * a_v * p.k_reaction * (concentration - p.C_eq)
        c_top, velocity = self.top_state(concentration[-1], epsilon[-1], d_eff[-1], height)

        liquid_flux = np.empty(n + 1)
        structural_density_birth = (
            p.rho_solution * p.epsilon_birth
            + p.rho_crystal * (1.0 - p.epsilon_birth)
        )
        liquid_flux[-1] = structural_density_birth * velocity + p.evaporation_mass_flux
        density_factor = 1.0 - p.rho_solution / p.rho_crystal
        for i in range(n - 1, -1, -1):
            liquid_flux[i] = liquid_flux[i + 1] + density_factor * reaction[i] * volumes[i]

        mesh_velocity = self.faces * velocity
        epsilon_face = np.empty(n + 1)
        epsilon_face[0] = epsilon[0]
        epsilon_face[-1] = p.epsilon_birth
        epsilon_face[1:-1] = 0.5 * (epsilon[:-1] + epsilon[1:])

        salt_flux = np.empty(n + 1)
        # Bottom reservoir: C=Ceq; include advective and diffusive salt flux.
        bottom_distance = height * (self.centers[0] - self.faces[0])
        bottom_diffusive = (
            -p.rho_solution * epsilon[0] * d_eff[0]
            * (concentration[0] - p.C_eq) / bottom_distance
        )
        relative_liquid_bottom = liquid_flux[0]
        salt_flux[0] = p.C_eq * relative_liquid_bottom + bottom_diffusive

        for face in range(1, n):
            left = face - 1
            right = face
            distance = height * (self.centers[right] - self.centers[left])
            k_left = epsilon[left] * d_eff[left]
            k_right = epsilon[right] * d_eff[right]
            half_left = height * (self.faces[face] - self.centers[left])
            half_right = height * (self.centers[right] - self.faces[face])
            k_face = distance / (half_left / k_left + half_right / k_right)
            diffusive = (
                -p.rho_solution * k_face
                * (concentration[right] - concentration[left]) / distance
            )
            relative_liquid = (
                liquid_flux[face]
                - p.rho_solution * epsilon_face[face] * mesh_velocity[face]
            )
            advected_c = concentration[left] if relative_liquid >= 0.0 else concentration[right]
            salt_flux[face] = relative_liquid * advected_c + diffusive

        # Moving-interface salt balance: dissolved salt becomes new-layer solid.
        salt_flux[-1] = p.rho_crystal * (1.0 - p.epsilon_birth) * velocity

        epsilon_mesh_flux = np.empty(n + 1)
        epsilon_mesh_flux[0] = 0.0
        # Relative material velocity is -mesh_velocity, so the upwind state is right.
        epsilon_mesh_flux[1:-1] = -mesh_velocity[1:-1] * epsilon[1:]
        epsilon_mesh_flux[-1] = -velocity * p.epsilon_birth
        return {
            "volumes": volumes,
            "a_v": a_v,
            "d_eff": d_eff,
            "reaction": reaction,
            "c_top": c_top,
            "velocity": velocity,
            "liquid_flux": liquid_flux,
            "salt_flux": salt_flux,
            "epsilon_mesh_flux": epsilon_mesh_flux,
        }

    def rhs(self, _time: float, y: np.ndarray) -> np.ndarray:
        p = self.p
        n = self.n.n_cells
        epsilon, concentration, height, _ = self.unpack(y)
        state = self.flux_state(epsilon, concentration, height)
        volumes = state["volumes"]
        reaction = state["reaction"]
        velocity = state["velocity"]

        d_e_integral = (
            state["epsilon_mesh_flux"][:-1]
            - state["epsilon_mesh_flux"][1:]
            - reaction * volumes / p.rho_crystal
        )
        d_q_integral = (
            state["salt_flux"][:-1]
            - state["salt_flux"][1:]
            - reaction * volumes
        )
        d_epsilon = (d_e_integral - epsilon * self.widths * velocity) / volumes
        pore_integral = epsilon * volumes
        d_concentration = (
            d_q_integral - p.rho_solution * concentration * d_e_integral
        ) / (p.rho_solution * pore_integral)

        derivative = np.empty(2 * n + 1 + self.N_INVENTORIES)
        derivative[:n] = d_epsilon
        derivative[n:2 * n] = d_concentration
        derivative[2 * n] = velocity
        inv = 2 * n + 1
        derivative[inv + 0] = p.rho_crystal * (1.0 - p.epsilon_birth) * velocity
        derivative[inv + 1] = float(np.sum(reaction * volumes))
        derivative[inv + 2] = state["salt_flux"][0]
        derivative[inv + 3] = state["liquid_flux"][0]
        derivative[inv + 4] = p.evaporation_mass_flux
        return derivative

    def initial_state(self) -> np.ndarray:
        n = self.n.n_cells
        y0 = np.zeros(2 * n + 1 + self.N_INVENTORIES)
        y0[:n] = self.p.epsilon_ref
        y0[n:2 * n] = self.p.C_eq
        y0[2 * n] = self.p.h0
        return y0

    def _events(self):
        n = self.n.n_cells
        num = self.n

        def eps_low(_t, y):
            return float(np.min(y[:n]) - num.porosity_stop_min)

        def eps_high(_t, y):
            return float(num.porosity_stop_max - np.max(y[:n]))

        def c_low(_t, y):
            return float(np.min(y[n:2 * n]) - num.concentration_stop_min)

        def c_high(_t, y):
            return float(num.concentration_stop_max - np.max(y[n:2 * n]))

        for event in (eps_low, eps_high, c_low, c_high):
            event.terminal = True
            event.direction = -1
        return [eps_low, eps_high, c_low, c_high]

    def run(self) -> DynamicResult:
        n = self.n.n_cells
        output_times = np.arange(0.0, self.n.t_end + 0.5 * self.n.output_interval,
                                 self.n.output_interval)
        if output_times[-1] > self.n.t_end:
            output_times[-1] = self.n.t_end
        atol = np.concatenate([
            np.full(2 * n, self.n.atol_field),
            np.array([self.n.atol_height]),
            np.full(self.N_INVENTORIES, self.n.atol_inventory),
        ])
        event_functions = self._events()
        solution = solve_ivp(
            self.rhs,
            (0.0, self.n.t_end),
            self.initial_state(),
            method="BDF",
            t_eval=output_times,
            rtol=self.n.rtol,
            atol=atol,
            max_step=self.n.max_step,
            events=event_functions,
        )
        time = solution.t
        state_matrix = solution.y
        message = solution.message
        event_names = ("porosity_min", "porosity_max", "concentration_min", "concentration_max")
        for event_index, event_times in enumerate(solution.t_events):
            if event_times.size:
                event_time = float(event_times[-1])
                event_state = solution.y_events[event_index][-1]
                message = f"{solution.message} event={event_names[event_index]} at t={event_time:.9g} s"
                if time.size == 0 or abs(time[-1] - event_time) > 1e-10 * max(1.0, event_time):
                    time = np.append(time, event_time)
                    state_matrix = np.column_stack([state_matrix, event_state])
                break
        epsilon = state_matrix[:n, :]
        concentration = state_matrix[n:2 * n, :]
        height = state_matrix[2 * n, :]
        inventory = state_matrix[2 * n + 1:, :]
        return DynamicResult(
            closure=self.closure,
            physical=self.p,
            numerical=self.n,
            faces=self.faces,
            centers=self.centers,
            widths=self.widths,
            time=time,
            epsilon=epsilon,
            concentration=concentration,
            height=height,
            top_solid=inventory[0],
            internal_solid=inventory[1],
            bottom_salt=inventory[2],
            bottom_solution=inventory[3],
            evaporated=inventory[4],
            status=solution.status,
            message=message,
        )


def snapshot_diagnostics(result: DynamicResult, index: int) -> dict[str, float | np.ndarray]:
    model = DynamicALEModel(result.physical, result.numerical, result.closure)
    eps = result.epsilon[:, index]
    concentration = result.concentration[:, index]
    height = float(result.height[index])
    state = model.flux_state(eps, concentration, height)
    p = result.physical
    volumes = height * result.widths
    dissolved = float(np.sum(p.rho_solution * eps * concentration * volumes))
    solid = float(np.sum(p.rho_crystal * (1.0 - eps) * volumes))
    bulk = float(np.sum((p.rho_solution * eps + p.rho_crystal * (1.0 - eps)) * volumes))

    eps0 = p.epsilon_ref
    initial_dissolved = p.rho_solution * eps0 * p.C_eq * p.h0
    initial_solid = p.rho_crystal * (1.0 - eps0) * p.h0
    initial_bulk = (p.rho_solution * eps0 + p.rho_crystal * (1.0 - eps0)) * p.h0
    salt_residual = dissolved + solid - (initial_dissolved + initial_solid + result.bottom_salt[index])
    bulk_residual = bulk - (initial_bulk + result.bottom_solution[index] - result.evaporated[index])

    reaction_positive = np.maximum(state["reaction"], 0.0)
    weights = reaction_positive * volumes
    total_positive = float(np.sum(weights))
    if total_positive > 0.0:
        cumulative = np.cumsum(weights[::-1])
        position = int(np.searchsorted(cumulative, 0.95 * total_positive, side="left"))
        cell = len(weights) - 1 - min(position, len(weights) - 1)
        l95 = height * (1.0 - result.faces[cell])
    else:
        l95 = 0.0

    a_top = float(state["a_v"][-1])
    d_top = float(state["d_eff"][-1])
    groups = dimensionless_groups(p, height, float(eps[-1]), a_top, d_top)
    return {
        "velocity": float(state["velocity"]),
        "c_top": float(state["c_top"]),
        "reaction": state["reaction"],
        "dissolved_salt": dissolved,
        "solid_salt": solid,
        "bulk_mass": bulk,
        "salt_residual": salt_residual,
        "bulk_residual": bulk_residual,
        "salt_residual_normalized": salt_residual / max(abs(dissolved + solid), 1e-30),
        "bulk_residual_normalized": bulk_residual / max(abs(bulk), 1e-30),
        "L95": l95,
        "Pe": groups.Pe,
        "Da": groups.Da,
        "Da_top": groups.Da_top,
        "Da_star": groups.Da_star,
        "delta_C": groups.delta_C,
        "cells_across_delta": cells_across_top_layer(result.faces, height, groups.delta_C),
        "epsilon_min": float(np.min(eps)),
        "epsilon_top_cell": float(eps[-1]),
        "concentration_min": float(np.min(concentration)),
        "concentration_max": float(np.max(concentration)),
    }
