"""Independent cell/face reconstruction, separate from flux_state and inventories."""
import numpy as np


def local_audit(model, y, dy=None):
    p, n = model.p, model.n.n_cells
    eps, c, h, _ = model.unpack(y)
    if dy is None:
        dy = model.rhs(0, y)
    v = dy[2*n]
    vol = h*model.widths
    dvol = v*model.widths
    dE = dy[:n]*vol+eps*dvol
    dQ = p.rho_solution*(eps*vol*dy[n:2*n]+c*dE)
    av = np.full(n, p.a_v_ref) if model.closure == "B1" else 36*(1-eps)/p.d_sub
    D = np.full(n, p.D_eff_ref) if model.closure == "B1" else np.sqrt(eps)*p.D_molecular
    R = p.rho_solution*av*p.k_reaction*(c-p.C_eq)
    ct = p.C_eq+v*p.rho_crystal/(p.beta*p.rho_solution*p.k_reaction)
    jf = p.rho_crystal*(1-p.epsilon_birth)*v
    F = np.empty(n+1)
    F[-1] = p.evaporation_mass_flux+(p.rho_solution*p.epsilon_birth+p.rho_crystal*(1-p.epsilon_birth))*v
    for i in range(n-1, -1, -1):
        F[i] = F[i+1]+(1-p.rho_solution/p.rho_crystal)*R[i]*vol[i]
    em = np.zeros(n+1)
    S = np.empty(n+1)
    L = np.empty(n+1)
    for face in range(n+1):
        if face == 0:
            ca, cb = p.C_eq, c[0]
            resistance = h*model.centers[0]/(eps[0]*D[0])
            cf = p.C_eq
        elif face == n:
            em[face] = -v*p.epsilon_birth
            ca, cb, cf = c[-1], ct, ct
            resistance = h*(1-model.centers[-1])/(eps[-1]*D[-1])
        else:
            w = v*model.faces[face]
            em[face] = -w*(eps[face] if w >= 0 else eps[face-1])
            ca, cb = c[face-1], c[face]
            resistance = h*((model.faces[face]-model.centers[face-1])/(eps[face-1]*D[face-1])
                            +(model.centers[face]-model.faces[face])/(eps[face]*D[face]))
            cf = ca if F[face]+p.rho_solution*em[face] >= 0 else cb
        L[face] = F[face]+p.rho_solution*em[face]
        delta_c = (v*p.rho_crystal/(p.beta*p.rho_solution*p.k_reaction)-(c[-1]-p.C_eq)) if face == n else cb-ca
        S[face] = cf*L[face]-p.rho_solution*delta_c/resistance
    # Kinetic/top diffusion evaluated with mass fractions has a subtraction floor.
    # Keep top residual separate; cell balances use physical Jf boundary.
    interface = S[-1]-jf
    S[-1] = jf
    liquid = p.rho_solution*dE+L[1:]-L[:-1]+R*vol
    salt = dQ+S[1:]-S[:-1]+R*vol
    solid_mesh_flux = -model.faces*v-em
    solid = p.rho_crystal*(dvol-dE+np.diff(solid_mesh_flux))-R*vol
    scale = max(abs(p.evaporation_mass_flux), abs(jf), 1e-30)
    result = {}
    for name, residual in (("liquid", liquid), ("salt", salt), ("solid", solid)):
        result[name+"_linf"] = float(max(abs(residual))/scale)
        result[name+"_weighted_l1"] = float(np.sum(abs(residual)*vol)/np.sum(vol)/scale)
    result["interface_salt_relative"] = float(abs(interface)/scale)
    result["interface_liquid_relative"] = float(abs(L[-1]-p.evaporation_mass_flux-jf)/scale)
    result["front_volume_relative"] = float(abs(jf-p.rho_crystal*(1-p.epsilon_birth)*v)/scale)
    return result
