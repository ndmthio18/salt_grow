"""Reading and writing a DuMuX / QSPNM ``.dgf`` pore-network file.  No physics: numbers in,
arrays out (and back)."""
import numpy as np


def parse_dgf_2d(path):
    lines = open(path).read().splitlines()
    iv = next(i for i, l in enumerate(lines) if l.startswith("SIMPLEX"))
    V = np.array([[float(x) for x in l.split()] for l in lines[5:iv] if l and l[0] not in "#%"])
    E = np.array([[float(x) for x in l.split()] for l in lines[iv + 2:] if l and l[0] not in "#%"
                  and len(l.split()) == 5])
    return V, E


def write_column_dgf(path, n_pores, spacing, pore_radius, pore_volume, throat_radius, throat_length):
    """Write a vertical column of ``n_pores`` identical free pores over one reservoir pore, in the
    layout ``parse_dgf_2d`` reads (five header lines, the second coordinate vertical).

    All lengths in m, the volume in m^3.  The reservoir pore (label 1) sits one ``spacing`` below
    the lowest free pore and is joined to it by a throat like every other; the top pore carries
    label 2.  The radius and the volume of a pore are written independently, as the format has
    them.  Numbers are written with ``repr``, so reading the file gives back the same floats.
    """
    if n_pores < 1:
        raise ValueError("a column needs at least one free pore")
    for name, value in (("spacing", spacing), ("pore_radius", pore_radius),
                        ("pore_volume", pore_volume), ("throat_radius", throat_radius),
                        ("throat_length", throat_length)):
        if not value > 0.0:
            raise ValueError(f"{name} must be positive")
    lines = ["DGF",
             "Vertex % Coordinates, volumes and boundary flags of the pore bodies",
             "parameters 3",
             "% Vertex parameters: PoreInscribedRadius PoreVolume PoreLabel",
             "% Element parameters: ThroatInscribedRadius ThroatLength ThroatLabel"]
    for k in range(n_pores + 1):                       # vertex 0 is the reservoir pore
        label = 1 if k == 0 else (2 if k == n_pores else -1)
        lines.append(f"0 {float(k * spacing)!r} {float(pore_radius)!r} {float(pore_volume)!r} {label}")
    lines += ["#", "SIMPLEX % Connections of the pore bodies (pore throats)", "parameters 3"]
    for k in range(n_pores):
        label = 1 if k == 0 else (2 if k == n_pores - 1 else -1)
        lines.append(f"{k} {k + 1} {float(throat_radius)!r} {float(throat_length)!r} {label}")
    lines.append("#")
    with open(path, "w") as f:
        f.write("\n".join(lines) + "\n")
    return path
