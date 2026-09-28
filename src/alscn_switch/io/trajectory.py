"""Memory-bounded extended-XYZ reader with a native C fast path."""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterator

import numpy as np

try:
    from .. import _xyzstream
except ImportError:  # source checkout before editable install
    _xyzstream = None


@dataclass(slots=True)
class Frame:
    index: int
    time_fs: float
    species: np.ndarray  # uint8: 1=Al, 2=Sc, 3=N
    positions: np.ndarray
    unwrapped_positions: np.ndarray
    lattice: np.ndarray
    charges: np.ndarray
    forces: np.ndarray
    bec: np.ndarray

    @property
    def natoms(self) -> int:
        return int(self.species.size)


def _frame_from_native(data: dict) -> Frame:
    n = data["natoms"]
    return Frame(
        index=data["frame"], time_fs=data["time_fs"],
        species=np.frombuffer(data["species"], dtype=np.uint8),
        positions=np.frombuffer(data["positions"], dtype=np.float64).reshape(n, 3),
        unwrapped_positions=np.frombuffer(data["unwrapped_positions"], dtype=np.float64).reshape(n, 3),
        lattice=np.asarray(data["lattice"], dtype=np.float64).reshape(3, 3),
        charges=np.frombuffer(data["charges"], dtype=np.float64),
        forces=np.frombuffer(data["forces"], dtype=np.float64).reshape(n, 3),
        bec=np.frombuffer(data["bec"], dtype=np.float64).reshape(n, 9),
    )


def _fallback_frames(path: Path, stride: int) -> Iterator[Frame]:
    """Strict pure-Python extxyz fallback, retaining only one frame at a time."""
    with path.open("rt", encoding="utf-8") as stream:
        index = 0
        while True:
            count_line = stream.readline()
            if not count_line:
                return
            if not count_line.strip():
                continue
            try:
                n = int(count_line)
            except ValueError as exc:
                raise ValueError(f"bad atom count at frame {index + 1}") from exc
            header = stream.readline()
            if not header:
                raise ValueError(f"truncated header at frame {index + 1}")
            props_start = header.find("Properties=")
            if props_start < 0:
                raise ValueError(f"missing Properties at frame {index + 1}")
            tail = header[props_start + len("Properties="):]
            quote = tail[0] if tail and tail[0] in "\"'" else ""
            prop_text = tail[1:].split(quote, 1)[0] if quote else tail.split()[0]
            fields = prop_text.split(":")
            cols: dict[str, tuple[int, int]] = {}
            col = 0
            for j in range(0, len(fields) - 2, 3):
                name, _kind, width = fields[j:j + 3]
                w = int(width)
                cols[name] = (col, w)
                col += w
            sp_key = "species" if "species" in cols else "element"
            if sp_key not in cols or "pos" not in cols or cols[sp_key][1] != 1 or cols["pos"][1] != 3:
                raise ValueError(f"unsupported Properties schema at frame {index + 1}")
            needed = max(c + w for c, w in cols.values())
            wanted = index % stride == 0
            species = np.empty(n, np.uint8) if wanted else None
            pos = np.empty((n, 3)) if wanted else None
            uw = np.empty((n, 3)) if wanted else None
            charge = np.full(n, np.nan) if wanted else None
            force = np.full((n, 3), np.nan) if wanted else None
            bec = np.full((n, 9), np.nan) if wanted else None
            codes = {"Al": 1, "Sc": 2, "N": 3}
            for atom in range(n):
                row = stream.readline()
                if not row:
                    raise ValueError(f"truncated frame {index + 1}, atom {atom}")
                if not wanted:
                    continue
                tokens = row.split()
                if len(tokens) < needed:
                    raise ValueError(f"short atom row at frame {index + 1}, atom {atom}")
                elem = tokens[cols[sp_key][0]]
                species[atom] = codes.get(elem, 0)
                c, _ = cols["pos"]
                pos[atom] = [float(x) for x in tokens[c:c + 3]]
                if "unwrapped_position" in cols and cols["unwrapped_position"][1] == 3:
                    c = cols["unwrapped_position"][0]
                    uw[atom] = [float(x) for x in tokens[c:c + 3]]
                else:
                    uw[atom] = pos[atom]
                if "charge" in cols and cols["charge"][1] == 1:
                    charge[atom] = float(tokens[cols["charge"][0]])
                if "forces" in cols and cols["forces"][1] == 3:
                    c = cols["forces"][0]
                    force[atom] = [float(x) for x in tokens[c:c + 3]]
                if "bec" in cols and cols["bec"][1] == 9:
                    c = cols["bec"][0]
                    bec[atom] = [float(x) for x in tokens[c:c + 9]]
            if wanted:
                lattice_text = header.split("Lattice=", 1)[1].split('"', 2)[1]
                lattice = np.fromstring(lattice_text, sep=" ").reshape(3, 3)
                time = header.split("Time=", 1)[1].split()[0] if "Time=" in header else "nan"
                yield Frame(index + 1, float(time), species, pos, uw, lattice, charge, force, bec)
            index += 1


def iter_extxyz(path: str | Path, stride: int = 1, backend: str = "auto") -> Iterator[Frame]:
    """Yield selected extxyz frames; memory use is O(one frame).

    `backend` is `auto`, `c`, or `python`. The C extension is built by pip.
    """
    path = Path(path)
    if stride < 1:
        raise ValueError("stride must be >= 1")
    if backend not in {"auto", "c", "python", "ase"}:
        raise ValueError("backend must be auto, c, python, or ase")
    if backend == "c" and _xyzstream is None:
        raise RuntimeError("native C reader is unavailable; install the package with pip")
    if backend == "ase":
        from ase.io import iread
        from ase.data import chemical_symbols
        for index, atoms in enumerate(iread(str(path), index=":")):
            if index % stride: continue
            symbols = atoms.get_chemical_symbols()
            codes = {"Al": 1, "Sc": 2, "N": 3}
            n = len(atoms); arr = atoms.arrays
            time_fs = atoms.info.get("Time", atoms.info.get("time", float("nan")))
            try: charges = atoms.get_charges()
            except Exception: charges = np.full(n, np.nan)
            try: forces = atoms.get_forces(apply_constraint=False)
            except Exception: forces = np.full((n, 3), np.nan)
            yield Frame(index + 1, float(time_fs), np.asarray([codes.get(x, 0) for x in symbols], np.uint8),
                        atoms.positions.copy(), arr.get("unwrapped_position", atoms.positions).copy(),
                        np.asarray(atoms.cell), np.asarray(charges).copy(), np.asarray(forces).copy(),
                        arr.get("bec", np.full((n, 9), np.nan)).copy())
        return
    if backend == "python" or _xyzstream is None:
        yield from _fallback_frames(path, stride)
        return
    cap = _xyzstream.open_stream(str(path))
    try:
        while True:
            data = _xyzstream.next_frame(cap)
            if data is None:
                return
            frame_no = data["frame"] - 1
            if frame_no % stride == 0:
                yield _frame_from_native(data)
    finally:
        _xyzstream.close_stream(cap)


def benchmark_readers(path: str | Path, frames: int = 10) -> dict[str, float]:
    """Compare reader throughput on the same leading frames."""
    import time
    result = {}
    for backend in ("ase", "c"):
        start = time.perf_counter()
        atoms = 0
        for frame in iter_extxyz(path, backend=backend):
            atoms += frame.natoms
            if frame.index >= frames:
                break
        elapsed = time.perf_counter() - start
        result[backend] = {"seconds": elapsed, "frames": min(frames, frame.index if atoms else 0),
                           "atoms_per_second": atoms / elapsed if elapsed else float("inf")}
    return result
