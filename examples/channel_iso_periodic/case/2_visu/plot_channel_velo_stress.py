"""
README
------
Plot CHAPSim2 channel statistics against MKM reference data.

Supported plot groups:
- velocity: ``u1``, ``u2``, ``u3``
- pressure: ``p``
- stress: ``uu``, ``uv``, ``uw``, ``vv``, ``vw``, ``ww``
- vorticity: ``omega_x``, ``omega_y``, ``omega_z``

The vorticity plots use CHAPSim gradient files ``domain1_tsp_avg_dudx**_<iter>.dat``
when available. If these are not present, the script falls back to the 1D channel
relations from the mean profiles:
- ``omega_x = d(u_z)/dy``
- ``omega_y = 0``
- ``omega_z = -d(u_x)/dy``

Reference mean-vorticity curves are derived from ``chan<Re_tau>.means``:
- ``omega_x^+ = (dW^+/d(y/h)) / Re_tau``
- ``omega_y^+ = 0``
- ``omega_z^+ = -(dU^+/d(y/h)) / Re_tau``

Note: ``chan<Re_tau>.vortvar`` contains vorticity fluctuation variances, not mean
vorticity, so it is not used here.

Examples
--------
python3 plot_channel_velo_stress.py --dns-time 85000 --re 2800
python3 plot_channel_velo_stress.py --dns-time 85000 --re 2800 --groups velocity pressure
python3 plot_channel_velo_stress.py --dns-time 85000 --re 2800 --groups stress vorticity --input-dir ../1_data
"""

from __future__ import annotations

import argparse
import math
from pathlib import Path

import matplotlib.cm as cm
import matplotlib.pyplot as plt
import numpy as np
from pylab import rcParams


FILEMAP_MEAN = {"u1": "u1", "u2": "u2", "u3": "u3", "p": "pr"}
FILEMAP_REY = {
    "uu": "uu11",
    "uv": "uu12",
    "uw": "uu13",
    "vv": "uu22",
    "vw": "uu23",
    "ww": "uu33",
}
DUDX_COMPONENTS = (
    "dudx11",
    "dudx12",
    "dudx13",
    "dudx21",
    "dudx22",
    "dudx23",
    "dudx31",
    "dudx32",
    "dudx33",
)

PLOT_GROUPS = {
    "velocity": ("u1", "u2", "u3"),
    "pressure": ("p",),
    "stress": ("uu", "uv", "uw", "vv", "vw", "ww"),
    "vorticity": ("omega_x", "omega_y", "omega_z"),
}

PARAMS = {
    "u1": {
        "group": "velocity",
        "ref_kind": "means",
        "ref_key": "umean",
        "scaling": "mean",
        "ylabel": r"$u_x^+$",
        "title": "Mean Streamwise Velocity",
    },
    "u2": {
        "group": "velocity",
        "ref_kind": "none",
        "ref_key": None,
        "scaling": "mean",
        "ylabel": r"$u_y^+$",
        "title": "Mean Wall-Normal Velocity",
    },
    "u3": {
        "group": "velocity",
        "ref_kind": "means",
        "ref_key": "wmean",
        "scaling": "mean",
        "ylabel": r"$u_z^+$",
        "title": "Mean Spanwise Velocity",
    },
    "p": {
        "group": "pressure",
        "ref_kind": "means",
        "ref_key": "pmean",
        "scaling": "pressure",
        "ylabel": r"$p^+$",
        "title": "Mean Pressure",
    },
    "uu": {
        "group": "stress",
        "ref_kind": "reystress",
        "ref_key": "uu",
        "scaling": "reynolds",
        "ylabel": r"$\overline{u^\prime u^\prime}^+$",
        "title": "Reynolds Stress $uu$",
    },
    "uv": {
        "group": "stress",
        "ref_kind": "reystress",
        "ref_key": "uv",
        "scaling": "reynolds",
        "ylabel": r"$\overline{u^\prime v^\prime}^+$",
        "title": "Reynolds Stress $uv$",
    },
    "uw": {
        "group": "stress",
        "ref_kind": "reystress",
        "ref_key": "uw",
        "scaling": "reynolds",
        "ylabel": r"$\overline{u^\prime w^\prime}^+$",
        "title": "Reynolds Stress $uw$",
    },
    "vv": {
        "group": "stress",
        "ref_kind": "reystress",
        "ref_key": "vv",
        "scaling": "reynolds",
        "ylabel": r"$\overline{v^\prime v^\prime}^+$",
        "title": "Reynolds Stress $vv$",
    },
    "vw": {
        "group": "stress",
        "ref_kind": "reystress",
        "ref_key": "vw",
        "scaling": "reynolds",
        "ylabel": r"$\overline{v^\prime w^\prime}^+$",
        "title": "Reynolds Stress $vw$",
    },
    "ww": {
        "group": "stress",
        "ref_kind": "reystress",
        "ref_key": "ww",
        "scaling": "reynolds",
        "ylabel": r"$\overline{w^\prime w^\prime}^+$",
        "title": "Reynolds Stress $ww$",
    },
    "omega_x": {
        "group": "vorticity",
        "ref_kind": "means",
        "ref_key": "wmeandy",
        "scaling": "vorticity",
        "ylabel": r"$\omega_x^+$",
        "title": "Vorticity $\\omega_x$",
    },
    "omega_y": {
        "group": "vorticity",
        "ref_kind": "means",
        "ref_key": None,
        "scaling": "vorticity",
        "ylabel": r"$\omega_y^+$",
        "title": "Vorticity $\\omega_y$",
    },
    "omega_z": {
        "group": "vorticity",
        "ref_kind": "means",
        "ref_key": "umeandy",
        "scaling": "vorticity",
        "ylabel": r"$\omega_z^+$",
        "title": "Vorticity $\\omega_z$",
    },
}

cbrg = cm.get_cmap("brg")
mlst = ["o", "<", "*", "v", "^", ">", "1", "2", "3", "4", "x", "s", "8", "+"]
plt.rc("figure", facecolor="white")
plt.rc("legend", fontsize=13)
rcParams["legend.loc"] = "best"

FIGSIZE = (9, 6)
DPI = 350
REF_HEADER_LINES = 25


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Plot CHAPSim2 channel velocity, pressure, Reynolds stress, and vorticity profiles."
    )
    parser.add_argument(
        "--dns-time",
        required=True,
        help="Time stamp for DNS data, e.g. 85000.",
    )
    parser.add_argument(
        "--re",
        type=float,
        default=2800.0,
        help="Channel Reynolds number used for wall scaling, e.g. 2800.",
    )
    parser.add_argument(
        "--ref-retau",
        type=int,
        default=None,
        help="Reference Re_tau to compare against, e.g. 180 or 395.",
    )
    parser.add_argument(
        "--groups",
        nargs="+",
        choices=("all", "velocity", "pressure", "stress", "vorticity"),
        default=("all",),
        help="Plot groups to generate. Default: all.",
    )
    parser.add_argument(
        "--input-dir",
        default="../1_data",
        help="Directory containing domain1_tsp_avg_*.dat files.",
    )
    parser.add_argument(
        "--output-dir",
        default=".",
        help="Directory for generated figures.",
    )
    return parser.parse_args()
def wall_distance(y: np.ndarray) -> np.ndarray:
    y_bottom = -1.0
    y_top = 1.0
    return np.minimum(y - y_bottom, y_top - y)


def derivative_y(values: np.ndarray, y: np.ndarray) -> np.ndarray:
    edge_order = 2 if len(y) >= 3 else 1
    return np.gradient(values, y, edge_order=edge_order)


class ChannelFlowPlotter:
    def __init__(self, args: argparse.Namespace):
        self.dns_time = args.dns_time
        self.re = args.re
        self.input_dir = (Path(__file__).resolve().parent / args.input_dir).resolve()
        self.output_dir = (Path(__file__).resolve().parent / args.output_dir).resolve()
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.ref_root = Path(__file__).resolve().parents[2]

        self.means_cache: dict[str, np.ndarray] = {}
        self.rey_cache: dict[str, np.ndarray] = {}
        self.dudx_cache: dict[str, np.ndarray] = {}

        self._load_wall_units()
        self.ref_retau = self._select_reference_retau(args.ref_retau)
        self.ref_cache: dict[str, np.ndarray] = {}

    def _load_ascii_column(self, stem: str) -> np.ndarray:
        path = self.input_dir / f"domain1_tsp_avg_{stem}_{self.dns_time}.dat"
        if not path.exists():
            raise FileNotFoundError(f"Missing required file: {path}")

        data = np.loadtxt(path)
        data = np.atleast_2d(data)
        if data.shape[1] < 3:
            raise ValueError(f"{path} must contain at least 3 columns: index, y, value")
        return data

    def _load_mean_profile(self, name: str) -> np.ndarray:
        if name not in self.means_cache:
            self.means_cache[name] = self._load_ascii_column(FILEMAP_MEAN[name])
        return self.means_cache[name]

    def _load_reynolds_profile(self, name: str) -> np.ndarray:
        if name not in self.rey_cache:
            self.rey_cache[name] = self._load_ascii_column(FILEMAP_REY[name])
        return self.rey_cache[name]

    def _load_dudx_profile(self, name: str) -> np.ndarray:
        if name not in self.dudx_cache:
            self.dudx_cache[name] = self._load_ascii_column(name)
        return self.dudx_cache[name]

    def _load_wall_units(self) -> None:
        ux = self._load_mean_profile("u1")
        self.grid_index = ux[:, 0]
        self.grid_y = ux[:, 1]

        y_wall = -1.0
        y1 = self.grid_y[0] - y_wall
        y2 = self.grid_y[1] - y_wall
        u1 = ux[0, 2]
        u2 = ux[1, 2]

        dudy_wall = (u1 * y2 * y2 - u2 * y1 * y1) / (y1 * y2 * (y2 - y1))
        self.tauw = dudy_wall / self.re
        self.utau = math.sqrt(abs(self.tauw))
        self.retau = self.re * self.utau
        self.yplus = self.retau * wall_distance(self.grid_y)

        print(f"Computed utau = {self.utau:.6f}, Re_tau = {self.retau:.2f}")

    def _available_reference_retaus(self) -> list[int]:
        available = []
        for entry in sorted(self.ref_root.iterdir()):
            if not (entry.is_dir() and entry.name.startswith("MKM") and entry.name.endswith("_profiles")):
                continue

            retau_str = entry.name[len("MKM") : -len("_profiles")]
            if not retau_str.isdigit():
                continue

            retau = int(retau_str)
            required = (
                entry / f"chan{retau}.means",
                entry / f"chan{retau}.reystress",
            )
            if all(path.exists() for path in required):
                available.append(retau)
        return available

    def _select_reference_retau(self, requested_retau: int | None) -> int:
        available = self._available_reference_retaus()
        if not available:
            raise FileNotFoundError(f"No MKM reference datasets found under {self.ref_root}")

        if requested_retau is not None:
            if requested_retau not in available:
                raise FileNotFoundError(
                    f"Requested reference Re_tau={requested_retau} is not available. Found: {available}"
                )
            selected = requested_retau
        else:
            selected = min(available, key=lambda retau: abs(retau - self.retau))

        print(
            f"Using reference Re_tau = {selected} "
            f"(simulation Re_tau = {self.retau:.2f})"
        )
        return selected

    def _load_reference(self, ref_kind: str) -> np.ndarray | None:
        if ref_kind == "none":
            return None
        if ref_kind in self.ref_cache:
            return self.ref_cache[ref_kind]

        folder = self.ref_root / f"MKM{self.ref_retau}_profiles"
        path = folder / f"chan{self.ref_retau}.{ref_kind}"

        if ref_kind == "means":
            names = ["y", "yplus", "umean", "umeandy", "wmean", "wmeandy", "pmean"]
        elif ref_kind == "reystress":
            names = ["y", "yplus", "uu", "vv", "ww", "uv", "uw", "vw"]
        else:
            raise ValueError(f"Unsupported reference kind: {ref_kind}")

        self.ref_cache[ref_kind] = np.genfromtxt(path, skip_header=REF_HEADER_LINES, names=names)
        return self.ref_cache[ref_kind]

    def _dns_series(self, name: str) -> tuple[np.ndarray, np.ndarray]:
        param = PARAMS[name]

        if param["group"] == "velocity":
            data = self._load_mean_profile(name)
            dns_val = data[:, 2] / self.utau
            return self.yplus, dns_val

        if param["group"] == "pressure":
            data = self._load_mean_profile("p")
            dns_val = data[:, 2] / self.tauw
            dns_val = dns_val - dns_val[0]
            return self.yplus, dns_val

        if param["group"] == "stress":
            stress = self._load_reynolds_profile(name)[:, 2]
            u1 = self._load_mean_profile("u1")[:, 2]
            u2 = self._load_mean_profile("u2")[:, 2]
            u3 = self._load_mean_profile("u3")[:, 2]
            mean_pair = {
                "uu": u1 * u1,
                "uv": u1 * u2,
                "uw": u1 * u3,
                "vv": u2 * u2,
                "vw": u2 * u3,
                "ww": u3 * u3,
            }[name]
            dns_val = (stress - mean_pair) / (self.utau * self.utau)
            return self.yplus, dns_val

        if param["group"] == "vorticity":
            omega = self._compute_dns_vorticity()[name]
            dns_val = omega / (self.re * self.utau * self.utau)
            return self.yplus, dns_val

        raise ValueError(name)

    def _compute_dns_vorticity(self) -> dict[str, np.ndarray]:
        try:
            dudx = {name: self._load_dudx_profile(name)[:, 2] for name in DUDX_COMPONENTS}
            return {
                "omega_x": dudx["dudx32"] - dudx["dudx23"],
                "omega_y": dudx["dudx13"] - dudx["dudx31"],
                "omega_z": dudx["dudx21"] - dudx["dudx12"],
            }
        except FileNotFoundError:
            u1 = self._load_mean_profile("u1")[:, 2]
            u3 = self._load_mean_profile("u3")[:, 2]
            return {
                "omega_x": derivative_y(u3, self.grid_y),
                "omega_y": np.zeros_like(self.grid_y),
                "omega_z": -derivative_y(u1, self.grid_y),
            }

    def _reference_series(self, name: str) -> tuple[np.ndarray, np.ndarray] | None:
        param = PARAMS[name]
        ref = self._load_reference(param["ref_kind"])
        if ref is None:
            return None

        if param["group"] == "vorticity":
            if name == "omega_x":
                ref_val = ref["wmeandy"] / self.ref_retau
            elif name == "omega_y":
                ref_val = np.zeros_like(ref["yplus"])
            elif name == "omega_z":
                ref_val = -ref["umeandy"] / self.ref_retau
            else:
                raise ValueError(name)
            return ref["yplus"], ref_val

        if param["ref_key"] is None:
            return None

        ref_val = ref[param["ref_key"]]

        if param["group"] == "pressure":
            ref_val = ref_val - ref_val[0]

        return ref["yplus"], ref_val

    def plot_quantity(self, name: str) -> None:
        param = PARAMS[name]
        fig, ax = plt.subplots(figsize=FIGSIZE, dpi=DPI)
        ax.set_xlabel(r"$y^+$", fontsize=18)
        ax.set_ylabel(param["ylabel"], fontsize=18)
        ax.set_title(param["title"], fontsize=16)
        ax.set_xscale("log")

        ref_series = self._reference_series(name)
        if ref_series is not None:
            ref_x, ref_y = ref_series
            ref_label = f"MKM{self.ref_retau}"
            ax.plot(
                ref_x,
                ref_y,
                marker=mlst[1],
                mfc="none",
                ms=4,
                color=cbrg(0.0),
                linestyle="None",
                label=ref_label,
            )

        dns_x, dns_y = self._dns_series(name)
        ax.plot(
            dns_x,
            dns_y,
            linestyle="--",
            color=cbrg(0.55),
            linewidth=1.6,
            label="CHAPSim2",
        )

        ax.grid(True, which="both", ls="-", alpha=0.2)
        ax.legend()
        ax.set_xlim(0.1, max(500.0, float(np.nanmax(dns_x))))

        outfile = self.output_dir / f"channel_{name}_{self.dns_time}.png"
        fig.savefig(outfile, bbox_inches="tight")
        plt.close(fig)
        print(f"Saved {outfile}")


def expand_groups(groups: tuple[str, ...] | list[str]) -> list[str]:
    if "all" in groups:
        ordered_groups = ("velocity", "pressure", "stress", "vorticity")
    else:
        ordered_groups = groups

    names: list[str] = []
    for group in ordered_groups:
        names.extend(PLOT_GROUPS[group])
    return names


def main() -> None:
    args = parse_args()
    names = expand_groups(args.groups)

    print(f"\nUsing DNS_TIME = {args.dns_time}")
    print(f"Input directory  = {Path(args.input_dir)}")
    print(f"Output directory = {Path(args.output_dir)}\n")

    plotter = ChannelFlowPlotter(args)

    for name in names:
        print(f"=== Processing {name} ===")
        plotter.plot_quantity(name)

    print("\nAll requested plots completed.\n")


if __name__ == "__main__":
    main()
