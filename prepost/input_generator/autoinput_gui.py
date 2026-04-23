#!/usr/bin/env python3
"""CHAPSim2 GUI Auto Input Generator - User-friendly configuration interface."""

import tkinter as tk
from tkinter import ttk, messagebox
import configparser
import math
from enum import Enum

# Constants
PI = round(math.pi, 6)
TWO_PI = 2.0 * PI
DEFAULT_FILENAME = "input_chapsim_gui.ini"


class Case(Enum):
    CHANNEL = 1
    PIPE = 2
    ANNULAR = 3
    TGV3D = 4
    DUCT = 5


class Init(Enum):
    RESTART = 0
    INTRPL = 1
    RANDOM = 2
    INLET = 3
    GIVEN = 4
    POISEUILLE = 5
    FUNCTION = 6
    GVBCLN = 7


class Stretching(Enum):
    NONE = 0
    CENTRE = 1
    SIDE2 = 2
    BOTTOM = 3
    TOP = 4


class BC(Enum):
    INTERIOR = 0
    PERIODIC = 1
    SYMM = 2
    ASYMM = 3
    DIRICHLET = 4
    NEUMANN = 5
    INTRPL = 6
    CONVOL = 7
    TURGEN = 8
    PROFL = 9
    DATABS = 10
    PARABOLIC = 11
    OTHERS = 12


class Drvfc(Enum):
    NONE = 0
    XMFLUX = 1
    XTAUW = 2
    XDPDX = 3
    ZMFLUX = 4
    ZTAUW = 5
    ZDPDZ = 6


def bool_to_string(value):
    """Converts 0/1 to Fortran boolean strings."""
    return ".true." if value else ".false."


class CustomConfigParser(configparser.ConfigParser):
    """Custom ConfigParser that formats output with space after '='."""

    def write(self, fp):
        for section in self.sections():
            fp.write(f"[{section}]\n")
            for key, value in self.items(section):
                fp.write(f"{key}= {value}\n")
            fp.write("\n")


class CHAPSimGUI:
    def __init__(self, root):
        self.root = root
        self.root.title("CHAPSim2 Input Generator")
        self.root.geometry("1100x850")

        # Create notebook for tabs
        self.notebook = ttk.Notebook(root)
        self.notebook.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)

        # State variables
        self.icase = tk.IntVar(value=1)
        self.ithermo = tk.BooleanVar(value=False)
        self.imhd = tk.BooleanVar(value=False)
        self.is_restart = tk.BooleanVar(value=False)
        self.is_write = tk.BooleanVar(value=False)
        self.iinlet = tk.IntVar(value=1)  # 1 = periodic (default)

        # Create tabs (order matters - mesh and bc must exist before domain triggers updates)
        self.create_process_tab()
        self.create_decomp_tab()
        self.create_flow_tab()
        self.create_thermo_tab()
        self.create_mhd_tab()
        self.create_mesh_tab()
        self.create_bc_tab()
        self.create_scheme_tab()
        self.create_simcontrol_tab()
        self.create_io_tab()
        self.create_probe_tab()
        # Create domain tab last so it can trigger updates to other tabs
        self.create_domain_tab()

        # Generate button
        button_frame = ttk.Frame(root)
        button_frame.pack(fill=tk.X, padx=10, pady=10)

        ttk.Button(
            button_frame, text="Generate INI File", command=self.generate_ini
        ).pack(side=tk.LEFT, padx=5)
        ttk.Button(button_frame, text="Exit", command=root.quit).pack(
            side=tk.RIGHT, padx=5
        )

    def create_labeled_input(
        self, parent, label_text, default_value, row, col=0, input_type="str", width=20
    ):
        """Helper to create labeled input fields."""
        if label_text:  # Only create label if text provided
            ttk.Label(parent, text=label_text).grid(
                row=row, column=col, sticky=tk.W, padx=5, pady=5
            )

        if input_type == "bool":
            var = tk.BooleanVar(value=int(default_value))
            widget = ttk.Checkbutton(parent, variable=var)
            widget.grid(row=row, column=col + 1, sticky=tk.W, padx=5, pady=5)
            return (var, widget)
        elif input_type == "choice":
            var = tk.StringVar(value=str(default_value))
            widget = ttk.Combobox(parent, textvariable=var, state="readonly", width=width)
            widget.grid(row=row, column=col + 1, sticky=tk.W, padx=5, pady=5)
            return (var, widget)
        else:
            var = tk.StringVar(value=str(default_value))
            widget = ttk.Entry(parent, textvariable=var, width=width)
            widget.grid(row=row, column=col + 1, sticky=tk.W, padx=5, pady=5)
            return (var, widget)

    def set_enabled(self, widget, enabled):
        """Enable or disable a widget."""
        if widget:
            state = "normal" if enabled else "disabled"
            if isinstance(widget, ttk.Entry):
                widget.config(state=state)
            elif isinstance(widget, ttk.Combobox):
                widget.config(state="readonly" if enabled else "disabled")
            elif isinstance(widget, ttk.Checkbutton):
                widget.config(state=state)

    def create_process_tab(self):
        """Process settings tab."""
        tab = ttk.Frame(self.notebook)
        self.notebook.add(tab, text="Process")

        self.is_prerun, _ = self.create_labeled_input(
            tab, "Enable prerun only?", 0, 0, 0, "bool"
        )
        self.is_postprocess, _ = self.create_labeled_input(
            tab, "Enable postprocess?", 0, 1, 0, "bool"
        )

    def create_decomp_tab(self):
        """Decomposition settings tab."""
        tab = ttk.Frame(self.notebook)
        self.notebook.add(tab, text="Decomposition")

        self.is_decomp = tk.BooleanVar(value=True)
        ttk.Label(tab, text="Auto domain decomposition?").grid(
            row=0, column=0, sticky=tk.W, padx=5, pady=5
        )
        self.is_decomp_check = ttk.Checkbutton(
            tab, variable=self.is_decomp, command=self.on_decomp_changed
        )
        self.is_decomp_check.grid(row=0, column=1, sticky=tk.W, padx=5, pady=5)

        self.p_row, self.p_row_w = self.create_labeled_input(
            tab, "Subdomain division Y", 0, 1
        )
        self.p_col, self.p_col_w = self.create_labeled_input(
            tab, "Subdomain division Z", 0, 2
        )

        self.on_decomp_changed()

    def on_decomp_changed(self, *args):
        """Handle decomposition checkbox changes."""
        enabled = not self.is_decomp.get()
        self.set_enabled(self.p_row_w, enabled)
        self.set_enabled(self.p_col_w, enabled)

    def create_domain_tab(self):
        """Domain settings tab."""
        tab = ttk.Frame(self.notebook)
        self.notebook.add(tab, text="Domain")

        cases = ["1:Channel", "2:Pipe", "3:Annular", "4:TGV3D", "5:DUCT"]
        self.icase_var, self.icase_combo = self.create_labeled_input(
            tab, "Simulation case", "1:Channel", 0, 0, "choice"
        )
        self.icase_combo["values"] = cases
        self.icase_combo.bind("<<ComboboxSelected>>", lambda e: self.on_case_changed())

        self.lxx, self.lxx_w = self.create_labeled_input(
            tab, "Streamwise length (Lx/h)", TWO_PI, 1
        )
        self.lzz, self.lzz_w = self.create_labeled_input(
            tab, "Spanwise length (Lz/h)", PI, 2
        )
        self.lyb, self.lyb_w = self.create_labeled_input(
            tab, "Bottom boundary", -1.0, 3
        )
        self.lyt, self.lyt_w = self.create_labeled_input(tab, "Top boundary", 1.0, 4)

        # Trigger initial update after all widgets are created
        self.root.after(100, self.on_case_changed)

    def on_case_changed(self):
        """Handle case selection changes - updates all dependent defaults."""
        case_str = self.icase_var.get()
        case_num = int(case_str.split(":")[0])
        self.icase.set(case_num)

        # Enable all first
        self.set_enabled(self.lxx_w, True)
        self.set_enabled(self.lzz_w, True)
        self.set_enabled(self.lyb_w, True)
        self.set_enabled(self.lyt_w, True)

        # Set defaults based on case
        if case_num == Case.CHANNEL.value:
            self.lxx.set(str(TWO_PI))
            self.lzz.set(str(PI))
            self.lyt.set("1.0")
            self.lyb.set("-1.0")
            self.set_enabled(self.lyb_w, False)
            self.set_enabled(self.lyt_w, False)
        elif case_num == Case.PIPE.value:
            self.lxx.set(str(TWO_PI))
            self.lzz.set(str(TWO_PI))
            self.lyt.set("1.0")
            self.lyb.set("0.0")
            self.set_enabled(self.lzz_w, False)
            self.set_enabled(self.lyb_w, False)
            self.set_enabled(self.lyt_w, False)
        elif case_num == Case.ANNULAR.value:
            self.lxx.set(str(TWO_PI))
            self.lzz.set(str(TWO_PI))
            self.lyt.set("1.0")
            self.set_enabled(self.lzz_w, False)
            self.set_enabled(self.lyt_w, False)
        elif case_num == Case.TGV3D.value:
            self.lxx.set(str(TWO_PI))
            self.lzz.set(str(TWO_PI))
            self.lyb.set(str(-PI))
            self.lyt.set(str(PI))
            self.set_enabled(self.lxx_w, False)
            self.set_enabled(self.lzz_w, False)
            self.set_enabled(self.lyb_w, False)
            self.set_enabled(self.lyt_w, False)
        elif case_num == Case.DUCT.value:
            self.lxx.set("2.0")
            self.lzz.set("12.0")
            self.lyb.set("-1.0")
            self.lyt.set("1.0")
            self.set_enabled(self.lyb_w, False)
            self.set_enabled(self.lyt_w, False)

        # Update mesh stretching defaults
        self.update_mesh_defaults()
        
        # Update BC defaults
        self.update_bc_defaults()
        
        # Update flow defaults
        self.update_flow_defaults()

    def create_flow_tab(self):
        """Flow settings tab."""
        tab = ttk.Frame(self.notebook)
        self.notebook.add(tab, text="Flow")

        ttk.Label(tab, text="Flow restart?").grid(
            row=0, column=0, sticky=tk.W, padx=5, pady=5
        )
        self.is_restart_check = ttk.Checkbutton(
            tab, variable=self.is_restart, command=self.on_restart_changed
        )
        self.is_restart_check.grid(row=0, column=1, sticky=tk.W, padx=5, pady=5)

        self.irestartfrom, self.irestartfrom_w = self.create_labeled_input(
            tab, "Restart from iteration", 2000, 1
        )

        self.velo1, self.velo1_w = self.create_labeled_input(
            tab, "Initial velocity X", 0.0, 2
        )
        self.velo2, self.velo2_w = self.create_labeled_input(
            tab, "Initial velocity Y", 0.0, 3
        )
        self.velo3, self.velo3_w = self.create_labeled_input(
            tab, "Initial velocity Z", 0.0, 4
        )
        self.noiselevel, self.noiselevel_w = self.create_labeled_input(
            tab, "Noise level (0-1)", 0.25, 5
        )
        self.ren, self.ren_w = self.create_labeled_input(
            tab, "Reynolds number", 2800, 6
        )
        self.reni, self.reni_w = self.create_labeled_input(
            tab, "Initial Reynolds number", 20000, 7
        )
        self.nreni, self.nreni_w = self.create_labeled_input(
            tab, "Iterations for initial Re", 10000, 8
        )

        self.on_restart_changed()

    def update_flow_defaults(self):
        """Update flow defaults based on case."""
        # Safety check - ensure flow widgets exist
        if not hasattr(self, 'noiselevel'):
            return
            
        case_num = self.icase.get()
        
        # Update noise level visibility
        if case_num == Case.TGV3D.value:
            self.noiselevel.set("0.0")
            self.set_enabled(self.noiselevel_w, False)
        else:
            self.noiselevel.set("0.25")
            self.set_enabled(self.noiselevel_w, not self.is_restart.get())
        
        # Update reni/nreni for TGV3D
        if case_num == Case.TGV3D.value:
            ren_val = self.ren.get()
            self.reni.set(ren_val)
            self.nreni.set("0")
            self.set_enabled(self.reni_w, False)
            self.set_enabled(self.nreni_w, False)
        else:
            if self.reni.get() == "0" or self.reni.get() == self.ren.get():
                self.reni.set("20000")
            if self.nreni.get() == "0":
                self.nreni.set("10000")
            self.set_enabled(self.reni_w, not self.is_restart.get())
            self.set_enabled(self.nreni_w, not self.is_restart.get())

    def on_restart_changed(self):
        """Handle restart checkbox changes."""
        enabled = self.is_restart.get()
        self.set_enabled(self.irestartfrom_w, enabled)
        self.set_enabled(self.velo1_w, not enabled)
        self.set_enabled(self.velo2_w, not enabled)
        self.set_enabled(self.velo3_w, not enabled)
        
        case_num = self.icase.get()
        if case_num != Case.TGV3D.value:
            self.set_enabled(self.noiselevel_w, not enabled)
            self.set_enabled(self.reni_w, not enabled)
            self.set_enabled(self.nreni_w, not enabled)

    def create_thermo_tab(self):
        """Thermal settings tab."""
        tab = ttk.Frame(self.notebook)
        self.notebook.add(tab, text="Thermo")

        ttk.Label(tab, text="Enable thermal field?").grid(
            row=0, column=0, sticky=tk.W, padx=5, pady=5
        )
        self.ithermo_check = ttk.Checkbutton(
            tab, variable=self.ithermo, command=self.on_thermo_changed
        )
        self.ithermo_check.grid(row=0, column=1, sticky=tk.W, padx=5, pady=5)

        self.icht, self.icht_w = self.create_labeled_input(
            tab, "Conjugate heat transfer?", 0, 1, 0, "bool"
        )

        gravity_options = ["0:None", "1:+X", "-1:-X", "2:+Y", "-2:-Y", "3:+Z", "-3:-Z"]
        self.igravity, self.igravity_w = self.create_labeled_input(
            tab, "Gravity direction", "0:None", 2, 0, "choice"
        )
        self.igravity_w["values"] = gravity_options

        fluid_options = [
            "1:scp-H2O",
            "2:scp-CO2",
            "3:sodium",
            "4:lead",
            "5:bismuth",
            "6:LBE",
        ]
        self.ifluid, self.ifluid_w = self.create_labeled_input(
            tab, "Fluid type", "1:scp-H2O", 3, 0, "choice"
        )
        self.ifluid_w["values"] = fluid_options

        self.refl0, self.refl0_w = self.create_labeled_input(
            tab, "Reference length (m)", 0.001, 4
        )
        self.refT0, self.refT0_w = self.create_labeled_input(
            tab, "Reference temperature (K)", 645.15, 5
        )

        inittm_options = [
            "0:Restart",
            "1:Intrpl",
            "2:Random",
            "3:Inlet",
            "4:Given",
            "5:Poiseuille",
            "6:Function",
            "7:GivenBCMix",
        ]
        self.inittm, self.inittm_w = self.create_labeled_input(
            tab, "Thermal initialization", "4:Given", 6, 0, "choice"
        )
        self.inittm_w["values"] = inittm_options

        self.Tini, self.Tini_w = self.create_labeled_input(
            tab, "Initial temperature (K)", 645.15, 7
        )

        self.on_thermo_changed()

    def on_thermo_changed(self):
        """Handle thermal checkbox changes."""
        enabled = self.ithermo.get()
        self.set_enabled(self.icht_w, enabled)
        self.set_enabled(self.igravity_w, enabled)
        self.set_enabled(self.ifluid_w, enabled)
        self.set_enabled(self.refl0_w, enabled)
        self.set_enabled(self.refT0_w, enabled)
        self.set_enabled(self.inittm_w, enabled)
        self.set_enabled(self.Tini_w, enabled)

    def create_mhd_tab(self):
        """MHD settings tab."""
        tab = ttk.Frame(self.notebook)
        self.notebook.add(tab, text="MHD")

        ttk.Label(tab, text="Enable MHD?").grid(
            row=0, column=0, sticky=tk.W, padx=5, pady=5
        )
        self.imhd_check = ttk.Checkbutton(
            tab, variable=self.imhd, command=self.on_mhd_changed
        )
        self.imhd_check.grid(row=0, column=1, sticky=tk.W, padx=5, pady=5)

        mhd_type = ["1:Stuart", "2:Hartmann"]
        self.mhd_type, self.mhd_type_w = self.create_labeled_input(
            tab, "MHD number type", "2:Hartmann", 1, 0, "choice"
        )
        self.mhd_type_w["values"] = mhd_type

        self.NS, self.NS_w = self.create_labeled_input(tab, "Stuart number", 10.0, 2)
        self.NH, self.NH_w = self.create_labeled_input(tab, "Hartmann number", 10.0, 3)
        self.b1, self.b1_w = self.create_labeled_input(tab, "Magnetic field X", 0.0, 4)
        self.b2, self.b2_w = self.create_labeled_input(tab, "Magnetic field Y", 1.0, 5)
        self.b3, self.b3_w = self.create_labeled_input(tab, "Magnetic field Z", 0.0, 6)

        self.on_mhd_changed()

    def on_mhd_changed(self):
        """Handle MHD checkbox changes."""
        enabled = self.imhd.get()
        self.set_enabled(self.mhd_type_w, enabled)
        self.set_enabled(self.NS_w, enabled)
        self.set_enabled(self.NH_w, enabled)
        self.set_enabled(self.b1_w, enabled)
        self.set_enabled(self.b2_w, enabled)
        self.set_enabled(self.b3_w, enabled)

    def create_mesh_tab(self):
        """Mesh settings tab."""
        tab = ttk.Frame(self.notebook)
        self.notebook.add(tab, text="Mesh")

        self.ncx, self.ncx_w = self.create_labeled_input(tab, "Cells in X", 64, 0)
        self.ncy, self.ncy_w = self.create_labeled_input(tab, "Cells in Y", 64, 1)
        self.ncz, self.ncz_w = self.create_labeled_input(tab, "Cells in Z", 64, 2)

        stretch_options = ["0:None", "1:Centre", "2:2-sides", "3:Bottom", "4:Top"]
        self.istret, self.istret_w = self.create_labeled_input(
            tab, "Grid clustering", "2:2-sides", 3, 0, "choice"
        )
        self.istret_w["values"] = stretch_options
        self.istret_w.bind(
            "<<ComboboxSelected>>", lambda e: self.on_stretching_changed()
        )

        self.rstret1, self.rstret1_w = self.create_labeled_input(
            tab, "Stretching method", 1, 4
        )
        self.rstret2, self.rstret2_w = self.create_labeled_input(
            tab, "Stretching factor", 0.12, 5
        )

        self.on_stretching_changed()

    def update_mesh_defaults(self):
        """Update mesh stretching defaults based on case."""
        # Safety check - ensure mesh widgets exist
        if not hasattr(self, 'istret'):
            return
            
        case_num = self.icase.get()
        
        if case_num in [Case.CHANNEL.value, Case.DUCT.value, Case.ANNULAR.value]:
            self.istret.set("2:2-sides")
            self.rstret1.set("1")
            self.rstret2.set("0.12")
        elif case_num == Case.PIPE.value:
            self.istret.set("4:Top")
            self.rstret1.set("2")
            self.rstret2.set("0.15")
        elif case_num == Case.TGV3D.value:
            self.istret.set("0:None")
            self.rstret1.set("0")
            self.rstret2.set("0.0")
        
        self.on_stretching_changed()

    def on_stretching_changed(self):
        """Handle stretching selection changes."""
        istret_val = int(self.istret.get().split(":")[0])
        enabled = istret_val != 0
        self.set_enabled(self.rstret1_w, enabled)
        self.set_enabled(self.rstret2_w, enabled)
        
        # Update stretching method defaults based on case
        if enabled:
            case_num = self.icase.get()
            if case_num in [Case.CHANNEL.value, Case.DUCT.value]:
                if self.rstret1.get() not in ["1", "2", "3"]:
                    self.rstret1.set("1")
                if float(self.rstret2.get()) < 0.01:
                    self.rstret2.set("0.12")
            elif case_num in [Case.PIPE.value, Case.ANNULAR.value]:
                if self.rstret1.get() not in ["1", "2", "3"]:
                    self.rstret1.set("2")
                if float(self.rstret2.get()) < 0.01:
                    self.rstret2.set("0.15")

    def create_bc_tab(self):
        """Boundary condition settings tab."""
        tab = ttk.Frame(self.notebook)
        self.notebook.add(tab, text="BC")

        # Create a scrollable frame
        canvas = tk.Canvas(tab)
        scrollbar = ttk.Scrollbar(tab, orient="vertical", command=canvas.yview)
        scrollable_frame = ttk.Frame(canvas)

        scrollable_frame.bind(
            "<Configure>", lambda e: canvas.configure(scrollregion=canvas.bbox("all"))
        )

        canvas.create_window((0, 0), window=scrollable_frame, anchor="nw")
        canvas.configure(yscrollcommand=scrollbar.set)

        # BC type options
        bc_options = [
            "0:Interior",
            "1:Periodic",
            "2:Symm",
            "3:Asymm",
            "4:Dirichlet",
            "5:Neumann",
            "6:Intrpl",
            "7:ConvOut",
            "8:TurGen",
            "9:Profile",
            "10:Database",
            "11:Parabolic",
            "12:Others"
        ]

        # Store all BC widgets in dictionary
        self.bc_widgets = {}
        row = 0

        # Headers
        ttk.Label(scrollable_frame, text="Boundary", font=("Arial", 10, "bold")).grid(
            row=row, column=0, padx=5, pady=5
        )
        ttk.Label(scrollable_frame, text="BC Type 1", font=("Arial", 10, "bold")).grid(
            row=row, column=1, padx=5, pady=5
        )
        ttk.Label(scrollable_frame, text="BC Type 2", font=("Arial", 10, "bold")).grid(
            row=row, column=2, padx=5, pady=5
        )
        ttk.Label(scrollable_frame, text="Value 1", font=("Arial", 10, "bold")).grid(
            row=row, column=3, padx=5, pady=5
        )
        ttk.Label(scrollable_frame, text="Value 2", font=("Arial", 10, "bold")).grid(
            row=row, column=4, padx=5, pady=5
        )
        row += 1

        # X-direction BC
        ttk.Label(scrollable_frame, text="X-direction (ifbcx_u/v/w)").grid(
            row=row, column=0, sticky=tk.W, padx=5, pady=5
        )
        self.bc_widgets["ifbcx_u_1"], _ = self._create_bc_combo(
            scrollable_frame, bc_options, "1", row, 1
        )
        self.bc_widgets["ifbcx_u_2"], _ = self._create_bc_combo(
            scrollable_frame, bc_options, "1", row, 2
        )
        self.bc_widgets["ifbcx_u_v1"], _ = self.create_labeled_input(
            scrollable_frame, "", "0.0", row, 3, "str", width=10
        )
        self.bc_widgets["ifbcx_u_v2"], _ = self.create_labeled_input(
            scrollable_frame, "", "0.0", row, 4, "str", width=10
        )
        row += 1

        ttk.Label(scrollable_frame, text="X-direction (ifbcx_p)").grid(
            row=row, column=0, sticky=tk.W, padx=5, pady=5
        )
        self.bc_widgets["ifbcx_p_1"], _ = self._create_bc_combo(
            scrollable_frame, bc_options, "1", row, 1
        )
        self.bc_widgets["ifbcx_p_2"], _ = self._create_bc_combo(
            scrollable_frame, bc_options, "1", row, 2
        )
        self.bc_widgets["ifbcx_p_v1"], _ = self.create_labeled_input(
            scrollable_frame, "", "0.0", row, 3, "str", width=10
        )
        self.bc_widgets["ifbcx_p_v2"], _ = self.create_labeled_input(
            scrollable_frame, "", "0.0", row, 4, "str", width=10
        )
        row += 1

        ttk.Label(scrollable_frame, text="X-direction (ifbcx_T)").grid(
            row=row, column=0, sticky=tk.W, padx=5, pady=5
        )
        self.bc_widgets["ifbcx_T_1"], _ = self._create_bc_combo(
            scrollable_frame, bc_options, "1", row, 1
        )
        self.bc_widgets["ifbcx_T_2"], _ = self._create_bc_combo(
            scrollable_frame, bc_options, "1", row, 2
        )
        self.bc_widgets["ifbcx_T_v1"], _ = self.create_labeled_input(
            scrollable_frame, "", "0.0", row, 3, "str", width=10
        )
        self.bc_widgets["ifbcx_T_v2"], _ = self.create_labeled_input(
            scrollable_frame, "", "0.0", row, 4, "str", width=10
        )
        row += 1

        # Y-direction BC
        ttk.Label(scrollable_frame, text="Y-direction (ifbcy_u/v/w)").grid(
            row=row, column=0, sticky=tk.W, padx=5, pady=5
        )
        self.bc_widgets["ifbcy_u_1"], _ = self._create_bc_combo(
            scrollable_frame, bc_options, "1", row, 1
        )
        self.bc_widgets["ifbcy_u_2"], _ = self._create_bc_combo(
            scrollable_frame, bc_options, "1", row, 2
        )
        self.bc_widgets["ifbcy_u_v1"], _ = self.create_labeled_input(
            scrollable_frame, "", "0.0", row, 3, "str", width=10
        )
        self.bc_widgets["ifbcy_u_v2"], _ = self.create_labeled_input(
            scrollable_frame, "", "0.0", row, 4, "str", width=10
        )
        row += 1

        ttk.Label(scrollable_frame, text="Y-direction (ifbcy_p)").grid(
            row=row, column=0, sticky=tk.W, padx=5, pady=5
        )
        self.bc_widgets["ifbcy_p_1"], _ = self._create_bc_combo(
            scrollable_frame, bc_options, "1", row, 1
        )
        self.bc_widgets["ifbcy_p_2"], _ = self._create_bc_combo(
            scrollable_frame, bc_options, "1", row, 2
        )
        self.bc_widgets["ifbcy_p_v1"], _ = self.create_labeled_input(
            scrollable_frame, "", "0.0", row, 3, "str", width=10
        )
        self.bc_widgets["ifbcy_p_v2"], _ = self.create_labeled_input(
            scrollable_frame, "", "0.0", row, 4, "str", width=10
        )
        row += 1

        ttk.Label(scrollable_frame, text="Y-direction (ifbcy_T)").grid(
            row=row, column=0, sticky=tk.W, padx=5, pady=5
        )
        self.bc_widgets["ifbcy_T_1"], _ = self._create_bc_combo(
            scrollable_frame, bc_options, "1", row, 1
        )
        self.bc_widgets["ifbcy_T_2"], _ = self._create_bc_combo(
            scrollable_frame, bc_options, "1", row, 2
        )
        self.bc_widgets["ifbcy_T_v1"], _ = self.create_labeled_input(
            scrollable_frame, "", "0.0", row, 3, "str", width=10
        )
        self.bc_widgets["ifbcy_T_v2"], _ = self.create_labeled_input(
            scrollable_frame, "", "0.0", row, 4, "str", width=10
        )
        row += 1

        # Z-direction BC
        ttk.Label(scrollable_frame, text="Z-direction (ifbcz_u/v/w)").grid(
            row=row, column=0, sticky=tk.W, padx=5, pady=5
        )
        self.bc_widgets["ifbcz_u_1"], _ = self._create_bc_combo(
            scrollable_frame, bc_options, "1", row, 1
        )
        self.bc_widgets["ifbcz_u_2"], _ = self._create_bc_combo(
            scrollable_frame, bc_options, "1", row, 2
        )
        self.bc_widgets["ifbcz_u_v1"], _ = self.create_labeled_input(
            scrollable_frame, "", "0.0", row, 3, "str", width=10
        )
        self.bc_widgets["ifbcz_u_v2"], _ = self.create_labeled_input(
            scrollable_frame, "", "0.0", row, 4, "str", width=10
        )
        row += 1

        ttk.Label(scrollable_frame, text="Z-direction (ifbcz_p)").grid(
            row=row, column=0, sticky=tk.W, padx=5, pady=5
        )
        self.bc_widgets["ifbcz_p_1"], _ = self._create_bc_combo(
            scrollable_frame, bc_options, "1", row, 1
        )
        self.bc_widgets["ifbcz_p_2"], _ = self._create_bc_combo(
            scrollable_frame, bc_options, "1", row, 2
        )
        self.bc_widgets["ifbcz_p_v1"], _ = self.create_labeled_input(
            scrollable_frame, "", "0.0", row, 3, "str", width=10
        )
        self.bc_widgets["ifbcz_p_v2"], _ = self.create_labeled_input(
            scrollable_frame, "", "0.0", row, 4, "str", width=10
        )
        row += 1

        ttk.Label(scrollable_frame, text="Z-direction (ifbcz_T)").grid(
            row=row, column=0, sticky=tk.W, padx=5, pady=5
        )
        self.bc_widgets["ifbcz_T_1"], _ = self._create_bc_combo(
            scrollable_frame, bc_options, "1", row, 1
        )
        self.bc_widgets["ifbcz_T_2"], _ = self._create_bc_combo(
            scrollable_frame, bc_options, "1", row, 2
        )
        self.bc_widgets["ifbcz_T_v1"], _ = self.create_labeled_input(
            scrollable_frame, "", "0.0", row, 3, "str", width=10
        )
        self.bc_widgets["ifbcz_T_v2"], _ = self.create_labeled_input(
            scrollable_frame, "", "0.0", row, 4, "str", width=10
        )
        row += 2

        # Inlet BC type
        ttk.Label(scrollable_frame, text="Streamwise BC type").grid(
            row=row, column=0, sticky=tk.W, padx=5, pady=5
        )
        inlet_options = ["1:Periodic", "4:Dirichlet", "9:Profile", "10:Database"]
        self.iinlet_var, self.iinlet_w = self._create_bc_combo(
            scrollable_frame, inlet_options, "1", row, 1
        )
        row += 1

        # Flow driving options
        ttk.Label(scrollable_frame, text="Flow driving method").grid(
            row=row, column=0, sticky=tk.W, padx=5, pady=5
        )
        driven_options = [
            "0:NONE",
            "1:XMFLUX",
            "2:XTAUW",
            "3:XDPDX",
            "4:ZMFLUX",
            "5:ZTAUW",
            "6:ZDPDZ"
        ]
        self.idriven, self.idriven_w = self._create_bc_combo(
            scrollable_frame, driven_options, "1:XMFLUX", row, 1
        )
        row += 1

        ttk.Label(scrollable_frame, text="Driving force magnitude").grid(
            row=row, column=0, sticky=tk.W, padx=5, pady=5
        )
        self.drivenCf, self.drivenCf_w = self.create_labeled_input(
            scrollable_frame, "", "0.0", row, 1, "str", width=20
        )

        canvas.pack(side="left", fill="both", expand=True)
        scrollbar.pack(side="right", fill="y")

    def update_bc_defaults(self):
        """Update BC defaults based on case type."""
        # Safety check - ensure BC widgets exist
        if not hasattr(self, 'bc_widgets') or not self.bc_widgets:
            return
            
        case_num = self.icase.get()
        
        # Y-direction defaults
        if case_num in [Case.CHANNEL.value, Case.DUCT.value, Case.ANNULAR.value]:
            self.bc_widgets["ifbcy_u_1"].set("4:Dirichlet")
            self.bc_widgets["ifbcy_u_2"].set("4:Dirichlet")
            self.bc_widgets["ifbcy_p_1"].set("5:Neumann")
            self.bc_widgets["ifbcy_p_2"].set("5:Neumann")
        elif case_num == Case.PIPE.value:
            self.bc_widgets["ifbcy_u_1"].set("0:Interior")
            self.bc_widgets["ifbcy_u_2"].set("4:Dirichlet")
            self.bc_widgets["ifbcy_p_1"].set("0:Interior")
            self.bc_widgets["ifbcy_p_2"].set("5:Neumann")
            self.bc_widgets["ifbcy_T_1"].set("0:Interior")
            self.bc_widgets["ifbcy_T_2"].set("1:Periodic")
        elif case_num == Case.TGV3D.value:
            self.bc_widgets["ifbcy_u_1"].set("1:Periodic")
            self.bc_widgets["ifbcy_u_2"].set("1:Periodic")
            self.bc_widgets["ifbcy_p_1"].set("1:Periodic")
            self.bc_widgets["ifbcy_p_2"].set("1:Periodic")
        
        # X-direction defaults
        if case_num == Case.DUCT.value:
            self.bc_widgets["ifbcx_u_1"].set("4:Dirichlet")
            self.bc_widgets["ifbcx_u_2"].set("4:Dirichlet")
            self.bc_widgets["ifbcx_p_1"].set("5:Neumann")
            self.bc_widgets["ifbcx_p_2"].set("5:Neumann")
            self.bc_widgets["ifbcx_T_1"].set("5:Neumann")
            self.bc_widgets["ifbcx_T_2"].set("5:Neumann")
            self.idriven.set("4:ZMFLUX")
        else:
            self.bc_widgets["ifbcx_u_1"].set("1:Periodic")
            self.bc_widgets["ifbcx_u_2"].set("1:Periodic")
            self.bc_widgets["ifbcx_p_1"].set("1:Periodic")
            self.bc_widgets["ifbcx_p_2"].set("1:Periodic")
            if case_num != Case.TGV3D.value:
                self.idriven.set("1:XMFLUX")
            else:
                self.idriven.set("0:NONE")
        
        # Z-direction defaults (mostly periodic)
        self.bc_widgets["ifbcz_u_1"].set("1:Periodic")
        self.bc_widgets["ifbcz_u_2"].set("1:Periodic")
        self.bc_widgets["ifbcz_p_1"].set("1:Periodic")
        self.bc_widgets["ifbcz_p_2"].set("1:Periodic")
        self.bc_widgets["ifbcz_T_1"].set("1:Periodic")
        self.bc_widgets["ifbcz_T_2"].set("1:Periodic")

    def _create_bc_combo(self, parent, options, default, row, col):
        """Helper to create BC combobox."""
        var = tk.StringVar(value=f"{default}:Periodic" if ":" not in str(default) else str(default))
        widget = ttk.Combobox(parent, textvariable=var, state="readonly", width=15)
        widget["values"] = options
        widget.grid(row=row, column=col, padx=5, pady=5)
        return var, widget

    def create_scheme_tab(self):
        """Numerical scheme settings tab."""
        tab = ttk.Frame(self.notebook)
        self.notebook.add(tab, text="Scheme")

        self.dt, _ = self.create_labeled_input(tab, "Time step size", 0.00001, 0)
        
        accuracy_options = ["1:2nd CD", "2:4th CD", "3:4th CP", "4:6th CP"]
        self.iAccuracy, _ = self.create_labeled_input(
            tab, "Spatial accuracy", "1:2nd CD", 1, 0, "choice"
        )
        self.iAccuracy.set("1:2nd CD")
        # Need to get the widget to set values
        for child in tab.winfo_children():
            if isinstance(child, ttk.Combobox) and child.cget("textvariable") == str(self.iAccuracy):
                child["values"] = accuracy_options
                break

    def create_simcontrol_tab(self):
        """Simulation control settings tab."""
        tab = ttk.Frame(self.notebook)
        self.notebook.add(tab, text="SimControl")

        self.nIterFlowFirst, _ = self.create_labeled_input(
            tab, "First iteration for flow", 1, 0
        )
        self.nIterFlowLast, _ = self.create_labeled_input(
            tab, "Last iteration for flow", 1000000, 1
        )
        self.nIterThermoFirst, _ = self.create_labeled_input(
            tab, "First iteration for thermal", 1, 2
        )
        self.nIterThermoLast, _ = self.create_labeled_input(
            tab, "Last iteration for thermal", 1000000, 3
        )

    def create_io_tab(self):
        """I/O settings tab."""
        tab = ttk.Frame(self.notebook)
        self.notebook.add(tab, text="I/O")

        self.cpu_nfre, _ = self.create_labeled_input(
            tab, "CPU info print frequency", 1, 0
        )
        self.ckpt_nfre, _ = self.create_labeled_input(
            tab, "Checkpoint save frequency", 1000, 1
        )
        self.visu_nfre, _ = self.create_labeled_input(
            tab, "Visualization frequency", 500, 2
        )
        self.stat_istart, _ = self.create_labeled_input(
            tab, "Start statistics from iter", 1000, 3
        )

        ttk.Label(tab, text="Write outlet plane data?").grid(
            row=4, column=0, sticky=tk.W, padx=5, pady=5
        )
        self.is_write_check = ttk.Checkbutton(
            tab, variable=self.is_write, command=self.on_write_changed
        )
        self.is_write_check.grid(row=4, column=1, sticky=tk.W, padx=5, pady=5)

        self.wrt_read_nfre1, self.wrt_read_nfre1_w = self.create_labeled_input(
            tab, "Plane data frequency", 1000, 5
        )
        self.wrt_read_nfre2, self.wrt_read_nfre2_w = self.create_labeled_input(
            tab, "Start saving from iter", 2001, 6
        )
        self.wrt_read_nfre3, self.wrt_read_nfre3_w = self.create_labeled_input(
            tab, "Stop saving at iter", 10000, 7
        )

        ttk.Label(tab, text="Read inlet plane data?").grid(
            row=8, column=0, sticky=tk.W, padx=5, pady=5
        )
        self.is_read = tk.BooleanVar(value=False)
        self.is_read_check = ttk.Checkbutton(tab, variable=self.is_read)
        self.is_read_check.grid(row=8, column=1, sticky=tk.W, padx=5, pady=5)

        self.on_write_changed()

    def on_write_changed(self):
        """Handle write checkbox changes."""
        enabled = self.is_write.get() or self.is_read.get()
        self.set_enabled(self.wrt_read_nfre1_w, enabled)
        self.set_enabled(self.wrt_read_nfre2_w, enabled)
        self.set_enabled(self.wrt_read_nfre3_w, enabled)

    def create_probe_tab(self):
        """Probe settings tab."""
        tab = ttk.Frame(self.notebook)
        self.notebook.add(tab, text="Probe")

        ttk.Label(tab, text="Auto-generate 5 probe points?").grid(
            row=0, column=0, sticky=tk.W, padx=5, pady=5
        )
        self.is_auto_probe = tk.BooleanVar(value=True)
        ttk.Checkbutton(tab, variable=self.is_auto_probe).grid(
            row=0, column=1, sticky=tk.W, padx=5, pady=5
        )
        
        ttk.Label(tab, text="(Probes will be auto-generated at runtime)", 
                 font=("Arial", 9, "italic")).grid(
            row=1, column=0, columnspan=2, sticky=tk.W, padx=5, pady=5
        )

    def generate_ini(self):
        """Generate the INI file from GUI inputs."""
        try:
            config = CustomConfigParser()
            case_num = self.icase.get()

            # Process
            config["process"] = {
                "is_prerun": bool_to_string(self.is_prerun.get()),
                "is_postprocess": bool_to_string(self.is_postprocess.get()),
            }

            # Decomposition
            is_decomp = self.is_decomp.get()
            config["decomposition"] = {
                "nxdomain": 1,
                "p_row": 0 if is_decomp else int(self.p_row.get()),
                "p_col": 0 if is_decomp else int(self.p_col.get()),
            }

            # Domain
            config["domain"] = {
                "icase": case_num,
                "lxx": float(self.lxx.get()),
                "lyt": float(self.lyt.get()),
                "lyb": float(self.lyb.get()),
                "lzz": float(self.lzz.get()),
            }

            # Flow - determine initfl based on restart and case
            if self.is_restart.get():
                initfl = Init.RESTART.value
                irestartfrom = int(self.irestartfrom.get())
            else:
                if case_num in [Case.CHANNEL.value, Case.DUCT.value, Case.PIPE.value, Case.ANNULAR.value]:
                    initfl = Init.POISEUILLE.value
                elif case_num == Case.TGV3D.value:
                    initfl = Init.FUNCTION.value
                else:
                    initfl = Init.POISEUILLE.value
                irestartfrom = 0

            config["flow"] = {
                "initfl": initfl,
                "irestartfrom": irestartfrom,
                "veloinit": f"{float(self.velo1.get())},{float(self.velo2.get())},{float(self.velo3.get())}",
                "noiselevel": float(self.noiselevel.get()),
                "reni": int(self.reni.get()),
                "nreni": int(self.nreni.get()),
                "ren": int(self.ren.get()),
            }

            # Thermo
            if self.ithermo.get():
                config["thermo"] = {
                    "ithermo": bool_to_string(self.ithermo.get()),
                    "icht": bool_to_string(self.icht.get()),
                    "igravity": int(self.igravity.get().split(":")[0]),
                    "ifluid": int(self.ifluid.get().split(":")[0]),
                    "ref_l0": float(self.refl0.get()),
                    "ref_T0": float(self.refT0.get()),
                    "inittm": int(self.inittm.get().split(":")[0]),
                    "irestartfrom": 0,
                    "Tini": float(self.Tini.get()),
                }

            # MHD
            if self.imhd.get():
                is_stuart = self.mhd_type.get().startswith("1")
                config["mhd"] = {
                    "imhd": bool_to_string(self.imhd.get()),
                    "NStuart": f"{bool_to_string(is_stuart)},{float(self.NS.get())}",
                    "NHartmn": f"{bool_to_string(not is_stuart)},{float(self.NH.get())}",
                    "B_static": f"{float(self.b1.get())},{float(self.b2.get())},{float(self.b3.get())}",
                }

            # Mesh
            istret_val = int(self.istret.get().split(":")[0])
            config["mesh"] = {
                "ncx": int(self.ncx.get()),
                "ncy": int(self.ncy.get()),
                "ncz": int(self.ncz.get()),
                "istret": istret_val,
                "rstret": f"{int(self.rstret1.get())},{float(self.rstret2.get())}"
                if istret_val != 0
                else "0,0.0",
            }

            # BC
            bc_dict = {}
            for key in self.bc_widgets:
                if "_1" in key or "_2" in key or "_v1" in key or "_v2" in key:
                    base_key = key.rsplit("_", 1)[0]
                    if base_key not in bc_dict:
                        # Build BC string for this boundary
                        bc1 = self.bc_widgets[f"{base_key}_1"].get().split(":")[0]
                        bc2 = self.bc_widgets[f"{base_key}_2"].get().split(":")[0]
                        v1 = self.bc_widgets[f"{base_key}_v1"].get()
                        v2 = self.bc_widgets[f"{base_key}_v2"].get()
                        bc_dict[base_key] = f"{bc1},{bc2},{v1},{v2}"

            config["bc"] = {
                "ifbcx_u": bc_dict.get("ifbcx_u", "1,1,0.0,0.0"),
                "ifbcx_v": bc_dict.get("ifbcx_u", "1,1,0.0,0.0"),
                "ifbcx_w": bc_dict.get("ifbcx_u", "1,1,0.0,0.0"),
                "ifbcx_p": bc_dict.get("ifbcx_p", "1,1,0.0,0.0"),
                "ifbcx_T": bc_dict.get("ifbcx_T", "1,1,0.0,0.0"),
                "ifbcy_u": bc_dict.get("ifbcy_u", "4,4,0.0,0.0"),
                "ifbcy_v": bc_dict.get("ifbcy_u", "4,4,0.0,0.0"),
                "ifbcy_w": bc_dict.get("ifbcy_u", "4,4,0.0,0.0"),
                "ifbcy_p": bc_dict.get("ifbcy_p", "5,5,0.0,0.0"),
                "ifbcy_T": bc_dict.get("ifbcy_T", "1,1,0.0,0.0"),
                "ifbcz_u": bc_dict.get("ifbcz_u", "1,1,0.0,0.0"),
                "ifbcz_v": bc_dict.get("ifbcz_u", "1,1,0.0,0.0"),
                "ifbcz_w": bc_dict.get("ifbcz_u", "1,1,0.0,0.0"),
                "ifbcz_p": bc_dict.get("ifbcz_p", "1,1,0.0,0.0"),
                "ifbcz_T": bc_dict.get("ifbcz_T", "1,1,0.0,0.0"),
                "idriven": int(self.idriven.get().split(":")[0]),
                "drivenfc": float(self.drivenCf.get()),
            }

            # Scheme
            config["scheme"] = {
                "dt": float(self.dt.get()),
                "iTimeScheme": 3,
                "iAccuracy": int(self.iAccuracy.get().split(":")[0]),
                "iviscous": 1,
            }

            # Simulation Control
            config["simcontrol"] = {
                "nIterFlowFirst": int(self.nIterFlowFirst.get()),
                "nIterFlowLast": int(self.nIterFlowLast.get()),
                "nIterThermoFirst": int(self.nIterThermoFirst.get()) if self.ithermo.get() else 0,
                "nIterThermoLast": int(self.nIterThermoLast.get()) if self.ithermo.get() else 0,
            }

            # I/O
            is_read = 1 if self.is_read.get() else 0
            is_write = 1 if self.is_write.get() else 0

            if is_write == 0 and is_read == 0:
                wrt_read_nfre = "0,0,0"
            else:
                wrt_read_nfre = f"{int(self.wrt_read_nfre1.get())},{int(self.wrt_read_nfre2.get())},{int(self.wrt_read_nfre3.get())}"

            config["io"] = {
                "cpu_nfre": int(self.cpu_nfre.get()),
                "ckpt_nfre": int(self.ckpt_nfre.get()),
                "visu_idim": 0,
                "visu_nfre": int(self.visu_nfre.get()),
                "visu_nskip": "1,1,1",
                "stat_istart": int(self.stat_istart.get()),
                "stat_nskip": "1,1,1",
                "is_wrt_read_bc": f"{bool_to_string(is_write)},{bool_to_string(is_read)}",
                "wrt_read_nfre": wrt_read_nfre,
            }

            # Probe - auto-generate 5 points
            probe_settings = self._generate_probe_settings()
            config["probe"] = probe_settings

            # Save file
            with open(DEFAULT_FILENAME, "w") as configfile:
                config.write(configfile)

            messagebox.showinfo("Success", f"Configuration saved to {DEFAULT_FILENAME}")

        except Exception as e:
            messagebox.showerror("Error", f"Failed to generate INI file:\n{str(e)}")

    def _generate_probe_settings(self):
        """Generate automatic 5 probe points."""
        npp = 5
        lxx = float(self.lxx.get())
        lzz = float(self.lzz.get())
        lyt = float(self.lyt.get())
        lyb = float(self.lyb.get())

        lxp = [lxx / 2.0] * npp
        lzp = [lzz / 2.0] * npp
        lyp = [lyb + (lyt - lyb) * (i + 1) / (npp + 1) for i in range(npp)]

        result = {"npp": npp}
        for i in range(npp):
            result[f"pt{i + 1}"] = f"{lxp[i]},{lyp[i]},{lzp[i]}"

        return result


if __name__ == "__main__":
    root = tk.Tk()
    gui = CHAPSimGUI(root)
    root.mainloop()
