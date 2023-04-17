'''Front-end code for HARMONI simulator
This handles the GUI and command line interfaces.
'''

import os
import sys
import argparse
import configparser
import multiprocessing as mp
import collections
import subprocess

from src.config import config_data


def get_version_number():
	try:
		result = subprocess.run(['git', 'describe'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
		version = result.stdout.strip().decode("utf-8")
		if version[0] == "v":
			version = version[1:]
		return version
	except:
		pass
	
	try:
		dir_path = os.path.dirname(os.path.realpath(__file__))
		with open(dir_path + '/PKG-INFO') as f:
			return f.readlines()[2][:-1].split('Version: ')[1]

	except:
		return "???"

def get_cpu_count():
	nprocs = mp.cpu_count() - 1
	if nprocs <= 0:
		nprocs = 1
	return nprocs

def get_grating_list():
	
	# Build grating list from the config file
	gr_low = []
	gr_mid = []
	gr_high = []
	for gr_name, gr_data in config_data["gratings"].items():
		if gr_data.R < 5000:
			gr_low.append((gr_name, gr_data.lmin))
		elif gr_data.R < 10000:
			gr_mid.append((gr_name, gr_data.lmin))
		else:
			gr_high.append((gr_name, gr_data.lmin))
	
	# sort by wavelenght
	gr_low = sorted(gr_low, key=lambda x: x[1])
	gr_mid = sorted(gr_mid, key=lambda x: x[1])
	gr_high = sorted(gr_high, key=lambda x: x[1])

	# join low, mid and high resolution gratings
	return [_[0] for _ in gr_low] + [_[0] for _ in gr_mid] + [_[0] for _ in gr_high]


if __name__ == "__main__":

	hsim_version = get_version_number()
	
	Parameter = collections.namedtuple("Parameter", "name,help,type,default,choices")
	Parameter.__new__.__defaults__ = (None, None, str, None, None)
	
	simulation_parameters = [Parameter("input_cube", "FITS input cube"),
				Parameter("output_dir", "Output directory"),
				Parameter("grating", "HARMONI grating", choices = get_grating_list()),
				Parameter("spaxel_scale", "Spaxel Scale", choices = ["4x4", "10x10", "20x20", "30x60", "60x60", "120x60"]),
				Parameter("exposure_time", "Exposure time [s]", type=int),
				Parameter("n_exposures", "Number of exposures", type=int),
				Parameter("ao_mode", "AO Mode", choices = ["LTAO", "SCAO", "HCAO", "noAO", "Airy", "User"]),
				Parameter("user_defined_psf", "User defined PSF FITS file when ao-mode is set to User", default="''"),
				Parameter("ao_star_hmag", "H magnitude of the LTAO AO star", type=float, default=17.5, choices = [15., 17.5, 19.]),
				Parameter("ao_star_distance", "Distance from the HARMONI FoV center to the LTAO AO star [arcsec]", type=int, default=30, choices = [15, 30, 45]),
				Parameter("hc_apodizer", "HCAO High-contrast apodizer", default="HSP1", choices = ["HSP1", "HSP2"]),
				Parameter("hc_fp_mask", "HCAO High-contrast focal plane mask", default="FPM1", choices = ["FPM1"]),
				Parameter("zenith_seeing", "Optical 500nm atmospheric seeing FWHM at zenith [arcsec]", type=float, choices = config_data["PSD_params"]["seeings"]),
				Parameter("air_mass", "Air mass of the observation", type=float, choices = config_data["PSD_params"]["air_masses"]),
				Parameter("moon_illumination", "Moon fractional illumination", type=float, default = 0., choices = [0.0, 0.5, 1.0]),
				Parameter("detector_systematics", "FITS input cube", default="False", choices = ["True", "False"]),
				Parameter("detector_tmp_path", "Directory to save interim detector files", default="''"),
				Parameter("adr", "Simulate atmospheric differential refraction", default="True", choices = ["True", "False"]),
				Parameter("mci", "Use minimum compliant instrument parameters", default="False", choices = ["True", "False"]),
				Parameter("telescope_temp", "Telescope temperature [K]", type=float, default = 280),
				Parameter("fprs_temp", "FPRS temperature [C]", type=float, default = +2),
				Parameter("scattered_sky", "Scattered sky fraction [%%]", type=float, default = 20),
				Parameter("extra_jitter", "Additional telescope PSF blur [mas]", type=str, default = "0"),
				Parameter("noise_seed", "Noise random number generator seed", type=int, default = 100),
				Parameter("n_cpus", "Number of processors", type=int, default = get_cpu_count()),
				Parameter("spectral_sampling", "Internal spectral oversampling factor", type=float, default = -1),
				Parameter("spatial_sampling", "Internal spatial oversampling factor", type=float, default = -1),
			  ]
	
	# Define argument parser
	parser = argparse.ArgumentParser(argument_default=argparse.SUPPRESS, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	
	parser.add_argument("-v", "--version", action="version", version="HARMONI Simulator version " + hsim_version)
	parser.add_argument("-b", dest="batch_mode", action="store_true", help="Batch mode. Do not show HSIM GUI")
	parser.add_argument("-c", dest="config_file", type=str, help="Simulation configuration file")
	
	parameter_actions = {}
	for param in simulation_parameters:
		parameter_actions[param.name] = parser.add_argument("--" + param.name.replace("_", "-"), dest=param.name, type=param.type, help=param.help if param.default is None else param.help + " (default: " + str(param.default) + ")", choices=param.choices)
	
	parameter_actions["debug"] = parser.add_argument("-d", "--debug", dest="debug", action="store_true", help="Produce debug outputs (default: False)")
		
	#
	args = parser.parse_args()
	
	# Get input parameters from command line and configuration file
	input_parameters = {}
	
	# Define default parameters
	for param in simulation_parameters:
		input_parameters[param.name.lower()] = param.default
	input_parameters["debug"] = False
	
	# Read ini configuration file if provided
	if hasattr(args, "config_file"):
		
		if not os.path.isfile(args.config_file):
			print("ERROR: Cannot open configuration file " + args.config_file)
			sys.exit()
		
		config = configparser.ConfigParser()
		config.read(args.config_file)
		
		if "HSIM" not in config:
			print("ERROR: [HSIM] section not found in configuration file " + args.config_file)
			sys.exit()
		
		for key in config["HSIM"]:
			value = config["HSIM"][key]
			key = key.lower()
			if key not in input_parameters:
				print("ERROR: Unknown option '" +  key + "' in configuration file " + args.config_file)
				sys.exit()
			
			action = parameter_actions[key]
			
			if action.type is not None:
				try:
					value = action.type(value)
				except:
					print("ERROR: '" +  str(value) + "' is not a valid " + str(action.type.__name__) + " value for '" +  key + "' in configuration file " + args.config_file)
					sys.exit()
			
			if action.choices is not None:
				if value not in action.choices:
					print("ERROR: Not valid value '" + str(value) + "' for '" +  key + "' in configuration file " + args.config_file)
					print("Valid values are: {" + ", ".join(map(str, action.choices)) + "}")
					sys.exit()
					
			input_parameters[key] = value
			
		input_parameters["config_file"] = args.config_file
	else:
		input_parameters["config_file"] = "command line"
		
		
	# Override with command line options
	for key, value in args.__dict__.items():
		key = key.lower()
		if key in input_parameters:
			input_parameters[key] = value
	
	
	input_parameters["version"] = get_version_number()
	
	if hasattr(args, "batch_mode") and args.batch_mode == True:
		# Batch mode
		
		# Check that all parameters are properly defined
		for key, value in input_parameters.items():
			if value is None:
				print("ERROR: '" + key + "' input parameter is not defined")
				sys.exit()
				
	
		from src.main import main
		main(input_parameters)
		
	else:
		# GUI mode
		from tkinter import *
		from tkinter import filedialog, messagebox
		from tkinter.ttk import *

		class panel_gui():
			def __init__(self, parent, title, column):
				self.parent = parent
				self.panel = Frame(parent)
				self.panel.grid(row=1, column=column, stick=N)
				Label(self.panel, text=title, font=(font_size_title)).grid(row=1, column=1, stick=W, padx=10)
				self.last_row = 1
				
				self.st_om = Style()
				self.st_om.configure("tm.TMenubutton", font=default_font)
				self.st_bt = Style()
				self.st_bt.configure("tm.TButton", font=default_font)

				
			def add_field(self, label, WidgetClass, **kw):
				self.last_row += 1
				Label(self.panel, text=label).grid(row=self.last_row, column=1, stick=W, padx=15)
				if WidgetClass is not None:
					var = StringVar()
					
					try:
						extra = kw["extra"].copy()
					except:
						extra = []
					
					extra_kw = {}
					if WidgetClass is Entry:
						extra_kw = {"textvariable":var, "width":10}
					elif WidgetClass is OptionMenu:
						extra_kw = {"style":"tm.TMenubutton"}
						var.set(str(extra[0]))
						extra.insert(0, "")
						extra.insert(0, var)
					elif WidgetClass is Button:
						extra_kw = {"textvariable":var, "command":kw["command"], "style":"tm.TButton"}
					else:
						extra_kw = {"variable":var}

					if "default" in kw:
						var.set(str(kw["default"]))
					
					if var.get() == "":
						var.set("0")
					
					w = WidgetClass(self.panel, *extra, **extra_kw)
					w.grid(row=self.last_row, column=2, stick=EW, padx=8, pady=5)
					
					return var
				
				return None

		class HSIM_GUI():
			def __init__(self, parent):
				
				# Menu 
				menubar = Menu(parent)
				parent.config(menu=menubar)
				
				file_menu = Menu(menubar)
				def OnAbout():
					messagebox.showinfo("HSIM", r"HSIM " + get_version_number() + "\nHARMONI simulation pipeline\nhttps://github.com/HARMONI-ELT/HSIM")

				
				file_menu.add_command(label="About", command=OnAbout)
				file_menu.add_command(label="Exit", command=parent.destroy)
				
				menubar.add_cascade(label="File", menu=file_menu)


				def create_field(name, panel_field):
					var = panel_field
					try:
						default_value = str(input_parameters[name]).strip()
						if name == "grating":
							for _ in grating_choices:
								if _.split("[")[0].strip() == default_value:
									var.set(_)
									break
						elif default_value != "None":
							if default_value == "True":
								default_value = "1"
							elif default_value == "False":
								default_value = "0"
							elif default_value == "''":
								default_value = "(None)"
							var.set(default_value)
					except:
						pass
					
					setattr(self, name, var)
				

				# Instrument frame
				panel_instrument = panel_gui(parent, "Instrument", 1)

				def browse_input_file(self):
					try:
						filename = filedialog.askopenfilename(filetypes = (("FITS files","*.fits"),("all files","*.*")))
						self.input_cube.set(os.path.relpath(filename))
					except:
						pass
				def browse_dir(self):
					try:
						filename = filedialog.askdirectory(initialdir = self.output_dir.get())
						self.output_dir.set(os.path.relpath(filename))
					except:
						pass
					
					
				create_field("input_cube", panel_instrument.add_field("Input cube", Button, command=lambda : browse_input_file(self), default="(None)"))
				create_field("output_dir", panel_instrument.add_field("Output dir", Button, command=lambda : browse_dir(self), default="./output_cubes"))
				create_field("exposure_time", panel_instrument.add_field("Exposure time [s]", Entry, default=600))
				create_field("n_exposures", panel_instrument.add_field("Number of exposures", Entry, default=3))
				spaxel_choices = list(parameter_actions["spaxel_scale"].choices)
				create_field("spaxel_scale", panel_instrument.add_field("Spaxel scale [mas]", OptionMenu, extra=spaxel_choices, default=spaxel_choices[1]))
	
				grating_choices = ["{name} [{info.lmin:.2f}-{info.lmax:.2f} um] (R={info.R:.0f})".format(name=_, info=config_data["gratings"][_]) for _ in get_grating_list()]
				create_field("grating", panel_instrument.add_field("Grating", OptionMenu, extra=grating_choices, default=grating_choices[6]))

				# Telescope frame
				def browse_psf_file(self):
					try:
						filename = filedialog.askopenfilename(filetypes = (("FITS files","*.fits"),("all files","*.*")))
						self.user_defined_psf.set(os.path.relpath(filename))
					except:
						pass

				panel_telescope = panel_gui(parent, "Telescope", 2)
				seeing_choices = list(map(str, parameter_actions["zenith_seeing"].choices))
				create_field("zenith_seeing", panel_telescope.add_field("Zenith seeing [arcsec]", OptionMenu, default=seeing_choices[2], extra=seeing_choices))
				ao_choices = list(parameter_actions["ao_mode"].choices)
				create_field("ao_mode", panel_telescope.add_field("AO mode", OptionMenu, extra=ao_choices))
				ao_star_hmag_choices = list(parameter_actions["ao_star_hmag"].choices)
				create_field("ao_star_hmag", panel_telescope.add_field("LTAO star H mag", OptionMenu, extra=ao_star_hmag_choices))
				ao_star_distance_choices = list(parameter_actions["ao_star_distance"].choices)
				create_field("ao_star_distance", panel_telescope.add_field("LTAO star distance [arcsec]", OptionMenu, extra=ao_star_distance_choices))
				hc_apodizer_choices = list(parameter_actions["hc_apodizer"].choices)
				create_field("hc_apodizer", panel_telescope.add_field("HCAO High-contrast apodizer", OptionMenu, extra=hc_apodizer_choices))
				hc_mask_choices = list(parameter_actions["hc_fp_mask"].choices)
				create_field("hc_fp_mask", panel_telescope.add_field("HCAO High-contrast mask", OptionMenu, extra=hc_mask_choices))
				create_field("user_defined_psf", panel_telescope.add_field("Used defined PSF", Button, command=lambda : browse_psf_file(self)))
				
				air_mass_choices = list(parameter_actions["air_mass"].choices)
				create_field("air_mass", panel_telescope.add_field("Air mass", OptionMenu, default=air_mass_choices[1], extra=air_mass_choices))
				create_field("moon_illumination", panel_telescope.add_field("Moon illumination", OptionMenu, extra=list(map(str, [0, 0.5, 1]))))
				

				# Misc frame
				panel_misc = panel_gui(parent, "Miscellaneous", 3)
				create_field("telescope_temp", panel_misc.add_field("Telescope temperature [K]", Entry))
				create_field("fprs_temp", panel_misc.add_field("FPRS temperature [C]", Entry))
				create_field("scattered_sky", panel_misc.add_field("Scattered sky fraction [%]", Entry))
				create_field("extra_jitter", panel_misc.add_field("Additional jitter [mas]", Entry))
				create_field("adr", panel_misc.add_field("ADR on/off", Checkbutton, default=1, height=1000))
				create_field("detector_systematics", panel_misc.add_field("Detector systematics", Checkbutton))
				create_field("detector_tmp_path", panel_misc.add_field("Detector tmp dir", Button, command=lambda : browse_dir(self)))
				create_field("n_cpus", panel_misc.add_field("No. of processors (1-" + str(mp.cpu_count())+")", Entry))
				create_field("noise_seed", panel_misc.add_field("Noise seed", Entry))
				panel_misc.add_field("Internal oversampling:", None)
				create_field("spectral_sampling", panel_misc.add_field("   Spectral (default = -1)", Entry))
				create_field("spatial_sampling", panel_misc.add_field("   Spatial (default = -1)", Entry))
				create_field("mci", panel_misc.add_field("Minimum compliant instrument", Checkbutton, default=1, height=1000))


				def OnClick():
					input_parameters["config_file"] = "GUI"
					for param in simulation_parameters:
						if param.choices == ["True", "False"]:
							input_parameters[param.name] = "True" if int(getattr(self, param.name).get()) == 1 else "False"
						elif param.name == "grating":
							input_parameters["grating"] = str(self.grating.get()).split("[")[0].strip()
						else:
							input_parameters[param.name] = param.type(getattr(self, param.name).get())
						
					from src.main import main
					main(input_parameters)

				Button(parent, text="Commence simulation", command=OnClick, style="tm.TButton").grid(row=2, column=2, pady=10)


		root = Tk()
		root.title("HARMONI Simulator Interface v" + get_version_number())
		font_size_title = "helvetica 20 bold"
		default_font = "helvetica 13"
		root.option_add("*Font", default_font)
		gui = HSIM_GUI(root)

		root.mainloop()
