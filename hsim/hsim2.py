'''Front-end code for HARMONI simulator
This handles the GUI and command line interfaces.
'''
import getopt
import sys
import os
import multiprocessing as mp
import os.path
import matplotlib
matplotlib.use('Agg')

from src.main import main
from src.config import config_data


if __name__ == "__main__":

	#Get version number
	try:
		dir_path = os.path.dirname(os.path.realpath(__file__))
		with open(dir_path + '/PKG-INFO') as f:
			ver = f.readlines()[2][:-1].split('Version: ')[1]

	except:
		ver = "???"

	optlist, args = getopt.getopt(sys.argv[1:], 'dhcp:o:', ['debug', 'help', 'cline', 'proc', 'odir'])

	for o, a in optlist:
		if o in ("-h", "--help"):
			print("")
			print('-'*19)
			print(' HARMONI Simulator')
			print('-'*19)
			print('Version: ', ver)
			print("---")
			print('TO RUN GUI')
			print('>>> python hsim2.py')
			print("")
			print('Command line')
			print('>>> python hsim2.py -c arg arg2 ...')
			print("")
			print("---")
			print('OPTIONS')
			print('-h or --help = display this message and exit')
			print('-c or --cline = use command line. Use: >>> python hsim2.py -c to display arguments list')
			print('-p or --proc = set the number of processors when using command line (1-'+str(mp.cpu_count())+')')
			print('-o or --odir = set the output file directory when using the command line (default: ./output_cubes)')
			print("")
			sys.exit()

	nprocs = mp.cpu_count() - 1
	for o, a in optlist:
		if o in ("-p", "--proc"):
			nprocs = int(a)

		if nprocs <= 0:
			nprocs = 1
			print("Using 1 CPU")

		if nprocs > mp.cpu_count():
			print('Only ' + str(mp.cpu_count()) + ' CPUs. Using ' + str(mp.cpu_count()))
			nprocs = mp.cpu_count()

	odir = 'output_cubes'
	for o, a in optlist:
		if o in ("-o", "--odit"):
			odir = a
			break

	debug = False
	for o, a in optlist:
		if o in ("-d"):
			debug = True
			break

	for o, a in optlist:
		if o in ("-c", "--cline") and len(args) != 14:
			print("")
			print('COMMAND LINE USAGE')
			print("")
			print('Enter command line arguments in following order:')
			print('1. datacube: Input datacube filepath')
			print('2. DIT: Detector Integration Time [s]')
			print('3. NDIT: No. of exposures')
			print('4. grating - V+R, Iz+J, H+K, Iz, J, H, K, z-high, J-high, H-high, K-short, K-long, J-short, J-long')
			print('5. spax: spatial pixel (spaxel) scale [mas] - 4x4, 10x10, 20x20, 30x60 ')
			print('6. seeing: Atmospheric seeing FWHM [arcsec] - 0.43, 0.57, 0.64, 0.72, 1.04')
			print('7. air mass - 1.1, 1.3, 1.5, 2.0')
			print('8. Moon fractional illumination - 0 0.5 1.0')
			print('9. jitter: Additional telescope PSF blur [mas]')
			print('10. site temp: Site/telescope temperature [K]')
			print('11. ADR on: (True/False) - atmospheric differential refraction')
			print('12. noise seed')
			print('13. AO mode [LTAO/SCAO//noAO/Airy/User defined PSF fits file]')
			print('14. Use IR detector systematics (True/False)')
			print("")

			sys.exit()
		elif o in ("-c", "--cline") and len(args) == 14:
			if not os.path.exists(odir) or not os.path.isdir(odir):
				print("Output directory '" + odir + "'  does not exist or is not a directory. Exiting.")
				sys.exit()
				
			datacube = str(args[0])
			DIT = int(args[1])
			NDIT = int(args[2])
			grat = str(args[3])
			spax = str(args[4])
			seeing = float(args[5])
			air_mass = float(args[6])
			moon = float(args[7])
			jitter = str(args[8])
			site_temp = float(args[9])
			adr = str(args[10])
			noise_seed = int(float(args[11]))
			ao_mode = str(args[12])
			systematics = str(args[13])

			#Start main function
			main(os.path.join(".", datacube), os.path.join(".", odir), DIT, NDIT, grat, spax, seeing, air_mass, ver,
				res_jitter=jitter, moon=moon, site_temp=site_temp, adr_switch=adr, det_switch=systematics,
				seednum=noise_seed, nprocs=nprocs, debug=debug, aoMode=ao_mode)

			sys.exit()

	#Use GUI interface if no command line option
	if len(optlist) == 0 and len(args) == 0:
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
						extra = kw["extra"]
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
					ver = "212"
					messagebox.showinfo("HSIM", r"HSIM " + ver + "\nHARMONI simulation pipeline\nhttps://github.com/HARMONI-ELT/HSIM")

				
				file_menu.add_command(label="About", command=OnAbout)
				file_menu.add_command(label="Exit", command=parent.destroy)
				
				menubar.add_cascade(label="File", menu=file_menu)


				# Instrument frame
				panel_instrument = panel_gui(parent, "Instrument", 1)

				def browse_file(self):
					try:
						filename = filedialog.askopenfilename(filetypes = (("FITS files","*.fits"),("all files","*.*")))
						self.input_cube.set(os.path.relpath(filename))
					except:
						pass

				
				def browse_dir(self):
					try:
						filename = filedialog.askdirectory(initialdir = self.outputdir.get())
						self.outputdir.set(os.path.relpath(filename))
					except:
						pass
					
					
				self.input_cube = panel_instrument.add_field("Input cube", Button, command=lambda : browse_file(self), default="(None)")
				self.outputdir = panel_instrument.add_field("Output dir", Button, command=lambda : browse_dir(self), default="./output_cubes")
				self.exp_time = panel_instrument.add_field("Exposure time [s]", Entry, default=600)
				self.n_exp = panel_instrument.add_field("Number of exposures", Entry, default=3)
				self.spax_scale = panel_instrument.add_field("Spaxel scale [mas]", OptionMenu, extra=["4x4", "10x10", "20x20", "30x60"])
				grating_list = ["V+R", "Iz+J", "H+K", "Iz", "J", "H", "K", "z-high", "J-high", "J-short", "J-long", "H-high", "K-short", "K-long"]
				grating_choices = ["{name} [{info.lmin:.2f}-{info.lmax:.2f} um] (R={info.R:.0f})".format(name=_, info=config_data["gratings"][_]) for _ in grating_list]
				self.grating = panel_instrument.add_field("Grating", OptionMenu, extra=grating_choices, default=grating_choices[6])


				# Telescope frame
				panel_telescope = panel_gui(parent, "Telescope", 2)
				seeing_choices = list(map(str, sorted(config_data["PSD_cube"]["seeings"])))
				self.seeing = panel_telescope.add_field("Zenith seeing [arcsec]", OptionMenu, default=seeing_choices[2], extra=seeing_choices)
				self.ao_mode = panel_telescope.add_field("AO mode", OptionMenu, extra=["LTAO", "SCAO", "noAO", "Airy"])
				air_mass_choices = list(map(str, sorted(config_data["PSD_cube"]["air_masses"])))
				self.air_mass = panel_telescope.add_field("Air mass", OptionMenu, default=air_mass_choices[1], extra=air_mass_choices)
				self.moon = panel_telescope.add_field("Moon illumination", OptionMenu, extra=list(map(str, [0, 0.5, 1])))
				self.tel_temp = panel_telescope.add_field("Telescope temperature [K]", Entry, default=280)


				# Misc frame
				panel_misc = panel_gui(parent, "Miscellaneous", 3)
				self.jitter = panel_misc.add_field("Additional jitter [mas]", Entry, default=3)
				self.adr_var = panel_misc.add_field("ADR on/off", Checkbutton, default=1, height=1000)
				self.det_var = panel_misc.add_field("Detector systematics", Checkbutton)
				self.ncpu = panel_misc.add_field("No. of processors (1-" + str(mp.cpu_count())+")", Entry, default=mp.cpu_count()-1)
				self.noise = panel_misc.add_field("Noise seed", Entry, default=100)


				def OnClick():
					cubefile = str(self.input_cube.get())
					ditval = float(self.exp_time.get())
					nditval = float(self.n_exp.get())
					spaxval = str(self.spax_scale.get())
					photoband = str(self.grating.get()).split(' ')[0]
					seeingval = float(self.seeing.get())
					aomode = str(self.ao_mode.get()).split(' ')[0]
					airmassaval = float(self.air_mass.get())
					moon = float(self.moon.get())
					resjitval = str(self.jitter.get())
					sitetempval = float(self.tel_temp.get())
					noiseseedval = int(self.noise.get())
					nprocs = int(self.ncpu.get())
					odir = str(self.outputdir.get())
					return_adrval = "True" if int(self.adr_var.get()) == 1 else "False"
					return_detval = "True" if int(self.det_var.get()) == 1 else "False"
					
					#start main program
					main(os.path.join(".", cubefile), os.path.join(".", odir), ditval, nditval, photoband, spaxval, seeingval, airmassaval, ver,
						res_jitter=resjitval, site_temp=sitetempval, adr_switch=return_adrval, det_switch=return_detval,
						seednum=noiseseedval, nprocs=nprocs, aoMode=aomode, moon=moon)

				Button(parent, text="Commence simulation", command=OnClick, style="tm.TButton").grid(row=2, column=2, pady=10)


		root = Tk()
		root.title("HARMONI Simulator Interface")
		font_size_title = "helvetica 20 bold"
		default_font = "helvetica 13"
		root.option_add("*Font", default_font)
		gui = HSIM_GUI(root)

		#lab.pack()
		root.mainloop()
