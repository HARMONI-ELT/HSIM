'''Front-end code for HARMONI simulator
This handles the GUI and command line interfaces.
'''
import matplotlib
matplotlib.use('Agg')

import getopt
import sys
import os

import multiprocessing as mp

import os.path

from src.main import main
from src.config import config_data


if __name__=="__main__":

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
			print ""
			print '-'*19
			print ' HARMONI Simulator'
			print '-'*19
			print 'Version: ', ver
			print "---"
			print 'TO RUN GUI'
			print '>>> python hsim2.py'
			print ""
			print 'Command line'
			print '>>> python hsim2.py -c arg arg2 ...'
			print ""
			print "---"
			print 'OPTIONS'
			print '-h or --help = display this message and exit'
			print '-c or --cline = use command line. Use: >>> python hsim2.py -c to display arguments list'
			print '-p or --proc = set the number of processors when using command line (1-'+str(mp.cpu_count())+')'
			print '-o or --odir = set the output file directory when using the command line (default: ./output_cubes)'
			print ""
			sys.exit()

	nprocs = mp.cpu_count() - 1
	for o, a in optlist:
		if o in ("-p", "--proc"):
			nprocs = int(a)
		
		if nprocs <= 0:
			nprocs = 1
			print "Using 1 CPU"
		
		if nprocs > mp.cpu_count():
			print 'Only ' + str(mp.cpu_count()) + ' CPUs. Using ' + str(mp.cpu_count())
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
		if o in ("-c", "--cline") and len(args) != 13:
			print ""
			print 'COMMAND LINE USAGE'
			print ""
			print 'Enter command line arguments in following order:'
			print '1. datacube: Input datacube filepath'
			print '2. DIT: Detector Integration Time [s]'
			print '3. NDIT: No. of exposures'
			print '4. grating - V+R, Iz+J, H+K, Iz, J, H, K, z-high, J-high, H-high, K-short, K-long, J-short, J-long'
			print '5. spax: spatial pixel (spaxel) scale [mas] - 4x4, 10x10, 20x20, 30x60 '
			print '6. seeing: Atmospheric seeing FWHM [arcsec] - 0.43", 0.57", 0.64", 0.72", 1.04"'
			print '7. air mass - 1.1, 1.3, 1.5, 2.0'
			print '8. Moon fractional illumination - 0 0.5 1.0'
			print '9. jitter: Additional telescope PSF blur [mas]'
			print '10. site temp: Site/telescope temperature [K]'
			print '11. ADR on/off: (True/False) - atmospheric differential refraction'
			print '12. noise seed'
			print '13. AO mode [LTAO/SCAO//noAO/Airy]'
			#print '14. Use detector systematics (True/False)'
			print ""
			
			sys.exit()
		elif o in ("-c", "--cline") and len(args) == 13:
			if not os.path.exists(odir) or not os.path.isdir(odir):
				print "Output directory '" + odir + "'  does not exist or is not a directory. Exiting."
				sys.exit()

			datacube = str(args[0])
			DIT = int(args[1])
			NDIT = int(args[2])
			grat = str(args[3])
			spax = str(args[4])
			seeing = float(args[5])
			air_mass = float(args[6])
			moon = float(args[7])
			jitter = float(args[8])
			site_temp = float(args[9])
			adr = str(args[10])
			noise_seed = int(args[11])
			ao_mode = str(args[12])
			#systematics = str(args[12])

			#Start main function
			main(os.path.join(".", datacube), os.path.join(".", odir), DIT, NDIT, grat, spax, seeing, air_mass, ver,
				res_jitter=jitter, moon=moon, site_temp=site_temp, adr_switch=adr, det_switch='False',
				seednum=noise_seed, nprocs=nprocs, debug=debug, aoMode=ao_mode)

			sys.exit()

	#Use GUI interface if no command line option
	wxfound = 0
	if len(optlist) == 0 and len(args) == 0:
		try:
			import wx
		except:
			print ""
			print 'HARMONI Simulator'
			print '    --------     '
			print 'NO WX MODULE'
			print 'wxPython can be downloaded from: http://www.wxpython.org/'
			print ""
			print 'COMMAND LINE USAGE ONLY'
			print ""
			print '>>> python hsim.py -c arg1 arg2 ...'
			print ""
			print 'To see help:'
			print ""
			print '>>> python hsim.py -c -h'
			print ""
			print 'To see required arguments:'
			print ""
			print '>>> python hsim.py -c'
			print ""
			print "---"
			sys.exit()


		class Form1(wx.Frame):
			def __init__(self, parent, title):
				super(Form1, self).__init__(parent, title=title)

				self.InitUI()
				self.Centre()
				self.Show()

			def InitUI(self):
				panel = wx.Panel(self)
				titlefont = wx.Font(18, wx.DEFAULT, wx.NORMAL, wx.BOLD)

				# Set up the menu.
				filemenu  = wx.Menu()
				menuAbout = filemenu.Append(wx.ID_ABOUT, "&About"," Information about this program")
				menuExit  = filemenu.Append(wx.ID_EXIT,"&Exit"," Terminate the program")

				# Creating the menubar.
				menuBar = wx.MenuBar()
				menuBar.Append(filemenu,"&File") # Adding the "filemenu" to the MenuBar
				self.SetMenuBar(menuBar)  # Adding the MenuBar to the Frame content.

				# Events.
				self.Bind(wx.EVT_MENU, self.OnExit, menuExit)
				self.Bind(wx.EVT_MENU, self.OnAbout, menuAbout)
				
				vbox = wx.BoxSizer(wx.VERTICAL)
				hbox = wx.BoxSizer(wx.HORIZONTAL)

				grid_size = (7, 2, 15, 10)
				fg = wx.FlexGridSizer(*grid_size)

				#Instrument parameters
				inst = wx.StaticText(panel, label='Instrument')
				inst.SetFont(titlefont)
				subinst = wx.StaticText(panel, label=' ')
				INPUTCUBE = wx.StaticText(panel, label="Input cube")
				self.INPUTCUBEVAL = wx.FilePickerCtrl(panel, path="")
				DIT = wx.StaticText(panel, label="Exposure time [s]")
				self.DITVAL = wx.TextCtrl(panel, value='600')
				NDIT = wx.StaticText(panel, label="Number of exposures")
				self.NDITVAL = wx.TextCtrl(panel, value='3')
				SPAX = wx.StaticText(panel, label="Spaxel scale [mas]")
				self.SPAXVAL = wx.Choice(panel, choices=['4x4', '10x10', '20x20', '30x60'])
				self.SPAXVAL.SetStringSelection('10x10')
				PHOTOBAND = wx.StaticText(panel, label='Grating')
				grating_list = ["V+R", "Iz+J", "H+K", "Iz", "J", "H", "K", "z-high", "J-high", "J-short", "J-long", "H-high", "K-short", "K-long"]
				grating_choices = ["{name} [{info.lmin:.2f}-{info.lmax:.2f} um] (R={info.R:.0f})".format(name=_, info=config_data['gratings'][_]) for _ in grating_list]
				self.PHOTOBANDVAL = wx.Choice(panel, choices=grating_choices)
				self.PHOTOBANDVAL.SetStringSelection(grating_choices[6])
				DIR = wx.StaticText(panel, label='Output dir')
				self.DIRVAL = wx.DirPickerCtrl(panel, path=os.path.join(".", "output_cubes/"))
				
				fg.AddMany([(inst), (subinst), (INPUTCUBE),
						(self.INPUTCUBEVAL, 1, wx.EXPAND),
						(DIR), (self.DIRVAL, 1, wx.EXPAND),
						(DIT), (self.DITVAL, 1, wx.EXPAND),
						(NDIT), (self.NDITVAL, 1, wx.EXPAND),
						(SPAX), (self.SPAXVAL, 1, wx.EXPAND),
						(PHOTOBAND), (self.PHOTOBANDVAL, 1, wx.EXPAND),])
				
				fg.AddGrowableCol(1,1)

				hbox.Add(fg, flag=wx.ALL|wx.EXPAND, border=10)
				
				#Telescope parameters
				fgs = wx.FlexGridSizer(*grid_size)

				tele = wx.StaticText(panel, label='Telescope')
				tele.SetFont(titlefont)
				subtele = wx.StaticText(panel, label=' ')
				SEEING = wx.StaticText(panel, label="Zenith seeing [arcsec]:")
				self.SEEINGVAL = wx.Choice(panel, choices=map(str, sorted(config_data["PSD_cube"]["seeings"])))
				self.SEEINGVAL.SetStringSelection(str(config_data["PSD_cube"]["seeings"][2]))
				AOMODE = wx.StaticText(panel, label="AO mode")
				self.AOMODEVAL = wx.Choice(panel, choices=['LTAO', 'SCAO', 'noAO', 'Airy'])
				self.AOMODEVAL.SetStringSelection('LTAO')
				MOON = wx.StaticText(panel, label="Moon illumination")
				self.MOONVAL = wx.Choice(panel, choices=map(str, [0, 0.5, 1]))
				self.MOONVAL.SetStringSelection("0")
				AIRMASS = wx.StaticText(panel, label="Air mass")
				self.AIRMASSVAL = wx.Choice(panel, choices=map(str, sorted(config_data["PSD_cube"]["air_masses"])))
				self.AIRMASSVAL.SetStringSelection(str(config_data["PSD_cube"]["air_masses"][1]))
				SITETEMP = wx.StaticText(panel, label="Telescope temperature [K]:")
				self.SITETEMPVAL = wx.TextCtrl(panel, value='280')

				fgs.AddMany([(tele), (subtele),
						(SEEING), (self.SEEINGVAL, 1, wx.EXPAND),
						(AOMODE), (self.AOMODEVAL, 1, wx.EXPAND),
						(AIRMASS), (self.AIRMASSVAL, 1, wx.EXPAND),
						(MOON), (self.MOONVAL, 1, wx.EXPAND),
						(SITETEMP), (self.SITETEMPVAL, 1, wx.EXPAND)])
				
				fgs.AddGrowableCol(1,1)

				hbox.Add(fgs, flag=wx.ALL|wx.EXPAND, border=10)
				
				#Misc. parameters
				fgss = wx.FlexGridSizer(*grid_size)

				misc = wx.StaticText(panel, label="Miscellaneous")
				misc.SetFont(titlefont)
				submisc = wx.StaticText(panel, label=' ')
				NOISESEED = wx.StaticText(panel, label="Noise seed")
				self.NOISESEEDVAL = wx.TextCtrl(panel, value='100')
				N_PROC = wx.StaticText(panel, label='No. of processors (1-'+str(mp.cpu_count())+')')
				self.N_PROCVAL = wx.TextCtrl(panel, value=str(mp.cpu_count()-1))
				RESJIT = wx.StaticText(panel, label="Additional jitter [mas]:")
				self.RESJITVAL = wx.TextCtrl(panel, value='3')
				ADR = wx.StaticText(panel, label="ADR on/off")
				self.ADRVAL = wx.CheckBox(panel)
				self.ADRVAL.SetValue(True)
				#DET = wx.StaticText(panel, label="Detector systematics")
				#self.DETVAL = wx.CheckBox(panel)
				#self.DETVAL.SetValue(False)

				fgss.AddMany([(misc), (submisc),
						(RESJIT), (self.RESJITVAL, 1, wx.EXPAND),
						(ADR), (self.ADRVAL, 1, wx.EXPAND),
						#(DET), (self.DETVAL, 1, wx.EXPAND),
						(N_PROC), (self.N_PROCVAL, 1, wx.EXPAND),
						(NOISESEED), (self.NOISESEEDVAL, 1, wx.EXPAND)])

				fgss.AddGrowableCol(1,1)

				hbox.Add(fgss, flag=wx.ALL|wx.EXPAND, border=10)
				
				vbox.Add(hbox, flag=wx.ALL|wx.EXPAND, border=10)

				button = wx.Button(panel, 10, "Commence simulation")
				wx.EVT_BUTTON(panel, 10, self.OnClick)
				vbox.Add(button, flag=wx.ALL|wx.CENTER, border=10)

				panel.SetSizer(vbox)
				
				vbox.SetSizeHints(self)


			def OnClick(self,event):
				#Extract parameters for simulation run:
				cubefile = str(self.INPUTCUBEVAL.GetPath())
				ditval = float(self.DITVAL.GetValue())
				nditval = float(self.NDITVAL.GetValue())
				spaxval = str(self.SPAXVAL.GetStringSelection())
				photoband = str(self.PHOTOBANDVAL.GetStringSelection()).split(' ')[0]
				seeingval = float(self.SEEINGVAL.GetStringSelection())
				aomode = str(self.AOMODEVAL.GetStringSelection()).split(' ')[0]
				airmassaval = float(self.AIRMASSVAL.GetStringSelection())
				moon = float(self.MOONVAL.GetStringSelection())
				resjitval = float(self.RESJITVAL.GetValue())
				sitetempval = float(self.SITETEMPVAL.GetValue())
				noiseseedval = int(self.NOISESEEDVAL.GetValue())
				nprocs = int(self.N_PROCVAL.GetValue())
				odir = str(self.DIRVAL.GetPath())
				return_adrval = str(self.ADRVAL.GetValue())
				#return_detval = str(self.DETVAL.GetValue())
				
				#start main program
				main(os.path.join(".", cubefile), os.path.join(".", odir), ditval, nditval, photoband, spaxval, seeingval, airmassaval, ver,
					res_jitter=resjitval, site_temp=sitetempval, adr_switch=return_adrval,
					seednum=noiseseedval, nprocs=nprocs, aoMode=aomode, moon=moon)


			def OnExit(self,e):
				self.Close(True)  # Close the frame.

			def OnAbout(self,e):
				info = wx.AboutDialogInfo()
				info.SetName('HSIM')
				info.SetVersion(ver)
				info.SetDescription('HARMONI simulation pipeline')
				info.SetWebSite('http://www-astro.physics.ox.ac.uk/instr/HARMONI/Simulator/simulator.html')
				wx.AboutBox(info)


		app = wx.App()
		Form1(None, title="HARMONI Simulator Interface")
		app.MainLoop()
		
