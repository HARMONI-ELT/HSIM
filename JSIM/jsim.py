'''Front-end code for JWST/NIRSpec simulator
This handles the GUI and command line interfaces.

Author: Simon Zieleniewski

Created: 16-10-15

Last updated: 03-12-15

'''


import getopt
import sys
import os
import numpy as n
import multiprocessing as mp
from src.JSIM_main import main
import commands
from src.modules.misc_utils import path_setup



if __name__=="__main__":

    #Get version number
    try:
        verpath = path_setup('../../')
        verinfo = open(verpath+'PKG-INFO')
        ver = verinfo.readlines()[2][:-1].split('Version: ')[1]
        print 'VERSION = ', ver
        verinfo.close()
    except:
        ver = str(commands.getoutput("svnversion -c ./ | sed -e 's/[MS]//g' -e 's/^[[:digit:]]*://'"))


    optlist, args = getopt.getopt(sys.argv[1:], 'hcp:o:', ['help', 'cline', 'proc', 'odir'])

    for o, a in optlist:
        if o in ("-h", "--help"):
            print""
            print '    --------     '
            print 'JWST/NIRSpec Simulator'
            print '    --------     '
            print 'SVN Version: ', ver
            print"---"
            print 'TO RUN'
            print 'GUI'
            print '>>> python jsim.py'
            print""
            print 'Command line'
            print '>>> python jsim.py -c arg arg2 ...'
            print"---"
            print 'OPTIONS'
            print '-h or --help = display this message and exit'
            print '-c or --cline = use command line. Use: >>> python jsim.py -c to display arguments list'
            print '-p or --proc = set the number of processors when using command line (1-'+str(mp.cpu_count())+')'
            print '-o or --odir = set the output file directory when using the command line (default: /JSIM-#/Output_cubes)'
            print""
            sys.exit()

    # nprocs = mp.cpu_count()-1
    nprocs = 1
    for o, a in optlist:
        if o in ("-p", "--proc"):
            nprocs = int(a)
            if nprocs > mp.cpu_count():
                print 'Only '+str(mp.cpu_count())+' CPUs. Using '+str(mp.cpu_count())
                nprocs = mp.cpu_count()

    odir = path_setup('../../Output_cubes/')
    for o, a in optlist:
        if o in ("-o", "--odit"):
            if os.path.exists(a) and os.path.isdir(a):
                odir = a
            elif not os.path.exists(a):
                print "OUTPUT DIRECTORY DOESN'T EXIST (OR IS NOT A DIRECTORY)! PLEASE CHOOSE AGAIN. EXITING"
                sys.exit()
    print 'LEN ARGS = ', len(args)
    for o, a in optlist:
        if o in ("-c", "--cline") and len(args) != 13:
            print""
            print 'COMMAND LINE USAGE'
            print ""
            print 'Enter command line arguments in following order:'
            print '1. datacube: Input datacube filepath'
            print '2. DIT: Exposure time [s]'
            print '3. NDIT: No. of exposures'
            print '4. grating: [F070-G140M, F100-G140M, G235M, G395M, F070-G140H, F100-G140H, G235H, G395H, Prism]'
            print '5. Ignore LSF convolution? (True/False) - turnes off LSF convolution in spectral dimension'
            print '6. PSF blur: Additional telescope PSF blur [mas]'
            print '7. Spec_Nyquist: (True/False) - Set spectral sampling to Nyquist sample spectral resolution element'
            print '8. Spec samp: Spectral sampling [A/pix] - Only used if Spec_Nyquist=False, but enter a value regardless!'
            print '9. Noise force seed: Force random number seed to take a set value (0=No, 1-10=yes and takes that value)'
            print '10. Remove background: (True/False) - Subtract background spectrum'
            print '11. Return object cube: (True/False) - Return object FITS cube'
            print '12. Return transmission cube: (True/False) - Return transmission FITS cube'
            print '13. Approximate Drizzle technique: (True/False) - Output sampling of 50 mas'
            print ""
            sys.exit()
        elif o in ("-c", "--cline") and len(args) == 13:
            print args
            datacube = str(args[0])
            DIT = int(args[1])
            NDIT = int(args[2])
            grat = str(args[3])
            ignoreLSFval = str(args[4])
            if ignoreLSFval == 'True':
                ignoreLSF = True
            elif ignoreLSFval == 'False':
                ignoreLSF = False
            res_jitter = float(args[5])
            # site_temp = float(args[6])
            site_temp = 1.0
            combine_ndits = True
            spec_Nyquist = args[6]
            if spec_Nyquist == 'True':
                spec_N = True
            elif spec_Nyquist == 'False':
                spec_N = False
            Spec_samp = float(args[7])
            noise_force_seed = int(float(args[8]))
            remove_background = args[9]
            return_obj = args[10]
            return_tra = args[11]
            drizzle = args[12]
            if drizzle == 'True':
                drizzleval = True
            elif drizzle == 'False':
                drizzleval = False
            #Start main function
            main(datacube, odir, DIT, NDIT, grat, ignoreLSF, res_jitter,
                 site_temp, combine_ndits, spec_N, Spec_samp, noise_force_seed,
                 remove_background, return_obj, return_tra, drizzleval, ver, nprocs)

    #Use GUI interface if no command line option
    wxfound = 0
    if len(optlist) == 0 and len(args) == 0:
        try:
            import wx
            wxfound = 1
            class Form1(wx.Frame):

                def __init__(self, parent, title):
                    super(Form1, self).__init__(parent, title=title,
                        size=(850, 390))

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

                    hbox = wx.BoxSizer(wx.HORIZONTAL)

                    fg = wx.FlexGridSizer(8, 2, 15, 10)

                    #Instrument parameters
                    inst = wx.StaticText(panel, label='Instrument')
                    inst.SetFont(titlefont)
                    subinst = wx.StaticText(panel, label=' ')
                    INPUTCUBE = wx.StaticText(panel, label="Input Cube")
                    self.INPUTCUBEVAL = wx.FilePickerCtrl(panel, path="")
                    DIT = wx.StaticText(panel, label="DIT [s]")
                    self.DITVAL = wx.TextCtrl(panel, value='900')
                    NDIT = wx.StaticText(panel, label="NINT")
                    self.NDITVAL = wx.TextCtrl(panel, value='1')
                    PHOTOBAND = wx.StaticText(panel, label='Grating')
                    self.PHOTOBANDVAL = wx.Choice(panel, choices=['F070-G140M [0.7-1.2 um] (R=1000)','F100-G140M [1.0-1.8 um] (R=1000)',
                                                                  'G235M [1.7-3.0 um]  (R=1000)', 'G395M [2.9-5.0 um] (R=1000)',
                                                                  'F070-G140H [0.7-1.2 um] (R=2700)', 'F100-G140H [1.0-1.8 um] (R=2700)',
                                                                  'G235H [1.7-3.0 um] (R=2700)', 'G395H [2.9-5.0 um] (R=2700)',
                                                                  'Prism [0.6-5.0 um] (R~100)'])
                    IGNORESPEC = wx.StaticText(panel, label='Ignore LSF convolution')
                    self.IGNORESPECVAL = wx.CheckBox(panel)
                    DRIZZLE = wx.StaticText(panel, label='Drizzle (50 mas sampling)')
                    self.DRIZZLEVAL = wx.CheckBox(panel)
                    DIR = wx.StaticText(panel, label='Output Dir')
                    self.DIRVAL = wx.DirPickerCtrl(panel, path=path_setup('../../Output_cubes/'))

                    fg.AddMany([(inst), (subinst), (INPUTCUBE),
                                (self.INPUTCUBEVAL, 1, wx.EXPAND),
                                (DIR), (self.DIRVAL, 1, wx.EXPAND),
                                (DIT), (self.DITVAL, 1, wx.EXPAND),
                                (NDIT), (self.NDITVAL, 1, wx.EXPAND),
                                (PHOTOBAND), (self.PHOTOBANDVAL, 1, wx.EXPAND),
                                (IGNORESPEC), (self.IGNORESPECVAL, 1, wx.EXPAND),
                                (DRIZZLE), (self.DRIZZLEVAL, 1, wx.EXPAND)])

                    fg.AddGrowableCol(1,1)

                    hbox.Add(fg, proportion=1, flag=wx.ALL|wx.EXPAND, border=10)

                    # #Telescope parameters
                    # fgs = wx.FlexGridSizer(2, 2, 15, 10)
                    # tele = wx.StaticText(panel, label='Telescope')
                    # tele.SetFont(titlefont)
                    # subtele = wx.StaticText(panel, label=' ')
                    # # SITETEMP = wx.StaticText(panel, label="Temperature [K]:")
                    # # self.SITETEMPVAL = wx.TextCtrl(panel, value='3.0')
                    # fgs.AddMany([(tele), (subtele),
                    #              (USERPSF), (self.USERPSFVAL, 1, wx.EXPAND)])
                    #             #  (SITETEMP), (self.SITETEMPVAL, 1, wx.EXPAND)])
                    # fgs.AddGrowableCol(1,1)
                    # hbox.Add(fgs, proportion=1, flag=wx.ALL|wx.EXPAND, border=10)

                    #Misc. parameters
                    fgss = wx.FlexGridSizer(9, 2, 15, 10)

                    misc = wx.StaticText(panel, label="Miscellaneous")
                    misc.SetFont(titlefont)
                    submisc = wx.StaticText(panel, label=' ')
                    REMOVE_BG = wx.StaticText(panel, label='Subtract Background')
                    self.REMOVE_BGVAL = wx.CheckBox(panel)
                    OUTPUT_OB = wx.StaticText(panel, label='Return Object Cube')
                    self.OUTPUT_OBVAL = wx.CheckBox(panel)
                    OUTPUT_TS = wx.StaticText(panel, label='Return Transmission Cube')
                    self.OUTPUT_TSVAL = wx.CheckBox(panel)
                    gapaa = wx.StaticText(panel, label=' ------------ ')
                    gapbb = wx.StaticText(panel, label=' ')
                    NOISESEED = wx.StaticText(panel, label="Noise Seed")
                    self.NOISESEEDVAL = wx.Choice(panel, choices=['Random', '1', '2'])
                    SETR = wx.StaticText(panel, label="Set Spec Samp [A/pix]")
                    self.SETRVAL = wx.TextCtrl(panel)
                    N_PROC = wx.StaticText(panel, label='No. of processors (1-'+str(mp.cpu_count())+')')
                    self.N_PROCVAL = wx.TextCtrl(panel, value=str(1))
                    RESJIT = wx.StaticText(panel, label="Additional PSF Blur [mas]:")
                    self.RESJITVAL = wx.TextCtrl(panel, value='0')

                    fgss.AddMany([(misc), (submisc),
                                  (REMOVE_BG), (self.REMOVE_BGVAL, 1, wx.EXPAND),
                                  (OUTPUT_OB), (self.OUTPUT_OBVAL, 1, wx.EXPAND),
                                  (OUTPUT_TS), (self.OUTPUT_TSVAL, 1, wx.EXPAND),
                                  (gapaa), (gapbb),
                                  (N_PROC), (self.N_PROCVAL, 1, wx.EXPAND),
                                  (NOISESEED), (self.NOISESEEDVAL, 1, wx.EXPAND),
                                  (SETR), (self.SETRVAL, 1, wx.EXPAND),
                                  (RESJIT), (self.RESJITVAL, 1, wx.EXPAND)])

                    fgss.AddGrowableCol(1,1)

                    hbox.Add(fgss, proportion=1, flag=wx.ALL|wx.EXPAND, border=10)

                    panel.SetSizer(hbox)

                    # A button
                    button =wx.Button(panel, 10, "Commence Simulation", wx.Point(360, 325))
                    wx.EVT_BUTTON(panel, 10, self.OnClick)

                def OnClick(self,event):
                    #Extract parameters for simulation run:
                    cubefile = str(self.INPUTCUBEVAL.GetPath())
                    ditval = float(self.DITVAL.GetValue())
                    nditval = float(self.NDITVAL.GetValue())
                    photoband = str(self.PHOTOBANDVAL.GetStringSelection()).split(' ')[0]
                    ignorespecval = str(self.IGNORESPECVAL.GetValue())
                    if ignorespecval == 'True':
                        ignoreLSF = True
                    elif ignorespecval == 'False':
                        ignoreLSF = False
                    resjitval = float(self.RESJITVAL.GetValue())
                    # sitetempval = float(self.SITETEMPVAL.GetValue())
                    sitetempval = 1.0
                    noiseseedval = self.NOISESEEDVAL.GetStringSelection()
                    if noiseseedval == 'Random':
                        noiseseedval = 0.
                    setrchoice = self.SETRVAL.GetValue()
                    if setrchoice:
                        spec_nyq = False
                        setrval = float(setrchoice)
                    elif not setrchoice:
                        spec_nyq = True
                        setrval = 1.
                    remove_bgval = str(self.REMOVE_BGVAL.GetValue())
                    return_obval = str(self.OUTPUT_OBVAL.GetValue())
                    return_tsval = str(self.OUTPUT_TSVAL.GetValue())
                    drizzle = str(self.DRIZZLEVAL.GetValue())
                    if drizzle == 'True':
                        drizzleval = True
                    elif drizzle == 'False':
                        drizzleval = False
                    nprocs = int(self.N_PROCVAL.GetValue())
                    odir = str(self.DIRVAL.GetPath())
                    combnditsval=True

                    print 'Filename: ', cubefile
                    print 'DIT = ', ditval
                    print 'NINT = ', nditval
                    print 'Grating = ', photoband
                    print 'Ignore LSF conv?', ignorespecval
                    print 'Residual jitter = ', resjitval
                    print 'Temperature = ', sitetempval
                    print 'Noise seed = ', noiseseedval
                    print 'Spectral Nyquist sampling? ', spec_nyq
                    print 'Subtract background? ', remove_bgval
                    print 'Return object cube? ', return_obval
                    print 'Return transmission? ', return_tsval
                    print 'Drizzle?', drizzleval
                    print 'No. of processors = ', nprocs
                    print 'Output directory = ', odir

                    #start main program
                    main(cubefile, odir, ditval, nditval, photoband,
                    ignoreLSF, resjitval, sitetempval,
                    combnditsval, spec_nyq, setrval, noiseseedval,
                    remove_bgval, return_obval, return_tsval, drizzleval,
                    ver, nprocs)

                def OnExit(self,e):
                    self.Close(True)  # Close the frame.

                def OnAbout(self,e):
                    info = wx.AboutDialogInfo()
                    #info.SetIcon(wx.Icon('hunter.png', wx.BITMAP_TYPE_PNG))
                    info.SetName('JSIM')
                    info.SetVersion(ver)
                    info.SetDescription('JWST/NIRSpec simulation pipeline')
                    #info.SetCopyright('(C) 2014')
                    #info.SetWebSite('http://www.url.co.uk')
                    #info.SetLicence(licence)
                    info.AddDeveloper('Simon Zieleniewski, Sarah Kendrew, Ryan Houghton & Niranjan Thatte')
                    info.AddDocWriter('Simon Zieleniewski')
                    #info.AddArtist('Artist')
                    #info.AddTranslator('Translator')
                    wx.AboutBox(info)

        except:
            print""
            print 'JWST/NIRSpec Simulator'
            print '    --------     '
            print 'NO WX MODULE'
            print 'wxPython can be downloaded from: http://www.wxpython.org/'
            print""
            print 'COMMAND LINE USAGE ONLY'
            print""
            print '>>> python jsim.py -c arg1 arg2 ...'
            print""
            print 'To see help:'
            print""
            print '>>> python jsim.py -c -h'
            print""
            print 'To see required arguments:'
            print""
            print '>>> python jsim.py -c'
            print""
            print"---"
            sys.exit()
        if wxfound:
            #GUI interface
            app = wx.App()
            Form1(None, title="JSIM/NIRSpec Interface")
            app.MainLoop()
