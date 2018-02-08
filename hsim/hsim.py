'''Front-end code for HARMONI simulator
This handles the GUI and command line interfaces.

Author: Simon Zieleniewski

Last updated: 09-09-15

'''


import getopt
import sys
import os
import numpy as n
import multiprocessing as mp
from src.TheSimulator import main
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
            print 'HARMONI Simulator'
            print '    --------     '
            print 'SVN Version: ', ver
            print"---"
            print 'TO RUN'
            print 'GUI'
            print '>>> python hsim.py'
            print""
            print 'Command line'
            print '>>> python hsim.py -c arg arg2 ...'
            print"---"
            print 'OPTIONS'
            print '-h or --help = display this message and exit'
            print '-c or --cline = use command line. Use: >>> python hsim.py -c to display arguments list'
            print '-p or --proc = set the number of processors when using command line (1-'+str(mp.cpu_count())+')'
            print '-o or --odir = set the output file directory when using the command line (default: /hsim-#/Output_cubes)'
            print""
            sys.exit()

    nprocs = mp.cpu_count()-1           
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
            
    for o, a in optlist:
        if o in ("-c", "--cline") and len(args) != 20:
            print""
            print 'COMMAND LINE USAGE'
            print ""
            print 'Enter command line arguments in following order:'
            print '1. datacube: Input datacube filepath'
            print '2. DIT: Exposure time [s]'
            print '3. NDIT: No. of exposures'
            print '4. grating: [V+R, Iz+J, H+K, Iz, J, H, K, z, J-high, H-high, K-high, None]'
            print '5. x-spax: x spatial pixel (spaxel) scale [mas]'
            print '6. y-spax: y spatial pixel (spaxel) scale [mas]'
            print '7. Telescope: Telescope type (E-ELT, VLT)'
            print '8. AO: AO mode (LTAO, SCAO, Gaussian)'
            print '9. Seeing: Atmospheric seeing FWHM [arcsec] - Between 0.67"-1.10"'
            print '10. Zenith angle: Zenith angle [deg]'
            print '11. User uploaded PSF: Path to user uploaded PSF FITS file - Enter None if not required'
            print '12. PSF blur: Additional telescope PSF blur [mas]'
            print '13. Site temp: Site/telescope temperature [K]'
            print '14. Spec_Nyquist: (True/False) - Set spectral sampling to Nyquist sample spectral resolution element'
            print '15. Spec samp: Spectral sampling [A/pix] - Only used if Spec_Nyquist=False, but enter a value regardless!'
            print '16. Noise force seed: Force random number seed to take a set value (0=No, 1-10=yes and takes that value)'
            print '17. Remove background: (True/False) - Subtract background spectrum'
            print '18. Return object cube: (True/False) - Return object FITS cube'
            print '19. Return transmission cube: (True/False) - Return transmission FITS cube'
            print '20. Turn ADR off: (True/False) - atmospheric differential refraction'
            print ""
            sys.exit()
        elif o in ("-c", "--cline") and len(args) == 20:
            print args
            datacube = str(args[0])
            DIT = int(args[1])
            NDIT = int(args[2])
            grat = str(args[3])
            spax = (float(args[4]), float(args[5]))
            telescope = str(args[6])
            AO = str(args[7])
            seeing = float(args[8])
            zenith_ang = float(args[9])
            user_PSF = str(args[10])
            res_jitter = float(args[11])
            site_temp = float(args[12])
            combine_ndits = True
            spec_Nyquist = args[13]
            if spec_Nyquist == 'True':
                spec_N = True
            elif spec_Nyquist == 'False':
                spec_N = False
            Spec_samp = float(args[14])
            noise_force_seed = int(float(args[15]))
            remove_background = args[16]
            return_obj = args[17]
            return_tra = args[18]
            return_adr = args[19]
            #Start main function
            main(datacube, odir, DIT, NDIT, grat, spax, seeing, zenith_ang, telescope,
                 user_PSF, AO, res_jitter, site_temp, combine_ndits, spec_N, Spec_samp,
                 noise_force_seed, remove_background, return_obj, return_tra, return_adr,
                 ver, nprocs)
            
##    else:
##        print""
##        print 'HARMONI Simulator'
##        print""
##        print "Unknown option keys. Use either -h for help, or -c for command line usage"
##        print""
##        sys.exit()
        

    #Use GUI interface if no command line option
    wxfound = 0
    if len(optlist) == 0 and len(args) == 0:
        try:
            import wx
            wxfound = 1
            class Form1(wx.Frame):
                
                def __init__(self, parent, title):
                    super(Form1, self).__init__(parent, title=title, 
                        size=(1120, 390))

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
                    XSPAX = wx.StaticText(panel, label="X Scale [mas]")
                    self.XSPAXVAL = wx.TextCtrl(panel, value='20')
                    YSPAX = wx.StaticText(panel, label="Y Scale [mas]")
                    self.YSPAXVAL = wx.TextCtrl(panel, value='20')
                    PHOTOBAND = wx.StaticText(panel, label='Grating')
                    self.PHOTOBANDVAL = wx.Choice(panel, choices=['V+R [0.47-0.81 um] (R=3500)', 'Iz+J [0.8-1.36 um]  (R=3500)', 'H+K [1.45-2.45 um] (R=3500)',
                                                                  'Iz [0.82-1.03 um] (R=7500)','J [1.08-1.36 um] (R=7500)', 'H [1.46-1.83 um] (R=7500)', 'K [1.95-2.45 um] (R=7500)',
                                                                  'z [0.82-0.91 um] (R=17000)','J-high [1.17-1.29 um] (R=17000)','H-high [1.545-1.715 um] (R=17000)','K-high [2.09-2.32 um] (R=17000)',
                                                                  'None'])
                    DIR = wx.StaticText(panel, label='Output Dir')
                    self.DIRVAL = wx.DirPickerCtrl(panel, path=path_setup('../../Output_cubes/'))
                    
                    fg.AddMany([(inst), (subinst), (INPUTCUBE),
                                (self.INPUTCUBEVAL, 1, wx.EXPAND),
                                (DIR), (self.DIRVAL, 1, wx.EXPAND),
                                (DIT), (self.DITVAL, 1, wx.EXPAND),
                                (NDIT), (self.NDITVAL, 1, wx.EXPAND),
                                (XSPAX), (self.XSPAXVAL, 1, wx.EXPAND),
                                (YSPAX), (self.YSPAXVAL, 1, wx.EXPAND),
                                (PHOTOBAND), (self.PHOTOBANDVAL, 1, wx.EXPAND),])
                    
                    fg.AddGrowableCol(1,1)

                    hbox.Add(fg, proportion=1, flag=wx.ALL|wx.EXPAND, border=10)
                       
                    #Telescope parameters(AO mode, seeing, zenith angle etc):
                    fgs = wx.FlexGridSizer(8, 2, 15, 10)

                    tele = wx.StaticText(panel, label='Telescope')
                    tele.SetFont(titlefont)
                    subtele = wx.StaticText(panel, label=' ')
                    TELESCOPE = wx.StaticText(panel, label="Telescope:")
                    self.TELESCOPEVAL = wx.Choice(panel, choices=['E-ELT', 'VLT'])
                    AOMODE = wx.StaticText(panel, label="AO Mode:")
                    self.AOMODEVAL = wx.Choice(panel, choices=['LTAO', 'SCAO', 'Gaussian'])
                    SEEING = wx.StaticText(panel, label="Zenith Seeing [arcsec]:")
                    self.SEEINGVAL = wx.TextCtrl(panel, value='0.67')
                    ZENITHA = wx.StaticText(panel, label="Zenith Angle [deg]")
                    self.ZENITHAVAL = wx.TextCtrl(panel, value='0')
                    USERPSF = wx.StaticText(panel, label='User PSF (replaces AO choice)')
                    self.USERPSFVAL = wx.FilePickerCtrl(panel, path="None")
                    gapa = wx.StaticText(panel, label=' ------------ ')
                    gapb = wx.StaticText(panel, label=' ')
                    SITETEMP = wx.StaticText(panel, label="Telescope Temperature [K]:")
                    self.SITETEMPVAL = wx.TextCtrl(panel, value='280.5')

                    fgs.AddMany([(tele), (subtele), (TELESCOPE),
                                 (self.TELESCOPEVAL, 1, wx.EXPAND),
                                 (AOMODE), (self.AOMODEVAL, 1, wx.EXPAND),
                                 (SEEING), (self.SEEINGVAL, 1, wx.EXPAND),
                                 (ZENITHA), (self.ZENITHAVAL, 1, wx.EXPAND),
                                 (gapa), (gapb),
                                 (USERPSF), (self.USERPSFVAL, 1, wx.EXPAND),
                                 (SITETEMP), (self.SITETEMPVAL, 1, wx.EXPAND)])
                    
                    fgs.AddGrowableCol(1,1)

                    hbox.Add(fgs, proportion=1, flag=wx.ALL|wx.EXPAND, border=10)
                    
                    #Misc. parameters
                    fgss = wx.FlexGridSizer(10, 2, 15, 10)

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
                    self.N_PROCVAL = wx.TextCtrl(panel, value=str(mp.cpu_count()-1))
                    RESJIT = wx.StaticText(panel, label="Additional PSF Blur [mas]:")
                    self.RESJITVAL = wx.TextCtrl(panel, value='0')
                    ADR = wx.StaticText(panel, label="Turn off ADR")
                    self.ADRVAL = wx.CheckBox(panel)

                    fgss.AddMany([(misc), (submisc),
                                  (REMOVE_BG), (self.REMOVE_BGVAL, 1, wx.EXPAND),
                                  (OUTPUT_OB), (self.OUTPUT_OBVAL, 1, wx.EXPAND),
                                  (OUTPUT_TS), (self.OUTPUT_TSVAL, 1, wx.EXPAND),
                                  (gapaa), (gapbb),
                                  (N_PROC), (self.N_PROCVAL, 1, wx.EXPAND),
                                  (NOISESEED), (self.NOISESEEDVAL, 1, wx.EXPAND),
                                  (SETR), (self.SETRVAL, 1, wx.EXPAND),
                                  (RESJIT), (self.RESJITVAL, 1, wx.EXPAND),
                                  (ADR), (self.ADRVAL, 1, wx.EXPAND)])

                    fgss.AddGrowableCol(1,1)

                    hbox.Add(fgss, proportion=1, flag=wx.ALL|wx.EXPAND, border=10)

                    panel.SetSizer(hbox)

                    # A button
                    button =wx.Button(panel, 10, "Commence Simulation", wx.Point(465, 325))
                    wx.EVT_BUTTON(panel, 10, self.OnClick)

                def OnClick(self,event):       
                    #Extract parameters for simulation run:
                    cubefile = str(self.INPUTCUBEVAL.GetPath())
                    ditval = float(self.DITVAL.GetValue())
                    nditval = float(self.NDITVAL.GetValue())
                    xspaxval = float(self.XSPAXVAL.GetValue())
                    yspaxval = float(self.YSPAXVAL.GetValue())
                    photoband = str(self.PHOTOBANDVAL.GetStringSelection()).split(' ')[0]
                    telescopeval = str(self.TELESCOPEVAL.GetStringSelection())
                    aomodeval = str(self.AOMODEVAL.GetStringSelection())
                    seeingval = float(self.SEEINGVAL.GetValue())
                    zenithaval = float(self.ZENITHAVAL.GetValue())
                    resjitval = float(self.RESJITVAL.GetValue())
                    sitetempval = float(self.SITETEMPVAL.GetValue())
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
                    return_adrval = str(self.ADRVAL.GetValue())
                    nprocs = int(self.N_PROCVAL.GetValue())
                    odir = str(self.DIRVAL.GetPath())
                    user_PSF = str(self.USERPSFVAL.GetPath())
                    combnditsval=True

                    print 'Filename: ', cubefile
                    print 'DIT = ', ditval
                    print 'NINT = ', nditval
                    print 'X spaxels = ', xspaxval
                    print 'Y spaxels = ', yspaxval
                    print 'Grating = ', photoband
                    print 'Telescope: ', telescopeval
                    print 'AO = ', aomodeval
                    print 'Seeing = ', seeingval
                    print 'Zenith angle = ', zenithaval
                    print 'Residual jitter = ', resjitval
                    print 'Temperature = ', sitetempval
                    print 'Noise seed = ', noiseseedval
                    print 'Spectral Nyquist sampling? ', spec_nyq
                    print 'Subtract background? ', remove_bgval
                    print 'User PSF? ' , user_PSF
                    print 'Return object cube? ', return_obval
                    print 'Return transmission? ', return_tsval
                    print 'ADR off? ', return_adrval
                    print 'No. of processors = ', nprocs
                    print 'Output directory = ', odir
                           
                    #start main program
                    main(cubefile, odir, ditval, nditval, photoband, (xspaxval, yspaxval),
                         seeingval, zenithaval, telescopeval, user_PSF, aomodeval, resjitval,
                         sitetempval, combnditsval, spec_nyq, setrval, noiseseedval, remove_bgval,
                         return_obval, return_tsval, return_adrval, ver, nprocs)

                def OnExit(self,e):
                    self.Close(True)  # Close the frame.

                def OnAbout(self,e):
                    info = wx.AboutDialogInfo()
                    #info.SetIcon(wx.Icon('hunter.png', wx.BITMAP_TYPE_PNG))
                    info.SetName('HSIM')
                    info.SetVersion(ver)
                    info.SetDescription('HARMONI simulation pipeline')
                    #info.SetCopyright('(C) 2014')
                    info.SetWebSite('http://www-astro.physics.ox.ac.uk/instr/HARMONI/Simulator/simulator.html')
                    #info.SetLicence(licence)
                    info.AddDeveloper('Simon Zieleniewski, Sarah Kendrew, Ryan Houghton & Niranjan Thatte')
                    info.AddDocWriter('Simon Zieleniewski')
                    #info.AddArtist('Artist')
                    #info.AddTranslator('Translator')
                    wx.AboutBox(info)

        except:
            print""
            print 'HARMONI Simulator'
            print '    --------     '
            print 'NO WX MODULE'
            print 'wxPython can be downloaded from: http://www.wxpython.org/'
            print""
            print 'COMMAND LINE USAGE ONLY'
            print""
            print '>>> python hsim.py -c arg1 arg2 ...'
            print""
            print 'To see help:'
            print""
            print '>>> python hsim.py -c -h'
            print""
            print 'To see required arguments:'
            print""
            print '>>> python hsim.py -c'
            print""
            print"---"
            sys.exit()
        if wxfound:
            #GUI interface
            app = wx.App()
            Form1(None, title="HARMONI Simulator Interface")
            app.MainLoop()



