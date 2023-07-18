# -*- coding: utf-8 -*-
"""
Created on Fri May 22 08:14:19 2020

@author: ERIC
"""



import yaml
import numpy as np
import nmrglue as ng
from matplotlib import pyplot as plt
from matplotlib.ticker import NullFormatter


class NMRexpt:
    
    vkys_yaml_str = """

seqfil:
    type: s
    description: Pulse sequence name
    units: 
        
ACQtime: 
    type: f
    description: Acquision time
    units: s
    
H1reffrq: 
    type: f
    description: Proton Larmor Frequncy for spectrometer
    units: MHz
    
acqdim: 
    type: i
    description: Number of dimension in experiment, could be 1, 2 or 3 ...
    units: None
    
at: 
    type: f
    description: Acquisition time in direct dimension, f2 in NMR tradition
    units: seconds
    
at1: 
    type: f     
    description: Acquisition time in indirect dimension f1 in NMR tradition    
    units: seconds  
dfrq: 
    type: f
    description: Transmitter frequency of first decoupler
    units: MHz
dn:
    type: s
    description: Nucleus for first decoupler
    units: 
dreffrq: 
    type: f
    description: Reference frequency for frst decoupler
    units: MHz
explabel: 
    type: s
    description: Name of pulse sequence
    units: 
fn: 
    type: i
    description: Fourier number in directly detected dimension
    units:
fn1: 
    type: i
    description: Fourier number in 1st indirectly detected dimension
    units:
lb: 
    type: f
    description: Line broadening in directly detected dimension
    units: Hz
lb1: 
    type: f
    description: Line broadening in 1st indirectly detected dimension
    units: Hz
lp: 
    type: f
    description: First-order phase in directly detected dimension
    units: Degrees
lp1: 
    type: f
    description: First-order phase in 1st indirectly detected dimension
    units: Degrees
lp2: 
    type: f
    description: First-order phase in 2nd indirectly detected dimension
    units: Degrees
np: 
    type: i
    description: Number of data points
    units:
nt: 
    type: i
    description: Number of transients
    units:
obsSW: 
    type: s
    description: Holds command to set specttral width limits
    units:
reffrq: 
    type: f
    description: Reference frequency of reference line
    units: MHz
reffrq1: 
    type: f
    description: Reference freq. of reference line in 1st indirect dimension
    units: MHz
rfl: 
    type: f
    description: Reference peak position in directly detected dimension
    units: Hz
rfl1: 
    type: f
    description: Reference peak position in 1st indirectly detected dimension
    units: Hz
rfp: 
    type: f
    description: Reference peak frequency in directly detected dimension
    units: ppm
rfp1: 
    type: f
    description: Reference peak frequency in in-directly detected dimension
    units: ppm
rl: 
    type: f
    description: Set reference line in directly detected dimension
    units:
rp: 
    type: f
    description: Zero-order phase in directly detected dimension
    units: degrees
rp1: 
    type: f
    description: Zero-order phase in 1st indirectly detected dimension
    units: degrees
rp2: 
    type: f
    description: Zero-order phase in 2nd indirectly detected dimension
    units: degrees
sample: 
    type: s
    description: Experiment directory name
    units:
samplename: 
    type: s
    description: Name of group the experiment belongs to
    units:
sfrq: 
    type: f
    description: Transmitter frequency of observe nucleus
    units: MHz
sp: 
    type: f
    description: Start of plot in directly detected dimension
    units: ppm
sreffrq: 
    type: f
    description: 
    units: MHz
sw: 
    type: f
    description: Spectral width in directly detected dimension
    units: Hz
sw1: 
    type: f
    description: Spectral width in 1st indirectly detected dimension
    units: Hz
tn: 
    type: s
    description: Nucleus for observe transmitter
    units:
tnref: 
    type: s
    description: Reference molecule and solvent
    units:
tof: 
    type: f
    description: Frequency offset for observe transmitter
    units: MHz
wp1: 
    type: f
    description: Width of plot in directly detected dimension
    units: ppm

"""   
    
    vkys = yaml.safe_load(vkys_yaml_str)
    
    
    def __init__(self ):
        
        pass
        
    def produce_spectrum(self,  **kwargs):
        
        self.update_udic(**kwargs)
        
        if self.udic['ndim'] == 1:
            
            self.produce_1D(**kwargs)
            
        elif self.udic['ndim'] == 2:
            
            self.produce_2D(**kwargs)
       
        
    def produce_1D(self, **kwargs):
        
        self.update_udic(**kwargs)
        
        H1 = 0
        
        zf = self.udic[H1]['zf']
        rp = self.udic[H1]['rp']
        lp = self.udic[H1]['lp']
        lb = self.udic[H1]['lb']
        # sw = self.udic[H1]['sw']
        
        fid = ng.proc_base.zf_size(self.data, zf)
        fid = ng.proc_base.ps(fid, p0=rp, p1=lp)
        # fid = ng.process.proc_base.em(fid,lb/sw)
        fid = ng.process.proc_base.em(fid,lb)
        sss1 = ng.process.proc_base.fft(fid)
                
        self.sss_0 = sss1/sss1.max()
        
        del fid
        del sss1

 
    def produce_2D(self, **kwargs):
        
        encoding = self.udic[0]['encoding']
        
        if encoding == 'states-tppi':
            
            # print('plot 2D states-tppi spectrum')
            
            self.zf_lb_fft_2D_states_tppi()
            
        elif encoding == 'states':
            
            # print('plot 2D states spectrum')
            pass
            
        elif encoding == 'magnitude':
            
            # print('plot 2D magnitude spectrum')
            self.zf_lb_fft_2D_magnitude()

        

    def plot_spectrum(self,  **kwargs):
        
        self.update_udic(**kwargs)
        
        if self.udic['ndim'] == 1:
            
            self.plot_1D(**kwargs)
            
        elif self.udic['ndim'] == 2:
            
            self.plot_2D(**kwargs)
            
            
        
    def plot_1D(self, fig=None, ax=None,  **kwargs):
        
        if isinstance(fig, type(None)) and isinstance(ax, type(None)):
            fig, ax = plt.subplots(figsize=(9, 5))
        else:
            # clear anything drawn to ax
            ax.clear()
            
        ax.set_ylim([-0.2, 1.2])
        fig.canvas.toolbar_visible = False
        fig.canvas.header_visible  = False
        fig.canvas.footer_visible  = False
        fig.canvas.resizable       = False
        
        ax.spines['top'].set_visible(False)
        ax.spines['left'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.set_xlabel('ppm', fontsize=14, gid='ppm')
        ax.set_yticks([])   
        
        h1_ppm_axis = self.udic[0]['unit_conversion'].ppm_scale()
        h1_ppm_limits = self.udic[0]['unit_conversion'].ppm_limits()
        
        full_spectrum = ax.plot(h1_ppm_axis, self.sss_0.real, color='black', lw=0.5)
        ax.set_xlim(h1_ppm_limits)
    
    
    def plot_2D(self, **kwargs):
        
        encoding = self.udic[0]['encoding']
        
        if encoding == 'states-tppi':
            
            # print('plot 2D states-tppi spectrum')
            
            self.zf_lb_fft_2D_states_tppi()
            self.matplotlib_plot2D()
            
        elif encoding == 'states':
            
            # print('plot 2D states spectrum')
            pass
            
        elif encoding == 'magnitude':
            
            # print('plot 2D magnitude spectrum')
            self.zf_lb_fft_2D_magnitude()
            self.matplotlib_plot2D(pos_and_neg = False)
            
     
    def zf_lb_fft_2D_magnitude(self):
        
        H1 = 1
        H1i = 0
        
        zf = self.udic[H1]['zf']
        # rp = self.udic[H1]['rp']
        # lp = self.udic[H1]['lp']
        lb = self.udic[H1]['lb']
        # sw = self.udic[H1]['sw']
        
        fid = ng.proc_base.zf_size(self.data, zf)
        # fid = ng.proc_base.ps(fid, p0=rp, p1=lp)
        # fid = ng.process.proc_base.em(fid,lb/sw)
        fid = ng.process.proc_base.sine(fid)
        # fid = ng.process.proc_base.em(fid,lb)
        fid1 = ng.process.proc_base.fft(fid)
        
        zf = self.udic[H1i]['zf']
        # rp = self.udic[H1i]['rp']
        # lp = self.udic[H1i]['lp']
        lb = self.udic[H1i]['lb']
        # sw = self.udic[H1i]['sw']
        
        fid1 = ng.proc_base.zf_size(fid1.T, zf)
        # fidri = ng.process.proc_base.em(fidri,lb/sw)
        fid1 = ng.process.proc_base.em(fid1,lb)
        # fid = ng.process.proc_base.sine(fid)
        sss1 = ng.process.proc_base.fft(fid1)
        
        sss1 = np.abs(sss1)
        
        self.sss_0 = sss1/sss1.max()
        
        del fid
        del fid1
        del sss1
        
        
        
    def zf_lb_fft_2D_states_tppi(self):
        
        H1=1
        C13=0
        
        zf = self.udic[H1]['zf']
        rp = self.udic[H1]['rp']
        lp = self.udic[H1]['lp']
        lb = self.udic[H1]['lb']
        sw = self.udic[H1]['sw']
    

        fid = ng.proc_base.zf_size(self.data, zf)
        # print("p0 =", rp)
        # print("p1 =", lp)
        fid = ng.proc_base.ps(fid, p0=rp, p1=lp)
        # fid = ng.process.proc_base.em(fid,lb/sw)
        fid = ng.process.proc_base.em(fid,lb)
        fid1 = ng.process.proc_base.fft(fid)
                
        fidr = fid1[::2]
        fidi = fid1[1::2]
        
        # fidrr = fidr.real + 1j*fidi.real
        fidri = fidr.real + 1j*fidi.imag
        fidir = fidr.imag + 1j*fidi.real
        # fidii = fidr.imag + 1j*fidi.imag
        
        # fidrr = fidrr.T
        fidri = fidri.T
        fidir = fidir.T
        # fidii = fidii.T
        
        zf = self.udic[C13]['zf']
        rp = self.udic[C13]['rp']
        lp = self.udic[C13]['lp']
        lb = self.udic[C13]['lb']
        sw = self.udic[C13]['sw']
                
        fidri = ng.proc_base.zf_size(fidri, zf)
        fidri = ng.proc_base.ps(fidri, p0=rp, p1=lp)
        # fidri = ng.process.proc_base.em(fidri,lb/sw)
        fidri = ng.process.proc_base.em(fidri,lb)
        sssri = ng.process.proc_base.fft(fidri)
        
        
        fidir = ng.proc_base.zf_size(fidir, zf)
        fidir = ng.proc_base.ps(fidir, p0=rp, p1=lp)
        # fidir = ng.process.proc_base.em(fidir,lb/sw)
        fidir = ng.process.proc_base.em(fidir, lb)
        sssir = ng.process.proc_base.fft(fidir)
        
        
        sss = sssri.real + sssir.imag
        
        if abs(sss.min()) < abs(sss.max()):
            smax = sss.max()
        else:
            smax = abs(sss.min())
            
        self.sss_0 = sss/smax
        
        del sss
        del fidir
        del fidri
        del fidr
        del fidi
        del fid1
        del fid
        

    def matplotlib_plot2D(self, pos_and_neg = True):  
        
        # print("sss_0",self.sss_0.shape, self.sss_0.max(), self.sss_0.min())
        

        dddCC = (self.sss_0.T).copy()
        dddA = (self.sss_0.T).copy()
        # dddB = dddA.copy()
        
        
        dddA[dddA<0.01]=0 
        
        if pos_and_neg:
            dddB = dddA.copy()
            dddB[dddB>-0.04]=0
            dddC = dddA + dddB
        else:
            dddC = dddA

        H1=1
        C13=0
        h1_ppm_axis = self.udic[H1]['unit_conversion'].ppm_scale()
        h1_ppm_limits = self.udic[H1]['unit_conversion'].ppm_limits()
        
        c13_ppm_axis = self.udic[C13]['unit_conversion'].ppm_scale()
        c13_ppm_limits = self.udic[C13]['unit_conversion'].ppm_limits()
        
        if pos_and_neg:
            ppp = ng.analysis.peakpick.pick(dddCC, pthres=0.1, nthres=-0.1, table=False )
        else:
            ppp = ng.analysis.peakpick.pick(dddCC, pthres=0.02, table=False )
        
        xy = (np.array(ppp[0])).T
        self.xy = xy
        
        c13_pk_pos = c13_ppm_axis[::-1][xy[0]]
        h1_pk_pos = h1_ppm_axis[xy[1]]
             
        nullfmt   = NullFormatter()         # no labels
        
        left, width = 0.1, 0.65
        bottom, height = 0.1, 0.65
        bottom_h = left_h = left+width+0.02
        
        rect_contour   = [left, bottom, width, height]
        rect_line_H1x  = [left, bottom_h, width, 0.2]
        rect_line_C13y = [left_h, bottom, 0.2, height]
        
        # start with a rectangular Figure
        plt.figure(2, figsize=(6,6))
        
        axContour = plt.axes(rect_contour)
        ax_line_xH1 = plt.axes(rect_line_H1x)
        ax_line_yC13 = plt.axes(rect_line_C13y)
        axContour.grid()
        
        axContour.contour(h1_ppm_axis, c13_ppm_axis[::-1], dddA, colors='blue')
        if pos_and_neg:
            axContour.contour(h1_ppm_axis, c13_ppm_axis[::-1], dddB, colors='red')
        
        axContour.set_xlim(h1_ppm_limits)
        axContour.set_ylim(c13_ppm_limits)
        
        ax_line_xH1.plot(h1_ppm_axis, dddC.max(axis=0), 'b', lw=0.66)
        if pos_and_neg:
            ax_line_xH1.plot(h1_ppm_axis, dddC.min(axis=0), 'r', lw=0.66)
        ax_line_xH1.set_xlim(h1_ppm_limits)
        
        ax_line_yC13.plot(dddC.max(axis=1), c13_ppm_axis[::-1], 'b', lw=0.66)
        if pos_and_neg:
            ax_line_yC13.plot(dddC.min(axis=1), c13_ppm_axis[::-1], 'r', lw=0.66)
        ax_line_yC13.set_ylim(c13_ppm_limits)
        
        # no labels
        ax_line_xH1.xaxis.set_major_formatter(nullfmt)
        ax_line_yC13.yaxis.set_major_formatter(nullfmt)
        
        axContour.tick_params(axis='x', labelsize=14)
        axContour.tick_params(axis='y', labelsize=14)
        axContour.set_xlabel('{} [ppm]'.format(self.udic[H1]['label']), fontsize=14)
        axContour.set_ylabel('{} [ppm]'.format(self.udic[C13]['label']), fontsize=14)
        
        ax_line_xH1.axis('off')
        ax_line_yC13.axis('off');
        
        axContour.scatter(h1_pk_pos,c13_pk_pos, alpha=0.5, s=20)
        
        plt.show() 


        
         
        
    def update_udic( self, **kwargs ):
        
        # print(kwargs)
        
        for ky, vals in kwargs.items():
            
            if ky in self.udic[0]:
                
                if isinstance(vals,list):
                    # print(ky, vals)
            
                    for n,v in enumerate(vals):
                        self.udic[n][ky]=v
                        
                        
        # update axis information 
        
        for n in range(self.udic['ndim']):
        
            self.udic[n]['unit_conversion'] = ng.fileiobase.unit_conversion(self.udic[n]['zf'],
                                                                     self.udic[n]['complex'], 
                                                                     self.udic[n]['sw'], 
                                                                     self.udic[n]['obs'], 
                                                                     self.udic[n]['car'])
        
        
        
    @classmethod    
    def from_varian( class_object, varian_expt_dirname='.', **kwargs):
        
        co =class_object()
        
        co.data_origin = "varian"
        
        co.dicv, co.data = ng.varian.read(varian_expt_dirname, **kwargs)
        
        co.fillin_varian_paramters()
        
        co.create_udic()
        
        return co
      
    
    def create_udic( self):

        dim_index = ["1", "",] + [str(i) for i in range(2,10)]
        
        self.udic = ng.varian.guess_udic(self.dicv, self.data )
        
        if self.udic['ndim'] == 1:
            dim_index = ["", "1",] + [str(i) for i in range(2,10)]
        else:
            dim_index = ["1", "",] + [str(i) for i in range(2,10)]
            
        
        for n in range(self.udic['ndim']):                
            
            dimension_char = dim_index[n]
            
            self.udic[n]['sw']  = self.varian_parameters['sw'+     dimension_char][0]
            self.udic[n]['obs'] = self.varian_parameters['reffrq'+ dimension_char][0]
            self.udic[n]['zf']  = self.varian_parameters['fn'+     dimension_char][0]
            
            self.udic[n]['rp'] = self.varian_parameters['rp'+     dimension_char][0]
            self.udic[n]['lp'] = self.varian_parameters['lp'+     dimension_char][0]
            
            self.udic[n]['lb'] = self.varian_parameters['lb'+     dimension_char][0]
            
            
            self.udic[n]['car'] = self.udic[n]['sw']/2 - self.varian_parameters['rfl'+dimension_char][0]
            
            self.udic[n]['unit_conversion'] = ng.fileiobase.unit_conversion(self.udic[n]['zf'],
                                                                             self.udic[n]['complex'], 
                                                                             self.udic[n]['sw'], 
                                                                             self.udic[n]['obs'], 
                                                                             self.udic[n]['car'])
        if self.udic['ndim'] > 1:  
            # Change label to nucleus rather than X,Y,Z     
            for n in range(self.udic['ndim']):                
                
                dimension_char = dim_index[n]
    
                if n==1: # direct dimension
                    
                    self.udic[n]['label'] = self.varian_parameters['tn'][0]
                    
                else:
                    if self.udic[1]['obs'] == self.udic[n]['obs']:
                        
                        self.udic[n]['label'] = self.varian_parameters['tn'][0]
                        
                    else:
                        # print(self.udic[n]['obs']//1, self.varian_parameters['dfrq'][0]//1)
                        if self.udic[n]['obs']//1 == self.varian_parameters['dfrq'][0]//1:
                        
                            self.udic[n]['label'] = self.varian_parameters['dn'][0]
                
            

    def fillin_varian_paramters(self ):
        
        self.varian_parameters = {}
        
        for k in self.dicv['procpar']:
            if k in NMRexpt.vkys.keys():

                if NMRexpt.vkys[k]['type'] == 'f':
                    values = [float(v) for v in self.dicv['procpar'][k]['values']]
                    self.varian_parameters[k]=values
                elif NMRexpt.vkys[k]['type'] == 'i':
                    values = [int(v) for v in self.dicv['procpar'][k]['values']]
                    self.varian_parameters[k]=values
                elif NMRexpt.vkys[k]['type'] == 's':
                    values = [v for v in self.dicv['procpar'][k]['values']]
                    self.varian_parameters[k]=values
                else:
                    raise NameError('Invalid type string')
                    
                    
if __name__ == "__main__":
    
    
    hsqc_varian_dir = r"C:\Users\ERIC\OneDrive\2-ethyl-1-indanone\Feb14154743\gHSQCAD_01.fid"
    cosy_varian_dir = r"C:\Users\ERIC\OneDrive\2-ethyl-1-indanone\Feb14154743\gCOSY_01.fid"
    hmbc_varian_dir = r"C:\Users\ERIC\Dropbox\projects\programming\2020\python\liquidNMRinterpretation\AlanKenwright\nmrData\varianData\Feb14154743\gHMBCAD_01.fid"
    proton1D_varian_dir = r"C:\Users\ERIC\Dropbox\projects\programming\2020\python\liquidNMRinterpretation\AlanKenwright\nmrData\varianData\Feb14154743\PROTON_01.fid"
   
    hsqcExpt = NMRexpt.from_varian(hsqc_varian_dir, read_blockhead=True)
    cosyExpt = NMRexpt.from_varian(cosy_varian_dir)
    hmbcExpt = NMRexpt.from_varian(hmbc_varian_dir)
    proton1dExpt = NMRexpt.from_varian(proton1D_varian_dir)
    
    proton1dExpt.produce_1D(lb = [0.001], rp=[190])
    
    proton1dExpt.plot_spectrum()
    
    # hsqcExpt.plot_spectrum( zf = [1024, 1024],
    #                         rp = [20, 210],
    #                         lp = [0,0],
    #                         lb = [0.001, 0.001],
    #                         encoding = ["states-tppi", "direct"])
    
    # cosyExpt.plot_spectrum( zf = [1024, 4096], 
    #                        lb = [0.001, 0.001], 
    #                        encoding = ["magnitude", "direct"])
    
    # hmbcExpt.plot_spectrum( zf = [1024,4096], 
    #                         lb = [0.001, 0.001], 
    #                         encoding = ["states-tppi", "direct"])
    
    
    
