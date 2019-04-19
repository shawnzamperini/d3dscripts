# -*- coding: utf-8 -*-
"""
Created on Mon Jul 02 16:02:16 2018

Plot routines for DIVIMP W impurity transport results, with collector probes

@author: Jacob
"""

import re
import numpy as np
import matplotlib as mpl
import oedge_omfit
from collections import OrderedDict
from customcolormaps import decadal2,decadal3,decadal4,decadal5,RdBkGr
from get_separatrix import get_sep_nc

class oedgeNCimp(oedge_omfit.oedgeNC):
    """
    subclass of oedgeNC class, which inherits oedge_omfit.oedgeNC methods.
    remember to call load() before running any oedgeNCimp methods.
    """
    
    def get_integrated_impurity_density(self):
        
        print 'nizs:', self.nizs
        print 'absfac:', self.absfac
        
        rawddlims=self['DDLIMS']['data']
        count = 0
        
        data = np.zeros(self.num_cells)
        
        for ir in range(self.nrs):
            for ik in range(self.nks[ir]):
                if self.area[ir,ik] != 0.0 :
                    for iz in range(self.nizs):
                        # Note: the results array indexing goes from -1:nizs for charge states. 
                        # -1 maps to 0 in the python array ... charge 0 is in index 1 - so add 
                        # one to charge to get correct indexing 
                        data[count]=data[count]+(rawddlims[iz+1][ir][ik]*self.absfac)
                    if data[count]<10**3:  #clip zero values so log plots will work
                        data[count]=10**3
                    count = count + 1
        
        print 'max n_iz:', np.max(data)
        return data
    
    def plot_integrated_impurity_density(self,zoom=None,vmin=None,vmax=None):
        
        plot_collector_probes=True    # set=True to mark collector probe locations
        plot_Wrings=True              # set=True to mark W ring locations
        plot_sep=True                # set=True to overplot separatrix
        
        data=self.get_integrated_impurity_density()
        
        fig = mpl.pyplot.figure(figsize=[6,8])
        ax = fig.add_subplot(111)
        
        if vmin==None:
            vmin=np.max(data)/10000
        else:
            vmin=vmin
        if vmax==None:
            vmax=np.max(data)
        else:
            vmax=vmax
        
        #colormap=decadal4()    # custom colormap to emphasize orders of magnitude
        colormap='nipy_spectral'
        #colormap='hsv'
        cNorm = mpl.colors.LogNorm(vmin=vmin, vmax=vmax)
        scalarMap = mpl.cm.ScalarMappable(norm=cNorm, cmap=mpl.pyplot.get_cmap(colormap,lut=21))
            
        coll = mpl.collections.PolyCollection(self.grid,array=data,cmap=scalarMap.cmap,norm=scalarMap.norm,edgecolors='none')
        ax.add_collection(coll)
        ax.plot(self.rvesm,self.zvesm,color='k',linewidth=1)
        ax.set_aspect("equal")
        ax.set_xlabel('R (m)',size=24)
        ax.set_ylabel('Z (m)',size=24)
        ax.tick_params(axis='both',which='major',labelsize=18)
        #coll.cmap.set_under('white')
        
        if zoom==None:
            ax.autoscale_view()
        elif type(zoom)==list and len(zoom)==4:
            ax.set_xlim(zoom[0],zoom[1])
            ax.set_ylim(zoom[2],zoom[3])
        else:
            print "ERROR: zoom must be either None or [x1,x2,y1,y2]"
            ax.autoscale_view()
        
        cbar=fig.colorbar(coll,ax=ax,extend='both')
        cbar.set_label('Impurity density (m^-3)',size=24)
        cbar.ax.tick_params(labelsize=18)
        
        if plot_collector_probes==True:
            # DIII-D midplane probe
            probe_mid=mpl.patches.Rectangle((2.21,-0.200),0.17,0.03)
            # DIII-D hypothetical upper probe
            probe_upper=mpl.patches.Rectangle((1.385,0.88),0.03,0.32)
            #probes=[probe_mid,probe_upper]
            probes=[probe_mid]
            pc=mpl.collections.PatchCollection(probes,facecolor='k',alpha=0.7,edgecolor=None)
            ax.add_collection(pc)
        
        if plot_Wrings:
            lines=[[(1.404,-1.25),(1.404,-1.3)],[(1.455,-1.25),(1.455,-1.3)],
                    [(1.32,-1.363),(1.32,-1.41)],[(1.37,-1.363),(1.37,-1.41)]]
            lc=mpl.collections.LineCollection(lines)
            ax.add_collection(lc)
        
        if plot_sep:
            sep=get_sep_nc(self)
            sc=mpl.collections.LineCollection(sep,color='k')
            ax.add_collection(sc)
        
        mpl.pyplot.show()
        
    def get_wall_fluxes(self,sputter_opt=7,cbomb_frac=1.0):
        """
        get plasma flux to wall, eroding impurity flux (if applicable), and eroded
        impurity flux. work in progress.
        """
        
        nds=self['NDS']['data']
        print 'wall elements:',self.nvesm
        print 'target elements:',nds
       
        idw=self['WALLPT']['data'][16]   # array of wall indices for wall elements
        idt=self['WALLPT']['data'][17]   # array of target indices for wall elements
        
        # define an along-wall coordinate s_wall
        wallcoords=self.wall
        print 's_wall starts at :',wallcoords[0]
        cumulative=0
        s_wall=np.array([0])  # first wall segment has s_wall=0
        for i in range(1,self.nvesm):
            disttonext=np.sqrt((wallcoords[i][0]-wallcoords[i-1][0])**2+(wallcoords[i][1]-wallcoords[i-1][1])**2)
            cumulative=cumulative+disttonext
            #print i,disttonext,cumulative
            s_wall=np.append(s_wall,cumulative)
        # get the distance from the final point to the first point
        distfinal=np.sqrt((wallcoords[0][0]-wallcoords[-1][0])**2+(wallcoords[0][1]-wallcoords[-1][1])**2)
        # total wall circumference
        self.circumference = s_wall[-1]+distfinal
        
        # shift s_wall such that zero is at the inner midplane (z=0)
        shift = 0.0 - wallcoords[0][1]
        s_wall = s_wall - shift
        s_wall[s_wall<0] += self.circumference
        # the conditional array s_wall_sort puts s_wall=0 first in all wall arrays, e.g. to avoid cross lines in plots
        s_wall_sort = np.argsort(s_wall)
        self.s_wall = s_wall
        self.s_wall_sort = s_wall_sort
        
        self.flxhw2=self['FLXHW2']['data'][0:self.nvesm]  # EIRENE atom+ion flux
        self.flxhw6=self['FLXHW6']['data'][0:self.nvesm]  # EIRENE atom flux
        self.flxhw8=self['FLXHW8']['data'][0:self.nvesm]  # EIRENE ion flux
        self.kflux=self['KFLUX']['data'][0:nds]
        self.kfy=self['KFY']['data'][0:nds]
        self.kyield=self['KYIELD']['data'][0:nds]
        
        # the eroding flux is kflux*cbomb_frac
        # for sputter opt 7, cbomb_frac=1, since kflux is the C flux from a previous divimp case
        # for sputter opt 1, cbomb_frac must be specified, since kflux is the D ion flux
        
        # the eroded flux is kfy
        
        self.cbombfrac=cbomb_frac
        self.cflux=self.kflux*self.cbombfrac
        
        # convert target-indexed quantities to wall-indexed
        self.kflux_w=np.zeros(self.nvesm)
        self.kfy_w=np.zeros(self.nvesm)
        self.kyield_w=np.zeros(self.nvesm)
        self.cflux_w=np.zeros(self.nvesm)
        for t in range(0,self.nvesm):
            self.kflux_w[t]=self.kflux[int(idt[t])]
            self.kfy_w[t]=self.kfy[int(idt[t])]
            self.kyield_w[t]=self.kyield[int(idt[t])]
            self.cflux_w[t]=self.cflux[int(idt[t])]
        
    def plot_wall_fluxes(self,cbomb_frac=1.0):
        """
        plot plasma flux to wall, eroding impurity flux (if applicable), and eroded
        impurity flux.
        """
        self.get_wall_fluxes(cbomb_frac=cbomb_frac)
       
        fig = mpl.pyplot.figure(figsize=[6,8])
        ax1 = fig.add_subplot(211)
        #ax2 = fig.add_subplot(212,sharex=ax1,sharey=ax1)
        
        ax1.axvline(x=self.circumference/2,color='k',ls='--')
        ax1.semilogy(self.s_wall[self.s_wall_sort],self.flxhw2[self.s_wall_sort],'r',lw=4,label='D')
        ax1.semilogy(self.s_wall[self.s_wall_sort],self.cflux_w[self.s_wall_sort],'b',lw=4,label='C')
        ax1.set_ylabel('Incident Flux (/m2/s)',size=16)
        ax1.tick_params(axis='both',which='major',labelsize=12)
        ax1.legend(loc='upper left',fontsize=14,edgecolor='w')
        ax1.set_xlabel('Poloidal dist along wall from IMP (m)',size=16)
        
        #ax2.axvline(x=self.circumference/2,color='k',ls='--')
        #ax2.semilogy(self.s_wall[self.s_wall_sort],self.kfy_w[self.s_wall_sort],'g',lw=4,label='W')
        #ax2.set_ylabel('Eroded Flux (/m2/s)',size=16)
        #ax2.set_xlabel('Poloidal dist along wall from IMP (m)',size=16)
        #ax2.tick_params(axis='both',which='major',labelsize=12)
        #ax2.legend(loc='upper left',fontsize=14,edgecolor='w')
        
        ax1.set_xlim([0,self.circumference])
        ax1.set_ylim([1.0e17,1.0e23])
        
        mpl.pyplot.subplots_adjust(hspace=0.1, right=0.95,left=0.15)
        mpl.pyplot.show()
        
    def get_deposition(self):
        
        self.wallpts=self['WALLPTS']['data']  # number of wall elements
        self.totwallse=self['WALLSE']['data'][-1] # total erosion counts
        self.totwallsi=self['WALLSI']['data'][-1] # total ion deposition counts
        self.totwallsn=self['WALLSN']['data'][-1] # total neutral deposition counts
        self.wallse=self['WALLSE']['data'][:self.wallpts]   # erosion count vs wall index
        self.wallsi=self['WALLSI']['data'][:self.wallpts]   # ion deposition count vs wall index
        self.wallsn=self['WALLSN']['data'][:self.wallpts]   # neutral deposition count vs wall index
        self.absfac_neut=self['ABSFAC_NEUT']['data'] # scaling factor
        self.pol_len=self['WALLPT']['data'][6,:self.wallpts] # poloidal length of each bin
        self.r_wall=self['RVESM']['data'][0][:self.nvesm]
        self.z_wall=self['ZVESM']['data'][0][:self.nvesm]
        self.wallse_flx=self.wallse/self.totwallse/self.pol_len*self.absfac_neut
        self.wallsi_flx=self.wallsi/(self.totwallsi+self.totwallsn)/self.pol_len*self.absfac_neut
        self.wallsn_flx=self.wallsn/(self.totwallsi+self.totwallsn)/self.pol_len*self.absfac_neut
        
    def plot_bin_deposition(self,zoombin=268):
        """
        plot ion and neutral deposition versus bin index
        """
        self.get_deposition()
        
        # plot raw counts
        maxval1=np.amax(self.wallse)
        fig1 = mpl.pyplot.figure(figsize=[6,8])
        ax1 = fig1.add_subplot(211)
        ax1.semilogy(self.wallse,'r-',label='source')
        ax1.semilogy((self.wallsi+self.wallsn),'b-',label='ion+neutral deposition')
        ax1.semilogy(self.wallsn,'k-',label='neutral deposition')
        ax1.set_ylim([maxval1/10**4,maxval1])
        ax1.set_ylabel('counts',fontsize=14)
        ax1.set_xlabel('bin index',fontsize=14)
        
        ax2 = fig1.add_subplot(212)
        ax2.semilogy(self.wallse,'r-',label='source')
        ax2.semilogy((self.wallsi+self.wallsn),'b-',label='ion+neutral deposition')
        ax2.semilogy(self.wallsn,'k-',label='neutral deposition')
        ax2.set_xlim([zoombin-8,zoombin+8])
        ax2.set_ylim([maxval1/10**4,maxval1])
        ax2.set_ylabel('counts',fontsize=14)
        ax2.set_xlabel('bin index',fontsize=14)
        
        ax1.legend(loc='upper left',fontsize=14,edgecolor='w')
        
        # plot fluxes
        maxval2=np.amax(self.wallse_flx)
        fig2 = mpl.pyplot.figure(figsize=[6,8])
        ax3 = fig2.add_subplot(211)
        ax3.semilogy(self.wallse_flx,'r-',label='source')
        ax3.semilogy((self.wallsi_flx+self.wallsn_flx),'b-',label='ion+neutral deposition')
        ax3.semilogy(self.wallsn_flx,'k-',label='neutral deposition')
        ax3.set_ylim([maxval2/10**4,maxval2])
        ax3.set_ylabel('particle flux',fontsize=14)
        ax3.set_xlabel('bin index',fontsize=14)
        
        ax4 = fig2.add_subplot(212)
        ax4.semilogy(self.wallse_flx,'r-',label='source')
        ax4.semilogy((self.wallsi_flx+self.wallsn_flx),'b-',label='ion+neutral deposition')
        ax4.semilogy(self.wallsn_flx,'k-',label='neutral deposition')
        ax4.set_xlim([zoombin-8,zoombin+8])
        ax4.set_ylim([maxval2/10**4,maxval2])
        ax4.set_ylabel('particle flux',fontsize=14)
        ax4.set_xlabel('bin index',fontsize=14)
        
        ax3.legend(loc='upper left',fontsize=14,edgecolor='w')
        
        mpl.pyplot.show()
        
    def plot_target_deposition(self,divertor='lower',vmax=None):
        """
        Plot erosion/deposition versus R for given hemisphere
        """
        plot_shelfring=True
        plot_dimes=True
        plot_OSP=True
        plot_neutraldep=False
        plot_legend=False
        
        self.get_deposition()
        
        if divertor=='lower':
            mask=np.logical_and((self.z_wall<0.0),np.logical_and((self.r_wall>1.2),(self.r_wall<1.8)))
        elif divertor=='upper':
            mask=np.logical_and((self.z_wall>0.0),np.logical_and((self.r_wall>1.2),(self.r_wall<1.8)))
        else:
            print 'Only divertor=lower or upper are supported'
            return
        
        # plot fluxes
        if not vmax:
            vmax=np.amax(self.wallse_flx)
        
        fig2 = mpl.pyplot.figure(figsize=[6,8])
        ax3 = fig2.add_subplot(211)
        ax3.plot(self.r_wall[mask],self.wallse_flx[mask],'r-',label='source')
        ax3.plot(self.r_wall[mask],(self.wallsi_flx[mask]+self.wallsn_flx[mask]),'b-',label='ion+neutral deposition')
        if plot_neutraldep:
            ax3.plot(self.r_wall[mask],self.wallsn_flx[mask],'k-',label='neutral deposition')
        ax3.set_ylim([0,vmax])
        ax3.set_ylabel('particle flux (#/m2/s)',fontsize=14)
        ax3.set_xlabel('R(m)',fontsize=14)
        
        ax4 = fig2.add_subplot(212)
        ax4.semilogy(self.r_wall[mask],self.wallse_flx[mask],'r-',label='source')
        ax4.semilogy(self.r_wall[mask],(self.wallsi_flx[mask]+self.wallsn_flx[mask]),'b-',label='ion+neutral deposition')
        if plot_neutraldep:
            ax4.semilogy(self.r_wall[mask],self.wallsn_flx[mask],'k-',label='neutral deposition')
        ax4.set_ylim([vmax/10**3,vmax])
        ax4.set_ylabel('particle flux (#/m2/s)',fontsize=14)
        ax4.set_xlabel('R(m)',fontsize=14)
        
        if plot_shelfring:
            ax3.axvline(x=1.404,color='k',ls='--',lw=1)
            ax3.axvline(x=1.455,color='k',ls='--',lw=1)
            ax4.axvline(x=1.404,color='k',ls='--',lw=1)
            ax4.axvline(x=1.455,color='k',ls='--',lw=1)
        if plot_dimes:
            ax3.axvline(x=1.466,color='k',ls='--',lw=1)
            ax3.axvline(x=1.514,color='k',ls='--',lw=1)
            ax4.axvline(x=1.466,color='k',ls='--',lw=1)
            ax4.axvline(x=1.514,color='k',ls='--',lw=1)
        if plot_OSP:
            ax3.axvline(x=1.423,color='g',lw=4,alpha=0.5)
            ax4.axvline(x=1.423,color='g',lw=4,alpha=0.5)
        
        ax3.set_xlim(1.38,1.58)
        ax4.set_xlim(1.38,1.58)
        
        if plot_legend:
            ax3.legend(loc='upper right',fontsize=14,edgecolor='w')
        
        mpl.pyplot.show()
        
    def plot_deposition_fluence(self,divertor='lower',times=[1,10,100]):
        """
        plot deposition flux multiplied by a range of times, for comparison to 
        measured deposition fluences
        """
        plot_shelfring=True
        plot_dimes=True
        plot_OSP=True
        plot_legend=True
        print_fluence=True
        
        if divertor=='lower':
            mask=np.logical_and((self.z_wall<0.0),np.logical_and((self.r_wall>1.2),(self.r_wall<1.8)))
        elif divertor=='upper':
            mask=np.logical_and((self.z_wall>0.0),np.logical_and((self.r_wall>1.2),(self.r_wall<1.8)))
        else:
            print 'Only divertor=lower or upper are supported'
            return
        
        # ion + neutral deposition flux (at/A2/s)
        depflx=(self.wallsi_flx[mask]+self.wallsn_flx[mask])/1.0e20
        
        fig = mpl.pyplot.figure(figsize=[6,5])
        ax = fig.add_subplot(111)
        for time in times:
            ax.semilogy(self.r_wall[mask],time*depflx,label='t={0}'.format(time))
        ax.set_ylim([1e-3,5e0])
        ax.set_xlim([1.4,1.6])
        ax.set_ylabel('W areal density (at/A2)',fontsize=14)
        ax.set_xlabel('Major Radius(m)',fontsize=14)
        
        if plot_shelfring:
            ax.axvline(x=1.404,color='k',ls='--',lw=1)
            ax.axvline(x=1.455,color='k',ls='--',lw=1)
        if plot_dimes:
            ax.axvline(x=1.466,color='k',ls='--',lw=1)
            ax.axvline(x=1.514,color='k',ls='--',lw=1)
        if plot_OSP:
            ax.axvline(x=1.423,color='g',lw=4,alpha=0.5)
        
        if plot_legend:
            ax.legend(loc='upper right',fontsize=14,edgecolor='w')
        
        if print_fluence:
            for time in times:
                print 't={0}'.format(time)
                print 'R(m) W(at/A2)'
                for j in range(len(depflx)):
                    print self.r_wall[mask][j], time*depflx[j]
                
        
        mpl.pyplot.show()
        
class collector_probe():
    """
    Data encoded in .collector_probe by OUT option 361 - DIVIMP collector probe model
    """
    
    def __init__(self,filename):
        self.filename=filename
        
    def find_token(self,curline,token):
        """
        Finds a token in the current line and sets the
        self.m_tok_end_idx variable to the end of the token
        and self.m_tok_start_idx to the start of the index
        """
        self.tok_start_idx = -1
        self.tok_end_idx = -1
        self.tok_start_idx = curline.find(token)
        if(self.tok_start_idx >= 0):
            self.tok_end_idx = self.tok_start_idx + len(token)
        else:
            self.tok_end_idx = -1;
        return self.tok_start_idx 
    
    def parse_file(self):
        """
        Quasi-intelligently reads contents of collector_probe file.
        Returns a list of dictionaries as probes attribute, e.g. self.probes[1]['IMPFLUX_IDF']
        accesses the IDF impurity flux profile for the 2nd probe in the file.
        """
        
        # define tokens which delimit specific parts of the collector_probe file
        TK_SUMMARY='COLLECTOR PROBE SUMMARY INFORMATION:'
        TK_LOCATION='LOCATION:'
        TK_DIAMETER_DPERP='DIAMETER,DPERP:'
        TK_ABSFAC='ABSFAC:'
        TK_HEADER='INDEX'
        
        try:
            fptr = open(self.filename,'r')
        except:
            print 'collector_probe.parse file: failed to open file: {0}'.format(self.filename)
            return
        
        # count occurences of tk_summary in order to count number of separate probes in file
        wholefile=fptr.read()
        result=re.findall(TK_SUMMARY,wholefile)
        self.num_probes=len(result)
        print 'number of probes found in file:',self.num_probes
        
        # go back to start of file, and loop over fixed structure for each probe
        fptr.seek(0)
        probes=[]
        curline=fptr.readline()
        try:
            while(len(curline)>0):
                if self.find_token(curline,TK_SUMMARY) >= 0:
                    # found a new probe, so initialize a new data dict
                    data = OrderedDict()
                if self.find_token(curline,TK_LOCATION) >= 0:
                    locline=curline.split()
                    data['R1P']=float(locline[1])
                    data['Z1P']=float(locline[2])
                    data['R2P']=float(locline[3])
                    data['Z2P']=float(locline[4])
                elif self.find_token(curline,TK_DIAMETER_DPERP) >= 0:
                    dialine=curline.split()
                    data['DIAMETER']=float(dialine[1])
                    data['DPERP']=float(dialine[2])
                elif self.find_token(curline,TK_ABSFAC) >= 0:
                    absline=curline.split()
                    data['ABSFAC']=float(absline[1])
                    data['ABSFAC_NEUT']=float(absline[2])
                elif self.find_token(curline,TK_HEADER) >= 0:
                    # found the main data block, which is of unknown length
                    header=curline.split()
                    numheader=len(header)
                    for value in header:   # initialize dictionary entries
                        data[value]=[]
                    idxline=header
                    while(len(idxline)==numheader):
                        idxline=fptr.readline().split()
                        if len(idxline)>0: 
                            for j,value in enumerate(header,0):
                                data[value].append(float(idxline[j]))
                        else:       # reached end of data block for this probe
                            for value in header: # convert all arrays to ndarrays
                                data[value]=np.array(data[value])
                            probes.append(data) # save data dict for this probe
                            break
                curline=fptr.readline()
        except Exception, Ex:
            print 'collector_probe.parse file: failed to read data, exception: {0}'.format(Ex)
            return
        fptr.close()
        
        # sanity check
        if len(probes)!=self.num_probes:
            print 'collector_probe.parse file: something went wrong in loading multiple probes'
        
        self.probes=probes
        self.absfac_applied=False
        
        # apply absfac to convert divimp units to physical fluxes/densities
        self.apply_absfac()
        
        return probes
    
    def apply_absfac(self,new_absfac=None):
        """
        multiplies reported fluxes and densities by absfac, in order to convert from divimp units
        of /part/m-tor to physical units. Absfac is typically reported by the divimp simulation.
        """
        for probe in self.probes:
            
            if new_absfac==None:
                absfac=probe['ABSFAC']
            else:
                absfac=new_absfac

            probe['IMPDENS_IDF']=probe['IMPDENS_IDF']*absfac
            probe['IMPDENS_ODF']=probe['IMPDENS_ODF']*absfac
            probe['IMPDENS_CENT']=probe['IMPDENS_CENT']*absfac
            probe['IMPFLUX_IDF']=probe['IMPFLUX_IDF']*absfac
            probe['IMPFLUX_ODF']=probe['IMPFLUX_ODF']*absfac
            probe['IMPFLUX_CENT']=probe['IMPFLUX_CENT']*absfac
            
        self.absfac_applied=True
        return
    
    def get_rminusrsep(self,probenum=0):
        """
        DIVIMP collector probe model returns PsiN and R-Rsep_OMP. sometimes it is useful to have
        R-Rsep at the probe location, which is directly comparable to the length scale of deposition
        observed on the probe itself. this function returns an array of R-Rsep_probe for the given
        probenum (index in the self.probes array), with the same length as the other probe data arrays.
        """
        psi=self.probes[probenum]['PSI']
        r=self.probes[probenum]['RSECT']
        
        # find indices of psin values that bracket psin=1
        psi_up=psi[psi > 1.00].min()
        psi_lo=psi[psi < 1.00].max()
        psi_up_index=np.argwhere(psi==psi_up)[0][0]
        psi_lo_index=np.argwhere(psi==psi_lo)[0][0]

        # find fractional index for linear interpolation
        psi_frac_index=(1.00 - psi_lo)/(psi_up - psi_lo)

        # find r values that correspond to these psin values
        r_lo=r[psi_lo_index]
        r_up=r[psi_up_index]
        r_frac=psi_frac_index * (r_up - r_lo)
        
        rsep= r_lo + r_frac
        print 'probe {0} rsep (m): {1}'.format(probenum,rsep)
        
        rminusrsep = r - rsep
        rminusrsep = rminusrsep * 100    # m to cm
        return rminusrsep
    
    def get_zminuszsep(self,probenum=0):
        """
        Same as get_rminusrsep, except z-zsep for vertical probes
        """
        psi=self.probes[probenum]['PSI']
        z=self.probes[probenum]['ZSECT']
        
        # find indices of psin values that bracket psin=1
        psi_up=psi[psi > 1.00].min()
        psi_lo=psi[psi < 1.00].max()
        psi_up_index=np.argwhere(psi==psi_up)[0][0]
        psi_lo_index=np.argwhere(psi==psi_lo)[0][0]

        # find fractional index for linear interpolation
        psi_frac_index=(1.00 - psi_lo)/(psi_up - psi_lo)

        # find z values that correspond to these psin values
        z_lo=z[psi_lo_index]
        z_up=z[psi_up_index]
        z_frac=psi_frac_index * (z_up - z_lo)
        
        zsep= z_lo + z_frac
        print 'probe {0} zsep (m): {1}'.format(probenum,zsep)
        
        zminuszsep = z - zsep
        zminuszsep = zminuszsep * 100    # m to cm
        return zminuszsep
    
    def plot_collector_probes(self,external_cp=None):
        """
        plots the ITF and OTF impurity flux to the two DIII-D collector probes (1 real, 1 virtual).
        can also specify an external collector_probe object to overplot another case.
        """
        probes=self.parse_file()
        
        rminusrsep0=self.get_rminusrsep(probenum=0)
        zminuszsep1=self.get_zminuszsep(probenum=1)
        
        fig = mpl.pyplot.figure(figsize=[6,8])
        ax1 = fig.add_subplot(211)
        ax2 = fig.add_subplot(212,sharex=ax1)
        
        # add dashed lines for another case if external_cp is specified
        if external_cp != None:
            ext = external_cp
            ext_probes = ext.parse_file()
            ext_rminusrsep0 = ext.get_rminusrsep(probenum=0)
            ext_zminuszsep1 = ext.get_zminuszsep(probenum=1)
            ax2.plot(ext_rminusrsep0,ext_probes[0]['IMPFLUX_IDF'],'r--',lw=4,label='ITF DIII-D')
            ax2.plot(ext_rminusrsep0,ext_probes[0]['IMPFLUX_ODF'],'b--',lw=4,label='OTF DIII-D')
            #ax1.plot(ext_zminuszsep1,ext_probes[1]['IMPFLUX_IDF'],'r--',lw=4,label='ITF DIII-D')
            #ax1.plot(ext_zminuszsep1,ext_probes[1]['IMPFLUX_ODF'],'b--',lw=4,label='OTF DIII-D')
            # plot midplane probe profile only here, since it represents a "good" signal
            ax1.plot(ext_rminusrsep0,ext_probes[0]['IMPFLUX_IDF'],'r--',lw=4,label='ITF DIII-D')
            ax1.plot(ext_rminusrsep0,ext_probes[0]['IMPFLUX_ODF'],'b--',lw=4,label='OTF DIII-D')
        
        # main plots
        ax2.plot(rminusrsep0,probes[0]['IMPFLUX_IDF'],'r',lw=4,label='ITF')
        ax2.plot(rminusrsep0,probes[0]['IMPFLUX_ODF'],'b',lw=4,label='OTF')
        ax2.set_ylabel('Calc W flux to probe (/m2/s)',size=16)
        ax2.set_xlabel('R-Rsep (cm)',size=16)
        ax2.tick_params(axis='both',which='major',labelsize=12)
        ax2.legend(loc='upper right',fontsize=14,edgecolor='w')
        ax2.set_title('Midplane probe: 3cm',fontsize=16)
        
        ax1.plot(zminuszsep1,probes[1]['IMPFLUX_IDF'],'r',lw=4,label='ITF')
        ax1.plot(zminuszsep1,probes[1]['IMPFLUX_ODF'],'b',lw=4,label='OTF')
        ax1.set_ylabel('Calc W flux to probe (/m2/s)',size=16)
        ax1.set_xlabel('Z-Zsep (cm)',size=16)
        ax1.tick_params(axis='both',which='major',labelsize=12)
        ax1.legend(loc='upper right',fontsize=14,edgecolor='w')
        ax1.set_title('Crown probe: 3cm',fontsize=16)
        
        #ax1.set_xlim([3,14])
        #ax1.set_xlim([2,10])
        ax1.set_xlim([0,14])
        ax1.set_ylim([0,2e21])
        ax2.set_ylim([0,2e21])
        
        mpl.pyplot.subplots_adjust(hspace=0.4)
        mpl.pyplot.show()
        
        
            
            
        