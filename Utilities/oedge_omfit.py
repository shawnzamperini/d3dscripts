# -*- coding: utf-8 -*-
"""
Created on Mon May 14 17:38:11 2018

Originally part of David Elder's OEDGE OMFIT module.
Slightly modified to work outside of OMFIT (changes marked by jhnmod).

Copied from source 5/14/18.

@author: Jacob Nichols
"""
import math
import numpy as np
import matplotlib as mpl
import netCDF4 as nc
from collections import OrderedDict
from customcolormaps import *
from get_separatrix import get_sep_nc


class oedgeNC(OrderedDict):
#class OMFIToedgeNC(OMFITnc):  jhnmod
    """OEDGE output NC data file"""

    plot_types_available = ['contour','along ring']
    # possible future plots to be added 
    # plot_types_available = ['contour','along ring','diagnostic','along wall','along target']

    # background plasma quantities
    # present for cases with a background plasma (basically everything)
    bg_plots = OrderedDict()

    bg_plots['ne']      = {'data':['KNBS'], 'targ':['KNDS'], 'label':'Density','units':'m-3',  'lineaxes':['P','S'], 'group':False,'notebook':False}
    bg_plots['Te']      = {'data':['KTEBS'],'targ':['KTEDS'],'label':'Te',     'units':'eV',   'lineaxes':['P','S'], 'group':False,'notebook':False}
    bg_plots['Ti']      = {'data':['KTIBS'],'targ':['KTIDS'],'label':'Ti',     'units':'eV',   'lineaxes':['P','S'], 'group':False,'notebook':False}
    bg_plots['Vb']      = {'data':['KVHS'], 'targ':['KVDS'], 'label':'Vpara',  'units':'m/s',  'lineaxes':['P','S'], 'group':False,'notebook':False, 'scale':'1/QTIM'}
    bg_plots['Plasma']  = {'data':['ne','Te','Ti','Vb'],     'label':'Background Plasma',      'lineaxes':['P','S'], 'group':True, 'notebook':False}

    bg_plots['E']       = {'data':['KES'],  'targ':['KEDS'], 'label':'Epara',  'units':'V/m',  'lineaxes':['P','S'], 'group':False,'notebook':False}
    bg_plots['ExB Pol'] = {'data':['E_POL'],'targ':None,     'label':'Epol',   'units':'V/m',  'lineaxes':['P','S'], 'group':False,'notebook':False}
    bg_plots['ExB Rad'] = {'data':['E_RAD'],'targ':None,     'label':'Erad',   'units':'V/m',  'lineaxes':['P','S'], 'group':False,'notebook':False}
    bg_plots['Efields'] = {'data':['E','ExB Pol','ExB Rad'], 'label':'Electric Fields',        'lineaxes':['P','S'], 'group':True, 'notebook':False}

    bg_plots['V_pol']   = {'data':['EXB_P'],'targ':None,     'label':'Vpol',   'units':'m/s',  'lineaxes':['P','S'], 'group':False,'notebook':False}
    bg_plots['V_rad']   = {'data':['EXB_R'],'targ':None,     'label':'Vrad',   'units':'m/s',  'lineaxes':['P','S'], 'group':False,'notebook':False}
    bg_plots['Drifts']  = {'data':['V_pol','V_rad'],         'label':'Drifts',                 'lineaxes':['P','S'], 'group':True, 'notebook':False}


    # Data present for cases that ran EIRENE to generate hydrogenic quantities
    h_plots=OrderedDict()
    h_plots['H Power']      = {'data':['HPOWLS'], 'targ':None, 'label':'Hydrogen Radiated Power','units':'W/m3',   'lineaxes':['P','S'], 'group':False,'notebook':False}
    h_plots['H Line']       = {'data':['HLINES'], 'targ':None, 'label':'Hydrogen Line Radiation','units':'W/m3',   'lineaxes':['P','S'], 'group':False,'notebook':False}
    h_plots['H Dalpha']     = {'data':['PINALP'], 'targ':None, 'label':'Hydrogen Dalpha Emission','units':'ph/m3/s','lineaxes':['P','S'], 'group':False,'notebook':False}
    h_plots['H Ionization'] = {'data':['PINION'], 'targ':None, 'label':'Hydrogen Ionization','units':'/m3/s','lineaxes':['P','S'], 'group':False,'notebook':False}
    h_plots['H Recomb']     = {'data':['PINREC'], 'targ':None, 'label':'Hydrogen Recombination','units':'/m3/s','lineaxes':['P','S'], 'group':False,'notebook':False}
    h_plots['H2 Density']   = {'data':['PINMOL'], 'targ':None, 'label':'Hydrogen Molecule Density','units':'/m3','lineaxes':['P','S'], 'group':False,'notebook':False}
    h_plots['H0 Density']   = {'data':['PINATO'], 'targ':None, 'label':'Hydrogen Atom Density','units':'/m3','lineaxes':['P','S'], 'group':False,'notebook':False}
    h_plots['H Atom Temp']  = {'data':['PINENA'], 'targ':None, 'label':'Hydrogen Atom Temperature','units':'eV','lineaxes':['P','S'], 'group':False,'notebook':False}
    h_plots['H Mol Temp']   = {'data':['PINENM'], 'targ':None, 'label':'Hydrogen Molecule Temperature','units':'eV','lineaxes':['P','S'], 'group':False,'notebook':False}
    h_plots['H Ion Energy Loss']  = {'data':['PINQI'], 'targ':None, 'label':'Hydrogen-Ion Energy Loss Term','units':'W/m3(?)','lineaxes':['P','S'], 'group':False,'notebook':False}
    h_plots['H Elec Energy Loss'] = {'data':['PINQE'], 'targ':None, 'label':'Hydrogen-Electron Energy Loss Term','units':'W/m3(?)','lineaxes':['P','S'], 'group':False,'notebook':False}
    h_plots['H Quantities'] = {'data':['H0 Density','H2 Density','H Ionization','H Recomb'], 'label':'Hydrogen Quantities','lineaxes':['P','S'], 'group':False,'notebook':False}

    # Data present for cases that ran impurities
    imp_plots = OrderedDict()
    imp_plots['Imp Density']     = {'data':['DDLIMS'], 'targ':None, 'label':'Impurity Density','units':'/m3','lineaxes':['P','S'], 'group':False,'notebook':False, 'scale':'ABSFAC'}
    imp_plots['Imp Temperature'] = {'data':['DDTS'], 'targ':None, 'label':'Impurity Temperature','units':'eV','lineaxes':['P','S'], 'group':False,'notebook':False}
    imp_plots['Imp Ionization']  = {'data':['TIZS'], 'targ':None, 'label':'Impurity Ionization','units':'/m3/s','lineaxes':['P','S'], 'group':False,'notebook':False, 'scale':'ABSFAC'}
    imp_plots['Imp Radiated Power']  = {'data':['POWLS'], 'targ':None, 'label':'Impurity Radiated Power','units':'W/m3','lineaxes':['P','S'], 'group':False,'notebook':False, 'scale':'ABSFAC'}
    imp_plots['Impurity Quantities']  = {'data':['Imp Density','Imp Temperature','Imp Ionization','Imp Radiated Power'], 'label':'Impurity Quantities','lineaxes':['P','S'], 'group':True,'notebook':False}

    # Plots along surfaces (?) separate category or make a different dict entry under appropriate sections? 
    surface_plots = OrderedDict()

    # Synthetic diagnostic plots are a separate category since they will require special treatment
    diagnostic_plots = OrderedDict()

    # available plot types - contour/false colour plots, line plots along ring, line plots along surfaces (use different plotting axes)
    plot_types=['contour-polygon','contour-interpolate','along ring','surface']
    #default_contour = 'contour-polygon'

    all_plots = OrderedDict()

    all_plots.update(bg_plots)
    all_plots.update(h_plots)
    all_plots.update(imp_plots)
    all_plots.update(surface_plots)

    plot_categories = ['background','hydrogen','impurity','surface','diagnostic']


    def __init__(self, filename,**kw):
#        OMFITnc.__init__(self,filename,**kw)
#        self.dynaLoad=True
        
        # jhnmod: initialize object as a blank ordered dictionary
        OrderedDict.__init__(self)
        self.filename=filename
        
#    @dynaLoad
    def load(self):

        # empty the dictionary in case this is a reload
        self.clear()    

        #print 'Load netcdf:'
#        OMFITnc.load(self)   jhnmod

        # jhnmod: the following section reproduces the basic functionality of OMFITnc
        with nc.Dataset(self.filename,'r') as f:
            # report data formats
            print 'file format: ', f.file_format
            print 'disk format: ', f.disk_format
            
            #load all netcdf variables into OMFIT-equivalent dictionaries
            for i in f.variables.keys():
                #print i
                var=f.variables[i][:]
                self[i]={'data':var}
                # TODO: load metadata like units, dims, etc.
        
        
        #
        # Do some post-processing after loading 
        # 
        # add anything else it should automatically do with the data on load ...
        # should these be done in the load routine or in the constructor?
        # Presumably the constructor calls load()

        #print 'Load self variables:'
        
        self.load_simulation_data()
        self.plots_available()        




    def load_simulation_data(self):
        # This routine labels specific geometry data from the NC file into local variables
        # This data is frequently used in the plots .. since it appears that one case object
        # handles all the plotting then it is worthwhile adding some efficiencies
        
        # Certain data must be present to make plots possible ... if it is not available for 
        # some reason then flag it as an error. 

        # load wall and create path object to define inside and outside the vessel
        # This is used in the tricontourf plots to eliminate triangles outside the vessel
        debug = True

        # vessel wall definition
        # jhnmod: switch to nvesm (full wall) from nves (wall w/o grid intersections)
        self.nvesm = self['NVESM']['data'] 
        self.rvesm = self['RVESM']['data'][:][:self.nvesm]
        self.zvesm = self['ZVESM']['data'][:][:self.nvesm]
        self.wall = [[self.rvesm[0][i],self.zvesm[0][i]] for i in range(self.nvesm)]
        self.wall_path = mpl.path.Path(self.wall)


        # load the computational mesh polygon data into a PolyCollection for plotting code results 
        # remove zero volume polygons from the list and also exclude when loading simulation data

        self.nvertp = self['NVERTP']['data']
        self.rvertp = self['RVERTP']['data']
        self.zvertp = self['ZVERTP']['data']
        self.korpg = self['KORPG']['data']

        # cell centers
        self.rs = self['RS']['data']
        self.zs = self['ZS']['data']

        # Significant surfaces in the computational mesh
        self.nrs = self['NRS']['data']
        self.nks = self['NKS']['data']
        self.irsep = self['IRSEP']['data']
        # no irsep2 in netcdf
        #self.irsep2 = self['IRSEP2']['data']
        self.irwall = self['IRWALL']['data']
        self.irwall2 = self['IRWALL2']['data']
        self.irtrap = self['IRTRAP']['data']
        self.irtrap2 = self['IRTRAP2']['data']

        # grid connection map data
        self.ikins = self['IKINS']['data']
        self.ikouts = self['IKOUTS']['data']
        self.irins = self['IRINS']['data']
        self.irouts = self['IROUTS']['data']

        # properties of the grid and target
        self.area= self['KAREAS']['data']
        
        # idds is an index array mapping from (targ,ir) -> (target index)
        # However, the index is from fortran arrays and numbers from 1
        # so when mapped to python arrays you need to subtract 1 which 
        # should work if idds is a numpy array
        self.idds= self['IDDS']['data'].copy() -1
        self.smax= self['KSMAXS']['data']
        self.pmax= self['KPMAXS']['data']
        self.kss= self['KSS']['data']
        self.kps= self['KPS']['data']


        mesh = []
        count = 0
        debug = True

        for ir in range(self.nrs):
            for ik in range(self.nks[ir]):
                index = self.korpg[ir,ik]-1
                if self.area[ir,ik] != 0.0:
                    mesh.append(zip(self.rvertp[index][0:4],self.zvertp[index][0:4]))
                    count = count +1

                # check to see that the cell center point is inside the set of vertices for the cell
                # this verifies that the grid is loading correctly
                if debug and self.area[ir,ik] != 0.0:
                    vert = zip(self.rvertp[index][0:4],self.zvertp[index][0:4])     
                    cell = mpl.path.Path(vert)
                    r = self['RS']['data'][ir,ik]
                    z = self['ZS']['data'][ir,ik]
                    if not cell.contains_point([r,z]):
                        print("grid loading error:")
                        print("vert:",vert)
                        print("test:",ir,ik,r,z,cell.contains_point([r,z]))

        # save the number of non-zero volume cells in the mesh
        # save the mesh polygons as in grid
        self.num_cells = count
        self.grid = mesh

        # The grid is used to create a polygon plot
        #coll = PolyCollection(self.grid,array=data,cmap=pyplot.cm.jet,edgecolors='none')

        #
        # The above created a mesh of polygons from the computational grid
        # Also needed is a mesh composed of the cell center points both with and without the 
        # target coordinates. These meshes have two uses:
        # 1) Interpolated triangular meshes which are based on cell center and target data
        # 2) Interpolation of the OEDGE results to generate an ERO plasma file. 
        #
        
        self.num_points = self.num_cells + (self.nrs -self.irsep +1) * 2

        self.r_cells = np.zeros(self.num_cells)
        self.z_cells = np.zeros(self.num_cells)

        self.r_points = np.zeros(self.num_points)
        self.z_points = np.zeros(self.num_points)

        cp = 0
        cc = 0

        for ir in range(self.nrs):
            for ik in range(self.nks[ir]):
        
                if ik == 0 and ir >= self.irsep-1 :
                    # add first target point 
                    id0 = self.idds[1,ir]
                    self.r_points[cp] = self['RP']['data'][id0]
                    self.z_points[cp] = self['ZP']['data'][id0]
                    cp = cp + 1

                if self.area[ir,ik] != 0.0 :
                    self.r_points[cp] = self.rs[ir][ik]
                    self.z_points[cp] = self.zs[ir][ik]
                    cp = cp + 1
                    self.r_cells[cc] = self.rs[ir][ik]
                    self.z_cells[cc] = self.zs[ir][ik]
                    cc = cc + 1

                if ik == self.nks[ir]-1 and ir >= self.irsep-1:
                    # add last target point
                    id1 = self.idds[0,ir]
                    self.r_points[cp] = self['RP']['data'][id1]
                    self.z_points[cp] = self['ZP']['data'][id1]
                    cp = cp + 1

        # Impurity scaling factor if available 
        # All impurity results are stored scaled to 1 part/m-tor/s entering the system
        # to get absolute values the results are multipled by absfac (absolute scaling factor)
        if 'ABSFAC' in self:
            self.absfac = self['ABSFAC']['data']
        else:
            self.absfac = 1.0

        # Simulation time steps
        # ion time step
        if 'QTIM' in self:
            self.qtim = self['QTIM']['data']
        else:
            self.qtim = 1.0

        # neutral time step
        if 'FSRATE' in self:
            self.fsrate = self['FSRATE']['data']
        else:
            self.fsrate = 1.0
        
        #print 'qtim:',self.qtim,self.fsrate


        if 'NIZS' in self:
            self.nizs = self['NIZS']['data']
        else:
            self.nizs = None

        # jhnmod: drift normalization
        self.driftfact=np.sqrt(self['KBFS']['data']**2-1) # geometric factor for ExB poloidal drift
        if 'EXB_P' in self:
            self['EXB_P_n']={'data':self['EXB_P']['data']/self.qtim/self.driftfact}
        if 'EXB_R' in self:
            self['EXB_R_n']={'data':self['EXB_R']['data']/self.qtim}

    def plots_available(self):
        # this routine returns a list of quantities that this object can plot
        # background plasma plots
        self.plots_avail=[]
        self.bg_plots_avail=[]
        self.h_plots_avail=[]
        self.imp_plots_avail=[]
        self.surf_plots_avail=[]
        
        for plotkey in oedgeNC.all_plots.keys():
            # Check to see that all the required datasets are available in the NC file for each class of plots
            avail = True                
            for datakey in oedgeNC.all_plots[plotkey]['data']:
                dataname,targname,label,scalef = self.get_names(plotkey,datakey)
                if dataname not in self:
                   avail = False
                if targname is not None:
                    if targname not in self:
                        avail = False
                
            # if all elements of a group plot are available then that plot is available
            if avail:
                self.plots_avail.append(plotkey)
                if plotkey in oedgeNC.bg_plots:
                    self.bg_plots_avail.append(plotkey)
                elif plotkey in oedgeNC.h_plots:
                    self.h_plots_avail.append(plotkey)
                elif plotkey in oedgeNC.imp_plots:
                    self.imp_plots_avail.append(plotkey)
                elif plotkey in oedgeNC.surface_plots:
                    self.surf_plots_avail.append(plotkey)



    def get_plot_categories(self):
        # return a list of the available plot categories
        return self.plot_categories



    def get_plots_available(self,kind=None):
        # return a list of the plots that are available in the current dataset
        # Since the plot list is getting long - allow for splitting by type so that the UI
        # can have more concise drop down lists
        #print "get_plots_available:",kind
        if kind is None:
            return self.plots_avail
        elif kind == 'background':
            return self.bg_plots_avail
        elif kind == 'hydrogen':
            return self.h_plots_avail
        elif kind == 'impurity':
            return self.imp_plots_avail
        else:
            return []


    def get_plot_types(self):
        # this routine returns a list of the supported plot types
        return self.plot_types


    def get_along_ring_axis_types(self,plot_select):
        # look up the selected plot and return the along ring or linear axis types available
        if plot_select in oedgeNC.all_plots:
            return oedgeNC.all_plots[plot_select]['lineaxes']
        # if this fails for some reason then return None ... but it shouldn't
        return None


    def need_ionization_state(self,selection):
        # Check to see if a charge state is needed for the specified plot
        res = False
        if selection in oedgeNC.imp_plots:
            res = True
        elif selection in {'H Power','H Lines'}:
            res=True
        return res


    def need_ring_number(self,selection):
        # Check to see if a ring number is needed for specified plot type
        res = False
        if selection == 'along ring':
            res = True
        return res



    def get_data_2d(self,dataname,targname=None,charge=None,scalef=1.0):
        # load a set of 2D data
        #

        rawdata=self[dataname]['data']
        count = 0
        #
        # Note: targname should always be None for cases where charge is specified
        # charge is a leading index into the data array mapping. 

        # if targname is specified then add the target data at the ends of the rings
        if targname is None:
            # leave out target data when plotting a polygon mesh
            data = np.zeros(self.num_cells)
        else:
            # add space to include target values in the mesh - used in triangulation interpolation plots
            data = np.zeros(self.num_points)

        # loop through data arrays
        # Note ... all of this is required because the rows of the array are of inconsistent lengths .. it isn't a rectangular grid
        # in terms of indexing ... core, sol and pfz rings may all have different numbers of cells though nsol_cells=ncore+npfz is usually true
        # along a ring
        count = 0
        for ir in range(self.nrs):
            for ik in range(self.nks[ir]):                
                if targname is not None and ik == 0 and ir >= self.irsep-1:
                    # add first target point 
                    id0 = self.idds[1,ir]
                    data[count] = self[targname]['data'][id0]
                    count = count + 1

                if self.area[ir,ik] != 0.0 :
                    # include index for charge state if specified - this requires an offset 
                    if charge is None:
                        data[count] = rawdata[ir][ik] * scalef
                    else:
                        # Note: the results array indexing goes from -1:nizs for charge states. 
                        # -1 maps to 0 in the python array ... charge 0 is in index 1 - so add one to charge
                        # to get correct indexing 
                        data[count] = rawdata[charge+1][ir][ik] * scalef

                    count = count + 1        

                if targname is not None and ik == self.nks[ir]-1 and ir >= self.irsep-1 :
                    # add last target point
                    id1 = self.idds[0,ir]
                    data[count] = self[targname]['data'][id1]
                    count = count + 1


        return data



    def plot(self,plot_select,plot_type,axis_type='S',ring_range=[],charge_range=[],zoom=None):
        # produce a contour plot 
        # Either interpolated or polygon

        #print 'plot inputs1:',plot_select,plot_type,axis_type,ring_range,charge_range,zoom


        # Set inputs not needed for specific plots to default values
        # Argument ELIMINATION doesn't propagate ... while argument object MANIPULATION would propagate back to calling routine
        # This is needed since I am not updating the GUI when updating ring or charge ranges ... so these variables keep
        # values from previous plot calls that may not be relevant
        # I could just reload the GUI to avoid this but I am not sure the overhead is worth it
        if not self.need_ring_number(plot_type):
            ring_range=[]

        if not self.need_ionization_state(plot_select):
            charge_range=[]
        else:   # plot needs ionization state
            if charge_range==[]:
                print("ERROR: Charge states not specified correctly for charge-resolved plots")
                return

        #print 'plot inputs2:',plot_select,plot_type,axis_type,ring_range,charge_range,zoom
        

        if plot_type == 'contour-polygon' or plot_type == 'contour-interpolate':
            self.plot_contour(plot_select,plot_type,charge_range,zoom)

        elif plot_type == 'along ring':
            if len(ring_range)==0 :
                print("ERROR: Rings not specified correctly for along ring plots")
            else:
                self.plot_along_ring(plot_select,axis_type,ring_range,charge_range,zoom)

        elif plot_type == 'surface':
            print("Plot not yet implemented")
            return

        elif plot_type == 'diagnostic':
            print("Plot not yet implemented")
            return

        else:
            print("ERROR: Specified plot type: `%s` not found."%plot_type)
            
        #TODO: plot range and log options, colormap options. Use option object?



    def plot_contour(self,plot_select,plot_type,charge_range=[],zoom=None):
        # This produces polygon contour plots matching the
        # computationals mesh of the data included in the plot list

        # If nplots is 1 then the 'data' reference contains a self[] reference otherwise it contains a list of all_plot dictionary entries        
        # set up subplots

        #
        # If charge states are specified for an impurity quantity two things happen
        # 1) The figures go into a figure notebook - one/charge state
        # 2) The get_data_2d routine in the specific contour plotting routines will require the charge state to be specified. 
        #
        if len(charge_range) == 0:
            # if no charge states then just plot the figure
            # Note: OMFITx.py Figure is intended to put a figure inside an OMFIT GUI element
            #fig = Figure(returnFigure=True)        
            fig = mpl.pyplot.figure()
            self.plot_contour_fig(fig,plot_select,plot_type,charge=None,zoom=zoom)
        else:
#            # use a figure notebook
#            # For lists of charge states 
#            fig_notebook = FigureNotebook(nfig=0,name='Contour plots by charge state')
#            # loop through charge states
#            for charge in charge_range:
#                fig = fig_notebook.add_figure(label='IZ={}'.format(charge))
#                self.plot_contour_fig(fig,plot_select,plot_type,charge,zoom=zoom)
#                fig.canvas.draw()

            # jhnmod: use multiple figures for charge states rather than OMFIT FigureNotebook
            for charge in charge_range:
                fig = mpl.pyplot.figure()
                self.plot_contour_fig(fig,plot_select,plot_type,charge=charge,zoom=zoom)


    def plot_contour_fig(self,fig,plot_select,plot_type,charge=None,zoom=None):
        # given the figure - add the appropriate contour plots

        plot_list = oedgeNC.all_plots[plot_select]['data']
        nplts = len(plot_list)

        # get figure layout for number of plots on the page
        nrows,ncols = self.calculate_layout(nplts)

        for cnt in range(len(plot_list)):
            ax = fig.add_subplot(nrows,ncols,cnt+1)
            dataname,targname,label,scalef = self.get_names(plot_select,plot_list[cnt])
            
            if plot_type == 'contour-polygon':
                #print "plot_contour_polygon_call:",scalef
                self.plot_contour_polygon(fig,ax,dataname,charge,zoom,scalef=scalef)
                ax.set_title(label)

            elif plot_type == 'contour-interpolate':
                self.plot_contour_interpolate(fig,ax,dataname,targname,charge,zoom,scalef=scalef)
                ax.set_title(label)
            # set the figure title
            if charge==None:
                fig.suptitle(plot_select)
            else:
                fig.suptitle('{0} Z={1}'.format(plot_select,charge))




    def plot_contour_polygon(self,fig,ax,dataname,charge=None,zoom=None,scalef=1.0):
        # polygon plots never use the target data
        # produce polygon contour plot
        if charge is None:
            data = self.get_data_2d(dataname,scalef=scalef)
        else:
            data = self.get_data_2d(dataname,targname=None,charge=charge,scalef=scalef)
            
        coll = mpl.collections.PolyCollection(self.grid,array=data,cmap=mpl.pyplot.cm.jet,edgecolors='none')
        ax.add_collection(coll)
        ax.plot(self.rvesm,self.zvesm,color='k',linewidth=1)
        #ax.plot(wall,linewidth=1)
        ax.autoscale_view()
        ax.set_aspect("equal")
        fig.colorbar(coll,ax=ax)



    def plot_contour_interpolate(self,fig,ax,dataname,targname=None,charge=None,zoom=None,scalef=1.0):
        # produce interpolated contour plot
        # not all 2D data has target data
        
        if charge is None:
            data = self.get_data_2d(dataname,targname,scalef=scalef)
        else:
            data = self.get_data_2d(dataname,targname,charge,scalef=scalef)

        if targname == None:
            r = self.r_cells
            z = self.z_cells
        else:
            r = self.r_points
            z = self.z_points

        triang = mpl.tri.Triangulation(r, z)

        # Mask off unwanted triangles by finding ones with center points outside of wall
        # This could be changed to ones with any corner outside wall but then you might get
        # missing triangle gaps inside the wall
        rmid = r[triang.triangles].mean(axis=1)
        zmid = z[triang.triangles].mean(axis=1)
        points = [[rmid[i],zmid[i]] for i in range(len(rmid))]

        res = self.wall_path.contains_points(points)
        mask = np.where(res==False,1,0)
        #print "mask",mask

        # triangles to be plotted should be set now
        triang.set_mask(mask)

        #ax.set_cmap("jet")
        ax.set_aspect('equal')
        ax.plot(self.rvesm,self.zvesm,color='k',linewidth=1)
        # contour levels are set to 20 but this could be pulled into a settable parameter
        cs=ax.tricontourf(triang,data,20,cmap=mpl.pyplot.cm.jet)
        fig.colorbar(cs,ax=ax)

    def plot_contour_standalone(self,dataname,charge=None,zoom=None,vmin=None,vmax=None,scalef=1.0,normtype='Lin'):
        # jhnmod: added this function to make contour-polygon plots independently of
        # primary plot dictionary. Also serves as a test bed for graphics option development.
        
        plot_Wrings=True
        plot_sep=True                # set=True to overplot separatrix
        iaea_plot=False  # special settings for Unterberg IAEA 2018
        
        if charge is None:
            data = self.get_data_2d(dataname,scalef=scalef)
        else:
            data = self.get_data_2d(dataname,targname=None,charge=charge,scalef=scalef)
        
        fig = mpl.pyplot.figure()
        ax = fig.add_subplot(111)
        
        maxval=np.max(data)
        minval=np.min(data)
        levels=None
        
        if normtype=='Lin':   # single-signed linear scaling
            if np.abs(maxval)>np.abs(minval):  # positive data
                newmax=maxval
            else:              # negative data
                newmax=minval
            if not vmax:
                vmax=newmax
            if not vmin:
                vmin=0
                
#            colormap=decadal5()
#            colormap=yesno()
#            colormap='viridis'
            colormap='nipy_spectral'
            
#            cNorm = mpl.colors.Normalize(vmin=0, vmax=newmax)
#            cNorm = mpl.colors.Normalize(vmin=minval, vmax=newmax)
#            cNorm = mpl.colors.Normalize(vmin=1.0e19, vmax=3.0e19)
            cNorm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
            scalarMap = mpl.cm.ScalarMappable(norm=cNorm, cmap=mpl.pyplot.get_cmap(colormap,lut=20))

        elif normtype=='Log':   # single-signed logarithmic scaling
            # no negative numbers allowed for now
            logmax=np.log10(maxval)
            newmax=maxval
 #           newmax=np.power(10,np.ceil(logmax))
 #           newmax=np.power(10,np.floor(logmax))           
            if not vmax:
                vmax=newmax
            if not vmin:
                vmin=newmax/1000
                
 #           colormap='viridis'
            colormap='nipy_spectral'
 #           colormap=decadal4()
#            colormap='seismic'
 
#            cNorm = mpl.colors.LogNorm(vmin=10**-2, vmax=10**0)
#            cNorm = mpl.colors.LogNorm(vmin=1E-6, vmax=1E-2)
  #          cNorm = mpl.colors.LogNorm(vmin=10**0, vmax=10**4)
#            cNorm = mpl.colors.LogNorm(vmin=10**-2, vmax=10**2)
            cNorm = mpl.colors.LogNorm(vmin=vmin, vmax=vmax)
            scalarMap = mpl.cm.ScalarMappable(norm=cNorm, cmap=mpl.pyplot.get_cmap(colormap,lut=20))
            
        elif normtype=='SymLin':   # diverging linear scaling
            absval=np.max(np.abs([maxval,minval]))
            if not vmax:
                vmax=absval
            if not vmin:
                vmin=-absval
#            colormap='seismic'
#            colormap=RdBkGr()
            colormap='coolwarm'
            
            cNorm = mpl.colors.Normalize(vmin=vmin,vmax=vmax)
            scalarMap = mpl.cm.ScalarMappable(norm=cNorm, cmap=mpl.pyplot.get_cmap(colormap,lut=21))
            
        elif normtype=='SymLog':   # diverging logarithmic scaling
            absval=np.max(np.abs([maxval,minval]))
            if not vmax:
                vmax=absval
            if not vmin:
                vmin=-absval
                
            if iaea_plot:
                thresh=1./20
                num=11
                colormap='coolwarm'
                levels=np.concatenate((np.geomspace(vmin,vmin*thresh,num/2),np.geomspace(vmax*thresh,vmax,num/2)))
                print levels
            else:
                thresh=1./100
                num=21
                colormap='coolwarm'

            cNorm = mpl.colors.SymLogNorm(linthresh=vmax*thresh,linscale=1,vmin=vmin,vmax=vmax)

            scalarMap = mpl.cm.ScalarMappable(norm=cNorm, cmap=mpl.pyplot.get_cmap(colormap,lut=num))
            
        else:
            print "ERROR: normtype unsupported"
            
        coll = mpl.collections.PolyCollection(self.grid,array=data,cmap=scalarMap.cmap,norm=scalarMap.norm,edgecolors='none')
#        coll = mpl.collections.PolyCollection(self.grid,array=data,cmap=mpl.pyplot.cm.jet,edgecolors='none')
        ax.add_collection(coll)
        ax.plot(self.rvesm,self.zvesm,color='k',linewidth=1)
        ax.set_aspect("equal")
        
        if zoom==None:
            ax.autoscale_view()
        elif type(zoom)==list and len(zoom)==4:
            ax.set_xlim(zoom[0],zoom[1])
            ax.set_ylim(zoom[2],zoom[3])
        else:
            print "ERROR: zoom must be either None or [x1,x2,y1,y2]"
            ax.autoscale_view()
        
        ax.set_xlabel('R (m)')
        ax.set_ylabel('Z (m)')
        
        if plot_Wrings:
            lines=[[(1.404,-1.25),(1.404,-1.3)],[(1.455,-1.25),(1.455,-1.3)],
                    [(1.32,-1.363),(1.32,-1.41)],[(1.37,-1.363),(1.37,-1.41)]]
            lc=mpl.collections.LineCollection(lines)
            ax.add_collection(lc)
        
        if plot_sep:
            sep=get_sep_nc(self)
            sc=mpl.collections.LineCollection(sep,color='k')
            ax.add_collection(sc)
        
        if levels is not None:
            fig.colorbar(coll,ax=ax,boundaries=levels,ticks=levels)
        else:
            fig.colorbar(coll,ax=ax,extend='both')

        
    def calculate_layout(self,nplts):
        #
        # Calculate a layout that will hold all the figures
        #
        # Handle special cases - then use generic algorithm
        # 
        if nplts == 3:
            nrows = 1
            ncols = 3
        elif nplts == 7:
            nrows = 2
            ncols = 4
        else:
            nrows = math.ceil(math.sqrt(nplts))
            ncols = math.ceil(math.sqrt(nplts))        
            # if number of plots will fit in a smaller grid remove a row
            if ncols*(nrows-1) >= nplts:
                nrows = nrows -1

        # return the layout
        return nrows,ncols


    def get_names(self,plot_select,plot_name):
        # This returns the names from the NC dictionary to be plotted based on the oedgeNC.all_plots summary
        # It also returns the label associated with the individual plot element
        dataname = None
        targname = None
        label = None
        scale_desc = None

        # each of the dictionary entries contains a list - this list is either a single reference to 
        # the NC datafile or a list of reference to other plots. Plot_name should always be a singular plot if it is present

        if plot_select in oedgeNC.all_plots:
            if plot_name in self:
                # Need to return plot_name as dataname and get targname from the plot select reference
                dataname = plot_name
                if 'targ' in oedgeNC.all_plots[plot_select]:
                    if oedgeNC.all_plots[plot_select]['targ'] is not None:
                        targname = oedgeNC.all_plots[plot_select]['targ'][0]
                # get label
                if 'label' in oedgeNC.all_plots[plot_select]:
                    if oedgeNC.all_plots[plot_select]['label'] is not None:
                        label = oedgeNC.all_plots[plot_select]['label']
                        if 'units' in  oedgeNC.all_plots[plot_select]:
                            if oedgeNC.all_plots[plot_select]['units'] is not None:
                                label = label + ' ' +  oedgeNC.all_plots[plot_select]['units']
                # get scaling factor if any
                if 'scale' in oedgeNC.all_plots[plot_select]:
                    if oedgeNC.all_plots[plot_select]['scale'] is not None:
                        scale_desc = oedgeNC.all_plots[plot_select]['scale']

            elif plot_name in oedgeNC.all_plots:
                # Need to load both dataname and targname from the plot list
                if 'data' in  oedgeNC.all_plots[plot_name]:
                    if oedgeNC.all_plots[plot_name]['data'] is not None:
                        dataname = oedgeNC.all_plots[plot_name]['data'][0]

                if 'targ' in oedgeNC.all_plots[plot_name]:
                    if oedgeNC.all_plots[plot_name]['targ'] is not None:
                        targname = oedgeNC.all_plots[plot_name]['targ'][0]
                # get data label
                if 'label' in  oedgeNC.all_plots[plot_name]:
                    if oedgeNC.all_plots[plot_name]['label'] is not None:
                        label = oedgeNC.all_plots[plot_name]['label']
                        if 'units' in  oedgeNC.all_plots[plot_name]:
                            if oedgeNC.all_plots[plot_name]['units'] is not None:
                                label = label + ' ' +  oedgeNC.all_plots[plot_name]['units']
                # get scaling factor if any
                if 'scale' in oedgeNC.all_plots[plot_name]:
                    if oedgeNC.all_plots[plot_name]['scale'] is not None:
                        scale_desc = oedgeNC.all_plots[plot_name]['scale']

#
#            Comment out - this code is used both to get plots AND to check whether they exist in the nc file
#
#            else:
#                print "ERROR: plot_name must be either a reference to self or a plot in oedgeNC.all_plots",plot_select,plot_name

        # set numerical scaling factor for the volumetric data
        # Scaling factor is set 1.0 by default ... only changed if some other quantity is specified
        if scale_desc == 'ABSFAC':
            scalef = self.absfac
        elif scale_desc == '1/QTIM':
            scalef = 1.0/self.qtim
        else:
            scalef = 1.0

        return dataname,targname,label,scalef



    def get_data_along_ring(self,ir,dataname=None,targname=None,charge=None,scalef=1.0):
        # loads a set of 2D data        
        # get along ring data
        #print "get_data_along_ring:",ir,dataname,targname,charge

        if dataname is not None:
            if charge is None:
                data = self[dataname]['data'][ir][:self.nks[ir]] * scalef
            else:
                # Note: the results array indexing goes from -1:nizs for charge states. 
                # -1 maps to 0 in the python array ... charge 0 is in index 1 - so add one to charge
                # to get correct indexing 
                data = self[dataname]['data'][charge+1][ir][:self.nks[ir]] * scalef

        #print "data:",dataname,data

        # Insert the target data at the beginning and end when it is available
        # Targt data is only available on sol and pfz rings
        if targname is not None and ir >= self.irsep-1:
            # Target data does not need the scaling
            # get target indices in the target data arrays
            id0 = self.idds[1,ir]
            id1 = self.idds[0,ir]

            #print 'id0,id1:',id0,id1,self.idds[1][ir],self.idds[0][ir]
            #print 'idds:',self.idds

            #print 'targ:',targname,self[targname]['data']
            #print 'targ2 id:',targname,self.idds[0][:]
            #print 'targ1 id:',targname,self.idds[1][:]
            #print 'targ data1:',targname,id0,self[targname]['data'][id0]
            

            data=np.insert(data,0,[self[targname]['data'][id0]])

            #print 'targ data2:',targname,id1,self[targname]['data'][id1]

            data=np.append(data,[self[targname]['data'][id1]])
            #print 'data2:',data


        # return the data
        return data

        

    def get_axis(self,ir,axis_type='S',targname=None):
        # This loads the particular axis for the along ring plot (either S or P)
        # It also adds the target elements if targname is not None
        if axis_type == 'P':
            axis = self.kps[ir][:self.nks[ir]]
            label = 'P (m)'
            if targname is not None and ir >self.irsep-1:
                axis = np.insert(axis,0,0.0)
                axis = np.append(axis,self.pmax[ir])
        else:    #  if axis_type == 'S':  # if someone forgets to set axis_type treat it as 'S'
            axis = self.kss[ir][:self.nks[ir]]
            label = 'S (m)'
            if targname is not None and ir >=self.irsep-1:
                axis = np.insert(axis,0,0.0)
                axis = np.append(axis,self.smax[ir])

        return axis,label



    def plot_ring(self,ir,fig,plot_select,axis_type='S',charge_range=[]):
        # plot selected data for specified ring on figure specified

        if ir < 1 or ir > self.nrs:
            print("ERROR: Ring specified is out of range available: %d"%ir)
            return
        else:
            # ring references are in 1 .. nrs from fortran while python arrays index from zero
            ir_ref = ir-1

        plot_list = oedgeNC.all_plots[plot_select]['data']
        nplts = len(plot_list)
        
        # If nplots is 1 then the 'data' reference contains a self[] reference otherwise it contains a list of all_plot dictionary entries
        
        # set up subplots

        nrows,ncols = self.calculate_layout(nplts)

        for cnt in range(len(plot_list)):
            ax=fig.add_subplot(nrows,ncols,cnt+1)
            if len(charge_range) == 0:
                dataname,targname,label,scalef = self.get_names(plot_select,plot_list[cnt])                
                axis,x_label=self.get_axis(ir_ref,axis_type,targname)
                data=self.get_data_along_ring(ir_ref,dataname,targname,scalef=scalef)                

                # add a marker for the end points if targname is not None
                if targname is not None:
                    markers_on = [0,len(axis)-1]
                else:
                    markers_on = []
                # plot the data
                ax.plot(axis,data,linestyle='-',marker='o',markevery=markers_on,linewidth=1)
                ax.set_xlabel(x_label)
                ax.set_ylabel(label)
                ax.set_xlim(0.0,axis[-1])
            else:
                for charge in charge_range:
                    dataname,targname,label,scalef = self.get_names(plot_select,plot_list[cnt])
                    axis,x_label=self.get_axis(ir_ref,axis_type,targname)
                    data=self.get_data_along_ring(ir_ref,dataname,targname,charge,scalef=scalef)
                    # add a marker for the end points if targname is not None
                    if targname is not None:
                        markers_on = [0,len(axis)-1]
                    else:
                        markers_on = []
                    # plot the data
                    ax.plot(axis,data,label='IZ={}'.format(charge),linestyle='-',marker='o',markevery=markers_on,linewidth=1)
                    ax.set_xlabel(x_label)
                    ax.set_ylabel(label)
                    ax.set_xlim(0.0,axis[-1])
                ax.legend()



    def plot_along_ring(self,plot_select,axis_type,ring_range,charge_range=[],zoom=None):
        # Routine plots the quantities along ring ... number of plots/page is controlled by number
        # of items in the plot list.
      
        # Single figures
        # Single ring/single quantity
        # single ring/multiple quantities
        #
        # Figure notebook spanning ring range
        # multi rings/single quantity
        # multi rings/multiple quantities
        #
        # If more than one ring is specified then a figure notebook is used
        #
        # If more than one charge - try to put them on the same axes
        #
        # plot
        
        nrings = len(ring_range)
        if nrings <= 0:
            # error condition
            print("ERROR: No rings specified for the plot")
            return

        minring = 1
        maxring = self.nrs
        
        #print 'rings:',nrings,ring_range

        res = [ring_range[i]<minring or ring_range[i]>maxring for i in range(len(ring_range))]

        #print 'res:',res
        
        out_of_bounds=np.where(res==True)
        #print 'out_of_bounds:',len(out_of_bounds),
        #print 'test:',ring_range[out_of_bounds]

    
        # Check to see if at least some rings are in range 
        if len(out_of_bounds) >= nrings:
            print("ERROR: No valid rings numbers specified: min=%d max=%d range=%s"%(minring,maxring,str(ring_range)))
            return

        plot_list = oedgeNC.all_plots[plot_select]['data']
        nplts = len(plot_list)

        # if nplts is 1 then the list contains a self reference ... if greater than 1 then it contains a list of plot references
        # Do everything in a figure notebook since there can be 1+ figures generated


        #print 'rings:',nrings,ring_range,plot_select,axis_type

#        # use figure notebook since it supports an unknown number of plots
#        fig_notebook = FigureNotebook(nfig=0,name='Along Ring Plots')
#
#        for ir in ring_range:
#            fig = fig_notebook.add_figure(label='IR={}'.format(ir))
#            self.plot_ring(ir,fig,plot_select,axis_type,charge_range)
#            fig.canvas.draw()
        
        # jhnmod: use multiple figures for rings rather than OMFIT FigureNotebook
        for ir in ring_range:
            psin=self['PSIFL']['data'][ir][0]
            fig=mpl.pyplot.figure()
            self.plot_ring(ir,fig,plot_select,axis_type,charge_range)
            fig.suptitle('IR={0}  PsiN={1:0.4f}'.format(ir,psin))

