
"""
Read and interpolate contents of ADAS adf11
Python 2.7
JH Nichols 3/12/18
"""

import numpy as np
import scipy.interpolate as spi
import matplotlib.pyplot as plt

#print('Remember that ADAS densities have unit 1/cm^3')

  ##------------------------------------------------------------------------
  # Read adf11 (iso-nuclear master files) from ADAS.
  # For now: Metastable unresolved files only.
  #
  # Coefficients include acd, scd, ccd, plt, plt, prb, prc
  #
  # input:
  #  filename  - Path and filename
  #  debug     - Triggers optional terminal output
  #  degree    - Optional higher order spline interpolation
  #  evaldens  - Density [cm^-3] at which to evaluate coefficient
  #  evaltemp  - Temperature [eV] at which to evaluate coefficient
  #  evaliz    - Charge state at which to evaluate coefficient

  # usage:
  #  a=adf11(filename,(debug=False))
  #  a.interpolate(,(degree=1))
  #  a.EvalInterpolate(evaldens,evaltemp,evaliz)
  #     ...do whatever with the result
  #
  # plots
  #  a.plot_adf11_contour(evaliz)
  ##------------------------------------------------------------------------
class adf11():

#most interesting:
#species, dens, temp, charge, coefficent[dens,temp,charge]

    def __init__(self,filename,debug=False):
        if debug:
            print('Read ' + filename)
        fp = open(filename,'r')
        header=[]
        header = fp.readline().split()
        if debug:
            print("header:" + str(header))

        self.izmax=int(header[0])   # nuclear charge
        self.n_dens=int(header[1])  # number of densities
        self.n_temp=int(header[2])  # number of temperatures
        self.iz1min=int(header[3])  # lowest ion charge +1
        self.iz1max=int(header[4])  # highest ion charge
        self.element=header[5][1:]
        self.source='-'.join(header[6:])[1:]

        if debug:
            print("element:",self.element)
            print("source:",self.source)

        dens=np.zeros(self.n_dens)
        temp=np.zeros(self.n_temp)
        data=np.zeros((self.n_dens,self.n_temp,self.izmax))

        fp.readline() #skip ---------

        #Read density vector
        line=[]
        for j in range(int(self.n_dens/8)):
            line.extend(fp.readline().split())
        if (self.n_dens%8!=0):
            line.extend(fp.readline().split())
        dens[0:len(line)]=[float(x) for x in line]
        if debug:
            print("log densities:" + str(dens))

        #Read temperature vector
        line=[]
        for j in range(int(self.n_temp/8)):
            line.extend(fp.readline().split())
        if (self.n_temp%8!=0):
            line.extend(fp.readline().split())
        temp[0:len(line)]=[float(x) for x in line]
        if debug:
            print("log temperatures:" + str(temp))

        #Read data array for each charge state
        igrd=[] # recombined ion metastable index
        iprt=[] # parent metastable index
        z1=[]
        date=[]
        for z in range(int(self.izmax)):
            header=fp.readline().split()
            if debug:
                print("header" + str(z) + ":" + str(header))
            for s in enumerate(header):
                if s[1].find('IGRD=')>=0:
                    igrd.append(int(header[s[0]+1]))
                if s[1].find('IPRT=')>=0:
                    iprt.append(int(header[s[0]+1]))
                if s[1].find('Z1=')>=0:
                    #z1.append(int(header[s[0]+1]))
                    pass
                if s[1].find('DATE=')>=0:
                    date.append(header[s[0]+1])
            for j in range(int(self.n_temp)):
                line=[]
                for k in range(int(self.n_dens/8)):
                    line.extend(fp.readline().split())
                if (self.n_dens%8!=0):
                    line.extend(fp.readline().split())
                #print line
                data[0:len(line),j,z]=[float(x) for x in line]

        self.comment=fp.read()

        fp.close()

        self.dens=np.asarray(dens.T)
        self.temp=np.asarray(temp.T)
        self.data=np.asarray(data)
        self.igrd=np.asarray(igrd)
        self.iprt=np.asarray(iprt)
        self.z1=np.asarray(z1)
        self.date=np.asarray(date)

        if debug:
            print("recombined metastable indices:" + str(self.igrd))
            print("parent metastable indices:" + str(self.iprt))
            print("ion charges:" + str(self.z1))
            print("dates:" + str(self.date))
            print(str(self.comment))

    def interpolate(self,degree=1):
        #create an interpolating spline to evaluate at arbitrary points
            #use first order as standard as this can cause no problems
            #better (higher order) interpolation possible by using 'degree'

        #It is recommended to use logarithmic values for interpolation
            #adf11 values are already written as log10
        logtemp=self.temp
        logdens=self.dens
        logdata=self.data

        #do actual interpolation
        self.logInterpolate=[]
        for i in range(self.izmax):
            self.logInterpolate.append(spi.RectBivariateSpline(logdens, \
                logtemp, logdata[:,:,i], kx=degree, ky=degree))
        print("Ionization coefficient interpolation done.")

    def EvalInterpolation(self,evaldens,evaltemp,evaliz):
        #evaluate coefficient for arbitrary density/temperature
        #must execute interpolate() method before calling EvalInterpolation

        #density in cm^-3, temperature in eV
        #evaliz=0 for ground state

        evaldens=np.float(evaldens)
        evaltemp=np.float(evaltemp)
        try:
            log10result=self.logInterpolate[evaliz](np.log10(evaldens),\
                                           np.log10(evaltemp))
        except AttributeError:
            print("the interpolated data cannot be found.")
            print("execute interpolate() before calling EvalInterpolation")
            return 0
        return np.power(10,log10result)[0]

    def plot_adf11_contour(self,evaliz=0):
        #plot density and temperature dependence of coefficient

        xlist = np.logspace(8, 15, 100)
        ylist = np.logspace(0, 4, 100)
        X, Y = np.meshgrid(xlist, ylist)
        vfunc=np.vectorize(self.EvalInterpolation)
        Z=vfunc(X,Y,evaliz)
        plt.figure()
        levels = [-15,-14.5,-14,-13.5,-13,-12.5,-12,-11.5,-11,-10.5,-10,-9.5,\
                  -9,-8.5,-8,-7.5,-7.0,-6.5,-6,-5.5,-5,-4.5,-4]
        contour = plt.contour(np.log10(X), np.log10(Y), \
                              np.log10(Z),levels,colors='k')
        plt.xlabel("Log10(Ne) [cm-3]")
        plt.ylabel("Log10(Te) [eV]")
        plt.title("your title here")
        plt.clabel(contour, colors = 'k', fmt = '%2.1f', fontsize=12)
        contour_filled = plt.contourf(np.log10(X), np.log10(Y), np.log10(Z),\
                                      100,cmap=plt.cm.jet)
        plt.colorbar(contour_filled)
        plt.show()
