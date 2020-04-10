# Script written by Eric Emdee, given to Shawn Zamperini on 2/27/20.

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colorsMPL
import math
import os.path

from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
import scipy
from scipy import interpolate
from matplotlib import ticker
from matplotlib.patches import Ellipse

mp=1.6726219236900001e-27
mn=mp
qe=1.6021766339999999E-019
ev=qe
eps0=8.8541878128000006E-012
me=9.10938356e-31

def stupidFunction():
	print('hey stupid')

def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False

def read_field(f,fieldname,dims,intField=False):
    line = f.readline().rstrip()
    # find the right field
    while fieldname not in line:
        line = f.readline().rstrip()
        if line == -1:
            print('read_ifield: EOF reached without finding '+str(fieldname))

    # Consistency check: number of elements specified
    # in the file should equal prod(dims)

    for i in range(len(line.split())):
        if is_number(line.split()[i]): numin = int(line.split()[i])

    if (numin != np.prod(dims)):
        print(line)
        print('read_ifield: inconsistent number of input elements.');

    fieldVal=[]
    # collect field values
    while (len(fieldVal) != numin):
        line = f.readline().rstrip()
        for i in range(len(line.split())):
            if (intField): fieldVal.append(int(line.split()[i]))
            else:
                fieldVal.append(float(line.split()[i]))
    fieldVal=np.array(fieldVal)

    if (np.size(dims) > 1): fieldVal = fieldVal.reshape(dims,order='F').copy()
    return fieldVal

def read_ft44_field(fid,ver,fieldname,dims,intField=False):
# Auxiliary routine to read fields from fort.44 file
# Verion 20160829: field label and size are specified in fort.44
# Do consistency check on the data
    if (ver >= 20160829):
        # Search the file until identifier 'fieldname' is found
        line = fid.readline().rstrip()
        while fieldname not in line:
            line = fid.readline().rstrip()
            if line == -1: print('EOF reached without finding '+str(fieldname))
        # Consistency check: number of elements specified in the file should equal
        # prod(dims)
        for i in range(len(line.split())):
            if is_number(line.split()[i]): numin = int(line.split()[i])

        if (numin != np.prod(dims) and 'wld' not in fieldname):
            print('issue with field '+fieldname)
            print("numin="+str(numin))
            print("np.prod(dims)="+str(np.prod(dims)))
            print('read_ft44_rfield: inconsistent number of input elements.')

    # Read the data
    fieldVal=[]
    # collect field values
    while (len(fieldVal) != numin):
        line = fid.readline().rstrip()
        if ('wld' in fieldname) and len(fieldVal)>=numin-1: break
        for i in range(len(line.split())):
            if ('wlpump' in fieldname):
                if not is_number(line.split()[i]): continue
            if (intField): fieldVal.append(int(line.split()[i]))
            else: fieldVal.append(float(line.split()[i]))
    fieldVal=np.array(fieldVal)

    if (np.size(dims) > 1): fieldVal = fieldVal.reshape(dims,order='F').copy()

    return fieldVal


def read_b2fgmtry(b2fgmtryLoc):
    fieldname = 'nx,ny'
    fid = open(b2fgmtryLoc)
    line = fid.readline().rstrip()#read the first line
    version = line[7:17]
    print('read_b2fgmtry -- file version '+version)#check the version
    dim = read_field(fid,'nx,ny',2,True)#check the grid size
    nx  = dim[0]
    ny  = dim[1]
    qcdim = [nx+2,ny+2]
    if (version >= '03.001.000'): qcdim = [nx+2,ny+2,2]

    #initialize class that will hold all the data of b2fgmtry
    class gmtryResults:
        def __init__(self):

            # Read symmetry information

            self.isymm = read_field(fid,'isymm',1,True)

            # Read gmtry variables
            self.crx  = read_field(fid,'crx' ,[nx+2,ny+2,4])
            self.cry  = read_field(fid,'cry' ,[nx+2,ny+2,4])
            self.fpsi = read_field(fid,'fpsi',[nx+2,ny+2,4])
            self.ffbz = read_field(fid,'ffbz',[nx+2,ny+2,4])
            self.bb   = read_field(fid,'bb'  ,[nx+2,ny+2,4])
            self.vol  = read_field(fid,'vol' ,[nx+2,ny+2])
            self.hx   = read_field(fid,'hx'  ,[nx+2,ny+2])
            self.hy   = read_field(fid,'hy'  ,[nx+2,ny+2])
            self.qz   = read_field(fid,'qz'  ,[nx+2,ny+2,2])
            self.qc   = read_field(fid,'qc'  ,qcdim)
            self.gs   = read_field(fid,'gs'  ,[nx+2,ny+2,3])

            # Some other geometrical parameters

            self.nlreg       = read_field(fid,'nlreg',1,True)
            self.nlxlo       = read_field(fid,'nlxlo',self.nlreg,True)
            self.nlxhi       = read_field(fid,'nlxhi',self.nlreg,True)
            self.nlylo       = read_field(fid,'nlylo',self.nlreg,True)
            self.nlyhi       = read_field(fid,'nlyhi',self.nlreg,True)
            self.nlloc       = read_field(fid,'nlloc',self.nlreg,True)

            self.nncut       = read_field(fid,'nncut'    ,1,True)
            self.leftcut     = read_field(fid,'leftcut'  ,self.nncut,True)
            self.rightcut    = read_field(fid,'rightcut' ,self.nncut,True)
            self.topcut      = read_field(fid,'topcut'   ,self.nncut,True)
            self.bottomcut   = read_field(fid,'bottomcut',self.nncut,True)

            self.leftix      = read_field(fid,'leftix'   ,[nx+2,ny+2],True)
            self.rightix     = read_field(fid,'rightix'  ,[nx+2,ny+2],True)
            self.topix       = read_field(fid,'topix'    ,[nx+2,ny+2],True)
            self.bottomix    = read_field(fid,'bottomix' ,[nx+2,ny+2],True)
            self.leftiy      = read_field(fid,'leftiy'   ,[nx+2,ny+2],True)
            self.rightiy     = read_field(fid,'rightiy'  ,[nx+2,ny+2],True)
            self.topiy       = read_field(fid,'topiy'    ,[nx+2,ny+2],True)
            self.bottomiy    = read_field(fid,'bottomiy' ,[nx+2,ny+2],True)

            self.region      = read_field(fid,'region'     ,[nx+2,ny+2,3],True)
            self.nnreg       = read_field(fid,'nnreg'      ,3,True)
            self.resignore   = read_field(fid,'resignore'  ,[nx+2,ny+2,2],True)
            self.periodic_bc = read_field(fid,'periodic_bc',1,True)

            self.pbs  = read_field(fid,'pbs' ,[nx+2,ny+2,2])

            self.parg = read_field(fid,'parg',100)
    gmtry = gmtryResults()#instantiate class
    # Close file
    fid.close()
    print('done reading geometry file')
    return gmtry

def read_b2fstate(b2fstateLoc):
    fieldname = 'nx,ny'
    fid = open(b2fstateLoc)
    line = fid.readline().rstrip()#read the first line
    version = line[7:17]
    print('read_b2fstate -- file version '+version)#check the version
    dim = read_field(fid,'nx,ny,ns',3,True)#check the grid size
    nx  = dim[0]
    ny  = dim[1]
    ns  = dim[2]
    fluxdim  = [nx+2,ny+2,2]
    fluxdimp = [nx+2,ny+2]
    fluxdims = [nx+2,ny+2,2,ns]
    if (version >= '03.001.000'):
        fluxdim  = [nx+2,ny+2,2,2]
        fluxdimp = fluxdim
        fluxdims = [nx+2,ny+2,2,2,ns]
    #initialize class that will hold all the data of b2fstate
    class stateResults:
        def __init__(self):
            # Read charges etc.

            self.zamin = read_field(fid,'zamin',ns)
            self.zamax = read_field(fid,'zamax',ns)
            self.zn    = read_field(fid,'zn',ns)
            self.am    = read_field(fid,'am',ns)

            # Read state variables

            self.na     = read_field(fid,'na',[nx+2,ny+2,ns])
            self.ne     = read_field(fid,'ne',[nx+2,ny+2])
            self.ua     = read_field(fid,'ua',[nx+2,ny+2,ns])
            self.uadia  = read_field(fid,'uadia',[nx+2,ny+2,2,ns])
            self.te     = read_field(fid,'te',[nx+2,ny+2])
            self.ti     = read_field(fid,'ti',[nx+2,ny+2])
            self.po     = read_field(fid,'po',[nx+2,ny+2])

            # Read fluxes

            self.fna    = read_field(fid,'fna',fluxdims)
            self.fhe    = read_field(fid,'fhe',fluxdim)
            self.fhi    = read_field(fid,'fhi',fluxdim)
            self.fch    = read_field(fid,'fch',fluxdim)
            self.fch_32 = read_field(fid,'fch_32',fluxdim)
            self.fch_52 = read_field(fid,'fch_52',fluxdim)
            self.kinrgy = read_field(fid,'kinrgy',[nx+2,ny+2,ns])
            self.time   = read_field(fid,'time',1)
            self.fch_p  = read_field(fid,'fch_p',fluxdimp)

    state = stateResults()#instantiate class
    # Close file
    fid.close()
    print('done reading state file')
    return state

def read_b2fplasmf(fileName,nx,ny,ns):
#
# Read formatted b2fplasmf file created by B2.5.

    if not (os.path.isfile(fileName)):
        print("b2fplasmf: Cannot find the filename")
        return 0
    fid = open(fileName)
    if (fid == -1): print("read_b2fplasmf: can't open file")


    # Get version of the b2fstate file

    line    = fid.readline().rstrip()
    version = line[7:16]

    print('read_b2fplasmf -- file version '+version)


    # Expected array sizes, gmtry
    qcdim = [nx+2,ny+2]
    if version >= '03.001.000': qcdim  = [nx+2,ny+2,2]


    # Expected array sizes, state
    fluxdim  = [nx+2,ny+2,2]
    fluxdims = [nx+2,ny+2,2,ns]
#     if version >= '03.001.000':
#         fluxdim  = [nx+2,ny+2,2,2]
#         fluxdims = [nx+2,ny+2,2,2,ns]

    class plasmfResults:
        def __init__(self):
            # Read gmtry variables

            self.crx  = read_field(fid,'crx' ,[nx+2,ny+2,4])
            self.cry  = read_field(fid,'cry' ,[nx+2,ny+2,4])

                        # Read state variables

            self.fch    = read_field(fid,'fch'   ,fluxdim)
            self.fch0   = read_field(fid,'fch0'  ,fluxdim)
            self.fchp   = read_field(fid,'fchp'  ,fluxdim)
            self.fhe    = read_field(fid,'fhe'   ,fluxdim)
            self.fhe0   = read_field(fid,'fhe0'  ,fluxdim)
            self.fhep   = read_field(fid,'fhep'  ,fluxdim)
            self.fhet   = read_field(fid,'fhet'  ,fluxdim)
            self.fhi    = read_field(fid,'fhi'   ,fluxdim)
            self.fhi0   = read_field(fid,'fhi0'  ,fluxdim)
            self.fhip   = read_field(fid,'fhip'  ,fluxdim)
            self.fhit   = read_field(fid,'fhit'  ,fluxdim)
            self.fna    = read_field(fid,'fna'   ,fluxdims)
            self.fna0   = read_field(fid,'fna0'  ,fluxdims)
            self.fnap   = read_field(fid,'fnap'  ,fluxdims)
            self.fne    = read_field(fid,'fne'   ,fluxdim)
            self.fni    = read_field(fid,'fni'   ,fluxdim)
            self.na     = read_field(fid,'na'    ,[nx+2,ny+2,ns])
            self.na0    = read_field(fid,'na0'   ,[nx+2,ny+2,ns])
            self.nap    = read_field(fid,'nap'   ,[nx+2,ny+2,ns])
            self.ne     = read_field(fid,'ne'    ,[nx+2,ny+2])
            self.ne0    = read_field(fid,'ne0'   ,[nx+2,ny+2])
            self.ne2    = read_field(fid,'ne2'   ,[nx+2,ny+2])
            self.nep    = read_field(fid,'nep'   ,[nx+2,ny+2])
            self.ni     = read_field(fid,'ni'    ,[nx+2,ny+2,2])
            self.ni0    = read_field(fid,'ni0'   ,[nx+2,ny+2,2])
            self.pb     = read_field(fid,'pb'    ,[nx+2,ny+2])
            self.po     = read_field(fid,'po'    ,[nx+2,ny+2])
            self.po0    = read_field(fid,'po0'   ,[nx+2,ny+2])
            self.pop    = read_field(fid,'pop'   ,[nx+2,ny+2])
            self.te     = read_field(fid,'te'    ,[nx+2,ny+2])
            self.te0    = read_field(fid,'te0'   ,[nx+2,ny+2])
            self.tep    = read_field(fid,'tep'   ,[nx+2,ny+2])
            self.ti     = read_field(fid,'ti'    ,[nx+2,ny+2])
            self.ti0    = read_field(fid,'ti0'   ,[nx+2,ny+2])
            self.tip    = read_field(fid,'tip'   ,[nx+2,ny+2])
            self.ua     = read_field(fid,'ua'    ,[nx+2,ny+2,ns])
            self.ua0    = read_field(fid,'ua0'   ,[nx+2,ny+2,ns])
            self.uap    = read_field(fid,'uap'   ,[nx+2,ny+2,ns])
            self.uadia  = read_field(fid,'uadia' ,[nx+2,ny+2,2,ns])
            self.fchdia = read_field(fid,'fchdia',fluxdim)
            self.fmo    = read_field(fid,'fmo'   ,fluxdims)
            self.fna_32 = read_field(fid,'fna_32',fluxdims)
            self.fna_52 = read_field(fid,'fna_52',fluxdims)
            self.fni_32 = read_field(fid,'fni_32',fluxdim)
            self.fni_52 = read_field(fid,'fni_52',fluxdim)
            self.fne_32 = read_field(fid,'fne_32',fluxdim)
            self.fne_52 = read_field(fid,'fne_52',fluxdim)
            self.wadia  = read_field(fid,'wadia' ,[nx+2,ny+2,2,ns])
            self.vaecrb = read_field(fid,'vaecrb',[nx+2,ny+2,2,ns])
            self.facdrift     = read_field(fid,'facdrift'    ,[nx+2,ny+2])
            self.fac_ExB      = read_field(fid,'fac_ExB'     ,[nx+2,ny+2])
            self.fchvispar    = read_field(fid,'fchvispar'   ,fluxdim)
            self.fchvisper    = read_field(fid,'fchvisper'   ,fluxdim)
            self.fchin        = read_field(fid,'fchin'       ,fluxdim)
            self.fna_nodrift  = read_field(fid,'fna_nodrift' ,fluxdims)
            self.fac_vis      = read_field(fid,'fac_vis'     ,[nx+2,ny+2])
            self.fna_mdf      = read_field(fid,'fna_mdf'     ,fluxdims)
            self.fhe_mdf      = read_field(fid,'fhe_mdf'     ,fluxdim)
            self.fhi_mdf      = read_field(fid,'fhi_mdf'     ,fluxdim)
            self.fnaPSch      = read_field(fid,'fnaPSch'     ,fluxdims)
            self.fhePSch      = read_field(fid,'fhePSch'     ,fluxdim)
            self.fhiPSch      = read_field(fid,'fhiPSch'     ,fluxdim)
            self.fna_fcor     = read_field(fid,'fna_fcor'    ,fluxdims)
            self.fna_he       = read_field(fid,'fna_he'      ,fluxdims)
            self.fchvisq      = read_field(fid,'fchvisq'     ,fluxdim)
            self.fchinert     = read_field(fid,'fchinert'    ,fluxdim)
            self.resco        = read_field(fid,'resco'       ,[nx+2,ny+2,ns])
            self.reshe        = read_field(fid,'reshe'       ,[nx+2,ny+2])
            self.reshi        = read_field(fid,'reshi'       ,[nx+2,ny+2])
            self.resmo        = read_field(fid,'resmo'       ,[nx+2,ny+2,ns])
            self.resmt        = read_field(fid,'resmt'       ,[nx+2,ny+2])
            self.respo        = read_field(fid,'respo'       ,[nx+2,ny+2])
            self.sch          = read_field(fid,'sch'         ,[nx+2,ny+2,4])
            self.she          = read_field(fid,'she'         ,[nx+2,ny+2,4])
            self.shi          = read_field(fid,'shi'         ,[nx+2,ny+2,4])
            self.smo          = read_field(fid,'smo'         ,[nx+2,ny+2,4,ns])
            self.smq          = read_field(fid,'smq'         ,[nx+2,ny+2,4,ns])
            self.sna          = read_field(fid,'sna'         ,[nx+2,ny+2,2,ns])
            self.sne          = read_field(fid,'sne'         ,[nx+2,ny+2,2])
            self.rsana        = read_field(fid,'rsana'       ,[nx+2,ny+2,ns])
            self.rsahi        = read_field(fid,'rsahi'       ,[nx+2,ny+2,ns])
            self.rsamo        = read_field(fid,'rsamo'       ,[nx+2,ny+2,ns])
            self.rrana        = read_field(fid,'rrana'       ,[nx+2,ny+2,ns])
            self.rrahi        = read_field(fid,'rrahi'       ,[nx+2,ny+2,ns])
            self.rramo        = read_field(fid,'rramo'       ,[nx+2,ny+2,ns])
            self.rqahe        = read_field(fid,'rqahe'       ,[nx+2,ny+2,ns])
            self.rqrad        = read_field(fid,'rqrad'       ,[nx+2,ny+2,ns])
            self.rqbrm        = read_field(fid,'rqbrm'       ,[nx+2,ny+2,ns])
            self.rcxna        = read_field(fid,'rcxna'       ,[nx+2,ny+2,ns])
            self.rcxhi        = read_field(fid,'rcxhi'       ,[nx+2,ny+2,ns])
            self.rcxmo        = read_field(fid,'rcxmo'       ,[nx+2,ny+2,ns])
            self.b2stbr_sna   = read_field(fid,'b2stbr_sna'  ,[nx+2,ny+2,ns])
            self.b2stbr_smo   = read_field(fid,'b2stbr_smo'  ,[nx+2,ny+2,ns])
            self.b2stbr_she   = read_field(fid,'b2stbr_she'  ,[nx+2,ny+2])
            self.b2stbr_shi   = read_field(fid,'b2stbr_shi'  ,[nx+2,ny+2])
            self.b2stbr_sch   = read_field(fid,'b2stbr_sch'  ,[nx+2,ny+2])
            self.b2stbr_sne   = read_field(fid,'b2stbr_sne'  ,[nx+2,ny+2])
            self.b2stbc_sna   = read_field(fid,'b2stbc_sna'  ,[nx+2,ny+2,ns])
            self.b2stbc_smo   = read_field(fid,'b2stbc_smo'  ,[nx+2,ny+2,ns])
            self.b2stbc_she   = read_field(fid,'b2stbc_she'  ,[nx+2,ny+2])
            self.b2stbc_shi   = read_field(fid,'b2stbc_shi'  ,[nx+2,ny+2])
            self.b2stbc_sch   = read_field(fid,'b2stbc_sch'  ,[nx+2,ny+2])
            self.b2stbc_sne   = read_field(fid,'b2stbc_sne'  ,[nx+2,ny+2])
            self.b2stbm_sna   = read_field(fid,'b2stbm_sna'  ,[nx+2,ny+2,ns])
            self.b2stbm_smo   = read_field(fid,'b2stbm_smo'  ,[nx+2,ny+2,ns])
            self.b2stbm_she   = read_field(fid,'b2stbm_she'  ,[nx+2,ny+2])
            self.b2stbm_shi   = read_field(fid,'b2stbm_shi'  ,[nx+2,ny+2])
            self.b2stbm_sch   = read_field(fid,'b2stbm_sch'  ,[nx+2,ny+2])
            self.b2stbm_sne   = read_field(fid,'b2stbm_sne'  ,[nx+2,ny+2])
            self.b2sihs_divue = read_field(fid,'b2sihs_divue',[nx+2,ny+2])
            self.b2sihs_divua = read_field(fid,'b2sihs_divua',[nx+2,ny+2])
            self.b2sihs_exbe  = read_field(fid,'b2sihs_exbe' ,[nx+2,ny+2])
            self.b2sihs_exba  = read_field(fid,'b2sihs_exba' ,[nx+2,ny+2])
            self.b2sihs_visa  = read_field(fid,'b2sihs_visa' ,[nx+2,ny+2])
            self.b2sihs_joule = read_field(fid,'b2sihs_joule',[nx+2,ny+2])
            self.b2sihs_fraa  = read_field(fid,'b2sihs_fraa' ,[nx+2,ny+2])
            self.b2sihs_str   = read_field(fid,'b2sihs_str' ,[nx+2,ny+2])
            self.b2npmo_smaf  = read_field(fid,'b2npmo_smaf' ,[nx+2,ny+2,4,ns])
            self.b2npmo_smag  = read_field(fid,'b2npmo_smag' ,[nx+2,ny+2,4,ns])
            self.b2npmo_smav  = read_field(fid,'b2npmo_smav' ,[nx+2,ny+2,4,ns])
            self.smpr         = read_field(fid,'smpr'        ,[nx+2,ny+2,ns])
            self.smpt         = read_field(fid,'smpt'        ,[nx+2,ny+2,ns])
            self.smfr         = read_field(fid,'smfr'        ,[nx+2,ny+2,ns])
            self.smcf         = read_field(fid,'smcf'        ,[nx+2,ny+2,ns])
            self.ext_sna      = read_field(fid,'ext_sna'     ,[nx+2,ny+2,ns])
            self.ext_smo      = read_field(fid,'ext_smo'     ,[nx+2,ny+2,ns])
            self.ext_she      = read_field(fid,'ext_she'     ,[nx+2,ny+2])
            self.ext_shi      = read_field(fid,'ext_shi'     ,[nx+2,ny+2])
            self.ext_sch      = read_field(fid,'ext_sch'     ,[nx+2,ny+2])
            self.ext_sne      = read_field(fid,'ext_sne'     ,[nx+2,ny+2])
            self.calf         = read_field(fid,'calf'        ,fluxdim)
            self.cdna         = read_field(fid,'cdna'        ,fluxdims)
            self.cdpa         = read_field(fid,'cdpa'        ,fluxdims)
            self.ceqp         = read_field(fid,'ceqp'        ,[nx+2,ny+2])
            self.chce         = read_field(fid,'chce'        ,fluxdim)
            self.chci         = read_field(fid,'chci'        ,fluxdim)
            self.chve         = read_field(fid,'chve'        ,fluxdim)
            self.chvemx       = read_field(fid,'chvemx'      ,[nx+2,ny+2])
            self.chvi         = read_field(fid,'chvi'        ,fluxdim)
            self.chvimx       = read_field(fid,'chvimx'      ,[nx+2,ny+2])
            self.csig         = read_field(fid,'csig'        ,fluxdim)
            self.cvla         = read_field(fid,'cvla'        ,fluxdims)
            self.cvsa         = read_field(fid,'cvsa'        ,fluxdims)
            self.cthe         = read_field(fid,'cthe'        ,[nx+2,ny+2,ns])
            self.cthi         = read_field(fid,'cthi'        ,[nx+2,ny+2,ns])
            self.csigin       = read_field(fid,'csigin'      ,[fluxdims[0],fluxdims[1],fluxdims[2],fluxdims[3],ns])
            self.cvsa_cl      = read_field(fid,'cvsa_cl'     ,fluxdims)
            self.fllime       = read_field(fid,'fllime'      ,[nx+2,ny+2])
            self.fllimi       = read_field(fid,'fllimi'      ,[nx+2,ny+2])
            self.fllim0fna    = read_field(fid,'fllim0fna'   ,fluxdims)
            self.fllim0fhi    = read_field(fid,'fllim0fhi'   ,fluxdims)
            self.fllimvisc    = read_field(fid,'fllimvisc'   ,[nx+2,ny+2,ns])
            self.sig0         = read_field(fid,'sig0'        ,[nx+2,ny+2])
            self.hce0         = read_field(fid,'hce0'        ,[nx+2,ny+2])
            self.alf0         = read_field(fid,'alf0'        ,[nx+2,ny+2])
            self.hci0         = read_field(fid,'hci0'        ,[nx+2,ny+2])
            self.hcib         = read_field(fid,'hcib'        ,[nx+2,ny+2,ns])
            self.dpa0         = read_field(fid,'dpa0'        ,[nx+2,ny+2,ns])
            self.dna0         = read_field(fid,'dna0'        ,[nx+2,ny+2,ns])
            self.vsa0         = read_field(fid,'vsa0'        ,[nx+2,ny+2,ns])
            self.vla0         = read_field(fid,'vla0'        ,[nx+2,ny+2,2,ns])
            self.csig_an      = read_field(fid,'csig_an'     ,fluxdim)
            self.calf_an      = read_field(fid,'calf_an'     ,fluxdim)
            nstra              = read_field(fid,'nstra'       ,[1],True)
            self.sclstra      = read_field(fid,'sclstra'     ,[ns+1,nstra[0]])
            self.sclrtio      = read_field(fid,'sclrtio'     ,[ns+1,nstra[0]])
    plasmf = plasmfResults()
    fid.close()
    print('done reading b2fplasmf')
    return plasmf


def read_ft44(fileName):
#
# Read fort.44 file
#
# For now
# - only fort.44 version 20081111 recognized
# - assuming nfla = 1 until a better fix
# - assuming nlwrmsh = 1 until a better fix
#

    print('read_ft44: assuming nlwrmsh = 1, nfla = 1.')
    nlwrmsh = 1
    nfla = 1

    fid = open(fileName)
    if (fid == -1): print("read_ft44: can't open file")

# Read dimensions

# nx, ny, version
    dims = fid.readline().rstrip().split()
    nx = int(dims[0])
    ny = int(dims[1])
    ver = int(dims[2])

    if (ver != 20081111 and ver != 20160829 and ver != 20170328):
        print('read_ft44: unknown format of fort.44 file')

# go to new line (skip reading a possible git-hash)
#     fid.readline().rstrip()

# natm, nmol, nion
    dims = fid.readline().rstrip().split()
    natm = int(dims[0])
    nmol = int(dims[1])
    nion = int(dims[2])
# for now, ignore reading species labels

    for i in range(natm): line = fid.readline().rstrip()
    for i in range(nmol): line = fid.readline().rstrip()
    for i in range(nion): line = fid.readline().rstrip()

# Read basic data, there is more, I might grab it if I find out I need it
    class ft44Results:
        def __init__(self):
            self.dab2     = read_ft44_field(fid,ver,'dab2',[nx,ny,natm]);
            self.tab2     = read_ft44_field(fid,ver,'tab2',[nx,ny,natm]);
            self.dmb2     = read_ft44_field(fid,ver,'dmb2',[nx,ny,nmol]);
            self.tmb2     = read_ft44_field(fid,ver,'tmb2',[nx,ny,nmol]);
            self.dib2     = read_ft44_field(fid,ver,'dib2',[nx,ny,nion]);
            self.tib2     = read_ft44_field(fid,ver,'tib2',[nx,ny,nion]);
            self.rfluxa   = read_ft44_field(fid,ver,'rfluxa',[nx,ny,natm]);
            self.rfluxm   = read_ft44_field(fid,ver,'rfluxm',[nx,ny,nmol]);
            self.pfluxa   = read_ft44_field(fid,ver,'pfluxa',[nx,ny,natm]);
            self.pfluxm   = read_ft44_field(fid,ver,'pfluxm',[nx,ny,nmol]);
            self.refluxa  = read_ft44_field(fid,ver,'refluxa',[nx,ny,natm]);
            self.refluxm  = read_ft44_field(fid,ver,'refluxm',[nx,ny,nmol]);
            self.pefluxa  = read_ft44_field(fid,ver,'pefluxa',[nx,ny,natm]);
            self.pefluxm  = read_ft44_field(fid,ver,'pefluxm',[nx,ny,nmol]);
            self.emiss    = read_ft44_field(fid,ver,'emiss',[nx,ny,1]);
            self.emissmol = read_ft44_field(fid,ver,'emissmol',[nx,ny,1]);
            self.srcml    = read_ft44_field(fid,ver,'srcml',[nx,ny,nmol]);
            self.edissml  = read_ft44_field(fid,ver,'edissml',[nx,ny,nmol]);
            self.wldnek   = read_ft44_field(fid,ver,'wldnek(0)',[88])
            self.wldspt0   = read_ft44_field(fid,ver,'wldspt(0)',[88])
            self.wldspt1   = read_ft44_field(fid,ver,'wldspt(  1)',[88])
            self.wldspt2   = read_ft44_field(fid,ver,'wldspt(  2)',[88])
            self.wldspt3   = read_ft44_field(fid,ver,'wldspt(  3)',[88])
            self.wldspt4   = read_ft44_field(fid,ver,'wldspt(  4)',[88])
            self.wldspt5   = read_ft44_field(fid,ver,'wldspt(  5)',[88])
            self.wldspt6   = read_ft44_field(fid,ver,'wldspt(  6)',[88])
 #           self.wlpumpA   = read_ft44_field(fid,ver,'wlpump(A)',[natm,88])
 #           self.wlpumpM   = read_ft44_field(fid,ver,'wlpump(M)',[nmol,88])
            self.eneutrad  = read_ft44_field(fid,ver,'eneutrad',[nx,ny,natm])
            self.emolrad  = read_ft44_field(fid,ver,'emolrad',[nx,ny,nmol])
            self.eionrad  = read_ft44_field(fid,ver,'eionrad',[nx,ny,nion])
    ft44 = ft44Results()
    fid.close()
    print('done reading ft44 file')
    return ft44

def read_ft46(fileName):
    # ft46 = read_ft46(file)
    #
    # Read fort.46 file. Convert to SI units.
    #
    # For now, only fort.46 version 20160513 recognized
    #

    fid = open(fileName)
    if (fid == -1): print("read_ft44: can't open file")

    # Read dimensions

    # ntri, version, avoid reading git-hash
    line = fid.readline().rstrip().split()
    ntri = int(line[0])
    ver  = int(line[1])

    if ver != 20160513 and ver != 20160829 and ver != 20170930:
        print('read_ft46: unknown format of fort.46 file')

    # natm, nmol, nion
    dims = fid.readline().rstrip().split()
    natm = int(dims[0])
    nmol = int(dims[1])
    nion = int(dims[2])

    # for now, ignore reading species labels
    for i in range(natm): fid.readline().rstrip()
    for i in range(nmol): fid.readline().rstrip()
    for i in range(nion): fid.readline().rstrip()

    eV   = 1.6021765650000000E-019
    # Read data
    class ft46Results:
        def __init__(self):
            self.pdena  = read_ft44_field(fid,ver,'pdena',[ntri,natm])*1e6# m^{-3}
            self.pdenm  = read_ft44_field(fid,ver,'pdenm',[ntri,nmol])*1e6
            self.pdeni  = read_ft44_field(fid,ver,'pdeni',[ntri,nion])*1e6

            self.edena  = read_ft44_field(fid,ver,'edena',[ntri,natm])*1e6*eV# J m^{-3}
            self.edenm  = read_ft44_field(fid,ver,'edenm',[ntri,nmol])*1e6*eV
            self.edeni  = read_ft44_field(fid,ver,'edeni',[ntri,nion])*1e6*eV

            self.vxdena = read_ft44_field(fid,ver,'vxdena',[ntri,natm])*1e1# kg s^{-1} m^{-2}
            self.vxdenm = read_ft44_field(fid,ver,'vxdenm',[ntri,nmol])*1e1
            self.vxdeni = read_ft44_field(fid,ver,'vxdeni',[ntri,nion])*1e1

            self.vydena = read_ft44_field(fid,ver,'vydena',[ntri,natm])*1e1# kg s^{-1} m^{-2}
            self.vydenm = read_ft44_field(fid,ver,'vydenm',[ntri,nmol])*1e1
            self.vydeni = read_ft44_field(fid,ver,'vydeni',[ntri,nion])*1e1

            self.vzdena = read_ft44_field(fid,ver,'vzdena',[ntri,natm])*1e1# kg s^{-1} m^{-2}
            self.vzdenm = read_ft44_field(fid,ver,'vzdenm',[ntri,nmol])*1e1
            self.vzdeni = read_ft44_field(fid,ver,'vzdeni',[ntri,nion])*1e1
            self.vol =    read_ft44_field(fid,ver,'volumes',[ntri])
            self.pux    = read_ft44_field(fid,ver,'pux',[ntri])
            self.puy    = read_ft44_field(fid,ver,'puy',[ntri])

    ft46 = ft46Results()
    # Close file
    fid.close()
    print('done reading ft46 file')
    return ft46

def read_ft33(fileName):
    #
    # Read fort.33-files (triangle nodes). Converts to SI units (m).
    #
    #

    fid = open(fileName);
    if (fid == -1): print("can't open fort.33 file")

    print('read_ft33: assuming ntrfrm = 0.')
    ntrfrm = 0


    # Read coordinates

    # number of nodes
    nnodes = int(fid.readline().rstrip())
    nodes  = [[],[]]

    if (ntrfrm==0):
            for line in fid:
                for j in range(len(line.split())):
                    nodes[0].append(float(line.split()[j]))
                if len(nodes[0])>=nnodes: break
            for line in fid:
                for j in range(len(line.split())):
                    nodes[1].append(float(line.split()[j]))
                if len(nodes[1])>=nnodes: break

    else: print('read_ft33: wrong ntrfrm.')

    # Convert from cm to m
    nodes = np.array(nodes)*1e-2
    # close file
    fid.close()
    return nodes

def read_ft34(fileName):
    # cells = read_ft34(file)
    #
    # Read fort.34-files (nodes composing each triangle).
    #

    fid = open(fileName)
    if (fid == -1): print("can't open fort.34 file")



    # Read data

    # number of triangels
    ntria = int(fid.readline().rstrip())

    cells = [[],[],[]]

    for i in range(ntria):
        line = fid.readline().rstrip().split()
        cells[0].append(int(line[1]))
        cells[1].append(int(line[2]))
        cells[2].append(int(line[3]))

    # close file
    fid.close()
    return cells

def read_ft35(fileName):
    # links = read_ft35(file)
    #
    # Read fort.35-files (triangle data).
    #

    fid = open(fileName)
    if (fid == -1): print("can't open fort.34 file")

    # Read data

    # number of triangles
    ntria = int(fid.readline().rstrip())

    class ft35Results():
        def __init__(self):
            self.nghbr = np.zeros([ntria,3]);
            self.side  = np.zeros([ntria,3]);
            self.cont  = np.zeros([ntria,3]);
            self.ixiy  = np.zeros([ntria,2]);

            for i in range (ntria):
                data = fid.readline().rstrip().split()
                data = [int(i) for i in data]
                self.nghbr[i,:] = data[1::3][0:3]
                self.side[i,:]  = data[2::3][0:3]
                self.cont[i,:]  = data[3::3]
                self.ixiy[i,:]  = data[10:12]
    links=ft35Results()
    # close file
    fid.close()
    return links

def read_triangle_mesh(fort33fn,fort34fn,fort35fn):
    # triangles = read_triangle_mesh(fort33,fort34,fort35)
    #
    # Wrapper routine to read all triangle data at once.
    #
    # Returns nodes, cells, nghbr, side and cont as fiels of triangles-struct.
    #

    class triangleResults:
        def __init__(self):
            self.nodes = np.array(read_ft33(fort33fn))#list of
            self.cells = np.array(read_ft34(fort34fn))
            centroidsX  = []
            centroidsY  = []
            nodeXs      = []
            nodeYs      = []
            for i in range(np.shape(self.cells)[1]):#loop through every triangle
                cntrX=0
                cntrY=0
                nodeX=[]
                nodeY=[]
                for j in range(3):#loop through each node on each triangle
                    cntrX=cntrX+self.nodes[:,self.cells[:,i][j]-1][0]
                    nodeX.append(self.nodes[:,self.cells[:,i][j]-1][0])
                    cntrY=cntrY+self.nodes[:,self.cells[:,i][j]-1][1]
                    nodeY.append(self.nodes[:,self.cells[:,i][j]-1][1])
                cntrX=cntrX/3#calculate centroid of triangle
                cntrY=cntrY/3
                nodeXs.append(nodeX)
                nodeYs.append(nodeY)
                centroidsX.append(cntrX)#make list of triangle centroid x-coordinate
                centroidsY.append(cntrY)#make list of triangle centroid y-coordinate
            self.nodeXs = np.array(nodeXs)
            self.nodeYs = np.array(nodeYs)
            self.triaX = np.array(centroidsX)
            self.triaY = np.array(centroidsY)
            links      = read_ft35(fort35fn)

            self.nghbr = links.nghbr
            self.side  = links.side
            self.cont  = links.cont
            self.ixiy  = links.ixiy
    triangles=triangleResults()
    return triangles

def readB2Plot(fileLoc):
    fid = open(fileLoc)
    title = fid.readline().rstrip()
    line  = fid.readline().rstrip().split()
    dataList =[[],[],[]]
    while (is_number(line[0])):
        dataList[1].append(float(line[0]))
        dataList[2].append(float(line[1]))
        line = fid.readline().rstrip().split()
    line = fid.readline().rstrip().split()
    while (is_number(line[0])):
        dataList[0].append(float(line[0]))
        line = fid.readline().rstrip().split()
        if not line: break
    fid.close()
    return np.array(dataList)

def read_tally_field(fid,fieldname):
    line = fid.readline().rstrip()
    while fieldname not in line: line = fid.readline().rstrip()
    line = fid.readline().rstrip()
    data=[]
    while line:
        if is_number(line.split()[0]): data.append(np.array(line.split()).astype(np.float))
        else: data.append(np.array(line.split()[1:]).astype(np.float))
        line = fid.readline().rstrip()
    if np.shape(data)[0]==1: data=data[0]
    return np.array(data)

#reads display.tallies which is the result of display_tallies > display.tallies
def readTallyDisplay(fileLoc,timestep=-1):
    fid = open(fileLoc)
    line = fid.readline().rstrip()
    while 'ITER' not in line: line = fid.readline().rstrip()
 #   maxIter=line.split()[5]
 #   if (timestep==-1): timestep=maxIter
 #   while timestep!=line.split()[1]: line = fid.readline().rstrip()
    class tallyResults:
        def __init__(self):
            self.rsanareg       = read_tally_field(fid,'rsanareg')
            self.rsahireg       = read_tally_field(fid,'rsahireg')
            self.rsamoreg       = read_tally_field(fid,'rsamoreg')
            self.rranareg       = read_tally_field(fid,'rranareg')
            self.rrahireg       = read_tally_field(fid,'rrahireg')
            self.rramoreg       = read_tally_field(fid,'rramoreg')
            self.rqahereg       = read_tally_field(fid,'rqahereg')
#            self.rqradreg       = read_tally_field(fid,'rqradreg')
#            self.rqbrmreg       = read_tally_field(fid,'rqbrmreg')
            self.rcxnareg       = read_tally_field(fid,'rcxnareg')
            self.rcxhireg       = read_tally_field(fid,'rcxhireg')
            self.rcxmoreg       = read_tally_field(fid,'rcxmoreg')
            self.fnaxreg        = read_tally_field(fid,'fnaxreg')
            self.fnayreg        = read_tally_field(fid,'fnayreg')
            self.fhixreg        = read_tally_field(fid,'fhixreg')
            self.fhiyreg        = read_tally_field(fid,'fhiyreg')
            self.fhexreg        = read_tally_field(fid,'fhexreg')
            self.fheyreg        = read_tally_field(fid,'fheyreg')
            self.fhpxreg        = read_tally_field(fid,'fhpxreg')
            self.fhpyreg        = read_tally_field(fid,'fhpyreg')
            self.fhmxreg        = read_tally_field(fid,'fhmxreg')
            self.fhmyreg        = read_tally_field(fid,'fhmyreg')
            self.fchxreg        = read_tally_field(fid,'fchxreg')
            self.fchyreg        = read_tally_field(fid,'fchyreg')
            self.fhtxreg        = read_tally_field(fid,'fhtxreg')
            self.fhtyreg        = read_tally_field(fid,'fhtyreg')
            self.fhjxreg        = read_tally_field(fid,'fhjxreg')
            self.fhjyreg        = read_tally_field(fid,'fhjyreg')
#            self.qconvixreg     = read_tally_field(fid,'qconvixreg')
#            self.qconviyreg     = read_tally_field(fid,'qconviyreg')
#            self.qconvexreg     = read_tally_field(fid,'qconvexreg')
#            self.qconveyreg     = read_tally_field(fid,'qconveyreg')
            self.b2stbr_sna_reg = read_tally_field(fid,'b2stbr_sna_reg')
#            self.b2stbr_sne_reg = read_tally_field(fid,'b2stbr_sne_reg')
            self.b2stbr_she_reg = read_tally_field(fid,'b2stbr_she_reg')
            self.b2stbr_shi_reg = read_tally_field(fid,'b2stbr_shi_reg')
            self.b2stbr_sch_reg = read_tally_field(fid,'b2stbr_sch_reg')
            self.b2stbc_sna_reg = read_tally_field(fid,'b2stbc_sna_reg')
            self.b2stbc_she_reg = read_tally_field(fid,'b2stbc_she_reg')
            self.b2stbc_shi_reg = read_tally_field(fid,'b2stbc_shi_reg')
            self.b2stbm_she_reg = read_tally_field(fid,'b2stbm_she_reg')
            self.b2stbm_shi_reg = read_tally_field(fid,'b2stbm_shi_reg')
            self.nareg          = read_tally_field(fid,'nareg')
            self.tereg          = read_tally_field(fid,'tereg')
            self.nereg          = read_tally_field(fid,'nereg')
            self.ne2reg         = read_tally_field(fid,'ne2reg')
            self.tireg          = read_tally_field(fid,'tireg')
            self.nireg          = read_tally_field(fid,'nireg')
            self.poreg          = read_tally_field(fid,'poreg')
            self.volreg         = read_tally_field(fid,'volreg')
            self.b2brem         = read_tally_field(fid,'b2brem')
            self.b2rad          = read_tally_field(fid,'b2rad')
            self.b2qie          = read_tally_field(fid,'b2qie')
            self.b2vdp          = read_tally_field(fid,'b2vdp')
            self.b2divue        = read_tally_field(fid,'b2divue')
            self.b2divua        = read_tally_field(fid,'b2divua')
            self.b2exbe         = read_tally_field(fid,'b2exbe')
            self.b2exba         = read_tally_field(fid,'b2exba')
            self.b2visa         = read_tally_field(fid,'b2visa')
            self.b2joule        = read_tally_field(fid,'b2joule')
            self.b2fraa         = read_tally_field(fid,'b2fraa')
            self.b2she          = read_tally_field(fid,'b2she')
            self.b2shi          = read_tally_field(fid,'b2shi')
            self.b2she0         = read_tally_field(fid,'b2she0')
            self.b2shi0         = read_tally_field(fid,'b2shi0')
    tally = tallyResults()
    print('finished reading tallies')
    return tally

def plotvar(xPts, yPts, var,minColor='none',maxColor='none', cbScale='linear',cbTitle=r'Density $m^{-3}$',title='SOLPS data',
            xlims=[1.2,2.4],ylims=[-1.15,0.75],colorBarOn=True,filename='NONE'):

    patches = []
    nx = np.shape(xPts)[0]
    ny = np.shape(xPts)[1]
    for iy in np.arange(0,ny):
        for ix in np.arange(0,nx):
            rcol = xPts[ix,iy,[0,1,3,2]]
            zcol = yPts[ix,iy,[0,1,3,2]]
            rcol.shape=(4,1)
            zcol.shape=(4,1)
            polygon = Polygon(np.column_stack((rcol,zcol)), True,linewidth=0)
            patches.append(polygon)

    vals=var.T.flatten()
    if (cbScale=='symlog'):
        p = PatchCollection(patches,True,cmap='bwr')
    else:
        p = PatchCollection(patches,True,cmap='viridis')
    p.set_array(np.array(vals))
    if (minColor!='none'):
        if (cbScale=='linear'):
#            p.vmin=10
#            p.vmin=20
            p.set_clim([minColor,maxColor])
        if (cbScale=='log'):
            p.norm=colorsMPL.LogNorm(vmin=minColor,vmax=maxColor)
        if (cbScale=='symlog'):
#            p.cmap='bwr'
            p.norm=colorsMPL.SymLogNorm(linthresh=maxColor/1800,linscale=0.0003,vmin=minColor,vmax=maxColor)



    fig,axs = plt.subplots(1,figsize=(6, 12))

    axs.add_collection(p)
    if (colorBarOn):
        cb = plt.colorbar(p,ax=axs,pad=0.01)
        cb.ax.tick_params(labelsize=20)
        cb.set_label(cbTitle,fontsize=25)

    axs.set_title(title,fontsize=20)
    axs.set_ylim(ylims)
    axs.set_xlim(xlims)
    axs.tick_params(axis='both',labelsize=15)
#     ax.set_ylim([-1.15,-0.5])
#     ax.set_xlim([1.3,1.8])
    plt.xlabel('R [m]',fontsize=20)
    plt.ylabel('Z [m]',fontsize=20)
    plt.grid(True)

    if filename != 'NONE':
        plt.savefig(filename)
    plt.show()


def read_b2wdat(filename):
    f = open(filename)
    line = f.readline().rstrip().split()
    fieldVal = []
    while (line!=[]):
        line = f.readline().rstrip().split()
        if line==[]: break
        fieldVal.append([float(i) for i in line][1:])
    return np.array(fieldVal[::-1]).T

def b2tlnl(nx, ny, te, ti, ne,icase=0):

    ev=1.6021766339999999E-019
#-----------------------------------------------------------------------
#
# purpose
#
#     B2TLNL computes the Coulomb logarithm according to Braginskii or
#     Wesson formulas.
#
#
#     lnlam is the log of the plasma parameter.
#
#-----------------------------------------------------------------------
#declarations
    lnlam=np.zeros(np.shape(te))

#.computation
    lamda=-5.0
    if (lamda<0):
        if icase==0:#Braginskii
            for iy in range(ny):
                for ix in range(nx):
                    if(te[ix][iy]/ev <= 50.0):
                        lnlam[ix][iy]=max(-lamda,23.4 - 1.15*math.log10(ne[ix][iy]/1.0e6) +
                                      3.45*math.log10(te[ix][iy]/ev))
                    else:
                        lnlam[ix][iy]=max(-lamda,25.3 - 1.15*math.log10(ne[ix][iy]/1.0e6) +
                                      2.30*math.log10(te[ix][iy]/ev))
    else:
        lnlam = lamda
    return lnlam

def is_neutral(a):
    if a==0 or a==2:
        return True
    else:
        return False

def fce1(z): return ((1.0+0.24*z)*(1.0+0.93*z))/((1.0+2.56*z)*(1.0+0.29*z))

def fce2(z): return ((1.0+1.40*z)*(1.0+0.52*z))/((1.0+2.56*z)*(1.0+0.29*z))*1.56

def fce2n(z): return fce2(z)/(z+math.sqrt(2.0)/2.0)

def fal_cen(z): return -fce2n(z)/fce1(z)

def zmffCalc(zamax,na,ns,ismain):
    zmff   = np.zeros(np.shape(na[:,:,0]))
    for sI in range(ns):
        if(sI!=ismain): zmff = zmff + zamax[sI]**2 * na[:,:,sI]
    zmff=zmff/(zamax[ismain]**2 * na[:,:,ismain])
    return zmff

def fkabvp(a,b,zamax,na):
    ismain = 1
    ns=len(zamax)
    zmff   = zmffCalc(zamax,na,ns,ismain)
    cimp1=fce1(zmff)
    if (a==ismain and (b!=a) and not is_neutral(a) and not is_neutral(b)):
        fkabvp=cimp1
    elif ((b==ismain) and (a!=b) and not is_neutral(a) and not is_neutral(b)):
        fkabvp=cimp1
    elif((a!=b and a!=ismain) and (b!=ismain) and not is_neutral(a) and not is_neutral(b)):
        fkabvp=np.ones(np.shape(na[:,:,0]))
    elif((a==b) and (not is_neutral(a)) and (not is_neutral(b))):
        fkabvp=np.zeros(np.shape(na[:,:,0]))
    else:
        fkabvp=np.zeros(np.shape(na[:,:,0]))
    return fkabvp

def fkabtf(a, b,zamax,na):
    ismain=1
    ns=len(zamax)
    zmff = zmffCalc(zamax,na,ns,ismain)
    cimp2=fce2(zmff)
    if ((b==ismain) and (b!=a) and (not is_neutral(a)) and (not is_neutral(b))):
        fkabtf=cimp2
    elif ((a!=ismain) and (b!=ismain) and (not is_neutral(a)) and (not is_neutral(b))):
        fkabtf=0.0
    elif ((a==ismain)):
        fkabtf=0.0
    else:
        fkabtf=0.0
    return fkabtf

def fka(a,zamax,na,am):
    ns = len(zamax)
    rz2 = zamax**2
    fka = np.zeros(np.shape(na[:,:,0]))
    for r in range(ns):
        fka = fka + rz2[r]*na[:,:,r]*math.sqrt(mp)*math.sqrt(am[a]*am[r]/(am[a]+am[r]))
    fka = fka*rz2[a]
    return fka


def b2xpne(ns, rza, na):# b2aux/b2xpne.F
#     ------------------------------------------------------------------
#     B2XPNE computes the electron density, ne:
#       ne(,) = (sum is :: rza(,,is)*na(,,is))
#     I'm using it to calculate ne2 which is used in some functions below
#     ------------------------------------------------------------------
    ne = np.zeros(np.shape(na[:,:,0]))
    for species in range(ns):
        ne = ne + rza[species]*na[:,:,species]
    return ne
