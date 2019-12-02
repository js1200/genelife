""" genelife_c_interface.py
    profect genelife

    Wrapping a C library function that does update of long unsigned int arrays gol, golg
    input using the numpy.ctypeslib.
    See https://docs.python.org/2/library/ctypes.html
    See https://docs.scipy.org/doc/numpy-1.15.0/reference/generated/numpy.dtype.html#numpy.dtype  for numpy.dtype specs
    Method gleaned from 
    http://www.scipy-lectures.org/advanced/interfacing_with_c/interfacing_with_c.html
    
------------------------------------------------------------ copyright ----------------------------------------------------------------------------------
  Written by John S. McCaskill and Norman H. Packard 2017-2019
  First created by John McCaskill on 14.07.2017. Last modified Nov 2019.

  MIT License
    Copyright (c) 2017,2018,2019 John S. McCaskill and Norman H. Packard

    Permission is hereby granted, free of charge, to any person obtaining a copy of
    this software and associated documentation files (the "Software"), to deal in
    the Software without restriction, including without limitation the rights to
    use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies
    of the Software, and to permit persons to whom the Software is furnished to do
    so, subject to the following conditions:

    The above copyright notice and this permission notice shall be included in all
    copies or substantial portions of the Software.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
    OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
    SOFTWARE.
----------------------------------------------------------------------------------------------------------------------------------------------------------
"""
import numpy as np
import numpy.ctypeslib as npct
from ctypes import c_int
from ctypes import c_int16
from ctypes import c_uint64
from ctypes import c_uint32
from ctypes import c_uint16
from ctypes import c_float

# input type for genelife c to python array functions
# must be a short/long/_ unsigned int array, with single dimension that is contiguous
uint16_array = npct.ndpointer(dtype=np.uint16, ndim=1, flags='CONTIGUOUS')
uint64_array = npct.ndpointer(dtype=np.uint64, ndim=1, flags='CONTIGUOUS')
uint_array = npct.ndpointer(dtype=np.uint32, ndim=1, flags='CONTIGUOUS')
int_array = npct.ndpointer(dtype=np.int32, ndim=1, flags='CONTIGUOUS')

# communicate component type for connected components if not using numpy
# from ctypes import *
# class COMPONENT(Structure):
#    _fields_ = [('N',c_uint16),('S',c_uint16),('W',c_uint16),('E',c_uint16),('log2n',c_uint16),('patt',c_uint16),('lastrc',c_uint16),('dummy',c_uint16),
#                ('label',c_uint32),('pixels',c_uint32),('quad',c_uint64),('gcolor',c_float),('reserve',c_uint32)]

# genedata type to retrieve hashed genes from C using numpy and ctypes
genedtype=[('popcount',c_uint32),('firsttime',c_uint16),('recenttime',c_uint16),('lasttime',c_uint16),('lastextinctiontime',c_int16),#    1 unsigned, 3 short unsigned, 1 short signed int
           ('activity',c_uint32),('nextinctions',c_uint32),('dummy',c_uint32),('gene',c_uint64),('firstancestor',c_uint64),('recentancestor',c_uint64)] #    2 unsigned 32 bit and 3 unsigned 64 bit int
gene_array = npct.ndpointer(dtype=genedtype, ndim=1, flags=['CONTIGUOUS','ALIGNED'])

# component type to retrieve connected components from C using numpy and ctypes
compdtype=[('N',c_uint16),('S',c_uint16),('W',c_uint16),('E',c_uint16),('log2n',c_uint16),('patt',c_uint16),('lastrc',c_uint16),('dummy',c_uint16),
           ('label',c_uint32),('pixels',c_uint32),('quad',c_uint64),('gcolor',c_float),('reserve',c_uint32)]
comp_array = npct.ndpointer(dtype=compdtype, ndim=1, flags=['CONTIGUOUS','ALIGNED'])

# quadnode type to retrieve quadnodes from C using numpy and ctypes
quaddtype=[('hashkey',c_uint64),('nw',c_uint64),('ne',c_uint64),('sw',c_uint64),('se',c_uint64),('isnode',c_uint16),('size',c_uint16),
           ('activity',c_uint32),('pop1s',c_uint32),('firsttime',c_uint32),('lasttime',c_uint32),('reserve',c_uint32)]
quad_array = npct.ndpointer(dtype=quaddtype, ndim=1, flags=['CONTIGUOUS','ALIGNED'])

# smallpatt type to retrieve smallpatts from C using numpy and ctypes
smallpattdtype=[('topactivity',c_uint32),('activity',c_uint32),('firsttime',c_uint32),('lasttime',c_uint32)]
smallpatt_array = npct.ndpointer(dtype=smallpattdtype, ndim=1, flags=['CONTIGUOUS','ALIGNED'])

# load the library, using numpy mechanisms
try:
    libcd = npct.load_library("libgenelife", ".")
except:
    try:
        libcd = npct.load_library("libgenelife", "..")
    except:
        try:
            libcd = npct.load_library("libgenelife", "../..")
        except:
            try:
                libcd = npct.load_library("libgenelife", "./genelifepy")
            except:
                print("Error: genelife library not found at . , .. , ../.. or ./genelifepy")
            
# setup the return types and argument types
libcd.get_log2N.restype = c_int
libcd.get_log2N.argtypes = []
libcd.genelife_update.restype = None
libcd.genelife_update.argtypes = [ c_int, c_int, c_int]
libcd.initialize.restype = None
libcd.initialize.argtypes = [int_array, c_int, int_array, c_int]
libcd.initialize_planes.restype = None
libcd.initialize_planes.argtypes = [int_array, c_int]
libcd.countspecies1.restype = None
libcd.countspecies1.argtypes = [uint64_array, uint64_array, c_int]
libcd.countspecies.restype = None
libcd.countspecies.argtypes = []
libcd.countspecieshash.restype = None
libcd.countspecieshash.argtypes = []
#libcd.print_gol.restype = None
#libcd.print_gol.argtypes = [uint64_array, c_int, c_int]
#libcd.printscreen.restype = None
#libcd.printscreen.argtypes = [uint64_array, uint64_array, c_int, c_int]
libcd.get_histo.restype = None
libcd.get_histo.argtypes = [int_array, c_int]
libcd.get_curgol.restype = None
libcd.get_curgol.argtypes = [uint64_array,c_int]
libcd.get_curgolg.restype = None
libcd.get_curgolg.argtypes = [uint64_array,c_int]
libcd.get_curgolbr.restype = None
libcd.get_curgolbr.argtypes = [uint64_array,uint64_array,c_int]
libcd.get_curgolgstats.restype = None
libcd.get_curgolgstats.argtypes = [uint64_array,c_int]
libcd.get_gliderinfo.restype = None
libcd.get_gliderinfo.argtypes = [uint64_array,c_int]
libcd.get_glider_count.restype = c_int
libcd.get_glider_count.argtypes = [uint_array,c_int]
libcd.stash.restype = None
libcd.stash.argtypes = []
libcd.label2stash.restype = None
libcd.label2stash.argtypes = [c_int]
libcd.unstash.restype = None
libcd.unstash.argtypes = None
libcd.get_nspecies.argtypes = None
libcd.get_nspecies.restype = c_int
libcd.get_nlive.argtypes = None
libcd.get_nlive.restype = c_int
libcd.get_genealogydepth.argtypes = None
libcd.get_genealogydepth.restype = c_int
libcd.get_curtime.argtypes = None
libcd.get_curtime.restype = c_int
libcd.get_stats.restype = None
libcd.get_stats.argtypes = [int_array, int_array, int_array, int_array, c_int]
libcd.get_activities.restype = c_int
libcd.get_activities.argtypes = [uint64_array, int_array, c_int]
libcd.get_all_activities.restype = c_int
libcd.get_all_activities.argtypes = [uint64_array, int_array, c_int]
libcd.get_small_activities.restype = c_int
libcd.get_small_activities.argtypes = [uint64_array, int_array, c_int]
libcd.get_quad_activities.restype = c_int
libcd.get_quad_activities.argtypes = [uint64_array, int_array, c_int]
libcd.get_all_small_activities.restype = c_int
libcd.get_all_small_activities.argtypes = [uint64_array, int_array, c_int]
libcd.get_all_quad_activities.restype = c_int
libcd.get_all_quad_activities.argtypes = [uint64_array, int_array, c_int]
libcd.get_acttrace.restype = None
libcd.get_acttrace.argtypes = [uint64_array,c_int]
libcd.get_poptrace.restype = None
libcd.get_poptrace.argtypes = [uint64_array,c_int]
libcd.get_genealogytrace.restype = None
libcd.get_genealogytrace.argtypes = [uint64_array, c_int]
libcd.get_nnovelcells.restype = None
libcd.get_nnovelcells.argtypes = [uint_array, c_int]
libcd.get_sorted_popln_act.restype = c_int
libcd.get_sorted_popln_act.argtypes = [int_array, uint64_array, int_array, int_array]
libcd.get_connected_comps.restype = c_int
libcd.get_connected_comps.argtypes = [uint_array, uint_array, c_int, c_int]
libcd.get_ncomponents.restype = c_int
libcd.get_ncomponents.argtypes = None
libcd.get_components.restype = c_int
libcd.get_components.argtypes = [comp_array, c_int]
libcd.get_smallpatts.restype = c_int
libcd.get_smallpatts.argtypes = [smallpatt_array, c_int]
libcd.get_quadnodes.restype = c_int
libcd.get_quadnodes.argtypes = [quad_array, c_int]
libcd.get_genes.restype = c_int
libcd.get_genes.argtypes = [gene_array, c_int]
libcd.get_genealogies.restype = c_int
libcd.get_genealogies.argtypes = [gene_array,c_int]
libcd.get_shist.restype = None
libcd.get_shist.argtypes = [int_array]
libcd.colorgenes.restype = None
libcd.colorgenes.argtypes = [int_array, c_int, c_int, c_int, c_int]
libcd.set_colorfunction.restype = None
libcd.set_colorfunction.argtypes = [c_int]
libcd.setget_act_ymax.restype = c_int
libcd.setget_act_ymax.argtypes = [c_int]
libcd.setget_act_ymaxq.restype = c_int
libcd.setget_act_ymaxq.argtypes = [c_int]
libcd.set_selectedgene.restype = None
libcd.set_selectedgene.argtypes = [c_uint64]
libcd.set_offsets.restype = None
libcd.set_offsets.argtypes = [c_int,c_int,c_int]
libcd.set_quadrant.restype = None
libcd.set_quadrant.argtypes = [c_int]
libcd.set_randominflux.restype = None
libcd.set_randominflux.argtypes = [c_int]
libcd.set_rbackground.restype = None
libcd.set_rbackground.argtypes = [c_int,c_int]
libcd.set_repscheme_bits.restype = c_uint32
libcd.set_repscheme_bits.argtypes = [c_int,c_int,c_int, uint_array]
libcd.set_repscheme.restype = None
libcd.set_repscheme.argtypes = [c_uint32]
libcd.set_rulemod.restype = None
libcd.set_rulemod.argtypes = [c_uint32]
libcd.set_surviveover64.restype = None
libcd.set_surviveover64.argtypes = [uint_array, c_int]
libcd.set_vscrolling.restype = None
libcd.set_vscrolling.argtypes = []
libcd.set_noveltyfilter.restype = None
libcd.set_noveltyfilter.argtypes = []
libcd.set_activity_size_colormode.restype = None
libcd.set_activity_size_colormode.argtypes = []
libcd.set_gcolors.restype = None
libcd.set_gcolors.argtypes = []
libcd.set_seed.restype = None
libcd.set_seed.argtypes = [c_int]
libcd.set_nbhist.restype = None
libcd.set_nbhist.argtypes = [c_int]
libcd.set_genealogycoldepth.restype = None
libcd.set_genealogycoldepth.argtypes = [c_int]
libcd.set_ancestortype.restype = None
libcd.set_ancestortype.argtypes = [c_int]
libcd.set_info_transfer_h.restype = None
libcd.set_info_transfer_h.argtypes = [c_int,c_int]
libcd.set_activityfnlut.restype = None
libcd.set_activityfnlut.argtypes = [c_int]
libcd.set_colorupdate1.restype = None
libcd.set_colorupdate1.argtypes = [c_int]
libcd.set_colorfunction2.restype = None
libcd.set_colorfunction2.argtypes = [c_int]

def get_log2N():
    return libcd.get_log2N()

def genelife_update(nsteps, nhist, nstat):
    return libcd.genelife_update(nsteps, nhist, nstat)

def get_curgol(gol):
    return libcd.get_curgol(gol, int(len(gol)))

def get_curgolg(golg):
    return libcd.get_curgolg(golg, int(len(golg)))

def get_curgolbr(golb,golr):
    return libcd.get_curgolbr(golb, golr, int(len(golr)))

def get_curgolgstats(golgstats):
    return libcd.get_curgolgstats(golgstats, int(len(golgstats)))

def get_gliderinfo(gliderinfo):
    return libcd.get_gliderinfo(gliderinfo, int(len(gliderinfo)))

def get_glider_count(golrstats):
    return libcd.get_glider_count(golrstats, int(len(golrstats)))

def get_nspecies():
    return libcd.get_nspecies()

def get_curtime():
    return libcd.get_curtime()

def initialize(runparams,simparams):
    return libcd.initialize(runparams, len(runparams), simparams, len(simparams))

def initialize_planes(offsets):
    return libcd.initialize_planes(offsets, int(len(offsets)))

def countspecies1(gol, golg):
    return libcd.countspecies1(gol, golg,  len(golg))

def countspecies():
    return libcd.countspecies()

def countspecieshash():
    return libcd.countspecieshash()

#def print_gol( gol, N):
#    return libcd.print_gol( gol, N, len(gol))

#def printscreen( gol, golg, N):
#    return libcd.printscreen( gol, golg, N, len(gol))

def get_histo(gol):
    return libcd.getconfigs( histo, len(histo))

def get_stats(livesites,genotypes,stepstats,configstats,nstats):
    return libcd.get_stats(livesites,genotypes,stepstats,configstats,nstats)

def get_activities(genes,activities):
    return libcd.get_activities(genes,activities,int(len(genes)))

def get_all_activities(genes,activities):
    return libcd.get_all_activities(genes,activities,int(len(genes)))

def get_quad_activities(quads,activities):
    return libcd.get_quad_activities(quads,activities,int(len(quads)))

def get_small_activities(smalls,activities):
    return libcd.get_small_activities(smalls,activities,int(len(smalls)))

def get_acttrace(acttrace):
    return libcd.get_acttrace(acttrace, int(len(acttrace)))

def get_poptrace(poptrace):
    return libcd.get_poptrace(poptrace, int(len(poptrace)))

def get_genealogytrace(genealogytrace):
    return libcd.get_genealogytrace(genealogytrace, int(len(genealogytrace)))

def get_nnovelcells(nnovelcells):
    return libcd.get_nnovelcells(nnovelcells,int(len(nnovelcells)))

def get_genealogies(genealogydat):
    return libcd.get_genealogies(genealogydat, int(len(genealogydat)))

def get_nspecies():
    return libcd.get_nspecies()

def get_genealogydepth():
    return libcd.get_genealogydepth()

def get_sorted_popln_act( gindices, genes, popcnts, activities):
    return libcd.get_sorted_popln_act(gindices, genes, popcnts, activities)

def get_connected_comps(connlabel,connlen,x,y):
    return libcd.get_connected_comps(connlabel,connlen,x,y)

def get_ncomponents():
    return libcd.get_ncomponents()

def get_components(components):
    return libcd.get_components(components, int(len(components)))

def get_smallpatts(smallpatts):
    return libcd.get_smallpatts(smallpatts, int(len(smallpatts)))

def get_quadnodes(quadnodes):
    return libcd.get_quadnodes(quadnodes, int(len(quadnodes)))

def get_genes(genelist):
    return libcd.get_genes(genelist, int(len(genelist)))

def unstash():
    return libcd.unstash()

def get_nlive():
    return libcd.get_nlive()

def get_shist(shist):
    return libcd.get_shist(shist)

def colorgenes(cgolg,colorfunction, winnr, nfrstep):
    return libcd.colorgenes( cgolg, len(cgolg), colorfunction, winnr, nfrstep)

def set_colorfunction(colorfunctionval):
    return libcd.set_colorfunction(colorfunctionval)

def setget_act_ymax(ymax):
    return libcd.setget_act_ymax(ymax)

def setget_act_ymaxq(ymaxq):
    return libcd.setget_act_ymaxq(ymaxq)

def set_selectedgene(gene):
    return libcd.set_selectedgene(gene)

def set_offsets(dx,dy,dt):
    return libcd.set_offsets(dx,dy,dt);

def set_quadrant(quadrant):
    return libcd.set_quadrant(quadrant)

def set_randominflux(randominflux):
    return libcd.set_randominflux(randominflux)

def set_rbackground(rbackground,randominflux):
    return libcd.set_rbackground(rbackground,randominflux)

def set_repscheme_bits(quadrant, x, y, surviveover):
    return libcd.set_repscheme_bits(quadrant, x, y, surviveover)

def set_repscheme(repscheme):
    return libcd.set_repscheme(repscheme)

def set_rulemod(rulemod):
    return libcd.set_rulemod(rulemod)

def set_surviveover64(surviveover):
    return libcd.set_surviveover64( surviveover, len(surviveover))

def set_vscrolling():
    return libcd.set_vscrolling()
    
def set_noveltyfilter():
    return libcd.set_noveltyfilter()

def set_activity_size_colormode():
    return libcd.set_activity_size_colormode()

def set_gcolors():
    return libcd.set_gcolors()

def set_seed(seed):
    return libcd.set_seed(seed)

def set_nbhist(nbhist):
    return libcd.set_nbhist(nbhist)

def set_genealogycoldepth(genealogycoldepth):
    return libcd.set_genealogycoldepth(genealogycoldepth)

def set_ancestortype(ancestortype):
    return libcd.set_ancestortype(ancestortype)

def stash():
    return libcd.stash();

def label2stash(cumul):
    return libcd.label2stash(cumul)

def set_info_transfer_h(do_info_transfer,nbhood):
    return libcd.set_info_transfer_h(do_info_transfer,nbhood);
    
def set_activityfnlut(activityfnlut):
    return libcd.set_activityfnlut(activityfnlut)

def set_colorupdate1(update1):
    return libcd.set_colorupdate1(update1)

def set_colorfunction2(colorfunction2):
    return libcd.set_colorfunction2(colorfunction2)
