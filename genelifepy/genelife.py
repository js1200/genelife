""" genelife.py
    profect genelife

------------------------------------------------------------- copyright ----------------------------------------------------------------------------------
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
---------------------------------------------------------------------------------------------------------------------------------------------------------
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mpcolors
from   matplotlib.colors import ListedColormap
import matplotlib.animation as animation
import matplotlib
import genelife_c_interface as genelife
import sdl2
import sdl2.ext
import sdl2.timer
import sdl2.sdlttf  # install library from https://www.libsdl.org/projects/SDL_ttf/
import ctypes
# import python_utilities

# converted to python 3 from python 2.7  in Mar 2019
# converted from pygame to PySDL2 with SDL2 in Mar 2019
version = sdl2.SDL_version()
sdl2.SDL_GetVersion(ctypes.byref(version))
print("Running with SDL version %d.%d.%d." % (version.major, version.minor, version.patch))

# Variables that are read-only in python notebook
log2N = genelife.get_log2N()                #  log2N=7 => N=128 get value from C library which is where it must be changed
N = 2**log2N
N2 = N*N
Nmask = N-1
Width = N
Height = N
nNhist = 20

gol = np.zeros(N2,np.uint64)
golg = np.zeros(N2,np.uint64)
golb = np.zeros(N2,np.uint64)
golr = np.zeros(N2,np.uint64)
golgstats = np.zeros(N2,np.uint64)

nspecies = 0
ncomponents = 0
connlabel = np.zeros(N2,np.uint32)
connlen = np.zeros(N2//4,np.uint32)
nnovelcells = np.zeros(N*nNhist,np.uint32)
                                            # graphics
cgrid = np.zeros((N,N),np.int32)
cgridt = np.zeros((N,N),np.int32)
cgolg =np.zeros(N2,np.int32)
colorfunction = 0
nfrstep = 0
surface = None
surfacex1 = None
surfacex2 = None
window = None
windowx1 = None
update1 = True
scalex2 = False
rescale = True
ncanon=[]
cancol=[]
caption = ""
dispinit = False
display = [window,surface,cgrid,caption,dispinit]

cgrid2 = np.zeros((N,N),np.int32)
cgolg2 =np.zeros(N2,np.int32)
window2 = None
window2x1 = None
surface2 = None
surface2x1 = None
surface2x2 = None
caption2 = ""
dispinit2 = False
display2 = [window2,surface2,cgrid2,caption2,dispinit2]

render2 = False                            # whether to use renderer on surface 2 : NB frees surface2 after defining texture2
update2 = True
renderer2 = None
texture2 = None
image2 = None
message2 = None
font2 = None
textColor2 = None
windowID2 = None


grect = sdl2.SDL_Rect()
grect.x = Width//2
grect.y = Height//2
grect.w = Width//2
grect.h = 20

grect1 = sdl2.SDL_Rect()
grect1.x = 0
grect1.y = Height
grect1.w = Width
grect1.h = 16

updatesenabled = True
mat = []

gogo = True
pixeldat = ""
paramdat = ""
mouseclicked = False
mouseclicked2 = False
pause = 0
ymax = 10000
oldymax = ymax
ymaxq = 10000
oldymaxq = ymaxq
maxPlane = 4
offdx = offdy = offdt = 0
quadrants = -1
gcolor = 0
                                         # counter and toggle initialization

cnt = 0
framenr = 0
mstime = 0
framerate=0.0
savecnt = 0                              # counter for saved images
randominflux = 0
vscrolling = 0
noveltyfilter = 0
activity_size_colormode = 0
ancestortype = 0
info_transfer_h = 0
activityfnlut = 0
                                         # parameter initialization
runparams = np.zeros(10,np.int32)        # 10 parameters passed to C
simparams = np.zeros(6,np.int32)         # 6 parameters passed to C
nrun=1; ndisp=1000; nskip=0; niter=1;    # simulation time stepping parameters: nrun CA updates per step, ndisp nr steps to display before skip,
                                         # nskip nr of CA updates to skip over display, niter nr of repeats of disp-skip cycle
nhist = 0                                # set to n to turn on histogram configurations every nth step
nbhist = -1                              # set block for display of traces of activity and population
nNhist = 20                              # number of additional blocks for trace memory
genealogycoldepth = 0                    # set depth of genealogical ancestor for gene display in colorfunction 11
nstat = 0                                # set to n to turn on statistics trace every nth step
rulemod = runparams[0] = 1               # 0,1 whether to allow GoL rule modifications
                                         # with rulemod 1 2-live-nb birth, 3-live-nb non-birth & non-survival possible
repscheme = runparams[1] = 8             # repscheme bit 3 (val 0x8) determines whether random choice of ancestor amongst live neighbours
                                         # repscheme mod 8 i.e. 0-7 determines selection scheme based on gene
                                         # 0 minimum gene as value  # 1 maximum gene as value
                                         # 2 minimum number of ones # 3 maximum number of ones
                                         # 4 neutral selection # 5 neutral but different selection
                                         # 6 penalty function -1 for a survival rule -2 for a birth rule  # 7 not allowed
selection = runparams[2] = 10            # fitness for 2 live neighbor rule : 0-15 see subgenelife.c code
overwritemask = runparams[3]= 0x3        # whether to overwrite existing genes and allow birth
survivalmask = runparams[4] = 0x06       # for selection=8-13 this is the GoL survival mask
birthmask = runparams[7] = 0x04          # for selection=8-13 this is the GoL birth mask
ancselectmask = runparams[8] = 0xff      # bit mask for enabling gene-selective choice of ancestor for different birth rules
colorfunction = runparams[5] = 0         # color function 0(hash), >=1(fnal), 2 nongolstate 3 report of golgstats in 3 colours and 3 levels each
                                         # 4 activities 5 populations 6 genealogy steps 7 genealogy temporal with activity scaled colors
                                         # 8 glider detection 9 connected component labelling 10 connected component activities
                                         # 11 genealogy depth 12 genetic glider detection
colorfunction2 = runparams[9] = -1       # colorfunction for 2nd window: -1 same as first window
initfield = runparams[6] = 100           # 1 init via 32x32 genepat.dat, n>1 init via nxn rand array
nlog2pmut = simparams[0] = 8             # log2 gene mutation probability (0 or >56 means no mutation)
initial1density = simparams[1] =  16384  # initial 1 density in GOL state
                                         # 16384 = nearest to half of guaranteed C rand max value 32767 = 2**15 - 1
initialrdensity = simparams[2] = 32768   # initial density of random genes, 0 if all initial gens are startgenes
ncoding = simparams[3] = 0               # for selection 10, non zero value means grow plane community from 0
                                         # otherwise (selection<10) no of bits used to encode valid connection functions 1-16
                                         # for selection==8, lut, ncoding 1,2,3 bits per lut entry : 0 implies 3.
startgenechoice = simparams[4] = 8       # initialize genes to startgene number 0-8 : 8 is random choice of 0-7
ranseed = simparams[5] = 1234            # initial seed for random number generator

                                         # offset initialization
"""offsets = [[-1, 0,-1],
           [ 1, 0,-1],
           [ 0,-1,-1],
           [ 0, 1,-1]]"""
offsets = [[ 1, 1, 0],
           [ 0, 0,-1],
           [ 0, 0,-2],
           [ 0, 0,-3],
           [ 0, 0,-4],
           [ 0, 0,-5],
           [ 0, 0,-6],
           [ 0, 0,-7]]
numHis = len(offsets)
histo=np.zeros(numHis,np.uint64)
flatoff =  [x for sublist in offsets for x in sublist]
npoffsets = np.array(flatoff,np.int32)

# setup of color map : black for 0, colors for 1 to LEN+1 or 257 for colormethod 0 or 1
#-----------------------------------------------------------------------------------------------------------
def rand_cmap(nlabels, type='bright', first_color_black=True, last_color_black=False, verbose=True):

    """
    Creates a random colormap to be used together with matplotlib. Useful for segmentation tasks
    :param nlabels: Number of labels (size of colormap)
    :param type: 'bright' for strong colors, 'soft' for pastel colors
    :param first_color_black: Option to use first color as black, True or False
    :param last_color_black: Option to use last color as black, True or False
    :param verbose: Prints the number of labels and shows the colormap. True or False
    :return: colormap for matplotlib
    Author of this color function: Delestro, stackoverflow or https://github.com/delestro/rand_cmap
    """
    from matplotlib.colors import LinearSegmentedColormap
    import colorsys
    import numpy as np

    if type not in ('bright', 'soft', 'grad'):
        print ('Please choose "bright" or "soft" or "grad" for type')
        return
    if verbose:
        print(('Number of labels: ' + str(nlabels)))
    np.random.seed(123456)
    # Generate color map for bright colors, based on hsv
    if type == 'bright':
        randHSVcolors = [(np.random.uniform(low=0.0, high=1),
                      np.random.uniform(low=0.2, high=1),
                      np.random.uniform(low=0.9, high=1)) for i in range(nlabels)]
        randRGBcolors = []
        for HSVcolor in randHSVcolors:  # Convert HSV list to RGB
            randRGBcolors.append(colorsys.hsv_to_rgb(HSVcolor[0], HSVcolor[1], HSVcolor[2]))
        if first_color_black:
            randRGBcolors[0] = [0, 0, 0]
#            randRGBcolors[1] = [0, 0, 1] # this color otherwise is black
        if last_color_black:
            randRGBcolors[-1] = [0, 0, 0]
        random_colormap = LinearSegmentedColormap.from_list('new_map', randRGBcolors, N=nlabels)

    if type == 'grad':
#        randRGBcolors = np.array((np.linspace(0.0, 1.0, num=256),
#                                  np.linspace(1.0, 0.0, num=256),
#                                  np.concatenate(np.linspace(1.0, 0.0, num=128),np.linspace(0.0, 1.0, num=128)))).T
        randRGBcolors =[((i-1)/255.,(255.-(i-1))/255.,(i-1)/128. if i<=128 else (255.-(i-1))/128.) for i in range(nlabels)]
        randRGBcolors[0] = [0, 0, 0]
        random_colormap = LinearSegmentedColormap.from_list('new_map', randRGBcolors, N=nlabels)
        
    # Generate soft pastel colors, by limiting the RGB spectrum
    if type == 'soft':
        low = 0.6
        high = 0.95
        randRGBcolors = [(np.random.uniform(low=low, high=high),
                      np.random.uniform(low=low, high=high),
                      np.random.uniform(low=low, high=high)) for i in range(nlabels)]
        if first_color_black:
            randRGBcolors[0] = [0, 0, 0]
        if last_color_black:
            randRGBcolors[-1] = [0, 0, 0]
        random_colormap = LinearSegmentedColormap.from_list('new_map', randRGBcolors, N=nlabels)

    # Display colorbar
    if verbose:
        from matplotlib import colors, colorbar
        from matplotlib import pyplot as plt
        fig, ax = plt.subplots(1, 1, figsize=(15, 0.5))
        bounds = np.linspace(0, nlabels, nlabels + 1)
        norm = colors.BoundaryNorm(bounds, nlabels)
        cb = colorbar.ColorbarBase(ax, cmap=random_colormap, norm=norm, spacing='proportional', ticks=None,
                               boundaries=bounds, format='%1i', orientation='horizontal')

    return random_colormap
#-----------------------------------------------------------------------------------------------------------

def colorgrid(colorfunction,cgolg,cgrid,winnr,nfrstep):
    """ colors array according to grid and genegrid using colormethod"""
    global N,cgridt
    genelife.colorgenes(cgolg,colorfunction,winnr,nfrstep)
    cgridt=np.reshape(cgolg,(N,N)).T
    cgrid[:,0:N] = cgridt   # is there a faster version of this copy that moves the data?
    return
#-----------------------------------------------------------------------------------------------------------

def packrepscheme(repscheme,survivalmask,overwritemask):
    if survivalmask<4 and overwritemask<4:
        repscheme = repscheme + (survivalmask<<24) + (overwritemask<<26)
    else:
        print("Error: can't pack masks, they are too large!")
    return(repscheme)

#-------------------------------------- SDL implemented graphics routines ----------------------------------
#  Draws filled rectangle = [x, y, width, height] on the passed surface
def draw_rect(surface, color, rectangle):
        sdl_color = sdl2.ext.Color(color[0],color[1],color[2])
        sdl2.ext.fill(surface, color, rectangle)
#  Get mouse coordinates
def mouse_get_pos():
    mouse32bitstate=sdl2.mouse.SDL_GetMouseState(ctypes.byref(x),ctypes.byref(y))
    return((x,y))
#  Set window caption
def set_caption(window, title):
    # sdl2.SDL_SetWindowTitle(window, title)
    window.title = title
    # sdl2.ext.Window.DEFAULTPOS = (1000, 32)

#-----------------------------------------------------------------------------------------------------------

def init_button_arrays():
    """ initialize information for button area at base of display"""
    global ncanon
                                                            # initialize lists to empty
    ncanon = []                                             # lists of number of successive buttons with the same color
    cancolors = []                                          # lists of button colors
    cancol =[]                                              # expanded lists of button colors, by button
                                                            # no of entries per color region
    ncanon.append([2,2,2,2,2,4,2,2,2])                      # selection 0-7
    ncanon.append([8,8,8,4,4,3,1])                          # selection 8,9
    ncanon.append([2,3,4,5,4,3,2,2,3,4,5,4,3,2,8,4,4,3,1])  # selection 10,11
    ncanon.append([4,7,10,7,4,4,7,10,7,4,8,4,4,3,1])        # selection 12,13
    ncanon.append([1,2,6,10,13,1,2,6,10,13,8,4,4,3,1])      # selection 14,15
                                                            # colors [R,G,B] for different color regions for buttons, colorvals must be < 128 to allow 2x highlight
    cancolors.append([[0,100,0],[0,50,100],[0,80,80],[0,100,50],[100,100,0],[50,100,0],[100,0,100],[0,0,127],[100,0,0]]) # selection 0-7
    cancolors.append([[0,0,127],[0,100,0],[100,0,0],[100,100,100],[50,100,50],[100,50,50],[100,100,100]]) # selection 8,9
    cancolors.append([[0,0,127],[50,0,127],[80,0,120],[100,0,120],[80,0,120],[50,0,127],[0,0,127],
                      [0,127,0],[50,127,0],[60,120,0],[100,120,0],[80,120,0],[50,120,0],[0,127,0],[100,0,0],[100,100,100],[50,100,50],[100,50,50],[100,100,100]]) # selection 10,11
    cancolors.append([[50,0,127],[80,0,120],[100,0,120],[80,0,120],[50,0,127],[50,127,0],[80,120,0],[100,120,0],[80,120,0],[50,127,0],[100,0,0],[100,100,100],[50,100,50],[100,50,50],[100,100,100]])  # selection 12,13
    cancolors.append([[50,0,127],[80,0,120],[100,0,120],[80,0,120],[50,0,127],[50,127,0],[80,120,0],[100,120,0],[80,120,0],[50,127,0],[100,0,0],[100,100,100],[50,100,50],[100,50,50],[100,100,100]])  # selection 14,15
                                                            # lists of colors for individual buttons expanded from above, first initialize to zero
    cancol.append(np.zeros((20,3),np.int32))
    cancol.append(np.zeros((3*8+12,3),np.int32))
    cancol.append(np.zeros((2*23+8+12,3),np.int32))
    cancol.append(np.zeros((2*32+8+12,3),np.int32))
    cancol.append(np.zeros((2*32+8+12,3),np.int32))
    
    for l in range(len(ncanon)):                            # buttons for different selection schemes
        k=0
        for j in range(len(ncanon[l])):
            for i in range(ncanon[l][j]):
                cancol[l][k]=cancolors[l][j]
                k = k+1
    return(cancol)
#-----------------------------------------------------------------------------------------------------------

def init_buttons(surface):    # initialize parameter buttons
    global repscheme,survivalmask,birthmask,overwritemask,ancselectmask,selection,ncoding
    global scalex2
    global Height,Width
    global log2N
    global ncanon

    if scalex2:
        sc = 1
    else:
        sc = 2
    cancol=init_button_arrays()
    draw_rect(surface,[50,50,50],[0,Height+4,Width,10*sc])
    # draw_rect(surface,[50,50,50],[0,Height+6,Width,7*sc])

    if selection<8:
        for k in range(20):
            if k<16:
                bit = (repscheme>>k)&0x1
            elif k<18:
                bit = (survivalmask>>(k-16))&0x1
            elif k<20:
                bit = (overwritemask>>(k-18))&0x1
            draw_rect(surface,cancol[0][k]*(1+bit),[k<<(log2N-6),Height+6,3*sc,3*sc])
        j = 0;
        for k in range(len(ncanon[0])):                       # white grouping markers
            draw_rect(surface,[200,200,200],[(j<<(log2N-6))-1 if j else 0,Height+4,sc,sc])
            j = j+ncanon[0][k]
    elif selection<10:
        for k in range(8):
            draw_rect(surface,cancol[1][k]*(1+((survivalmask>>k)&0x1)),[k<<(log2N-6),Height+6,3*sc,3*sc])
            draw_rect(surface,cancol[1][k+8]*(1+((birthmask>>(k))&0x1)),[(k+8)<<(log2N-6),Height+6,3*sc,3*sc])
            draw_rect(surface,cancol[1][k+16]*(1+((overwritemask>>(k))&0x1)),[(k+16)<<(log2N-6),Height+6,3*sc,3*sc])
        for k in range(12):
            draw_rect(surface,cancol[1][k+24]*(1+((repscheme>>k)&0x1)),[(k+24)<<(log2N-6),Height+6,3*sc,3*sc])
        j = 0;
        for k in range(len(ncanon[1])):                       # white grouping markers
            draw_rect(surface,[200,200,200],[(j<<(log2N-6))-1 if j else 0,Height+4,sc,sc])
            j = j+ncanon[1][k]
    elif selection<12:
        # draw_rect(surface,[200,200,200],[(23<<(log2N-6))-1,Height+6,1,9])
        for k in range(23):
            draw_rect(surface,cancol[2][k]*(1+((survivalmask>>k)&0x1)),[k<<(log2N-6),Height+6,3*sc,3*sc])
            draw_rect(surface,cancol[2][k+23]*(1+((birthmask>>(k))&0x1)),[(k+23)<<(log2N-6),Height+6,3*sc,3*sc])
        for k in range(8):
            draw_rect(surface,cancol[2][k+46]*(1+((overwritemask>>(k))&0x1)),[(k+46)<<(log2N-6),Height+6,3*sc,3*sc])
        for k in range(12):
            draw_rect(surface,cancol[2][k+54]*(1+((repscheme>>k)&0x1)),[k<<(log2N-6),Height+8+3*sc,3*sc,3*sc])
        j = 0;
        for k in range(len(ncanon[2])):                       # white grouping markers
            jj = (j if j<54 else j-54);
            draw_rect(surface,[200,200,200],[(jj<<(log2N-6))-1 if jj else 0,Height+4 if j<54 else Height+6+3*sc,sc,sc])
            j = j+ncanon[2][k]
    elif selection<14:
        # draw_rect(surface,[200,200,200],[(32<<(log2N-6))-1,Height+6,1,9])
        for k in range(32):
            draw_rect(surface,cancol[3][k]*(1+((survivalmask>>k)&0x1)),[k<<(log2N-6),Height+6,3*sc,3*sc])
            draw_rect(surface,cancol[3][k+32]*(1+((birthmask>>(k))&0x1)),[(k+32)<<(log2N-6),Height+6,3*sc,3*sc])
        for k in range(8):
            draw_rect(surface,cancol[3][k+64]*(1+((overwritemask>>(k))&0x1)),[(k+32)<<(log2N-6),Height+8+3*sc,3*sc,3*sc])
        for k in range(12):
            draw_rect(surface,cancol[3][k+72]*(1+((repscheme>>k)&0x1)),[k<<(log2N-6),Height+8+3*sc,3*sc,3*sc])
        j = 0;
        for k in range(len(ncanon[3])):                       # white grouping markers
            if j<64:
                jj = j;
            elif j<72:
                jj = j-32;
            else:
                jj = j-72;
            draw_rect(surface,[200,200,200],[(jj<<(log2N-6))-1 if jj else 0,Height+4 if j<64 else Height+6+3*sc,sc,sc])
            j = j+ncanon[3][k]
    elif selection<16:
        # draw_rect(surface,[200,200,200],[(32<<(log2N-6))-1,Height+6,1,9])
        for k in range(32):
            draw_rect(surface,cancol[4][k]*(1+((survivalmask>>k)&0x1)),[k<<(log2N-6),Height+6,3*sc,3*sc])
            draw_rect(surface,cancol[4][k+32]*(1+((birthmask>>(k))&0x1)),[(k+32)<<(log2N-6),Height+6,3*sc,3*sc])
        for k in range(8):
            draw_rect(surface,cancol[4][k+64]*(1+((overwritemask>>(k))&0x1)),[(k+32)<<(log2N-6),Height+8+3*sc,3*sc,3*sc])
        for k in range(12):
            draw_rect(surface,cancol[4][k+72]*(1+((repscheme>>k)&0x1)),[k<<(log2N-6),Height+8+3*sc,3*sc,3*sc])
        j = 0;
        for k in range(len(ncanon[4])):                       # white grouping markers
            if j<64:
                jj = j;
            elif j<72:
                jj = j-32;
            else:
                jj = j-72;
            draw_rect(surface,[200,200,200],[(jj<<(log2N-6))-1 if jj else 0,Height+4 if j<64 else Height+6+3*sc,sc,sc])
            j = j+ncanon[4][k]
    return(cancol)
#-----------------------------------------------------------------------------------------------------------
def display_init():
    global scalex2,Width,Height,rescale,cnt
    global window,windowx1,surface,surfacex1,surfacex2,cgrid,caption,dispinit

    if not dispinit:
        sdl2.ext.init()
    dispinit = True
    # clock = sdl2.ext.time.Clock()  # only works in modified version of pysdl2 see https://lukems.github.io/py-sdl2/modules/sdl2ext_time.html
    caption = "Gene Life at iteration %d" % cnt
    if (Height <= 512 and rescale):            # Not yet ported well to SDL2
        scalex2 = True
        #window = sdl2.ext.Window(caption,(2*Width, 2*(Height+16)),(1000,60),sdl2.SDL_PIXELFORMAT_BGRA8888)     # opens sdl2 window, add flags for last parameter
        windowx1 = sdl2.ext.Window(caption,(Width, Height+16),(500,60),sdl2.SDL_WINDOW_HIDDEN)     # opens sdl2 window but does not show it
        window = sdl2.ext.Window(caption,(2*Width, 2*(Height+16)),(500,60),
                                          sdl2.SDL_WINDOW_SHOWN|sdl2.SDL_WINDOW_INPUT_FOCUS|sdl2.SDL_WINDOW_MOUSE_FOCUS)     # opens sdl2 window x2
        surfacex1 = sdl2.ext.Window.get_surface(windowx1)
        surface = surfacex2 = sdl2.ext.Window.get_surface(window)
        sdl2.SDL_SetSurfaceBlendMode(surfacex1 , sdl2.SDL_BLENDMODE_NONE);
        sdl2.SDL_SetSurfaceBlendMode(surfacex2 , sdl2.SDL_BLENDMODE_NONE);
    else:
        scalex2 = False
        window = sdl2.ext.Window(caption,(Width, Height+16),(500,60),
                                          sdl2.SDL_WINDOW_SHOWN|sdl2.SDL_WINDOW_INPUT_FOCUS|sdl2.SDL_WINDOW_MOUSE_FOCUS)     # opens sdl2 window
        surface = surfacex1 = sdl2.ext.Window.get_surface(window)  # ARGB format pixels, use SDL_ConvertSurface if need to convert surface efficiently
        sdl2.SDL_SetSurfaceBlendMode(surfacex1 , sdl2.SDL_BLENDMODE_NONE);

    # pf = sdl2.SDL_GetWindowPixelFormat(window.window)   # https://stackoverflow.com/questions/24576570/updating-window-position-in-pysdl2-help
    # pfname = sdl2.SDL_GetPixelFormatName(pf)
    # print("pixel format name is %s" % pfname)

    cnt = 0
    # window.show()
    cgrid=sdl2.ext.pixels2d(surfacex1)
    sdl2.SDL_RaiseWindow(window.window)
    if scalex2:
        # sdl2.ext.fill(surfacex2, 0)
        sdl2.SDL_BlitScaled(surfacex1,None,surfacex2,None)
    sdl2.ext.Window.refresh(window)

#-----------------------------------------------------------------------------------------------------------
def display_init2():
    global scalex2,Width,Height,rescale,cnt
    global window2,window2x1,surface2,surface2x1,surface2x2,cgrid2,caption2,dispinit2

    if not dispinit and not dispinit2:
        sdl2.ext.init()
    dispinit2 = True
    # clock = sdl2.ext.time.Clock()  # only works in modified version of pysdl2 see https://lukems.github.io/py-sdl2/modules/sdl2ext_time.html
    caption2 = "Gene Life at iteration %d" % cnt
    if (Height <= 512 and rescale):            # Not yet ported well to SDL2
        scalex2 = True
        
        window2x1 = sdl2.ext.Window(caption2,(Width, Height+16),(0,0),sdl2.SDL_WINDOW_HIDDEN)                   # opens sdl2 window x1 never displayed
        window2 = sdl2.ext.Window(caption2,(2*Width, 2*(Height+16)),(500,620),sdl2.SDL_WINDOW_HIDDEN)           # opens sdl2 window 2, initially not displayed
        surface2x1 = sdl2.ext.Window.get_surface(window2x1)
        surface2 = surface2x2 = sdl2.ext.Window.get_surface(window2)
        sdl2.SDL_SetSurfaceBlendMode(surface2x1 , sdl2.SDL_BLENDMODE_NONE);
        sdl2.SDL_SetSurfaceBlendMode(surface2x2 , sdl2.SDL_BLENDMODE_NONE);
    else:
        scalex2 = False
        window2 = sdl2.ext.Window(caption2,(Width, Height+16),(1000,620),
                                          sdl2.SDL_WINDOW_SHOWN|sdl2.SDL_WINDOW_INPUT_FOCUS|sdl2.SDL_WINDOW_MOUSE_FOCUS)     # opens sdl2 window
        surface2 = surface2x1 = sdl2.ext.Window.get_surface(window2)  # ARGB format pixels, use SDL_ConvertSurface if need to convert surface efficiently
        sdl2.SDL_SetSurfaceBlendMode(surface2 , sdl2.SDL_BLENDMODE_NONE);

    pf = sdl2.SDL_GetWindowPixelFormat(window2.window)   # https://stackoverflow.com/questions/24576570/updating-window-position-in-pysdl2-help
    pfname = sdl2.SDL_GetPixelFormatName(pf)
    # print("pixel format name is %s" % pfname)

    cnt = 0
    # window.show()
    cgrid2=sdl2.ext.pixels2d(surface2x1)
    sdl2.SDL_RaiseWindow(window2.window)
    if scalex2:
        sdl2.SDL_BlitScaled(surface2x1,None,surface2x2,None)
    sdl2.ext.Window.refresh(window2)
#-----------------------------------------------------------------------------------------------------------

def display_init2F():           # window with rendering e.g. for text, currently no longer used
    global Width,Height,cnt
    global window2,windowID2,surface2,cgrid2,caption2,dispinit2,renderer2,sdlrenderer2,factory2,font2,textColor2,grect,image2,message2
    
    dispinit2 = True
    caption2 = "Gene Life Magnified 2X"
    window2 = sdl2.ext.Window(caption2,(2*Width, 2*(Height+16)),(800,360))     # opens sdl2 window
    windowID2 = sdl2.SDL_GetWindowID(window2.window)
    window2.show()
    if render2:                                  # see this tutorial https://dev.to/noah11012/using-sdl2-2d-accelerated-renderering-1kcb
        # renderer2 = sdl2.SDL_CreateRenderer( window.window, -1, sdl2.SDL_RENDERER_ACCELERATED | sdl2.SDL_RENDERER_PRESENTVSYNC );
        renderer2 = sdl2.ext.Renderer(window2)
        if (not renderer2):
            print("Failed to create renderer for window")
            print("SDL2 Error: ", sdl2.SDL_GetError())
        sdlrenderer2 = renderer2.sdlrenderer
        sdl2.SDL_SetRenderDrawColor( sdlrenderer2, 0x00, 0x00, 0x00, 0x00 )
        sdl2.SDL_RenderClear(sdlrenderer2)
        
        if sdl2.sdlttf.TTF_Init() == -1:
            print("TTF_Init: %s" % sdl2.sdlttf.TTF_GetError())
        # font = sdl2.sdlttf.TTF_OpenFont( str.encode("lazy.ttf"), 28 )
        font2 = sdl2.sdlttf.TTF_OpenFont( str.encode("Arial.ttf"), 10 )
        textColor2=sdl2.SDL_Color()
        textColor2.r = 255
        textColor2.a = 100
        textColor2.g = 255
        textColor2.b = 255
        
        message2 = sdl2.sdlttf.TTF_RenderText_Solid( font2, str.encode("Genelife frame %5d" % cnt), textColor2 )
        grect.x = Width*2-message2.contents.w-20
        grect.y = 20
        grect.w = message2.contents.w
        grect.h = message2.contents.h
        if not message2:
            print("error rendering message")

        # font_file = sysfont.get_font("freesans")
        # print(font_file)
        # font_manager = sdl2.ext.FontManager(font_file, size=24)
        # factory2 = sdl2.ext.SpriteFactory(sdl2.ext.TEXTURE, renderer=renderer2,fontmanager=font_manager)
        
        factory2 = sdl2.ext.SpriteFactory(sdl2.ext.TEXTURE, renderer=renderer2)
        # surface2 = sdl2.ext.Window.get_surface(window2)
        # surface2 = sdl2.SDL_CreateRGBSurface( 0, Width*2, (Height+16)*2, 32, 0x00FF0000, 0x0000FF00, 0x000000FF, 0xFF000000 ) # format not recog by fill
        surface2=sdl2.SDL_CreateRGBSurfaceWithFormat( 0, Width*2, (Height+16)*2, 32, sdl2.SDL_PIXELFORMAT_ARGB32)
        # print(surface2.contents.format.contents)
        sdl2.ext.fill(surface2.contents, 0)
        image2=factory.from_surface(surface2)
        sdl2.SDL_RenderCopy(sdlrenderer2, image2.texture, None, None)
        message2 = factory2.from_surface(message2)
        sdl2.SDL_RenderCopy(sdlrenderer2, message2.texture, None, grect)
        sdl2.SDL_RenderPresent(sdlrenderer2)

        # texture2 = sdl2.SDL_CreateTextureFromSurface(renderer2, surface2);    # destroy with sdl2.SDL_DestroyTexture(texture)
        # texture2 = sdl2.SDL_CreateTexture( renderer2, sdl2.SDL_PIXELFORMAT_ARGB8888, sdl2.SDL_TEXTUREACCESS_STREAMING, Width, Height+16 );
        # if(not texture2):
        #    print("Failed to convert surface into a texture")
        #    print("SDL2 Error: ", sdl2.SDL_GetError())
        # surface2 = sdl2.SDL_CreateRGBSurface( 0, Width, Height, 32, 0x00FF0000, 0x0000FF00, 0x000000FF, 0xFF000000 )
        # sdl2.SDL_RenderClear(renderer2)
        # sdl2.SDL_RenderCopy(renderer2, texture2, None, None)
        # sdl2.SDL_RenderPresent(renderer2)
        # factory2 = sdl2.ext.SpriteFactory(sdl2.ext.TEXTURE, renderer=renderer2)
    else:
        surface2 = sdl2.ext.Window.get_surface(window2)  # ARGB format pixels, use SDL_ConvertSurface if need to convert surface efficiently
        sdl2.ext.Window.refresh(window2)
    cgrid2=sdl2.ext.pixels2d(surface2)

#-----------------------------------------------------------------------------------------------------------

def update_sim(nrun, ndisp, nskip, niter, nhist, nstat, count=True):
    global gol, cgolg, cgolg2, cgrid, cgrid2, colorfunction, colorfunction2
    global golg
    global log2N
    global runparams
    global cnt,framenr
    global update1,update2

    cnt = cnt+nrun
    if cnt % ndisp == 0 and nrun:  # insert the non-displayed iterations & count species : NB nrun must divide ndisp
        genelife.genelife_update(nskip, nhist, nstat)
        framenr = framenr + nskip
        if(count): genelife.countspecieshash()
    genelife.genelife_update(nrun, nhist, nstat)
    framenr = framenr+nrun
    if update1: colorgrid(colorfunction,cgolg,cgrid,0,0)  # sets  cgrid
    if update2 and (colorfunction2 != -1):
            colorgrid(colorfunction2,cgolg2,cgrid2,1,0)
    return

#-----------------------------------------------------------------------------------------------------------

def pr_params():
    print("runparams[0] = rulemod = ",rulemod)                  
    print("runparams[1] = repscheme = %x"%repscheme)
    print("runparams[2] = selection = ",selection)              
    print("runparams[3] = overwritemask = %x"%overwritemask)
    print("runparams[4] = survivalmask = %x"%survivalmask)      
    print("runparams[7] = birthmask = %x"%birthmask)
    print("runparams[8] = ancselectmask = %x"%ancselectmask)
    print("runparams[5] = colorfunction = ",colorfunction)
    print("runparams[9] = colorfunction2 = ",colorfunction2)
    print("runparams[6] = initfield = ",initfield)                
    print("simparams[0] = nlog2pmut = ",nlog2pmut)                
    print("simparams[1] = initial1density = ",initial1density)          
    print("simparams[2] = initialrdensity = ",initialrdensity)          
    print("simparams[3] = ncoding = ",ncoding)                  
    print("simparams[4] = startgenechoice = ",startgenechoice)
    print("simparams[5] = ranseed = ",ranseed)
    
#-----------------------------------------------------------------------------------------------------------

def set_params():
    global rulemod,repscheme,survivalmask,birthmask,overwritemask,ancselectmask,selection,ncoding
    global startgenechoice,initialrdensity,initial1density,nlog2pmut,initfield
    global colorfunction
    global runparams, simparams
    
    runparams[0] = rulemod                   # 0,1 whether to allow GoL rule modifications
                                                 # with rulemod 1 2-live-nb birth, 3-live-nb non-birth & non-survival possible
    runparams[1] = repscheme                 # 0-7 20 control bits for repscheme 8-15 <12 different control bits for repscheme 16-19 nplanes
    runparams[2] = selection                 # 0-7 fitness for 2 live neighbor rule, 8-15 LUT symmetry model and gene coding modes 16-19 16 planes 20-23 64 planes 24 matching
    runparams[3] = overwritemask             # mask of bits to overwrite sum or lut entry with birth instead of survival
    runparams[4] = survivalmask              # 8-32 bit survival mask for allowing genes to modify LUTs
    runparams[7] = birthmask                 # 8-32 bit birth mask for allowing genes to modify LUTs
    runparams[8] = ancselectmask             # 8-32 bit ancesor selection mask to allow genetic selection to determine ancesotr for LUT rule
    runparams[5] = colorfunction             # color function; 0(hash), 1-3 (functional), 2 nongulstate or color gol planes, 3 notgolrul yellow
                                             # 4 activities 5 popln 6 genealogy steps 7 genealogy temporal 8 glider detection
                                             # 9 connected component labels and novelty (n) 10 connected component activities
                                             # 11 genealogy depth 12 genetic glider detection
    runparams[9] = colorfunction2            # -1 for same as clorfunction, otherwise values as in colorfunction
    runparams[6] = initfield                 # 0 full field random or start depending on initialrdensity, 1 init via 32x32 genepat.dat, n>1 init via nxn rand array
    simparams[0] = nlog2pmut                 # log2 gene mutation probability (0 or >56 means no mutation)
    simparams[1] = initial1density           # initial 1 density in GOL state
                                                 # 16384 = nearest to half of guaranteed C rand max value 32767 = 2**15 - 1
    simparams[2] = initialrdensity           # initial density of random genes
    simparams[3] = ncoding                   # for selection 10, non zero value means grow plane community from 0
                                                 # otherwise (selection<10) no of bits used to encode valid connection functions 1-16
                                                 # for selection==8, lut, ncoding 1,2,3 bits per lut entry : 0 implies 3.
    simparams[4] = startgenechoice           # initialize genes to startgene number 0-8 : 8 is random choice of 0-7
    simparams[5] = ranseed                   # initialize ranseed
    
    pr_params()
    
#-----------------------------------------------------------------------------------------------------------
#
# activity: run for N generations and plot the quantiles vs time (semilog)
# assumes runparams, simparams are set

def all_activity(N=1000,nquant=10,acttype="live", # all = all genes, live = live genes only,  quad = quad patterns, small = small patterns
             maxnum=100000,init=False):        # for quads maxnum needs to be bumped up (e.g. > 2*10^6), which slows things down
    if acttype=="live":
        doact = genelife.get_activities
        doquad = genelife.get_quad_activities
        dosmall = genelife.get_small_activities
    elif acttype=='all':
        doact = genelife.get_all_activities
        doquad = genelife.get_all_quad_activities
        dosmall = genelife.get_all_small_activities
    else:
        print("unknown acttype:  ",acttype)

    activities = np.zeros(maxnum,np.int32)
    actsmall = np.zeros(65536,np.int32)
    actquad = np.zeros(maxnum,np.int32)

    actdata = np.zeros(maxnum,np.uint64)
    smalldata = np.zeros(65536,np.uint64)
    quaddata = np.zeros(maxnum,np.uint64)


    qactivity = [[None]*nquant for _ in range(N)]
    qsmall = [[None]*nquant for _ in range(N)]
    qquad = [[None]*nquant for _ in range(N)]
    qq = [x/nquant for x in range(nquant)] # same number of quantiles for all
    nspecies = [None]*N
    nquad = [None]*N
    nsmall = [None]*N

    if init:
        genelife.initialize_planes(npoffsets)
        genelife.initialize(runparams,simparams)
        genelife.set_seed(ranseed)
    genelife.genelife_update(1,0,0)
    nspecies[0]=doact(actdata,activities)
    nquad[0] = doquad(quaddata,actquad)
    nsmall[0] = dosmall(smalldata,actsmall)

    for j in range(1,N):
        # stash the stats
        ac = [activities[i] for i in range(len(activities)) if activities[i]>1]
        qactivity[j] = np.quantile(ac,qq)
        ac = [actquad[i] for i in range(len(activities)) if activities[i]>1]
        qquad[j] = np.quantile(ac,qq)
        ac = [actsmall[i] for i in range(len(activities)) if activities[i]>1]
        qsmall[j] = np.quantile(ac,qq)
        # iterate one time step
        genelife.genelife_update(1,0,0)
        # compute the stats
        nspecies[j]=doact(actdata,activities)
        nsmall[j] = dosmall(smalldata,actsmall)
        nquad[j] = doquad(quaddata,actquad)

    rtn = {}
    rtn['actquantiles'] = qactivity
    rtn['smallquantiles'] = qsmall
    rtn['quadquantiles'] = qquad
    rtn['nspecies'] = nspecies
    rtn['nsmall'] = nsmall
    rtn['nquad'] = nquad
    return(rtn)

def activity(N=1000,nquant=10,acttype="live", # all = all genes, live = live genes only,  quad = quad patterns, small = small patterns
             maxnum=100000,init=True):        # for quads maxnum needs to be bumped up (e.g. > 2*10^6), which slows things down
    if acttype=="live":
        doact = genelife.get_activities
    elif acttype=="all":
        doact = genelife.get_all_activities
    elif acttype=="quad":
        doact = genelife.get_quad_activities
    elif acttype=="small":
        doact = genelife.get_small_activities
    else:
        print("unknown acttype:  ",acttype)

    activities = np.zeros(maxnum,np.int32)

    data=np.zeros(maxnum,np.uint64)
    qqq = [[None]*nquant for _ in range(N)]
    qq = [x/nquant for x in range(nquant)]
    nspecies = [None]*N

    if init:
        genelife.initialize_planes(npoffsets)
        genelife.initialize(runparams,simparams)
        genelife.set_seed(ranseed)
    genelife.genelife_update(1,0,0)
    nspecies[0]=doact(data,activities)

    for j in range(1,N):
        ac = [activities[i] for i in range(len(activities)) if activities[i]>1]
        qqq[j] = np.quantile(ac,qq)
        genelife.genelife_update(1,0,0)
        nspecies[j]=doact(data,activities)

    rtn = {}
    rtn['quantiles']=qqq
    rtn['nspecies'] = nspecies
    return(rtn)

    # for j in range(10):
    #     foo = [qqq[i][j] for i in range(N)]
    #     plt.semilogy(foo)
    # plt.show()
    # plt.plot(nspecies);

def plot_activity(dat):
    qqq = dat['quantiles']
    pop = dat['nspecies']
    N = len(qqq)
    nq = len(qqq[0])
    for j in range(nq):
        foo = [qqq[i][j] for i in range(N)]
        plt.semilogy(foo)
    plt.title("Activity quantiles")
    plt.xlabel('time')
    plt.ylabel('activity')
    plt.show()
    plt.plot(pop)
    plt.title('Number of species')
    plt.xlabel('time')

#-----------------------------------------------------------------------------------------------------------

def show0(count=True):
# display initial population and count species
    global framenr
    global surface, surfacex1, window, scalex2, caption, dispinit, colorfunction, colorfunction2
    global surface2, surface2x1, window2, caption2, dispinit2
    # global repscheme,survivalmask,overwritemask,ancselectmask,selection
    global cancol,cgolg,cgolg2,cgrid,cgrid2

    if not dispinit:
        display_init()
    if not dispinit2:
        display_init2()
        # display_init2F()
    caption = "Gene Life at iteration %d" % framenr
    set_caption(window,caption)

    cancol=init_buttons(surfacex1)                           # initialize parameter buttons

    colorgrid(colorfunction,cgolg,cgrid,0,0)
    
    if colorfunction2 != -1:
        colorgrid(colorfunction2,cgolg2,cgrid2,1,0)
        # sdl2.ext.fill(surface2, 0)
        if sdl2.SDL_BlitScaled(surfacex1,grect1,surface2x1,grect1): print("BlitScaled failed")
        if sdl2.SDL_BlitScaled(surface2x1,None,surface2,None): print("BlitScaled failed")

    else:
        if sdl2.SDL_BlitScaled(surface,None,surface2,None): print("BlitScaled failed")

    sdl2.ext.Window.refresh(window)
    sdl2.ext.Window.refresh(window2)
    
    if(count):
        genelife.countspecieshash()

#-----------------------------------------------------------------------------------------------------------

def step(count=True):
    """single step and update display and species counts"""
    global framenr
    #global gol,golg,golgstats
    global surface, surfacex1, window, scalex2, dispinit, grect1
    global surface2, surface2x1, window2, caption2, dispinit2
    
    if not dispinit:
        display_init()
    if not dispinit2:
        display_init2()
        # display_init2F()

    update_sim(1, 1, 0, 1, 0, 0, count)
    caption = "Gene Life at iteration %d" % framenr
    
    set_caption(window,caption)
    set_caption(window2,caption)

    if colorfunction2 != -1:
        if sdl2.SDL_BlitScaled(surfacex1,grect1,surface2x1,grect1): print("BlitScaled failed")
        if sdl2.SDL_BlitScaled(surface2x1,None,surface2,None): print("BlitScaled failed")
    else:
        if sdl2.SDL_BlitScaled(surface,None,surface2,None): print("BlitScaled failed")

    sdl2.ext.Window.refresh(window)
    sdl2.ext.Window.refresh(window2)
    
    if (count):
        genelife.countspecieshash()
#-----------------------------------------------------------------------------------------------------------
def construct_caption(colorfunction1or2,pixeldat,buttonhelp,win):
    """ construct window caption
    """
    global colorfunction,framenr,nspecies,selection
    global quadrants,ymax,ymaxq,offdx,offdy,offdt,ncomponents,genealogycoldepth,ancestortype
    
    selectiontext0007 = ["largest value","most ones","scissors-well-stone-paper","not well ordered","two target","predator prey","cooperative","neutral"];
    selectiontext0815 = ["sum fixed","sum variable","edge fixed","edge variable","canonical fixed","canonical variable","2D sym fixed","2D sym variable"];

    caption = "Gene Life at step %d coloring %d nspecies %d " % (framenr,colorfunction1or2,nspecies)
    if(win == 2): caption = "2nd " + caption
    if selection < 8:
        caption = caption + "pairwise selection " + selectiontext0007[selection] + " " + buttonhelp
    elif selection<16:
        caption = caption + "LUT encoding " + selectiontext0815[selection-8] + " " + buttonhelp

    if quadrants >= 0:
        paramdat = "repscheme %06x surv. %01x overw. %01x ncoding %06x" % (repscheme,survivalmask,overwritemask,ncoding)
        caption = caption + ("q%1d " % quadrants) + paramdat
    if colorfunction1or2 == 4: caption = caption + ("ymax %d " % ymax)
    elif colorfunction1or2 == 6 or colorfunction1or2 == 7 or colorfunction1or2 == 11:
        if colorfunction1or2 == 11: caption = caption + ("gcdepth %d " % genealogycoldepth)
        if   ancestortype == 0: caption = caption + "anc first "
        elif ancestortype == 1: caption = caption + "anc clonal "
        elif ancestortype == 2: caption = caption + "anc first & clonal"
    elif colorfunction1or2 == 8: caption = caption + ("offsets (%d,%d,%d) " % (offdx,offdy,offdt))
    elif colorfunction1or2 == 9:
        ncomponents=genelife.get_ncomponents()
        caption = caption + ("ncomponents %d " % (ncomponents))
    elif colorfunction1or2 == 10: caption = caption + ("ymaxq %d " % ymaxq)

    if pixeldat: caption = caption + pixeldat
    return caption
#-----------------------------------------------------------------------------------------------------------
# infinite loop of display updates
# cmd click in graphics window to stop, click for pixel details or quadrant selection,
# alt-click or arrow keys for recolor,
# +/- keys reserved for activity ymax : actually the crossover value in N* act/(ymax+act)
# keys lower case - decrement, upper case - increment, alt - input value: y,Y ymax q,Q quadrant
# misc. keys save image
def run(nrun, ndisp, nskip, niter, nhist, nstat, count=True, maxsteps=100000):
    global mstime,framenr,framerate
    global window, surface, surfacex1, surfacex2, scalex2, caption, dispinit, update1, grect1
    global window2, surface2, surface2x1, surface2x2, caption2, dispinit2, grect, render2, update2, renderer2, factory2, image2, message2, font2, textColor2, windowID2
    global N,nNhist
    global gol,golg,golb,golr,golgstats
    global connlabel,connlen,ncomponents,nnovelcells
    global colorfunction,colorfunction2,gcolor,nfrstep,genealogycoldepth,ancestortype
    global nspecies,ymax,ymaxq,oldymax,oldymaxq,nbhist,nNhist
    global updatesenabled
    global rulemod,repscheme,survivalmask,birthmask,overwritemask,ancselectmask,selection,ncoding
    global savecnt
    global cancol
    global Height,Width
    global randominflux,vscrolling,noveltyfilter,activity_size_colormode,info_transfer_h,activityfnlut
    global gogo,pause,mouseclicked,mouseclicked2,pixeldat,paramdat
    global maxPlane,offdx,offdy,offdt,quadrants,displayoneplane
    global parhelp

    mstime = sdl2.timer.SDL_GetTicks()
    actsizecoltxt= [" color from hashkey"," color log2 of size"," color # pixels"," color sqrt # pixels"]
    
    buttonhelp0007 =    ["0. selective birth for 3-live-nbs ","1. selective birth for 2-live-nbs ",
                         "2. canonical 0 position vs difft  ","3. bypass selection for 2-live-nbs ",
                         "4. enforce birth for 3-live-nbs ","5. enforce birth for 2-live-nbs ",
                         "6. 2nd nb genetic modulation ","7. 1st nb genetic masking ",
                         "8. enforce GoL if non GoL rule ","9. enforce GoL last change by non GoL ",
                         "10. allow 2-nb birth for canonical config 0 ","11. allow 2-nb birth for canonical config 1 ",
                         "12. allow 2-nb birth for canonical config 2 ","13. allow 2-nb birth for canonical config 3 ",
                         "14. Parentdies is forced ","15. Dummy, not yet used ",
                         "16. Survival for 3-live-nbs ","17. Survival for 2-live-nbs ",
                         "18. Gene overwrite for 3-live-nbs ","19. Gene overwrite for 2-live-nbs "
                         ]
    buttonhelp0815 =    ["0. survival gene central/nbs ","1. OR/AND of LUTs of nbs ",
                         "2. canonical 0 position vs difft  ","3. Parentdies is forced ",
                         "4. bit 0 of selection mode ","5. bit 1 of selection mode ",
                         "6. bit 2 of selection mode ","7. bit 3 of selection mode (golr) ",
                         "8. bit 0 of disambig mode 0-7 ","9. bit 1 of disambig mode 0-7 ",
                         "10. bit 2 of disambig mode 0-7 ","11. random ancestor choice "
                         ]
    buttonhelp = ""
    if scalex2:
        sc = 1
    else:
        sc = 2

    if not dispinit:
        display_init()
    if not dispinit2:
        display_init2()
        # display_init2F()
    # if scalex2:
    cancol=init_buttons(surfacex1)
    
    surviveover = np.array([survivalmask,birthmask,overwritemask],dtype=np.uint32)
    gogo = True
    pixeldat = ""
    paramdat = ""
    mouseclicked = False
    mouseclicked2 = False
    pause = 0
    ymax = ymaxq = 10000
    maxPlane = 4
    offdx = offdy=offdt=0
    quadrants = -1
    oldymax = genelife.setget_act_ymax(ymax)
    oldymaxq = genelife.setget_act_ymaxq(ymaxq)
    nfrstep = 0

    gcolor=0
    
    sdl2.SDL_ShowWindow(window.window)
    sdl2.SDL_RaiseWindow(window.window)
    if colorfunction2 != -1:
        sdl2.SDL_ShowWindow(window2.window)
        sdl2.SDL_RaiseWindow(window2.window)
    else:
        sdl2.SDL_HideWindow(window2.window)
    event = sdl2.SDL_Event()
    while (gogo):
        for event in sdl2.ext.get_events():
        # if (sdl2.SDL_PollEvent(ctypes.byref(event))):
            if event.type == sdl2.SDL_QUIT:
                mouseclicked = False
                gogo = False
                # sdl2.ext.quit()                  # check that quitting SDL here is OK
            if event.type == sdl2.SDL_MOUSEBUTTONDOWN and update1:
                if   (sdl2.SDL_GetModState() & sdl2.KMOD_ALT) | (event.button.button == sdl2.SDL_BUTTON_MIDDLE): # quit event loop on middle mouse button (option-click)
                    mouseclicked = False
                    gogo = False
                    # sdl2.ext.quit()              # check that quitting SDL here is OK
                elif event.button.button == sdl2.SDL_BUTTON_LEFT:          # get mouse coords on mouse event
                    mouseclicked = True
                    if scalex2:                    # or (event.window.windowID == windowID2):
                        x = (int) (event.button.x//2)
                        y = (int) (event.button.y//2)
                    else:
                        x = event.button.x
                        y = event.button.y
                    if y >= N:
                        k=x>>(log2N-6)
                        if selection<8:
                            if k<20:
                                if k<16:
                                    repscheme = repscheme ^ (1<<k)
                                    bit = (repscheme>>k)&0x1
                                    print(("step %d repscheme changed to %x" % (framenr,repscheme)))
                                elif k<18:
                                    survivalmask = survivalmask ^ (1<<(k-16))
                                    bit = (survivalmask>>(k-16))&0x1
                                    print(("step %d survivalmask changed to %x" % (framenr,survivalmask)))
                                else:
                                    overwritemask = overwritemask ^ (1<<(k-18))
                                    bit = (overwritemask>>(k-18))&0x1
                                    print(("step %d overwritemask changed to %x" % (framenr,overwritemask)))
                                draw_rect(surfacex1,cancol[0][k]*(1+bit),[k<<(log2N-6),Height+6,3*sc,3*sc])
                                surviveover[0],surviveover[1]= survivalmask,overwritemask      # 2nd elt only picked up in C as overwrite for selection<8
                                genelife.set_surviveover64(surviveover)
                                genelife.set_repscheme(repscheme)
                        elif selection < 10:
                            if k<24:
                                if k<8:
                                    survivalmask = survivalmask ^ (1<<k)
                                    print(("step %d survivalmask changed to %x" % (framenr,survivalmask)))
                                    draw_rect(surfacex1,cancol[1][k]*(1+((survivalmask>>k)&0x1)),[k<<(log2N-6),Height+6,3*sc,3*sc])
                                elif k<16:
                                    birthmask = birthmask ^ (1<<(k-8))
                                    print(("step %d birthmask changed to %x" % (framenr,birthmask)))
                                    draw_rect(surfacex1,cancol[1][k]*(1+((birthmask>>(k-8))&0x1)),[k<<(log2N-6),Height+6,3*sc,3*sc])
                                elif k<24:
                                    overwritemask = overwritemask ^ (1<<(k-16))
                                    print(("step %d overwritemask changed to %x" % (framenr,overwritemask)))
                                    draw_rect(surfacex1,cancol[1][k]*(1+((overwritemask>>(k-16))&0x1)),[k<<(log2N-6),Height+6,3*sc,3*sc])
                                surviveover[0],surviveover[1],surviveover[2]= survivalmask,birthmask,overwritemask
                                genelife.set_surviveover64(surviveover)
                            elif k<24+12:
                                repscheme = repscheme ^ (1<<(k-24))
                                print(("step %d repscheme changed to %x" % (framenr,repscheme)))
                                draw_rect(surfacex1,cancol[1][k]*(1+((repscheme>>(k-24))&0x1)),[k<<(log2N-6),Height+6,3*sc,3*sc])
                                genelife.set_repscheme(repscheme)
                        elif selection < 12:
                            if k<12 and y>=N+12:
                                repscheme = repscheme ^ (1<<k)
                                print(("step %d repscheme changed to %x" % (framenr,repscheme)))
                                draw_rect(surfacex1,cancol[2][54+k]*(1+((repscheme>>k)&0x1)),[k<<(log2N-6),Height+8+3*sc,3*sc,3*sc])
                                genelife.set_repscheme(repscheme)
                            elif k<54:
                                if k<23:
                                    survivalmask = survivalmask ^ (1<<k)
                                    print(("step %d survivalmask changed to %x" % (framenr,survivalmask)))
                                    draw_rect(surfacex1,cancol[2][k]*(1+((survivalmask>>k)&0x1)),[k<<(log2N-6),Height+6,3*sc,3*sc])
                                elif k<46:
                                    birthmask = birthmask ^ (1<<(k-23))
                                    print(("step %d birthmask changed to %x" % (framenr,birthmask)))
                                    draw_rect(surfacex1,cancol[2][k]*(1+((birthmask>>(k-23))&0x1)),[k<<(log2N-6),Height+6,3*sc,3*sc])
                                else:
                                    overwritemask = overwritemask ^ (1<<(k-46))
                                    print(("step %d overwritemask changed to %x" % (framenr,overwritemask)))
                                    draw_rect(surfacex1,cancol[2][k]*(1+((overwritemask>>(k-46))&0x1)),[k<<(log2N-6),Height+6,3*sc,3*sc])
                                surviveover[0],surviveover[1],surviveover[2]= survivalmask,birthmask,overwritemask
                                genelife.set_surviveover64(surviveover)
                        elif selection<14:
                            if k<64 and y<N+12:
                                if k<32:
                                    survivalmask = survivalmask ^ (1<<k)
                                    print(("step %d survivalmask changed to %x" % (framenr,survivalmask)))
                                    draw_rect(surfacex1,cancol[3][k]*(1+((survivalmask>>k)&0x1)),[k<<(log2N-6),Height+6,3*sc,3*sc])
                                else:
                                    birthmask = birthmask ^ (1<<(k-32))
                                    print(("step %d birthmask changed to %x" % (framenr,birthmask)))
                                    draw_rect(surfacex1,cancol[3][k]*(1+((birthmask>>(k-32))&0x1)),[k<<(log2N-6),Height+6,3*sc,3*sc])
                                surviveover[0],surviveover[1],surviveover[2]= survivalmask,birthmask,overwritemask
                                genelife.set_surviveover64(surviveover)
                            elif k<12:
                                repscheme = repscheme ^ (1<<k)
                                print(("step %d repscheme changed to %x" % (framenr,repscheme)))
                                draw_rect(surfacex1,cancol[3][k+72]*(1+((repscheme>>k)&0x1)),[k<<(log2N-6),Height+8+3*sc,3*sc,3*sc])
                                genelife.set_repscheme(repscheme)
                            elif k>=32 and k<40:
                                overwritemask = overwritemask ^ (1<<(k-32))
                                print(("step %d overwritemask changed to %x" % (framenr,overwritemask)))
                                draw_rect(surfacex1,cancol[3][k+32]*(1+((overwritemask>>(k-32))&0x1)),[k<<(log2N-6),Height+8+3*sc,3*sc,3*sc])
                                surviveover[0],surviveover[1],surviveover[2]= survivalmask,birthmask,overwritemask
                                genelife.set_surviveover64(surviveover)
                        elif selection<16:
                            if k<64 and y<N+12:
                                if k<32:
                                    survivalmask = survivalmask ^ (1<<k)
                                    print(("step %d survivalmask changed to %x" % (framenr,survivalmask)))
                                    draw_rect(surfacex1,cancol[4][k]*(1+((survivalmask>>k)&0x1)),[k<<(log2N-6),Height+6,3*sc,3*sc])
                                else:
                                    birthmask = birthmask ^ (1<<(k-32))
                                    print(("step %d birthmask changed to %x" % (framenr,birthmask)))
                                    draw_rect(surfacex1,cancol[4][k]*(1+((birthmask>>(k-32))&0x1)),[k<<(log2N-6),Height+6,3*sc,3*sc])
                                surviveover[0],surviveover[1],surviveover[2]= survivalmask,birthmask,overwritemask
                                genelife.set_surviveover64(surviveover)
                            elif k<12:
                                repscheme = repscheme ^ (1<<k)
                                print(("step %d repscheme changed to %x" % (framenr,repscheme)))
                                draw_rect(surfacex1,cancol[4][k+72]*(1+((repscheme>>k)&0x1)),[k<<(log2N-6),Height+8+3*sc,3*sc,3*sc])
                                genelife.set_repscheme(repscheme)
                            elif k>=32 and k<40:
                                overwritemask = overwritemask ^ (1<<(k-32))
                                print(("step %d overwritemask changed to %x" % (framenr,overwritemask)))
                                draw_rect(surfacex1,cancol[4][k+32]*(1+((overwritemask>>(k-32))&0x1)),[k<<(log2N-6),Height+8+3*sc,3*sc,3*sc])
                                surviveover[0],surviveover[1],surviveover[2]= survivalmask,birthmask,overwritemask
                                genelife.set_surviveover64(surviveover)
                    else: # y<N
                        if colorfunction < 4 or colorfunction == 8 or colorfunction >10 :
                            genelife.get_curgol(gol)    # get current gol,golg,golgstats arrays
                            genelife.get_curgolg(golg)
                            if colorfunction > 10:
                                genelife.get_curgolbr(golb,golr)
                            else:
                                genelife.get_curgolgstats(golgstats)
                            if quadrants >= 0 and selection<8:   # set the two bits in repscheme corresponding to quadrant
                                repscheme=genelife.set_repscheme_bits(quadrants,x,y,surviveover)
                                survivalmask  = surviveover[0]
                                overwritemask = surviveover[1]
                                repscheme=packrepscheme(repscheme,survivalmask,overwritemask)
                                print(("step %d repscheme changed to %x" % (framenr,repscheme)))
                                quadrants = -1
                                pixeldat = ""
                            else:
                                if colorfunction < 10:
                                    genelife.get_curgolgstats(golgstats)
                                    if colorfunction == 2:
                                        pixeldat = "(%d,%d) gol %01x gene %016x s %03d" % (x,y,gol[x+y*N],golg[x+y*N],(golgstats[x+y*N] % 8)+1)
                                    else:
                                        pixeldat = "(%d,%d) gol %01x gene %016x status %016x" % (x,y,gol[x+y*N],golg[x+y*N],golgstats[x+y*N])
                                elif colorfunction == 11:
                                    pixeldat = "(%d,%d) gol %01x cloneid %016x gene %016x" % (x,y,gol[x+y*N],golb[x+y*N],golg[x+y*N])
                                elif colorfunction == 12:
                                    pixeldat = "(%d,%d) gol %01x golr %016x golg %016x" % (x,y,gol[x+y*N],golr[x+y*N],golg[x+y*N])
                                print(("step %d pixel data %s" % (framenr,pixeldat)))
                                if selection == 8:                              # color rule table rectangles at base by rule derived from gene at current pixel
                                    for k in range(16):
                                        draw_rect(surfacex1,cancol[1][k]*(1+(np.right_shift(np.uint64(golg[x+y*N]),np.uint64(k))&np.uint64(0x1))),[k<<(log2N-6),Height+6,3*sc,3*sc])
                                elif selection == 10:
                                    for k in range(46):
                                        draw_rect(surfacex1,cancol[2][k]*(1+(np.right_shift(np.uint64(golg[x+y*N]),np.uint64(k))&np.uint64(0x1))),[k<<(log2N-6),Height+6,3*sc,3*sc])
                                elif selection == 12:
                                    for k in range(64):
                                        draw_rect(surfacex1,cancol[3][k]*(1+(np.right_shift(np.uint64(golg[x+y*N]),np.uint64(k))&np.uint64(0x1))),[k<<(log2N-6),Height+6,3*sc,3*sc])
                                elif selection == 14:
                                    for k in range(64):
                                        draw_rect(surfacex1,cancol[4][k]*(1+(np.right_shift(np.uint64(golg[x+y*N]),np.uint64(k))&np.uint64(0x1))),[k<<(log2N-6),Height+6,3*sc,3*sc])
                        elif colorfunction == 4:
                            genelife.get_acttrace(golg)
                            pixeldat = "(%d,%d) gene %016x" % (x,y,golg[x+y*N])
                        elif colorfunction == 5:
                            genelife.get_poptrace(golg)
                            pixeldat = "(%d,%d) gene %016x" % (x,y,golg[x+y*N])
                        elif colorfunction <= 7:
                            genelife.get_genealogytrace(golg)
                            pixeldat = "(%d,%d) gene %016x" % (x,y,golg[x+y*N])
                            genelife.set_selectedgene(golg[x+y*N])
                            print(("step %d pixel data %s" % (framenr,pixeldat)))
                        elif colorfunction == 9:
                            ncomponents=genelife.get_connected_comps(connlabel,connlen,x,y)
                            colorgrid(colorfunction,cgolg,cgrid,0,0)
                            pixeldat = "(%d,%d) label %4d nrconn %d" % (x,y,connlabel[y*N+x],connlen[connlabel[y*N+x]])
                        elif colorfunction == 10:
                            ncomponents=genelife.get_connected_comps(connlabel,connlen,x,y)
                            colorgrid(colorfunction,cgolg,cgrid,0,0)
                            pixeldat = "(%d,%d)" % (x,y)
                elif event.button.button ==  sdl2.SDL_BUTTON_RIGHT:          # info on button or selection<8 model - right mouse button (-click)
                    # print("right mouse button pressed")
                    if scalex2 or event.window.windowID == windowID2:
                        x = (int) (event.button.x//2)
                        y = (int) (event.button.y//2)
                    else:
                        x = event.button.x
                        y = event.button.y
                    if selection<16:
                        if y >= N:
                            k=x>>(log2N-6)
                            if selection < 8:
                                if k<20:
                                    buttonhelp = buttonhelp0007[k]
                            elif selection < 10:
                                if k >= 24 and k < 24 + 12:
                                    buttonhelp = buttonhelp0815[k-24]
                            elif selection < 16:
                                if y >= N+8 and k<12:
                                        buttonhelp = buttonhelp0815[k]
                    mouseclicked2 = True
            elif event.type==sdl2.SDL_MOUSEBUTTONUP:
                mouseclicked = False
                mouseclicked2 = False
                buttonhelp = ""
                if selection == 8:                                  # reset mask control buttons to survivalmask and birthmask control colours
                    for k in range(16):
                        if k<8: draw_rect(surfacex1,cancol[1][k]*(1+((survivalmask>>k)&0x1)),[k<<(log2N-6),Height+6,3*sc,3*sc])
                        else: draw_rect(surfacex1,cancol[1][k]*(1+((birthmask>>(k-8))&0x1)),[k<<(log2N-6),Height+6,3*sc,3*sc])
                elif selection == 10:
                    for k in range(46):
                        if k<23: draw_rect(surfacex1,cancol[2][k]*(1+((survivalmask>>k)&0x1)),[k<<(log2N-6),Height+6,3*sc,3*sc])
                        else: draw_rect(surfacex1,cancol[2][k]*(1+((birthmask>>(k-23))&0x1)),[k<<(log2N-6),Height+6,3*sc,3*sc])
                elif selection==12:
                    for k in range(64):
                        if k<32: draw_rect(surfacex1,cancol[3][k]*(1+((survivalmask>>k)&0x1)),[k<<(log2N-6),Height+6,3*sc,3*sc])
                        else: draw_rect(surfacex1,cancol[3][k]*(1+((birthmask>>(k-32))&0x1)),[k<<(log2N-6),Height+6,3*sc,3*sc])
                elif selection==14:
                    for k in range(64):
                        if k<32: draw_rect(surfacex1,cancol[4][k]*(1+((survivalmask>>k)&0x1)),[k<<(log2N-6),Height+6,3*sc,3*sc])
                        else: draw_rect(surfacex1,cancol[4][k]*(1+((birthmask>>(k-32))&0x1)),[k<<(log2N-6),Height+6,3*sc,3*sc])
                if colorfunction==9 or colorfunction==10:
                    ncomponents=genelife.get_connected_comps(connlabel,connlen,-1,-1)
                    colorgrid(colorfunction,cgolg,cgrid,0,0)
                pixeldat = ""
            elif event.type==sdl2.SDL_MOUSEMOTION:
                if mouseclicked:
                    if scalex2 or event.window.windowID == windowID2:
                        x = (int) (event.motion.x//2)
                        y = (int) (event.motion.y//2)
                    else:
                        x = event.motion.x
                        y = event.motion.y
                    if x < N and y < N:
                        if colorfunction < 4 or colorfunction == 8 or colorfunction > 10:
                            genelife.get_curgol(gol)    # get current gol,golg,golgstats arrays
                            genelife.get_curgolg(golg)
                            if colorfunction < 10:
                                genelife.get_curgolgstats(golgstats)
                                if colorfunction == 2:
                                    pixeldat = "(%d,%d) gol %01x gene %016x s %03d" % (x,y,gol[x+y*N],golg[x+y*N],(golgstats[x+y*N] % 8)+1)
                                else:
                                    pixeldat = "(%d,%d) gol %01x gene %016x status %016x" % (x,y,gol[x+y*N],golg[x+y*N],golgstats[x+y*N])
                            elif colorfunction == 11:
                                    genelife.get_curgolbr(golb,golr)
                                    pixeldat = "(%d,%d) gol %01x cloneid %016x gene %016x" % (x,y,gol[x+y*N],golb[x+y*N],golg[x+y*N])
                            elif colorfunction == 12:
                                    genelife.get_curgolbr(golb,golr)
                                    pixeldat = "(%d,%d) gol %01x golr %016x golg %016x" % (x,y,gol[x+y*N],golr[x+y*N],golg[x+y*N])
                            if selection == 8:
                                for k in range(16):
                                    draw_rect(surfacex1,cancol[1][k]*(1+(np.right_shift(np.uint64(golg[x+y*N]),np.uint64(k))&np.uint64(0x1))),[k<<(log2N-6),Height+6,3*sc,3*sc])
                            elif selection ==10:
                                for k in range(46):
                                    draw_rect(surfacex1,cancol[2][k]*(1+(np.right_shift(np.uint64(golg[x+y*N]),np.uint64(k))&np.uint64(0x1))),[k<<(log2N-6),Height+6,3*sc,3*sc])
                            elif selection ==12:
                                for k in range(64):
                                    draw_rect(surfacex1,cancol[3][k]*(1+(np.right_shift(np.uint64(golg[x+y*N]),np.uint64(k))&np.uint64(0x1))),[k<<(log2N-6),Height+6,3*sc,3*sc])
                            elif selection ==14:
                                for k in range(64):
                                    draw_rect(surfacex1,cancol[4][k]*(1+(np.right_shift(np.uint64(golg[x+y*N]),np.uint64(k))&np.uint64(0x1))),[k<<(log2N-6),Height+6,3*sc,3*sc])
                        elif colorfunction == 4:
                            genelife.get_acttrace(golg)
                            pixeldat = "(%d,%d) gene %016x" % (x,y,golg[x+y*N])
                        elif colorfunction == 5:
                            genelife.get_poptrace(golg)
                            pixeldat = "(%d,%d) gene %016x" % (x,y,golg[x+y*N])
                        elif colorfunction <= 7:
                            genelife.get_genealogytrace(golg)
                            pixeldat = "(%d,%d) gene %016x" % (x,y,golg[x+y*N])
                            genelife.set_selectedgene(golg[x+y*N])
                        elif colorfunction == 9:
                            ncomponents=genelife.get_connected_comps(connlabel,connlen,x,y)
                            colorgrid(colorfunction,cgolg,cgrid,0,0)
                            pixeldat = "(%d,%d) label %4d nr.conn %d" % (x,y,connlabel[y*N+x],connlen[connlabel[y*N+x]])
                        elif colorfunction == 10:
                            ncomponents=genelife.get_connected_comps(connlabel,connlen,x,y)
                            colorgrid(colorfunction,cgolg,cgrid,0,0)
                            pixeldat = "(%d,%d)" % (x,y)
                elif mouseclicked2:
                    if colorfunction == 2:
                        if scalex2:                 # or event.window.windowID == windowID2:
                            x = (int) (event.motion.x//2)
                            y = (int) (event.motion.y//2)
                        else:
                            x = event.motion.x
                            y = event.motion.y
                    if selection<16:
                        if scalex2 or event.window.windowID == windowID2:
                            x = (int) (event.motion.x//2)
                            y = (int) (event.motion.y//2)
                        else:
                            x = event.motion.x
                            y = event.motion.y
                        if y >= N:
                            k=x>>(log2N-6)
                            if selection < 8:
                                if k<20:
                                    buttonhelp = buttonhelp0007[k]
                            elif selection < 10:
                                if k >= 24 and k < 24 + 12:
                                    buttonhelp = buttonhelp0815[k-24]
                            elif selection < 16:
                                if y >= N+8 and k<12:
                                        buttonhelp = buttonhelp0815[k]
            elif event.type == sdl2.SDL_KEYDOWN:
                # keystatus = sdl2.SDL_GetKeyboardState(None) # keystatus should also reveal if pressed if event structure doesn't work
                # if keystatus[sdl2.SDL_SCANCODE_H]:
                # print("key pressed")
                if event.key.keysym.scancode == sdl2.SDL_SCANCODE_H:                   # alternatively use:
                    if   sdl2.SDL_GetModState() &  sdl2.KMOD_SHIFT:
                        rulemod = rulemod ^ 2   # horizon mode with GoL in upper half toggled on/off
                        print(("step %d rulemod changed to %x (Horizon mode)" % (framenr,rulemod)))
                        genelife.set_rulemod(rulemod)
                    else:
                        parhelp()
                elif event.key.keysym.scancode == sdl2.SDL_SCANCODE_SPACE:
                    pause = 1-pause
                    if not pause:
                        nfrstep = 0                                                    # number of frames behind viewing returned to zero
                elif event.key.keysym.scancode == sdl2.SDL_SCANCODE_RIGHT:
                    colorfunction = (colorfunction + 1) % 13
                    genelife.set_colorfunction(colorfunction)
                    print('step',framenr,'colorfunction changed to',colorfunction)
                elif event.key.keysym.scancode == sdl2.SDL_SCANCODE_LEFT:
                    colorfunction = (colorfunction - 1) % 13
                    genelife.set_colorfunction(colorfunction)
                    print('step',framenr,'colorfunction changed to',colorfunction)
                elif event.key.keysym.scancode == sdl2.SDL_SCANCODE_UP:
                    colorfunction2 = (colorfunction2 + 1)
                    if colorfunction2 == 13: colorfunction2 = -1
                    genelife.set_colorfunction2(colorfunction2)
                    print('step',framenr,'colorfunction2 changed to',colorfunction2)
                elif event.key.keysym.scancode == sdl2.SDL_SCANCODE_DOWN:
                    colorfunction2 = (colorfunction2 - 1)
                    if colorfunction2 == -2: colorfunction2 = 12
                    genelife.set_colorfunction2(colorfunction2)
                    print('step',framenr,'colorfunction2 changed to',colorfunction2)
                elif event.key.keysym.scancode == sdl2.SDL_SCANCODE_EQUALS or event.key == sdl2.SDL_SCANCODE_KP_PLUS:
                    if (colorfunction == 4) or (colorfunction == 5):
                        ymax = ymax * 2
                        oldymax = genelife.setget_act_ymax(ymax)
                        print('step',framenr,'new ymax =',ymax)
                    elif colorfunction == 10:
                        ymaxq = ymaxq * 2
                        oldymaxq = genelife.setget_act_ymaxq(ymaxq)
                        print('step',framenr,'new ymaxq =',ymaxq)
                    elif colorfunction == 11:
                        genealogycoldepth = genealogycoldepth + 1
                        genelife.set_genealogycoldepth(genealogycoldepth)
                        print('step',framenr,'new genealogycoldepth =',genealogycoldepth)
                elif event.key.keysym.scancode == sdl2.SDL_SCANCODE_MINUS or event.key == sdl2.SDL_SCANCODE_KP_MINUS:
                    if (colorfunction == 4) or (colorfunction == 5):
                        ymax = ymax // 2
                        oldymax = genelife.setget_act_ymax(ymax)
                        print('step',framenr,'new ymax =',ymax)
                    elif colorfunction == 10:
                        ymaxq = ymaxq // 2
                        oldymaxq = genelife.setget_act_ymaxq(ymaxq)
                        print('step',framenr,'new ymaxq =',ymaxq)
                    elif colorfunction == 11:
                        if genealogycoldepth > 0:
                            genealogycoldepth = genealogycoldepth - 1
                            genelife.set_genealogycoldepth(genealogycoldepth)
                            print('step',framenr,'new genealogycoldepth =',genealogycoldepth)
                elif event.key.keysym.scancode == sdl2.SDL_SCANCODE_1:
                    update1=1-update1
                    genelife.set_colorupdate1(update1)
                    if update1: print('step',framenr,'first window updates turned on')
                    else:       print('step',framenr,'first window updates turned off')
                elif event.key.keysym.scancode == sdl2.SDL_SCANCODE_2:
                    update2=1-update2
                    if update2: print('step',framenr,'second window updates turned on')
                    else:       print('step',framenr,'second window updates turned off')
                elif event.key.keysym.scancode == sdl2.SDL_SCANCODE_B:
                    nbhistmax=framenr//(N//2)
                    nbhistold = nbhist
                    if nbhistmax >= nNhist*2-1: nbhistmax=nNhist*2-2;
                    if sdl2.SDL_GetModState() & (sdl2.KMOD_LSHIFT|sdl2.KMOD_RSHIFT):
                        nbhist = nbhist+1
                        if nbhist > nbhistmax:
                            nbhist = -1
                    else:
                        nbhist = nbhist-1
                        if nbhist < -1:
                            nbhist=nbhistmax
                    if nbhist!=nbhistold:
                        print('step',framenr,"nbhist changed to ",nbhist)
                        genelife.set_nbhist(nbhist)
                elif event.key.keysym.scancode == sdl2.SDL_SCANCODE_C:
                    if colorfunction in [0,1,2,3,11,12]:
                        if sdl2.SDL_GetModState() & (sdl2.KMOD_LSHIFT|sdl2.KMOD_RSHIFT):
                            nfrstep = nfrstep + 1
                            if nfrstep == 8: nfrstep = 0
                        else:
                            nfrstep = nfrstep - 1
                            if nfrstep < 0: nfrstep = 8-1
                        print('step',framenr,"nfrstep changed to ",nfrstep)
                        # colorgrid(colorfunction,cgolg,cgrid,0,nfrstep)
                elif event.key.keysym.scancode == sdl2.SDL_SCANCODE_D:
                    print('step',framenr,"simulation continues with stashed patterns")
                    genelife.unstash()
                elif event.key.keysym.scancode == sdl2.SDL_SCANCODE_E:
                    if colorfunction in [9,10]:
                        if sdl2.SDL_GetModState() & (sdl2.KMOD_LSHIFT|sdl2.KMOD_RSHIFT):
                            genelife.label2stash(1)       # stash pattern cumulatively
                            print('step',framenr,"selected pattern added to stash")
                        else:
                            genelife.label2stash(0)       # clear stash and add pattern
                            print('step',framenr,"selected pattern isolated to stash")
                elif event.key.keysym.scancode == sdl2.SDL_SCANCODE_F:
                    print("entering key F")
                    if   sdl2.SDL_GetModState() &  (sdl2.KMOD_LSHIFT|sdl2.KMOD_RSHIFT):
                        print("entering shift key F")
                        if scalex2 or event.window.windowID == windowID2:
                            windowsize=(2*Width, 2*(Height+16))
                        else:
                            windowsize=(Width, Height+16)
                        if window.get_flags() & pg.FULLSCREEN:
                            window = pg.display.set_mode(windowsize)
                            # screen = pg.display.set_mode(screensize,pg.DOUBLEBUF|pg.OPENGL,32)
                        else:
                            sdl2.SDL_SetWindowFullscreen(window)
                    else:
                        print("entering key f")
                        print("no of frames per second (av. last 10) = %f" % framerate)
                elif event.key.keysym.scancode == sdl2.SDL_SCANCODE_G:
                    if colorfunction == 9:
                        gcolor = (gcolor+1)%10;
                        genelife.set_gcolors()
                        print('step',framenr,'new gcolor =',gcolor)
                    elif colorfunction == 6 or colorfunction == 7 or colorfunction == 11:
                        ancestortype= ancestortype + 1
                        if ancestortype > 2:
                            ancestortype = 0
                        genelife.set_ancestortype(ancestortype)
                        if ancestortype == 2:
                            print('step',framenr,'new ancestor choice clonal for win 2 first for win 1')
                        elif ancestortype == 1:
                            print('step',framenr,'new ancestor choice clonal')
                        else:
                            print('step',framenr,'new ancestor choice first')
                elif event.key.keysym.scancode == sdl2.SDL_SCANCODE_I:
                    nbhoods=[7,3,5,7]
                    info_transfer_h = info_transfer_h + 1;
                    if (info_transfer_h == 4):  info_transfer_h = 0;
                    if (info_transfer_h != 0): do_info_transfer = 1;
                    else: do_info_transfer = 0;
                    genelife.set_info_transfer_h(do_info_transfer,nbhoods[info_transfer_h])
                    print('step',framenr,'info_transfer_h =',info_transfer_h)
                elif event.key.keysym.scancode == sdl2.SDL_SCANCODE_L:
                    activityfnlut = 1 - activityfnlut
                    genelife.set_activityfnlut(activityfnlut)
                    print('step',framenr,'activityfnlut =',activityfnlut)
                elif event.key.keysym.scancode == sdl2.SDL_SCANCODE_N:
                    if   sdl2.SDL_GetModState() &  (sdl2.KMOD_LSHIFT|sdl2.KMOD_RSHIFT):
                        genelife.get_nnovelcells(nnovelcells)
                        print('step',framenr,"collected trace of nnovelcells")
                    else:
                        noveltyfilter=1-noveltyfilter
                        print('step',framenr,"noveltyfilter changed to ",noveltyfilter)
                        genelife.set_noveltyfilter()
                elif event.key.keysym.scancode == sdl2.SDL_SCANCODE_P:
                    activity_size_colormode=(activity_size_colormode+1)%4
                    print('step',framenr,"activity_size_colormode changed to ",activity_size_colormode)
                    genelife.set_activity_size_colormode()
                    pixeldat=actsizecoltxt[activity_size_colormode]
                elif event.key.keysym.scancode == sdl2.SDL_SCANCODE_Q:
                    if   sdl2.SDL_GetModState() & sdl2.KMOD_LALT:
                        quadrants = eval(input("Enter an integer between -1 and 6: "))
                    elif   sdl2.SDL_GetModState() &  sdl2.KMOD_SHIFT:
                        if quadrants < 7: quadrants = quadrants+1
                    else:
                        if quadrants >= 0: quadrants = quadrants-1
                    print('step',framenr,"quadrants changed to ",quadrants)
                    genelife.set_quadrant(quadrants)
                elif event.key.keysym.scancode == sdl2.SDL_SCANCODE_R:
                    if   sdl2.SDL_GetModState() & sdl2.KMOD_LALT:
                        rbackground,randominflux = eval(input("Enter rbackground [0-32768] and randominflux (3 deletions only, 2 GoL gene, 1 random gene:"))
                        print('step',framenr,"rbackground changed to ",rbackground,"with random mode ",randominflux)
                        genelife.set_rbackground(rbackground,randominflux)
                    elif   sdl2.SDL_GetModState() &  sdl2.KMOD_LSHIFT:
                        randominflux = 2 if randominflux !=2 else 0
                        print('step',framenr,"randominflux changed to ",randominflux)
                        genelife.set_randominflux(randominflux)
                    elif   sdl2.SDL_GetModState() &  sdl2.KMOD_RSHIFT:
                        randominflux = 3 if randominflux !=3 else 0
                        print('step',framenr,"randominflux changed to ",randominflux)
                        genelife.set_randominflux(randominflux)
                    else:
                        randominflux = 1 if randominflux !=1 else 0
                        print('step',framenr,"randominflux changed to ",randominflux)
                        genelife.set_randominflux(randominflux)
                elif event.key.keysym.scancode == sdl2.SDL_SCANCODE_S:
                    if   (sdl2.SDL_GetModState() &  (sdl2.KMOD_LSHIFT | sdl2.KMOD_RSHIFT)) and (colorfunction2 != -1):
                        fname = "./images/genelife_sel%02d_t%03d_r%08x_s%03d_c%02d.jpeg" % (selection,framenr,repscheme,savecnt,colorfunction2)
                        err = sdl2.sdlimage.IMG_SavePNG(surface2,bytes(fname, encoding="ascii"))
                    else:
                        fname = "./images/genelife_sel%02d_t%03d_r%08x_s%03d_c%02d.jpeg" % (selection,framenr,repscheme,savecnt,colorfunction)
                        # err = sdl2.SDL_SaveBMP(ctypes.byref(surface),bytes(fname, encoding="ascii")) # does not work because of surface type error
                        err = sdl2.sdlimage.IMG_SavePNG(surface,bytes(fname, encoding="ascii"))
                    if err:
                        print("error %d file not saved" % err)
                    else:
                        print("image saved "+fname)
                    savecnt = savecnt + 1
                elif event.key.keysym.scancode == sdl2.SDL_SCANCODE_T:
                    if   sdl2.SDL_GetModState() &  (sdl2.KMOD_LSHIFT|sdl2.KMOD_RSHIFT):
                        if(offdt<0): offdt = offdt+1
                    elif offdt>-maxPlane+1: offdt = offdt-1
                    print('step',framenr,"offset dt changed to ",offdt)
                    genelife.set_offsets(offdx,offdy,offdt)
                elif event.key.keysym.scancode == sdl2.SDL_SCANCODE_V:
                    vscrolling=1-vscrolling
                    print('step',framenr,"vscrolling changed to ",vscrolling)
                    genelife.set_vscrolling()
                elif event.key.keysym.scancode == sdl2.SDL_SCANCODE_X:
                    if   sdl2.SDL_GetModState() &  (sdl2.KMOD_LSHIFT|sdl2.KMOD_RSHIFT): offdx = offdx+1
                    else: offdx = offdx-1
                    print('step',framenr,"offset dx changed to ",offdx)
                    genelife.set_offsets(offdx,offdx,offdt)
                elif event.key.keysym.scancode == sdl2.SDL_SCANCODE_Y:
                    if   sdl2.SDL_GetModState() &  (sdl2.KMOD_LSHIFT|sdl2.KMOD_RSHIFT): offdy = offdy+1
                    else: offdy = offdy-1
                    print('step',framenr,"offset dy changed to ",offdy)
                    genelife.set_offsets(offdx,offdy,offdt)
        if (not mouseclicked):
            if updatesenabled and not pause and (framenr < maxsteps):
                update_sim(nrun, ndisp, nskip, niter, nhist, nstat, count)
            else:
                if update1:
                    colorgrid(colorfunction,cgolg,cgrid,0,nfrstep)
                if update2:
                    if colorfunction2 != -1:
                        colorgrid(colorfunction2,cgolg2,cgrid2,1,0)

        nspecies=genelife.get_nspecies()
        caption=construct_caption(colorfunction,pixeldat,buttonhelp,1)
        if colorfunction2 == -1:
            caption2 = "2nd "+caption
        else:
            caption2=construct_caption(colorfunction2,pixeldat,buttonhelp,2)
        
        if update1:                                         # window 1 updates switched on
            set_caption(window, caption)                    # sets the window caption with current status data
            if scalex2:
                sdl2.SDL_BlitScaled(surfacex1,None,surfacex2,None)
            sdl2.ext.Window.refresh(window)                 # copies the window to the display
        elif framenr % 10 == 0: set_caption(window, caption)# only less frequent caption update and no graphics update of window 1

        if update2:                                         # if update of 2nd window on
            if render2:
                if colorfunction2 == -1:
                    sdlrenderer2 = renderer2.sdlrenderer
                    # sdl2.SDL_SetRenderDrawColor( sdlrenderer2, 0xFF, 0xFF, 0xFF, 0xFF )
                    sdl2.SDL_SetRenderDrawColor( sdlrenderer2, 0x00, 0x00, 0x00, 0x00 )
                    sdl2.SDL_RenderClear( sdlrenderer2 )
                    sdl2.ext.fill(surface2.contents, 0)
                    sdl2.SDL_BlitScaled(surface,None,surface2,None)
                    image2=factory2.from_surface(surface2)
                    sdl2.SDL_RenderCopy(sdlrenderer2, image2.texture, None, None)  # last two parameters are source and dest rect (e.g. grect)
                    message = sdl2.sdlttf.TTF_RenderText_Solid( font2, str.encode("Genelife frame %5d" % cnt), textColor2 )
                    message2 = factory2.from_surface(message)
                    grect.w = message.contents.w
                    grect.h = message.contents.h
                    set_caption(window2, caption2)
                    sdl2.SDL_RenderCopy(sdlrenderer2, message2.texture, None, grect)
                    sdl2.SDL_RenderPresent(sdlrenderer2)
                    
                    """
                    rect = sdl2.SDL_Rect()
                    rect.x = 0
                    rect.y = 0
                    rect.w = Width
                    rect.h = Height+16
                    # https://stackoverflow.com/questions/21651976/how-to-pass-sdl-surface-to-sdl-locktexture-with-pysdl2
                    sdl2.SDL_LockTexture(texture2, rect, ctypes.byref(ctypes.c_void_p(surface2.contents.pixels)), ctypes.byref(ctypes.c_int(surface2.contents.pitch)))
                    # https://stackoverflow.com/questions/4355524/getting-data-from-ctypes-array-into-numpy
                    # buffer = np.core.multiarray.int_asbuffer(surface2.contents.pixels, 8*Width*Height)  # at runtime says xxx not implemented
                    buffer_from_memory = ctypes.pythonapi.PyMemoryView_FromObject
                    buffer_from_memory.restype = ctypes.py_object
                    #buffer = buffer_from_memory(ctypes.c_void_p(surface2.contents.pixels), 8*Width*Height)
                    buffer = buffer_from_memory(surface2.contents.pixels, 8*Width*Height)
                    # nppixels = typeslib.as_array((ctypes.c_uint32 * (Height*Width)).from_address(ctypes.c_void_p(surface2.contents.pixels)))
                    nppixels = np.frombuffer(buffer, ctypes.c_uint32)
                    nppixels[:] = cgolg[:]
                    sdl2.SDL_UnlockTexture(texture2)
                    sdl2.SDL_RenderCopy( renderer2, texture2, rect, None);
                    sdl2.SDL_RenderPresent( renderer2 )
                    """
                else:
                    winflags = sdl2.SDL_GetWindowFlags(window2.window)
                    if winflags & sdl2.SDL_WINDOW_HIDDEN:
                        sdl2.SDL_ShowWindow(window2.window)
                        sdl2.SDL_RaiseWindow(window2.window)
                    set_caption(window2, caption2)
                    if scalex2:
                        if sdl2.SDL_BlitScaled(surface2x1,None,surface2x2,None):
                            print("error",sdl2.SDL_GetError(),"executing BlitScaled at step",framenr)
                    sdl2.ext.Window.refresh(window2)                # copies the 2nd window to the display
            else:
                set_caption(window2, caption2)
                if colorfunction2 == -1:
                    sdl2.SDL_HideWindow(window2.window)
                    # if sdl2.SDL_BlitScaled(surfacex1,None,surface2,None):
                    #    print("error",sdl2.SDL_GetError(),"executing BlitScaled at step",framenr)
                else:
                    winflags = sdl2.SDL_GetWindowFlags(window2.window)
                    if winflags & sdl2.SDL_WINDOW_HIDDEN:
                        sdl2.SDL_ShowWindow(window2.window)
                        sdl2.SDL_RaiseWindow(window2.window)
                    if sdl2.SDL_BlitScaled(surfacex1,grect1,surface2x1,grect1): print("BlitScaled failed")
                    if scalex2:
                        if sdl2.SDL_BlitScaled(surface2x1,None,surface2,None):
                            print("error",sdl2.SDL_GetError(),"executing BlitScaled at step",framenr)
                sdl2.ext.Window.refresh(window2)                # copies the 2nd window to the display
                # print("executing BlitScaled at step",framenr)
        
        if framenr % 10 == 0:
            mslasttime = mstime
            mstime = sdl2.timer.SDL_GetTicks()
            framerate = 1./((mstime-mslasttime)/10000.0)
        # sdl2.ext.time.Clock.tick()                   # requires modified pysdl2 see sdl2.ext.time.Clock above
        # framerate = sdl2.ext.time.Clock.get_fps()    # requires modified pysdl2 see sdl2.ext.time.Clock above

#-----------------------------------------------------------------------------------------------------------

def parhelp():
    """ definition of parameters"""
    if selection < 8:
        print("Control bits (left to right, also in caption on right mouse click): ")
        print("____________________________")
        print("Green    ","0. selective birth for 3-live-nbs  ","1. selective birth for 2-live-nbs")
        print("Mid Blue ","2. canonical 0 position vs difft   ","3. bypass selection for 2-live-nbs")
        print("Teal blue","4. enforce birth for 3-live-nbs    ","5. enforce birth for 2-live-nbs")
        print("Green    ","6. 2nd neighbour genetic modulation","7. 1st neighbour genetic masking")
        print("Yellow   ","8. enforce GoL rule if non GoL rule","9. enforce GoL rule last change by non GoL rule */")
        print("Green    ","10-13. allow 2-nb birth only for active subset of 4 canonical configs")
        print("Blue     ","14. Survival for 3-live-nbs        ","15. Survival for 2-live-nbs")
        print("Red      ","16. Gene overwrite for 3-live-nbs  ","17. Gene overwrite for 2-live-nbs")
        print("Purple   ","18. Parentdies forced")
    elif selection < 16:
        print("")
        if selection < 10:
            print("Control bits for masks to enable gene encoded LUTs (left to right):")
            print("___________________________________________________________________")
            print("Masks for LUT rules for sums 2-6 for survival /left half) and birth (right half)")
            print("blue         ","survival  for sum=1-8 separate buttons for the 8 non-zero sums")
            print("green        ","birth     for sum=1-8 separate buttons for the 8 non-zero sums")
        elif selection < 12:
            print("Control bits for masks to enable gene encoded LUTs (left to right):")
            print("___________________________________________________________________")
            print("Distance dept (2 classes) LUT rules for s,se s 1-7 for survival /left half) and birth (right half)")
            print("blue-purple ","survival sum=1-7 separate buttons for s = 1-7 se = 0-1 0-2 0-3 0-4 1-4 2-4 3-4")
            print("green-yellow","birth    sum=1-7 separate buttons for s = 1-7 se = 0-1 0-2 0-3 0-4 1-4 2-4 3-4")
        elif selection < 14:
            print("Control bits for masks to enable gene encoded LUTs (left to right):")
            print("___________________________________________________________________")
            print("Canonical rotation LUT rules for sums 2-6 for survival /left half) and birth (right half)")
            print("blue-purple ","survival sum=2-6 separate buttons for s = 2-6 canonical rotns = 4 7 10 7 4")
            print("green-yellow","birth    sum=2-6 separate buttons for s = 2-6 canonical rotns = 4 7 10 7 4")
        elif selection < 16:
            print("Control bits for masks to enable gene encoded LUTs (left to right):")
            print("___________________________________________________________________")
            print("Full 2D symmetry LUT rules for sums 0-4 for survival /left half) and birth (right half)")
            print("blue-purple ","survival sum=0-4 separate buttons for s = 0-4 pattern offset = 1 2 6 10 13")
            print("green-yellow","birth    sum=0-4 separate buttons for s = 0-4 pattern offset = 1 2 6 10 13")
        print("red          ","overwrite for sum=1-8 separate buttons for the 8 non-zero sums")
        print("none         ","ancselectmask for s=1-8 whether to do gene selection on ancestor")
        print("grey         ","repscheme 0: central/nb survival, 1: OR/AND gene LUTs, 2: anc 0-canonical/most difft, 3: parentdies forced/not")
        print("teal         ","repscheme 4-7: one of 16 ancestor selection schemes in selectone_of_s, 7: choose golr selection modes ")
        print("brown        ","repscheme 8-10: one of 8 disambiguation schemes, 11: random choice of ancestor from all live")
        print("")
        print("Control bits for selection are in repscheme")
        print("______________________________________________________")
        print("repscheme bits 0-3 : ")
        print("   bit 0 determines choice of gene(s) for survival rule : OR of live neighbors (0) central gene (1)")
        print("   bit 1 determines combination of nb gene function for even LUT selection 8-14: AND (0) OR (1)")
        print("   bit 2 determines choice of neighbour in canonical rotation : most central/different (0) or first bit (1)")
        print("   bit 3 determines whether parentdies set: parents forced to die (1) or not (0)")
        print("repscheme bits 4-7 determine selection scheme based on gene")
        print("   # 0 minimum gene as value (penalizes proliferation) # 1 maximum gene as value (rewards proliferation)")
        print("   # 2 minimum number of ones (penalizes proliferation) # 3 maximum number of ones (rewards proliferation)")
        print("   # 4 neutral selection # 5 neutral but different selection (will not survive homogeneous gene start)")
        print("   # 6 penalty function -1 for a survival rule -2 for a birth rule  # 7 reserved, currently same as 6")
        print("   # 8 selection for minimum period glider # 9 selection for maximum period glider ")
        print("   # 10 NYI same as 8 #11 NYI same as 9")
        print("   # 12 NYI same as 8 #13 NYI same as 9")
        print("   # 14 NYI same as 8 #15 NYI same as 9")
        print("repscheme bits 8-10 determine disambiguation method for symmetric cases sum=2,crot=3 and sum=4,crot=2,9")
        print("   0 random choice amongst as yet unresolved: this involves a departure from determinism for these cases")
        print("   1 ignore problem and choose selected bit of canonical configuration : live with minimal asymmetry")
        print("   2 disallow birth : effectively modifies the rules and is like excluding this LUT entry from the table")
        print("   3 choose a default gene such as the gene coding for the Game of Life LUT to resolve ambiguity")
        print("   4 choose lesser in value of genes if different (otherwise it makes no difference)")
        print("   5 choose first gene which is most likely to not survive with s<2 or s>3")
        print("   6 choose first gene which is most likely to not survive or to be overwritten with (s<2 or s>3) xor overwrite")
        print("   7 generate a random gene to give birth to for these ambiguous instances")
        print("repscheme bit 11 overrides bits 4-7 with random choice of ancestor amongst live neighbours")
    print("")
    print("Other controls:")
    print("_______________")
    print("<space>     ","pause simulation, allowing ongoing display control")
    print("middle mouse","stop simulation [data is retained for possible run() for run/analysis with updatesenabled=True/False]")
    print("left mouse  ","extract information about local state inside the array, or control buttons below")
    print("right mouse ","choose single plane for GoL display in colorfunction 2 for selection 16-19")
    print("<- , ->     ","decrement or increment the colorfunction analysis type mod 12")
    print("Down , Up   ","decrement or increment the colorfunction2 analysis type mod 12 : -1 for same as colorfunction")
    print("1, 2        ","toggle first or second window updates off/on")
    print("b , B       ","decrement or increment the half block for trace display: in range -1,0 to nNhist*2-2=38")
    print("c , C       ","decrement or increment the playback counter to replay last up to 8 timesteps while paused (with space), colorfn 0-3,11,12 only")
    print("d           ","continue run with stashed pattern(s) only")
    print("e , E       ","extract current selected pattern to stash e: singly E:cumulatively")
    print("f           ","print frame rate in fps (average of last 10 frames NYI")
    print("F           ","toggle to fullscreen NYI")
    print("g           ","toggle on/off inherited coloring of connected cpts from overlapping cpts (colfn 9) or cycle first/clonal/both ancestors (6,7,11)")
    print("h           ","print this help")
    print("H           ","toggle horizon mode on or off: upper half of array obeys unmodified GoL rule")
    print("i           ","toggle display and calculation of info_transfer_histogram on/off")
    print("l           ","toggle activity of aggregate functional genes activityfnlut on/off")
    print("n , N       ","n: toggle novelty filter on/off for connected component color function 9 N: get trace of nnovecells")
    print("p           ","rotate activity_size_colormode 0,1,2,3 for (no,log2n,pixels,sqrt(pixels)) size display of activities in color function 10")
    print("q , Q       ","incr or decr quadrant parameter choice : -1 = no quadrants, 0-4 are first 5 bit pairs of repscheme, 5,6 surv and overwrite")
    print("r           ","toggle random influx domain on or off")
    print("R leftshift ","toggle intermittent feathered random influx domain on or off")
    print("R rightshift","toggle random deletion perturbation at rbackground rate on or off")
    print("alt-r       ","input background rate of perturbation for random influx")
    print("s , S       ","save current image from window 1 (s) 2 (S) to file in image subdirectory")
    print("x,X y,Y t,T ","lower (lc) or raise (uc) the (dx,dy,dt) offsets for glider tracking (colorfn 8) (0,0,0)=(all 8 nnb dt=-1)")
    print("v           ","toggle vertical scroll tracking mode : following top most objects and losing lowest objects in contact with 0 row")
    print("+ , -       ","increase or decrease ymax or ymaxq for activity display scaled as act/(ymax+act) by a factor of 2")
    print("+ , -       ","increase or decrease genealogycoldepth if colorfunction==11 for gene ancestor display")
#-----------------------------------------------------------------------------------------------------------

if __name__ == '__main__':
    """ main program to run with default example parameters (see global parameters above) """
    
    genelife.initialize_planes(npoffsets)
    genelife.initialize(runparams,simparams)
    framenr = 0
    cnt=0
    show0()
    # step()
    run(nrun, ndisp, nskip, niter, nhist, nstat, cnt, 1000)
    
