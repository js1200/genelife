//
//  external_set.c
//  project genelife
//
//  external setting of various quantities (typically called from python via genelife_c_interface.py)
//---------------------------------------------------------- copyright ----------------------------------------------------------------------------------
//  Written by John S. McCaskill and Norman H. Packard 2017-2019
//  First created by John McCaskill on 14.07.2017. Last modified Nov 2019.
//
/*  MIT License
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
*/
//-------------------------------------------------------------------------------------------------------------------------------------------------------
//
#include "genelife.h"
//
// set_colorfunction    set color function integer from GUI for use in patterning and coloring display
// setget_act_ymax      set activity ymax for scaling of gene activity plot
// setget_act_ymaxq     set activity ymax for scaling of quad activity plot
// set_selectedgene     set selected gene for highlighting from current mouse selection in graphics window
// set_offsets          set offsets for detection of glider structures in display for color function 8
// set_quadrant         set the pair of bits in repscheme (or survivalmask or overwritemask) used for quadrant variation 0-6
// set_randominflux     change the randominflux activation for rbackground if nonzero or continual updating of init field with random states and genes : 2,1,0
// set_rbackground      set the backround random live gene input rate per frame and site to rbackground/32768
// set_repscheme_bits   set the two of the repscheme (or survivalmask or overwritemask) bits corresponding to the selected quadrant
// set_repscheme        set repscheme from python
// set_rulemod          set rulemod from python
// set_surviveover      set the two masks for survival and overwrite from python (survivalmask, overwritemask)
// set_vscrolling       set vertical scrolling to track fronts of growth in vertical upwards direction
// set_noveltyfilter    set novelty filter for darkening already encountered components in connected component display (colorfunction 9)
// set_activity_size_colormode set colormode by size for colorfunction 10 : 0 by ID  1 log2 enclosing square size 2 use #pixels 3 use sqrt(#pixels)
// set_gcolors          set connected component colors as inherited colors from colliding connected components with random drift
// set_seed             set random number seed
// set_nbhist           set nbhist N-block of time points for trace from GUI for use in activity and population display traces
// set_genealogycoldepth set genealogycoldepth for colorfunction=11 display
// set_ancestortype     set ancestortype for genealogy display and return of first (0), clonal (1) or first & clonal in 2 windows (2)
// set_info_transfer_h  set information transfer histogram display value (0,1) from python
// set_activityfnlut    set collection of functional activity statistics corresponding to functional aggregate of genes by non-neutral bits
// set_colorupdate1     control update of colorgenes and regular print statements via flag colorupdate1
// set_colorfunction2   choice of colorfunction for window 2
//-------------------------------------------------------------------- set ...--------------------------------------------------------------------------
void set_colorfunction(int colorfunctionin) {
    if((colorfunctionin>12) || (colorfunctionin<0)) fprintf(stderr,"error colorfunction value passed %d out of range\n",colorfunctionin);
    else     colorfunction = colorfunctionin;
}
//.......................................................................................................................................................
void set_colorfunction2(int colorfunctionin) {
    if((colorfunctionin>12) || (colorfunctionin<-1)) fprintf(stderr,"error colorfunction value passed %d out of range\n",colorfunctionin);
    else     colorfunction2 = colorfunctionin;
}
//.......................................................................................................................................................
int setget_act_ymax(int actymax) {                  // sets ymax for activities only if argument nonzero, reads old value
    int ymaxold;
    ymaxold = ymax;
    ymax = actymax;
    return(ymaxold);
}
//.......................................................................................................................................................
int setget_act_ymaxq(int actymaxq) {                  // sets ymax for activities only if argument nonzero, reads old value
    int ymaxqold;
    ymaxqold = ymaxq;
    ymaxq = actymaxq;
    return(ymaxqold);
}
//.......................................................................................................................................................
void set_selectedgene(uint64_t gene) {
    selectedgene=gene;
    fprintf(stderr,"selected gene set to %"PRIx64"\n",selectedgene);
}
//.......................................................................................................................................................
void set_offsets(int dx,int dy,int dt) {
    offdx =dx;
    offdy = dy;
    if(dt>0) {
        dt=0;
        fprintf(stderr,"positive time offsets not allowed, looking into the future not possible\n");
    }
    if(dt<=-maxPlane) {
        dt=-maxPlane+1;
        fprintf(stderr,"not enough planes set for this time offset, recompile software with larger maxPlane value\n");
    }
    offdt = dt;
}
//.......................................................................................................................................................
void set_quadrant(int quadrant) {
    if (quadrant >= -1 && quadrant < 7) quadrants = quadrant;
    repscheme &= ~R_quadrant;                                           // remove all quadrant bits : only one set at a time in interactive version
    if(quadrant >= 0 && quadrant < 7) {
        repscheme |= R_16_quadrant_sele<<quadrant;                          // assumes quadrant selectors are 7 successive bits following R_16_...
    }
}
//.......................................................................................................................................................
void set_randominflux(int randominfluxin) {
    randominflux=randominfluxin;
}
//.......................................................................................................................................................
void set_rbackground(int rbackgroundin, int randominfluxin) {
    rbackground=rbackgroundin;
    randominflux=randominfluxin;
    if(!rbackground) randominflux=0;   // do not leave patch randominflux variable active when turning off random background
}
//.......................................................................................................................................................
unsigned int set_repscheme_bits(int quadrant, int x, int y, unsigned int surviveover[]) {
    unsigned int quadrantval;

    quadrantval=(x<(Nmask>>1)? 0 : 1) + (y<(Nmask>>1)? 0 : 2);              // determine selected quadrant
    if(quadrants >= 0 && quadrants < 5) {                                   // assumes repscheme bits in pairs starting from bit 0 matching quadrants
        repscheme &=  ~(0x3llu<<(quadrants<<1));
        repscheme |=  quadrantval<<(quadrant<<1);
    }
    else if (quadrant < 6) {
        survivalmask &= ~0x3u;
        survivalmask|= quadrantval;
    }
    else if (quadrant<7) {
        overwritemask &= ~0x3u;
        overwritemask|= quadrantval;
    }
    repscheme &= ~R_quadrant;                                                // remove all quadrant bits : only one set at a time in interactive version
    quadrants = -1;                                                          // reset internal quadrants choice so that full display is shown
    surviveover[0]=survivalmask;
    surviveover[1]=overwritemask;
    return(repscheme);
}
//.......................................................................................................................................................
void set_repscheme(unsigned int repscheme_in) {
    repscheme = repscheme_in;
}
//.......................................................................................................................................................
void set_rulemod(unsigned int rulemod_in) {
    rulemod = rulemod_in;
}
//.......................................................................................................................................................
void set_surviveover64(unsigned int surviveover[], int len ) {
    if (len==3) {
        survivalmask = surviveover[0];
        if(selection<8) overwritemask = surviveover[1];
        else {
            birthmask = surviveover[1];
            overwritemask = surviveover[2];
        }
    }
    else fprintf(stderr,"surviveover64 needs three parameters, %d provided\n",len);
}
//.......................................................................................................................................................
void set_vscrolling() {
    vscrolling=1-vscrolling;
    if(!vscrolling) last_scrolled = 0;
}
//.......................................................................................................................................................
void set_noveltyfilter() {
    noveltyfilter=1-noveltyfilter;
}
//.......................................................................................................................................................
void set_activity_size_colormode() {
    activity_size_colormode = (activity_size_colormode+1) % 4;
}
//.......................................................................................................................................................
void set_gcolors() {
    gcolors = (gcolors+1)%10;
}
//.......................................................................................................................................................
void set_seed(int seed) {
    ranseed = seed;
}
//.......................................................................................................................................................
void set_nbhist(int nbhistin) {
    if(nbhist<nNhist*2) nbhist=nbhistin;
    else fprintf(stderr,"nbhist out of range %d > %d\n",nbhistin,nNhist*2-1);
}
//.......................................................................................................................................................
void set_genealogycoldepth(int genealogycoldepthin) {
    genealogycoldepth = genealogycoldepthin;
}
//.......................................................................................................................................................
void set_ancestortype(int ancestortypein) {
    if(ancestortypein <3) ancestortype = ancestortypein;
    else fprintf(stderr,"ancestor type %d out of range [0..2]\n",ancestortypein);
}
//.......................................................................................................................................................
void set_info_transfer_h(int do_info_transfer, int nbhood) {
    info_transfer_h = do_info_transfer;
    if(nbhood == 3 || nbhood == 5 || nbhood == 7)
        it_nbhood = nbhood;
    else fprintf(stderr,"error in nbhood value %d, allowed values are 3,5,7\n",nbhood);
}
//.......................................................................................................................................................
void set_activityfnlut(int activityfnlutin) {
    activityfnlut = activityfnlutin;
}
//.......................................................................................................................................................
void set_colorupdate1(int update1) {
    colorupdate1 = update1;
}
