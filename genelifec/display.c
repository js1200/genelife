//
//  display.c
//  genelife
//
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
// set_color            function to assign rainbow colors
// label_color          map label (quasi-uniformly) into larger unsigned colour space
// rgba                 converts hue (0..1) to rgb+alpha
// mix_color            mix colors from overlapping components, add random drift of colour
// delay                time delay in ms for graphics
// printxy              simple terminal print of array
// printscreen          terminal screen print of array on xterm, moving cursor with esc codes */
// colorgenes           colour display of genes in one of 12 colorfunction modes, including activities, pattern analysis, genealogies and glider detection
//--------------------------------------------------------------- colorgenes ----------------------------------------------------------------------------
extern INLINE void setcolor(unsigned int *color,int n) { // for coloring by quad size...
    // rainbow colors from running from fields::tim.colors(11) in R by Tim Hoar : assuming log2N <= 10 (limited by display size)
    // "#00008F" "#0000F5" "#005AFF" "#00BDFF" "#23FFDC" "#87FF78" "#ECFF13" "#FFAD00" "#FF4A00" "#E40000" "#800000"
    static const unsigned int rainbow[]={0x800000FF, 0xE40000FF, 0xFF4A00FF, 0xFFAD00FF, 0xECFF13FF, 0x87FF78FF, 0x23FFDCFF, 0x00BDFFFF, 0x005AFFFF, 0x0000F5FF, 0x00008FFF};

    unsigned int mycol;
    mycol = rainbow[n];
    color[2] = (mycol >> 24) & 0xff; // B
    color[1] = (mycol >> 16) & 0xff; // G
    color[0] = (mycol >> 8) & 0xff;  // R
}
//.......................................................................................................................................................
extern INLINE unsigned int labelcolor( unsigned int label) {
    uint64_t mask;
    unsigned int color;
    mask = (uint64_t) label;
    mask = mask * 11400714819323198549ull;          // map label (quasi-uniformly) into larger unsigned colour space
    color = mask >> (64 - 32);                      // hash with optimal prime multiplicator down to 32 bits
    color |= 0x080808ffull;                         // ensure visible (slightly more pastel) color at risk of improbable redundancy, make alpha opaque
    return color;
}
//.......................................................................................................................................................
extern INLINE unsigned int rgba( float hue) {                            // converts hue (0..1) to rgb+alpha
    unsigned int r,g,b;
    unsigned int k;
    float hue6 = hue*6.0;

    k = (unsigned int) (6.*hue); if (k==6) k = 0;  // deal with possible round off error
    r=g=b=0x00;
    switch (k) {
        case 0: r = 0xff; g = (unsigned int)((hue6)*256.0);   g = g>0xff?0xff:g;break;
        case 1: g = 0xff; r = (unsigned int)((2.-hue6)*256.0);r = r>0xff?0xff:r;break;
        case 2: g = 0xff; b = (unsigned int)((hue6-2.)*256.0);b = b>0xff?0xff:b;break;
        case 3: b = 0xff; g = (unsigned int)((4.-hue6)*256.0);g = g>0xff?0xff:g;break;
        case 4: b = 0xff; r = (unsigned int)((hue6-4.)*256.0);r = r>0xff?0xff:r;break;
        case 5: r = 0xff; b = (unsigned int)((6.-hue6)*256.0);b = b>0xff?0xff:b;break;
        default: fprintf(stderr,"step %d Error in hut to rgb color switch - case %d out of bounds\n",totsteps,k);
    }
    return ((r<<8) | (g<<16) | (b<<24) | 0xff);
}
//.......................................................................................................................................................
extern INLINE float mixcolor( unsigned int label, uint64_t rand) { // mix colors from overlapping components, add random drift of colour
    unsigned int conn;
    float color,color1,x1,y1,x,y;
    float eps = 0.0001*gcolors;
    #define PI 3.14159265
    const float i2pi = 1./(2.*PI);

    x = y = 0.0;
    conn = connlists[label];                                            // NB conn is not a label but an index in the array of possible connections
    while(conn) {
        color1= oldcomplist[connections[conn].oldlab].gcolor;           // mix colors based on aoverlap weights
        x1 = cosf(2.*PI*color1);                                        // averaging of circular variable requires converting to x,y vector and average
        y1 = sinf(2.*PI*color1);
        /* color2 = atan2f(y1,x1)*i2pi;
        color2 = color2 < 0. ? color2 + 1. : color2;
        if (fabsf(color2-color1) > eps) fprintf(stderr,"step %d label %d non inverted arctan color1 %f color2 %f\n",totsteps,label,color1,color2);*/
        x += x1*connections[conn].aoverlap;
        y += y1*connections[conn].aoverlap;
        conn=connections[conn].next;
    }
    color = atan2f(y,x)*i2pi;                                           // atan2 returns principal value in range [-PI, PI]
    color = color < 0. ? color + 1. : color;                            // for negative values add 2*PI to get back to range [0,2PI]
    color += eps*(((float)(rand&0xff)) - 127.5);                        // random drift of color
    color = color>1.0 ? color-1.0 : color;
    color = color<0.0 ? color+1.0 : color;
    return color;
}
//.......................................................................................................................................................
void delay(int milliseconds) {
    long pause;
    clock_t now,then;
    pause = milliseconds*(CLOCKS_PER_SEC/1000);
    now = then = clock();
    while( (now-then) < pause )
        now = clock();
}
//.......................................................................................................................................................
void printxy (uint64_t gol[],uint64_t golg[]) {                         // print the game of life configuration
    int    ij, col, X, Y;
    // https://stackoverflow.com/questions/27159322/rgb-values-of-the-colors-in-the-ansi-extended-colors-index-17-255
    for (ij=0; ij<N2; ij++) {
        if(gol[ij]>0){
            col = 32+((golg[ij]>>57)&0x7f);
            X = ij % N;
            Y = ij / N;
            printf("%d %d %d ",col,X,Y);
      }
    }
    printf("\n");
}
//.......................................................................................................................................................
void printscreen (uint64_t gol[], uint64_t golg[]) {   /* print basic genelife on xterm, moving cursor with esc codes */
//  The simulation uses the terminal text output as a colour display for the GoL.
//  In order to fit N=128 GoL display on cinema display, use terminal preferences to change font spacings column 1.3 and line 0.65 at 12 pt
//  In order to fit N=128 GoL display on smaller displays, use terminal preferences to change font spacings column 1.0 and line 0.5 at 10 pt
//  To allow escape codes to move cursor, in terminal profiles select "Allow Vt100 application keypad mode"
    int    ij, col;
    // https://stackoverflow.com/questions/27159322/rgb-values-of-the-colors-in-the-ansi-extended-colors-index-17-255
    printf("\e[38;5;255;48;5;238m");
    for (ij=0; ij<N2; ij++) {
        col = 32+((golg[ij]>>57)&0x7f);
        printf ("\e[38;5;%dm%c", col, gol[ij] ? '*' : ' ');
        if ((ij % N) == N -1) printf ("\n");
    }
    printf("\e[38;5;238;48;5;255m");
    fflush(stdout);
}
//.......................................................................................................................................................
void colorgenes( int cgolg[], int NN2, int colorfunction, int winnr, int nfrstep) {
    uint64_t gene, gdiff, g2c, mask, quad, clone;
    uint64_t *ggol,*ggolg,*ggolb,*ggolr,*ggolgstats;
    int ij,k,j,nfrPlane,nbeven,activity,popcount,labelxy;
    unsigned int d,d0,d1,d2,jper;
    unsigned int color[3],colormax;
    double rescalecolor;
    uint64_t *traceptr;
    int *npopptr;
    const int colpopmax = 225;                                              // it would be useful to bring this parameter up to python
    const double colpoplow = 50.0;
    const double logcolpopmax= log(colpoplow+(double) colpopmax);
    int actmax = 0;                                                         // it would be useful to bring this parameter up to python
    int ancestortypec;
    int it_nbhood2;
    short unsigned int dscale[16] = {0xff,0xcf,0x7f,0x4f,0x2f,0x27,0x1f,0x1d,0x1b,0x19,0x17,0x15,0x14,0x13,0x12,0x11};
    int dx,dy;
    quadnode *q;
    uint64_t qid,ancestor,root;
    // static int numones[16]={0,1,1,2,1,2,2,3,1,2,2,3,2,3,3,4};
    int labelimage(uint64_t hashkeypatt, unsigned int labelimg[], unsigned int label, int offset);
    extern INLINE int log2size(const short unsigned int golpw);
    void get_gliderinfo(uint64_t outgliderinfo[], int narraysize);
    
    if (nfrstep==0) {
        ggol  = gol;
        ggolg = golg;
        ggolb = golb;
        ggolr = golr;
        ggolgstats = golgstats;
    }
    else {
        nfrPlane = (curPlane+nfrstep) % numPlane;
        ggol = planes[nfrPlane];                                            // get planes of gol,golg,golb,golr,golgstats data
        ggolg = planesg[nfrPlane];
        ggolgstats = planesgs[nfrPlane];
        ggolb = planesb[nfrPlane];
        ggolr = planesr[nfrPlane];
    }

    switch (colorfunction) {
      case 0:                                                               // colorfunction based on multiplicative hash
        // see https://stackoverflow.com/questions/6943493/hash-table-with-64-bit-values-as-key/33871291
        for (ij=0; ij<N2; ij++) {
            if (ggol[ij]) {
                gene = ggolg[ij];
                if (gene == 0ull) gene = 11778ull;                          // random color for gene==0
                // mask = (gene * 11400714819323198549ul) >> (64 - 8);      // hash with optimal prime multiplicator down to 8 bits
                // mask = (gene * 11400714819323198549ul) >> (64 - 32);     // hash with optimal prime multiplicator down to 32 bits
                mask = gene * 11400714819323198549ull;
                mask = mask >> (64 - 32);                                   // hash with optimal prime multiplicator down to 32 bits
                if(rulemod & 0x4) mask = (((ij>>log2N)==((N>>1)-(initfield>>1)-1) || ((ij>>log2N)==((N>>1)-(initfield>>1)-2))) && (ij & 0x1) ? 0x00ffffffull : mask);
                mask |= 0x080808ffull;                                      // ensure visible (slightly more pastel) color at risk of improbable redundancy, make alpha opaque
                cgolg[ij] = (int) mask;
            }
            else {
                if(rulemod & 0x4) mask = (((ij>>log2N)==((N>>1)-(initfield>>1)-1) || ((ij>>log2N)==((N>>1)-(initfield>>1)-2))) && (ij & 0x1) ? 0x00ffffffull : 0xffull);
                else mask = 0xffull;
                cgolg[ij] = (int) mask;
            }
        }
        break;
      case 1:                                                               // functional color assignment based on selection model
        for (ij=0; ij<N2; ij++) {
            if (ggol[ij]) {
                gene = ggolg[ij];
                POPCOUNT64C(gene,d);                                        // assigns number of ones in gene to d
                switch (selection) {
                        case 0 : gdiff=(gene>>40); mask = ((((gdiff&0xff)<<16)+(((gdiff>>8)&0xff)<<8)+(gdiff>>16))<<8) + 0xff; break;  // MSByte red shade, next 8 green, next 8 blue: highest=white
                        case 1 : mask = d==64? 0x0000ffff : ((((63-d)<<18)+(abs((int)d-32)<<10)+(d<<2))<<8) + 0xff; break;  // number of ones determine color gradient blue to red (all 1s red)
                        case 2 :
                        case 3 : d = d & 0x3; mask = d==3 ? 0xf0f0f0ff : ((0xff<<(d<<3))<<8)+0xff; break;  // scissors-stone-well-paper: red-green-blue-white
                        case 4 : mask = d < ncoding ? ((0x3f^d)<<19)+0xff : ((64-d < ncoding) ? ((0x3f^(64-d))<<11)+0xff : 0xf0f0f0ff); break; // near 0 green, near 1 red, others white
                        case 5 : mask = d >= 32 ? ((0x3f^(64-d))<<11)+0xff : ((0x3f^d)<<19)+0xff; break;  //predators green, prey red
                        case 6 : g2c = (1ull<<ncoding)-1ull;gdiff = gene^g2c; POPCOUNT64C(gdiff,d2);
                                 mask = d<d2 ? (d<<26)+0xff : (d2<<10)+0xff; break;
                        case 7 : g2c = (gene>>8)&((1ull<<ncoding)-1ull);
                                 gdiff = gene&0xff;POPCOUNT64C(gdiff,d2);d = d2>7? 7 : d2;
                                 mask = g2c ? 0xf0f0f0ff : ((0x1f+(d<<5))<<8)+(((gdiff>>4)&0xf)<<27)+(((gdiff&0xf)<<4)<<16)+0xff; break;
                        case 8 :                                            // colors according to nr of non zero LUT entries d from blue to red in increasing d, except for repselect 7
                                 if (((repscheme>>4) & 0x7) == 7) {d = d & 0x3; mask = d==3 ? 0xf0f0f0ff : ((0xff<<(d<<3))<<8)+0xff;} // repselect 7 : sc-st-we-pa: red-green-blue-white
                                 else if (ncoding==1) {gdiff=gene&0xffff;POPCOUNT64C(gdiff,d);mask = d==16 ? 0x0000ffff : ((((15-d)<<20)+(abs((int)d-8)<<12)+(d<<4))<<8) + 0xff;}
                                 else if (ncoding==2) {PATTERN2_32(gene,0x3,d);mask = d==16 ? 0x00ffffff : ((((15-d)<<20)+(abs((int)d-8)<<12)+(d<<4))<<8) + 0xff;}
                                 else {PATTERN4(gene,0xf,d);mask = d==16 ? 0x00ffffff : ((((15-d)<<20)+(abs((int)d-8)<<12)+(d<<4))<<8) + 0xff;}
                                 break;
                        case 9 : PATTERN4(gene,0x0,d0);PATTERN4(~gene&0x8888888888888888ull,0x8,d1);PATTERN4(gene&0x8888888888888888ull,0x8,d2);
                                 mask = (d2<<26)+((d1-d0)<<10)+0xff; break;
                        case 10: POPCOUNT64C(gene&0x7ffffull,d1); POPCOUNT64C((gene>>32)&0x7ffffull,d2);
                                 mask = (d2<<27)+(d1<<11)+0xff; break;
                        case 11: PATTERN8(gene,0x27,d0);d0=d0>1?2:d0;PATTERN8(gene,0x3f,d1);d1=d1>1?2:d1;d0=d0+d1;d0=d0>2?3:d0; // saturated copy nrs for survival GoL LUTs
                                 PATTERN8(gene,0xbf,d1);d1=d1>2?3:d1; PATTERN8(gene&0x8080808080808080ull,0x80,d2);d2=d2>7?7:d2; // saturated copy nrs for birth GoL and birth LUTs
                                 mask =(((d0<<22)+(d1<<14)+(d2<<5)) << 8) + 0xff; break;
                        case 12: POPCOUNT64C(gene&0xffffffffull,d1); POPCOUNT64C((gene>>32)&0xffffffffull,d2);
                                 mask = (d2<<27)+(d1<<11)+0xff; break;
                        case 13: PATTERN8(gene,0x2f,d0);d0=d0>1?2:d0;PATTERN8(gene,0x3f,d1);d1=d1>1?2:d1;d0=d0+d1;d0=d0>2?3:d0; // saturated copy nrs for survival GoL LUTs
                                 PATTERN8(gene,0xbf,d1);d1=d1>2?3:d1; PATTERN8(gene&0x8080808080808080ull,0x80,d2);d2=d2>7?7:d2; // saturated copy nrs for birth GoL and birth LUTs
                                 mask =(((d0<<22)+(d1<<14)+(d2<<5)) << 8) + 0xff; break;  //Check NYI
                        case 14: POPCOUNT64C(gene&0xffffffffull,d1); POPCOUNT64C((gene>>32)&0xffffffffull,d2);
                                 mask = (d2<<27)+(d1<<11)+0xff; break;
                        case 15: mask = 0xffffffff;  break; // NYI
                        default  : mask = ((d+(d<<6)+(d<<12)+(d<<18))<<8) + 0xff;
                }
                if(ggolgstats[ij]&F_nongolchg) mask = 0x00c0c0ff;       // color states changed by non GoL rule are shown in yellow
                // if(ggolgstats[ij]&F_survmut) mask = 0xff00ffff;         // color states for surviving fresh mutants (non-replicated) purple/pink
                // else if(ggolgstats[ij]&F_notgolrul) mask = 0x00ffffff;  // color states for not GoL rule yellow

                cgolg[ij] = (int) mask;
            }
            else cgolg[ij] = 0xff;
        }
        break;
      case 2:                                                               // colouring by s state and central state (and pattern)
        for (ij=0; ij<N2; ij++) {
            gene = ggolgstats[ij];
            if (ggol[ij]) {
                mask = 0xffull + ((F_s_live&gene) << (16+5));
                cgolg[ij] = (int) mask;
            }
            else if (F_s_live&gene) {
                mask = 0xffull + ((F_s_live&gene) << (24+5));
                cgolg[ij] = (int) mask;
            }
            else cgolg[ij] =0xff;
        }
        break;
      case 3:                                                               // color states report golgstats status of processing
                                                                            // principle is to map bits to one of 3 colours
                                                                            // green: 7: mutation 6:  birth 5: parent
                                                                            // blue:  7: survival or disambig used 6: parentdies 5: normal death
                                                                            // red:   7: parentaldeath 6: nongolchange 5: not gol rule
        for (ij=0; ij<N2; ij++) {
            gene = ggolgstats[ij];
            if (ggol[ij]||((F_death|F_parentaldeath)&gene)) {
                mask = 0xffull;                                             // opaque
                mask |= (((F_parentaldeath&gene) ? 0x80ull : 0ull) + ((F_nongolchg&gene) ? 0x40ull : 0ull) + ((F_notgolrul&gene) ? 0x20ull : 0ull))<<8;  // red
                mask |= (((F_mutation&gene) ? 0x80ull : 0ull) + ((F_birth&gene) ? 0x40ull : 0ull) + ((F_parent&gene) ? 0x20ull : 0ull))<<16;             // green
                mask |= ((((F_survival|F_disambig)&gene) ? 0x80ull : 0ull) + ((F_dummy&gene) ? 0x40ull : 0ull) + ((F_death&gene) ? 0x20ull : 0ull))<<24;           // blue
                cgolg[ij] = (int) (mask|0xffull);
            }
            else cgolg[ij] = 0xff;
        }
        break;
      case 4:                                                               // activities
        if(nbhist==-1) {
            traceptr =&acttrace[0];
            npopptr = &npopulation[0];
        }
        else {
            traceptr = &acttrace1[N2*(nbhist>>1)];
            npopptr = &npopulation1[N*(nbhist>>1)];
        }
        nbeven=(nbhist==-1)?1:1-(nbhist&0x1);
        for (ij=0; ij<N2; ij++) {
            int i,i1;
            i=ij&Nmask;
            i1=nbeven?0:(i<(N>>1)?N>>1:N2-(N>>1));
            gene=traceptr[ij+i1];
            if (gene == rootgene) mask = 0x3f3f3fff;            // grey color for background, all root genes
            else {
                if (gene == 0ull) gene = 11778L;                // random color for gene==0
                mask = gene * 11400714819323198549ul;
                mask = mask >> (64 - 32);                       // hash with optimal prime multiplicator down to 32 bits
                mask |= 0x080808ffull;                          // ensure visible (slightly more pastel) color at risk of improbable redundancy, make alpha opaque
                if(colpopmax && (diagnostics & diag_hash_genes)) {
                    popcount=0;
                    if((genedataptr = (genedata *) hashtable_find(&genetable, gene)) != NULL) popcount = genedataptr->popcount;
                    else fprintf(stderr,"gene not found in colorfunction for activities\n");
                    if(popcount>colpopmax) popcount=colpopmax;
                    colormax=0;
                    for(d=0;d<3;d++) if((color[d]=( (mask>>(8+(d<<3))) & 0xff))>colormax) colormax=color[d];
                    if(popcount)
                        rescalecolor=(log(colpoplow+(double)popcount)/logcolpopmax)*((double)0xff/(double)colormax);
                    else
                        rescalecolor=0.2;
                    for(d=0;d<3;d++) color[d]=(unsigned int) (((double) color[d])*rescalecolor);
                    for(d=0,mask=0xff;d<3;d++) mask |= color[d]<<((d<<3)+8);
                }
            }
            i1=nbeven?0:(i<(N>>1)?N>>1:N-(N>>1));
            i=nbeven?0:(i<(N>>1)?0:N);
            if ((npopptr[i+((ij+i1)&Nmask)]>>log2N)==(N-1-(ij>>log2N))) mask = 0xffffffff;  // overlay plot with trace of density in white (except if pop=N2)
            cgolg[ij]= (int) mask;
        }
        break;
      case 5:                                                   // populations of genes or clones
        if(nbhist==-1) {
            traceptr =&poptrace[0];
            npopptr = &npopulation[0];
        }
        else {
            traceptr = &poptrace1[N2*(nbhist>>1)];
            npopptr = &npopulation1[N*(nbhist>>1)];
        }
        nbeven=(nbhist==-1)?1:1-(nbhist&0x1);
        for (ij=0; ij<N2; ij++) {
            int i,i1;
            i=ij&Nmask;
            i1=nbeven?0:(i<(N>>1)?N>>1:N2-(N>>1));
            gene=traceptr[ij+i1];
            if (gene == rootgene) mask = 0x3f3f3fff;            // grey color for background, all root genes
            else {
                if (gene == 0ull) gene = 11778L;                // random color for gene==0
                mask = gene * 11400714819323198549ul;
                mask = mask >> (64 - 32);                       // hash with optimal prime multiplicator down to 32 bits
                mask |= 0x080808ffull;                          // ensure visible (slightly more pastel) color at risk of improbable redundancy, make alpha opaque
                if(actmax && (diagnostics & diag_hash_genes)) {
                    popcount=0;
                    if((genedataptr = (genedata *) hashtable_find(&genetable, gene)) != NULL) popcount = genedataptr->popcount;
                    else fprintf(stderr,"gene not found in colorfunction for population\n");
                    if(popcount>actmax) popcount=actmax;
                    colormax=0;
                    for(d=0;d<3;d++) if((color[d]=( (mask>>(8+(d<<3))) & 0xff))>colormax) colormax=color[d];
                    rescalecolor=(log((double)popcount)/log((double)actmax))*((double)0xff/(double)colormax);
                    for(d=0;d<3;d++) color[d]=(unsigned int) (((double) color[d])*rescalecolor);
                    for(d=0,mask=0xff;d<3;d++) mask |= color[d]<<((d<<3)+8);
                }
            }
            i1=nbeven?0:(i<(N>>1)?N>>1:N-(N>>1));
            i=nbeven?0:(i<(N>>1)?0:N);
            if ((npopptr[i+((ij+i1)&Nmask)]>>log2N)==(N-1-(ij>>log2N))) mask = 0xffffffff;  // overlay plot with trace of density in white (except if pop=N2)
            cgolg[ij]= (int) mask;
        }
        break;
      case 6:                                                   // genealogies
      case 7:                                                   // genealogies with temporal resolution
        if(colorfunction==6) k=2; else k=0;                             // double row for each ancestral step to make display more readable in colorfunction 6 mode
        if((colorfunction2==-1) || (ancestortype !=2)) ancestortypec = ancestortype;
        else if (winnr) ancestortypec = 1;
        else ancestortypec = 0;
        for (ij=0; ij<N2; ij++) {
            int i=ij&Nmask; int j=ij>>(log2N+k);
            if((ancestortypec==1) && (diagnostics &diag_hash_clones)) {
                clone=clonealogytrace[i+(j<<log2N)];
                if (clone == rootclone)   mask = 0x000000ff;            // black color for root
                else {
                    if (clone == 0ull) clone = 11778L;                  // random color for clone==0
                    mask = clone * 11400714819323198549ul;
                    mask = mask >> (64 - 32);                           // hash with optimal prime multiplicator down to 32 bits
                    mask |= 0x080808ffull; // ensure visible (slightly more pastel) color at risk of improbable redundancy, make alpha opaque
                    if((colorfunction==7) && (diagnostics & diag_hash_clones)) {                         // rescale color brightness by activity/activitymax
                        activity = 0;
                        if((clonedataptr = (clonedata *) hashtable_find(&clonetable, clone)) != NULL) activity = clonedataptr->activity;
                        else fprintf(stderr,"clone not found in colorfunction for clonealogy\n");
                        colormax=0;
                        for(d=0;d<3;d++) if((color[d]=( (mask>>(8+(d<<3))) & 0xff))>colormax) colormax=color[d];
                        rescalecolor=0.25+0.75*((double)(activity*255))/((double)(activitymax*colormax));         // integer version doesn't work
                        for(d=0;d<3;d++) color[d]=(unsigned int) (((double) color[d])*rescalecolor);     // rescale colors by activity/activitymax
                        for(d=0,mask=0xff;d<3;d++) mask |= color[d]<<((d<<3)+8);
                    }
                    else {};    // no rescaling of colours for case 6 (formerly cases 5 and 6)
                }
            }
            else {
                gene=genealogytrace[i+(j<<log2N)];
                if (gene == selectedgene)    mask = 0xffffffff;         // white color for selected gene
                else if (gene == rootgene)   mask = 0x000000ff;         // black color for root
                else if (gene == generepeat) mask = 0x3f3f3fff;         // grey color for repeated gene
                else {
                    if (gene == 0ull) gene = 11778L;                    // random color for gene==0
                    mask = gene * 11400714819323198549ul;
                    mask = mask >> (64 - 32);                           // hash with optimal prime multiplicator down to 32 bits
                    mask |= 0x080808ffull; // ensure visible (slightly more pastel) color at risk of improbable redundancy, make alpha opaque
                    if((colorfunction==7) && (diagnostics & diag_hash_genes)) {                           // rescale color brightness by activity/activitymax
                        activity = 0;
                        if((genedataptr = (genedata *) hashtable_find(&genetable, gene)) != NULL) activity = genedataptr->activity;
                        else fprintf(stderr,"gene not found in colorfunction for genealogy\n");
                        colormax=0;
                        for(d=0;d<3;d++) if((color[d]=( (mask>>(8+(d<<3))) & 0xff))>colormax) colormax=color[d];
                        rescalecolor=0.25+0.75*((double)(activity*255))/((double)(activitymax*colormax));         // integer version doesn't work
                        for(d=0;d<3;d++) color[d]=(unsigned int) (((double) color[d])*rescalecolor);     // rescale colors by activity/activitymax
                        for(d=0,mask=0xff;d<3;d++) mask |= color[d]<<((d<<3)+8);
                    }
                    else {};    // no rescaling of colours for case 6 (formerly cases 5 and 6)
                }
            }
            cgolg[ij]=(int) mask;
        }
        break;
      case 8:                                                           // colorfunction based on packed bit pattern with multiplicative hash, compare with offset
        it_nbhood2 = it_nbhood*it_nbhood+2;
        for (ij=0; ij<N2; ij++) {
                gene = golmix[ij];
                if(offdx==0 && offdy==0 && offdt==0) {
                    for (mask=0,k=0;k<8;k++) {
                        d1 = (gene>>(k<<3))&0xff;                   // no of mismatches : 0:0 black space perfect match is coded as 0xff=255
                        d1 = (d1 >= it_nbhood2-1) ? 0 : (it_nbhood2-1-d1);               // no of matches : 0:0 match set to 0 (no match) here
                        d1 = (d1 < it_nbhood2-1-15) ? 0 : d1-(it_nbhood2-1-15);          // 0 to 15 : perfect match is 15  (4 bits)
                        d1 = (d1==0xf) ? 0x1f : d1;                 // perfect match separated to value 31 (5 bits) for better contrast
                        if(k<3) mask+=d1<<(3+(k<<3));               // perfect match has full intensity colour
                        else if (k==3 && d1==0x1f) mask = (d1<<3)+(d1<<11)+(d1<<19); // the fourth channel has white colour : no others shown
                        else if (k<7 && d1==0x1f) mask+= (d1<<(3+((k-4)<<3)))+(d1<<(3+((k<6?k-3:0)<<3))); // mixed colours for NE SE SW
                        else if(d1==0x1f) mask+= (d1<<3)+(d1<<10)+(d1<<18); // mixed colour for NW
                    }
                    mask = (mask<<8)+0xff;
                }
                else {
                    POPCOUNT64C(gene,d);                            // assigns number of ones in gene to d. This 3 line version for one offset comparison
                    d=(d==64)?0:63-d;
                    mask = (d==63) ? 0xffffffff : ((((d&3)<<22)+(((d>>2)&3)<<14)+(((d>>4)&3)<<6))<<8) + 0xff;
                }

                cgolg[ij] = (int) mask;
        }
        if(info_transfer_h) {                                       // display histograms of glider matching in eight directions N E S W NE SE SW NW
            
            uint64_t *binomialp;
            uint64_t binomial9[11] = {1, 9, 36, 84, 126, 126, 84, 36, 9, 1, 1};
            uint64_t binomial25[27] = {1, 25, 300, 2300, 12650, 53130, 177100, 480700, 1081575, 2042975, 3268760, 4457400, 5200300, 5200300, 4457400, 3268760, 2042975, 1081575, 480700, 177100, 53130, 12650, 2300, 300, 25, 1, 1};
            uint64_t binomial49[51] = {1, 49, 1176, 18424, 211876, 1906884, 13983816, 85900584, 450978066, 2054455634, 8217822536, 29135916264, 92263734836, 262596783764, 675248872536, 1575580702584, 3348108992991, 6499270398159, 11554258485616, 18851684897584, 28277527346376, 39049918716424, 49699896548176, 58343356817424, 63205303218876, 63205303218876, 58343356817424, 49699896548176, 39049918716424, 28277527346376, 18851684897584, 11554258485616, 6499270398159, 3348108992991, 1575580702584, 675248872536, 262596783764, 92263734836, 29135916264, 8217822536, 2054455634, 450978066, 85900584, 13983816, 1906884, 211876, 18424, 1176, 49, 1, 1};
            double ratiomax=0.,ratio;

            if(it_nbhood == 3) {
                binomialp = binomial9;
            }
            else if (it_nbhood == 5) {
                binomialp=binomial25;
            }
            else {
                binomialp=binomial49;
            }

            get_gliderinfo(gliderinfo, it_nbhood2*8);

            for (int i=0; i<it_nbhood2*8; i++) {
                if((i%it_nbhood2) == it_nbhood2-1) {
                    gliderinfo[i]=0;
                    ratio = 0.;
                }
                else ratio = (double) gliderinfo[i] / (double) binomialp[i%it_nbhood2];
                if(ratio > ratiomax) {
                    ratiomax = ratio;
                }
            }
            
            // fprintf(stderr,"ratiomax %g\n",ratiomax);
            
            for (ij=N2>>1;ij<N2;ij++) cgolg[ij] = 0;
            for (int i=0; i<it_nbhood2*8; i++) {
                ratio = (double) gliderinfo[i] / (double) binomialp[i%it_nbhood2];
                if((i%it_nbhood2) == it_nbhood2-1) {
                    for (int j=0;j< (N>>1);j++)
                        cgolg[((N-j)<<log2N)+i] = 0x7f7f7fff;
                }
                else {
                    for (int j=0;  j< (int) ((ratio*(double)(N>>1))/ratiomax)  ;j++)
                        cgolg[((N-j)<<log2N)+i] = 0xff00ff + (((i*0x1f/it_nbhood2))<<24)-(((i*0x1f/it_nbhood2))<<16);
                }
            }
            
        }
        break;
      case 9:                                                       // colorfunction based on unique labelling of separate components in image
          if (diagnostics & diag_component_labels) {
            for (ij=0; ij<NN2; ij++) labelcc[ij]=relabel[label[ij]];    // transfer labels to interactive working area
            if(xdisplay>=0 && ydisplay>=0) {
                if ((labelxy=label[xdisplay+ydisplay*N])) {
                    for (ij=0; ij<NN2; ij++) if (label[ij]==complist[labelxy].label) labelcc[ij]=0xffffffff;    // label chosen component white at its current location
                    d=complist[labelxy].log2n;
                    unsigned int conn = connlists[labelxy];
                    while(conn) {
                        for (ij=0; ij<NN2; ij++) if (oldlabel[ij]==connections[conn].oldlab) {
                            if (labelcc[ij]==0xffffffff) labelcc[ij]=0xfffffffd; // color old connected components overlapping pink
                            else labelcc[ij]=0xfffffffe;                         // color old connected components not overlapping red
                        }
                        conn=connections[conn].next;
                    }
                    for (ij=0; ij<(1<<(d<<1)); ij++) labelcc[(ij&((1<<d)-1)) + (ij>>d)*N] = 0;  // initialize component drawing area to zero in top left corner
                    if (complist[labelxy].quad) labelimage(complist[labelxy].quad, labelcc, 0xfffffffc, 0); // extract this component and label it white
                    else labelimage(complist[labelxy].patt, labelcc, 0xfffffffc, 0);  // component involves patt not quad
                }
            }
     
            if (!gcolors) {
              for (ij=0; ij<NN2; ij++) {
                if (labelcc[ij]) {
                    // quad = (uint64_t) relabel[label[ij]];
                    quad = (uint64_t) labelcc[ij];
                    mask = quad * 11400714819323198549ull;              // map label (quasi-uniformly) into larger unsigned colour space
                    mask = mask >> (64 - 32);                           // hash with optimal prime multiplicator down to 32 bits
                    mask |= 0x080808ffull;                              // ensure visible (slightly more pastel) color at risk of improbable redundancy, make alpha opaque
                    if (noveltyfilter && (diagnostics & diag_hash_patterns)) {
                        quad = complist[label[ij]].quad;
                        if((d=complist[label[ij]].patt)) popcount=smallpatts[d].activity; // number of times small (<= 4x4) pattern encountered previously
                        else if((q = (quadnode *) hashtable_find(&quadtable, quad)) != NULL) // if we reach here, quad should have been stored in hash table
                            popcount=q->activity;                       // number of times large pattern encountered previously (poss. also as part of larger patt)
                        else popcount=1;                                // should never occur, but just in case, assume novel
                        if (popcount>1) mask &= (mask&0x3f3f3fff);      // darken non novel components (currently a little too much)
                    }
                    if (labelcc[ij] == 0xffffffff) mask = 0xffffffffull;      // recolor selected component white
                    else if (labelcc[ij] == 0xfffffffc) mask = 0x00ffffffull; // recolor selected component against black background in corner yellow
                    else if (labelcc[ij] == 0xfffffffd) mask = 0xc0c0ffffull; // recolor connected old components overlapping with selected component pink
                    else if (labelcc[ij] == 0xfffffffe) mask = 0x0000ffffull; // recolor connected old components not overlapping with selected component red
                    cgolg[ij] = (int) mask;
                }
                else cgolg[ij] = 0xff;
              }
            }
            else {                                                      // with gcolors set, color according to computed inherited label colours
              for (ij=0; ij<NN2; ij++) {
                if (label[ij]) cgolg[ij] = (int) rgba(complist[label[ij]].gcolor);
                else cgolg[ij] = 0xff;
              }
            }
          } // diag_component_labels
          break;
      case 10:                                                          //activities for patterns with size weighted colours
        for (ij=0; ij<NN2; ij++) labelcc[ij]=0;
        if(xdisplay>=0 && ydisplay>=0) {
            if ((qid=acttraceq[xdisplay+ydisplay*N])) {
                labelimage(qid, labelcc, 0xfffffffc, 0);                    // extract this component and label it yellow
            }
        }
        for (ij=0; ij<N2; ij++) {
            quad=acttraceq[ij];
            if(labelcc[ij]) mask = 0x00ffffff;
            else if (quad == rootgene) mask = 0x3f3f3fff;               // grey color for background, all root genes
            else {
                if (activity_size_colormode == 0) {
                    mask = quad * 11400714819323198549ul;
                    mask = mask >> (64 - 32);                           // hash with optimal prime multiplicator down to 32 bits
                    mask |= 0x080808ffull;                              // ensure visible (slightly more pastel) color at risk of improbable redundancy, make alpha opaque
                }
                else if (activity_size_colormode == 1) {                // color by log2n, enclosing square size
                    if(acttraceqt[ij] && (q = (quadnode *) hashtable_find(&quadtable, quad)) != NULL) d = log2upper((unsigned int) q->size);
                    else if (!acttraceqt[ij] && quad<65536ull && smallpatts[quad].activity) {
                        d = log2size((short unsigned int) quad);
                    }
                    else {fprintf(stderr,"quad pattern not found in colorfunction for activities\n");d=0;}
                    if(d>log2N){fprintf(stderr,"Error in colorfunction 10, size error %d\n",d);d=0;}
                    setcolor(color,d);
                    mask = (color[0]<<8) | (color[1]<<16) |  (color[2]<<24) | 0xff;
                }
                else {                                                  // color by sqrt of nr of live pixels (up to max value of 255)
                    int popmax = 255;
                    if(acttraceqt[ij] && (q = (quadnode *) hashtable_find(&quadtable, quad)) != NULL) popcount = q->pop1s;
                    else if (!acttraceqt[ij] && quad<65536ull && smallpatts[quad].activity) {POPCOUNT64C((uint64_t) quad,popcount);}
                    else {fprintf(stderr,"quad pattern not found in colorfunction for activities\n");popcount=0;}
                    if (activity_size_colormode == 3) popcount = (int) sqrtupper(popcount);
                    // fprintf(stderr,"step %d ij %d popcount %d\n",totsteps,ij,popcount);
                    if(popcount>popmax) popcount=popmax;
                    color[0]=popcount;color[0]= color[0]>255 ? 255: color[0];
                    color[1]=popcount<<3;color[1]=color[1]>255 ? 255: color[1];
                    color[2]=popcount<<6;color[2]=color[2]>255 ? 255: color[2];
                    for(d=0,mask=0xff;d<3;d++) mask |= color[d]<<((d<<3)+8);
                }
            }
            if ((npopulation[ij&Nmask]>>log2N)==(ij>>log2N)) mask = 0xffffffff;  // overlay plot with trace of density in white (except if pop=N2)
            cgolg[ij]= (int) mask;
        }
        break;
      case 11:                                                          //genealogy based colours of ancestors at genealogycoldepth
        if((colorfunction2==-1) || (ancestortype !=2)) ancestortypec = ancestortype;
        else if (winnr) ancestortypec = 1;                              // ancestortype set to 3: use values 0 and 2 in two windows
        else ancestortypec = 0;
        if (ancestortypec==1) root = rootclone;
        else root = rootgene;

        for (ij=0; ij<N2; ij++) {
            if (gol[ij] && (diagnostics & diag_hash_genes)) {
                gene = (ancestortypec==2) ? ggolb[ij] : ggolg[ij];    // variable gene holds either birthid (clones) or gene at ij;
                ancestor=gene;
                for (j=1;j<=genealogycoldepth;j++) {
                    if (ancestor==root) break;                       // reached root, exit j loop
                    else {
                        gene = ancestor;
                        if (ancestortypec==1) {
                            if((clonedataptr = (clonedata *) hashtable_find(&clonetable, gene)) != NULL) {
                                ancestor=clonedataptr->parentid;
                            }
                            else fprintf(stderr,"ancestor not found in clonealogies\n");
                        }
                        else {
                            if((genedataptr = (genedata *) hashtable_find(&genetable, gene)) != NULL) {
                                ancestor=genedataptr->firstancestor;
                            }
                            else fprintf(stderr,"ancestor not found in genealogies\n");
                        }
                    }
                }
                if (gene == 0ull) gene = 11778ull;                  // random color for gene==0
                mask = gene * 11400714819323198549ull;
                mask = mask >> (64 - 32);                           // hash with optimal prime multiplicator down to 32 bits
                mask |= 0x080808ffull;                              // ensure visible (slightly more pastel) color at risk of improbable redundancy, make alpha opaque
                if((ancestortypec!=1) && (ancestor==rootgene)) mask = 0x3f3f3fff;                // grey color for rootclone
                if((ancestortypec==1) && (ancestor&rootclone)) mask = 0x3f3f3fff;                // grey color for rootclone

                cgolg[ij] = (int) mask;
            }
            else cgolg[ij] = 0xff;
        }
        break;
      case 12:                                                      // colouring based on periodicity of dynamic record of 16 last states in golr for live genes
        for (ij=0; ij<N2; ij++) {
            if (ggol[ij] && (diagnostics & diag_hash_genes)) {
                gene = ggolr[ij];                                   // variable gene holds dynamical record stored in golr
                golr_digest (gene, &d0, &d1, &jper, &dx, &dy);
                if (d0 > 7) mask = 0x080808ffull;                   // dark grey color for no significant periodic match found
                else {
                    mask = 0xffull;
                    if ((d0==d1) && (~gene&0x8ull)) {               // survival with constant mismatch: ie static pattern : colour dark red
                        mask |= 0x3f<<8;
                    }
                    else if (abs(dx)< 0.25*(16-jper) && abs(dy) < 0.25*(16-jper)) {
                        mask |= 0x7f<<8;                            // stationary pattern: slightly less dark red
                    }
                    else {
                        mask |= (jper*(8-d0))<<25;
                        mask |= ((15-jper)*(8-d0))<<17;
                        mask |= dscale[(d0<<1)+1]<<8;
                    }
                }
                cgolg[ij] = (int) mask;
            }
            else cgolg[ij] = 0xff;
        }
        break;
      default:
        fprintf(stderr,"Error: colorfunction %d not supported\n",colorfunction);
        break;
    }
    for (ij=0; ij<N2; ij++) {                                       // convert BGRA format (pygame) to ARGB (PySDL2)
        uint32_t c = cgolg[ij];
        cgolg[ij]=((c&0xff)<<24) | ((c&0xff00)<<8) | ((c&0xff0000)>>8) | ((c&0xff000000)>>24);
    }
}
