//
//  spatial_control.c
//  genelife
//
//  spatial control including vertical scrolling and localized random influx options
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
//
// vscroll              vertical scroll of array in addition to dynamics with wrapping block at top/bottom barrier
// random_influx        random influx of live genes in initialization containing square, in continuous or intermittent modes, with optiional boundary feathering
//----------------------------------------------------------------- spatial control --------------------------------------------------------------------
void v_scroll(uint64_t newgol[],uint64_t newgolg[],uint64_t newgolb[],uint64_t newgolr[]) {
    int ij,scroll_needed;

    scroll_needed = 0;
    for (ij=N2-N;ij<N2;ij++) {                              // clear top row only if needed
        if(newgol[ij]) {
            scroll_needed = 1;
            break;
        }
    }
    last_scrolled = scroll_needed;                          // global variable needed for correct glider detection if scrolling

    for (ij=0;ij<N;ij++) {                                  // delete genes in bottom buffer row
        if(newgol[ij]) {
            newgol[ij]=0ull;
            if(diagnostics & diag_hash_genes)
                hashdeletegene(newgolg[ij],newgolb[ij],"error in v_scroll hashdeletegene call for step %d with gene %"PRIx64"\n");
            newgolg[ij]=gene0;
            newgolb[ij]=0ull;
            newgolr[ij]=0ull;
        }
    }

    if(!scroll_needed) return;

    for (ij=0;ij<N2-N;ij++) {                               // scroll rows down 1 leaving top row intact
        newgol[ij]=newgol[ij+N];
        newgolg[ij]=newgolg[ij+N];
        newgolb[ij]=newgolb[ij+N];
        newgolr[ij]=newgolr[ij+N];
    }
    for (ij=0;ij<N;ij++) {                                  // delete all states and genes in new bottom buffer row
        if(newgol[ij]) {
            newgol[ij]=0ull;
            if(diagnostics & diag_hash_genes)
                hashdeletegene(newgolg[ij],newgolb[ij],"error in v_scroll hashdeletegene call for step %d with gene %"PRIx64"\n");
            newgolg[ij]=gene0;
            newgolb[ij]=0ull;
            newgolr[ij]=0ull;
        }
    }
    for (ij=N2-N;ij<N2;ij++) {                              // clear top row
        if(newgol[ij]) {
            newgol[ij]=0ull;
            if(diagnostics & diag_hash_genes)
                hashdeletegene(newgolg[ij],newgolb[ij],"error in v_scroll hashdeletegene call for step %d with gene %"PRIx64"\n");
            newgolg[ij]=gene0;
            newgolb[ij]=0ull;
            newgolr[ij]=0ull;
        }
    }
}
//.......................................................................................................................................................
void random_influx(uint64_t newgol[],uint64_t newgolg[],uint64_t newgolb[],uint64_t newgolr[]) {
    int Nf,i,j,ij,i0,j0,i1,j1,d;
    uint64_t randnr,mask,parentid;
    static unsigned int rmask = (1 << 15) - 1;
    unsigned int density;
    parentid = ((uint64_t) totsteps)<<32;
    
    if(rbackground) {                                           // homogeneous random background at rate rbackground/32768 per site per frame
        if(randominflux<3) {                                    // only create new genes as perturbation if randominflux<3
            if(randominflux < 2 || selection < 8) {             // random gene if randominflux<2 or selection < 8
                for(ij=0;ij<N2;ij++) {
                    if((rand()&rmask) < rbackground) {          // random event
                        if (!newgol[ij]) {                      // if not live, random genome
                            newgol[ij] = 1ull;
                            RAND128P(randnr);
                            newgolg[ij] = randnr;
                            if(diagnostics & diag_hash_genes) hashaddgene(ij,newgolg[ij],rootgene,newgolb+ij,(parentid+rootclone+ij),0x1ull);
                            newgolr[ij]=0ull;
                        }
                    }
                }
            }
            else if(randominflux==2) {                          // randominflux==2 and selection>=8
                for(ij=0;ij<N2;ij++) {
                    if((rand()&rmask) < rbackground) {
                        if (!newgol[ij]) {                      // if not live, fill with game of life genome
                            newgol[ij] = 1ull;
                            newgolg[ij] = genegol[selection-8]; // GoL gene for particular selection model coding
                            if(diagnostics & diag_hash_genes) hashaddgene(ij,newgolg[ij],rootgene,newgolb+ij,(parentid+rootclone+ij),0x1ull);
                            newgolr[ij]=0ull;
                        }
                    }
                }
            }
        }
        else if (randominflux==3) {                             // deletion perturbations only
            for(ij=0;ij<N2;ij++) {
                if((rand()&rmask) < rbackground) {              // random event
                    if (newgol[ij]) {                           // if live cell, delete gene
                        newgol[ij]=0ull;
                        if(diagnostics & diag_hash_genes) hashdeletegene(newgolg[ij],newgolb[ij],"error in randominflux=3 hashdeletegene call for step %d with gene %"PRIx64"\n");
                        newgolg[ij]=gene0;
                        newgolb[ij]=0ull;
                        newgolr[ij]=0ull;
                    }
                }
            }
        }
        else if (randominflux==4) {                             // deletion perturbations s-dependent
            int s, se, k;
            int nb[8], ij, i, j, jp1, jm1, ip1, im1;
            uint64_t gols, nb1i;
            for(ij=0;ij<N2;ij++) {
                if((rand()&rmask) < rbackground) {              // random event
                    if (newgol[ij]) {                           // if live cell, delete gene
                        // compute s
                        i = ij & Nmask;  j = ij >> log2N;                                          // row & column
                        jp1 = ((j+1) & Nmask)*N; jm1 = ((j-1) & Nmask)*N;                          // toroidal (j+1)*N and (j-1)*N
                        ip1 =  (i+1) & Nmask; im1 =  (i-1) & Nmask;                                // toroidal i+1, i-1
                        nb[0]=jm1+im1; nb[1]=jm1+i; nb[2]=jm1+ip1; nb[3]=j*N+ip1;                  // new order of nbs
                        nb[4]=jp1+ip1; nb[5]=jp1+i; nb[6]=jp1+im1; nb[7]=j*N+im1;
                        for (s=se=0,nb1i=0ull,k=0;k<8;k++) {                                       // packs non-zero nb indices in first up to 8*4 bits
                            gols=newgol[nb[k]];                                                    // whether neighbor is alive
                            s += gols;                                                             // s is number of live nbs
                            se += k&0x1&gols;                                                      // se is number of edge-centred live neighbours (odd k)
                            nb1i = (nb1i << (gols<<2)) + (gols*k);                                 // nb1i is packed list of live neighbour indices
                        }
                        // compute s done
                        if(s>1){
                            newgol[ij]=0ull;
                            if(diagnostics & diag_hash_genes) hashdeletegene(newgolg[ij],newgolb[ij],"error in randominflux=4 hashdeletegene call for step %d with gene %"PRIx64"\n");
                            newgolg[ij]=gene0;
                            newgolb[ij]=0ull;
                            newgolr[ij]=0ull;
                        }
                    }
                }
            }
        }
        return;
    }

    if(randominflux>=2)
        if ((totsteps & 0xf) || randominflux == 3) return;      // only execute remaining code once every 16 time steps and if randominflux!=3

    density = initial1density;
    mask = ~0ull;
    Nf = initfield;
    if (Nf==0 || Nf>N) Nf=N;
    i0 = j0 = (N>>1)-(Nf>>1);

    for (i=0; i<Nf; i++) {
        for (j=0; j<Nf; j++) {
            ij=i0+i+N*(j0+j);
            if(randominflux==2) {                                 // border feathering as well as intermittent every 16 steps
                i1 = i<j ? i : j;                               // swap so that i1<=j1
                j1 = i<j ? j : i;
                d= j1< (Nf>>1) ? i1 : (i1 < Nf-j1 ? i1 : Nf-j1);// find Manhatten distance to border ij1
                density = (d <= 8 ? (initial1density >> (16-(d<<1))) : initial1density);
            }
            if(!newgol[ij]) {
                newgol[ij] = ((rand() & rmask) < density)?1ull:0ull;
 
                if (newgol[ij]) {                               // if live cell, fill with game of life genome or random genome
                    if (selection<8) {
                        RAND128P(randnr);
                        newgolg[ij] = gene0^(randnr&mask);
                    }
                    else {
                        newgolg[ij] = genegol[selection-8];
                    }
                    if(diagnostics & diag_hash_genes) hashaddgene(ij,newgolg[ij],rootgene,newgolb+ij,(parentid+rootclone+ij),0x1ull);
                }
            }
        }
    }
}
