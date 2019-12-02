//
//  pack.c
//  project genelife
//
//  functions to pack, unpack information from the local arrays into integer words to allow parallel processing
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
// pack012neighbors     pack all up to 2nd neighbours in single word
// pack0123neighbors    pack all up to 3rd neighbours in single word
// pack49neighbors      fast routine to pack all up to 1st,2nd,3rd neighbours in single word : order of bits dictated by hierarchical assembly
// pack16neighbors      pack 4x4 blocks in single uint64_t word using 16 bits
// unpack16neighbors    unpack 16 bit word to 4x4 block at offset in full array of labels, marking with chosen label
// log2size             return log2 of linear size of pattern in integer power of 2 for small patterns 0-65535
// pack64neighbors      pack 8x8 blocks in single uint64_t (long) word with 64 bits
// unpack64neighbors    unpack 64 bit word to 8x8 block at offset in full array of labels, marking with chosen label
// compare_neighbors    compare packed pack neighbours with one given x,y shift of arbitrary size
// compare_all_neighbors compare packed pack neighbours with all nearest neighbour x,y shifts
// packandcompare       pack and compare either all 1-shifted 3-neighbourhoods with t=-1 or chosen (dx,dy,dt) 3-neighbourhoods
// golr_digest          digest information in golr (displacement record), extracting min and max mismatch, period, and x,y displacement of period
//----------------------------------------------------------- pack012,3neighbors ------------------------------------------------------------------------
extern INLINE void pack012neighbors(uint64_t gol[],uint64_t golp[]) {              // routine to pack all up to 2nd neighbours in single word
    unsigned int ij,k;
    uint64_t gs;
    int nbx[24] = {-1,0,1,1,1,0,-1,-1,-2,-1,0,1,2,2,2,2,2,1,0,-1,-2,-2,-2,-2};
    int nby[24] = {-1,-1,-1,0,1,1,1,0,-2,-2,-2,-2,-2,-1,0,1,2,2,2,2,2,1,0,-1};

    for (ij=0;ij<N2;ij++)  golp[ij] = gol[ij];       // copy 1 bit gol to golp

    for (ij=0;ij<N2;ij++) {                          // copy up to 2nd neighbours to golp in upper 32-bit word
        for(k=0,gs=0ull;k<24;k++) {
            gs=gol[DELTAXY(ij,nbx[k],nby[k])];
            golp[ij] |= gs<<(k+32);
        }
    }
}
//.......................................................................................................................................................
extern INLINE void pack0123neighbors(uint64_t gol[],uint64_t golp[]) {              // routine to pack all up to 3rd neighbours in single word
    unsigned int ij,k;
    uint64_t gs;
    int nbx[48] = {-1,0,1,1,1,0,-1,-1,-2,-1,0,1,2,2,2,2,2,1,0,-1,-2,-2,-2,-2,-3,-2,-1,0,1,2,3,3,3,3,3,3,3,2,1,0,-1,-2,-3,-3,-3,-3,-3,-3};
    int nby[48] = {-1,-1,-1,0,1,1,1,0,-2,-2,-2,-2,-2,-1,0,1,2,2,2,2,2,1,0,-1,-3,-3,-3,-3,-3,-3,-3,-2,-1,0,1,2,3,3,3,3,3,3,3,2,1,0,-1,-2};

    for (ij=0;ij<N2;ij++) golp[ij] = gol[ij];       // copy 1 bit gol to golp

    for (ij=0;ij<N2;ij++) {
        for(k=0,gs=0ull;k<48;k++) {
            gs=gol[DELTAXY(ij,nbx[k],nby[k])];
            golp[ij] |= gs<<(k+8);                  // copy up to 2nd neighbours to golp in upper 7 bytes of word
        }
    }
}
//.......................................................................................................................................................
extern INLINE void pack49neighbors(uint64_t gol[],uint64_t golp[], int nbhood) {              // routine to pack all up to 3rd neighbours in single word
    unsigned int ij,k;
    int nbx[6] = {1,0,2,0,-4,-4};
    int nby[6] = {0,1,0,2,0,-4};
    /* 48 49 52 53 32 33 36 37                                                   // indexing order of bits in packed 64-bit word for neighborhoods around 0 bit
       50 51 54 55 34 35 38 39                                                   // unused bits are masked out after six step recursive packing
       56 57 60 61 40 41 44 45
       58 59 62 63 42 43 46 47
       16 17 20 21  0  1  4  5
       18 19 22 23  2  3  6  7
       24 25 28 29  8  9 12 13
       26 27 30 31 10 11 14 15 */

    for (ij=0;ij<N2;ij++) golp[ij] = gol[ij];                                     // copy 1 bit gol to golp
    for(k=0;k<6;k++)                                                              // hierarchical bit copy and swap
        for (ij=0;ij<N2;ij++)
             golp[ij] |= golp[DELTAXY(ij,nbx[k],nby[k])]<<(1<<k);                 // 8x8 packed arrays
    if(nbhood == 7)                                                               // masks out 15 values in top row and left column to give 7x7 neighbourhoods
        for (ij=0;ij<N2;ij++) golp[ij] = golp[ij]&0xfac8ffccfafaffffull;          // mask removes bit numbers 16,18,24,26,32,33,36,37,48,49,50,52,53,56,58
    else if (nbhood == 5)                                                         // masks in 25 values to give 5x5 neighbourhoods
        for (ij=0;ij<N2;ij++) golp[ij] = golp[ij]&0xf0005f0030f0135full;          // mask in bit nrs 0,1,2,3,4,6,8,9,12,20,21,22,23,28,29,40,41,42,43,44,46,60,61,62,63
    else                                                                          // masks in 9 values to give 3x3 neighbourhoods
        for (ij=0;ij<N2;ij++) golp[ij] = golp[ij]&0x80000c0000a0000full;          // mask in bit nrs 0,1,2,3,21,23,42,43,63

}
//.......................................................................................................................................................
extern INLINE short unsigned int pack16neighbors(uint64_t wgol[], int log2n) {    // routine to pack up to 4x4 subarray of binary square array wgol (nxn) into single word
    uint64_t golp;                                                      // assuming side length of square is n (power of 2)
    if (log2n==0) return((short unsigned int) wgol[0]);
    else if (log2n==1) return((short unsigned int) (wgol[0]+(wgol[1]<<1)+(wgol[2]<<2)+(wgol[3]<<3)));
    else if (log2n==2) {
        golp=0ull;
        golp |=  (wgol[0]+(wgol[1]<<1)+(wgol[4]<<2)+(wgol[5]<<3));
        golp |= ((wgol[2]+(wgol[3]<<1)+(wgol[6]<<2)+(wgol[7]<<3))<<4);
        golp |= ((wgol[8]+(wgol[9]<<1)+(wgol[12]<<2)+(wgol[13]<<3))<<8);
        golp |= ((wgol[10]+(wgol[11]<<1)+(wgol[14]<<2)+(wgol[15]<<3))<<12);
        return((short unsigned int) golp);
    }
    else {fprintf(stderr,"pack16neighbours called with not permitted value of log2n %d\n",log2n);return(0);}
}
//.......................................................................................................................................................
extern INLINE void unpack16neighbors(const short unsigned golpw, unsigned int labelimg[],const unsigned int label,const int offset){
    int k,ij;
    if (golpw < 2) labelimg[0] = golpw ? label : 0;
    else if (golpw < 16) {
        for(k=0;k<4;k++) {
            ij = DELTAXY(offset,k&0x1,k>>1);
            labelimg[ij] = (golpw>>k)&0x1 ? label : 0;
        }
    }
    else {
        for(k=0;k<16;k++) {
            int k1 = k&0x1; int k2 = (k>>1)&0x1;
            ij = DELTAXY(offset,k1+(((k>>2)&0x1)<<1),k2+((k>>3)<<1));
            labelimg[ij] = (golpw>>k)&0x1 ? label : 0;
        }
    }
}
//.......................................................................................................................................................
extern INLINE int log2size(const short unsigned int golpw) {
    if (golpw < 2) return 0;
    else if (golpw < 16) return 1;
    else return 2;
}
//.......................................................................................................................................................
extern INLINE void pack64neighbors(uint64_t gol[],uint64_t golp[],int log2n) {    // routine to pack 8x8 subarrays of full binary array gol into single words
    int n = (1 << log2n);
    int ij,ij1,k;                                                                 // assuming golp length >= (n*n)>>2^6
    int n2 = n*n;

    if (log2n<3) {
        fprintf(stderr,"Error trying to pack too small an array into 64 bit words: need >= 8x8 have %d x %d\n",n,n);
        golp[0]=0;
        return;
    }
    for (ij1=0;ij1<(n2>>6);ij1++) golp[ij1]=0ull;
    for(k=0;k<64;k++) {
        int kx = ((k&0xf)&3)+(((k>>4)&1)<<2); int ky = ((k&0xf)>>2)+(((k>>4)>>1)<<2);
        for (ij1=0,ij=0;ij<n2;ij+=8) {
            // golp[ij1++] |= gol[(ij &(n-1))+(k&0x7)+n*((ij>>log2n)+(k>>3))]<<k;// blocked as 8*8 not compatible with smallpatts
            golp[ij1++] |= gol[(ij &(n-1))+ kx +n*((ij>>log2n)+ky)]<<k;          // blocked as 4*4*4
            ij+= ((ij+8)&(n-1)) ? 0 : (8-1)*n;
        }
    }
}
//.......................................................................................................................................................
extern INLINE void unpack64neighbors(const uint64_t golpw, unsigned int labelimg[], const unsigned int label, const int offset){
    int k,ij;                                                                   // only unpacks one word
    for(k=0;k<64;k++) {                                                         // bits blocked as 4*4*4 so that 1st 16 bits are nw 4*4 block
        int k1 = k&0xf; int k2 = k>>4;
        ij = DELTAXY(offset,(k1&3)+((k2&1)<<2),(k1>>2)+((k2>>1)<<2));
        labelimg[ij] = (golpw>>k)&0x1 ? label : 0;
    }
}
//.......................................................................................................................................................
extern INLINE void compare_neighbors(uint64_t a[],uint64_t b[], int dx, int dy) {  // routine to compare packed pack neighbours with shift, result in a
    unsigned int ij,ijs,scrollN;
    uint64_t bij;

    if (last_scrolled) scrollN = N;
    else scrollN = 0;

    for (ij=0;ij<N2;ij++) {
        ijs=(ij-scrollN)&N2mask;
        bij=b[DELTAXY(ijs,dx,dy)];
        a[ij] = (a[ij]|bij) ? a[ij]^bij : rootgene;
    }
}
//.......................................................................................................................................................
extern INLINE void compare_all_neighbors(uint64_t a[],uint64_t b[]) {  // routine to compare packed pack neighbours with shift, result in a
    unsigned int ij,ijs,scrollN;
    unsigned int d, k;
    uint64_t aij,bijk;
    int nbx[8] = {0,1,0,-1,1,1,-1,-1};   // N E S W NE SE SW NW
    int nby[8] = {-1,0,1,0,-1,1,1,-1};

    if (last_scrolled) scrollN = N;
    else scrollN = 0;
    
    for (ij=0;ij<N2;ij++) {
        ijs=(ij-scrollN)&N2mask;
        aij = a[ij];
        for (a[ij]=0ull,k=0;k<8;k++) {
            bijk=b[DELTAXY(ijs,nbx[k],nby[k])];
            POPCOUNT64C((aij^bijk),d);
            d = (aij&&bijk) ? d : 0xff;
            a[ij]|=((uint64_t) d)<<(k<<3);
        }
    }
}
//.......................................................................................................................................................
extern INLINE void packandcompare(uint64_t newgol[],uint64_t working[],uint64_t golmix[]) {
    pack49neighbors(newgol,working,it_nbhood);                      // 3x3,5x5 or 7x7 packed newgol values in working
    if(offdx==0 && offdy==0 && offdt==0) {
        pack49neighbors(gol,golmix,it_nbhood);                      // 3x3,5x5 or 7x7 packed gol values in golmix
        compare_all_neighbors(golmix,working);                      // compare all 8 directions N E S W NE SE SW NW;
    }                                                               // output=golmix will contain packed numbers of 7x7 differences for all 8 directions
    else {
        if (offdt<=-maxPlane) offdt=-maxPlane;
        if(offdt>0) offdt = 0;
        pack49neighbors(planesg[(newPlane-offdt)%maxPlane],golmix,it_nbhood);
        compare_neighbors(golmix,working,offdx,offdy);              // compare with a single direction (north) for gliders
    }
}
//--------------------------------------------------------------digest golr information ------------------------------------------------------------------
extern INLINE void golr_digest (uint64_t golr, unsigned int *mismatchmin, unsigned int *mismatchmax, unsigned int *period, int *pershx, int *pershy) {
                                                        // extract optimal period match: min, max numbers of mismatches,period and x,y offsets for opt. match period
    uint64_t gdiff;
    int j,k,d,d0,d1,dx,dy;
    unsigned int jper;
    int nbx[8] = {-1,0,1,1,1,0,-1,-1};
    int nby[8] = {-1,-1,-1,0,1,1,1,0};
    int dsx,dsy;

    
    gdiff = golr;                                       // variable golr holds dynamical record
    d0 = 16; d1 = 0;                                    // min,max number of mismatches

    jper = 0;
    for (j=0;j<15;j++) {
        gdiff = (gdiff>>4)|((gdiff&0xfull)<<60);        // rotate record cyclically by one time step
        POP4COUNT64C((golr^gdiff),d);                   // number of difference positions between gene and gdiff
        if(d<d0) {
            d0=d;
            jper = j;
        }
        if(d>d1) d1=d;
    }
    *mismatchmin = d0;
    *mismatchmax = d1;
    *period = jper;
    
    dsx = dsy = dx = dy = 0;
    for (j=0;j<16;j++) {
        gdiff = (golr>>(j<<2))&0xfull;
        if(gdiff&0x8ull) {
            k = ((unsigned int) gdiff) & 0x7;
            dx += nbx[k];
            dy += nby[k];
        }
        if (j>=jper) {
            dsx+=dx;
            dsy+=dy;
            dx = dy = 0;
        }
    }
                                                         // average periodic displacements are now dsx/(16-jper), dsy/(16-jper)
    *pershx = dsx;
    *pershy = dsy;
}
