//
//  select_nbs.c
//  project genelife
//  routines associated with the selection of a particular neighbouring gene
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
// selectone_of_2       select one (or none) of two genes based on selection model parameter selection :  returns birth and newgene
// selectone_of_s       select one (or none) of s genes based on selection model parameter selection :  returns birth and newgene
// selectone_nbs        select one of two genes based on pattern of their live 2nd shell neighbours and their genetic encoding
// selectdifft1         select the gene at the single active neighbour position : algorithm could be optimized
// selectdifft2         select the right or left of two genes bunched with least number of empty genes between them
// selectdifft3         select the unique most different (by symmetry) of three live neighbours or first one in canonical rotation
// selectdifft4         select the most central (left) of four live neighbours or first one in canonical rotation
// selectdifft5         select the most central (left) of five live neighbours or first one in canonical rotation
// selectdifft6         select the most central (left) of six live neighbours or first one in canonical rotation
// selectdifft7         select the most central (left) of seven live neighbours or first one in canonical rotation
// selectdifft          select the most central (left) of sum live neighbours or first one in canonical rotation : calls 1 of selectdifft1-7
// analyze_nbs          analyze neighbour states, calculating s and packed set of live indices
// disambiguate         disambiguate the cases where the canonical rotation does not uniquely identify a pattern start point : 1 of 8 methods
//------------------------------------------------------------- selectone -------------------------------------------------------------------------------
extern INLINE void selectone_of_2(int s, uint64_t nb2i, int nb[], uint64_t golg[], uint64_t golb[],uint64_t * birth, uint64_t *newgene, uint64_t *parentid, unsigned int *kch) {
// birth is returned 1 if ancestors satisfy selection condition. Selection of which of two genes to copy is newgene. Non-random result.
    unsigned int k,d0,d1,d2,d3,dd,swap;                  // number of ones in various gene combinations
    uint64_t livegenes[2],gdiff,gdiff0,gdiff1;           // various gene combinations
    uint64_t gene2centre;                                // gene function centres in sequence space
    int g0011,g0110,prey,prey2,ijanc[2];

    for(k=0;k<2;k++) livegenes[k] = golg[ijanc[k]=nb[(nb2i>>(k<<2))&0x7]];
    POPCOUNT64C(livegenes[0],d0);
    POPCOUNT64C(livegenes[1],d1);
    switch (selection) {
        case 0:                                          // integer value of sequence as fitness
            *birth = (livegenes[0]^livegenes[1]) ? 1ull: 0ull; // birth condition is two genes different
            *newgene = livegenes[0]>livegenes[1] ?  livegenes[0] : livegenes[1]; // choose one with larger gene to replicate
            *parentid = (livegenes[0]>livegenes[1]) ?  golb[ijanc[0]] : golb[ijanc[1]];
            break;

        case 1:                                          // number of ones in sequence as fitness
            *birth = (d0^d1) ? 1ull: 0ull;               // birth condition is two genes different in number of ones
            *newgene = (d0>d1) ? livegenes[0] : livegenes[1];
            *parentid = (d0>d1) ? golb[ijanc[0]] : golb[ijanc[1]];
            break;
        case 2:                                          // scissors-stone-well-paper game on number ones mod 4
                                                         // scissors 0 stone 1 well 2 paper 3
                                                         // exception to numerical order: sc>pa
            d0=d0&0x3;d1=d1&0x3;
            *birth = (d0^d1) ? 1ull: 0ull;               // birth if 2 genes differ mod 4 in number of ones
            if(*birth) {
                swap = 0;
                if (d0>d1) {                             // swap d0 and d1 so that smaller one comes first
                    dd = d0; d0 = d1; d1 = dd; swap = 1;
                }
                *newgene = (d0==0 && d1==3) ? livegenes[swap^0] : livegenes[swap^1];
                *parentid = (d0==0 && d1==3) ? golb[ijanc[swap^0]] : golb[ijanc[swap^1]];
            }
            else *newgene = 0ull;
            break;

        case 3:                                          // 4 color game (next color wins) on number of ones mod 4
                                                         // red 0 green 1 blue 2 white 3
                                                         // exception to numerical order: 0>3 birth only if diff=1
            d0=d0&0x3;
            d1=d1&0x3;
            *birth = ((d0^d1)==1ull) ? 1ull: 0ull;       // birth if 2 genes differ by 1 mod 4 in number of ones
            if(*birth) {
                swap = 0;
                if (d0>d1) {                             // swap d0 and d1 so that smaller one comes first
                    dd = d0;
                    d0 = d1;
                    d1 = dd;
                    swap = 1;
                }
                *newgene = (d0==0 && d1==3) ? livegenes[swap^0] : livegenes[swap^1];
                *parentid = (d0==0 && d1==3) ? golb[ijanc[swap^0]] : golb[ijanc[swap^1]];
            }
            else *newgene = 0ull;
            break;
        case 4:                                          // birth if 2 genes cooperate : closer to all 0, all 1 targets than ncoding and closer to each other than 64-ncoding
            gdiff=livegenes[0]^livegenes[1];
            POPCOUNT64C(gdiff,dd);
            *birth = dd>0 && dd<(64-ncoding) && ((d0<ncoding && d1>64-ncoding) || (d1<ncoding && d0>64-ncoding))  ? 1ull: 0ull; // birth if 2 genes close enough to targets
            if (d0<ncoding) {if(d0>64-d1) swap=1;else swap=0;}                  // need 64-d1 != d0 to avoid asymmetry in direction               && (64-d1 != d0)
            else {if(64-d0>=d1) swap=1; else swap=0;}
            *newgene = livegenes[swap];
            *parentid = golb[ijanc[swap]];
            break;
        case 5:                                          // predator prey model: prey-prey evolves to all 0, predator to complement of prey
            gdiff=livegenes[0]^livegenes[1];
            gdiff1=livegenes[0]^(~livegenes[1]);
            POPCOUNT64C(gdiff1,dd);
            prey = (d0<32) || (d1<32);                   // prey present : newgene is one with less ones, 1 prey : predator wins
            prey2 = (d0<32) && (d1<32);                  // 2 prey : newgene is one with less ones, 1 prey : predator wins
            // *birth = (gdiff && prey && dd<ncoding) ? 1ull: 0ull;  // birth if different and >=1 prey and close enough match)
            *birth = (gdiff && prey2) || (prey && (!prey2) && (dd<ncoding)) ? 1ull: 0ull; // birth if different and >=1 prey and close enough match)
            *newgene = (prey2 ? ((d0<d1) ? livegenes[0] : livegenes[1]) : ((d0<32) ? livegenes[1] : livegenes[0]));
            *parentid = (prey2 ? ((d0<d1) ? golb[ijanc[0]] : golb[ijanc[1]]) : ((d0<32) ? golb[ijanc[1]] : golb[ijanc[0]]));
            break;
        case 6:                                         // birth if 2 genes differently functional (Not Yet Working)
                                                                            // 1st target is WLOG the all 0 sequence, d0 and d1 gives distances of two genes
            if(ncoding<64) gene2centre = (1ull<<ncoding)-1ull;              // first ncoding 1s in this sequence
            else gene2centre = ~0ull;                                       // otherwise 2nd target is all 64 ones
            gdiff  = livegenes[0]^livegenes[1];                             // difference between two genes
            gdiff0 = livegenes[0]^gene2centre;                              // difference of gene 0 to 2nd target
            gdiff1 = livegenes[1]^gene2centre;                              // difference of gene 1 to 2nd target
            POPCOUNT64C(gdiff,dd);                                          // dd is distance between genes
            POPCOUNT64C(gdiff0,d2);                                         // d2 is distance of gene 0 from 2nd target
            POPCOUNT64C(gdiff1,d3);                                         // d3 is distance of gene 1 from 2nd target
            g0011 = d0<dd && d3<dd;                      // gene 0 closer to 1st target and gene 1 closer to 2nd target than they are from each other
            g0110 = d2<dd && d1<dd;                      // gene 0 closer to 2nd target and gene 1 closer to 1st target than they are from each other
            *birth = (g0011 && (d0!=d3)) != (g0110 && (d2!=d1))  ? 1ull: 0ull; // birth if 2 genes closer to two different targets than to each other
            *newgene = (g0011 && (d0!=d3)) ? ((d0<d3) ? livegenes[0] : livegenes[1]) : ((d2<d1) ? livegenes[0] : livegenes[1]);
            *parentid = (g0011 && (d0!=d3)) ? ((d0<d3) ? golb[ijanc[0]] : golb[ijanc[1]]) : ((d2<d1) ? golb[ijanc[0]] : golb[ijanc[1]]);
            break;
        case 7:                                          // neutral selection but selective birth only occurs if two chosen sequences are same (different) (NB uses RNG)
                                                         // note that birth is already suppressed for 3 identical live nbs unless enforcebirth-3 bit on
                                                         // we want to use this routine only for two live nbs and only when genes not same
            *birth = livegenes[0]^livegenes[1] ? 1ull: 0ull;
            *newgene = golg[nb[*kch]];
            *parentid = golb[nb[*kch]];
            break;
        // cases 8+: this subroutine is not used
        default:
            fprintf(stderr,"Error: two live gene fitness value %d is not implemented\n",selection);
            exit(1);
    }
}
//.......................................................................................................................................................
extern INLINE int selectone_of_s(unsigned int *kch, int s, uint64_t nb1i, int nb[], uint64_t golg[], uint64_t golb[], uint64_t golr[], uint64_t *birth, uint64_t *newgene, uint64_t *parentid, uint64_t *nbmask, int ij) {
// result is number of equally fit best neighbours that could be the ancestor (0 if no birth allowed, 1 if unique) and list of these neighbour indices
// birth is returned 1 if ancestors satisfy selection condition. Selection of which of genes to copy is newgene. Non-random result.
    unsigned int k,nbest,ijanc[8],kchs[8];                // index for neighbours and number in best fitness category and ij, kch for up to 8 possible ancestors
    unsigned int d[8],dmax[8],p[8],dS,dB,d0,d1,d2;         // d number of ones or mismatches, dmx max mismatches,p period, + other variables
    int psx[8],psy[8];                                    // periodic shift in x and y for optimal period of recorded displacements
    unsigned int scores[8];                               // cumulative scores for pairwise games of individual livegenes (used case repselect == 7)
    uint64_t livegenes[8],gdiff,extremval,bestnbmask,birthid;
    int extreme1,extreme2;
    unsigned int repselect = (repscheme & R_47_repselect)>>4; //

    birthid = (uint64_t) totsteps;
    birthid = (birthid << 32)+rootclone+ij;
    if(s==0) {*birth = 1ull;  *newgene = genegol[selection-8]; *parentid = birthid; *kch = 0; return(1);}

    for(k=0;k<s;k++) {
        kchs[k]=(nb1i>>(k<<2))&0x7;
        ijanc[k] = nb[kchs[k]];
        livegenes[k] = golg[ijanc[k]];
        if(repselect<8) {POPCOUNT64C(livegenes[k],d[k]);}
        else golr_digest (golr[ijanc[k]], d+k, dmax+k, p+k, psx+k, psy+k);
    }

    switch (repselect) {
        case 0:                                          // integer value of sequence as fitness : smallest (case 0) or largest (case 1)
        case 1:
            if(repselect&0x1) for(extremval= 0ull,k=0;k<s;k++) extremval = (livegenes[k] >= extremval ? livegenes[k] : extremval); // find value of fittest gene max value
            else              for(extremval=~0ull,k=0;k<s;k++) extremval = (livegenes[k] <= extremval ? livegenes[k] : extremval); // find value of fittest gene min value
            for(bestnbmask=0ull,nbest=0,k=0;k<s;k++) bestnbmask |= ((livegenes[k]==extremval) ? 1ull<< (k+0*nbest++) : 0ull); // find set of genes with equal best value
            *birth = (nbest>0) ? 1ull: 0ull;             // birth condition may include later that genes not all same
            for(k=0;k<s;k++) if((bestnbmask>>k)&0x1) break;   // execute loop until first optimal live neighbour found at k (a one)
            if (k==s) {k=0;*birth = 0ull;}               // in case no genes with best value, no birth, avoid k being out of bounds below
            *newgene = livegenes[k&0x7];                 // choose first of selected set to replicate (can make positional dependent choice instead externally)
            *parentid=golb[ijanc[k&0x7]];
            *kch = kchs[k];
            break;
        case 2:
        case 3:                                          // number of ones in sequence as fitness : smallest (case 2) or largest (case 3)
            if(repselect&0x1) for(extremval= 0ull,k=0;k<s;k++) extremval = (d[k]>= extremval ? d[k] : extremval); // find value of fittest gene max nr ones
            else              for(extremval=~0ull,k=0;k<s;k++) extremval = (d[k]<= extremval ? d[k] : extremval); // find value of fittest gene min nr ones
            for(bestnbmask=0ull,nbest=0,k=0;k<s;k++) bestnbmask |= (d[k]==extremval ? 1ull<<(k+0*nbest++) : 0ull); // find set of genes with equal best value
            *birth = ((nbest>0) ? 1ull: 0ull);           // birth condition may include later that genes not all same
            for(k=0;k<s;k++) if((bestnbmask>>k)&0x1) break;
            if (k==s) {k=0;*birth = 0ull;}               // in case no genes with best value, no birth, avoid k being out of bounds below
            *newgene = livegenes[k&0x7];                 // choose first of selected set to replicate (can make positional dependent choice instead externally)
            *parentid=golb[ijanc[k&0x7]];
            *kch = kchs[k];
            break;
        case 4:                                          // neutral selection
            for(bestnbmask=0ull,nbest=0,k=0;k<s;k++) bestnbmask |= 1ull<<(k+0*nbest++); // find set of genes with equal best value : in this case all of them
            *birth = 1ull;                               // birth condition is always true
            for(k=0;k<s;k++) if((bestnbmask>>k)&0x1) break;
            if (k==s) {k=0;*birth = 0ull;}               // in case no genes with best value, no birth, avoid k being out of bounds below
            *newgene = livegenes[k&0x7];                 // choose first of selected set to replicate (can make positional dependent choice instead externally)
            *parentid=golb[ijanc[k&0x7]];
            *kch = kchs[k];
            break;
        case 5:                                          // neutral selection but birth only occurs if some sequences are different
            for(gdiff=0ull,k=1;k<s;k++) if ((gdiff=livegenes[0]^livegenes[k])) break; // test whether all genes the same
            if (gdiff) for(bestnbmask=0ull,nbest=0,k=0;k<s;k++) bestnbmask |= 1ull<<(k+0*nbest++); // find set of genes with equal best value : in this case all of them
            else { nbest = 0; bestnbmask = 0ull;}
            *birth = gdiff ? 1ull : 0ull;                // birth condition is genes not all same
            for(k=0;k<s;k++) if((bestnbmask>>k)&0x1) break;
            if (k==s) {k=0;*birth = 0ull;}               // in case no genes with best value, no birth, avoid k being out of bounds below
            *newgene = livegenes[k&0x7];                 // choose first of selected set to replicate (can make positional dependent choice instead externally)
            *parentid=golb[ijanc[k&0x7]];
            *kch = kchs[k];
            break;
        case 6:                                          // penalize genes by a cost for the numebr of LUT entries realized for survival and/or birth
          switch(selection) {
            case 8:                                      // totalistic lut penalty of gene in fixed length encoding : first 8 bits survival, next 8 bits birth
                for(k=0;k<s;k++) {
                    // livegenes[k] = golg[nb[(nb1i>>(k<<2))&0x7]];
                    POPCOUNT64C(livegenes[k]&0xffull,dS);
                    POPCOUNT64C((livegenes[k]>>8)&0xffull,dB);
                    d[k]=24-dS-(dB<<1);
                }
                break;
            case 9:                                      // totalistic lut penalty of gene in variable length encoding : 4 bit patterns S 0xxx  B 1xxx
                for(k=0;k<s;k++) {
                    //livegenes[k] = golg[nb[(nb1i>>(k<<2))&0x7]];
                    PATTERN4(livegenes[k],0x0,d0)
                    PATTERN4(~livegenes[k]&0x8888888888888888ull,0x8,d1)
                    PATTERN4(livegenes[k]&0x8888888888888888ull,0x8,d2)
                    d[k]=48-(d1-d0)-2*d2;
                }
                break;
            case 10:
                for(k=0;k<s;k++) {
                    // livegenes[k] = golg[nb[(nb1i>>(k<<2))&0x7]];
                    POPCOUNT64C(livegenes[k]&0x7fffffull,dS);        // 23 coding bits for survival
                    POPCOUNT64C((livegenes[k]>>32)&0x7fffffull,dB);  // 23 coding bits for birth
                    d[k]=57-dS-(dB<<1);
                }
                break;
            case 11:                                     // dist. lut penalty of gene in variable length encoding : 8 bit patterns S 0xxxrrrr  B 1xxxrrrr
                for(k=0;k<s;k++) {
                    // livegenes[k] = golg[nb[(nb1i>>(k<<2))&0x7]];
                    PATTERN8(livegenes[k],0x00,d0)
                    PATTERN8(~livegenes[k]&0x8080808080808080ull,0x80,d1)
                    PATTERN8(livegenes[k]&0x8080808080808080ull,0x80,d2)
                    d[k]=24-(d1-d0)-2*d2;
                }
                break;
            case 12:                                     // canon. lut penalty in fixed length encoding 32 bit survival 32 bit birth
            case 14:
                 for(k=0;k<s;k++) {
                    livegenes[k] = golg[nb[(nb1i>>(k<<2))&0x7]];
                    POPCOUNT64C(livegenes[k]&0xffffffffull,dS);        // 32 coding bits for survival
                    POPCOUNT64C((livegenes[k]>>32)&0xffffffffull,dB);  // 32 coding bits for birth
                    d[k]=96-dS-(dB<<1);
                }
                break;
            case 13:                                    // canon. lut penalty of gene in variable length encoding : 8 bit patterns S 0xxxwrrr  B 1xxxwrrr
                for(k=0;k<s;k++) {                      // could be made more accurate by accounting for upper pair bit states for rest of 5-subset too
                    // livegenes[k] = golg[nb[(nb1i>>(k<<2))&0x7]];
                    PATTERN8(livegenes[k]&0xffffffffffff,0x00,d0)
                    PATTERN8(~livegenes[k]&0x808080808080ull,0x80,d1)
                    PATTERN8(livegenes[k]&0x808080808080ull,0x80,d2)
                    d[k]=24-(d1-d0)-2*d2;
                }
                break;
            case 15:                                    // dist. or canon. lut penalty of gene in variable length encoding : 8 bit patterns S 0xxxwwrr  B 1xxxwwrr
                for(k=0;k<s;k++) {                      // could be made more accurate by accounting for upper quartet bit states for rest of 6-subset too
                    // livegenes[k] = golg[nb[(nb1i>>(k<<2))&0x7]];
                    PATTERN8(livegenes[k]&0xffffffffff,0x00,d0)
                    PATTERN8(~livegenes[k]&0x8080808080ull,0x80,d1)
                    PATTERN8(livegenes[k]&0x8080808080ull,0x80,d2)
                    d[k]=24-(d1-d0)-2*d2;
                }
                break;
            default:
                fprintf(stderr,"Error: selectone of s for selection %d is not implemented\n",selection);
                exit(1);
          }
          for(extremval=0ull,k=0;k<s;k++) extremval = d[k]>= extremval ? d[k] : extremval; // find value of fittest gene
          for(bestnbmask=0ull,nbest=0,k=0;k<s;k++) bestnbmask |= (d[k]==extremval) ? 1ull<<(k+0*nbest++) : 0ull; // find set of genes with equal best value
          *birth = 1ull;                                // birth condition is always met since extremval set is always > 0
          for(k=0;k<s;k++) if((bestnbmask>>k)&0x1) break;
          if (k==s) {k=0;*birth = 0ull;}                // in case no genes with best value, no birth, avoid k being out of bounds below
          // if(k==s) fprintf(stderr,"Error in selectone of s case 6: k>=s (%d > %d)\n",k,s);
          *newgene = livegenes[k&0x7];                  // choose first of selected set to replicate (can make positional dependent choice instead externally)
          *parentid=golb[ijanc[k&0x7]];
          *kch = kchs[k];
          break;
        case 7:                                         // scissors-stone-well-paper game on number ones mod 4, scissors 0 stone 1 well 2 paper 3
            for (k=0;k<s;k++) {                         // exception to numerical order: sc>pa. Do tournament with all others for each livegene and score
                scores[k]=0; d0 = d[k]&0x3;
                for (int k1=0;k1<s;k1++) {
                    d1 = d[k1]&0x3;                     // if k==k1 there is zero score in next line so allow it
                    scores[k] += ( (d0 > d1 && !(d0 == 3 && d1 == 0)) || (d0 == 0 && d1 == 3) ) ? 1 : 0;
                }
            }
            for(extremval= 0ull,k=0;k<s;k++) extremval = (scores[k]>= extremval ? scores[k] : extremval); // find value of fittest gene, best score in tournament
            if (extremval) {
                for(bestnbmask=0ull,nbest=0,k=0;k<s;k++) bestnbmask |= (scores[k]==extremval ? 1ull<<(k+0*nbest++) : 0ull); // find set of genes with equal best value
            }
            else {   // no positive scores, no birth
                bestnbmask=0ull;nbest=0;
            }
            *birth = nbest ? 1ull: 0ull;               // birth condition
            if (s==3) {                                // included for compatibility with update_23() for enforcebirth off
                for(d0=k=0;k<s;k++) for (int k1=0;k1<k;k1++) d0 |= livegenes[k]^livegenes[k1] ? 1 : 0; //  d0 = 1 if  genes are not all the same
                *birth = d0 ? 1ull: 0ull;              // birth condition modified to only be true if some genes difft
                nbest = d0 ? 3 : 0;
                bestnbmask = d0 ? 0x7ull : 0;
            }
            for(k=0;k<s;k++) if((bestnbmask>>k)&0x1) break;
            if (k==s) {k=0;*birth = 0ull;}             // in case no genes with best value, no birth, avoid k being out of bounds below

            *newgene = livegenes[k&0x7];               // choose first of selected set to replicate (can make positional dependent choice instead externally)
            *parentid=golb[ijanc[k&0x7]];
            *kch = kchs[k];
            break;
                                                       // cases 8-15 are intended for golr selection modes
        case 8:                                        // 8,9 simple optimization for period
        case 9:                                        // 8 is min period 9 is max period
            if(repselect&0x1) {
                for(extreme1= 0,k=0;k<s;k++) extreme1 = ((p[k]>= extreme1)&&(d[k]<=4) ? p[k] : extreme1); // find value of fittest gene with max period
            }
            else  {
                for(extreme1=16,k=0;k<s;k++) extreme1 = ((p[k]<= extreme1)&&(d[k]<=4)) ? p[k] : extreme1; // find value of fittest gene with min period
            }
            for(bestnbmask=0ull,nbest=0,k=0;k<s;k++) bestnbmask |= ((p[k]==extreme1)? 1ull<<(k+0*nbest++) : 0ull); // find set of genes with equal best value
            *birth = ((nbest>0) ? 1ull: 0ull);         // birth condition may include later that genes not all same
            for(k=0;k<s;k++) if((bestnbmask>>k)&0x1) break;
            if (k==s) {k=0;*birth = 0ull;}             // in case no genes with best value, no birth, avoid k being out of bounds below
            *newgene = livegenes[k&0x7];               // choose first of selected set to replicate (can make positional dependent choice instead externally)
            *parentid=golb[ijanc[k&0x7]];
            *kch = kchs[k];
            break;
        case 10:                                       // optimization for period and then most accurate periodicity (least mismatches)
        case 11:                                       // 10 is min period 11 is max period : both with secondary selection for precision of periodicity
            if(repselect&0x1) {
                for(extreme1= 0,k=0;k<s;k++) extreme1 = ((p[k]>= extreme1)&&(d[k]<=4) ? p[k] : extreme1); // find value of fittest gene with max period
                for(extreme2=16,k=0;k<s;k++) extreme2 = ((p[k]==extreme1)&&(d[k]<extreme2)) ? d[k] : extreme2; // find genes with most precise periodicity at this period
            }
            else  {
                for(extreme1=16,k=0;k<s;k++) extreme1 = ((p[k]<= extreme1)&&(d[k]<=4)) ? p[k] : extreme1; // find value of fittest gene with min period
                for(extreme2=16,k=0;k<s;k++) extreme2 = ((p[k]==extreme1)&&(d[k]<extreme2)) ? d[k] : extreme2; // find genes with most precise periodicity at this period
            }
            for(bestnbmask=0ull,nbest=0,k=0;k<s;k++) bestnbmask |= (((p[k]==extreme1)&&(d[k]==extreme2))? 1ull<<(k+0*nbest++) : 0ull); // find set of genes with equal best value
            *birth = ((nbest>0) ? 1ull: 0ull);         // birth condition may include later that genes not all same
            for(k=0;k<s;k++) if((bestnbmask>>k)&0x1) break;
            if (k==s) {k=0;*birth = 0ull;}             // in case no genes with best value, no birth, avoid k being out of bounds below
            *newgene = livegenes[k&0x7];               // choose first of selected set to replicate (can make positional dependent choice instead externally)
            *parentid=golb[ijanc[k&0x7]];
            *kch = kchs[k];
            break;
        case 12:                                       // optimization for large displacement and then for small period
        case 13:                                       // optimization for large displacement and then large period
            for(extreme1=0,k=0;k<s;k++) {
                d0=abs(psx[k])+abs(psy[k]);
                extreme1 = ((d0 >= extreme1)&&(d[k]<=4) ? d0 : extreme1); // find value of fittest gene with max displacement for clear periodicity
            }
            if(repselect&0x1) {
                for(extreme2=16,k=0;k<s;k++) {
                    d0=abs(psx[k])+abs(psy[k]);
                    extreme2 = ((d0==extreme1)&&(p[k]<extreme2)&&(d0>=4)) ? p[k] : extreme2; // find genes with shortest period at max displacement
                }
            }
            else  {
                for(extreme2=0,k=0;k<s;k++) {
                    d0=abs(psx[k])+abs(psy[k]);
                    extreme2 = ((d0==extreme1)&&(p[k]>extreme2)&&(d0>=4)) ? p[k] : extreme2; // find genes with longest period at max displacement
                }
            }
            for(bestnbmask=0ull,nbest=0,k=0;k<s;k++) {
                d0=abs(psx[k])+abs(psy[k]);
                bestnbmask |= (((d0==extreme1)&&(p[k]==extreme2))? 1ull<<(k+0*nbest++) : 0ull); // find set of genes with equal best value
            }
            *birth = ((nbest>0) ? 1ull: 0ull);         // birth condition may include later that genes not all same
            for(k=0;k<s;k++) if((bestnbmask>>k)&0x1) break;
            if (k==s) {k=0;*birth = 0ull;}             // in case no genes with best value, no birth, avoid k being out of bounds below
            *newgene = livegenes[k&0x7];               // choose first of selected set to replicate (can make positional dependent choice instead externally)
            *parentid=golb[ijanc[k&0x7]];
            *kch = kchs[k];
            break;
        case 14:                                        // optimization for diagonal displacement and then for small period : currently same as 12
        case 15:                                        // optimization for diagonal displacement and then for large period : currently same as 13
            for(extreme1=0,k=0;k<s;k++) {
                d0=abs(psx[k])+abs(psy[k]);
                extreme1 = ((d0 >= extreme1)&&(d[k]<=4) ? d0 : extreme1); // find value of fittest gene with max displacement for clear periodicity
            }
            if(repselect&0x1) {
                for(extreme2=16,k=0;k<s;k++) {
                    d0=abs(psx[k])+abs(psy[k]);
                    extreme2 = ((d0==extreme1)&&(p[k]<extreme2)&&(d0>=4)) ? p[k] : extreme2; // find genes with shortest period at max displacement
                }
            }
            else  {
                for(extreme2=0,k=0;k<s;k++) {
                    d0=abs(psx[k])+abs(psy[k]);
                    extreme2 = ((d0==extreme1)&&(p[k]>extreme2)&&(d0>=4)) ? p[k] : extreme2; // find genes with longest period at max displacement
                }
            }
            for(bestnbmask=0ull,nbest=0,k=0;k<s;k++) {
                d0=abs(psx[k])+abs(psy[k]);
                bestnbmask |= (((d0==extreme1)&&(p[k]==extreme2))? 1ull<<(k+0*nbest++) : 0ull); // find set of genes with equal best value
            }
            *birth = ((nbest>0) ? 1ull: 0ull);         // birth condition may include later that genes not all same
            for(k=0;k<s;k++) if((bestnbmask>>k)&0x1) break;
            if (k==s) {k=0;*birth = 0ull;}             // in case no genes with best value, no birth, avoid k being out of bounds below
            *newgene = livegenes[k&0x7];               // choose first of selected set to replicate (can make positional dependent choice instead externally)
            *parentid=golb[ijanc[k&0x7]];
            *kch = kchs[k];
            break;
        default:
            fprintf(stderr,"Error: s = %d live gene repselect %d is not implemented\n",s,repselect);
            exit(1);
    }
    for (*nbmask=0ull,k=0;k<s;k++)
        *nbmask |= ((bestnbmask>>k)&0x1ull)<<((nb1i>>(k<<2))&0x7);

    return(nbest);
}
//------------------------------------------------------------- selectone_nbs ---------------------------------------------------------------------------
extern INLINE void selectone_nbs(int s, uint64_t nb2i, int nb[], uint64_t gol[], uint64_t golg[],uint64_t golb[], uint64_t * birth, uint64_t *newgene, uint64_t *parentid, unsigned int *kch) {
// birth is returned 1 if ancestors satisfy selection condition. Selection of which of two genes to copy is newgene.
    int k, kanc, l, nb1[8], ij1, i1, j1, j1p1, j1m1, i1p1, i1m1;
    uint64_t nbi, nbil, gene, genelink, s2, sl;
    unsigned int cmask;

    for (k=kanc=0,s2=0;k<s;k++) {                           // loop only over live neigbours, s2 is number of live nbs in connected 2nd ring
        nbi = (nb2i>>(k<<2))&0x7;                           // kth live neighbour index
        ij1 = nb[nbi];                                      // neighbour site ij index
        gene = golg[ij1];                                   // gene at neighbour site
        cmask = 0;                                          // connection mask initialized to 0
        sl = 0;                                             // initialize number of connected live neighbours of this neighbour
        for(l=0;l<2;l++) {                                  // 2 possible connections encoded in 2 16-bit gene words
            genelink = (gene >> ((l+1)<<4)) & codingmask;   // ncoding bit sequences describing possible links: ncoding <=16
            if (genelink == codingmask) cmask = cmask|(1<<l);// set mask only if connection encoded (all ones), later use probs
        }
        if(cmask) {                                                  // only if there are some connections
            i1 = ij1 & Nmask;  j1 = ij1 >> log2N;                                   // row & column
            j1p1 = ((j1+1) & Nmask)*N; j1m1 = ((j1-1) & Nmask)*N;                   // toroidal (j+1)*N and (j-1)*N
            i1p1 =  (i1+1) & Nmask; i1m1 =  (i1-1) & Nmask;                         // toroidal i+1, i-1
            nb1[0]=j1m1+i1m1; nb1[1]=j1m1+i1; nb1[2]=j1m1+i1p1; nb1[3]=j1*N+i1p1;   //next nbs  0 to 3
            nb1[4]=j1p1+i1p1; nb1[5]=j1p1+i1; nb1[6]=j1p1+i1m1; nb1[7]=j1*N+i1m1;   //next nbs  4 to 7
            for(l=0;l<2;l++) {
                if ((cmask>>l)&0x1) {
                    nbil = (nbi+l)&0x7; // on the 2nd ring, wrt nb in direction k,the 2 nbs at & after (l=0,1)
                    if(gol[nb1[nbil]]) {
                        s2++;sl++;
                    } // if
                }  // if
            }  // for
        } // if
        if (sl&0x2) kanc = k;   // neighbour contributing 2 live next shell neighbours serves as ancestor, only works for repscheme==2
    } // for
    if(s2==3) {                 // 3 live neighbours in 2nd shell pointed to by live first shell neighbours
        *birth = 1ull;
        k = (nb2i>>(kanc<<2))&0x7;
        ij1 = nb[k];
        *parentid=golb[ij1];
        *newgene = golg[ij1];
        *kch = k;
    }
    // else no change to birth or kch
}
//-------------------------------------------------------------- selectdifftx ---------------------------------------------------------------------------
extern INLINE unsigned int selectdifft0(uint64_t nbmask, int *crot, int *kodd) {
    *kodd = 0;
    *crot = 0;
    return(0);                                                 // replication of live nb in bit 0 of canonical rotation
}
//.......................................................................................................................................................
extern INLINE unsigned int selectdifft1(uint64_t nbmask, int *crot, int *kodd) {
// selection based on canonical rotation
    int k,kmin;
    uint64_t nbmaskr, nbmaskrm;

    for (k=1,nbmaskrm=nbmaskr=nbmask,kmin=0;k<8;k++) {         // compute canonical rotation (minimum) of this mask
        nbmaskr = ((nbmaskr & 0x1ull)<<7) | (nbmaskr>>1);      // 8 bit rotate right
        if (nbmaskr < nbmaskrm) {                              // choose minimal value of mask rotation
            nbmaskrm = nbmaskr;                                // neighbor mask rotate min is current rotation
            kmin = k;                                          // no of times rotated to right
        }
    }
    *kodd = kmin & 0x1;
    *crot = 0;
    return(kmin);                                              // replication of live nb in bit 0 of canonical rotation
}
//.......................................................................................................................................................
extern INLINE unsigned int selectdifft2(uint64_t nbmask, int *crot, int *kodd) {
// selection based on canonical rotation to bunched pair, choose clockwise or anti-clockwise one (for R_2_canonical_nb)
    int k,kmin;
    uint64_t nbmaskr, nbmaskrm;

    for (k=1,nbmaskrm=nbmaskr=nbmask,kmin=0;k<8;k++) {         // compute canonical rotation (minimum) of this mask
        nbmaskr = ((nbmaskr & 0x1ull)<<7) | (nbmaskr>>1);      // 8 bit rotate right
        if (nbmaskr < nbmaskrm) {                              // choose minimal value of mask rotation
            nbmaskrm = nbmaskr;                                // neighbor mask rotate min is current rotation
            kmin = k;                                          // no of times rotated to right
        }
    }
    *kodd = kmin & 0x1;                                        // if not canonical, replication of live neighbour in other (non zero k) position
    switch (nbmaskrm) {                                        //              x03    x05    x09    x11
        case 0x03ull : k = 1; *crot = 0; break;                // 00000011    |01.|  <-
        case 0x05ull : k = 2; *crot = 1; break;                // 00000101    |...|  |0.2|  <-
        case 0x09ull : k = 3; *crot = 2; break;                // 00001001    |...|  |...|  |0..|   <-
        case 0x11ull : k = 4; *crot = 3; break;                // 00010001           |...|  |..3|  |0..|   <-
        default  : {                                           //                           |...|  |...|
                                                               //                                  |..4|
            fprintf(stderr,"Error in canonical rotation for two live neighbours nbmaskrm = %"PRIx64" for mask %"PRIx64"\n",nbmaskrm,nbmask); k = 0;
        } //default case
    } //switch

    if (canonical) return(kmin);                               // replication of live neigbour in bit 0 of canonical rotation
    else return((kmin+k)&0x7);                                 // rotate unique nb k left (kmin) back to orig nb pat
}
//.......................................................................................................................................................
extern INLINE unsigned int selectdifft3(uint64_t nbmask, int *crot, int *kodd) {
    unsigned int k,kmin;
    uint64_t nbmaskr,nbmaskrm;

    for (k=1,nbmaskrm=nbmaskr=nbmask,kmin=0;k<8;k++) {         // compute canonical rotation (minimum) of this mask
        nbmaskr = ((nbmaskr & 1ull)<<7) | (nbmaskr>>1);        // 8 bit rotate right
        if (nbmaskr < nbmaskrm) {                              // choose minimal value of mask rotation
            nbmaskrm = nbmaskr;                                // neighbor mask rotate min is current rotation
            kmin = k;                                          // no of times rotated to right
        }
    }
    *kodd = kmin & 0x1;                                        // replication of live neighbour in most different position
    switch (nbmaskrm) {                                        //              x07    x0b    x0d    x13    x15    x19    x25
        case 0x07ull : k = 1; *crot = 0; break;                // 00000111    |012|  <-
        case 0x0bull : k = 0; *crot = 1; break;                // 00001011    |...|  |01.|  <-
        case 0x0dull : k = 3; *crot = 2; break;                // 00001101    |...|  |..3|  |0.2|   <-
        case 0x13ull : k = 1; *crot = 3; break;                // 00010011           |...|  |..3|  |01.|   <-
        case 0x15ull : k = 2; *crot = 4; break;                // 00010101                  |...|  |...|  |0.2|   <-
        case 0x19ull : k = 0; *crot = 5; break;                // 00011001                         |..4|  |...|  |0..|   <-
        case 0x25ull : k = 5; *crot = 6; break;                // 00100101                                |..4|  |..3|  |0.2|  <-
        default  : {                                           //                                                |..4|  |...|
                                                               //                                                       |.5.|
            fprintf(stderr,"Error in canonical rotation for three live neighbours nbmaskrm = %"PRIx64" for mask %"PRIx64"\n",nbmaskrm,nbmask); k = 0;
        } //default case
    } //switch
    if (canonical) return(kmin);                               // replication of live neigbour in bit 0 of canonical rotation
    else return((kmin+k)&0x7);                                 // rotate unique nb k left (kmin) back to orig nb pat

}
//.......................................................................................................................................................
extern INLINE unsigned int selectdifft4(uint64_t nbmask, int *crot, int *kodd) {
    int k,kmin;
    uint64_t nbmaskr,nbmaskrm;

    for (k=1,nbmaskrm=nbmaskr=nbmask,kmin=0;k<8;k++) {         // compute canonical rotation (minimum) of this mask
        nbmaskr = ((nbmaskr & 0x1ull)<<7) | (nbmaskr>>1);      // 8 bit rotate right
        if (nbmaskr < nbmaskrm) {                              // choose minimal value of mask rotation
            nbmaskrm = nbmaskr;                                // neighbor mask rotate min is current rotation
            kmin = k;                                          // no of times rotated to right
        }
    }
    *kodd = kmin&0x1;                                          // replication of live neighbour in most central position (left disambiguation)
    switch (nbmaskrm) {                                        //              x07f    x17    x1b    x1d    x27    x2b    x2d   x33   x35    x55
        case 0x0full : k = 1; *crot = 0; break;                // 00001111    |012|  <-
        case 0x17ull : k = 2; *crot = 1; break;                // 00010111    |..3|  |012|  <-
        case 0x1bull : k = 1; *crot = 2; break;                // 00011011    |...|  |...|  |01.|   <-
        case 0x1dull : k = 2; *crot = 3; break;                // 00011101           |..4|  |..3|  |0.2|   <-
        case 0x27ull : k = 2; *crot = 4; break;                // 00100111                  |..4|  |..3|  |012|   <-
        case 0x2bull : k = 3; *crot = 5; break;                // 00101011                         |..4|  |...|  |01.|   <-
        case 0x2dull : k = 2; *crot = 6; break;                // 00101101                                |.5.|  |..3|  |0.2|  <-
        case 0x33ull : k = 1; *crot = 7; break;                // 00110011                                       |.5.|  |..3|  |01.|  <-
        case 0x35ull : k = 2; *crot = 8; break;                // 00110101                                              |.5.|  |...|  |0.2|  <-
        case 0x55ull : k = 2; *crot = 9; break;                // 01010101                                                     |.54|  |...|  |0.2|  <-
        default  : {                                           //                                                                     |.54|  |...|
                                                               //                                                                            |6.4|
            fprintf(stderr,"Error in canonical rotation for four live neighbours nbmaskrm = %"PRIx64" for mask %"PRIx64"\n",nbmaskrm,nbmask); k = 0;
        } //default case
    } //switch
    if (canonical) return(kmin);                               // replication of live neigbour in bit 0 of canonical rotation
    else return((kmin+k)&0x7);                                 // rotate unique nb k left (kmin) back to orig nb pat
}
//.......................................................................................................................................................
extern INLINE unsigned int selectdifft5(uint64_t nbmask, int *crot, int *kodd) {
    unsigned int k,kmin;
    uint64_t nbmaskr,nbmaskrm;

    for (k=1,nbmaskrm=nbmaskr=nbmask,kmin=0;k<8;k++) {         // compute canonical rotation (minimum) of this mask
        nbmaskr = ((nbmaskr & 1ull)<<7) | (nbmaskr>>1);        // 8 bit rotate right
        if (nbmaskr < nbmaskrm) {                              // choose minimal value of mask rotation
            nbmaskrm = nbmaskr;                                // neighbor mask rotate min is current rotation
            kmin = k;                                          // no of times rotated to right
        }
    }
    *kodd = kmin & 0x1;
    switch (nbmaskrm) {                                        //              x1f    x2f    x3d    x3b    x57    x37    x5b
        case 0x1full : k = 2; *crot = 0; break;                // 00011111    |012|  <-
        case 0x2full : k = 2; *crot = 1; break;                // 00101111    |..3|  |012|  <-
        case 0x3dull : k = 3; *crot = 2; break;                // 00111101    |..4|  |..3|  |0.2|   <-
        case 0x3bull : k = 3; *crot = 3; break;                // 00111011           |.5.|  |..3|  |01.|   <-
        case 0x57ull : k = 2; *crot = 4; break;                // 01010111                  |.54|  |..3|  |012|   <-
        case 0x37ull : k = 2; *crot = 5; break;                // 00110111                         |.54|  |...|  |012|   <-
        case 0x5bull : k = 3; *crot = 6; break;                // 01011011                                |6.4|  |...|  |01.|  <-
        default  : {                                           //                                                |.54|  |..3|
                                                               //                                                       |6.4|
            fprintf(stderr,"Error in canonical rotation for five live neighbours nbmaskrm = %"PRIx64" for mask %"PRIx64"\n",nbmaskrm,nbmask); k = 0;
        } //default case
    } //switch
    if (canonical) return(kmin);                               // replication of live neigbour in bit 0 of canonical rotation
    else return((kmin+k)&0x7);                                 // rotate unique nb k left (kmin) back to orig nb pat
}
//.......................................................................................................................................................
extern INLINE unsigned int selectdifft6(uint64_t nbmask, int *crot, int *kodd) {
// selection based on canonical rotation to bunched pair, choose clockwise or anti-clockwise one (for R_2_canonical_nb)
    int k,kmin;
    uint64_t nbmaskr, nbmaskrm;

    for (k=1,nbmaskrm=nbmaskr=nbmask,kmin=0;k<8;k++) {         // compute canonical rotation (minimum) of this mask
        nbmaskr = ((nbmaskr & 0x1ull)<<7) | (nbmaskr>>1);      // 8 bit rotate right
        if (nbmaskr < nbmaskrm) {                              // choose minimal value of mask rotation
            nbmaskrm = nbmaskr;                                // neighbor mask rotate min is current rotation
            kmin = k;                                          // no of times rotated to right
        }
    }
    *kodd = kmin & 0x1;
    switch (nbmaskrm) {                                        //              x3f    x5f    x6f    x77
        case 0x3full : k = 2; *crot = 0; break;                // 00111111    |012|  <-
        case 0x5full : k = 3; *crot = 1; break;                // 01011111    |..3|  |012|  <-
        case 0x6full : k = 3; *crot = 2; break;                // 01101111    |.54|  |..3|  |012|   <-
        case 0x77ull : k = 4; *crot = 3; break;                // 01110111           |6.4|  |..3|  |012|   <-
        default  : {                                           //                           |65.|  |...|
                                                               //                                  |654|
            fprintf(stderr,"Error in canonical rotation for six live neighbours nbmaskrm = %"PRIx64" for mask %"PRIx64"\n",nbmaskrm,nbmask); k = 0;
        } //default case
    } //switch
    if (canonical) return(kmin);                               // replication of live neigbour in bit 0 of canonical rotation
    else return((kmin+k)&0x7);                                 // rotate unique nb k left (kmin) back to orig nb pat
}
//.......................................................................................................................................................
extern INLINE unsigned int selectdifft7(uint64_t nbmask, int *crot, int *kodd) {
// selection based on canonical rotation
    int k,kmin;
    uint64_t nbmaskr, nbmaskrm;

    for (k=1,nbmaskrm=nbmaskr=nbmask,kmin=0;k<8;k++) {         // compute canonical rotation (minimum) of this mask
        nbmaskr = ((nbmaskr & 0x1ull)<<7) | (nbmaskr>>1);      // 8 bit rotate right
        if (nbmaskr < nbmaskrm) {                              // choose minimal value of mask rotation
            nbmaskrm = nbmaskr;                                // neighbor mask rotate min is current rotation
            kmin = k;                                          // no of times rotated to right
        }
    }
    *kodd = kmin & 0x1;
    *crot = 0;
    if (canonical) return(kmin);                               // replication of live nb in bit 0 of canonical rotation
    else  return((kmin+3)&0x7);                                // replication of live nb in bit 3 (middle) of canonical rotation
}
//.......................................................................................................................................................
extern INLINE unsigned int selectdifft(int sum, uint64_t nbmask, int *crot, int *kodd, int *nsame) {
        int kch;
        *nsame = 0;

        switch(sum) {
                    case 0:  return(selectdifft0(nbmask, crot, kodd));
                    case 1:  return(selectdifft1(nbmask, crot, kodd));
                    case 2:  kch=selectdifft2(nbmask, crot, kodd);
                             if (*crot==3) *nsame = 2;
                             return(kch);
                    case 3:  return(selectdifft3(nbmask, crot, kodd));
                    case 4:  kch=selectdifft4(nbmask, crot, kodd);
                             // if (*crot==2) *nsame = 2;      // this case can be resolved with good rotation symmetry
                             if (*crot==7) *nsame = 2;
                             else if (*crot==9) *nsame = 4;
                             return(kch);
                    case 5:  return(selectdifft5(nbmask, crot, kodd));
                    case 6:  kch=selectdifft6(nbmask, crot, kodd);
                             if (*crot==3) *nsame = 2;
                             return(kch);
                    case 7:  return(selectdifft7(nbmask, crot, kodd));
                    default: return(0);
        }
}
//.......................................................................................................................................................
extern INLINE void analyze_nbs(int ij, uint64_t gol[], int nnb[], uint64_t *nnb1i, int *ns) {
    int i,j,ip1,im1,jp1,jm1,s,k;
    uint64_t nb1i,gols;
    
    i = ij & Nmask;  j = ij >> log2N;                                       // row & column
    jp1 = ((j+1) & Nmask)*N; jm1 = ((j-1) & Nmask)*N;                       // toroidal (j+1)*N and (j-1)*N
    ip1 =  (i+1) & Nmask; im1 =  (i-1) & Nmask;                             // toroidal i+1, i-1
    nnb[0]=jm1+im1; nnb[1]=jm1+i; nnb[2]=jm1+ip1; nnb[3]=j*N+ip1;           // new order of nbs
    nnb[4]=jp1+ip1; nnb[5]=jp1+i; nnb[6]=jp1+im1; nnb[7]=j*N+im1;
    for (s=0,nb1i=0ull,k=0;k<8;k++) {                                       // packs non-zero nb indices in first up to 8*4 bits
        gols=gol[nnb[k]];                                                   // whether neighbor is alive
        s += gols;                                                          // s is number of live nbs
        nb1i = (nb1i << (gols<<2)) + (gols*k);                              // nb1i is packed list of live neighbour indices
    }
    *nnb1i = nb1i;
    *ns = s;
}
//.......................................................................................................................................................
extern INLINE uint64_t disambiguate(unsigned int *kchx, uint64_t nb1i, int nb[],  uint64_t gol[], uint64_t golg[], uint64_t golb[], int nsame, uint64_t *birth, uint64_t *parentid, uint64_t *ancestor, int ij) {
    uint64_t gene,newgene,randnr;
    int k,ijanc,nnb[8],ns,discase,deathlikely;
    unsigned int kch,kch1;
    uint64_t nnb1i;
    
    kch = *kchx;
    discase=(repscheme>>8)&0x7;
    switch (discase) {
        case 0:  RAND128P(randnr);                                                      // random choice
                 kch += ((nsame-1)&(randnr>>32))<< (nsame ==4 ? 1 : 2);
                 kch &= 0x7;*kchx = kch;
                 ijanc = nb[kch];
                 *parentid = golb[ijanc];
                 *ancestor = golg[ijanc];
                 return( golg[ijanc]);
        case 1:  ijanc = nb[kch];                                                        // ignore asymmetry issue, continue regardless;
                 *parentid = golb[ijanc];
                 // if(*parentid == 0ull) fprintf(stderr,"error in disambiguate case 1: parentid set to golb[%d] which is 0 kch %d\n",ijanc,kch);
                 *ancestor = golg[ijanc];
                 return( golg[ijanc]);
        case 2:  *birth = 0ull; return(0ull);                                            // abandom birth attempt
        case 3:  *parentid = (((uint64_t) totsteps) <<32) + rootclone + ij;              // choose one GoL input gene, default ancestor for input genes
                 *ancestor = rootgene;
                 return(genegol[selection-8]);
        case 4:  for (newgene=golg[nb[kch]],kch1=kch,k=1;k<nsame;k++) {                   // choose minimum value gene
                     kch1+=k*(nsame==4 ? 2 : 4);kch1 &= 0x7;
                     ijanc=nb[kch1];
                     gene=golg[ijanc];
                     if (gene<newgene) {newgene = gene;kch=kch1;}
                 };
                 *kchx = kch;
                 ijanc=nb[kch];
                 *parentid = golb[ijanc];
                 *ancestor = newgene;
                 return(newgene);
        case 5:
        case 6:  for (newgene=golg[nb[kch]],kch1=kch,k=0;k<nsame;k++) {                        // choose first gene which is likely to die
                     kch1+=k*(nsame==4 ? 2 : 4);
                     kch1 &= 0x7;
                     ijanc=nb[kch1];
                     analyze_nbs(ijanc, gol, nnb, &nnb1i, &ns);
                     if(discase==5) deathlikely = ns<2 || ns>3;
                     else           deathlikely = ns<2 || ns>3 || ((overwritemask>>(ns-1))&0x1ull);
                     if(deathlikely) {                           // approx. substitute for death/overwrite calculation
                        newgene = golg[ijanc];
                        kch=kch1;
                        break;
                     }
                 };
                 *kchx = kch;
                 ijanc=nb[kch];
                 *parentid = golb[ijanc];
                 *ancestor = newgene;
                 return(newgene);
        case 7:  *parentid = (uint64_t) totsteps; *parentid = (*parentid <<32) + rootclone + ij;// default ancestor for input genes
                 RAND128P(randnr);                                                      // random choice
                 // *kchx = kch;    no change, retains kch unaltered as in case 1
                 *ancestor = rootgene;
                 return(randnr);                                                         // choose random gene
        default: fprintf(stderr,"Error in switch of ambiguous rotation resolution, should never reach here\n");
                 return(0ull);
    }
}
