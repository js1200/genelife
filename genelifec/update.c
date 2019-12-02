//
//  update.c
//  project genelife
//
//  update the genelife array for various model symmetries: s=2 or 3, lut_sum, lut_dist, lut_rot, lut_2D_sym
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
// update_23            update gol, golg, golgstats for a single synchronous time step : for selection 0-7 with fixed GoL rule departures in repscheme
// finish_update_ij     inline subroutine inside main array ij loop that is in common for all LUT based update rules
// finish_update        inline subroutine after main array loop that is in common for all LUT based update rules
// update_lut_sum       update version for gene encoding look up table for totalistic survival and birth (disallowing 0 live neighbour entries) sel 8,9
// update_lut_dist      update version for gene encoding look up table for survival and birth based on corner & edge sums (2*19 states disallowing s=0,1,7,8): sel 10,11
// update_lut_canon_rot update version for gene encoding look up table for canonical rotation survival and birth (2*32 states, disallowing 0,1,7,8 entries) : sel 12,13
// update_lut_2D_sym    update version all different configurations under the standard 2D 4-rotation and 4-reflection symmetries are distinguished: sel 14,15
// genelife_update      master routine to call specific model update, collect statistics if required and rotate planes
//---------------------------------------------------------------- update_23 ----------------------------------------------------------------------------
void update_23(uint64_t gol[], uint64_t golg[], uint64_t golgstats[], uint64_t golb[],uint64_t golr[],uint64_t newgol[], uint64_t newgolg[], uint64_t newgolgstats[], uint64_t newgolb[],uint64_t newgolr[]){
    // update for dissection of genetic rule variants within nearest neighbor sum s=2 or 3 only for survival and birth
    // update GoL for toroidal field which has side length which is a binary power of 2
    // encode without if structures for optimal vector treatment
    int s, s0, k, k1, kodd, nmut, crot;
    unsigned int kch,rulemodij,add2nd,mask1st;
    unsigned int select23live,pos_canon_neutral,survival,overwrite,enforcebirth,add2ndmask1st,nongolnottwice;
    int nb[8], nbc, nbch, ij, i, j, jp1, jm1, ip1, im1, ij1;
    uint64_t g, gs, nb1i, nb2i, randnr, r2;
    uint64_t nbmask, nbmaskr;
    uint64_t newgene, ancestor, livegenes[3], parentid;
    uint64_t s2or3, birth, statflag, nextgolstate;
    // short unsigned int patt; // needed if quadimage uncommented

    canonical = repscheme & R_2_canonical_nb;
    survival = survivalmask;
    overwrite = overwritemask;
    select23live = ((repscheme & R_0_2sel_3live)?1:0)+((repscheme & R_1_2sel_2live)?2:0);
    pos_canon_neutral = ((repscheme & R_2_canonical_nb)?1:0)+((repscheme & R_3_neutral_pos)?2:0);
    enforcebirth =((repscheme & R_4_enforce3birth)?1:0)+((repscheme & R_5_enforce2birth)?2:0);
    add2ndmask1st = ((repscheme & R_6_2ndnb_genes)?1:0)+((repscheme & R_7_1stnb_masks)?2:0);
    nongolnottwice = ((repscheme & R_8_nongolstat)?1:0)+((repscheme & R_9_nongolstatnbs)?2:0);
    add2nd = add2ndmask1st&0x1;
    parentdies = (repscheme & R_14_parentdies_23) ? 1 : 0;

    if(parentdies) for (ij=0; ij<N2; ij++) newgolgstats[ij] = 0ull;         // need to update statistics of neighbours with parenting information, so init required

    for (ij=0; ij<N2; ij++) {                                               // loop over all sites of 2D torus with side length N
        i = ij & Nmask;  j = ij >> log2N;                                   // row & column
        jp1 = ((j+1) & Nmask)*N; jm1 = ((j-1) & Nmask)*N;                   // toroidal (j+1)*N and (j-1)*N
        ip1 =  (i+1) & Nmask; im1 =  (i-1) & Nmask;                         // toroidal i+1, i-1
        nb[0]=jm1+im1; nb[1]=jm1+i; nb[2]=jm1+ip1; nb[3]=j*N+ip1;           // new order of nbs
        nb[4]=jp1+ip1; nb[5]=jp1+i; nb[6]=jp1+im1; nb[7]=j*N+im1;
        for (k=0,s=0,nb1i=0;k<8;k++) {                                      // packs non-zero nb indices in first up to 8*4 bits
            gs=gol[nb[k]];                                                  // whether neighbor is alive
            s += gs;                                                        // s is number of live nbs
            nb1i = (nb1i << (gs<<2)) + (gs*k);                              // nb1i is packed list of live neighbour indices
        }
        s0 = s;                                                             // record unmodified value of s for comparison
        s2or3 = (s>>2) ? 0ull : (s>>1);                                     // s == 2 or s ==3 : checked by bits 2+ are zero and bit 1 is 1
        nextgolstate = s2or3 ? (gol[ij] ? 1ull : (s&0x1ull ? 1ull : 0ull )) : 0ull;   // GoL standard calculation next state

        statflag = 0ull;
        rulemodij = (rulemod&0x2) ? (ij>=(N2>>1) ? 1 : 0) : (rulemod&0x1);   // if rulemod bit 1 is on then split into half planes with/without mod
        nbmask = 0;
        if(s>1) {
          if (repscheme & R_19_quadrant_2nb1)  add2ndmask1st =     (ij > (N2>>1) ? 0x2 : 0x0) + ((ij&Nmask)>(N>>1) ? 0x1 : 0x0);
          if (repscheme & R_20_quadrant_ngol)  nongolnottwice =    (ij > (N2>>1) ? 0x2 : 0x0) + ((ij&Nmask)>(N>>1) ? 0x1 : 0x0);
          add2nd = add2ndmask1st&0x1; mask1st=(add2ndmask1st>>1)&0x1;

          if(nongolnottwice&0x1) {                                          // check for non GoL changed states
            if((nongolnottwice>>1)&0x1) {                                   // check in neighborhood
                int sng;                                                    // sum of nbs with non GoL rule change bit set
                for (k=0,sng=0;k<8;k++) {                                   // calc number of neighbours resulting from nongol rule change
                    sng += (golgstats[nb[k]]&F_nongolchg)?1:0;              // neighbors with state set by a non GoL change
                }
                if(sng) rulemodij = 0ull;
            }
            if(golgstats[ij]&F_nongolchg) rulemodij = 0ull;                 // if central state the result of a non GoL rule change
          }
          if(mask1st&&rulemodij) {                                         // recalculate effective new s as less than original sum sm
            for (k=0,nbmaskr=0;k<8;k++) {                                   // depending on rotated overlay of masks in live neighbour genes
                if(gol[nb[k]]) {
                    g= golg[nb[k]];                                         // fnal gene has ncoding 0s then 8 bit mask
                    //if(!((g>>8) & codingmask)) nbmaskr |= g&0xff;           // tried |=, &=, ^= .
                    if(!((g>>8) & codingmask)) nbmaskr |= (0x2<<(g&0x7L))&0xff;  // only one bit on for gene masks in this version: not self, so 0x2L

                }
                nbmaskr = ((nbmaskr & 0x1ull)<<7) | (nbmaskr>>1);           // 8 bit rotate right
            }
            for (k=0,s=0,nb1i=0ull;k<8;k++) {                               // recalculate sum and nb1i using combined mask
                gs =gol[nb[k]] & (((~nbmaskr)>>k)&0x1ull);                  // if mask bit set, count as if dead
                nbmask |= gs<<k;                                            // also calculate nbmask for use below
                s += gs;
                nb1i = (nb1i << (gs<<2)) + (gs*k);
            }
            if(s!=s0) s2or3 = (s>>2) ? 0ull : (s>>1);                       // redo s2or3 calculation for modified s
          }
          else for (k=0,nbmask=0;k<8;k++) nbmask |= (gol[nb[k]]<<k);        // 8-bit mask of GoL states of 8 nbs, clockwise from top left
          statflag |= F_livenbs & (nbmask<<16);                             // record live neighbour pattern
        } // end if s>1
        if (s2or3) {                                                        // if 2 or 3 neighbours alive
            if (repscheme & R_quadrant) {                                   // quarter the plane with 4 different parameter values
                if (repscheme & R_16_quadrant_sele)  select23live =      (ij > (N2>>1) ? 0x2 : 0x0) + ((ij&Nmask)>(N>>1) ? 0x1 : 0x0);
                if (repscheme & R_17_quadrant_posn)  pos_canon_neutral = (ij > (N2>>1) ? 0x2 : 0x0) + ((ij&Nmask)>(N>>1) ? 0x1 : 0x0);
                if (repscheme & R_18_quadrant_enfb)  enforcebirth =      (ij > (N2>>1) ? 0x2 : 0x0) + ((ij&Nmask)>(N>>1) ? 0x1 : 0x0);
                if (repscheme & R_21_quadrant_surv)  survival =          (ij > (N2>>1) ? 0x2 : 0x0) + ((ij&Nmask)>(N>>1) ? 0x1 : 0x0);
                if (repscheme & R_22_quadrant_over)  overwrite =         (ij > (N2>>1) ? 0x2 : 0x0) + ((ij&Nmask)>(N>>1) ? 0x1 : 0x0);
                canonical = pos_canon_neutral&0x1;       // global value since needed in ...difft2-6 subroutines
            }
            birth = 0ull;
            newgene = 0ull;
            parentid = 0ull;
            kch = 0;

            if (s&0x1ull) {  // s==3                                        // allow birth (with possible overwrite)
              statflag |= F_3_live;                                         // record instance of 3 live nbs
              if ((0x1ull&overwrite)||!gol[ij] ) {                          // central site empty or overwrite mode
                birth = 1ull;                                               // birth flag
                for(k=0;k<s;k++) livegenes[k] = golg[nb[(nb1i>>(k<<2))&0x7]]; // live neighbour genes
                kch=selectdifft3(nbmask, &crot, &kodd);
                parentid=golb[nb[kch&0x7]];                                 // default parent is for kch unless selected below
                if((livegenes[0]^livegenes[1])|(livegenes[0]^livegenes[2])) { // genes not all same, need ancestor calculation
                  nbch=nb[kch];
                  if (select23live&0x1) {                                   // execute selective replication of one of two otherwise unchosen live genes
                      nb2i = 0ull;
                      for(k1=k=0;k<3;k++) {                                 // choice of two other live genes for possible ancestor
                          nbc=(nb1i>>(k<<2))&0x7;
                          if(nb[nbc]!=nbch) nb2i = (nbc<<(k1++<<2))+nb2i;
                      }
                      if (add2nd) selectone_nbs(s,nb2i,nb,gol,golg,golb,&birth,&newgene,&parentid,&kch);  //2nd nb modulation
                      else selectone_of_2(s,nb2i,nb,golg,golb,&birth,&newgene,&parentid,&kch);
                      if (birth==0ull) {                                    // optional reset of ancestor & birth if no ancestors chosen in selectone
                        if((enforcebirth&0x1)||rulemodij)  {                // birth cannot fail or genes don't matter or no modification to gol rules
                            newgene = golg[nbch];
                            parentid=golb[nbch];
                            birth = 1ull;
                        }
                      }
                      else statflag |= F_2select;                           // ancestor has been chosen in selectone_of_2
                  }
                  else {
                      newgene = golg[nbch];
                      parentid = golb[nbch];
                  }
                } // end if not all live neighbors the same
                else {
                    statflag |= F_3g_same;
                    newgene = livegenes[0];                                 // genes all the same : copy first one, parentid default
                    if((~enforcebirth&0x1) && rulemodij) birth = 0ull;      // no birth for 3 identical genes if not enforcebirth3 and rulemod
                    kch=selectdifft3(nbmask, &crot, &kodd);                 // need a single deterministic ancestor for genealogy and golr
                }
              } // end central site empty or overwrite mode
            }  // end if s==3
            else {  // s==2                                                 // possible birth as exception to GoL rule
                statflag |= F_2_live;
                if (((select23live>>1)&0x1)&&(rulemodij||gol[ij])) {        // rule departure from GOL allowed or possible overwrite
                    if ((0x1ull&(overwrite>>1))||!gol[ij]) {                // either overwrite on for s==2 or central site is empty
                        if (add2nd) {
                            selectone_nbs(s,nb1i,nb,gol,golg,golb,&birth,&newgene,&parentid,&kch); //2nd nb modulation
                        }
                        else {
                            nbmask = (0x1ull<<(nb1i&0x7)) + (0x1ull<<((nb1i>>4)&0x7));
                            kch=selectdifft2(nbmask, &crot, &kodd);
                            if ((pos_canon_neutral>>1)&0x1) {               // enforce gene independent birth for s = 2 (corrected 31.1.2019, remove repscheme&)
                                newgene = golg[nb[kch]];
                                parentid= golb[nb[kch]];
                                birth = 1ull;
                            }
                            else {
                                selectone_of_2(s,nb1i,nb,golg,golb,&birth,&newgene,&parentid,&kch);
                            }
                            if(repscheme & R_10_13_2birth_k4) {             // birth only on the active subset of 4 canonical 2-live nb configs if any are active
                                if(~repscheme & (R_10_2birth_k0<<(kch-1))) birth = 0ull;   // cancel birth if chosen configuration not active
                            }
                        }

                        if(!birth && (enforcebirth&0x2)) {
                            nbmask = (0x1ull<<(nb1i&0x7)) + (0x1ull<<((nb1i>>4)&0x7));
                            kch=selectdifft2(nbmask, &crot, &kodd);
                            newgene = golg[nb[kch]];
                            parentid= golb[nb[kch]];
                            if(repscheme & R_10_13_2birth_k4) {
                                if(repscheme & (R_10_2birth_k0<<(kch-1))) birth = 1ull;
                                else birth = 0ull;
                            }
                            else birth = 1ull;
                        }
                        if (birth) statflag |= F_2select;
                    }
                }
            }

            if(birth){
                // compute random events for multiple single bit mutations, as well as mutation position nmut
                r2=1ull;
                ancestor = newgene;
                while (r2) {
                    RAND128P(randnr);                                           // inline exp so compiler recognizes auto-vec,
                    r2 = randprob(pmutmask,(unsigned int) randnr);
                    nmut = (randnr >> 56) & 0x3f;                               // choose mutation position for length 64 gene : from bits 56:61 of randnr
                    newgene = newgene ^ (r2<<nmut);                             // introduce single mutation with probability pmut = probmut
                    if(r2) {
                        statflag = statflag | F_mutation;
                        statflag = statflag | F_survmut;
                    }
                }
                if(diagnostics & diag_hash_genes) {
                    if(gol[ij]) {                                               // central old gene present: overwritten
                        hashdeletegene(golg[ij],golb[ij],"step %d hash delete error 1 in update, gene %"PRIx64" not stored\n");
                    }
                    hashaddgene(ij,newgene,ancestor,newgolb+ij,parentid,statflag & F_mutation);
                }
                newgol[ij]  =  1ull;                                        // new game of life cell value: alive
                newgolg[ij] =  newgene;                                     // if birth then newgene
                newgolr[ij] = (golr[nb[kch]]<<4) | kch | 0x8;               // register ancestor offset index in record_of_dynamics_gene golr along with 0x8 for birth event
                statflag = statflag | F_birth;
                if (parentdies) {
                    ij1 = nb[kch];
                    newgolgstats[ij1] = newgolgstats[ij1] | F_parent;
                }
            } // end birth
            else {
                if ((survival&s&0x1ull)|((survival>>1)&(~s)&0x1ull)|((~rulemodij)&0x1ull)) {// (surv bit 0 and s==3) or (surv bit 1 and s==2) or not rulemod1ij
                // if ((survival&s&0x1ull)|((survival>>1)&(~s)&0x1ull)) {   // survival bit 0 and s==3, or (survival bit 1 and s==2)
                    newgol[ij]  = gol[ij];                                  // new game of life cell value same as old
                    newgolg[ij] = golg[ij];                                 // gene stays as before, live or not
                    newgolb[ij] = golb[ij];
                    newgolr[ij] = (golr[ij]<<4) | (s-1);                    // register s-1 in record_of_dynamics_gene golr along with 0 value of bit3 for survival event
                    if(gol[ij]) {
                        statflag |= F_survival;
                        if (golgstats[ij]&F_survmut) statflag |= F_survmut;
                    }
                }
                else {
                    if(diagnostics & diag_hash_genes) {
                        if(gol[ij]) {                                       // death : need to update hash table
                            hashdeletegene(golg[ij],golb[ij],"step %d hash delete error 2 in update, gene %"PRIx64" not stored\n");
                        }
                    }
                    newgol[ij]  = 0ull;                                     // new game of life cell value dead
                    newgolg[ij] = 0ull;                                     // gene dies or stays dead
                    newgolb[ij] = 0ull;                                     // clone removed from site
                    newgolr[ij] = 0ull;                                     // displacement history erased on death
                    if(gol[ij]) statflag |= F_death;
                }
            } // end no birth
        }  // end if s2or3
        else {                                                              // else not birth or survival, 0 values for gol and gene
            if(diagnostics & diag_hash_genes) {
                if(gol[ij]) {                                               // death : need to update hash table
                    hashdeletegene(golg[ij],golb[ij],"step %d hash delete error 3 in update, gene %"PRIx64" not stored\n");
                }
            }
            newgol[ij]  = 0ull;                                             // new game of life cell value
            newgolg[ij] = 0ull;                                             // gene dies
            newgolb[ij] = 0ull;                                             // clone removed from site
            newgolr[ij] = 0ull;
            if(gol[ij]) statflag |= F_death;
        }
        if(gol[ij]) statflag |= F_golstate;
        if(newgol[ij]^nextgolstate) statflag |= F_notgolrul;
        if(gol[ij]^newgol[ij]) {
            statflag |= F_golchange;
            if(statflag&F_notgolrul) statflag |= F_nongolchg;
        }
        else if (golgstats[ij]&F_nongolchg) statflag |= F_nongolchg;        // maintain non-GoL chg status until state changed by GoL rule
        if (parentdies) newgolgstats[ij] |= statflag;                       // newgolgstats may already contain updated parenthood info F_parent
        else newgolgstats[ij] = statflag;
    }  // end for ij

    if(parentdies) {
        for (ij=0; ij<N2; ij++) {
            statflag = newgolgstats[ij];
            if(gol[ij] && !(statflag&F_birth) && !(statflag&F_death) && (statflag&F_parent)) {
                newgol[ij]  = 0ull;                                         // new game of life cell value dead
                newgolg[ij] = 0ull;                                         // gene dies
                newgolb[ij] = 0ull;
                newgolr[ij] = 0ull;
                if(diagnostics & diag_hash_genes)
                    hashdeletegene(golg[ij],golb[ij],"step %d hash delete error 2 in update, gene %"PRIx64" not stored\n");
                newgolgstats[ij] |= F_parentaldeath;
            }
        }
    }

    if(randominflux) random_influx(newgol,newgolg,newgolb,newgolr);
    if(vscrolling) v_scroll(newgol,newgolg,newgolb,newgolr);
    if ((colorfunction == 8) || (colorfunction2 == 8)) packandcompare(newgol,working,golmix);
    if(diagnostics & diag_component_labels) ncomponents=extract_components(newgol,newgolg);

    for (ij=0; ij<N2; ij++) {       // complete missing hash table records of extinction and activities
        if(gol[ij]) hashgeneextinction(golg[ij],"hash extinction storage error %d in update at step %d, gene %"PRIx64" not stored\n");
        if(newgol[ij]) hashgeneactivity(newgolg[ij],"hash activity storage error in update, gene %"PRIx64" not stored\n");
        if(newgolb[ij]) hashcloneactivity(newgolb[ij],"hash activity storage error in update, clone %"PRIx64" not stored\n");
    }
    // if(diagnostics & diag_hash_patterns) qimage = quadimage(newgol,&patt,log2N); // quadtree hash of entire image
}

//---------------------------------------------------------------- finish_update ------------------------------------------------------------------------
// two inline subroutines that are in common for all LUT based update rules
extern INLINE void finish_update_ij(int ij,int s,uint64_t golij,uint64_t gols,uint64_t nb1i,uint64_t nbmask,int nb[],uint64_t survive,uint64_t birth,
                                    uint64_t gol[],uint64_t golg[],uint64_t golgstats[],uint64_t golb[],uint64_t golr[],
                                    uint64_t newgol[],uint64_t newgolg[],uint64_t newgolgstats[],uint64_t newgolb[],uint64_t newgolr[]) {
        uint64_t randnr,r2,newgene,parentid,ancestor,statflag;
        int k, nmut, kodd, crot, nbest, nsame;
        unsigned int kch;
        statflag = F_livenbs & (nbmask<<16);                                    // record live neighbour pattern
        statflag |= F_s_live & (s&0x7);                                         // requires F_s_live to be in lowest 3 bits : val 8 mapped to 0
        if(birth) {                                                         // birth allowed by rules encoded in local genes (may still fail by selection)
            if (repscheme & R_11_random_resln) {
                RAND128P(randnr);                                           // inline exp so compiler recognizes auto-vec,
                kch = ((randnr>>32)&0xffff) % s;                            // choose random in this option only
                ancestor = newgene = golg[(nb[(nb1i>>(kch<<2))&0x7])];
                parentid = golb[(nb[(nb1i>>(kch<<2))&0x7])];
            }
            else {
                if ((ancselectmask>>(s-1)) &0x1) {                          // use genes to select ancestor
                    nbest=selectone_of_s(&kch,s,nb1i,nb,golg,golb,golr,&birth,&newgene,&parentid,&nbmask,ij);// selection scheme depends on repscheme parameter, selection depends on genes
                    ancestor = newgene;
                }
                else {                                                      // use positional information to select ancestor (leave birth on)
                    nbest = s;
                    ancestor = newgene = parentid = 0ull;kch = 0;           // ancestor, newgene, parentid, kch initialized here to avoid warning below (not needed though)
                    for (nbmask=0ull,k=0;k<s;k++) nbmask |= 0x1ull<<((nb1i>>(k<<2))&0x7);    // check whether this really needed here
                    if (!(nbest>1)) {                                       // for s==1 we define newgene ancestor immediately (s==0 does not reach here)
                        kch = nb1i&0x7;
                        ancestor = newgene = golg[nb[kch]];
                        parentid = golb[nb[kch]];
                    }
                }
                if(nbest>1 ) {
                    kch=selectdifft(nbest,nbmask,&crot,&kodd,&nsame);       // kch is chosen nb in range 0-7, nsame gives the number of undistinguished positions in canonical rotation
                    if(nsame) {
                        newgene = disambiguate(&kch, nb1i, nb, gol, golg, golb, nsame, &birth, &parentid, &ancestor, ij); // restore symmetry via one of 8 repscheme options
                        if(birth) statflag |= F_disambig;
                    }
                    else {
                        ancestor = newgene = golg[nb[kch]];
                        parentid = golb[nb[kch]];
                    }
                 }
            }
            if (birth) {                                                    // ask again because disambiguate may turn off birth
                statflag |= F_birth;
                r2=1ull;                                                    // compute random events for single bit mutation, as well as mutation position nmut
                while (r2) {
                    RAND128P(randnr);                                       // inline exp so compiler recognizes auto-vec,
                    r2 = randprob(pmutmask,(unsigned int) randnr);
                    nmut = (randnr >> 56) & 0x3f;                           // choose mutation position for length 64 gene : from bits 56:61 of randnr
                    newgene = newgene ^ (r2<<nmut);                         // introduce single mutation with probability pmut = probmut
                    if(r2) {
                        statflag = statflag | F_mutation;
                        statflag = statflag | F_survmut;
                    }                }
                newgol[ij]  =  1ull;                                        // new game of life cell value: alive
                newgolg[ij] =  newgene;                                     // if birth then newgene
                newgolr[ij] = (golr[nb[kch]]<<4) | kch | 0x8ull;            // register ancestor offset index in record_of_dynamics_gene golr along with 0x8 for birth event
                if(diagnostics & diag_hash_genes) {
                    if(golij) hashdeletegene(golg[ij],golb[ij],"step %d hash delete error 1 in update_lut_sum, gene %"PRIx64" not stored\n");
                    hashaddgene(ij,newgene,ancestor,newgolb+ij,parentid,statflag & F_mutation);
                }
                if (parentdies) {
                    unsigned int ij1 = nb[kch];
                    newgolgstats[ij1] = newgolgstats[ij1] | F_parent;
                }
            }
        }
        if(!birth) {                                                       // need instead of else because if(birth) section may change value of birth
            if(golij) {                                                    // death/survival
                if(survive) {                                              // survival coded
                    statflag |= F_survival;
                    if (golgstats[ij]&F_survmut) statflag |= F_survmut;    // gene is non-replicated mutant survivor
                    newgol[ij]  = golij;                                   // new game of life cell value same as old
                    newgolg[ij] = golg[ij];                                // gene stays same
                    newgolb[ij] = golb[ij];
                    newgolr[ij] = (golr[ij]<<4) | (s-1);                   // register s-1 in record_of_dynamics_gene golr along with 0 value of bit3 for survival event
                }
                else {                                                     // death
                    statflag |= F_death;
                    newgol[ij]  = 0ull;                                    // new game of life cell value dead
                    newgolg[ij] = 0ull;                                    // gene dies
                    newgolb[ij] = 0ull;
                    newgolr[ij] = 0ull;
                    if(diagnostics & diag_hash_genes)
                        hashdeletegene(golg[ij],golb[ij],"step %d hash delete error 2 in update, gene %"PRIx64" not stored\n");
                }
            }
            else {                                                         // empty and no birth, stays empty
                newgol[ij]  = golij;
                newgolg[ij] = golg[ij];
                newgolb[ij] = golb[ij];
                newgolr[ij] = golr[ij];
            }
        }
        if(newgol[ij]!=gols) {
            statflag |= F_notgolrul;
            if(newgol[ij]) statflag |= F_nongolchg;
        }
        if(golij) statflag |= F_golstate;                                   // this is the last gol state, not the new state
        if (parentdies) newgolgstats[ij] = newgolgstats[ij] | statflag;     // newgolgstats may already contain updated parenthood info F_parent
        else newgolgstats[ij] = statflag;
    
}
//........................................................................................................................................................
extern INLINE void finish_update(uint64_t newgol[], uint64_t newgolg[],uint64_t newgolgstats[],uint64_t newgolb[], uint64_t newgolr[], int nbshist[]) {
    int ij,k;
    uint64_t statflag;
    if(parentdies) {
        for (ij=0; ij<N2; ij++) {
            statflag = newgolgstats[ij];
            if(gol[ij] && !(statflag&F_birth) && !(statflag&F_death) && (statflag&F_parent)) {
                newgol[ij]  = 0ull;                                    // new game of life cell value dead
                newgolg[ij] = 0ull;                                    // gene dies
                newgolb[ij] = 0ull;
                newgolr[ij] = 0ull;
                if(diagnostics & diag_hash_genes)
                    hashdeletegene(golg[ij],golb[ij],"step %d hash delete error 2 in update, gene %"PRIx64" not stored\n");
                newgolgstats[ij] |= F_parentaldeath;
            }
        }
    }
    for (ambigsum=0,ij=0; ij<N2; ij++) ambigsum += (newgolgstats[ij]&F_disambig) ? 1 : 0;
    
    for(k=0;k<8;k++) nbshist[k]=0;
    for (ij=0; ij<N2; ij++) {
        if(newgolgstats[ij]&F_birth)
            for (k=0; k<8; k++) nbshist[k] += (newgolgstats[ij]&(0x1ull<<(16+k))) ? 1: 0;
    }
    
    if(randominflux) random_influx(newgol,newgolg,newgolb,newgolr);
    if(vscrolling) v_scroll(newgol,newgolg,newgolb,newgolr);
    if ((colorfunction == 8) || (colorfunction2 == 8)) packandcompare(newgol,working,golmix);
    if(diagnostics & diag_component_labels) ncomponents=extract_components(newgol,newgolg);
    if(diagnostics & diag_hash_genes) {
        for (ij=0; ij<N2; ij++) {       // complete missing hash table records of extinction and activities
            if(gol[ij]) hashgeneextinction(golg[ij],"hash extinction storage error %d in update at step %d, gene %"PRIx64" not stored\n");     // [**gol**]
            if(newgol[ij]) hashgeneactivity(newgolg[ij],"hash activity storage error in update, gene %"PRIx64" not stored\n");   // [**gol**]
        }
    }
    // if(diagnostics & diag_hash_patterns) qimage = quadimage(newgol,&patt,log2N); // quadtree hash of entire image
}
//---------------------------------------------------------------- update_lut_sum -----------------------------------------------------------------------
void update_lut_sum(uint64_t gol[], uint64_t golg[], uint64_t golgstats[], uint64_t golb[],uint64_t golr[],uint64_t newgol[], uint64_t newgolg[], uint64_t newgolgstats[], uint64_t newgolb[],uint64_t newgolr[]){    // selection models 8,9
// this version should work even if extra information is packed in the higher bits of gol: previous version relabelled to update_lut_sumx
// update GoL for toroidal field which has side length which is a binary power of 2
// encode without if structures for optimal vector treatment
/*
        0 <= s <= 8
        genome =
        8 bits for each possible s value for birth (center site 0) : exclude s=0 no birth of isolated live cells
        8 bits for each possible s value for survival/death (center site 1) : exclude s=0 no survival of isolated cells
        1. fixed length encoding for selection=8: using b0 for s=1, b1 for s=2, etc
        s-1   76543210    76543210
        GOL = 00000100(0)|00000110(0)
            = 0000 0100 0000 0110
            = 0x0406
        2. variable length encoding for selection=9: using 4 bit patterns anywhere on 4-bit raster in genome
        note that one entry suffices for any lut rule: duplicate rules may be used to encode robustness
        of four bits 0-3, bit 3 encodes birth/survival as 1/0. Bits 0-2 encode the value of s (0-7, s=8 excluded).
        GOL = 0xb32 or 0x23b2b for example (up to 64 bits)
                                                                                                                */
    int s, s1, s2or3, k, kch, nrfound_b, nrfound_s;
    unsigned int rulemodij;
    int nb[8],  ij, i, j, jp1, jm1, ip1, im1;
    uint64_t genecode, genecodeb, genecodes, gols, golij, nb1i, nbmask, found;
    uint64_t  survive, birth, overwrite, survivalgene, smask, bmask, ncodingmask, allcoding;

    canonical = repscheme & R_2_canonical_nb;                                   // set global choice of canonical rotation bit choice for selectdifftx
    survivalgene = repscheme & R_12_survivalgene;                                // gene determining survival is 1: central gene 0: determined by neighbours
    smask = (uint64_t) survivalmask;                                            // convert to 64 bit mask for efficient usage here
    bmask = (uint64_t) birthmask;
    parentdies = (repscheme & R_3_parentdies) ? 1 : 0;
    if(parentdies) for (ij=0; ij<N2; ij++) newgolgstats[ij] = 0ull;             // need to update statistics of neighbours with parenting information, so init required
    
    ncodingmask = (1ull<<ncoding)-1ull;                                         // mask for number of bits coding for each lut rule: <=4 for birth and survival case
    if (ncoding==4)         allcoding = 0xffffffffffffffff;                     // mask for total gene coding region
    else if (ncoding ==2)   allcoding = 0xffffffff;
    else                    allcoding = 0xffff;
  
    for (ij=0; ij<N2; ij++) {                                                   // loop over all sites of 2D torus with side length N
        i = ij & Nmask;  j = ij >> log2N;                                       // row & column
        jp1 = ((j+1) & Nmask)*N; jm1 = ((j-1) & Nmask)*N;                       // toroidal (j+1)*N and (j-1)*N
        ip1 =  (i+1) & Nmask; im1 =  (i-1) & Nmask;                             // toroidal i+1, i-1
        nb[0]=jm1+im1; nb[1]=jm1+i; nb[2]=jm1+ip1; nb[3]=j*N+ip1;               // new order of nbs
        nb[4]=jp1+ip1; nb[5]=jp1+i; nb[6]=jp1+im1; nb[7]=j*N+im1;
        for (s=0,nb1i=0ull,nbmask=0ull,k=0;k<8;k++) {                           // packs non-zero nb indices in first up to 8*4 bits
            gols=gol[nb[k]];                                                    // whether neighbor is alive
            s += gols;                                                          // s is number of live nbs
            nb1i = (nb1i << (gols<<2)) + (gols*k);                              // nb1i is packed list of live neighbour indices (each in range 0-7)
            nbmask |= (gols << k);
        }
        birth = 0ull;
        golij=gol[ij];
        if (s) {
            s1=s-1;
            s2or3 = (s>>2) ? 0ull : (s>>1);                                     // s == 2 or s ==3 : checked by bits 2+ are zero and bit 1 is 1
            gols = s2or3 ? (golij ? 1ull : (s&1ull ? 1ull : 0ull )) : 0ull;     // GoL calculation next state for non-genetic gol plane
            rulemodij = (rulemod&0x4) ? MEMBRANE : ((rulemod&0x2) ? (ij>=(N2>>1) ? 1 : 0) : (rulemod&0x1)); // if rulemod bit 2 then activate membrane of death
                                                                                // else if rulemod bit 1 is on then split into half planes with/without mod
            if(rulemodij==1) {
                overwrite = overwritemask&(0x1ull<<s1);
                if (selection==9) {                                             // selection == 9 : NB selection 9 does not use ncoding to determine coding length
                    nrfound_b = nrfound_s = 0;
                    for (k=0;k<s;k++) {                                         // decodes genes using variable position encoding: fixed length genecode from nbs variable position encoding
                        kch = (nb1i>>(k<<2))&0x7;
                        if (!survivalgene && golij) {
                            PATTERN4(golg[nb[kch]], (s&0x7), found);            //survival?
                            if(found) nrfound_s ++;
                        }
                        if (overwrite || !golij) {
                            PATTERN4(golg[nb[kch]], (s|0x8), found);            //birth?
                            if(found) nrfound_b ++;
                        }
                    }
                    if (survivalgene && golij) {                                // survival determined by central gene in this case
                            PATTERN4(golg[ij], (s&0x7), found);                 // survival?
                            if(found) nrfound_s = s;
                    }

                    if(repscheme&R_0_nb_majority) {
                        if((repscheme&R_1_nb_OR_AND)||(s&0x1)) found = (nrfound_s >  (s>>1)) ? 1 : 0;  // > s/2
                        else                                   found = (nrfound_s >= (s>>1)) ? 1 : 0;  // >= s/2
                    }
                    else {
                        if(repscheme&R_1_nb_OR_AND) found = nrfound_s ? 1 : 0;          // OR
                        else                        found = (nrfound_s == s) ? 1 : 0;   // AND
                    }
                    genecodes = found? (1ull << s1) : 0ull;
                    
                    if(repscheme&R_0_nb_majority) {
                        if((repscheme&R_1_nb_OR_AND)||(s&0x1)) found = (nrfound_b >  (s>>1)) ? 1 : 0;   // > s/2
                        else                                   found = (nrfound_b >= (s>>1)) ? 1 : 0;   // >= s/2
                    }
                    else {
                        if(repscheme&R_1_nb_OR_AND) found = nrfound_b ? 1 : 0;
                        else                        found = (nrfound_b == s) ? 1 : 0;
                    }
                    genecodeb = found ? (1ull << (s1+8)) : 0ull;
                    
                    genecode=genecodes|genecodeb;
                    survive = ((genecode&smask)>>s1) & 0x1ull;
                    genecode>>=8;
                    if (overwrite || !golij) birth   = ((genecode&bmask)>>s1) & 0x1ull;
                }
                else {                                                          // selection == 8
                    if(repscheme&R_0_nb_majority) {                                 // MAJORITY resolution of neighbours
                        for (nrfound_b=nrfound_s=0,k=0;k<s;k++) {                   // decodes genes with fixed length encoding length ncoding 1,2,4 bits
                            kch = (nb1i>>(k<<2))&0x7;
                            genecode = golg[nb[kch]];
                            nrfound_s += (!survivalgene && golij) && (((genecode>>((s-1)*ncoding)) & ncodingmask) == ncodingmask) ? 1 : 0;
                            nrfound_b +=  (overwrite || !golij) && (((genecode>>((8+s1)*ncoding)) & ncodingmask) == ncodingmask) ? 1 : 0;
                        }
                        if(survivalgene) {
                            genecode = golg[ij];
                            nrfound_s = golij && (((genecode>>((s-1)*ncoding)) & ncodingmask) == ncodingmask) ? s : 0;
                        }
                
                        if((repscheme&R_1_nb_OR_AND)||(s&0x1)) found = (nrfound_s >  (s>>1)) ? 1 : 0;  // > s/2
                        else                                   found = (nrfound_s >= (s>>1)) ? 1 : 0;  // >= s/2
                        survive = (smask>>s1) & (found ? 1ull : 0ull);
                        
                        if((repscheme&R_1_nb_OR_AND)||(s&0x1)) found = (nrfound_b >  (s>>1)) ? 1 : 0;
                        else                                   found = (nrfound_b >= (s>>1)) ? 1 : 0;
                        birth = (bmask>>s1) & (found ? 1ull : 0ull);
                    }
                    else {                                                          // this code is more efficient than above approach using AND or OR of neighbors
                        if(repscheme&R_1_nb_OR_AND)
                            for (genecode=0ull,k=0;k<s;k++) {                       // decodes genes with fixed length encoding by OR
                                kch = (nb1i>>(k<<2))&0x7;
                                genecode |= golg[nb[kch]];                          // OR of live neighbours encodes birth rule & survival rule
                            }
                        else {
                            for (genecode=allcoding,k=0;k<s;k++) {                  // decodes genes with fixed length encoding by AND
                                kch = (nb1i>>(k<<2))&0x7;
                                genecode &= golg[nb[kch]];                          // AND of live neighbours encodes birth rule & survival rule
                            }
                        }
                        if(survivalgene) genecode = (genecode&(0xffffffffull<<32)) | (golg[ij]&(8ull*ncoding-1ull));   // if central gene codetermines survival
                        survive=(((genecode>>((s-1)*ncoding)) & ncodingmask) == ncodingmask) && ((smask>>s1)&1ull) ? 1ull : 0ull;
                        if (overwrite || !golij) birth=(((genecode>>((8+s1)*ncoding)) & ncodingmask) == ncodingmask) && ((bmask>>s1)&1ull) ? 1ull : 0ull;
                    }
                }
            }
            else if (rulemodij==2){                                          // hard death on membrane defined above via macro "MEMBRANE"
                survive = 0ull;
                birth = 0ull;
            }
            else {
                    survive = s2or3;
                    birth = s2or3&s&0x1ull&~golij;
            }
        }
        else survive = birth = s2or3 = gols = 0ull;
    finish_update_ij(ij,s,golij,gols,nb1i,nbmask,nb,survive,birth,gol,golg,golgstats,golb,golr,newgol,newgolg,newgolgstats,newgolb,newgolr);
    }  // end for ij

    finish_update(newgol, newgolg, newgolgstats, newgolb, newgolr, nbshist);
}
//---------------------------------------------------------------- update_lut_dist ----------------------------------------------------------------------
void update_lut_dist(uint64_t gol[], uint64_t golg[], uint64_t golgstats[], uint64_t golb[],uint64_t golr[],uint64_t newgol[], uint64_t newgolg[], uint64_t newgolgstats[], uint64_t newgolb[],uint64_t newgolr[]) {     // selection models 10,11
// update GoL for toroidal field which has side length which is a binary power of 2
// encode without if structures for optimal vector treatment
/*
    Configurations are distinguished by the number of ones in two classes (corner,edge-centred) of the 8-bit live neighbour pattern
    - for the different s values                 0  1  2  3  4  5  6  7  8
    - there are a nr n of partitions             1  2  3  4  5  4  3  2  1   total 25   with se the number in edge centred sites N,S,E,W
    - chosen nr of partitions for exploration    0  2  3  4  5  4  3  2  0   total 23
    i.e. if we exclude the cases s = 0,8, then there are 23 bits required to distinguish cases of either birth (movement) or survival.
    This results in the following two encodings of genomes: (selected by 1. selection=10 or 2. selection=11):
    1. Fixed length encoding has 1 bit per LUT entry: 23 survival bits and then 23 birth bits
        B/S 1111111111111111111111100000000000000000000000      1=Birth  0=Survival
        s   7766655554444433332221177666555544444333322211
        se  1021032104321032102101010210321043210321021010
    2. Variable length encoding uses two 4-bit bytes : upper byte for B/S bit and 3-bit 0<s<8 and lower byte for bitmask of which se 0-3 are on
        The case s=4,se=4 is treated as a special case x0000001 where x is 1 for B and 0 for S. Note that s=0 is not allowed for survival or birth.
        bit 7    6    5    4    3    2    1    0        of rule aligned with 8-bit genes in genome
            B/S  s_2  s_1  s_0  se=3 se=2 se=1 se=0     where three bits of s are s_2,s_1,s_0
            B/S  0    0    0    0    0    0    1        exceptional encodings for s=4,se=4 case  (not otherwise used)                */
 
    int s, s1, se, s0, s2or3, k, kch, nrfound_b, nrfound_s;
    unsigned int rulemodij;
    int nb[8], ij, i, j, jp1, jm1, ip1, im1;
    uint64_t genecode, gols, golij, nb1i, nbmask, found;
    uint64_t survive, birth, overwrite, survivalgene, smask, bmask, gene;
    
    static uint64_t summasks[7] = {0x3ull,0x7ull,0xfull,0x1full,0xfull,0x7ull,0x3ull};
    static int sumoffs[7] = {0,2,5,9,14,18,21};
    // static int first = 1;

    canonical = repscheme & R_2_canonical_nb;                                      // set global choice of canonical rotation bit choice for selectdifftx
    survivalgene = repscheme & R_12_survivalgene;                                   // gene determining survival is 1: central gene 2: determined by neighbours
    smask = (uint64_t) survivalmask;                                               // convert to 64 bit mask for efficient usage here
    bmask = (uint64_t) birthmask;
    parentdies = (repscheme & R_3_parentdies) ? 1 : 0;
    if(parentdies) for (ij=0; ij<N2; ij++) newgolgstats[ij] = 0ull;                // need to update statistics of neighbours with parenting information, so init required
    
    for (ij=0; ij<N2; ij++) {                                                      // loop over all sites of 2D torus with side length N
        i = ij & Nmask;  j = ij >> log2N;                                          // row & column
        jp1 = ((j+1) & Nmask)*N; jm1 = ((j-1) & Nmask)*N;                          // toroidal (j+1)*N and (j-1)*N
        ip1 =  (i+1) & Nmask; im1 =  (i-1) & Nmask;                                // toroidal i+1, i-1
        nb[0]=jm1+im1; nb[1]=jm1+i; nb[2]=jm1+ip1; nb[3]=j*N+ip1;                  // new order of nbs
        nb[4]=jp1+ip1; nb[5]=jp1+i; nb[6]=jp1+im1; nb[7]=j*N+im1;
        for (s=se=0,nb1i=0ull,nbmask=0ull,k=0;k<8;k++) {                           // packs non-zero nb indices in first up to 8*4 bits
            gols=gol[nb[k]];                                                       // whether neighbor is alive
            s += gols;                                                             // s is number of live nbs
            se += k&0x1&gols;                                                      // se is number of edge-centred live neighbours (odd k)
            nb1i = (nb1i << (gols<<2)) + (gols*k);                                 // nb1i is packed list of live neighbour indices
            nbmask |= (gols << k);
        }
        survive = birth = 0ull;
        golij = gol[ij];
        
        if (s>0 && s<8) {
            s1 = s-1;
            s2or3 = (s>>2) ? 0ull : (s>>1);                                        // s == 2 or s ==3 : checked by bits 2+ are zero and bit 1 is 1
            gols = s2or3 ? (golij ? 1ull : (s&1ull ? 1ull : 0ull )) : 0ull;        // GoL calculation next state for non-genetic gol plane
            rulemodij = (rulemod&0x4) ? MEMBRANE : ((rulemod&0x2) ? (ij>=(N2>>1) ? 1 : 0) : (rulemod&0x1)); // if rulemod bit 2 then activate membrane of death
                                                                                   // else if rulemod bit 1 is on then split into half planes with/without mod
            if (rulemodij==1) {                // NB need to put gene calculation outside so that we can do genetic propagation with GoL rulemod off
                overwrite = overwritemask&(0x1ull<<s1);
                overwrite = (overwrite || !golij) ? 1ull : 0ull;
                if (golij) survive = (smask>>sumoffs[s1])&summasks[s1] ? 1ull : 0ull;
                if (overwrite) birth = (bmask>>sumoffs[s1])&summasks[s1] ? 1ull : 0ull;
                if (survive|birth) {
                    s0 = (s-4 > 0 ? s-4 : 0);
                    if (golij)   survive = (smask>>(sumoffs[s1]+se-s0))&0x1ull;     // refine decisions for specific combination of s and se
                    if (overwrite) birth = (bmask>>(sumoffs[s1]+se-s0))&0x1ull;     // only allowed if birth,survivalmask permits (ask this before consulting genes)
                    if (survive|birth) {                                            // complete determination of birth or survival
                        if (selection==11) {
                            nrfound_b = nrfound_s = 0;
                            for (k=0;k<s;k++)  {                                    // decodes genes with variable position encoding only for current s,se
                                kch = (nb1i>>(k<<2))&0x7;
                                gene = golg[nb[kch]] & (0xf0f0f0f0f0f0f0f0ull | (0x0101010101010101ull<<(se-s0))); // focus gene down to specific required se bit for exact matching
                                if (!survivalgene && golij) {                 // coding is 8 bits [(b/s) (s1 2 1 0) (se subset mask 3 2 1 0)]  (exception case s==4,se==4 is x0000001)
                                    PATTERN8(gene, (s==4 && se==4) ? 0x01 : ((s<<4)|(1<<(se-s0))), found);    //survival rule found? final decision for survival
                                    if(found) nrfound_s++;
                                }
                                if (overwrite) {
                                    PATTERN8(gene, (s==4 && se==4) ? 0x81 :(((8|s)<<4)|(1<<(se-s0))), found); //birth rule found? final decision for birth
                                    if(found) nrfound_b++;
                                }
                            }
                            if (survivalgene && golij) {                            // survival determined by central gene in this case
                                gene = golg[ij] & (0xf0f0f0f0f0f0f0f0ull | (0x0101010101010101ull<<(se-s0))); // focus gene down to specific required se bit for exact matching
                                PATTERN8(gene, (s==4 && se==4) ? 0x01 : ((s<<4)|(1ull<<(se-s0))), found);   // survival rule found? final decision for survival
                                if(found) nrfound_s = s;
                            }

                            if(repscheme&R_0_nb_majority) {
                                if((repscheme&R_1_nb_OR_AND)||(s&0x1)) {
                                    survive &= (nrfound_s >  (s>>1)) ? 1ull : 0ull;  // > s/2
                                    birth   &= (nrfound_b >  (s>>1)) ? 1ull : 0ull;  // > s/2
                                }
                                else  {
                                    survive &= (nrfound_s >= (s>>1)) ? 1ull : 0ull;  // >= s/2
                                    birth   &= (nrfound_b >= (s>>1)) ? 1ull : 0ull;  // >= s/2
                                }
                            }
                            else {
                                if(repscheme&R_1_nb_OR_AND) {
                                    survive &= (nrfound_s     ) ? 1ull : 0ull;   // OR
                                    birth   &= (nrfound_b     ) ? 1ull : 0ull;   // OR
                                }
                                else  {
                                    survive &= (nrfound_s == s) ? 1ull : 0ull;   // AND
                                    birth   &= (nrfound_b == s) ? 1ull : 0ull;   // AND
                                }
                            }
                        }
                        else {                                                      // selection == 10
                            if(repscheme&R_0_nb_majority) {                             // MAJORITY resolution of neighbours
                                for (nrfound_b=nrfound_s=0,k=0;k<s;k++) {               // decodes genes with fixed length encoding
                                    kch = (nb1i>>(k<<2))&0x7;
                                    genecode = golg[nb[kch]];
                                    nrfound_s += (!survivalgene && golij) && ((genecode>>(sumoffs[s1]+(se-s0)))&0x1ull) ? 1 : 0;
                                    nrfound_b +=  (overwrite || !golij) && ((genecode>>(32+sumoffs[s1]+(se-s0)))&0x1ull) ? 1 : 0;
                                }
                                
                                if(survivalgene) {
                                    genecode = golg[ij];
                                    nrfound_s = golij && ((genecode>>(sumoffs[s1]+(se-s0)))&0x1ull) ? s : 0;
                                }
                        
                                if((repscheme&R_1_nb_OR_AND)||(s&0x1)) found = (nrfound_s >  (s>>1)) ? 1 : 0;  // > s/2
                                else                                   found = (nrfound_s >= (s>>1)) ? 1 : 0;  // >= s/2
                                survive &= found ? 1ull : 0ull;
                                
                                if((repscheme&R_1_nb_OR_AND)||(s&0x1)) found = (nrfound_b >  (s>>1)) ? 1 : 0;  // > s/2
                                else                                   found = (nrfound_b >= (s>>1)) ? 1 : 0;  // >= s/2
                                birth &= found ? 1ull : 0ull;
                            }
                            else {
                                if(repscheme & R_1_nb_OR_AND)
                                    for (genecode=0ull,k=0;k<s;k++) {                   // decodes genes with fixed length encoding by OR
                                        kch = (nb1i>>(k<<2))&0x7;
                                        genecode |= golg[nb[kch]];                      // OR of live neighbours encodes birth rule & survival rule
                                    }
                                else
                                    for (genecode=~0ull,k=0;k<s;k++) {                  // decodes genes with fixed length encoding by AND
                                        kch = (nb1i>>(k<<2))&0x7;
                                        genecode &= golg[nb[kch]];                      // AND of live neighbours encodes birth rule & survival rule
                                    }
                                if(survivalgene) genecode = (genecode&(0xffffffffull<32)) | (golg[ij]&0xffffffffull);     // if central gene determines survival
                                genecode &= (bmask << 32)|smask;                        // just to be clear, not required as now already tested above
                                if (golij) survive &= (genecode>>(sumoffs[s1]+(se-s0)))&0x1ull;        // changed from |=
                                if (overwrite) birth &= (genecode>>(32+sumoffs[s1]+(se-s0)))&0x1ull;
                            }
                        }
                    }
                }
            }
            else if (rulemodij==2){                                                 // hard death on membrane defined above via macro "MEMBRANE"
                survive = 0ull;
                birth = 0ull;
            }
            else {                                                                  // GOL rule
                survive = s2or3;
                birth = s2or3&s&0x1&~golij;
            }
            // if(s>2)  survive=0ull;
        }
        else survive = birth = s2or3 = gols = 0ull;

    finish_update_ij(ij,s,golij,gols,nb1i,nbmask,nb,survive,birth,gol,golg,golgstats,golb,golr,newgol,newgolg,newgolgstats,newgolb,newgolr);
    }  // end for ij

    finish_update(newgol, newgolg, newgolgstats, newgolb, newgolr, nbshist);
}
//---------------------------------------------------------------- update_lut_canon_rot -----------------------------------------------------------------
void update_lut_canon_rot(uint64_t gol[], uint64_t golg[], uint64_t golgstats[], uint64_t golb[],uint64_t golr[],uint64_t newgol[], uint64_t newgolg[], uint64_t newgolgstats[], uint64_t newgolb[],uint64_t newgolr[]) {     // selection models 12,13
/*
    Configurations are distinguished by the number of ones and the canonical rotation of the 8-bit live neighbour pattern.
    The canonical rotation is the rotation with the minimum numerical value as 8-bit number (bits are numbered clockwise from 0-7):
    - for the different s values                 0  1  2  3  4  5  6  7  8
    - there are a nr of canonical rotations      1  1  4  7 10  7  4  1  1   total 36
    - if we exclude the cases s=0,1,7,8 which can create growth artefacts, then there are 32 bits
    This results in the following two encodings of genomes: (selected by 1. selection=12 or 2. selection=13):
    1. Fixed length encoding: Survival: s=2 bits 0-3, s=3 bits 4-10, s=4 bits 11-20, s=5 bits 21-27, s=6 bits 28-31; Birth: same + 32
        B/S  1111111111111111111111111111111100000000000000000000000000000000      1=Birth  0=Survival
        s    6666555555544444444443333333222266665555555444444444433333332222
        crot 3210654321098765432106543210321032106543210987654321065432103210
    2. Modular encoding: Up to 6 10-bit modules. Lower 8 bits stored in first 6*8=48 bits. Upper 2-bits stored pairwise in 12 bits 48-59.
       bit 63  xxxxxx989898989898765432107654321076543210765432107654321076543210   bit 0
       Module: 9     8    |  7     6     5     4     3     2     1     0       we abbreviate possible values of crot mod 5 as cr0-4 in this description
               cr4   cr3     s3    s2    s1    s0    crot5 cr2   cr1   cr0     8-bit pattern matching via remapping cr0-4 to lowest 3 bits: 0->001 1->010 2->100 3->001 4->010
                             s3    s2    s1    s0    crot5 0/cr2 cr4/1 cr3/0   remapping of bits depending on whether crot mod 5 is >2 / <=2
                                                                                                                                                                                    */
    int s, smid, s2, s2or3, k, k1, kodd, crot, crot5, crotmod5, pat, nsame, nrfound_b, nrfound_s;
    uint64_t survive, birth, overwrite, survivalgene, smask, bmask, found, rulemodij, golij;
    static uint64_t summasks[5] = {0xfull,0x7full,0x3ffull,0x7full,0xfull};
    static int sumoffs[5] = {0,4,11,21,28};
    int nb[8], ij, i, j, jp1, jm1, ip1, im1;
    unsigned int kch=0;
    uint64_t genecode, genecode1, gols, nb1i, nbmask;

    canonical = repscheme & R_2_canonical_nb;                                       // set global choice of canonical rotation bit choice for selectdifftx
    survivalgene = repscheme & R_12_survivalgene ? 1ull : 0ull;                      // gene determining survival is 1: central gene 2: determined by neighbours
    smask = (uint64_t) survivalmask;                                                // 32 bits of survivalmask used to limit space of rules, convert to 64 bit masks for efficient usage here
    bmask = (uint64_t) birthmask;                                                   // 32 bits of birthmask
    parentdies = (repscheme & R_3_parentdies) ? 1 : 0;
    
    if(parentdies) for (ij=0; ij<N2; ij++) newgolgstats[ij] = 0ull;                 // need to update statistics of neighbours with parenting information, so init required

    for (ij=0; ij<N2; ij++) {                                                       // loop over all sites of 2D torus with side length N
        i = ij & Nmask;  j = ij >> log2N;                                           // row & column
        jp1 = ((j+1) & Nmask)*N; jm1 = ((j-1) & Nmask)*N;                           // toroidal (j+1)*N and (j-1)*N
        ip1 =  (i+1) & Nmask; im1 =  (i-1) & Nmask;                                 // toroidal i+1, i-1
        nb[0]=jm1+im1; nb[1]=jm1+i; nb[2]=jm1+ip1; nb[3]=j*N+ip1;                   // new order of nbs
        nb[4]=jp1+ip1; nb[5]=jp1+i; nb[6]=jp1+im1; nb[7]=j*N+im1;
        for (s=0,nb1i=0ull,nbmask=0ull,k=0;k<8;k++) {                               // packs non-zero nb indices in first up to 8*4 bits
            gols=gol[nb[k]];                                                        // whether neighbor is alive
            s += gols;                                                              // s is number of live nbs
            nb1i = (nb1i << (gols<<2)) + (gols*k);                                  // nb1i is packed list of live neighbour indices
            nbmask |= (gols << k);
        }
        smid = s>1 && s<7; s2 = s-2;                                                // s in mid-range for possible lut rule
        survive = birth = 0ull;
        golij = gol[ij];

        if (smid) {
            s2or3 = (s>>2) ? 0ull : (s>>1);                                         // s == 2 or s ==3 : checked by bits 2+ are zero and bit 1 is 1
            gols = s2or3 ? (gol[ij] ? 1ull : (s&1ull ? 1ull : 0ull )) : 0ull;       // GoL calculation next state for non-genetic gol plane
            
            rulemodij = (rulemod&0x4) ? MEMBRANE : ((rulemod&0x2) ? (ij>=(N2>>1) ? 1 : 0) : (rulemod&0x1)); // if rulemod bit 2 then activate membrane of death
                                                                                    // else if rulemod bit 1 is on then split into half planes with/without mod
            if(rulemodij==1) {
                overwrite = s ? (overwritemask>>(s-1))&0x1ull : 0ull;               // allow birth to overwrite occupied cell = survival in GoL
                overwrite = overwrite | (~gol[ij] & 0x1ull);                        // either central cell is empty or overwrite bit set is required for birth
                if (gol[ij]) survive = (smask>>sumoffs[s2])&summasks[s2] ? 1ull : 0ull;
                if (overwrite) birth = (bmask>>sumoffs[s2])&summasks[s2] ? 1ull : 0ull;
                if (survive|birth) {
                    for (k=0,nbmask=0;k<8;k++) nbmask |= (gol[nb[k]]<<k);           // constuct mask of live bits for 8 neighbours
                    kch = selectdifft(s, nbmask, &crot, &kodd, &nsame);             // find the canonical rotation index of the live neighbour configuration
                    survive&= (smask>>(sumoffs[s2]+crot))&0x1ull;                   // refine decisions for specific canonical rotation configuration
                    birth  &= (bmask>>(sumoffs[s2]+crot))&0x1ull;                   // only allowed if birth/survivemask permits (ask this before consulting genes)
                    if (survive|birth) {                                            // complete determination of birth or survival
                        if (selection==13) {                                        // modular gene encoding one of up to two 5-bit subsets of crot (e.g. for s=3,4 5+2 and 5+5 resp.)
                            nrfound_b = nrfound_s = 0;
                            crot5 = crot > 4 ? 1 : 0;                               // crot is in range 0-9
                            crotmod5 = crot-5*crot5;                                // the 5 bits indexed by crotmod5 are stored in two places 0-2 in 8-bit word and 3-4 in pairs bits 48+
                            pat = (s<<4)|(crot5<<3)|(0x1<<((crotmod5>2)?crotmod5-3:crotmod5));
                            genecode = 0ull;
                            for (k=0;k<s;k++) {                                     // decodes genes with variable position encoding only for current s,crot
                                if ((!survivalgene && golij) || overwrite) {       // replace lowest 3 bits of 8-bit part of 10-bit modules with appropriate crot subset bit
                                    kch = (nb1i>>(k<<2))&0x7;
                                    genecode = golg[nb[kch]];
                                    if(crotmod5 > 2) for (genecode1=0ull,k1=0;k1<6;k1++) genecode1 |= ( (genecode>>(48+(k1<<1))) & (0x1ull<<(crotmod5-3)) ) << (k1<<3);
                                    else genecode1 = genecode & (0x010101010101 << crotmod5); // prepare lookup in all 6 modules on gene : in this case crot subset bits in right place
                                    genecode &= 0xf8f8f8f8f8f8;                    // accept upper 5 bits of lower 8-bits of all 6 coding modules
                                    genecode |= genecode1;                         // combine with lower 3 bits assembled from relocated appropriate crot subset bit
                                }
                                if (!survivalgene && golij) {
                                    PATTERN8(genecode, pat, found);
                                    if(found) nrfound_s++;
                                }
                                if (overwrite) {
                                    PATTERN8(genecode, (0x80|pat), found);
                                    if(found) nrfound_b++;
                                }
                            }
                            if (survivalgene & golij) {                           // survival determined by central gene in this case
                                genecode = golg[ij];
                                if(crotmod5 > 2) for (genecode1=0ull,k1=0;k1<6;k1++) genecode1 |= ((genecode>>(48+(k1<<1)))&(0x1<<(crotmod5-3)))<<(k1<<3);
                                else genecode1 = genecode & (0x010101010101 << crotmod5);
                                genecode &= 0xf8f8f8f8f8f8;
                                genecode |= genecode1;
                                PATTERN8(genecode, pat, found);
                                if(found) nrfound_s = s;
                            }
                            
                            if(repscheme&R_0_nb_majority) {
                                if((repscheme&R_1_nb_OR_AND)||(s&0x1)) {
                                    survive &= (nrfound_s >  (s>>1)) ? 1ull : 0ull;  // > s/2
                                    birth   &= (nrfound_b >  (s>>1)) ? 1ull : 0ull;  // > s/2
                                }
                                else  {
                                    survive &= (nrfound_s >= (s>>1)) ? 1ull : 0ull;  // >= s/2
                                    birth   &= (nrfound_b >= (s>>1)) ? 1ull : 0ull;  // >= s/2
                                }
                            }
                            else {
                                if(repscheme&R_1_nb_OR_AND) {
                                    survive &= (nrfound_s     ) ? 1ull : 0ull;   // OR
                                    birth   &= (nrfound_b     ) ? 1ull : 0ull;   // OR
                                }
                                else  {
                                    survive &= (nrfound_s == s) ? 1ull : 0ull;   // AND
                                    birth   &= (nrfound_b == s) ? 1ull : 0ull;   // AND
                                }
                            }
                        }
                        else {                                                      // selection == 12
                            if(repscheme&R_0_nb_majority) {                             // MAJORITY resolution of neighbours
                                for (nrfound_b=nrfound_s=0,k=0;k<s;k++) {               // decodes genes with fixed length encoding
                                    kch = (nb1i>>(k<<2))&0x7;
                                    genecode = golg[nb[kch]];
                                    nrfound_s += (!survivalgene && golij) && ((genecode>>(sumoffs[s2]+crot))&0x1ull) ? 1 : 0;
                                    nrfound_b +=  (overwrite || !golij) && ((genecode>>(32+sumoffs[s2]+crot))&0x1ull) ? 1 : 0;
                                }
                                
                                if(survivalgene) {
                                    genecode = golg[ij];
                                    nrfound_s = golij && ((genecode>>(sumoffs[s2]+crot))&0x1ull) ? s : 0;
                                }
                        
                                if((repscheme&R_1_nb_OR_AND)||(s&0x1)) found = (nrfound_s >  (s>>1)) ? 1ull : 0ull;  // > s/2
                                else                                   found = (nrfound_s >= (s>>1)) ? 1ull : 0ull;  // >= s/2
                                survive &= found ? 1ull : 0ull;
                                
                                if((repscheme&R_1_nb_OR_AND)||(s&0x1)) found = (nrfound_b >  (s>>1)) ? 1ull : 0ull;  // > s/2
                                else                                   found = (nrfound_b >= (s>>1)) ? 1ull : 0ull;  // >= s/2
                                birth &= found ? 1ull : 0ull;
                            }
                            else {
                                if(repscheme&R_1_nb_OR_AND)
                                    for (genecode=0ull,k=0;k<8;k++)                     // decodes genes with fixed length encoding by OR
                                        genecode |= (gol[nb[k]]?golg[nb[k]]:0ull);      // OR of live neighbours encodes birth rule & survival rule
                                else
                                    for (genecode=~0ull,k=0;k<8;k++)                    // decodes genes with fixed length encoding by AND
                                        genecode &= (gol[nb[k]]?golg[nb[k]]:~0ull);     // AND of live neighbours encodes birth rule & survival rule
                                if(survivalgene) genecode = (genecode&(0xffffffffull<<32)) | (golg[ij]&0xffffffffull);     // if central gene determines survival
                                genecode&=(bmask<<32)|smask;                            // actually no longer needed since test done above
                                if (gol[ij]) survive &= (genecode>>(sumoffs[s2]+crot))&0x1ull;         // changed from |=
                                if (overwrite) birth &= (genecode>>(32+sumoffs[s2]+crot))&0x1ull;
                            }
                        }
                    }
                }
            }
            else if (rulemodij==2){                                                 // hard death on membrane defined above via macro "MEMBRANE"
                survive = 0ull;
                birth = 0ull;
            }
            else {
                survive = s2or3;
                birth = s2or3&s&0x1&~gol[ij];
            }
        }
        else survive = birth = s2or3 = gols = 0ull;
        
    finish_update_ij(ij,s,golij,gols,nb1i,nbmask,nb,survive,birth,gol,golg,golgstats,golb,golr,newgol,newgolg,newgolgstats,newgolb,newgolr);
    }  // end for ij

    finish_update(newgol, newgolg, newgolgstats, newgolb, newgolr, nbshist);
}
//---------------------------------------------------------------- update_lut_2Dsym ---------------------------------------------------------------------
void update_lut_2D_sym(uint64_t gol[], uint64_t golg[], uint64_t golgstats[], uint64_t golb[], uint64_t golr[], uint64_t newgol[], uint64_t newgolg[], uint64_t newgolgstats[], uint64_t newgolb[], uint64_t newgolr[]) {     // selection models 14,15
/*
    All different configurations under the standard 2D 4-rotation and 4-reflection symmetries are distinguished
    i.e. by number of ones and edge-corner differences and additional distinctions in arrangement
    - for the different s values                     0  1  2  3  4  5  6  7  8
    - there are a nr of distinct configurations      1  2  6 10 13 10  6  2  1   total 51
    - if we exclude the cases s=0,1,7,8 which can create growth artefacts, then there are still 45 bits
    - if we focus on cases s=3,4,5 there are still one too many, i.e. 33 cases, for investigation in a single integer genome
    - instead, we investigate here the most interesting lower s domain s = 0,1,2,3,4 in full detail
 
    [Alternative approaches would be to:
        (i) investigate separately the configurations s=4-8 which requires 32 bits
        (ii) investigate s=2,3,5,6, excluding 4, with 32 bits for birth and survival
        (iii) code up the s=3,4,5 minus one : with the software specifying one omitted configuration
        (iv) include undifferentiated s=5 instead of 0 : 32 bits
        (v) include undifferentiated s=1,5,6 instead of 0,1  ie s=1-6 : 32 bits
        (vi) split the s-range allowed for birth and survival: S s=1-4 B s=3-5 ie 31 bits for survival and 33 for birth.]
 
    This results in the following two encodings of genomes: (selected by 1. selection=14 or 2. selection=15):
    1. Fixed length encoding: Survival: s=0 bit 0, s=1 bits 1,2, s=2 bits 3-8, s=3 bits 9-18, s=4 bits 19-31; Birth: same + 32
        B/S  1111111111111111111111111111111100000000000000000000000000000000      1=Birth  0=Survival
        s    4444444444444333333333322222211044444444444443333333333222222110
        crot 32106543210987654321065432103210cba98765432109876543210543210100
    2. Modular encoding: Up to 5 12-bit modules. Lower 8 bits in up to 5 modules in bits 0-39, upper 4 bits in up to 5 modules in bits 40-59
                                                                                                                */
    int s, se, slow, s2or3, k, k1, kodd, coff, crot, coffdiv6, coffmod6, pat, nsame, nrfound_b, nrfound_s;
    uint64_t survive, birth, overwrite, survivalgene, found, smask, bmask, rulemodij, golij;
    static uint64_t summasks[5] = {0x1ull,0x3ull,0x3full,0x3ffull,0x1fffull};  // masks for s= 0,1,2,3,4 with nr cases 1,2,6,10,13
    static int sumoffs[9] = {0,1,3,9,19,32,42,48,50};                          // cumulative offsets to start of coding region for s = 0,1,2,3,4,5,6,7,8 - only s<5 used
    static int csumoffs[9] = {0,2,4,12,26,46,60,68,2};                         // start of indexing in confoffs for crot,kodd lookup for s = 0,1,2,3,4,5,6,7,8
    static unsigned char confoffs[2+2+8+14+20+14+8+2+2] = {0,0, 0,0, 0,0,1,2,3,3,4,5, 0,1,2,3,3,2,4,5,6,7,4,5,8,9, 0,0,1,2,3,4,1,2,5,6,7,8,9,9,10,10,8,7,11,12,
                                        0,1,2,3,3,2,4,5,6,7,4,5,8,9, 0,0,1,2,3,3,4,5, 0,0, 0,0}; // look up for crot*2+kodd to gene bit offset
    int nb[8], ij, i, j, jp1, jm1, ip1, im1;
    unsigned int kch=0;
    uint64_t genecode, genecode1, gols, nb1i, nbmask;

/*
    static int first = 1;
    static unsigned char lut[256];
    if (first) {                                                                // call selectdifft once at start for all possible nbmask vals to accelerate processing
        first = 0;
        for (nbmask=0;nbmask<256;nbmask++) {
            POPCOUNT64C(nbmask,s);
            kch = selectdifft(s, nbmask, &crot, &kodd, &nsame);
            coff = confoffs[csumoffs[s]+(crot<<1)+kodd];
            lut[nbmask]=coff;
            // fprintf(stderr,"LUT nbmask %2llx s %1d goff %2d soff %2d coff %2d crot %2d kodd %1d nsame %2d\n",nbmask,s,sumoffs[s]+coff,sumoffs[s],coff,crot,kodd,nsame);
        }
    }
*/

    canonical = repscheme & R_2_canonical_nb;                                  // set global choice of canonical rotation bit choice for selectdifftx
    survivalgene = repscheme & R_12_survivalgene;                               // gene determining survival is 1: central gene 2: determined by neighbours
    smask = (uint64_t) survivalmask;                                           // 32 bits of survivalmask used to limit space of rules, convert to 64 bit masks for efficient usage here
    bmask = (uint64_t) birthmask;                                              // 32 bits of birthmask
    parentdies = (repscheme & R_3_parentdies) ? 1 : 0;

    if(parentdies) for (ij=0; ij<N2; ij++) newgolgstats[ij] = 0ull;            // need to update statistics of neighbours with parenting information, so init required

    for (ij=0; ij<N2; ij++) {                                                  // loop over all sites of 2D torus with side length N
        i = ij & Nmask;  j = ij >> log2N;                                      // row & column
        jp1 = ((j+1) & Nmask)*N; jm1 = ((j-1) & Nmask)*N;                      // toroidal (j+1)*N and (j-1)*N
        ip1 =  (i+1) & Nmask; im1 =  (i-1) & Nmask;                            // toroidal i+1, i-1
        nb[0]=jm1+im1; nb[1]=jm1+i; nb[2]=jm1+ip1; nb[3]=j*N+ip1;              // new order of nbs
        nb[4]=jp1+ip1; nb[5]=jp1+i; nb[6]=jp1+im1; nb[7]=j*N+im1;
        for (s=se=0,nb1i=0ull,nbmask=0ull,k=0;k<8;k++) {                       // packs non-zero nb indices in first up to 8*4 bits
            gols=gol[nb[k]];                                                   // whether neighbor is alive
            s += gols;                                                         // s is number of live nbs
            se += k&0x1&gols;                                                  // se is number of edge-centred live neighbours (odd k)
            nb1i = (nb1i << (gols<<2)) + (gols*k);                             // nb1i is packed list of live neighbour indices
            nbmask |= (gols << k);
        }
        slow = s<5;                                                            // s in low-range 0-4 for possible lut rule
        survive = birth = 0ull;
        golij = gol[ij];

        if (slow) {
            s2or3 = (s>>2) ? 0ull : (s>>1);                                    // s == 2 or s ==3 : checked by bits 2+ are zero and bit 1 is 1
            gols = s2or3 ? (gol[ij] ? 1ull : (s&1ull ? 1ull : 0ull )) : 0ull;  // GoL calculation next state for non-genetic gol plane
            rulemodij = (rulemod&0x4) ? MEMBRANE : ((rulemod&0x2) ? (ij>=(N2>>1) ? 1 : 0) : (rulemod&0x1)); // if rulemod bit 2 then activate membrane of death
                                                                               // else if rulemod bit 1 is on then split into half planes with/without mod
            if(rulemodij==1) {
                overwrite = s ? overwritemask&(0x1ull<<(s-1)) : 0;
                overwrite = (overwrite || !gol[ij]) ? 1ull : 0ull;
                if (gol[ij]) survive = (smask>>sumoffs[s])&summasks[s] ? 1ull : 0ull;
                if (overwrite) birth = (bmask>>sumoffs[s])&summasks[s] ? 1ull : 0ull;
                if (survive|birth) {
                    for (k=0,nbmask=0;k<8;k++) nbmask |= (gol[nb[k]]<<k);      // constuct mask of live bits for 8 neighbours
                    kch = selectdifft(s, nbmask, &crot, &kodd, &nsame);        // find the canonical rotation index and odd/even offset of the live neighbour configuration
                    coff = confoffs[csumoffs[s]+(crot<<1)+kodd];
                    coffdiv6 = coff >= 12 ? 2 : (coff >=6 ? 1 : 0);            // coff is in range 0-12 coffdiv6 0,1,2
                    coffmod6 = coff-6*coffdiv6;                                // the 6 bits indexed by coffmod6 are stored in two places 0-1 in 8-bit word and 2-5 in quartet bits 40+
                    pat = (s<<4)|(coffdiv6<<2)|0x1;
                    genecode = 0ull;
                    survive= (smask>>(sumoffs[s]+coff))&0x1ull;                // refine decisions for specific canonical rotation configuration
                    birth  = (bmask>>(sumoffs[s]+coff))&0x1ull;                // only allowed if birth/survivemask permits (ask this before consulting genes)
                    if (survive|birth) {                                       // complete the determination of birth or survival
                        if (selection==15) {                                   // selection == 15 variable length encoding
                            nrfound_b = nrfound_s = 0;
                            for (k=0;k<s;k++) {                                // decodes genes with variable position encoding only for current s,crot
                                if ((!survivalgene && gol[ij]) || overwrite) {
                                    kch = (nb1i>>(k<<2))&0x7;
                                    genecode = golg[nb[kch]];
                                    if(coffmod6 > 1) for (genecode1=0ull,k1=0;k1<5;k1++) genecode1 |= ( (genecode>>(40+(k1<<2)+coffmod6-2)) & 0x1ull ) << (k1<<3);
                                    else genecode1 = (genecode >> coffmod6) & 0x0101010101; // prepare lookup in all 5 modules on gene : in this case coff subset bits in right place
                                    genecode &= 0xfcfcfcfcfc;                   // accept upper 6 bits of lower 8-bits of all 5 coding modules
                                    genecode |= genecode1;                      // replace lowest 2 bits of 8-bit part of 12-bit modules with appropriate coff subset bits
                                }
                                if (!survivalgene && golij) {
                                    PATTERN8(genecode, pat, found);
                                    if(found) nrfound_s++;
                                }
                                if (overwrite) {
                                    PATTERN8(genecode, (0x80|pat), found);
                                    if(found) nrfound_b++;
                                }
                            }
                            
                            if (survivalgene & golij) {                     // survival determined by central gene in this case
                                genecode = golg[ij];
                                if(coffmod6 > 1) for (genecode1=0ull,k1=0;k1<5;k1++) genecode1 |= ( (genecode>>(40+(k1<<2)+coffmod6-2)) & 0x1ull ) << (k1<<3);
                                else genecode1 = (genecode >> coffmod6) & 0x0101010101; // prepare lookup in all 5 modules on gene : in this case coff subset bits in right place
                                genecode &= 0xfcfcfcfcfc; //
                                genecode |= genecode1;                         // replace lowest 2 bits of 8-bit part of 12-bit modules with appropriate coff subset bits
                                PATTERN8(genecode, pat, found);                // final decision for survival?  NB all 0 seq encodes survival with 0 live nbs
                                if(found) nrfound_s = s;
                            }
                            
                            if(repscheme&R_0_nb_majority) {
                                if((repscheme&R_1_nb_OR_AND)||(s&0x1)) {
                                    survive &= (nrfound_s >  (s>>1)) ? 1ull : 0ull;  // > s/2
                                    birth   &= (nrfound_b >  (s>>1)) ? 1ull : 0ull;  // > s/2
                                }
                                else  {
                                    survive &= (nrfound_s >= (s>>1)) ? 1ull : 0ull;  // >= s/2
                                    birth   &= (nrfound_b >= (s>>1)) ? 1ull : 0ull;  // >= s/2
                                }
                            }
                            else {
                                if(repscheme&R_1_nb_OR_AND) {
                                    survive &= (nrfound_s     ) ? 1ull : 0ull;   // OR
                                    birth   &= (nrfound_b     ) ? 1ull : 0ull;   // OR
                                }
                                else  {
                                    survive &= (nrfound_s == s) ? 1ull : 0ull;   // AND
                                    birth   &= (nrfound_b == s) ? 1ull : 0ull;   // AND
                                }
                            }
                        }
                        else {                                                 // selection == 14
                            if(repscheme&R_0_nb_majority) {                        // MAJORITY resolution of neighbours
                                for (nrfound_b=nrfound_s=0,k=0;k<s;k++) {          // decodes genes with fixed length encoding
                                    kch = (nb1i>>(k<<2))&0x7;
                                    genecode = golg[nb[kch]];
                                    nrfound_s += (!survivalgene && golij) && ((genecode>>(sumoffs[s]+coff))&0x1ull) ? 1 : 0;
                                    nrfound_b +=  (overwrite || !golij) && ((genecode>>(32+sumoffs[s]+coff))&0x1ull) ? 1 : 0;
                                }
                                
                                if(survivalgene) {
                                    genecode = golg[ij];
                                    nrfound_s = golij && ((genecode>>(sumoffs[s]+coff))&0x1ull) ? s : 0;
                                }
                        
                                if((repscheme&R_1_nb_OR_AND)||(s&0x1)) found = (nrfound_s >  (s>>1)) ? 1 : 0;  // > s/2
                                else                                   found = (nrfound_s >= (s>>1)) ? 1 : 0;  // >= s/2
                                survive &= found ? 1ull : 0ull;
                                
                                if((repscheme&R_1_nb_OR_AND)||(s&0x1)) found = (nrfound_b >  (s>>1)) ? 1 : 0;  // > s/2
                                else                                   found = (nrfound_b >= (s>>1)) ? 1 : 0;  // >= s/2
                                birth &= found ? 1ull : 0ull;
                            }
                            else {
                                if(repscheme&R_1_nb_OR_AND)
                                    for (genecode=0ull,k=0;k<8;k++)                // decodes genes with fixed length encoding by OR
                                        genecode |= (gol[nb[k]]?golg[nb[k]]:0ull); // OR of live neighbours encodes birth rule & survival rule
                                else
                                    for (genecode=~0ull,k=0;k<8;k++)               // decodes genes with fixed length encoding by AND
                                        genecode &= (gol[nb[k]]?golg[nb[k]]:~0ull);// AND of live neighbours encodes birth rule & survival rule
                                if(survivalgene) genecode = (genecode&(0xffffffffull<<32)) | (golg[ij]&0xffffffffull); // if central gene determines survival
                                genecode&=(bmask<<32)|smask;                       // actually no longer needed since test done above
                                if (gol[ij]) survive &= (genecode>>(sumoffs[s]+coff))&0x1ull;        // changed from |=
                                if (overwrite) birth &= (genecode>>(32+sumoffs[s]+coff))&0x1ull;
                            }
                        }
                    }
                }
            }
            else if (rulemodij==2){                                             // hard death on membrane defined above via macro "MEMBRANE"
                survive = 0ull;
                birth = 0ull;
            }
            else {
                survive = s2or3;
                birth = s2or3&s&0x1&~gol[ij];
            }
        }
        else survive = birth = s2or3 = gols = 0ull;

    finish_update_ij(ij,s,golij,gols,nb1i,nbmask,nb,survive,birth,gol,golg,golgstats,golb,golr,newgol,newgolg,newgolgstats,newgolb,newgolr);
    }  // end for ij

    finish_update(newgol, newgolg, newgolgstats, newgolb, newgolr, nbshist);
}
//------------------------------------------------------------------ genelife_update --------------------------------------------------------------------
void genelife_update (int nsteps, int nhist, int nstat) {
    /* update GoL and gene arrays for toroidal field which has side length which is a binary power of 2 */
    /* encode as much as possible without if structures (use ? : instead) in update routines for optimal vector treatment */
    int k,t;
    uint64_t *newgol, *newgolg, *newgolgstats, *newgolb, *newgolr;
    genedata *genedatap = NULL;
    int totalpoptrace(uint64_t gol[]);                                        // calculate total current population and store in scrolling trace npopulation
    int activitieshash(void);                                                 // count activities of all currently active gene species
    int activitieshashquad(void);                                             // count activities of all currently active quad pattern species
    int get_genealogies(genedata genealogydat[], int narraysize);             // genealogies of all currently active species
    int clonealogies(void);                                                   // clonealogies of all currently active clones
    void tracestats(uint64_t gol[],uint64_t golg[], uint64_t golgstats[], int NN2); // trace statistics based on gol,golg
    void countconfigs(void);
    void countspecies1(uint64_t gol[], uint64_t golg[], int N2);              // count species

    nhistG = nhist;                                                           // intervals for collecting histograms
    nstatG = nstat;
    
    for (t=0; t<nsteps; t++) {                                                // main iteration loop for nsteps
        newgol = planes[newPlane];
        newgolg = planesg[newPlane];
        newgolgstats = planesgs[newPlane];
        newgolb = planesb[newPlane];
        newgolr = planesr[newPlane];
        
        totsteps++;                                                           // simulation step counter
        totdisp++;                                                            // currently every step is counted for display in activities

        if (selection<8)        update_23(gol,golg,golgstats,golb,golr,newgol,newgolg,newgolgstats,newgolb,newgolr);           // calculate next iteration with detailed variants of version s=2-3
        else if (selection<10)  update_lut_sum(gol,golg,golgstats,golb,golr,newgol,newgolg,newgolgstats,newgolb,newgolr);      // calculate next iteration for lut sum (gene coded)   version s=1-8
        else if (selection<12)  update_lut_dist(gol,golg,golgstats,golb,golr,newgol,newgolg,newgolgstats,newgolb,newgolr);     // calculate next iteration for lut dist (corner/edge) version s=1-7
        else if (selection<14)  update_lut_canon_rot(gol,golg,golgstats,golb,golr,newgol,newgolg,newgolgstats,newgolb,newgolr);// calculate next iteration for lut canonical rotation version s=2-6
        else if (selection<16)  update_lut_2D_sym(gol,golg,golgstats,golb,golr,newgol,newgolg,newgolgstats,newgolb,newgolr);   // calculate next iteration for lut fully 2D symmetric version s=0-4

        if((diagnostics & diag_offset_statistics) && nhist && (totsteps%nhist == 0)) countconfigs(); // count configurations
        if((diagnostics & diag_general_statistics) && nstat && (totsteps%nstat == 0)) tracestats(gol,golg,golgstats,N2); // time trace point

        curPlane = (curPlane+1) % numPlane;                                   // update plane pointers to next cyclic position
        newPlane = (newPlane+1) % numPlane;
        gol = planes[curPlane];                                               // get planes of gol,golg,golb,golr,golgstats data
        golg = planesg[curPlane];
        golgstats = planesgs[curPlane];
        golb = planesb[curPlane];
        golr = planesr[curPlane];
        
        if (diagnostics & diag_scrolling_trace) totalpoptrace(gol);    // calculate total current population and store in scrolling trace npopulation
        
        if ((diagnostics & diag_activities) && (diagnostics & diag_hash_genes)) {
            nspeciesgene=activitieshash();                                    // colors acttrace and sets current population arrays, need to run always for continuity
            if(nspeciesgene<0) fprintf(stderr,"error returned from activitieshash\n");
        }
        if ((diagnostics & diag_activities ) && (diagnostics & diag_hash_patterns)) {
            nspeciesquad=activitieshashquad();                                 // colors acttraceq and sets current population arrays, need to run always for continuity
            if(nspeciesquad<0) fprintf(stderr,"error returned from activitieshashquad\n");
        }
        if(colorfunction==6 || colorfunction==7 || colorfunction2==6 || colorfunction2==7) { // genealogies
            ngenealogydeep=get_genealogies(genedatap,0);                       // calculates and colors genealogytrace
            if(ngenealogydeep<0) fprintf(stderr,"error returned from genealogies\n");
            nclonealogydeep=clonealogies();                                    // calculates and colors clonealogytrace
            if(nclonealogydeep<0) fprintf(stderr,"error returned from clonealogies\n");
        }
        if(!(totsteps%10) && colorupdate1) {
            if (diagnostics & diag_hash_genes)    nallspecies     = hashtable_count(&genetable);
            if (diagnostics & diag_hash_patterns) nallspeciesquad = hashtable_count(&quadtable);
            if (diagnostics & diag_hash_clones)   nallclones = hashtable_count(&clonetable);
            fprintf(stderr,"step %6d:",totsteps);
            if(diagnostics & diag_hash_genes) {
                countspecies1(gol, golg, N2);
                fprintf(stderr," genes %d/%d (extant/all)",nspeciesgene,nallspecies);
            }
            if(diagnostics & diag_hash_patterns) {
                fprintf(stderr," patterns %d/%d (extant/all)",nspeciesquad+nspeciessmall,nallspeciesquad+nallspeciessmall);
            }
            if (diagnostics & diag_hash_clones) {
                fprintf(stderr," clones %d (all)",nallclones);
            }
            fprintf(stderr," ambiguous state s=4 frequency %d",ambigsum);
            fprintf(stderr,"\n");
            fprintf(stderr," neighbour occupation statistics at birth (NW,N,NE,E,SE,S,SW,W):"); for(k=0;k<8;k++) fprintf(stderr," %d",nbshist[k]); fprintf(stderr,"\n");
            fprintf(stderr,"__________________________________________________________________________________________________________________________________________\n");
        }

    }
}
