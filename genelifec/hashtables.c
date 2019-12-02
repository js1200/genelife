//
//  hashtables.c
//  project genelife
//  management of hash tables for genelife, involving, genes, patterns and clones
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
//...................................................... hash table management for genes and clones .....................................................
// hashaddgene          add new gene to hash table, increment popln and ancestor information for genes already encountered. If mutation calls hashaddclone
// hashdeletegene       decrements population count for gene, printing error message if not found or already zero. Calls hashdeletefromclone
// hashgeneextinction   record gene extinctions in hash gene table, counting number of extinctions
// hashgeneactivity     update activity of gene
// hashaddclone         add new clone to hash table, increment popln information for existing clones
// hashdeletefromclone  decrements population count for clone, printing error message if not found or already zero
// hashcloneactivity    update activity of clone
//......................................................  pattern storage and analysis future ..........................................................
// rotate16             rotate bits in 4x4 pattern for 90 deg clockwise rotation
// rotate64             rotate bits in 8x8 pattern for 90 deg clockwise rotation
// rotate4x64           rotate bits in 16x16 pattern for 90 deg clockwise rotation
// rotatequad           rotate bits in quad pattern for 90 deg clockwise rotation
//...................................................... hash table pattern storage and analysis current ...............................................
// hash_patt8_store     8x8 bit pattern store
// hash_patt4_find      4x4 bit pattern lookup
// hash_patt8_find      8x8 bit pattern lookup
// patt_hash            hash function for a pattern specified by 4 64-bit (8x8) patterns
// node_hash            hash function for a node specified by 4 64-bit pointers
// hash_patt16_store    store new pattern in small pattern table
// hash_patt16_find     find quadtree hash for pattern (leaf of quadtree consists of 4 64bit integers defining a 16x16 pixel array)
// hash_node_store      store new node with hashkey h and subnodes nw,ne,sw,se in hash table
// hash_node_find       find quadtree hash for node (node is specified by its four quadrant pointers (64 bit))
// quadimage            construct quadtree for an entire image or connected component, reporting if the image has been found previously, returning hashkey
// labelimage           rebuild image in a chosen label array from quadimage at chosen offset with chosen label
//------------------------------------------------------------ hash gene inline fns ---------------------------------------------------------------------
extern INLINE void hashaddgene(int ij,uint64_t gene,uint64_t ancestor,uint64_t *golb,uint64_t parentid,uint64_t mutation) {
    genedata gdata;
    uint64_t birthid;
    extern INLINE void hashaddclone(uint64_t birthid, uint64_t parentid, uint64_t gene);
    
    if((genedataptr = (genedata *) hashtable_find(&genetable, gene)) != NULL) {
        genedataptr->lasttime = totsteps;
        if(mutation) {
            genedataptr->recenttime = (short unsigned int) totsteps;
        }
        genedataptr->popcount++;
    }
    else {
        gdata=ginitdata;
        gdata.gene = gene;
        gdata.firsttime = gdata.recenttime = gdata.lasttime = (short unsigned int) totsteps;
        gdata.firstancestor = ancestor;
        hashtable_insert(&genetable, gene,(genedata *) &gdata);
    }
    if(ancestor != rootgene) {
        if((genedataptr = (genedata *) hashtable_find(&genetable, ancestor)) == NULL) {
            fprintf(stderr,"error in hashaddgene step %d ij %d, the ancestor %"PRIx64" of gene %"PRIx64" to be stored is not stored (mutation %"PRIx64")\n",totsteps,ij,ancestor,gene,mutation);
        }
    }
    if(diagnostics & diag_hash_clones) {
        if (mutation || (parentid&rootclone)) {                   // when s=0 birth involves the creation of a new clone from root and rootclone bit set in parentid
            birthid = ((uint64_t) totsteps)<<32;
            birthid |= ij;
            hashaddclone(birthid,parentid,gene);
            *golb = birthid;
        }
        else {
            *golb = parentid;
            if((clonedataptr = (clonedata *) hashtable_find(&clonetable, parentid)) != NULL) clonedataptr->popln++;
            else {
                fprintf(stderr,"step %d error in hashclone update, %"PRIx64" clone not saved\n",totsteps,parentid);
                fprintf(stderr,"ij %"PRIx64" rootclone %d tstep %"PRId64"\n",parentid&N2mask,(parentid&rootclone)?1:0,parentid>>32);
            }
        }
    }
}
//.......................................................................................................................................................
extern INLINE void hashdeletegene(uint64_t gene,uint64_t birthid,const char errorformat[]) {
    extern INLINE void hashdeletefromclone(uint64_t birthid);
    
    if((genedataptr = (genedata *) hashtable_find(&genetable, gene)) != NULL) {genedataptr->popcount--;}
    else fprintf(stderr,errorformat,totsteps,gene);     // errorformat must contain %d and %"PRIx64" format codes in this order
    
    if(diagnostics & diag_hash_clones) hashdeletefromclone(birthid);
}
//.......................................................................................................................................................
extern INLINE void hashgeneextinction(uint64_t gene,const char errorformat[]) {
    if((genedataptr = (genedata *) hashtable_find(&genetable, gene)) != NULL) {
        if(genedataptr->popcount == 0) {
            genedataptr->lastextinctiontime = (short int) totsteps;
            genedataptr->nextinctions++;
        }
    }
    else fprintf(stderr,errorformat,3,totsteps,gene);
}
//.......................................................................................................................................................
extern INLINE void hashgeneactivity(uint64_t gene, const char errorformat[]) {
        if((genedataptr = (genedata *) hashtable_find(&genetable, gene)) != NULL) genedataptr->activity++;
        else fprintf(stderr,errorformat,4,totsteps,gene);
}
//------------------------------------------------------------ hash clone inline fns --------------------------------------------------------------------
extern INLINE void hashaddclone(uint64_t birthid, uint64_t parentid, uint64_t gene) {
    clonedata cdata;

    if((clonedataptr = (clonedata *) hashtable_find(&clonetable, birthid)) != NULL) {
        fprintf(stderr,"error in hashaddclone, %"PRIx64" already present\n",birthid);
    }
    else {
        cdata=cinitdata;
        cdata.birthid=birthid;
        cdata.parentid=parentid;
        cdata.gene = gene;
        hashtable_insert(&clonetable, birthid,(clonedata *) &cdata);
    }
}
//.......................................................................................................................................................
extern INLINE void hashdeletefromclone(uint64_t birthid) {
    if((clonedataptr = (clonedata *) hashtable_find(&clonetable, birthid)) != NULL) {
        if(clonedataptr->popln) clonedataptr->popln--;
        else fprintf(stderr,"error in deleting individual from clone: popln already zero\n");
        // if (!clonedataptr->popln) hashtable_remove( &clonetable, birthid );
    }
    else fprintf(stderr,"error deleting from clone with birthid %"PRIx64" at time %d : not found in hash table\n",birthid,totsteps);
}
//.......................................................................................................................................................
extern INLINE void hashcloneactivity(uint64_t birthid, const char errorformat[]) {
        if((clonedataptr = (clonedata *) hashtable_find(&clonetable, birthid)) != NULL) clonedataptr->activity++;
        else fprintf(stderr,errorformat,6,totsteps,birthid);
}
//------------------------------------------------------- hash quadtree inline fns ----------------------------------------------------------------------
extern INLINE uint16_t rotate16(uint16_t patt) {                                        // rotate bits in 4x4 pattern for 90 deg clockwise rotation
    uint16_t rotate4[16]={0,2,8,10,1,3,9,11,4,6,12,14,5,7,13,15};
    uint16_t nw,ne,sw,se;
    
    nw=patt&0xf;patt>>=4; nw=rotate4[nw];
    ne=patt&0xf;patt>>=4; ne=rotate4[ne];
    sw=patt&0xf;patt>>=4; sw=rotate4[sw];
    se=patt&0xf;          se=rotate4[se];
    return sw | (nw<<4) | (se<<8) | (ne<<12) ;
}
//.......................................................................................................................................................
extern INLINE uint64_t rotate64(uint64_t patt) {                                        // rotate bits in 8x8 pattern for 90 deg clockwise rotation
    uint64_t nw,ne,sw,se;
    
    nw=patt&0xff;patt>>=16; nw=(uint64_t) rotate16((uint16_t) nw);
    ne=patt&0xff;patt>>=16; ne=(uint64_t) rotate16((uint16_t) ne);
    sw=patt&0xff;patt>>=16; sw=(uint64_t) rotate16((uint16_t) sw);
    se=patt&0xff;           se=(uint64_t) rotate16((uint16_t) se);
    return sw | (nw<<16) | (se<<32) | (ne<<48) ;
}
//.......................................................................................................................................................
extern INLINE void rotate4x64(uint64_t *nw, uint64_t *ne, uint64_t *sw, uint64_t *se) { // rotate bits in 16x16 pattern for 90 deg clockwise rotation
    *nw=rotate64(*sw); *ne=rotate64(*sw); *sw=rotate64(*se); *se=rotate64(*ne);
}
//.......................................................................................................................................................
extern INLINE void rotatequad(uint64_t *nw, uint64_t *ne, uint64_t *sw, uint64_t *se) { // rotate bits in quad pattern for 90 deg clockwise rotation
// NYI 1. hash_node(*nw,*ne,*sw,*se) 2. lookup hash entry of hash & check id 3. if patt then call rotate4x64 else do 4 calls to rotatequad with subnodes
}
//.......................................................................................................................................................
extern INLINE quadnode * hash_patt8_store(const uint64_t h, const uint64_t patt) { //8x8 bit pattern store
    int nr1;
    quadinit.hashkey = h;
    quadinit.isnode=0;
    quadinit.nw=patt;
    quadinit.ne=0;
    quadinit.sw=0;
    quadinit.se=0;
    quadinit.size= 8;
    quadinit.firsttime=totsteps;
    quadinit.lasttime=totsteps;
    POPCOUNT64C(patt,nr1);
    quadinit.pop1s =nr1;
    hashtable_insert(&quadtable, h,(quadnode *) &quadinit);
    return (quadnode *) hashtable_find(&quadtable, h);
}
//.......................................................................................................................................................
extern INLINE void hash_patt4_find(const short unsigned int patt) {  // 4x4 bit pattern lookup

    if(smallpatts[patt].activity) {                                         // pattern found
        smallpatts[patt].activity++;
        smallpatts[patt].lasttime=totsteps;
    }
    else {                                                                  // store new pattern
        // smallpatts[patt].size=log2upper(patt);
        smallpatts[patt].activity++;
        smallpatts[patt].firsttime=totsteps;
        smallpatts[patt].lasttime=totsteps;
    }
}
//.......................................................................................................................................................
extern INLINE quadnode * hash_patt8_find(const uint64_t patt) {  // 8x8 bit pattern lookup
        quadnode *q;
        uint64_t h,npatt;
        short unsigned int nw,ne,sw,se;
        const uint64_t randomizer = 11400714819323198549ull;

        h = patt_hash(patt,patt+1,patt+2,patt+3);
        if((q = (quadnode *) hashtable_find(&quadtable, h)) != NULL) {
                                                        // check if pattern found in hash table is correct
            if( patt == q->nw &&  0ull == q->ne && 0ull == q->sw && 0ull ==q->se && !q->isnode) {
                q->activity++;q->lasttime=totsteps;
            }
            else {                                      // collision in hash table at leaf level
                npatt=patt*randomizer;
                h = patt_hash(npatt,npatt+1,npatt+2,npatt+3);
                if((q = (quadnode *) hashtable_find(&quadtable, h)) != NULL) {
                                                        // check if pattern found in hash table is correct
                    if( patt == q->nw &&  0ull == q->ne && 0ull == q->sw && 0ull ==q->se && !q->isnode) {
                        q->activity++;q->lasttime=totsteps;
                    }
                    else {                                      // collision in hash table at 8-leaf level
                        quadcollisions++;
                        fprintf(stderr,"at %d quadhash 2ndary 8-pattern collision %"PRIx64" hash %"PRIx64" collides with %"PRIx64" %"PRIx64" %"PRIx64" %"PRIx64"\n",
                                totsteps,patt,h,q->nw,q->ne,q->sw,q->se);
                    }
                }
                else {                                           // new node or pattern, save in hash table
                    q = hash_patt8_store(h,patt);
                }
            }
        }
        else {                                           // new node or pattern, save in hash table
            q = hash_patt8_store(h,patt);
        }
        nw = patt & 0xffff;ne = (patt>>16) & 0xffff;sw = (patt>>32) & 0xffff; se = (patt>>48) & 0xffff;
        if(nw) hash_patt4_find(nw);   // find or store 8x8 64-bit subpatterns, updating activities and lasttime
        if(ne) hash_patt4_find(ne);   // store if new, otherwise update
        if(sw) hash_patt4_find(sw);
        if(se) hash_patt4_find(se);
        return(q);
}
//.......................................................................................................................................................
extern INLINE uint64_t patt_hash(const uint64_t a, const uint64_t b, const uint64_t c, const uint64_t d) {
                                                        // this hash function seems to work much better than that used in Golly for example
    uint64_t a1,b1,c1,d1,r;
    a1 = a^0x7f0e1d2c3b4a5968ull;
    b1 = b^0xf0e1d2c3b4a59687ull;
    c1 = c^0xba9876543210fedcull;
    d1 = d^0x456789abcdef0123ull;
    r =  (a1>>13)+(b1>>23)+(c1>>29)+(d1>>31);
    r += ((d1<<16)+(c1<<8)+(b1<<4)+(a1<<2)) + (a1+b1+c1+d1);
    r += (r >> 11);
    // r=r*11400714819323198549ull;
    return r ;
}
//.......................................................................................................................................................
extern INLINE uint64_t node_hash(const uint64_t a, const uint64_t b, const uint64_t c, const uint64_t d) {
// now not used : could use instead of patt_hash for hash_node_find
   uint64_t r = (65537*(d)+257*(c)+17*(b)+5*(a));
   r += (r >> 11);
   return r ;
}
//.......................................................................................................................................................
extern INLINE quadnode * hash_patt16_store(const uint64_t h, const uint64_t nw, const uint64_t ne, const uint64_t sw, const uint64_t se) {
    int nr1,nr1s;
    quadinit.hashkey = h;
    quadinit.isnode=0;
    quadinit.nw=nw;
    quadinit.ne=ne;
    quadinit.sw=sw;
    quadinit.se=se;
    quadinit.size= 16;
    quadinit.firsttime=totsteps;
    quadinit.lasttime=totsteps;
    POPCOUNT64C(nw,nr1);nr1s=nr1;
    POPCOUNT64C(ne,nr1);nr1s+=nr1;
    POPCOUNT64C(sw,nr1);nr1s+=nr1;
    POPCOUNT64C(se,nr1);nr1s+=nr1;
    quadinit.pop1s =nr1s;
    hashtable_insert(&quadtable, h,(quadnode *) &quadinit);
    return (quadnode *) hashtable_find(&quadtable, h);
}
//.......................................................................................................................................................
extern INLINE quadnode * hash_patt16_find(const uint64_t nw, const uint64_t ne, const uint64_t sw, const uint64_t se) {
        quadnode *q;
        uint64_t h,nnw,nne,nsw,nse;
        const uint64_t randomizer = 11400714819323198549ull;

        h = patt_hash(nw,ne,sw,se);
        if((q = (quadnode *) hashtable_find(&quadtable, h)) != NULL) {
                                                        // check if pattern found in hash table is correct
            if( nw == q->nw &&  ne == q->ne && sw == q->sw && se ==q->se && !q->isnode) {
                q->activity++;q->lasttime=totsteps;
            }
            else {                                      // collision in hash table at leaf level
                // quadcollisions++;
                // fprintf(stderr,"at %d quadhash pattern collision %"PRIx64" %"PRIx64" %"PRIx64" %"PRIx64" hash %"PRIx64" collides %"PRIx64" %"PRIx64" %"PRIx64" %"PRIx64"\n",
                //    totsteps,nw,ne,sw,se,h,q->nw.patt,q->ne.patt,q->sw.patt,q->se.patt);
                nnw=nw*randomizer;
                nne=ne*randomizer;
                nsw=sw*randomizer;
                nse=se*randomizer;
                h = patt_hash(nnw,nne,nsw,nse);
                if((q = (quadnode *) hashtable_find(&quadtable, h)) != NULL) {
                                                        // check if pattern found in hash table is correct
                    if( nw == q->nw &&  ne == q->ne && sw == q->sw && se ==q->se && !q->isnode) {
                        q->activity++;q->lasttime=totsteps;
                    }
                    else {                                      // collision in hash table at leaf level
                        quadcollisions++;
                        fprintf(stderr,"at %d quadhash 2ndary pattern collision %"PRIx64" %"PRIx64" %"PRIx64" %"PRIx64" hash %"PRIx64" collides %"PRIx64" %"PRIx64" %"PRIx64" %"PRIx64"\n",
                                totsteps,nw,ne,sw,se,h,q->nw,q->ne,q->sw,q->se);
                    }
                }
                else {                                          // new node or pattern, save in hash table
                    q = hash_patt16_store(h,nw,ne,sw,se);
                }
            }
        }
        else {                                                  // new node or pattern, save in hash table
            q = hash_patt16_store(h,nw,ne,sw,se);
        }
        if(nw) (void) hash_patt8_find(nw);                        // find or store 8x8 64-bit subpatterns, updating activities and lasttime
        if(ne) (void) hash_patt8_find(ne);                        // store if new, otherwise update.
        if(sw) (void) hash_patt8_find(sw);
        if(se) (void) hash_patt8_find(se);
        return(q);
}
//.......................................................................................................................................................
extern INLINE quadnode * hash_node_store(uint64_t h, uint64_t nw, uint64_t ne, uint64_t sw, uint64_t se) {
    quadnode *q;
    quadinit.hashkey = h;
    quadinit.isnode=1;
    quadinit.nw=nw;quadinit.ne=ne;quadinit.sw=sw;quadinit.se=se;
    quadinit.firsttime=totsteps;quadinit.lasttime=totsteps;
    quadinit.pop1s=0;quadinit.size=0;                           // incremented below
    quadinit.activity=1;quadinit.topactivity=0;                 // incremented only in quadimage at top level
    
    if (nw) {
        if ((q = (quadnode *) hashtable_find(&quadtable, nw)) == NULL) fprintf(stderr,"Error nw node not found in hashtable.\n");
        else {quadinit.pop1s+=q->pop1s;quadinit.size=q->size<<1;}
    }
    if (ne) {
        if ((q = (quadnode *) hashtable_find(&quadtable, ne)) == NULL) fprintf(stderr,"Error ne node not found in hashtable.\n");
        else {quadinit.pop1s+=q->pop1s;quadinit.size=q->size<<1;}
    }
    if (sw) {
        if((q = (quadnode *) hashtable_find(&quadtable, sw)) == NULL) fprintf(stderr,"Error sw node not found in hashtable.\n");
        else {quadinit.pop1s+=q->pop1s;quadinit.size=q->size<<1;}
    }
    if (se) {
        if((q = (quadnode *) hashtable_find(&quadtable, se)) == NULL) fprintf(stderr,"Error se node not found in hashtable.\n");
        else {quadinit.pop1s+=q->pop1s;quadinit.size=q->size<<1;}
    }
    
    hashtable_insert(&quadtable, h,(quadnode *) &quadinit);
    return (quadnode *) hashtable_find(&quadtable, h);
}
//.......................................................................................................................................................
extern INLINE quadnode * hash_node_find(const uint64_t nw, const uint64_t ne, const uint64_t sw, const uint64_t se) {
                                                                // should only be called if arguments are hashtable keys, not bit patterns
        quadnode *q;
        uint64_t h,nnw,nne,nsw,nse;
    
        h = node_hash(nw,ne,sw,se);                             // alternatively use patt_hash : but this is optimized for 64 bit keys (unlikely to have 0 or low values)
        if((q = (quadnode *) hashtable_find(&quadtable, h)) != NULL) {
            if(nw == q->nw && ne == q->ne && sw == q->sw && se == q->se && q->isnode) { // node found in hash table
                q->activity++;q->lasttime=totsteps;
                // Still need to update activities of 4 subnodes if non-zero !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            }
            else {                                              // collision in hash table at node level
                nnw=nw * 11400714819323198549ull;
                nne= ne * 11400714819323198549ull;
                nsw=(uint64_t) sw * 11400714819323198549ull;
                nse=(uint64_t) se * 11400714819323198549ull;
                h = node_hash(nnw,nne,nsw,nse);                 // alternatively use patt_hash
                if((q = (quadnode *) hashtable_find(&quadtable, h)) != NULL) {
                                                                // check if pattern found in hash table is correct
                    if(nw == q->nw && ne == q->ne && sw == q->sw && se == q->se && q->isnode) { // node found in hash table
                        q->activity++;q->lasttime=totsteps;
                    }
                    else {                                      // collision in hash table at leaf level
                        quadcollisions++;
                        fprintf(stderr,"at %d quadhash node 2ndary collision %"PRIx64" %"PRIx64" %"PRIx64" %"PRIx64" hash %"PRIx64" collides %"PRIx64" %"PRIx64" %"PRIx64" %"PRIx64"\n", totsteps,
                                 nw, ne, sw, se, h, q->nw, q->ne,  q->sw, q->se);  // simple recording for now, later do further chaining or whatever
                    }
                }
                else q=hash_node_store(h,nw,ne,sw,se);          // new node, save in hash table
            }
        }
        else q=hash_node_store(h,nw,ne,sw,se);                  // new node, save in hash table

        return(q);
}
//.......................................................................................................................................................
uint64_t quadimage(uint64_t gol[], short unsigned int *patt, int log2n) {           // routine to generate a quadtree for an entire binary image of long words
                                                                                    // makes use of global linear and quadratic size variable N2 for sizing internal arrays golp,golq
                                                                                    // assumes that n is power of 2
    unsigned int ij,ij1,n3;
    uint64_t golp[N2>>6];
    quadnode * golq[N2>>8];
    int n = 1 << log2n;
    extern INLINE short unsigned int pack16neighbors(uint64_t gol[],int log2n);
    extern INLINE void pack64neighbors(uint64_t gol[],uint64_t golp[],int log2n);
    
    if(n<16) {                                                                      // n < 16
        if (n==8) {                                                                 // n == 8
            pack64neighbors(gol,golp,log2n);
            golq[0]=hash_patt8_find(golp[0]);
            golq[0]->topactivity++;
            return(golq[0]->hashkey);
        }
        else {                                                                      // n == 4,2,1   use smallpatts array to store patterns (more efficient than continued quadtree)
            *patt=pack16neighbors(gol,log2n);
            hash_patt4_find(*patt);
            smallpatts[*patt].topactivity++;
            return 0ull;                                                            // in this case image key is returned in *patt rather than quad key
        }
    }
    else  {                                                                         // n >= 16
        pack64neighbors(gol,golp,log2n);                                            // 8x8 blocks of gol pixels packed into single 64bit words in golp
        n3=n>>3;                                                                    // n3=n/8 is number of such 8x8 blocks along each side of square : n3 is at least 2 here (n>=16)
        // for(ij=0;ij<n3*n3;ij++) { if (ij%8 == 0) fprintf(stderr,"\n step %d ij %d",totsteps,ij);fprintf(stderr," %"PRIx64" ",golp[ij]);} fprintf(stderr,"\n");
        quadtable.expansion_frozen = 1;                                             // freeze quad hash table against expansion ( to ensure valid pointers during array ops)
        for (ij=ij1=0;ij<n3*n3;ij+=2,ij1++) {                                       //  hash all 16x16 patterns (2x2 of golp words) found as leaves of the quadtree
            golq[ij1]=hash_patt16_find(golp[ij],golp[ij+1],golp[(ij+n3)],golp[(ij+n3)+1]); // hash_patt16_find(nw,ne,sw,se) adds quad leaf entry if pattern not found
            ij+= ((ij+2)&(n3-1)) ? 0 : n3;                                          // skip odd rows since these are the northern parts of quads generated on even rows
        }

        for(n3 >>= 1; n3>1; n3 >>= 1) {                                             // proceed up the hierarchy amalgamating 2x2 blocks to next level until reach top
            for (ij=ij1=0;ij<n3*n3;ij+=2,ij1++) {                                   // hash_node_find(nw,ne,sw,se) adds quad node entry if node not found & updates activities
                golq[ij1]=hash_node_find(golq[ij]->hashkey,golq[ij+1]->hashkey,golq[ (ij+n3)]->hashkey,golq[(ij+n3)+1]->hashkey);
                ij+= ((ij+2)&(n3-1)) ? 0 : n3;                                      // skip odd rows since these are the southern parts of quads generated on even rows
            }
        }
        // if(golq[0]!=NULL) if(golq[0]->activity > 1) fprintf(stderr,"step %d image already found at t = %d activity %d\n",totsteps,golq[0]->firsttime,golq[0]->activity);
        golq[0]->topactivity++;
        quadtable.expansion_frozen = 0;                                                         // unfreeze quad hash table
        return(golq[0]->hashkey);
    }
}
//.......................................................................................................................................................
int labelimage(uint64_t hashkeypatt, unsigned int labelimg[], unsigned int label, int offset) { // rebuild image from quadimage at offset with label
    short unsigned int patt;
    int n;
    quadnode * q;
    // uint64_t nw,ne,sw,se;
    extern INLINE void unpack16neighbors(const short unsigned golpw, unsigned int labelimg[], const unsigned int label, const int offset);
    extern INLINE void unpack64neighbors(const uint64_t golpw, unsigned int labelimg[], const unsigned int label, const int offset);
    if(hashkeypatt<65536) {patt = hashkeypatt; unpack16neighbors(patt,labelimg,label,offset);}
    else if((q = (quadnode *) hashtable_find(&quadtable, hashkeypatt)) != NULL) {
        if(q->isnode) {                                     // hashed item is regular quadnode
            n = q->size >> 1;
            if (q->nw) labelimage(q->nw, labelimg, label, offset);
            if (q->ne) labelimage(q->ne, labelimg, label, offset-(offset&Nmask)+((offset+n)&Nmask));
            if (q->sw) labelimage(q->sw, labelimg, label, (offset+n*N)&N2mask);
            if (q->se) labelimage(q->se, labelimg, label, (offset-(offset&Nmask)+((offset+n)&Nmask)+n*N)&N2mask);
        }
        else {                                              // hashed item is pattern
            n=8;
            if (q->nw) {if (q->nw <65536) {patt = q->nw; unpack16neighbors(patt,labelimg,label,offset);}
                        else unpack64neighbors(q->nw,labelimg,label,offset);}
            if (q->ne) {if (q->ne< 65536) {patt = q->ne; unpack16neighbors(patt,labelimg,label,offset-(offset&Nmask)+((offset+n)&Nmask));}
                        else unpack64neighbors(q->ne,labelimg,label,offset-(offset&Nmask)+((offset+n)&Nmask));}
            if (q->sw) {if (q->sw< 65536) {patt = q->sw; unpack16neighbors(patt,labelimg,label,(offset+n*N)&N2mask);}
                        else unpack64neighbors(q->sw,labelimg,label,(offset+n*N)&N2mask);}
            if (q->se) {if (q->se< 65536) {patt = q->se; unpack16neighbors(patt,labelimg,label,(offset-(offset&Nmask)+((offset+n)&Nmask)+n*N)&N2mask);}
                        else unpack64neighbors(q->se,labelimg,label,(offset-(offset&Nmask)+((offset+n)&Nmask)+n*N)&N2mask);}
        }
    }
    return label;
}

