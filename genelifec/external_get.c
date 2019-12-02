//
//  external_get.c
//  project genelife
//
//  external retrieval of various quantities (typically called from python via genelife_c_interface.py)
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
// get_log2N            get the current log2N value from C to python
// get_curgol           get current gol array from C to python
// get_curgolg          get current golg array from C to python
// get_curgolbr         get current golb and golr arrays from C to python
// get_stats            get the traced statistics from C to python
// get_acttrace         get current acttrace array C to python
// get_poptrace         get current poptrace array C to python
// get genealogytrace   get current trace of genealogies to python
// get_nnovelcells      get N*nhist first counts of number of novel cells (live cells in novel components)
// get_nspecies         get number of species from C to python
// get_nlive            get number of live cells
// get_genealogydepth   get depth of genealogies returned by get_genealogies()
// get_genealogies      get genealogies: not in this file, but in genealogy.c
// get_hist             get the histogram from C to python
// get_activities       get the current activity statistics of genes from C to python
// get_all_activities   get all activity statistics of genes (since t=0) from C to python
// get_quad_activities  get current activity statistics of quads (since t=0) from C to python
// get_small_activities  get current activity statistics of smallpats (since t=0) from C to python
// get_all_quad_activities  get all activity statistics of quads (since t=0) from C to python
// get_all_small_activities  get all activity statistics of smallpats (since t=0) from C to python
// get_components       get all current connected component data structures
// get_smallpatts       get array of small pattern data structures including sizes and activities
// get_quadnodes        get all hashed quadnodes including hashkey, sizes and activities
// get_genes            get all hashed genes with data structures including activity counts, extinctions etc
// get_curgolgstats     get current golgstats array C to python
// get_sorted_popln_act return sorted population and activities (sorted by current population numbers)
// get_gliderinfo       get information about gliders from packed array representation
//------------------------------------------------------------------- get ... ---------------------------------------------------------------------------
// compare with colorfunction==12 code and selectone_of_s code above
int get_glider_count(uint32_t golrstats[], int NN2){
    int ij,rtn=0,d0,inc;
    unsigned int d,dmax,jper;
    int psx,psy;
    extern INLINE void golr_digest (uint64_t golr, unsigned int *mismatchmin, unsigned int *mismatchmax, unsigned int *period, int *pershx, int *pershy);
    
    if (NN2!=N2) {
        fprintf(stderr,"error, wrong number of entries in array, need %d\n",N2);
        return(-1);
    }
    for (ij=0; ij<N2; ij++) {
        if (gol[ij] && (diagnostics & diag_hash_genes)) {
            golr_digest (golr[ij], &d, &dmax, &jper, &psx, &psy); // dynamical record stored in golr
            golrstats[ij]=(d&0xf)|((dmax&0xf)<<4)|((jper&0xf)<<8)|((psx+16)<<16)|((psy+16)<<24); // export 32 bit integer with digested data about golr
            // 0:3 d, 4:7 dmax, 8:15 jper, 16:23 psx+16 24:31 psy+16  --- quality of match d, period = 1+jper, mean displacements ((psx+16)-16)/(16-jper)
            d0 = (psx<0 ? -psx : psx) + (psy<0 ? -psy : psy);
            if (d==dmax || d>4) inc = 0;                          // threshold for periodicity match not reached
            else if (2*d0<=(15-jper))  inc = 0;                   // threshold for mean mobility not reached
            else inc = jper + 1;                                   // glider cell counted with value equal to period
            rtn += inc;
        }
    }
    return(rtn);
}
    
void  get_shist(int outshist[]){
    int ij,idx;
    for(ij=0; ij<N2; ij++){
        idx = golgstats[ij] & F_s_live;
        if(idx>8){
            fprintf(stderr,"bad S value from golstats.\n");
            break;
        }
        outshist[idx]++;
    }
}

int get_log2N() {
    return(log2N);
}
//.......................................................................................................................................................
void get_gliderinfo(uint64_t outgliderinfo[], int narraysize){               // put 7x7 pattern averaged match counts into outgliderinfo array
    uint64_t *gitmp, gene;
    int ij,k,nbhood,sum=0;
    unsigned int d1;
    if(narraysize!=51*8 && narraysize!=27*8 && narraysize!=11*8){
        fprintf(stderr,"get_gliderinfo():  wrong data size (should be array of 8*11,8*27 or 8*51)\n");
    }
    nbhood=narraysize>>3;                                     // divide by eight
    if (nbhood != it_nbhood*it_nbhood+2) {
        fprintf(stderr,"error: mismatch between C value of it_nbhood %d and python expectation %d\n",it_nbhood,nbhood);
        return;
    }
    gitmp=outgliderinfo;
    for (ij=0; ij<N2; ij++) {
        gene = golmix[ij];
        for (k=0;k<8;k++) {                                   // each direction: N E S W NE SE SW NW
            gitmp = outgliderinfo + k*nbhood;
            d1 = (int) ((gene>>(k<<3))&0xffull);              // differences for this direction
            if(d1==0xff)
                gitmp[nbhood-1]++;
            else{
                if(d1<0 || d1>nbhood-2){
                    fprintf(stderr, "get_gliderinfo:  bad difference count value:  %d.", d1);
                    return;
                }
                d1 = nbhood-2-d1;                             // change to matches
                gitmp[d1]++;
            }
        }
    }
    for(d1=0;d1<narraysize;d1++) sum += gitmp[d1];
    // fprintf(stderr,"in get_gliderinfo sum is %d\n",sum);
}
//.......................................................................................................................................................
void get_curgol(uint64_t outgol[], int NN) {
    int ij;
    for (ij=0; ij<NN; ij++) {
        outgol[ij] = planes[curPlane][ij];
    }
}
//.......................................................................................................................................................
void get_curgolg(uint64_t outgolg[], int NN) {
    int ij;
    for (ij=0; ij<NN; ij++) {
        outgolg[ij] = planesg[curPlane][ij];
    }
}
//.......................................................................................................................................................
void get_curgolbr(uint64_t outgolb[], uint64_t outgolr[], int NN) {
    int ij;
    for (ij=0; ij<NN; ij++) {
        outgolb[ij] = planesb[curPlane][ij];
    }
    for (ij=0; ij<NN; ij++) {
        outgolr[ij] = planesr[curPlane][ij];
    }
}
//.......................................................................................................................................................
void get_stats(int outstats[], int outgtypes[], int outstepstats[], int outconfigstats[], int numStats ){
    int i;
    
    if (!(diagnostics & diag_general_statistics)) {
        fprintf(stderr,"statistics collection not enabled in C\n");
        return;
    }
    if(numStats > arraysize){
        fprintf(stderr,"Ack! numStats = %d  > arraysize = %d\n",numStats,arraysize);
        exit(1);
    }
    for(i=0; i<numStats; i++) outstats[i] = livesites[i];
    for(i=0; i<4*numStats; i++) outgtypes[i] = genestats[i];
    for(i=0; i<10*numStats; i++) outstepstats[i] = stepstats[i];
    if (nhistG==nstatG) for(i=0; i<Noff*numStats; i++) outconfigstats[i] = configstats[i];
}
//.......................................................................................................................................................
void get_acttrace(uint64_t outgolg[], int NN) {
    int ij;
    for (ij=0; ij<NN; ij++) {
        outgolg[ij] = acttrace[ij];
    }
}
//.......................................................................................................................................................
void get_acttraceq(uint64_t outgolg[], int NN) {
    int ij;
    for (ij=0; ij<NN; ij++) {
        outgolg[ij] = acttraceq[ij];
    }
}
//.......................................................................................................................................................
void get_poptrace(uint64_t outgolg[], int NN) {
    int ij;
    for (ij=0; ij<NN; ij++) {
        outgolg[ij] = poptrace[ij];
    }
}
//.......................................................................................................................................................
void get_genealogytrace(uint64_t outgolg[], int NN) {
    int ij;
    for (ij=0; ij<NN; ij++) {
        outgolg[ij] = genealogytrace[ij];
    }
}
//.......................................................................................................................................................
void get_nnovelcells(unsigned int outnnovelcells[], int Nh) {
    int i;
    if(Nh>N*nNhist) fprintf(stderr,"error, request for too many entries in get_nnovelcells_trace: %d > %d\n",Nh,N*nNhist);
    else
        for (i=0; i<Nh; i++)  outnnovelcells[i] = nnovelcells[i];
}
//.......................................................................................................................................................
int get_nspecies() {
    int k,nspecies,nspeciesnow;
    nspecies = hashtable_count(&genetable);
    geneitems = (genedata*) hashtable_items( &genetable );

    for (k=0,nspeciesnow=0; k<nspecies; k++)
        if(geneitems[k].popcount) nspeciesnow++;
    
    return(nspeciesnow);
}
//.......................................................................................................................................................
int get_nlive() {
    int ij,nlive;
    for(nlive=0,ij=0;ij<N2;ij++) nlive+= (gol[ij]>0) ? 1 : 0;
    return(nlive);
}
//.......................................................................................................................................................
int get_genealogydepth() {
    int j, jmax, i, nspecies, nspeciesnow;
    uint64_t gene, ancgene;
    int *gindices;
    uint64_t *genes;

    nspecies = hashtable_count(&genetable);
    genotypes = hashtable_keys(&genetable);
    geneitems = (genedata*) hashtable_items( &genetable );

    for (i=nspeciesnow=0; i<nspecies; i++)
        if(geneitems[i].popcount) nspeciesnow++;


    gindices = (int *) malloc(nspecies*sizeof(int));
    for (i=j=0; i<nspecies; i++) {
        if(geneitems[i].popcount) {
            gindices[j++]=i;
        }
        else gindices[nspeciesnow+i-j]=i;
    }

    genes = (uint64_t *) malloc(nspeciesnow*sizeof(uint64_t));

    for (i=0; i<nspeciesnow; i++) {
        genes[i]=genotypes[gindices[i]];
    }
    
    if(ancestortype>0) fprintf(stderr,"Warning: get_genealogydepth currently only implemented for ancestortypes 0 called with %d\n",ancestortype);
    for (i=jmax=0; i<nspeciesnow; i++) {                            // calculate max depth in genealogy jmax
        ancgene=geneitems[gindices[i]].firstancestor;
        for (j=1;;j++) {
            gene=ancgene;
            if(gene==rootgene) break;                               // reached root, exit j loop
            else {
                if((genedataptr = (genedata *) hashtable_find(&genetable, gene)) != NULL) {
                    ancgene=genedataptr->firstancestor;
                }
                else fprintf(stderr,"ancestor not found in genealogies\n");
            }
        }
        if (j>jmax) jmax=j;
    }
    genealogydepth = jmax+1;
    free(genes); free(gindices);
    return(genealogydepth);
}
//.......................................................................................................................................................
int get_curtime(){
    return(totsteps);
}
//.......................................................................................................................................................
void get_histo(int outhisto[],int numHistoC){
    int i;
   
    if (!(diagnostics & diag_offset_statistics)) {
        fprintf(stderr,"histogram of offsets not enabled, activate diag_offset_statistics\n");
        return;
    }
    if(numHistoC != numHisto){
        fprintf(stderr,"Ack! numHisto = %d  != numHistoC = %d\n",numHisto,numHistoC);
        exit(1);
    }
    for(i=0; i<numHisto; i++) outhisto[i] = histo[i];
}
//.......................................................................................................................................................
int get_activities(uint64_t actgenes[], int activities[], int narraysize) {
    int k, nlivegenes, nspecies;

    nspecies = hashtable_count(&genetable);
    genotypes = hashtable_keys(&genetable);
    geneitems = (genedata*) hashtable_items( &genetable );
    // fprintf(stderr,"The number of different species that have ever existed is %d\n",nspecies);

    for (k=nlivegenes=0; k<nspecies; k++) {
        if((genedataptr = (genedata *) hashtable_find(&genetable, genotypes[k])) != NULL) {
            if(genedataptr->popcount) {
                if (nlivegenes <= narraysize) {
                    actgenes[nlivegenes] = genotypes[k];
                    activities[nlivegenes] = genedataptr->activity;
                }
                nlivegenes++;
            }
        }
        else fprintf(stderr,"get_activities error, no entry for gene %"PRIx64" in hash table\n", genotypes[k]);
    }
    if (nlivegenes > narraysize) fprintf(stderr,"Error: array size %d to small to hold live activities %d, increase it\n",narraysize,nlivegenes);

    return nlivegenes;
}
//.......................................................................................................................................................
int get_all_activities(uint64_t genes[], int activities[], int narraysize) {
// get_all_activities   get all activity statistics of genes (since t=0) from C to python
    int k, nspecies;

    nspecies = hashtable_count(&genetable);
    genotypes = hashtable_keys(&genetable);
    geneitems = (genedata *) hashtable_items( &genetable );
    // fprintf(stderr,"The number of different species that have ever existed is %d\n",nspecies);
    if (nspecies > narraysize) {
        fprintf(stderr,"Error: array size %d to small to hold all activities %d, increase it\n",narraysize,nspecies);
        return nspecies;
    }

    for (k=0; k<nspecies; k++) {
        if((genedataptr = (genedata *) hashtable_find(&genetable, genotypes[k])) != NULL) {
            genes[k] = genotypes[k];
            activities[k] = genedataptr->activity;
        }
        else fprintf(stderr,"get_all_activities error, no entry for gene %"PRIx64" in hash table\n", genotypes[k]);
    }
    return nspecies;
}
//.......................................................................................................................................................
int get_quad_activities(uint64_t quads[], int activities[], int narraysize) {
// get_quad_activities  get *live* activity statistics of quads (since t=0) from C to python
    int k, nspecies, livecnt;
    quadnode *q;

    nspecies = hashtable_count(&quadtable);
    quadkeys = hashtable_keys(&quadtable);
    quaditems = (quadnode *) hashtable_items( &genetable );

    for (k=0,livecnt = 0; k<nspecies; k++)
        if((q = (quadnode *) hashtable_find(&quadtable, quadkeys[k])) != NULL){
            if(q->lasttime == totsteps) // test for currently live
                livecnt++;
        } else {
            fprintf(stderr,"get_quad_activities error, no entry for quad %"PRIx64" in hash table\n", quadkeys[k]);
        }
    // fprintf(stderr,"The number of different species that have ever existed is %d\n",nspecies);
    if (livecnt > narraysize) {
        fprintf(stderr,"Error: array size %d to small to hold all quad activities %d, increase it\n",narraysize,livecnt);
        return livecnt;
    }

    for (k=0; k<nspecies; k++) {
        if((q = (quadnode *) hashtable_find(&quadtable, quadkeys[k])) != NULL) {
            if(q->lasttime == totsteps){ // test for currently live
                quads[k] = quadkeys[k];
                activities[k] = q->activity;
            }
        }
        else fprintf(stderr,"get_quad_activities error, no entry for quad %"PRIx64" in hash table\n", quadkeys[k]);
    }
    return nspecies;
}
//.......................................................................................................................................................
int get_all_quad_activities(uint64_t quads[], int activities[], int narraysize) {
// get_all_quad_activities  get all activity statistics of quads (since t=0) from C to python
    int k, nspecies;
    quadnode *q;

    nspecies = hashtable_count(&quadtable);
    quadkeys = hashtable_keys(&quadtable);
    quaditems = (quadnode *) hashtable_items( &genetable );
    // fprintf(stderr,"The number of different species that have ever existed is %d\n",nspecies);
    if (nspecies > narraysize) {
        fprintf(stderr,"Error: array size %d to small to hold all quad activities %d, increase it\n",narraysize,nspecies);
        return nspecies;
    }

    for (k=0; k<nspecies; k++) {
        if((q = (quadnode *) hashtable_find(&quadtable, quadkeys[k])) != NULL) {
            quads[k] = quadkeys[k];
            activities[k] = q->activity;
        }
        else fprintf(stderr,"get_quad_activities error, no entry for quad %"PRIx64" in hash table\n", quadkeys[k]);
    }
    return nspecies;
}
//.......................................................................................................................................................
int get_small_activities(uint64_t smalls[], int activities[], int narraysize) {
// get_small_activities  get *live only* activity statistics of smallpatts (since t=0) from C to python
    int k, nspecies;

    if (narraysize<65536) {
        fprintf(stderr,"Error in get_small_activities : called with insufficent smallpatt holding array size %d < %d\n",narraysize,65536);
        return -1;
    }
    nspecies =0;
    for (k=0; k<65536; k++) {
        if(smallpatts[k].lasttime == totsteps){
            nspecies += smallpatts[k].activity ? 1 : 0;
            activities[k] = smallpatts[k].activity;
        }
    }
    return nspecies;
}
//.......................................................................................................................................................
int get_all_small_activities(uint64_t smalls[], int activities[], int narraysize) {
// get_all_small_activities  get all activity statistics of quads (since t=0) from C to python
    int k, nspecies;

    if (narraysize<65536) {
        fprintf(stderr,"Error in get_small_activities : called with insufficent smallpatt holding array size %d < %d\n",narraysize,65536);
        return -1;
    }
    nspecies =0;
    for (k=0; k<65536; k++) {
        nspecies += smallpatts[k].activity ? 1 : 0;
        activities[k] = smallpatts[k].activity;
    }
    return nspecies;
}
//.......................................................................................................................................................
int get_connected_comps(unsigned int outlabel[], unsigned int outconnlen[], int x, int y) {
    int i,ij;
    for (ij=0; ij<N2; ij++) {
        outlabel[ij] = (unsigned int) label[ij];
    }
    for (i=1;i<ncomponents+1;i++) {
        outconnlen[i] = (unsigned int) connlen[i];
    }
    xdisplay = x;
    ydisplay = y;
    return ncomponents;
}
//.......................................................................................................................................................
int get_ncomponents() {
    return(ncomponents);
}
//.......................................................................................................................................................
int get_components(component components[],int narraysize) {
    int i;
    if (narraysize<ncomponents) {
        fprintf(stderr,"Error in get_components : called with insufficent component holding array size %d < %d\n",narraysize,ncomponents);
        return -1;
    }
    for (i=1;i<=ncomponents;i++) {
        components[i-1]=complist[i];
    }

    return ncomponents;
}
//.......................................................................................................................................................
int get_smallpatts(smallpatt smallpattsout[],int narraysize) {
    int i,count;
    if (narraysize<65536) {
        fprintf(stderr,"Error in get_smallpatts : called with insufficent smallpatt holding array size %d < %d\n",narraysize,65536);
        return -1;
    }
    for (count=i=0;i<65536;i++) {
        // smallpattsout[i].size= smallpatts[i].size;
        smallpattsout[i].topactivity= smallpatts[i].topactivity;
        smallpattsout[i].activity= smallpatts[i].activity;
        count += smallpatts[i].activity ? 1 : 0;
        smallpattsout[i].firsttime= smallpatts[i].firsttime;
        smallpattsout[i].lasttime= smallpatts[i].lasttime;
    }
    return count;
}
//.......................................................................................................................................................
int get_quadnodes(quadnode quadnodes[],int narraysize) {
    int i;
    // these three calls executed through hashactivityquad if colorfunction 9 or 10
    // nallspeciesquad = hashtable_count(&quadtable);
    // quadkeys = hashtable_keys(&quadtable);
    // quaditems = (quadnode*) hashtable_items( &quadtable );

    if (narraysize<nallspeciesquad) {
        fprintf(stderr,"Error in get_quadnodes : called with insufficent quadnode holding array size %d < %d\n",narraysize,nallspeciesquad);
        return -1;
    }
    for (i=0;i<nallspeciesquad;i++) {
        quadnodes[i]=quaditems[i];
    }
    return nallspeciesquad;
}
//.......................................................................................................................................................
int get_genes(genedata genelist[],int narraysize) {
    int i;
    // these three calls executed already through hashactivity
    // nallspecies = hashtable_count(&genetable);
    // genotypes = hashtable_keys(&genetable);
    // geneitems = (genedata*) hashtable_items( &genetable );

    if (narraysize<nallspecies) {
        fprintf(stderr,"Error in get_genes : called with insufficent genedata holding array size %d < %d\n",narraysize,nallspecies);
        return -1;
    }
    for (i=0;i<nallspecies;i++) {
        genelist[i]=geneitems[i];          /* shallow copy OK if no pointers being copied, otherwise the structures they point to will not be copied */
    }
    return nallspecies;
}
//.......................................................................................................................................................
void get_curgolgstats(uint64_t outgolgstats[], int NN) {
    int ij;
    for (ij=0; ij<NN; ij++) {
        outgolgstats[ij] = golgstats[ij];                       // Note that golgstats is not dealt with in planes !
    }
}
//.......................................................................................................................................................
int get_sorted_popln_act( int gindices[], uint64_t genes[], int popln[], int activities[]) {
    int nspecies;
    int activitieshashx(int gindices[], uint64_t genes[], int popln[], int activities[]);
    nspecies=activitieshashx(gindices, genes, popln, activities);        // sets acttrace and returns current population arrays
    return(nspecies);
}
