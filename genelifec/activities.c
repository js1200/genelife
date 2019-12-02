//
//  activities.c
//  project genelife
//
//  calculation of activites for genes, clones and patterns
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
// activitieshash       calculate array of current gene activities and update acttrace array of genes in activity plot format
// activitieshashquad   calculate array of current quad activities and update acttraceq array of patterns in activity plot format
//------------------------------------------------------------ activitieshash ---------------------------------------------------------------------------
extern INLINE unsigned int genefnindex( uint64_t gene, uint64_t mask, int indexoff[]) {
    int k,d;
    uint64_t g,gf;
    uint64_t indexmask = 0xffffff;
    POPCOUNT64C(mask, d)
    if (d > 24) d=24;                                           // max of first 24 of active lut sites in mask allowed into genefnactivities array
    g = gene&mask;
    gf = 0ull;
    for (k=d;k>=0;k--) {
        gf |= (g>>indexoff[k])&0x1ull;
        gf <<= 1;
    }
    gf &= indexmask;
    return((unsigned int) gf);
}
//.......................................................................................................................................................
int activitieshash() {  /* count activities of all currently active gene species */
    int i, j, jmax, k, ij, ij1, x, nchist, nrhist, nspecies, nspeciesnow, cnt0, cnt1;
    int *gindices,*popln,*activities;
    double act,pop;
    uint64_t *genes, *traceptr;
    uint64_t gene,sbmask;
    const int maxact = 10000;
    int indexoff[24];
    extern int cmpfunc3 (const void * pa, const void * pb);
    
    sbmask = 0ull;
    if (activityfnlut) {
        for (j=0;j<24;j++) indexoff[j] = 0;
        sbmask = (((uint64_t) birthmask) << 32) | (uint64_t) survivalmask;
        for (j=k=0;k<64;k++) {
            if (j>=24) break;
            if ((sbmask>>k)&0x1ull) indexoff[j++]=k;
        }
        // fprintf(stderr,"indexoff:");for (j=0;j<24;j++) fprintf(stderr," %4d",indexoff[j]);fprintf(stderr,"\n");
    }

    nspecies = hashtable_count(&genetable);
    genotypes = hashtable_keys(&genetable);
    geneitems = (genedata *) hashtable_items( &genetable );

    // if(gindices != NULL && col) free(gindices);
    for (i=0,nspeciesnow=0; i<nspecies; i++)
        nspeciesnow+=geneitems[i].popcount ? 1 : 0;

    gindices = (int *) malloc(nspeciesnow*sizeof(int));

    for (i=j=0; i<nspecies; i++) {
        if(geneitems[i].popcount) {
            gindices[j++]=i;                                           // if col is 0 then the array gindices must be passed with sufficient length
        }
    }

    if(activityfnlut) {
        genes = (uint64_t *) malloc(nspeciesnow*sizeof(uint64_t));
        popln = (int *) malloc(nspeciesnow*sizeof(int));
        activities = (int *) malloc(nspeciesnow*sizeof(int));
        memset(popln,0,sizeof(int)*nspeciesnow);
        memset(activities,0,sizeof(int)*nspeciesnow);
        memset(genefnindices,0,sizeof(int)*(1<<24));
        for (i=0,jmax=1; i<nspeciesnow; i++) {
            k=genefnindex( genotypes[gindices[i]], sbmask, indexoff);
            j=genefnindices[k];
            if(!j) j=genefnindices[k]=jmax++;
            // genes[j-1]=genotypes[gindices[i]]&sbmask;
            genes[j-1]=genotypes[gindices[i]];
            popln[j-1]+=geneitems[gindices[i]].popcount;
            activities[j-1]+=geneitems[gindices[i]].activity;
        }
        nspeciesnow = jmax-1;
    }
    else {
        if (nspeciesnow > maxact) {                                      //sort in order of decreasing population
            qsort(gindices, nspeciesnow, sizeof(int), cmpfunc3);         // sort in decreasing count order
            nspeciesnow = maxact;
        }
        genes = (uint64_t *) malloc(nspeciesnow*sizeof(uint64_t));
        popln = (int *) malloc(nspeciesnow*sizeof(int));
        activities = (int *) malloc(nspeciesnow*sizeof(int));
        for (i=0; i<nspeciesnow; i++) {
            genes[i]=genotypes[gindices[i]];
            popln[i]=geneitems[gindices[i]].popcount;
            activities[i]=geneitems[gindices[i]].activity;
        }
    }

    if (totdisp>=N) {                                               // 1 pixel to left scroll when full
        for(ij=0;ij<N2;ij++) {
            ij1 = ((ij+1)&Nmask)+((ij>>log2N)<<log2N);              // (i+1)%N+j*N;
            // if(ij1>=N2) fprintf(stderr,"error in scroll of acttrace\n");
            acttrace[ij]=acttrace[ij1];
            poptrace[ij]=poptrace[ij1];
        }
        x=N-1;
    }
    else x=totdisp;

    for(i=0;i<N;i++) acttrace[x+i*N]=rootgene;                      // set column gray
    for(i=0;i<N;i++) poptrace[x+i*N]=rootgene;                      // set column gray
    //for(i=ymax1=0;i<nspeciesnow;i++)                              // ymax1 is current maximum of activities
    //    ymax1 = (activities[i]>ymax1) ? activities[i] : ymax1;
    // if (ymax1>ymax) ymax = ymax*2;                               // autoscale of activities
    // if (ymax1<ymax/2) ymax = ymax/2;                             // autoscale of activities
    for(j=0;j<nspeciesnow;j++) {
        gene = genes[j];
                                                                            // rescale populations and activities with saturation
        act = (double) activities[j];
        pop = (double) popln[j];
        // activities[j] = N-1 - (activities[j] * (N-1)) / ymax;    // linear scale, needs truncation if ymax superceded
        // activities[j] = (N-1) - (int) ((N-1)*log2(act)/log2ymax);// logarithmic scale, suffers from discrete steps at bottom
        activities[j] = (N-1) - (int) ((N-1)*act/(act+(double)ymax));
        popln[j] = (N-1) - (int) ((N-1)*pop/(pop+(double)ymax/10.));
        
        ij = (x&Nmask)+activities[j]*N;
        if(acttrace[ij]==rootgene)                                  // only one genotype to plot
            acttrace[ij] = gene;
        else {                                                      // plot species color with largest current population size if choice of multiple
            if(activityfnlut) {
                if((k=genefnindices[genefnindex( acttrace[ij], sbmask, indexoff)])) cnt0 = popln[k-1];
                else cnt0 = 0;
                cnt1 = popln[j];
                if(cnt1 >= cnt0) acttrace[ij]=gene;
            }
            else {
                if((genedataptr = (genedata *) hashtable_find(&genetable, acttrace[ij])) != NULL) cnt0 = genedataptr->popcount;
                else cnt0 = 0;
                if((genedataptr = (genedata *) hashtable_find(&genetable, gene)) != NULL) cnt1 = genedataptr->popcount;
                else cnt1 = 0;
                if(cnt1 >= cnt0) acttrace[ij]=gene;
            }
        }
        
        ij = (x&Nmask)+popln[j]*N;
        if(poptrace[ij]==rootgene)                                  // only one genotype to plot
            poptrace[ij] = gene;
        else {                                                      // plot species color with largest current activity if choice of multiple
            if(activityfnlut) {
                if((k=genefnindices[genefnindex( acttrace[ij], sbmask, indexoff)])) cnt0 = activities[k-1];
                else cnt0 = 0;
                cnt1 = activities[j];
                if(cnt1 >= cnt0) acttrace[ij]=gene;
            }
            else {
                if((genedataptr = (genedata *) hashtable_find(&genetable, poptrace[ij])) != NULL)
                   cnt0 = genedataptr->activity;
                else cnt0 = 0;
                if((genedataptr = (genedata *) hashtable_find(&genetable, gene)) != NULL)
                   cnt1 = genedataptr->activity;
                else cnt1 = 0;
                if(cnt1 >= cnt0) poptrace[ij]=gene;
            }
        }
    }
    if(totdisp<N*nNhist) {
        nchist=totdisp/N; nrhist=totdisp-nchist*N;
        traceptr=&acttrace1[N2*nchist];
        for(j=0;j<N;j++) traceptr[nrhist+j*N]=acttrace[x+j*N];
        traceptr=&poptrace1[N2*nchist];
        for(j=0;j<N;j++) traceptr[nrhist+j*N]=poptrace[x+j*N];
        npopulation1[totdisp]=npopulation[x];
        nnovelcells[totdisp]=novelcells();
    }
    free(gindices);free(activities);free(genes);free(popln);
    return(nspeciesnow);
}
//.......................................................................................................................................................
int activitieshashx(int gindices[], uint64_t genes[], int popln[], int activities[]) {  /* python interface to count activities of all currently active species, no display */
    int i, j, nspecies, nspeciesnow;
    const int maxact = 10000;
    extern int cmpfunc3 (const void * pa, const void * pb);

    nspecies = hashtable_count(&genetable);
    genotypes = hashtable_keys(&genetable);
    geneitems = (genedata*) hashtable_items( &genetable );

    // if(gindices != NULL && col) free(gindices);
    for (i=0,nspeciesnow=0; i<nspecies; i++)
        nspeciesnow+= geneitems[i].popcount ? 1 : 0;

    if (nspeciesnow>10000) return(-1);                              // exit with error need to allocate more space in python
    for (i=j=0; i<nspecies; i++) {
        gindices[j] = (geneitems[i].popcount) ? i : gindices[j];      // if col is 0 then the array gindices must be passed with sufficient length
        j += (geneitems[i].popcount) ? 1 : 0;
    }
    qsort(gindices, nspeciesnow, sizeof(int), cmpfunc3);            // sort in decreasing count order

    if (nspeciesnow > maxact) nspeciesnow = maxact;

    for (i=0; i<nspeciesnow; i++) {
        genes[i]=genotypes[gindices[i]];
        popln[i]=geneitems[gindices[i]].popcount;
        activities[i]=geneitems[gindices[i]].activity;
    }

    return(nspeciesnow);                                            // exit here without doing display

}
//.......................................................................................................................................................
int activitieshashquad() {  /* count activities of all currently active quad images of connected components */
    int i, j, ij, ij1, x, nspecies, nspeciesnow, popcnt0, popcnt1;
    int *qindices,*popln,*activities;
    int qsindices[65536]; // ,qsallindices[65536];
    quadnode *q;
    double act;
    uint64_t *qids;                                                   // 64 bit ids for components used for colouring and ID
    uint64_t qid;
    const int maxact = 10000;
    extern int cmpfunc3q (const void * pa, const void * pb);
    extern int cmpfunc3qs (const void * pa, const void * pb);

    nspecies = hashtable_count(&quadtable);
    quadkeys = hashtable_keys(&quadtable);
    quaditems = (quadnode*) hashtable_items( &quadtable );


    if(nhistG && (totsteps%nhistG == 0)) {                            // collect cumulative pattern size histograms of entire hash table
        for(i=0;i<=log2N;i++) histcumlogpattsize[i]=0;
        for(i=0;i<=N;i++) histcumpixelssqrt[i]=0;
        for (i=0; i<nspecies; i++) {
            histcumlogpattsize[log2upper(quaditems[i].size)]++;
            histcumpixelssqrt[sqrtupper(quaditems[i].pop1s)]++;
        }
    }
    
    if (totdisp>=N) {                                                 // 1 pixel to left scroll when full
        for(ij=0;ij<N2;ij++) {
            ij1 = ((ij+1)&Nmask)+((ij>>log2N)<<log2N);                // (i+1)%N+j*N;
            // if(ij1>=N2) fprintf(stderr,"error in scroll of acttraceq\n");
            acttraceq[ij]=acttraceq[ij1];
            acttraceqt[ij]=acttraceqt[ij1];
        }
        x=N-1;
    }
    else x=totdisp;
    for(i=0;i<N;i++) acttraceq[x+i*N]=rootgene;                       // set column gray, rootgene is used as unique pattern mapped to gray as for gene activities
    
    for (i=0,nspeciesnow=0; i<nspecies; i++)
        nspeciesnow+=quaditems[i].lasttime==totsteps ? 1 : 0;

    if (nspeciesnow) {
        qindices = (int *) malloc(nspeciesnow*sizeof(int));

        for (i=j=0; i<nspecies; i++) {
            if(quaditems[i].lasttime==totsteps) qindices[j++]=i;      // if col is 0 then the array qindices must be passed with sufficient length
        }
        if (nspeciesnow > maxact) {                                   // sort in order of decreasing pixel count
            qsort(qindices, nspeciesnow, sizeof(int), cmpfunc3q);
        }
        if (nspeciesnow > maxact) nspeciesnow = maxact;

        qids = (uint64_t *) malloc(nspeciesnow*sizeof(uint64_t));     // allocate arrays
        popln = (int *) malloc(nspeciesnow*sizeof(int));
        activities = (int *) malloc(nspeciesnow*sizeof(int));

        for (i=0; i<nspeciesnow; i++) {                               // set arrays of ids, popln (nr 1 pixels), and activities from hash table
            qids[i]=quadkeys[qindices[i]];
            popln[i]=quaditems[qindices[i]].pop1s;
            activities[i]=quaditems[qindices[i]].topactivity;
        }

        for(j=0;j<nspeciesnow;j++) {                                 // main loop to construct new display column for activities
            act = (double) activities[j];
            activities[j] = (N-1) - (int) ((N-1)*act/(act+(double)ymaxq));
            qid = qids[j];
            ij = (x&Nmask)+activities[j]*N;
            if(acttraceq[ij]==rootgene) {                            // first quadtype to plot at this position
                acttraceq[ij] = qid;
                acttraceqt[ij] = 1;                                  // tye of entry is quadtree (not small pattern)
            }
            else {                                                   // plot species color with largest current pop1s size if choice of multiple
                if((q = (quadnode *) hashtable_find(&quadtable, acttraceq[ij])) != NULL)
                   popcnt0 = q->pop1s;
                else popcnt0 = 0;
                if((q = (quadnode *) hashtable_find(&quadtable, qid)) != NULL)
                   popcnt1 = q->pop1s;
                else popcnt1 = 0;
                if(popcnt1 >= popcnt0) {
                    acttraceq[ij]=qid;
                    acttraceqt[ij] = 1;
                }
            }
        }
        free(qindices);free(activities);free(qids);free(popln);
    }
    
    if (nspeciesnow<maxact) {                                        // overlay activities of smallpatts up to maxact
        for (nallspeciessmall=nspeciessmall=i=0;i<65536; i++) {
            if (smallpatts[i].topactivity) {
                // qsallindices[nallspeciessmall++]=i;
                if (smallpatts[i].lasttime == totsteps) qsindices[nspeciessmall++]=i;   // indices of current patterns
            }
        }
        if (nspeciesnow+nspeciessmall > maxact) {                    //sort in order of decreasing pixel count
            qsort(qsindices, nspeciessmall, sizeof(int), cmpfunc3qs);
            nspeciessmall = maxact-nspeciesnow;
        }
        
        activities = (int *) malloc(nspeciessmall*sizeof(int));
        for (i=0; i<nspeciessmall; i++) {
            activities[i]=smallpatts[qsindices[i]].topactivity;
        }
        
        for(j=0;j<nspeciessmall;j++) {
            act = (double) activities[j];
            activities[j] = (N-1) - (int) ((N-1)*act/(act+(double)ymaxq));
            qid = qsindices[j];
            ij = (x&Nmask)+activities[j]*N;
            if(acttraceq[ij]==rootgene) {                          // first quadtype to plot at this position (prefer quad over smallpatts)
                acttraceq[ij] = qid;
                acttraceqt[ij] = 0;
            }
            else if (!acttraceqt[ij]) {                            // plot species color with largest current pop1s size if choice of multiple
                POPCOUNT64C(((uint64_t) acttraceq[ij]), popcnt0);
                POPCOUNT64C(((uint64_t) qid), popcnt1);
                if(popcnt1 >= popcnt0) {
                    acttraceq[ij]=qid;
                    acttraceqt[ij] = 0;
                }
            }
        }
        free(activities);
        if(totdisp<N*nNhist) {
            for(i=0;i<N;i++) acttraceq1[totdisp+i*N]=acttraceq[x+i*N];
            for(i=0;i<N;i++) acttraceqt1[totdisp+i*N]=acttraceqt[x+i*N];
        }
    }
    
    return(nspeciesnow);
}
