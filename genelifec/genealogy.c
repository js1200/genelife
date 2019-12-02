//
//  genealogy.c
//  project genelife
//
//  rapid construction of genealogies for genelife using stored hash tables, callable at each time step
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
// get_genealogies      calculate and retrieve or display genealogies, depending on size of array passed (0: display only, >0: retrieve only)
//--------------------------------------------------------------- genealogies ---------------------------------------------------------------------------
int get_genealogies(genedata genealogydat[], int narraysize) {  /* genealogies of all currently active species */
    int j, jmax, i, ij, k, nspecies, nspeciesnow;
    unsigned int birthstep;
    int j1, j2, j3, activity, gorder[N];
    uint64_t gene, ancgene, nextgene, genealogy1[N];
    uint64_t *curgen,*curgenealogy;
    int *gindices,*popln,*activities;
    genedata genedummy = {0,0,0,0,0,0,0,0ull,rootgene,rootgene};  // default data structure for gene data: note I put 4th member lastextinction -1 which is not a valid timestep
    extern int cmpfunc3 (const void * pa, const void * pb);
    extern int cmpfunc5 (const void * pa, const void * pb);

    nspecies = hashtable_count(&genetable);
    genotypes = hashtable_keys(&genetable);
    geneitems = (genedata*) hashtable_items( &genetable );

    for (i=nspeciesnow=0; i<nspecies; i++)
        if(geneitems[i].popcount) nspeciesnow++;

    gindices = (int *) malloc(nspecies*sizeof(int));
    for (i=j=0; i<nspecies; i++) {
        if(geneitems[i].popcount) {
            gindices[j]=i;
            j++;
        }
        else gindices[nspeciesnow+i-j]=i;
    }

    if(!narraysize) {
        qsort(gindices, nspeciesnow, sizeof(int), cmpfunc3);// sort in decreasing population size order
        // qsort(gindices, nspeciesnow, sizeof(int), cmpfunc4);// alternatively, sort in increasing birthstep order
        if (nspeciesnow>N) nspeciesnow=N;                           // can only display at most N species, chose oldest
        curgenealogy=genealogy1;
    }
    else curgenealogy=working;

    popln = (int *) malloc(nspeciesnow*sizeof(int));
    activities = (int *) malloc(nspeciesnow*sizeof(int));

    for (i=0; i<nspeciesnow; i++) {
        popln[i]=geneitems[gindices[i]].popcount;
        activities[i]=geneitems[gindices[i]].activity;
    }

    if(narraysize) {                                                // need to allocate data for genealogy array in python
        for (i=jmax=0; i<nspeciesnow; i++) {                        // calculate max depth in genealogy jmax
            gene=genotypes[gindices[i]];
            ancgene=geneitems[gindices[i]].firstancestor;
            working[0]=gene;                                        // use working instead of genealogy1 array here as length needed may be larger than N
            for (j=1;;j++) {
                gene=ancgene;
                if(gene==rootgene) break;                           // reached root, exit j loop
                else {
                    for (k=0;k<j;k++) if (gene==working[k]) {gene=generepeat;break;};   // if gene already in ancestry, break with generepeat
                    if(gene==generepeat) { j=j+1;break;}
                    if((genedataptr = (genedata *) hashtable_find(&genetable, gene)) != NULL) {
                        ancgene=genedataptr->firstancestor;
                    }
                    else fprintf(stderr,"ancestor not found in genealogies\n");
                }
                working[j]=gene;
            }
            if (j>jmax) jmax=j;
        }
        
        if(narraysize < (jmax+1)*nspeciesnow){
            fprintf(stderr,"get_genealogies(): narraysize not large enough.  Must be at least %d\n",nspeciesnow*(jmax+1));
            free(gindices);free(activities);free(popln);
            return(-1);
        }
        curgen = (uint64_t *) calloc(jmax,sizeof(uint64_t)); // current genealogy array
    }
    else {                                                        // display version
        for(ij=0;ij<N2;ij++) working[ij]=rootgene;                // set field to rootgene as background
        jmax=N-1;
        curgen = NULL;
    }
    
    activitymax=0;
    for (i=ij=genealogydepth=0; i<nspeciesnow; i++) {
        gene=genotypes[gindices[i]];                              // do not need to copy array to genes since only needed here
        if(narraysize) curgen[0]=gene;
        ancgene=geneitems[gindices[i]].firstancestor;
        activity=geneitems[gindices[i]].activity;
        if(activity>activitymax) activitymax=activity;
        curgenealogy[0]=gene;
        if(!narraysize) working[i]=gene;                          // ij = i for j=0
        for (j=k=1;j<=jmax;j++) {                                 // go back at most jmax links in genealogy
            gene=ancgene;
            if(gene==rootgene) break;                             // reached root, exit j loop
            else {
                for (k=0;k<j;k++) if (gene==curgenealogy[k]) {gene=generepeat;break;};  // if gene already in ancestry, break with generepeat
                if(gene==generepeat) { if(narraysize) curgen[j]=gene;else working[i+j*N]=gene;j=j+1;break;}
                if((genedataptr = (genedata *) hashtable_find(&genetable, gene)) != NULL) {
                    ancgene=genedataptr->firstancestor;
                    activity = genedataptr->activity;
                    if(activity>activitymax) activitymax=activity;
                }
                else fprintf(stderr,"ancestor not found in genealogies\n");
            }
            curgenealogy[j]=gene;
            if(narraysize) curgen[j]=gene; else working[i+j*N]=gene;
        }
        if (j>genealogydepth) genealogydepth=j;
        if(narraysize) {                                        // copy genealogy to python
            for(k=0; k<j; k++){
                if((genedataptr = (genedata *) hashtable_find(&genetable, curgen[k])) != NULL) {
                    genealogydat[ij+k] = *genedataptr; // cf comment on shallow copy in get_genes()
                }
                else fprintf(stderr,"ancestor not found in genealogies\n");
            }
            for(k=j; k<=jmax; k++){                             // fill list up to uniform jmax with dummy data
                genealogydat[ij+k] = genedummy;
            }
            ij += (jmax+1);
        }
    }
    if(narraysize) {                                            // if narraysize then return from routine here
        free(curgen);
        genealogydepth = jmax+1;
        free(gindices);free(activities);free(popln);
        return(genealogydepth);
    }
                                                                // reverse ancestries to allow comparison at same number of speciations
    for (i=0; i<nspeciesnow; i++) {
        for(j=0;j<N;j++) {
            if (working[i+j*N]==rootgene ) break;  // || working[i+j*N]==generepeat
        }
        for(j1=0;j1<(j>>1);j1++) {
            gene=working[i+(j-j1-1)*N];
            working[i+(j-j1-1)*N]=working[i+j1*N];
            working[i+j1*N]=gene;
        }
    }
    for (i=0; i<N; i++) gorder[i]=i;
    qsort(gorder, nspeciesnow, sizeof(int), cmpfunc5);          // sort according to ancestral lines - use cmpfunc5 to sorting genes laterally via gene value
    //qsort(gorder, nspeciesnow, sizeof(int), cmpfunc6);        // sort according to ancestral lines - use cmpfunc6 to sorting genes laterally via activity
    //qsort(gorder, nspeciesnow, sizeof(int), cmpfunc7);        // sort according to ancestral lines - use cmpfunc7 to sort genes laterally via population size

    for (i=0;i<N;i++) if((gorder[i]<0)||(gorder[i]>=N)) fprintf(stderr,"step %d error in gorder out of bounds at i = %d with value %d\n",totsteps,i,gorder[i]);

    for(ij=0;ij<N2;ij++) genealogytrace[ij]=rootgene;           // initialize genealogytrace to root gene before drawing part of it

    if(colorfunction==7 || colorfunction2==7) {                 // time trace of genealogies: NB simultaneus display of 6 and 7 not possible :> 2x 7 or 2x 6
      birthstep=0;
      for(i=0;i<nspeciesnow;i++) {
        for(j=0,j1=0;j<jmax;j++) {
            if(gorder[i]>=nspeciesnow) fprintf(stderr,"error in genealogies gorder at i=%d, order value %d out of range\n",i,gorder[i]);
            ij = gorder[i]+j*N;
            gene = working[ij];
            ij+=N;
            if(ij<N2) {
                nextgene = working[ij];
                if(nextgene==rootgene) birthstep=totsteps;
                else {
                    if((genedataptr = (genedata *) hashtable_find(&genetable, nextgene)) != NULL) {
                        birthstep = (unsigned int) genedataptr->firsttime;
                    }
                    else fprintf(stderr,"ancestor %"PRIx64" not found at (%d,%d) in genealogies during birthstep extraction\n",nextgene,ij&Nmask,ij>>log2N);
                }
            }
            else birthstep=totsteps;
            j2 = birthstep*N/totsteps;
            for (j3=j1;j3<j2;j3++) {
                ij = i+j3*N;
                genealogytrace[ij]=gene;
            }
            j1 = j2;
        }
      }
      // for(i=nspeciesnow;i<N;i++) for(j=0;j<N;j++) genealogytrace[gorder[i]+j*N]=rootgene;
    }
    else {                                                      // species changes only trace (colorfunction == 6)
      for(i=0;i<nspeciesnow;i++) {
        for(j=0;j<jmax;j++) {
            ij=i+j*N;
            genealogytrace[ij]=working[gorder[i]+j*N];
        }
      }
    }
    free(gindices);free(activities);free(popln);
    genealogydepth = jmax+1;

    return(genealogydepth);
}
