//
//  clonealogy.c
//  project genelife
//
//  rapid construction of clonealogies for genelife using stored hash tables, callable at each time step
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
// clonealogies         calculate and display clonealogies: genealogies of clones
//--------------------------------------------------------------- clonealogies --------------------------------------------------------------------------
int clonealogies() {                                            // genealogies of all clones
    int j, jmax, i, ij, nclones, nclonesnow, birthstep;
    int j1, j2, j3, activity, gorder[N];
    uint64_t birthid, parentid, nextclone;
    // uint64_t gene, ancgene, nextgene;
    int *gindices,*popln,*activities,*birthsteps;
    extern int cmpfunc3c (const void * pa, const void * pb);
    extern int cmpfunc5c (const void * pa, const void * pb);
    

    nclones = hashtable_count(&clonetable);
    clones = hashtable_keys(&clonetable);
    cloneitems = (clonedata*) hashtable_items( &clonetable );

    for (i=nclonesnow=0; i<nclones; i++)
        if(cloneitems[i].popln) nclonesnow++;

    gindices = (int *) malloc(nclones*sizeof(int));
    for (i=j=0; i<nclones; i++) {
        if(cloneitems[i].popln) {
            gindices[j]=i;
            j++;
        }
        else gindices[nclonesnow+i-j]=i;
    }

    qsort(gindices, nclonesnow, sizeof(int), cmpfunc3c);        // sort in decreasing population size order

    if (nclonesnow>N) nclonesnow=N;                             // can only display at most N clones, chose oldest

    popln = (int *) malloc(nclonesnow*sizeof(int));
    activities = (int *) malloc(nclonesnow*sizeof(int));
    birthsteps = (int *) malloc(nclonesnow*sizeof(int));

    for (i=0; i<nclonesnow; i++) {
        popln[i]=cloneitems[gindices[i]].popln;
        activities[i]=cloneitems[gindices[i]].activity;
        birthsteps[i]=(unsigned int) ((cloneitems[gindices[i]].birthid)>>32);
    }

    for(ij=0;ij<N2;ij++) working[ij]=rootclone;                 // set field to rootclone as background
    activitymax=0;
    for (i=jmax=0; i<nclonesnow; i++) {
        birthid=clones[gindices[i]];                            // do not need to copy array since only needed here
        parentid=cloneitems[gindices[i]].parentid;
        activity=cloneitems[gindices[i]].activity;
        if(activity>activitymax) activitymax=activity;
        working[i]=birthid;                                     // ij = i for j=0
        for (j=1;j<N;j++) {                                     // go back at most N links in clonealogy
            birthid=parentid;
            if(birthid&rootclone) break;                       // reached root, exit j loop
            else {
                if((clonedataptr = (clonedata *) hashtable_find(&clonetable, birthid)) != NULL) {
                    parentid=clonedataptr->parentid;
                    activity = clonedataptr->activity;
                    if(activity>activitymax) activitymax=activity;
                }
                else fprintf(stderr,"ancestor not found in clonealogies\n");
            }
            ij = i+j*N;
            working[ij]=birthid;
        }
        if (j>jmax) jmax=j;
    }
    clonealogydepth = jmax;

                                                                //reverse ancestries to allow comparison at same number of clonal speciations
    for (i=0; i<nclonesnow; i++) {
        for(j=0;j<N;j++) {
            if (working[i+j*N]&rootclone) break;
        }
        for(j1=0;j1<(j>>1);j1++) {
            birthid=working[i+(j-j1-1)*N];
            working[i+(j-j1-1)*N]=working[i+j1*N];
            working[i+j1*N]=birthid;
        }
    }
    for (i=0; i<N; i++) gorder[i]=i;
    qsort(gorder, nclonesnow, sizeof(int), cmpfunc5c);          // sort according to ancestral lines - use cmpfunc5 to sorting clones laterally via clone value

    for (i=0;i<N;i++) if((gorder[i]<0)||(gorder[i]>=N)) fprintf(stderr,"step %d error in gorder out of bounds at i = %d with value %d\n",totsteps,i,gorder[i]);

    for(ij=0;ij<N2;ij++) clonealogytrace[ij]=rootclone;          // initialize clonealogytrace to root clone before drawing part of it

    if(colorfunction==7 || colorfunction2==7) {                 // time trace of clonealogies (takes precedence with two difft displays
      birthstep=0;
      for(i=0;i<nclonesnow;i++) {
        for(j=0,j1=0;j<jmax;j++) {
            if(gorder[i]>=nclonesnow) fprintf(stderr,"error in clonealogies gorder at i=%d, order value %d out of range\n",i,gorder[i]);
            ij = gorder[i]+j*N;
            birthid = working[ij];
            ij+=N;
            if(ij<N2) {
                nextclone = working[ij];
                if(nextclone&rootclone) birthstep=totsteps;
                else {
                    if((clonedataptr = (clonedata *) hashtable_find(&clonetable, nextclone)) != NULL) birthstep = (unsigned int) (clonedataptr->birthid>>32);
                    else fprintf(stderr,"ancestor %"PRIx64" not found at (%d,%d) in clonealogies during birthstep extraction\n",nextclone,ij&Nmask,ij>>log2N);
                }
            }
            else birthstep=totsteps;
            j2 = birthstep*N/totsteps;
            for (j3=j1;j3<j2;j3++) {
                ij = i+j3*N;
                clonealogytrace[ij]=birthid;
            }
            j1 = j2;
        }
      }
      // for(i=nclonesnow;i<N;i++) for(j=0;j<N;j++) clonealogytrace[gorder[i]+j*N]=rootclone;
    }
    else {                                                      // species changes only trace (colorfunction == 6)
      for(i=0;i<nclonesnow;i++) {
        for(j=0;j<jmax;j++) {
            ij=i+j*N;
            clonealogytrace[ij]=working[gorder[i]+j*N];
        }
      }
    }
    free(gindices);free(activities);free(popln);free(birthsteps);
    return(jmax);
}
