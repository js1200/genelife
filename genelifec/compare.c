//
//  compare.c
//  project genelife
//
//  custom comparison functions, primarily for use in sorting
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
// cmpfunc              compare gene values as numerical unsigned numbers
// cmpfunc1             compare gene counts in population
// cmpfunc2             compare gene values corresponding to given number index in hash table
// cmpfunc3             compare population counts of hash stored genes
// cmpfunc3c            compare population counts of hash stored clones
// cmpfunc3q            compare pixel counts (pop1s) of hash stored quad patterns
// cmpfunc3qs           compare pixel counts (pop1s) of hash stored small patterns
// cmpfunc4             compare birth times of hash stored genes
// cmpfunc5             compare common genealogy level gene values of hash stored genes
// cmpfunc5c            compare common clonealogy level clone values of hash stored clones
// cmpfunct6            compare according to ancestry in genealogytrace using activity ordering
// cmpfunc7             compare according to ancestry in genealogytrace using population size ordering
//-------------------------------------------------------------- comparison fns -------------------------------------------------------------------------
int cmpfunc (const void * pa, const void * pb) {
   // return ( *(int*)pa - *(int*)pb );
   return ((*(const uint64_t *)pa > *(const uint64_t *)pb)  ? 1 : -1);
}
//.......................................................................................................................................................
int cmpfunc1 ( const void *pa, const void *pb ) {
    const uint64_t *a = (const uint64_t *) pa;
    const uint64_t *b = (const uint64_t *) pb;
    if(a[1] == b[1])
        return a[0] > b[0] ? 1 : -1;
    else
        return (int) (b[1] - a[1]);
}
//.......................................................................................................................................................
int cmpfunc2 (const void * pa, const void * pb) {
    return ( genotypes[*(const int*)pa] > genotypes[*(const int*)pb] ? 1 : -1);
}
//.......................................................................................................................................................
int cmpfunc3 (const void * pa, const void * pb) {
    return ( geneitems[*(const int*)pa].popcount < geneitems[*(const int*)pb].popcount ? 1 : -1);
}
//.......................................................................................................................................................
int cmpfunc3c (const void * pa, const void * pb) {
    return ( cloneitems[*(const int*)pa].popln < cloneitems[*(const int*)pb].popln ? 1 : -1);
}
//.......................................................................................................................................................
int cmpfunc3q (const void * pa, const void * pb) {
    return ( quaditems[*(const int*)pa].pop1s < quaditems[*(const int*)pb].pop1s ? 1 : -1);
}
//.......................................................................................................................................................
int cmpfunc3qs (const void * pa, const void * pb) {
    uint64_t a,b;
    int na,nb;
    a = (uint64_t) *(const int*)pa;
    b = (uint64_t) *(const int*)pb;
    POPCOUNT64C(a,na);
    POPCOUNT64C(b,nb);
    return ( na < nb ? 1 : -1);
}
//.......................................................................................................................................................
int cmpfunc4 (const void * pa, const void * pb) {
   return ( geneitems[*(const int *)pa].firsttime > geneitems[*(const int *)pb].firsttime ? 1 : -1);
}
//.......................................................................................................................................................
int cmpfunc5 (const void * pa, const void * pb) {               // sort according to ancestry in genealogytrace

   int i1,i2,ij1,ij2,j;
   uint64_t gene1,gene2;
   i1=*(const int *)pa; i2=*(const int *)pb;

   for (j=0;j<genealogydepth;j++) {
        ij1 = i1+j*N; ij2 = i2+j*N;
        gene1=working[ij1]; gene2=working[ij2];
        if(gene1!=gene2) return((((gene1 > gene2) && (gene1!=rootgene)) || (gene2==rootgene)) ? 1 : -1);
    }
    return(0);
}
//.......................................................................................................................................................
int cmpfunc5c (const void * pa, const void * pb) {               // sort according to ancestry in clonealogytrace

   int i1,i2,ij1,ij2,j;
   uint64_t birthid1,birthid2;
   i1=*(const int *)pa; i2=*(const int *)pb;

   for (j=0;j<clonealogydepth;j++) {
        ij1 = i1+j*N; ij2 = i2+j*N;
        birthid1=working[ij1]; birthid2=working[ij2];
        if(birthid1!=birthid2) return((((birthid1 > birthid2) && (birthid1!=rootclone)) || (birthid2==rootclone)) ? 1 : -1);
    }
    return(0);
}
//.......................................................................................................................................................
int cmpfunc6 (const void * pa, const void * pb) {               // sort according to ancestry in genealogytrace using activity ordering
   int i1,i2,ij1,ij2,j;
   uint64_t gene1,gene2;
   i1=*(const int *)pa; i2=*(const int *)pb;
   int act1,act2;

   for (j=0;j<genealogydepth;j++) {
        ij1 = i1+j*N; ij2 = i2+j*N;
        gene1=working[ij1]; gene2=working[ij2];
        if(gene1!=gene2)  {
            if((gene1!=rootgene) && (gene2 != rootgene)) {
                if((genedataptr = (genedata *) hashtable_find(&genetable, gene1)) != NULL) act1 = genedataptr->activity; else act1 = 0;
                if((genedataptr = (genedata *) hashtable_find(&genetable, gene2)) != NULL) act2 = genedataptr->activity; else act2 = 0;
                return(act1 > act2 ? 1 : (act1==act2 ? (gene1 > gene2 ? 1 : -1) : -1));
            }
            else return((gene2==rootgene) ? 1 : -1);
        }
    }
    return(0);
}
//.......................................................................................................................................................
int cmpfunc7 (const void * pa, const void * pb) {               // sort according to ancestry in genealogytrace using popln size ordering
   int i1,i2,ij1,ij2,j;
   uint64_t gene1,gene2;
   i1=*(const int *)pa; i2=*(const int *)pb;
   int pop1,pop2;

   for (j=0;j<genealogydepth;j++) {
        ij1 = i1+j*N; ij2 = i2+j*N;
        gene1=working[ij1]; gene2=working[ij2];
        if(gene1!=gene2)  {
            if((gene1!=rootgene) && (gene2 != rootgene)) {
                if((genedataptr = (genedata *) hashtable_find(&genetable, gene1)) != NULL) pop1 = genedataptr->popcount; else pop1 = 0;
                if((genedataptr = (genedata *) hashtable_find(&genetable, gene2)) != NULL) pop2 = genedataptr->popcount; else pop2 = 0;
                return(pop1 > pop2 ? 1 : (pop1==pop2 ? (gene1 > gene2 ? 1 : -1) : -1));
            }
            else return((gene2==rootgene) ? 1 : -1);
        }
    }
    return(0);
}
