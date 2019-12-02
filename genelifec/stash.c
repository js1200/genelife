//
//  stash.c
//  project genelife
//
//  miscellaneous functions associated with stashing, saving and retrieving data
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
//---------------------------------------------------------------- save and retrieve gol... data --------------------------------------------------------
// stash                stash current gol,golg, golb, golr, golgstats in stashgol, stshgolg, stashgolb, stashgolr, stashgolgstats
// label2stash          stash current gol,golg, golb, golr, golgstats from selected labelled component (either cumulatively or individually)
// unstash              retrieve current gol,golg, golb, golr, golgstats from stashed values
// savegols             save current arrays to file
// retrievegols         retrieve arrays to current from file
//.......................................................................................................................................................
void stash(){               // stash current gol,golg
    int ij;
    for (ij=0; ij<N2; ij++) {
        stashgol[ij] = planes[curPlane][ij];
        stashgolg[ij] = planesg[curPlane][ij];
        stashgolb[ij] = planesb[curPlane][ij];
        stashgolr[ij] = planesr[curPlane][ij];
        stashgolgstats[ij] = planesgs[curPlane][ij];
    }
}
//.......................................................................................................................................................
void label2stash(int cumul) { // stash current gol,golg, golb, golr, golgstats from selected labelled component (either cumulatively or individually)
    int ij;
    for (ij=0; ij<N2; ij++) {
        if (labelcc[ij]==0xffffffff) {
            stashgol[ij] = planes[curPlane][ij];
            stashgolg[ij] = planesg[curPlane][ij];
            stashgolb[ij] = planesb[curPlane][ij];
            stashgolr[ij] = planesr[curPlane][ij];
            stashgolgstats[ij] = planesgs[curPlane][ij];
        }
        else if (!cumul) {
            stashgol[ij] = 0ull;
            stashgolg[ij] = 0ull;
            stashgolb[ij] = 0ull;
            stashgolr[ij] = 0ull;
            stashgolgstats[ij] = 0ull;
        }
    }
}
//.......................................................................................................................................................
void unstash(){               // retrieve current gol,golg from stashed values
    int ij;
    for (ij=0; ij<N2; ij++) {
        planes[curPlane][ij] = stashgol[ij];
        planesg[curPlane][ij] = stashgolg[ij];
        planesb[curPlane][ij] = stashgolb[ij];
        planesr[curPlane][ij] = stashgolr[ij];
        planesgs[curPlane][ij] = stashgolgstats[ij];
    }
}
//.......................................................................................................................................................
int savegols( int step, uint64_t gol[], uint64_t golg[], uint64_t golgstats[], uint64_t golb[],uint64_t golr[]) {
    FILE *fp;
    char fname[30];

    sprintf(fname,"gol_gsbr_data%d.ext",step);
    fp = fopen( fname , "wb" );
    
    fwrite(gol  , sizeof(uint64_t) , N2 , fp );
    fwrite(golg , sizeof(uint64_t) , N2 , fp );
    fwrite(golgstats, sizeof(uint64_t) , N2 , fp );
    fwrite(golb , sizeof(uint64_t) , N2 , fp );
    fwrite(golr , sizeof(uint64_t) , N2 , fp );

    fclose(fp);
    return 0;
}
//---------------------------------------------------------------- save data ----------------------------------------------------------------------------
int retrievegols( int step, uint64_t gol[], uint64_t golg[], uint64_t golgstats[], uint64_t golb[],uint64_t golr[]) {
    FILE *fp;
    char fname[30];

    sprintf(fname,"gol_gsbr_data%d.ext",step);
    fp = fopen( fname , "rb" );
    
    fread(gol ,  sizeof(uint64_t) , N2, fp );
    fread(golg ,  sizeof(uint64_t) , N2, fp );
    fread(golgstats, sizeof(uint64_t) , N2, fp );
    fread(golb , sizeof(uint64_t) , N2, fp );
    fread(golr , sizeof(uint64_t) , N2, fp );

    fclose(fp);
    return 0;
}


