//
//  init.c
//  project genelife
//
//  initialization of planes, arrays, hash tables and other structures
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
// initialize_planes    initialize periodic sequence of planes to record rolling time window of up to maxPlanes time points (≤8)
// readFile             read file of gol/golg array (32x32) data
// writeFile            write file of gol/golg array (32x32) data
// testmacros           test macros used to accelerate processing (usually not called): FIRST1INDEX, PATTERN4, PATTERN8
// initialize           initialize simulation parameters and arrays
//----------------------------------------------------------------- initialize_planes -------------------------------------------------------------------
void initialize_planes(int offs[],  int Noffsets) {
    int i,j,idx;
    static int notfirst = 0;

    curPlane = 0;
    newPlane = 1;
    if (notfirst)   return;     // need to fix memory free at two levels unless this fix: no changes in planes structure during run allowed
    notfirst = 1;

    planes[0]  = plane0;  planes[1]  = plane1;     // initialize plane pointers:
    planesg[0] = planeg0; planesg[1] = planeg1;
    planesgs[0] = planegs0; planesgs[1] = planegs1;
    planesb[0]  = planeb0;  planesb[1]  = planeb1;
    planesr[0]  = planer0;  planesr[1]  = planer1;
#if maxPlane > 2
    planes[2]  = plane2;  planes[3]  = plane3;
    planesg[2] = planeg2; planesg[3] = planeg3;
    planesgs[2] = planegs2; planesgs[3] = planegs3;
    planesb[2]  = planeb2;  planesb[3]  = planeb3;
    planesr[2]  = planer2;  planesr[3]  = planer3;
#endif
#if maxPlane > 4
    planes[4]  = plane4;  planes[5]  = plane5;  planes[6]  = plane6;  planes[7]  = plane7;
    planesg[4] = planeg4; planesg[5] = planeg5; planesg[6] = planeg6; planesg[7] = planeg7;
    planesgs[4] = planegs4; planesgs[5] = planegs5; planesgs[6] = planegs6; planesgs[7] = planegs7;
    planesb[4]  = planeb4;  planesb[5]  = planeb5;  planesb[6]  = planeb6;  planesb[7]  = planeb7;
    planesr[4]  = planer4;  planesr[5]  = planer5;  planesr[6]  = planer6;  planesr[7]  = planer7;
#endif

    if (!(diagnostics & diag_offset_statistics)) return;
    
    if(Noffsets%3 !=0) fprintf(stderr,"Size of offsets array not a multiple of 3 as expected.");
    Noff = Noffsets/3;        // Noff global
    if(Noff>24){
        fprintf(stderr,"Too many offsets!  Max=24");
        exit(1);
    }
    numHisto = Noff;
    histo = (int *) calloc(numHisto,sizeof(int));

    // install offsets:
    offsets = (int **) calloc(Noff,sizeof(int *));
    for(i=0; i<Noff; i++) offsets[i] = (int *) calloc(3,sizeof(int)); // each containing xoff, yoff, toff.
    for(idx=0,i=0; i<Noff; i++){
        for(j=0; j<3; j++){
            offsets[i][j] = offs[idx];
            idx++;
        }
    }

    // compute number of planes from toff = 3rd element of each offest vec:
    int tall,tmx = 0;
    int toff, tmn = Noffsets;
    for(i=0; i<Noff; i++){
        toff = offsets[i][2];
        if(toff>tmx) tmx = toff; // it does not make sense to have +ve values (this would look into future) so tmx = 0
        if(toff<tmn) tmn = toff;
    }
    if(tmx>0)   {                // exit if positive values of z in (x,y,z) offsets have been entered : these look into future
        fprintf(stderr,"Error: offsets looking into future not allowed ------- tmx = %d, tmn = %d",tmx,tmn);
        exit(1);
    }    tall = tmx-tmn;
    // numPlane = 2 + tall;    // numPlane >= 2
    if(tall>4) numPlane = 8;                            // only numPlane values which are powers of two are allowed
    else if (tall > 2) numPlane = 4;
    else numPlane = 2;
    if (numPlane>maxPlane) {
        fprintf(stderr,"Not enough planes defined by maxPlane for given offsets (need > %d)\n",numPlane);
        exit(1);
    }
}
//---------------------------------------------------------------- readFile and writeFile ---------------------------------------------------------------
int readFile(char * code, const char *fileName) {
  FILE *file;
  char tmp;
  int cnt;
  file = fopen(fileName, "r");
  cnt = -1;
  do
  {
    *code++ = tmp = (char)fgetc(file);
    cnt++;
  } while(tmp != EOF);
  fclose(file);
  return cnt;
}
//.......................................................................................................................................................
int writeFile(char *fileName)  {            // initialize 32x32 genepat file with all empty sites
    FILE *file;
    int ij,error;
    file = fopen(fileName, "w");
    error = 0;
    for(ij=0;ij<32*32;ij++) {
        error=fputc(0,file);
        if (error) break;
    }
    fclose(file);
    return error;
}
//---------------------------------------------------------------- initialize ---------------------------------------------------------------------------
void testmacros() {
    int j,k;
    uint64_t g;
    g=1ull;                                                              // test of FIRST1INDEX
    for(j=0;j<10;j++) {
        FIRST1INDEX(g,k);
        fprintf(stderr,"test of first1index cnt %d val %"PRIx64" index %d\n",j,g,k);
        g*=42;
    }

    g=1ull;                                                              // test of PATTERN4
    for(k=0;k<10;k++) {
        int found;
        for(j=0;j<16;j++) {
            PATTERN4(g,j,found);
            fprintf(stderr,"test of pattern4 pat %x val %"PRIx64" found? %d\n",j,g,found);
        }
        g*=42;
    }

    g=1ull;                                                              // test of PATTERN8
    for(k=0;k<10;k++) {
        int found;
        for(j=0;j<64;j++) {
            PATTERN8(g,j,found);
            fprintf(stderr,"test of pattern8 pat %x val %"PRIx64" found? %d\n",j,g,found);
        }
        g*=42;
    }
    
    for (j=0;j<257;j++) {
        fprintf(stderr,"test of patterns j %5d log2upper(j) %5d log2upper(sqrtupper(j)) %5d\n",j,log2upper(j),log2upper(sqrtupper(j)));
    }
}
//........................................................................................................................................................
void initialize(int runparams[], int nrunparams, int simparams[], int nsimparams) {
    int hcnt;
    int ij,ij1,i0,j0,i,j,Nf,k,cnt,icf,nstartgenes;
    unsigned int ncodingin;
    uint64_t g;
    static unsigned int rmask = (1 << 15) - 1;
    static int notfirst = 0;
    uint64_t startgenes[16];
    char *golgin;

    rulemod = runparams[0];
    repscheme = runparams[1];
    selection = runparams[2];
    overwritemask = runparams[3];
    survivalmask = runparams[4];
    colorfunction = runparams[5];
    initfield = runparams[6];
    birthmask=runparams[7];
    ancselectmask=runparams[8];
    colorfunction2 = runparams[9];

    pmutmask = (unsigned int) simparams[0];                                      // low values of pmutmask <32 are interpreted as -nlog2pmut
    if(pmutmask<32) pmutmask = (0x1 << (32-pmutmask)) - 0x1ull;                  // NB if pmut==0, pmutmask==zero, no mutation.
    initial1density = simparams[1];
    initialrdensity = simparams[2];
    
    ncodingin = simparams[3];                                                    // used in selection  4,5,6,8,
    ncoding = ncodingin & 0xff;
    ncoding2 = (ncodingin>>8) & 0xff;
    if (ncoding > 64) { fprintf(stderr,"value %d of ncoding is out of range\n",ncoding);ncoding = 64;}
    if (selection==8 || selection == 9) if (ncoding<1 || ncoding>4) ncoding = 4; // ncoding range restriction for selection = 8 ie 16*ncoding bits of gene used
    if (selection<8) codingmask = (1ull<<ncoding2)-1ull;                         // bit mask corresponding to ncoding2 bits, only used in connection with add2ndmask1st
    else if (selection<10) codingmask = (1ull<<ncoding)-1ull;                    // coding mask used to encode number of bits per LUT (1-4)
    
    startgenechoice = simparams[4];

    if(nsimparams > 5) ranseed = simparams[5];
    
    srand(ranseed); // Range: rand returns numbers in the range of [0, RAND_MAX ), and RAND_MAX is specified with a minimum value of 32,767. i.e. 15 bit
    randstate[0] = rand();randstate[1] = rand();                                 // state vector for dedicated 64-bit random number generator macro RAND128P
    cnt = 0;
    totsteps = 0;
    totdisp = 0;
    statcnts = 0;
    quadrants = -1;
    rbackground = 0;
    quadcollisions = 0;
    randominflux = 0;
    vscrolling = last_scrolled = vscrolly = 0;
    quadrants = -1;
    gene0=0ull;                                                                 // normally default gene is 0ull : unused when gol state not live
    nstartgenes = 8;
    
    // testmacros();                                                             // test macros used to accelerate processing
    // writeFile("genepat.dat");                                                 // can be used to initialize formatted template for gene input of 32x32 array


    fprintf(stderr,"___________________________________________________________________________________________\n");
    fprintf(stderr,"_________________________________ genelife simulation _____________________________________\n");
    fprintf(stderr,"runparams %d %d %d %d %d %d %d %d %d %d\n",runparams[0],runparams[1],runparams[2],
                    runparams[3],runparams[4],runparams[5],runparams[6],runparams[7],runparams[8],runparams[9]);
    fprintf(stderr,"simparams %d %d %d %d %d %d\n",simparams[0],simparams[1],simparams[2],simparams[3],simparams[4],ranseed);
    fprintf(stderr,"pmutmask %x (NB 0 means no mutation)\n",pmutmask);

    switch (selection) {                                                         // initialize starting genes depending on selection model, encoding and symmetry
        case 0:  for (k=0;k<4;k++) { startgenes[k] = 0xf0f0f0f0f0f0f0f0; startgenes[k+4] = 0x0f0f0f0f0f0f0f0f;} break;
        case 1:  for (k=0;k<8;k++)   startgenes[k] = ((0x1ull<<k*3)-1ull)<<20;break;
        case 2:  for (k=0;k<8;k++)   startgenes[k] = ((0x7ull>>(k&3))<<61) | 0x606ull;break;
        case 3:  for (k=0;k<8;k++)   startgenes[k] =(((0x1ull<<20)-1ull)<<20)+((0x1ull<<k)-0x1ull);break;
        case 4:  for (k=0;k<8;k++)  {g = 0xff0ull; startgenes[k] = k<4 ? g+1 : ~((g<<16));} break;
        case 5:  for (k=0;k<8;k++)  {g = 0xf0ull + k; startgenes[k] = k<4 ? g : (~g)|(0xfull<<16);} break;
        case 6:
        case 7:  for (k=0;k<8;k++) startgenes[k]=(0x1ull<<(4+k*8))-1ull; break;

        case 8:  if (((repscheme>>4)&0x7)==7)
                      genegol[selection-8] = (codingmask<<((8+3-1)*ncoding))|(codingmask<<((8+2-1)*ncoding))|(codingmask<<((2-1)*ncoding))|(codingmask<<((3-1)*ncoding));
                 else genegol[selection-8] = (codingmask<<((8+3-1)*ncoding))|(codingmask<<((2-1)*ncoding))|(codingmask<<((3-1)*ncoding));
                 for (k=0;k<8;k++)   startgenes[k] = ((0x7ull>>(k&3))<<61) | genegol[selection-8];          // put up to 3 extra bits at top to ensure all nr 1s values occupied
                 break;
        case 9:  if (((repscheme>>4)&0x7)==7) genegol[selection-8] = 0xba32ull;                             // extended rule for scissors-paper-stone-well gliders
                 else                         genegol[selection-8] = 0xb32ull;                              //GoL rule for survival in totalistic LUT case, variable length encoding
                 for (k=0;k<8;k++)   startgenes[k] = ((0x7ull>>(k&3))<<61) | genegol[selection-8];          // put up to 3 extra bits at top to ensure all nr 1s values occupied
                 break;
        case 10: if (((repscheme>>4)&0x7)==7) genegol[selection-8] = (0x7ull<<2)|(0xfull<<5)|(0x7ull<<34)|(0xfull<<37); //GoL rule for S23 B23 in corner/edge dist LUT case, fixed length encoding
                 else                         genegol[selection-8] = (0x7ull<<2)|(0xfull<<5)|(0xfull<<37);  //GoL rule for S23 B3 in corner/edge dist LUT case, fixed length encoding
                 for (k=0;k<8;k++)   startgenes[k] = ((0x7ull>>(k&3))<<61) | genegol[selection-8];
                 break;
        case 11: if (((repscheme>>4)&0x7)==7) genegol[selection-8] = 0xbfa73f27ull;                         //GoL rule for S23 B23 in corner/edge dist LUT case, variable length encoding
                 else                         genegol[selection-8] = 0xbf3f27ull;                           //GoL rule for S23 B3 in corner/edge dist LUT case, variable length encoding
                 for (k=0;k<8;k++)   startgenes[k] = ((0x7ull>>(k&3))<<61) | genegol[selection-8];
                 break;
        case 12: if (((repscheme>>4)&0x7)==7) genegol[selection-8] = 0xfull|(0x7full<<4)|(0xfull<<32)|(0x7full<<36); //GoL rule for S23 B23 in canonical rotation case, fixed length encoding
                 else                         genegol[selection-8] = 0xfull|(0x7full<<4)|(0x7full<<36);     //GoL rule for 2,3 s and 3 b in canonical rotation case:  4,7,10,7,4 configs s=2,3,4,5,6
                 for (k=0;k<8;k++)   startgenes[k] = ((0x7ull>>(k&3))<<61) | genegol[selection-8];
                 break;
        case 13: if (((repscheme>>4)&0x7)==7) genegol[selection-8] = 0x04cda7bbb73b3727;                    //GoL rule for S23 B23 with modular encoding for canonical rotation case
                 else                         genegol[selection-8] = 0x00cd00bbb73b3727;                    //GoL rule for S23 B3 with modular encoding for canon rot'n case: 4,7,10,7,4 configs s=2,3,4,5,6
                 for (k=0;k<8;k++)   startgenes[k] = ((0x7ull>>(k&3))<<61) | genegol[selection-8];
                 break;
        case 14: if (((repscheme>>4)&0x7)==7) genegol[selection-8] = (0x3full<<3)|(0x3ffull<<9)|(0x3full<<(32+3))|(0x3ffull<<(32+9)); //GoL rule for S23 B23 for 2D_sym case, fixed length encoding
                 else                         genegol[selection-8] = (0x3full<<3)|(0x3ffull<<9)|(0x3ffull<<(32+9)); //GoL rule for S23 B3 for 2D_sym case with 1,2,6,10,13 configs for s=0,1,2,3,4
                 for (k=0;k<8;k++)   startgenes[k] = ((0x7ull>>(k&3))<<61) | genegol[selection-8];
                 break;
        case 15: if (((repscheme>>4)&0x7)==7) genegol[selection-8] = 0x03ffffb7b3a33323;                    //GoL rule for S2(all)3(1st 6) B23 with modular encoding for 2D_sym case (NB smask 0 anyway)
                 else                         genegol[selection-8] = 0x03f3ffb7b3373323;                    //GoL rule for S23 B3 with modular encoding for 2D_sym case for s=0,1,2,3,4
                 for (k=0;k<8;k++)   startgenes[k] = ((0x7ull>>(k&3))<<61) | genegol[selection-8];
                 break;

        default: for (k=0;k<8;k++) startgenes[k]=(0x1ull<<(4+k*8))-1ull;
    }
    
    if (diagnostics & diag_general_statistics) {                                  // general statistics
        if (livesites !=NULL) {free(livesites);livesites = NULL;}
        if (genestats !=NULL) {free(genestats);genestats = NULL;}
        if (stepstats !=NULL) {free(stepstats);stepstats = NULL;}
        if (configstats != NULL) {free(configstats);configstats = NULL;}
    
        arraysize = startarraysize;
        livesites = (int *) malloc(arraysize * sizeof(int));
        genestats = (int *) malloc(arraysize * 4 * sizeof(int));
        stepstats = (int *) malloc(arraysize * 10 * sizeof(int));
        if (nhistG==nstatG) configstats = (int *) malloc(arraysize * Noff * sizeof(int));
    }
    
    curPlane = 0;                                                                 // if we rerun initialize, we want to restart plane cycling from zero
    newPlane = 1;
    gol = planes[curPlane];
    golg = planesg[curPlane];
    golgstats = planesgs[curPlane];
    golb = planesb[curPlane];
    golr = planesr[curPlane];
    
    if(notfirst) {
        if(diagnostics & diag_hash_genes) hashtable_term(&genetable);
        if(diagnostics & diag_hash_patterns) {
            hashtable_term(&quadtable);
            memset(smallpatts,0,sizeof(smallpatt)*65536);
            // for (ij=0;ij<65536;ij++) smallpatts[ij].activity = 0;              // initialize small pattern table to no patterns hit, already done
        }
        if(diagnostics & diag_hash_clones) hashtable_term(&clonetable);
    }
    
    if(diagnostics & diag_hash_genes) hashtable_init(&genetable,sizeof(genedata),N2<<2,0);   // initialize dictionary for genes
    if(diagnostics & diag_hash_patterns) hashtable_init(&quadtable,sizeof(quadnode),N2<<2,0);// initialize dictionary for quadtree patterns
    if(diagnostics & diag_hash_clones) hashtable_init(&clonetable,sizeof(clonedata),N2<<4,0);// initialize dictionary for clones
    
    notfirst = 1;
    if (initfield==1) {                              // input from file genepat.dat with max size of 32*32 characters
        golgin = (char *) malloc(32* 32 * sizeof(char));
        icf=readFile(golgin, "genepat.dat");
        if (icf != 32*32) {
            icf = 0;
            fprintf(stderr,"error reading file, %d not 32*32 chars\n",icf);
        }
        for (ij=0; ij<N2; ij++) {
            gol[ij] = 0ull;
            golg[ij] = gene0;
            golgstats[ij] = 0ull;
            golb[ij] = 0ull;
            golr[ij] = 0ull;
        }
        for (ij1=0; ij1<32*32; ij1++) {
            if (N<(32)) {fprintf(stderr,"Error, array dim %d too small for file array dim %d\n",N,32);break;};
            ij=(N>>1)-16+(ij1&0x1f)+ N*((N>>1)-16+(ij1>>5));
            if (golgin[ij1] > 0)    {                // if live cell
                gol[ij] = 1ull;
                if(golgin[ij1] <= 8 ) golg[ij] = startgenes[golgin[ij1]-1];
                else if (golgin[ij1]>='0' && golgin[ij1]<'8') golg[ij] = startgenes[golgin[ij1]-'0'];
                else golg[ij] = startgenes[7];
                cnt++;
                golb[ij] = rootclone + ij;               // initialize clone to new clone at time 0
            }
            // if (golg[ij] == 0 && gol[ij] != 0) fprintf(stderr,"zero gene at %d\n",ij);
        }

    }
    else if (initfield>=0) {                         // initfield gives linear size of random block for initialization (0 => full frame, as before)
        Nf = initfield;
        if (Nf==0 || Nf>N) Nf=N;
        for (ij=0; ij<N2; ij++) {
            gol[ij] = 0ull;
            golg[ij] = gene0;
            golgstats[ij] = 0ull;
            golb[ij] = 0ull;
            golr[ij] = 0ull;
        }
        i0 = j0 = (N>>1)-(Nf>>1);
        for (i=0; i<Nf; i++) {
            for (j=0; j<Nf; j++) {
                ij=i0+i+N*(j0+j);
                gol[ij] |= ((rand() & rmask) < initial1density)?1ull:0ull;
            }
        }
        for (ij=0; ij<N2; ij++) {
            g = 0ull;
            if (gol[ij]) {                          //  fill with random genome g or randomly chosen startgene depending on initialrdensity
                if (((unsigned) rand() & rmask) < initialrdensity) {for (k=0; k<64; k++) g = (g << 1) | (rand() & 0x1);g=gene0^g;}
                else if (startgenechoice == nstartgenes) g = startgenes[0xf & rand() & (nstartgenes-1)];
                else if (startgenechoice > nstartgenes) fprintf(stderr,"startgenechoice %d out of range\n",startgenechoice);
                else g = startgenes[0xf & startgenechoice & (nstartgenes-1)];
                cnt++;
            }
            golg[ij] = g;
            golb[ij] = rootclone + ij;
        }
    }
    else {                                          // initfield < 0, use array values from stashed
        for (ij=0; ij<N2; ij++) {
            gol[ij] = stashgol[ij];
            golg[ij] = stashgolg[ij];
            golb[ij] = stashgolb[ij];
            golr[ij] = stashgolr[ij];
            golgstats[ij] = stashgolgstats[ij];
        }
    }

    if(diagnostics & diag_scrolling_trace) {
        for (i=0;i<N;i++) npopulation[i]=0;         // initialize scrolling total population trace
        for (ij=0; ij<N2; ij++) {
            poptrace[ij]=rootgene;                  // initialize population traces to root gene
            acttrace[ij]=rootgene;                  // initialize activity traces to root gene
            acttraceq[ij]=rootgene;                 // initialize activity traces of patterns to root gene
            acttraceqt[ij]=0;                       // initialize activity traces of types of patterns to 0
            genealogytrace[ij] = rootgene;          // initialize genealogy traces to root gene
        }
    }

    if(diagnostics & diag_longtime_trace) {
        for (i=0;i<N*nNhist;i++) {
            npopulation1[i]=0; // initialize longer beginning total population trace
            nnovelcells[i]=0;  // initialize longer beginning total novel cell count trace
        }
        for (ij=0; ij<N2*nNhist; ij++) {
            poptrace1[ij]=rootgene;                 // initialize long population traces to root gene
            acttrace1[ij]=rootgene;                 // initialize long activity traces to root gene
            acttraceq1[ij]=rootgene;                // initialize long activity traces of patterns to root gene
            acttraceqt1[ij]=0;                      // initialize long activity traces of types of patterns to 0
        }
    }
    
    if(diagnostics & diag_hash_genes) {
        for (ij=0; ij<N2; ij++) {
            if(gol[ij]) {
                hashaddgene(ij,golg[ij],rootgene,golb+ij,rootclone+ij,0x1ull);   // totsteps=0 so parentid is spatial ij + flag rootclone
            }
        }
                                                // enumerate gene hash keys and values
        hcnt=hashtable_count(&genetable);
        genotypes = hashtable_keys( &genetable );
        fprintf(stderr,"population size %d with %d different genes\n",cnt,hcnt);
    }

    // if(diagnostics & diag_hash_patterns) qimage = quadimage(newgol,&patt,log2N); // quadtree hash of entire image
    if(diagnostics & diag_component_labels) ncomponents=extract_components(gol,golg);
}
