//
//  genelife.h
//  project genelife
//
// declaration of global constants and variables
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
//------------------------------------------------------- acknowledgements ------------------------------------------------------------------------------
// Thanks to external contributions:
//
// Sean Eron Anderson
// Bit Twiddling Hacks © 1997-2005 (all snippets in public domain)
// https://graphics.stanford.edu/~seander/bithacks.html
//
// Mattias Gustavsson
// hashtable.h, in:
// https://github.com/mattiasgustavsson/libs
// dual license, MIT / Unlicense
// choosing Unlicense license (http://www.unlicense.org)
//
// Stackoverflow:
// license: creative commons https://creativecommons.org/licenses/by-sa/3.0/
//
// Anthony Bourdain
// https://stackoverflow.com/questions/27159322/rgb-values-of-the-colors-in-the-ansi-extended-colors-index-17-255
//
// MetallicPriest
// https://stackoverflow.com/questions/6943493/hash-table-with-64-bit-values-as-key/33871291
//
// Wikipedia:
//
// "Xorshift" algorithm rewritten here as inline macro
// https://en.wikipedia.org/wiki/Xorshift
// Vigna, Sebastiano. "xorshift*/xorshift+ generators and the PRNG shootout". Retrieved 2014-10-25.
// https://prng.di.unimi.it
// "Hamming Weight" popcount4c algorithm
// https://en.wikipedia.org/wiki/Hamming_weight
// "Disjoint-set data structure" (n.d.) Retrieved Oct 16, 2019
// https://en.wikipedia.org/wiki/Disjoint-set_data_structure
// The quadtree code for storing patterns benefited from but is different from that employed in hashlife
// https://en.wikipedia.org/wiki/Hashlife
//------------------------------------------------------------------------------------------------------------------------------------------------------
#ifndef genelife_h
#define genelife_h

#ifndef OEX
#define OEX extern
#define NOHASH
#define INIT(x,...) x /* __VA_ARGS__ */
#else
#define OEX
#define INIT(x,...) x = __VA_ARGS__
#endif  /* OEX */

#define INLINE

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#define __STDC_FORMAT_MACROS
#include <inttypes.h>
#include <time.h>
#include <math.h>
#include <string.h>

#include "genelife_size.h"              // definition of model size via log2N and derived parameters N,N2,NLM,NLC,Nmask,N2mask: N is length of array side
//--------------------------------------------------------- main parameters of model --------------------------------------------------------------------
OEX unsigned int INIT(rulemod,1);       // determine whether to modify GoL rules
OEX int INIT(selection,1);              // 0-7 fitness of 2 live neighbours, 8-15 lut symmetry and encoding scheme
                                        // 0.  integer value                        1. number of ones         2. scissors-stone-well-paper: wins over left 1-1-2-2
                                        // 3.  scissors-stone-well-paper 1-1-1-1    4. 5. predator prey       6. two target coding     7. selection always fails (no result)
                                        // 8.  genes encode lut based on sum s for survival/birth:  s 1-8 (8 configs): ncoding determines nr of fixed position bits per LUT entry
                                        // 9.  like 8 but position-independent encoding of lut entries: ncoding 1 bit only
                                        // 10. genes encode lut based on sum and edge-centred count (s,se) for survival/birth: s 1-7 se 0-4 (2+3+4+5+4+3+2=23 configs): ncoding 1-bit only
                                        // 11. like 10 but position-independent encoding of lut entries
                                        // 12. genes encode lut based on sum s and canonical rotations survival/birth sum 2-6 (4,7,10,7,4 canon.rotns) : ncoding 1 bit only
                                        // 13. like 12 but but position-independent encoding of lut entries
                                        // 14. genes encode lut based on 2D rotation/reflection symmetries for sum 0-4 only (1,2,6,10,13 cann.rotns) 2x32 bits: ncoding 1 bit only
                                        // 15. like 14 but position-independent encoding of lut entries
OEX unsigned int INIT(repscheme,1);     // replication scheme with separate meaning for selection 0-7 (update_23) and 8-15 other update fns
                                        // for selection 0-7: see #define R_... 1st section below
                                        // --  lowest 10 bits define 5 pairs of bit options for:
                                        // --  0,1 select birth 2,3 neighbour choice 4,5 enforce birth 6,7 2nd,1st neighbour genes 8,9 no successive nonGoL
                                        // --  bits 10-13 specify 4-bit mask for 2-live neighbour birth onl for those of the 4 canonical configurations set
                                        // --  bits 14-20 activate quadrant exploration of specific subset of the first 5 pairs of repscheme and survival and overwrite masks
                                        // --  bit 21 parentdies
                                        // for selection 8-15: see #define R... 2nd section below
OEX unsigned int INIT(survivalmask,0x3);// for selection 0-7 survive mask for two (bit 1) and three (bit 0) live neighbours
                                        // for selection=8,9 it is 16-bit and for 10-11 it is 19-bit and 12-13 32-bit survival mask to restrict luts
OEX unsigned int INIT(birthmask,0x0);   // for selection 0-7 unused, birth control is done via bits 0,1,4,5 of repscheme
                                        // for selection=8,9 it is 16-bit and for 10-11 it is 19-bit and 12-13 32-bit birthmask to restrict luts
OEX unsigned int INIT(overwritemask,0x3);// for selection 0-7 bit mask for 4 cases of overwrite: bit 0. s==3  bit 1. special birth s==2
                                        // for selection=8,9 it is 16-bit and for 10-11 it is 19-bit and 12-13 32-bit birthmask to restrict luts
OEX unsigned int INIT(ancselectmask,0xff);// whether to use selection between genes to determine ancestor: for each birth rule (0 use positional choice via selectdifft)
OEX int INIT(ncoding,1);                // byte 0 of python ncoding : number of coding bits per gene function
OEX int INIT(ncoding2,0);               // byte 1 of python ncoding: number of coding bits per gene function for masks in connection with repscheme add2ndmask1st R_6,7
OEX unsigned int pmutmask;              // binary mask so that prob of choosing zero is pmut = pmutmask/2^32. If value<32 interpret as integer -log2(prob)
//-----------------------------------------------------------initialization and color parameters---------------------------------------------------------
OEX int INIT(initial1density,(1<<15)>>1);// initial density of ones in gol as integer value, divide by 2^15 for true density
OEX int INIT(initialrdensity,(1<<15)>>1);// initial density of random genes in live sites, divide by 2^15 for true density
OEX int INIT(startgenechoice,8);        // selection for defined starting genes 0-8 (8 is random 0-7) otherwise choose only particular startgene
OEX int INIT(initfield,0);              // 0 input from random field of random or start genes, 1 input from file genepat.dat of 32x32 indexes to start genes
                                        // value i>1: initialized to random field on central ixi block only, outside this zero.
                                        // value<0: initialize with gol and golg values
OEX int INIT(colorfunction,0);          // color function choice of 0: hash or 1: functional (color classes depends on selection parameter) & if last step was non GoL change yellow
                                        // 2: use numbner of live neighbours to colour  3: colour with 3 colours and 3 levels using information in golgstats
                                        // 4: activities 5: populations 6: genealogies without time 7: genealogies with time brightness by activity 8: gliders
                                        // 9: connected components 10 : connected component activities 11: genealogy based individual colors 12: genetic glider det
OEX int INIT(colorfunction2,-1);        // colorfunction for second window: as above, but -1 means same value as colorfunction : only one fn call made
OEX int INIT(colorupdate1,1);           // flag to enable routine print statements to terminal during run : linked to colorfunction display in python
OEX int INIT(connectedprints,0);        // flag to enable printouts of statistics of connected components (every 10 steps)
OEX int INIT(ancestortype,0);           // display and return genealogies via first ancestor (0), clonal ancestry (1) or first & clonal in 2 win (2)
OEX int INIT(parentdies,0);             // model variant enhancing interpretation of non-proliferative birth as movement (1) or default (0): set in repscheme
//------------------------------------------------------- masks for named bits of control words ---------------------------------------------------------
//.................................................. masks for named repscheme bits (selection 0-7) .....................................................
#define R_0_2sel_3live     0x1          /* 1: for 3-live-n birth, employ selection on two least different live neighbours for ancestor */
#define R_1_2sel_2live     0x2          /* 1: allow 2-live-n birth, employ selection on 2 live neighbours for ancestor */
#define R_2_canonical_nb   0x4          /* 1: choose live neighbour at zero bit in canonical rotation 0: choose most difft position */
#define R_3_neutral_pos    0x8          /* 1: for 2-live-nb birth, choose canonical position as specified by R_2_canonical_nb rather than doing selection */
#define R_4_enforce3birth  0x10         /* 1: enforce birth for 3 live nbs (with ancestor most difft) in case 2-select fails (req. R_0_2sel_3live==1 for effect) */
#define R_5_enforce2birth  0x20         /* 1: enforce birth for 2 live nbs : choose 1st live nb from top left of canonical config as ancestor (asym!) */
#define R_6_2ndnb_genes    0x40         /* 1: execute genetically encoded 2nd neighbours of live neighbours to determine birth & 1st nb ancestor */
#define R_7_1stnb_masks    0x80         /* 1: check genetically encoded masks to mask out certain live nb pos's: rule as if these not there */
#define R_8_nongolstat     0x100        /* 1: enforce GoL rule if state of central cell was last changed by a non GoL rule */
#define R_9_nongolstatnbs  0x200        /* 1: enforce GoL rule if state of any cell in nbs was last changed by a non GoL rule */
#define R_10_2birth_k0     0x400        /* 1: bit position at start of k1-4 mask for selective subset of 2-births */
#define R_10_13_2birth_k4  0x3c00       /* 1: enforce birth for 2 live nbs canonical config for one of k= 1,2,3,4, next 4 bits: choose 1st live nb from TL (asym!) */
#define R_14_parentdies_23 0x4000       /* 1. parent is forced to die on birth 0: not: only used for selection 0-7 ie update_23 */
#define R_15_dummy         0x8000       /* not yet used but connected up to graphical interface */
#define R_quadrant         0x7f0000     /* 1: quarter the spatial domain with one or more of 7 pairs of repscheme bits ie 4 different values */
#define R_16_quadrant_sele 0x10000      /* q0 1: quarter the spatial domain with selection enable values for 2,3 live nbs: only in update ie for selection<8 */
#define R_17_quadrant_posn 0x20000      /* q1 1: quarter the spatial domain with selection enable values for 2,3 live nbs: only in update ie for selection<8 */
#define R_18_quadrant_enfb 0x40000      /* q2 1: quarter the spatial domain with enforce birth values for 2,3 live nbs: only in update ie for selection<8 */
#define R_19_quadrant_2nb1 0x80000      /* q3 1: quarter the spatial domain with 1st nb masks and/or 2nd nb addition: only in update ie for selection<8 */
#define R_20_quadrant_ngol 0x100000     /* q4 1: quarter the spatial domain with last non gol rule and/or non gol created state: only in update ie for selection<8 */
#define R_21_quadrant_surv 0x200000     /* q5 1: quarter the spatial domain with survival values for 2,3 live nbs: only in update ie for selection<8 */
#define R_22_quadrant_over 0x400000     /* q6 1: quarter the spatial domain with overwrite values for 2,3 live nbs: only in update ie for selection<8 */
//.................................................. masks for named repscheme bits LUT (selection 8-15) .................................................
#define R_0_nb_majority   0x1           /* 1: Majority rule (reinterprets next bit nb_OR_AND as 1: >= 0: >) 0: No majority rule (i.e. And or Or) */
#define R_1_nb_OR_AND     0x2           /* If Majority=0 1: OR of neighbours determines genetic LUT in selection 8,10,12,14 0: AND of neighbours Otherwise 1:>= 0:>*/
#define R_2_canonical_nb  0x4           /* 1: choose live neighbour at zero bit in canonical rotation 0: choose most difft position */
#define R_3_parentdies    0x8           /* 1: parent is forced to die on birth 0 not. Only used for selection 8-15 */
#define R_46_repselect    0x70          /* 0-7 choice of selection mechanism for LUT genes : 0: min 1: max 2: min 1s 3: max 1s 4: neutral 5: neutral difft 6,7: c-S-2B */
#define R_47_repselect    0xf0          /* +8-15 choice of golr selection mechanism : 8: 9: 10: 11: 12: 13: 14: 15:  */
#define R_7_golr_select   0x80          /* 1: choose 8-15 above 0: choose 0-7 above */
#define R_810_disambig    0x700         /* 0-7 choice of different disambiguation mechanisms for symmetric canonical rotations */
#define R_11_random_resln 0x800         /* 1: random choice amongst selected live neighbours 0: deterministic choice based on gene content and position */
#define R_12_survivalgene  0x1000       /* 1: survival gene chosen from central existing gene 0: survival gene taken from neighbours as in birth */
//........................................ status flag bits for recording site status in golgstats array...................................................
#define F_s_live      0x7ull            /* s value mod 8 (number of live neighbors) for selection 8-15  and separate bits below for selection 0 to 7 */
#define F_1_live      0x1ull            /* bit is bit0 of s for selection 8-15 or 1 if exactly 1 live neighbours for selection 0-7 : currently not set */
#define F_2_live      0x2ull            /* bit is bit1 of s for selection 8-15 or 1 if exactly 2 live neighbours for selection 0-7 */
#define F_3_live      0x4ull            /* bit is bit2 of s for selection 8-15 or 1 if exactly 3 live neighbours */
#define F_birth       0x8ull            /* bit is 1 if birth (includes overwriting of genes for some parameter values) */
#define F_mutation    0x10ull           /* bit is 1 if a mutation event occured */
#define F_disambig    0x80ull           /* bit is 1 if the disambiguate routine is used: nbest>1, nsame > 1 */
#define F_survival    0x20ull           /* bit is 1 if last step was a 1->1 gol survival */
#define F_death       0x40ull           /* bit is 1 if last step involved the death of a gene ie 1->0 gol transition */
#define F_golstate    0x100ull          /* bit is 1 if gol state is 1 */
#define F_golchange   0x200ull          /* bit is 1 if state changed at last step */
#define F_nongolchg   0x400ull          /* bit is 1 if state when produced (ie changed to) was made by a non GoL rule */
#define F_notgolrul   0x800ull          /* bit is 1 if last step not a GoL rule */
#define F_survmut     0x1000ull         /* bit is 1 if mutation or survival from non-replicated mutant: mutation(t) or mutation(t-1)&survival(t) */
#define F_parent      0x2000ull         /* bit is 1 if individual that was at this site was parent/ancestor/genetic donor of a new individual born in last step */
#define F_parentaldeath 0x4000ull       /* bit is 1 if individual died at this step only because it was a parent under parentdies=1 option */
#define F_dummy       0x8000ull         /* free: not yet used */
#define F_livenbs     0xff0000ull       /* mask for storing configuration of live neighbours : clockwise from top-left neighbour (NW) (sel 0-7 for s>1, sel 8-15 for s>0) */
#define F_2select     0x1000000ull      /* bit is 1 if the 2 live neighbour selection routine was employed : selection 0-7 only */
#define F_3g_same     0x2000000ull      /* bit is 1 if exactly 3 live nbs and all 3 have same gene : only for selection 0-7 */
//...........................................................diagnostic control .........................................................................
#define diag_all                0xffff  /* all diagnostics active, remove bits from this mask to run without them */
#define diag_hash_genes         0x1     /* enable hash storage of all genes encountered in simulation */
#define diag_hash_patterns      0x2     /* enable hash storage of all patterns encountered in simulation */
#define diag_hash_clones        0x4     /* enable hash storage of all clones encountered in simulation */
#define diag_activities         0x8     /* enable activity statistics recording for genes and patterns */
#define diag_simp_genealogies   0x10    /* enable simple genealogies with first or most recent ancestor */
#define diag_clone_genealogies  0x20    /* enable genealogies by hashed clones */
#define diag_component_labels   0x40    /* enable spatial connected component labelling, mapping and coloring */
#define diag_offset_statistics  0x80    /* enable collection of offset histogram statistics: histo,numHisto,offsets,Noff */
#define diag_scrolling_trace    0x100   /* enable scrolling time tracing of activities for genes and patterns,populations,genealogies */
#define diag_longtime_trace     0x200   /* enable longer time tracing of activities for genes and patterns,populations (poss.genealogies) */
#define diag_general_statistics 0x400   /* enable collection of general statistics: livesites,genestats,stepstats,configstats */
#define diag_info_transfer_hist 0x800   /* enable collection of glider information transfer histogram in 8 directions  N E S W NE SE SW NW */
OEX unsigned int INIT(diagnostics,diag_all);// bit mask for all diagnostics as defined by following constants
//----------------------------------------------------------hash table implementation of python style dictionary---------------------------------------
#ifndef NOHASH
#define HASHTABLE_IMPLEMENTATION        /* uses Mattias Gustavsson's hashtable (github) for unsigned 64 bit key dictionary */
#endif
#define HASHTABLE_U64 uint64_t          /* define the hashtable 64bit unsigned int type as uint64_t */
#define HASHTABLE_U32 uint32_t          /* define the hashtable 32bit unsigned int type as uint32_t */
#define HASHTABLE_SIZE_T uint64_t       /* use 64 bit unsigned key type consistent with this file */
#include "hashtable.h"                  // Gustavsson's file was modified because the 64bit to 32bit key compression code produces 0 for some values
OEX hashtable_t genetable;              // genetable contains the hash table for genes : the hash table keys are the 64bit genes themselves
typedef struct genedata {               // value of keys stored for each gene encountered in simulation
    unsigned int popcount;              // initialized to 1
    short unsigned int firsttime;       // first time the gene was created from ancestor: initialized to 0
    short unsigned int recenttime;      // most recent time the gene was created by mutation from ancestor or randomly generated : initialized to 0
    short unsigned int lasttime;        // last time this gene was seen : enables genes to be mapped to their last epoch (including present)
    short int lastextinctiontime;       // this is initialized to -1, meaning no extinctions yet
    unsigned int activity;              // initialized to 0
    unsigned int nextinctions;          // initialized to 0
    unsigned int dummy;                 // pad to 64 bit word boundary
    uint64_t gene;                      // stored gene : note that two difft 64-bit genes may be stored at same location, so need check
    uint64_t firstancestor;             // this is initialized to a special gene seq (rootgene) not likely ever to occur for starting genes
} genedata;
#define rootgene 0xfedcba9876543210     /* initial special gene as root for genealogies */
// OEX const uint64_t INIT(rootgene,0xfedcba9876543210);// initial special gene as root for genealogies
OEX const uint64_t INIT(generepeat,0x0123456789abcdef);// special gene sequence to denote repeated genes found in genealogy
OEX genedata INIT(ginitdata,{1,0,0,0,-1,0,0,0,0ull,rootgene});// initialization data structure for gene data
OEX genedata *genedataptr;              // pointer to a genedata instance
OEX HASHTABLE_SIZE_T const* genotypes;  // pointer to stored hash table keys (which are the genotypes)
OEX genedata* geneitems;                // list of genedata structured items stored in hash table
OEX int genefnindices[1<<24];           // table of activities for functional gene indices calculated by genefnindex
//.......................................................................................................................................................
OEX hashtable_t quadtable;              // hash table for quad tree of spatial live state patterns
typedef struct quadnode {               // stored quadtree binary pattern nodes for population over time (currently only for analysis not computation)
    uint64_t hashkey;                   // hash table look up key for node : enables tree exploration more directly than construction from nw,ne,sw,se including collision avoidance
    uint64_t nw, ne, sw, se;            // constant keys to hashed quadnodes or 64-bit patterns : we terminate one level higher than Gosper & golly
    unsigned short int isnode;          // 1 if this is a node not a pattern
    unsigned short int size;            // side length of square image corresponding to quadtree pattern
    unsigned int activity;              // number of references to finding this node
    unsigned int pop1s;                 // 32 bit number of 1s in quadnode
    unsigned int firsttime;             // first time node was identified
    unsigned int lasttime;              // last time node was identified : used to determine if part of current timestep
    unsigned int topactivity;           // activity of pattern as top of connected component (last field : giving size of record an even number of 64bit words)
} quadnode;
OEX quadnode INIT(quadinit,{0ull,0ull,0ull,0ull,0ull,0,0,1,0,0,0,0});
OEX uint64_t qimage;                    // key for quadimage datastructure of entire image
OEX int INIT(quadcollisions,0);
OEX HASHTABLE_SIZE_T const* quadkeys;   // pointer to stored hash table keys (which are the quadkeys)
OEX quadnode* quaditems;                // list of quadnode structured items stored in hash table
typedef struct smallpatt {              // stored binary patterns for 4*4 subarrays or smaller (16bit)
    unsigned int topactivity;           // number of references to finding this pattern as top level of connected component
    unsigned int activity;              // number of references to finding this pattern
    unsigned int firsttime;             // first time pattern was identified
    unsigned int lasttime;              // last time pattern was identified
} smallpatt;
OEX smallpatt smallpatts[65536];
//.......................................................................................................................................................
OEX hashtable_t clonetable;             // hash table for clone ancestry
typedef struct clonedata {
    uint64_t birthid;                   // birth time and place of clone (upper 32 bits: time, lower 32 bits: place ij) plus 1 bit (N2) 0: ancesotr 1: no ancestor
    uint64_t parentid;                  // ancestor's id : i.e. birth time and place (upper 32 bits: time, lower 32 bits: place ij)
    uint64_t gene;                      // gene of this individual
    unsigned int popln;                 // number of live individuals in this clone (32 bit)
    unsigned int activity;              // activity of clone (32 bit)
    // unsigned int subclones;          // number of mutant clones stemming from clone (32 bit): not used currently
} clonedata;
OEX clonedata INIT(cinitdata,{0ull,0ull,0ull,1,1}); // initial values of clonedata: popln set to 1
OEX clonedata *clonedataptr;            // pointer to clonedata instance
OEX HASHTABLE_SIZE_T const* clones;     // pointer to stored hash table keys (which are the clones, live or with live descendants)
OEX clonedata* cloneitems;              // list of clonedata structured items stored in hash table
OEX const uint64_t INIT(rootclone,N2);  // single bit for mask enquiries for clones with root heritage
//-------------------------------------------------------------------------------------------------------------------------------------------------------
OEX int INIT(totsteps,0);               // total number of simulation steps
OEX int INIT(totdisp,0);                // total number of displayed steps
OEX int INIT(statcnts,0);               // total number of statistic timepts counted
OEX uint64_t codingmask;                // ncoding derived mask for ncoding bits
OEX int INIT(nhistG,0);                 // interval for collecting config histogram data : 0 no collection, nstatG collection with time
OEX int INIT(nstatG,0);                 // interval for collecting other statistical trace data : 0 no collection
OEX int INIT(genealogydepth,0);         // depth of genealogies in current population
OEX int INIT(genealogycoldepth,0);      // genes coloured by colour of ancestor at this depth in colorfunction=11
OEX int ngenealogydeep;                 // depth of genealogy
OEX int INIT(clonealogydepth,0);        // depth of clonealogies in current population
OEX int nclonealogydeep;                // depth of clonealogy
OEX int ambigsum;                       // temporary statistic of ambiguous resolution cases
OEX int nbshist[8];                     // global histogram to count neighbourhoods to check for assymetries
//.........................................................resource management NYI.......................................................................
OEX int INIT(rmax,1);                   // max number of resources per cell : 0 also turns off resource processing
OEX int INIT(rthresh,3);                // minimum resource number per neighborhood to allow birth process
OEX int INIT(rbirth,2);                 // minimum resource number per neighborhood to allow birth replenishment
//---------------------------------------------------------main arrays in simulation---------------------------------------------------------------------
OEX uint64_t *gol, *golg, *golb, *golr; // pointers to one plane of gol, golg, golb & golr arrays: live/dead, gene, cloneid (birth t,x,y), resource
OEX uint64_t *golgstats;                // pointer to 64 bit masks for different events during processing at one of the plane cycle locations
OEX uint64_t stashgol[N2];              // for stashing state and recovering with initfield = -1 or for component extraction and running
OEX uint64_t stashgolg[N2];             // for stashing genes and recovering with initfield = -1 or for component extraction and running
OEX uint64_t stashgolb[N2];             // for stashing birthid and recovering with initfield = -1 or for component extraction and running
OEX uint64_t stashgolr[N2];             // for stashing record of displacements and recovering with initfield = -1 or for component extraction and running
OEX uint64_t stashgolgstats[N2];        // for stashing stats and recovering with initfield = -1 or for component extraction and running
OEX uint64_t golmix[N2];                // array for packing configs
OEX uint64_t gene0;                     // uncoupled planes background gene, non zero for selection==16,17
OEX uint64_t selectedgene;              // gene currently selected interactively in graphics window
OEX uint64_t genegol[16];               // genes encoding for GoL in various LUT selection models indexed as selection-8
//---------------------------------------------------------additional variables set in simulation--------------------------------------------------------
OEX unsigned int canonical;             // current value of choice of canonical position repscheme bit 2 : needed globally in ...difft2-6 routines
OEX int INIT(quadrants,-1);             // integer choice of bit pair from repscheme/survivalmask/overwritemask for quadrant division of array (-1 none)
OEX int INIT(randominflux,0);           // 1,2 steady or intermittent random input of genes and gol states into the square central region defined by initfield
                                        // 3 steady deletion of genes at random rate rbackground
OEX int INIT(rbackground,0);            // integer background rate of random gene input per frame per site : rate is rbackground/32768 (nonzero overides influx)
                                        // the gene input depends on randominflux value: 2 GoL 1 random
OEX int INIT(vscrolling,0);             // whether to do vertical scrolling to track upwards growth (losing all states that fall off downward cliff)
OEX int INIT(vscrolly,0);               // cumulative extent of scrolling (mod N) : to allow clone birthids to be parsed correctly
OEX int INIT(last_scrolled,0);          // whether vscrolling applied on last time step (needed for correct glider detection)
OEX int INIT(ymax,2000);                // gene activity scale max for plotting : will be adjusted dynamically or by keys
OEX int INIT(ymaxq,2000);               // quad pattern activity scale max for plotting : will be adjusted dynamically or by keys
// double log2ymax = 25.0;              // activity scale max 2^25 = 33.5 * 10^6 : suffers from discrete steps at bottom, not used
OEX int activitymax;                    // max of activity in genealogical record of current population
OEX int INIT(activityfnlut,0);          // whether to lump activities of genes which have the same sequence at functionally masked in active positions of LUT
OEX int INIT(noveltyfilter,0);          // novelty filter for colorfunction 9 : if on (key "n"), darkens non-novel components (activity>1) in display
OEX int INIT(activity_size_colormode,0);// color by size for colorfunction 10 : if on (key "p")  1 log2 enclosing square size 2 use #pixels 3 use sqrt(#pixels)
OEX int xdisplay,INIT(ydisplay,-1);     // display x and y coordinates selected by mouse in python
OEX int INIT(info_transfer_h,0);        // whether to display histogram on glider information transfer counts (non zero)
OEX int INIT(it_nbhood,7);              // size of neighborhood for collecting glider characterization histogram : default 7x7 nbhood
OEX uint64_t gliderinfo[408];           // histogram of counts for glider detection by match quality in eight directions N E S W NE SE SW NW
//------------------------------------------------ arrays for time tracing, activity and genealogies ----------------------------------------------------
#define startarraysize 1024             /* starting array size (used when initializing second run) */
// OEX const int INIT(startarraysize,1024);// starting array size (used when initializing second run)
OEX int INIT(arraysize,startarraysize); // size of trace array (grows dynamically)
OEX int INIT(*livesites,NULL);          // dynamic array pointer for statistics of number of live sites over time
OEX int INIT(*genestats,NULL);          // dynamic array pointer for statistics of number of 4 genotype classes over time
OEX int INIT(*stepstats,NULL);          // dynamic array pointer for statistics of site update types over time
OEX int INIT(*configstats,NULL);        // dynamic array pointer for statistics of gol site configurations (x,y,t) offsets
OEX uint64_t poptrace[N2];              // scrolled trace of last N time points of population of N most frequent genes
OEX uint64_t acttrace[N2];              // scrolled trace of last N time points of activity of N most frequent genes
OEX uint64_t acttraceq[N2];             // scrolled trace of last N time points of activity of N most frequent quad patterns
OEX unsigned char acttraceqt[N2];       // type of entry in acttraceq : 1 quad, 0 smallpatt (<65536) i.e. corresponding to isnode
OEX uint64_t genealogytrace[N2];        // image trace of genealogies for N most frequently populated genes
OEX uint64_t clonealogytrace[N2];       // image trace of clonealogies for N most frequently populated clones
OEX int INIT(nbhist,-1);                // current block for trace
enum {nNhist = 20,                      // maximum number of older blocks for trace
      N2h = N2*nNhist,
      Nh = N*nNhist,
      log2N1 = log2N+1,
      N1 = N+1};
OEX uint64_t poptrace1[N2h];            // trace of first N*nNhist time points of population of N most frequent genes
OEX uint64_t acttrace1[N2h];            // trace of first N*nNhist time points of activity of N most frequent genes
OEX uint64_t acttraceq1[N2h];           // trace of first N*nNhist time points of activity of N most frequent quad patterns
OEX unsigned char acttraceqt1[N2h];     // type of entry in acttraceq : 1 quad, 0 smallpatt (<65536) i.e. corresponding to isnode
OEX uint64_t working[N2];               // working space array for calculating genealogies and doing neighbour bit packing
OEX int npopulation[N];                 // number of live sites (gol 1s) in last N time steps up to current population
OEX int npopulation1[Nh];               // number of live sites (gol 1s) in first N*nNhist time steps
OEX unsigned int nnovelcells[Nh];       // number of live sites that are part of novel components in first N*nNhist time steps
OEX int nspeciesgene,nallspecies;       // number of gene species in current population, and that have ever existed
OEX int nallclones;                     // number of clones stored in hash table
OEX int nspeciesquad,nallspeciesquad;   // number of quad species in current population, and that have ever existed
OEX int nspeciessmall,nallspeciessmall; // number of small pattern species now, and that have ever existed
OEX int histcumlogpattsize[log2N1];     // histogram of patterns binned on log scale according to power of two side enclosing square
OEX int histcumpixelssqrt[N1];          // histogram of patterns binned on an integer sqrt scale according to number of pixels
//------------------------------------------------ arrays for connected component labelling and tracking ------------------------------------------------
// enum { NLM = N2};                    // maximum number of discrete components possible N*N
// enum { NLC = N2<<2};                 // maximum number of connections N*N*4
OEX unsigned int label[N2];             // labels for pixels in connected component labelling
OEX unsigned int oldlabel[N2];          // previous time step labels for connected component labelling
OEX unsigned int labelcc[N2];           // label array to reassemble and display individual components
typedef struct equivrec {               // equivalence record
  unsigned int pt;                      // parent of record : integer index into eqv array at lower integer location
  unsigned int rank;                    // rank of equivalence class tree, may be used to accelerate union-find
  unsigned int size;                    // number of cells (lattice sites) with this label in current array
  unsigned int overlaps;                // reserved for future use, filling record overall size to single long int (64bit)
} equivrec;
OEX equivrec eqv[NLM];                  // equivalences between labels
typedef struct component {              // data structure for identified connected components in gol array
    short unsigned int N,S,W,E;         // rectangular bounds of rectangle containing the component
    short unsigned int log2n;           // label index of component, enclosing square size n=2^log2n
    short unsigned int patt;            // for small components of size 4x4 pixels or less, the image is encoded directly in the pattern patt
    short unsigned int lastrc;          // last occupied row/col used to trace components wrapping across the N-1 to 0 border
    short unsigned int dummy;           // not used: just to maintain long word boundaries
    unsigned int label;                 // label of component
    unsigned int pixels;                // number of pixels on
    uint64_t quad;                      // hashkey for quadtree node of subimage for component
    float gcolor;                       // inherited drifting color hue mixed from connected comp's (structure needs whole nr of 64-bit words for ndarray python comm.)
    unsigned int reserve;               // reserve position to align components on 64 bit boundary for efficiency
} component;
OEX component complist[NLM];            // current array of components
OEX component oldcomplist[NLM];         // old (previous time step) list of components
OEX int INIT(ncomponents,0);                // current number of components (= current number of labels)
OEX int INIT(oldncomponents,0);             // number of components in previous time step
typedef struct connection {             // connection of connected component at time t to another component at t-1
    unsigned int oldlab;                // label at time t-1
    unsigned int newlab;                // label at time t
    unsigned int next;                  // next connection index in list of backward connections for a given labelled component at time t (newlab)
    unsigned int nextf;                 // next connection index in list of forward connections for a given labelled component at time t-1 (oldlab)
    unsigned int overlap;               // sum of neighboring pixels (9-nbhd) between old and new component
    unsigned int reserve;               // unused element to preserve 64 bit alignment after shift from short unsigned int to unsigned int
    float woverlap;                     // weighted overlap out of all forward connections from old component i.e. overlap/sum(overlaps))
    float aoverlap;                     // weighted woverlap out of all backward connections from new component i.e. woverlap/sum(woverlaps))
} connection;
OEX connection connections[NLC];        // open memory reservoir of connection nodes to use in connection lists
OEX unsigned int connlists[NLM];        // entry points for connection list corresponding to each labelled component
OEX unsigned int connlistsf[NLM];       // entry points for forward connection list corresponding to each old labelled component
OEX unsigned int connlen[NLM];          // lengths of connection lists for connected components to previous time ie from t to t-1
OEX unsigned int connlenf[NLM];         // lengths of connection lists for connected components to current time i.e. from t-1 to t
OEX unsigned int connpref[NLM];         // preferred backward connected component at time t-1 for each component at time t
OEX unsigned int connpreff[NLM];        // preferred forward connected component at time t for each component at time t-1
OEX int INIT(connused,0);               // used connection nodes
OEX int INIT(gcolors,0);                // use inherited colors to color connected components
OEX int INIT(conn_genetic,1);           // 0: if gol state only used to determine connected components 1: if golg state used as well
//............................................... optimal linear assignment t-1 to t ....................................................................
// using maxmatch.c                     // arrays for matching components from one time step to the next, using Hopcroft Karp maxmatch algorithm in maxmatch.c
// #include "lapjv.h"                   // no longer using lapmod.c, modified from Tomas Kazmar python interfaced implementation of Jonker-Volgenant LAPMOD algorithm
OEX unsigned int iilap[NLM];            // indices of start of each variable length row in sparse cost matrix : first entry 0, last entry nclap
OEX unsigned int cclap[N2];             // sparse cost matrix for mapping connected components at t-1 to t (overlaps) containing nclap entries
OEX unsigned int recolor[NLM];          // recolor connected component based on colors of connected components and random drift
OEX unsigned int kklap[N2];             // column indices for successive entries in cost matrix
OEX unsigned int xlap[NLM];             // returned list of assignments: columns assigned to rows
OEX unsigned int ylap[NLM];             // returned list of assignments: rows assigned to columns
OEX unsigned int dist[NLM];             // distance along augmented paths for Hopcroft Karp matching algorithm: maxmatch
OEX unsigned int relabel[NLM];          // array to relabel connected components matching to be compatible with previous step
OEX unsigned int oldrelabel[NLM];       // old relabel array at previous time step
OEX unsigned int queue_array[NLM];      // array for queue used in Hopcroft Karp matching algorithm: maxmatch
OEX int  nlap;                          // number of connected components at t-1 entering into the assignment, i.e. n for LAPMOD
OEX int  nclap;                         // number of edges between connected comp's t-1 to t entering into the assignment, == no. of cost matrix entries
OEX int  nmatched;                      // number of matched old labels in current label set
//------------------------------------------------ planes and configuration offsets----------------------------------------------------------------------
OEX int INIT(offdx,0);                  // display chosen offsets for glider analysis with colorfunction 8
OEX int INIT(offdy,0);
OEX int INIT(offdt,0);
OEX int INIT(Noff,9);                   // number of offsets
OEX int **offsets;                      // array of offsets (2D + time) for planes
OEX int *histo;
OEX int numHisto;
// initialize planes:
#define maxPlane 8                      /* maximum number of planes allowed : values 2,4,8 allowed */
OEX int INIT(curPlane,0);               // current plane index
OEX int INIT(newPlane,1);               // new plane index
OEX int INIT(numPlane,maxPlane);        // number of planes must be power of 2 to allow efficient modulo plane
OEX uint64_t *planes[maxPlane];         // ring buffer planes of gol array states
OEX uint64_t *planesg[maxPlane];        // ring buffer planes of golg genes
OEX uint64_t *planesgs[maxPlane];       // ring buffer planes of golgstatus bits
OEX uint64_t *planesb[maxPlane];        // ring buffer planes of birth id for clone ids
OEX uint64_t *planesr[maxPlane];        // ring buffer planes of birth id for resources
OEX uint64_t plane0[N2];                // gol   0
OEX uint64_t plane1[N2];                // gol   1
OEX uint64_t planeg0[N2];               // golg  0
OEX uint64_t planeg1[N2];               // golg  1
OEX uint64_t planegs0[N2];              // golgs 0
OEX uint64_t planegs1[N2];              // golgs 1
OEX uint64_t planeb0[N2];               // golb  0
OEX uint64_t planeb1[N2];               // golb  1
OEX uint64_t planer0[N2];               // golr  0
OEX uint64_t planer1[N2];               // golr  1
#if maxPlane > 2
OEX uint64_t plane2[N2];                // gol   2
OEX uint64_t plane3[N2];                // gol   3
OEX uint64_t planeg2[N2];               // golg  2
OEX uint64_t planeg3[N2];               // golg  3
OEX uint64_t planegs2[N2];              // golgs 2
OEX uint64_t planegs3[N2];              // golgs 3
OEX uint64_t planeb2[N2];               // golb  2
OEX uint64_t planeb3[N2];               // golb  3
OEX uint64_t planer2[N2];               // golr  2
OEX uint64_t planer3[N2];               // golr  3
#endif
#if maxPlane > 4
OEX uint64_t plane4[N2];                // gol   4
OEX uint64_t plane5[N2];                // gol   5
OEX uint64_t plane6[N2];                // gol   6
OEX uint64_t plane7[N2];                // gol   7
OEX uint64_t planeg4[N2];               // golg  4
OEX uint64_t planeg5[N2];               // golg  5
OEX uint64_t planeg6[N2];               // golg  6
OEX uint64_t planeg7[N2];               // golg  7
OEX uint64_t planegs4[N2];              // golgs 4
OEX uint64_t planegs5[N2];              // golgs 5
OEX uint64_t planegs6[N2];              // golgs 6
OEX uint64_t planegs7[N2];              // golgs 7
OEX uint64_t planeb4[N2];               // golb  4
OEX uint64_t planeb5[N2];               // golb  5
OEX uint64_t planeb6[N2];               // golb  6
OEX uint64_t planeb7[N2];               // golb  7
OEX uint64_t planer4[N2];               // golr  4
OEX uint64_t planer5[N2];               // golr  5
OEX uint64_t planer6[N2];               // golr  6
OEX uint64_t planer7[N2];               // golr  7
#endif
//------------------------------------------------------- fast macros for pattern counting and random number generator ---------------------------------
OEX uint64_t randstate[2];              // State for xorshift pseudorandom number generation. The state must be seeded so that it is not zero
#define RAND128P(val) {                 /* Adapted from Wikipedia "Xorshift" rewritten here as inline macro & Vigna, Sebastiano, see acknowledgements */    \
    uint64_t x = randstate[0]; uint64_t y = randstate[1];                       \
    randstate[0] = y;    x ^= x << 23;  randstate[1] = x ^ y ^ (x >> 17) ^ (y >> 26);  \
    val = randstate[1] + y;}
OEX int INIT(ranseed,1234);
//.......................................................................................................................................................
OEX const uint64_t INIT(m1 ,0x5555555555555555); //binary: 0101...           Constants for Hamming distance macro POPCOUNT64C
OEX const uint64_t INIT(m2 ,0x3333333333333333); //binary: 00110011..
OEX const uint64_t INIT(m4 ,0x0f0f0f0f0f0f0f0f); //binary:  4 zeros,  4 ones ...
OEX const uint64_t INIT(h01,0x0101010101010101); //the sum of 256 to the power of 0,1,2,3...
OEX const uint64_t INIT(m10,0x1111111111111111); //binary: 00010001...       Additional constant for 4-bit byte distance macro POP4COUNT64C
#define POPCOUNT64C(x, val) {                  /* Wikipedia "Hamming Weight" popcount4c alg */  \
    uint64_t xxxx;                             /* define copy of x argument so that we do not change it */ \
    xxxx = x;                                  /* copy x argument */ \
    xxxx -= (xxxx >> 1) & m1;                  /* put count of each 2 bits into those 2 bits */ \
    xxxx = (xxxx & m2) + ((xxxx >> 2) & m2);   /* put count of each 4 bits into those 4 bits */ \
    xxxx = (xxxx + (xxxx >> 4)) & m4;          /* put count of each 8 bits into those 8 bits */ \
    val = (xxxx * h01) >> 56;}                 /* left 8 bits of x + (x<<8) + (x<<16) + (x<<24) + ... */
#define POP4COUNT64C(x, val) {                 /* Count number of 4-bit bytes which are non zero */  \
    uint64_t xxxx;                             /* define copy of x argument so that we do not change it */ \
    xxxx = x;                                  /* copy x argument */ \
    xxxx = (xxxx | (xxxx >> 1)) & m1;          /* put or of each 2 bits into right bit of pair */ \
    xxxx = (xxxx | (xxxx >> 2)) & m10;         /* put or of each quartet bits into right bit of quartet */ \
    xxxx = (xxxx + (xxxx >> 4)) & m4;          /* put count of each 8 bits into those 8 bits */ \
    val = (xxxx * h01) >> 56;}                 /* left 8 bits of x + (x<<8) + (x<<16) + (x<<24) + ... */
//.......................................................................................................................................................
#define PATTERN2_32(x, pat, found) {           /* find number of 2-bit aligned 2-bit copies of pattern pat in 32-bits of 64-bit integer x */ \
    const uint64_t r1_32 = 0x11111111ull;       \
    const uint64_t r3_32 = 0x33333333ull;       \
    const uint64_t r5_32 = 0x55555555ull;       \
    const uint64_t rf_32 = 0xffffffffull;       \
    uint64_t xxxx;                             /* define internal variable */ \
    xxxx=(uint64_t) pat;                       /* copy pat to ensure it is not assumed to be a variable that can be changed */ \
    xxxx|=xxxx<<2;                             /* doubles the pattern to the left */ \
    xxxx|=xxxx<<4;                             /* doubles the patterns to the left : now 4 copies */ \
    xxxx|=xxxx<<8;                             /* doubles the patterns to the left : now 8 copies */ \
    xxxx|=xxxx<<16;                            /* doubles the patterns to the left : now 16 copies */ \
    xxxx = (xxxx ^ (uint64_t) x) & rf_32;      /* xor x argument with 16 copies of pat */ \
    if (xxxx) {                                /* need to exclude the 16 match case as the fast count fails here */ \
        xxxx = ~xxxx&rf_32;                    /* invert difference map to yield identity map, 2 ones if match at a position */ \
        xxxx = (xxxx&(xxxx>>1))&r5_32;         /* convert 2-bit set of identity patterns to 1/0 decision at bit 0 of 2 bit pattern */ \
        xxxx = (xxxx+(xxxx>>2))&r3_32;         /* do first sum pairwise so that set of 8 2-bit sums spaced by 4-bits */ \
        found = ((xxxx&r1_32) * r1_32) >> 28;} /* found is returned as the number of patterns found at any of the 16 positions : <16 */ \
    else found = 16;}                          /* match at all positions */
//.......................................................................................................................................................
OEX const uint64_t INIT(r1,0x1111111111111111ull);
#define PATTERN4(x, pat, found) {              /* find number of 4-bit aligned 4-bit copies of pattern pat in 64-bit x */ \
    uint64_t xxxx;                             /* define internal variable */ \
    xxxx=(uint64_t) pat;                       /* copy pat to ensure it is not assumed to be a variable that can be changed */ \
    xxxx|=xxxx<<4;                             /* doubles the pattern to the left */ \
    xxxx|=xxxx<<8;                             /* doubles the patterns to the left : now 4 copies */ \
    xxxx|=xxxx<<16;                            /* doubles the patterns to the left : now 8 copies */ \
    xxxx|=xxxx<<32;                            /* doubles the patterns to the left : now 16 copies */ \
    xxxx ^= (uint64_t) x;                      /* xor x argument with 16 copies of pat */ \
    if (xxxx) {                                /* need to exclude the 16 match case as the fast count fails here */ \
        xxxx = ~xxxx;                          /* invert difference map to yield identity map, 4 ones if match at a position */ \
        xxxx &= (xxxx>>1)&(xxxx>>2)&(xxxx>>3); /* convert 4-bit set of identity patterns to 1/0 decision at bit 0 of 4 bit pattern */ \
        found = ((xxxx&r1) * r1) >> 60;}       /* found is returned as the number of patterns found at any of the 16 positions : <16 */ \
    else found = 16;}                          /* match at all positions */
//.......................................................................................................................................................
// const uint64_t h01 = 0x0101010101010101;    /* multiplicand defined already for POPCOUNT24: used to sum up all 8 8-bit bytes */
#define PATTERN8(x, pat, found) {              /* find number of 8-bit aligned 8-bit copies of pattern pat in 64-bit x */ \
    uint64_t xxxx;                             /* define internal variable */ \
    xxxx=(uint64_t) pat;                       /* copy pat to ensure it is not assumed to be a variable that can be changed */ \
    xxxx|=xxxx<<8;                             /* doubles the pattern to the left : now 2 copies */ \
    xxxx|=xxxx<<16;                            /* doubles the patterns to the left : now 4 copies */ \
    xxxx|=xxxx<<32;                            /* doubles the patterns to the left : now 8 copies */ \
    xxxx ^= (uint64_t) x;                      /* xor x argument with 8 copies of pat */ \
    xxxx = ~xxxx;                              /* invert difference map to yield identity map, 8 ones if match at a position */ \
    xxxx &= xxxx>>1;                           /* convert 8-bit set of identity patterns to 1/0 decision */ \
    xxxx &= xxxx>>2;                           /*    at bit 0 of 8 bit pattern */ \
    xxxx &= xxxx>>4;                           /*    in 3 steps */  \
    found = ((xxxx&h01) * h01) >> 56;}         /* found is returned as the number of patterns found at any of the 8 positions */
//.......................................................................................................................................................
#define FIRST1INDEX(v, c) {                    /* starting point 64bit from Sean Eron Anderson https://graphics.stanford.edu/~seander/bithacks.html#ZerosOnRightParallel */  \
    uint64_t mmmm,mmmq;                        /* calculates position of rightmost (lsb) 1 : arguments must be of types uint64_t and int respectivley */ \
    int cccc;                                  /* takes on successive integer values 32,16,84,2,1 */ \
    int tttt;                                  /* logical to integer variable true=one false=zero : if a 1 in v under mask mmmm */ \
    c=v?0:1;                                   /* c will contain count of number of zeros on right of last one, here if v is all zeros then start from 1 */ \
    mmmm=~0ull;                                /* initially all ones, this is the mask from previous stage in loop below */ \
    for (cccc=1<<5;cccc>0;cccc>>=1) {          /* loop over cccc goes 32,16,8,4,2,1 : the amount of shift used in mask construction */ \
        mmmq = mmmm;                           /* query mask mmmq is to be the mask used to query if a one is under it at this stage, start with old mask */ \
        mmmq &= mmmm^(mmmm<<cccc);             /* mmmq: 0xffffffff00000000, 0xffff0000ffff0000, 0xff00ff00ff00ff00, 0xf0f0f0f0f0f0f0f0, 0xccc..., 0xaaa...*/ \
        tttt = v&mmmq?0:1;                     /* tttt is zero if a one under the query mask, one otherwise */ \
        mmmm=mmmq^(tttt*mmmm);                 /* the new mask for next stage is the query mask if a one is under it, otherwise the other half of mmmm */ \
        c+=tttt*cccc;                          /* the right zero counter is incremented by the length of the current interval cccc if a one was not under mask */ \
    }                                          /* note that Anderson's algorithm was incorrect, see also profile comparison in standalone lsb64.c */ \
}                                              /* this macro calculates the LSB 1 (ie from the bottom) not the MSB 1 (ie from the top) that the integer log function finds. */
//.......................................................................................................................................................
#define DELTAXY(ij,x,y)  ((ij - (ij&Nmask) + (((ij+(x))&Nmask) + (y)*N)) & N2mask)  /* calculate offset (x,y) from index ij into periodic 2D array of size N */
//.......................................................................................................................................................
#define MEMBRANE (((ij>>log2N)==((N>>1)-(initfield>>1)-1) || ((ij>>log2N)==((N>>1)-(initfield>>1)-2))) && (ij & 0x1) ? 2 : 1) /* formula for membrane of death */

//------------------------------------------------------ function prototypes -----------------------------------------------------------------------------
// ......................................................... mathext.c ...................................................................................
extern INLINE int integerSqrt(int n);
extern INLINE unsigned int log2r(unsigned int v);
extern INLINE unsigned int log2lower(unsigned int v);
extern INLINE unsigned int log2upper(unsigned int v);
extern INLINE unsigned int sqrtupper(unsigned int v);
extern INLINE uint64_t randprob(unsigned int uprob, unsigned int randnr);
// .......................................................... pack.c ......................................................................................
extern INLINE void pack012neighbors(uint64_t gol[],uint64_t golp[]);
extern INLINE void pack0123neighbors(uint64_t gol[],uint64_t golp[]);
extern INLINE void pack49neighbors(uint64_t gol[],uint64_t golp[], int nbhood);
extern INLINE short unsigned int pack16neighbors(uint64_t wgol[], int log2n);
extern INLINE void unpack16neighbors(const short unsigned golpw, unsigned int labelimg[],const unsigned int label,const int offset);
extern INLINE int log2size(const short unsigned int golpw);
extern INLINE void pack64neighbors(uint64_t gol[],uint64_t golp[],int log2n);
extern INLINE void unpack64neighbors(const uint64_t golpw, unsigned int labelimg[], const unsigned int label, const int offset);
extern INLINE void compare_neighbors(uint64_t a[],uint64_t b[], int dx, int dy);
extern INLINE void compare_all_neighbors(uint64_t a[],uint64_t b[]);
extern INLINE void packandcompare(uint64_t newgol[],uint64_t working[],uint64_t golmix[]);
extern INLINE void golr_digest (uint64_t golr, unsigned int *mismatchmin, unsigned int *mismatchmax, unsigned int *period, int *pershx, int *pershy);
// ....................................................... hashtables.c ....................................................................................
extern INLINE void hashaddgene(int ij,uint64_t gene,uint64_t ancestor,uint64_t *golb,uint64_t parentid,uint64_t mutation);
extern INLINE void hashdeletegene(uint64_t gene,uint64_t birthid,const char errorformat[]);
extern INLINE void hashgeneextinction(uint64_t gene,const char errorformat[]);
extern INLINE void hashgeneactivity(uint64_t gene, const char errorformat[]);
extern INLINE void hashaddclone(uint64_t birthid, uint64_t parentid, uint64_t gene);
extern INLINE void hashdeletefromclone(uint64_t birthid);
extern INLINE void hashcloneactivity(uint64_t birthid, const char errorformat[]);
extern INLINE uint16_t rotate16(uint16_t patt);
extern INLINE uint64_t rotate64(uint64_t patt);
extern INLINE void rotate4x64(uint64_t *nw, uint64_t *ne, uint64_t *sw, uint64_t *se);
extern INLINE void rotatequad(uint64_t *nw, uint64_t *ne, uint64_t *sw, uint64_t *se);
extern INLINE quadnode * hash_patt8_store(const uint64_t h, const uint64_t patt);
extern INLINE void hash_patt4_find(const short unsigned int patt);
extern INLINE quadnode * hash_patt8_find(const uint64_t patt);
extern INLINE uint64_t patt_hash(const uint64_t a, const uint64_t b, const uint64_t c, const uint64_t d);
extern INLINE uint64_t node_hash(const uint64_t a, const uint64_t b, const uint64_t c, const uint64_t d);
extern INLINE quadnode * hash_patt16_store(const uint64_t h, const uint64_t nw, const uint64_t ne, const uint64_t sw, const uint64_t se);
extern INLINE quadnode * hash_patt16_find(const uint64_t nw, const uint64_t ne, const uint64_t sw, const uint64_t se);
extern INLINE quadnode * hash_node_store(uint64_t h, uint64_t nw, uint64_t ne, uint64_t sw, uint64_t se);
extern INLINE quadnode * hash_node_find(const uint64_t nw, const uint64_t ne, const uint64_t sw, const uint64_t se);
extern uint64_t quadimage(uint64_t gol[], short unsigned int *patt, int log2n);
extern int labelimage(uint64_t hashkeypatt, unsigned int labelimg[], unsigned int label, int offset);
// ......................................................... label.c ......................................................................................
extern INLINE unsigned int lab_union(equivrec eqv[], unsigned int i, unsigned int j);
extern INLINE unsigned int label_cell(int ij,unsigned int *nlabel);
extern INLINE unsigned int label_cell_genetic(int ij,unsigned int *nlabel,uint64_t golg[]);
extern void checklabels(equivrec eqv[],unsigned int *nlabel);
extern void flattenlabels(equivrec eqv[],unsigned int *nlabel);
extern unsigned int label_components(uint64_t gol[],uint64_t golg[]);
extern unsigned int extract_components(uint64_t gol[],uint64_t golg[]);
extern int novelcells(void);
//...................................................... mapping of connected components .................................................................
int maxmatch(int m, unsigned int kk[], unsigned int ii[], unsigned int pairU[], unsigned int pairV[], unsigned int dist[]);
// ......................................................... display.c ....................................................................................
extern INLINE void setcolor(unsigned int *color,int n);
extern INLINE unsigned int labelcolor( unsigned int label);
extern INLINE unsigned int rgba( float hue);
extern INLINE float mixcolor( unsigned int label, uint64_t rand);
extern void delay(int milliseconds);
extern void printxy (uint64_t gol[],uint64_t golg[]);
extern void printscreen (uint64_t gol[], uint64_t golg[]);
extern void colorgenes( int cgolg[], int NN2, int colorfunction, int winnr, int nfrstep);
// ....................................................... select_nbs.c ....................................................................................
extern INLINE void selectone_of_2(int s, uint64_t nb2i, int nb[], uint64_t golg[], uint64_t golb[],uint64_t * birth, uint64_t *newgene, uint64_t *parentid, unsigned int *kch);
extern INLINE int selectone_of_s(unsigned int *kch, int s, uint64_t nb1i, int nb[], uint64_t golg[], uint64_t golb[], uint64_t golr[], uint64_t *birth, uint64_t *newgene, uint64_t *parentid, uint64_t *nbmask, int ij);
extern INLINE void selectone_nbs(int s, uint64_t nb2i, int nb[], uint64_t gol[], uint64_t golg[],uint64_t golb[], uint64_t * birth, uint64_t *newgene, uint64_t *parentid, unsigned int *kch);
extern INLINE unsigned int selectdifft0(uint64_t nbmask, int *crot, int *kodd);
extern INLINE unsigned int selectdifft1(uint64_t nbmask, int *crot, int *kodd);
extern INLINE unsigned int selectdifft2(uint64_t nbmask, int *crot, int *kodd);
extern INLINE unsigned int selectdifft3(uint64_t nbmask, int *crot, int *kodd);
extern INLINE unsigned int selectdifft4(uint64_t nbmask, int *crot, int *kodd);
extern INLINE unsigned int selectdifft5(uint64_t nbmask, int *crot, int *kodd);
extern INLINE unsigned int selectdifft6(uint64_t nbmask, int *crot, int *kodd);
extern INLINE unsigned int selectdifft7(uint64_t nbmask, int *crot, int *kodd);
extern INLINE unsigned int selectdifft(int sum, uint64_t nbmask, int *crot, int *kodd, int *nsame);
extern INLINE void analyze_nbs(int ij, uint64_t gol[], int nnb[], uint64_t *nnb1i, int *ns);
extern INLINE uint64_t disambiguate(unsigned int *kchx, uint64_t nb1i, int nb[],  uint64_t gol[], uint64_t golg[], uint64_t golb[], int nsame, uint64_t *birth, uint64_t *parentid, uint64_t *ancestor, int ij);
// ....................................................... spatial_control.c ..............................................................................
extern void v_scroll(uint64_t newgol[],uint64_t newgolg[],uint64_t newgolb[],uint64_t newgolr[]);
extern void random_influx(uint64_t newgol[],uint64_t newgolg[],uint64_t newgolb[],uint64_t newgolr[]);
//............................................ many other routines not called externally from C, see separate files .......................................
//---------------------------------------------------------------------------------------------------------------------------------------------------------
#endif /* genelife_h */
