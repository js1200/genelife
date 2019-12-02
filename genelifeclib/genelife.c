//
// "subgenelife.h"
// project genelife
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
//
//----------------------------------------------------- list of subroutines -----------------------------------------------------------------------------
//......................................................    macros    ...................................................................................
// RAND128P(val)               64-bit random number algorithm adapted from Vigna Sebastiano, see acknowledgements
// POPCOUNT64C(x, val)         Wikipedia "Hamming Weight" popcount4c alg
// POP4COUNT64C(x, val)        Count number of 4-bit bytes which are non zero
// PATTERN2_32(x, pat, found)  find number of 2-bit aligned 2-bit copies of pattern pat in 32-bits of 64-bit integer x
// PATTERN4(x, pat, found)     find number of 4-bit aligned 4-bit copies of pattern pat in 64-bit x
// PATTERN8(x, pat, found)     find number of 8-bit aligned 8-bit copies of pattern pat in 64-bit x
// FIRST1INDEX(v, c)           index of first one in 64-bit binary pattern
// DELTAXY(ij,x,y)             calculate offset (x,y) from index ij into periodic 2D array of size N
// MEMBRANE                    formula for membrane of death, optional places in array where cell death is forced
//......................................................  fast integer processing   .....................................................................
// integerSqrt          direct bit processing algorithm to implement integer sqrt (largest integer smaller than sqrt) : but floating point sqrt is faster
// log2r                fast integer logarithm working only for arguments which are powers of 2 (not used but slightly faster than log2a when applicable)
// log2lower            fast integer logarithm working for all integers : largest integer smaller than or equal to logarithm base 2 of argument
// log2upper            fast integer logarithm working for all integers : smallest integer larger than or equal to logarithm base 2 of argument
// sqrtupper            fast integer sqrt working for all integers : smallest integer larger than or equal to sqrt of argument
// randprob             random event with probability determined by a 32 bit unsigned integer iprob as iprob / 2^32 using RAND128 and uint64_t
//......................................................  color functions   .............................................................................
// set_color            function to assign rainbow colors
// label_color          map label (quasi-uniformly) into larger unsigned colour space
// rgba                 converts hue (0..1) to rgb+alpha
// mix_color            mix colors from overlapping components, add random drift of colour
// delay                time delay in ms for graphics
// printxy              simple terminal print of array
// printscreen          terminal screen print of array on xterm, moving cursor with esc codes */
// golr_digest          digest information in golr (displacement record), extracting min and max mismatch, period, and x,y displacement of period
// colorgenes           colour display of genes in one of 12 colorfunction modes, including activities, pattern analysis, genealogies and glider detection
//......................................................  selection of genes for birth  .................................................................
// selectone_of_2       select one (or none) of two genes based on selection model parameter selection :  returns birth and newgene
// selectone_of_s       select one (or none) of s genes based on selection model parameter selection :  returns birth and newgene
// selectone_nbs        select one of two genes based on pattern of their live 2nd shell neighbours and their genetic encoding
// selectdifft1         select the gene at the single active neighbour position : algorithm could be optimized
// selectdifft2         select the right or left of two genes bunched with least number of empty genes between them
// selectdifft3         select the unique most different (by symmetry) of three live neighbours or first one in canonical rotation
// selectdifft4         select the most central (left) of four live neighbours or first one in canonical rotation
// selectdifft5         select the most central (left) of five live neighbours or first one in canonical rotation
// selectdifft6         select the most central (left) of six live neighbours or first one in canonical rotation
// selectdifft7         select the most central (left) of seven live neighbours or first one in canonical rotation
// selectdifft          select the most central (left) of sum live neighbours or first one in canonical rotation : calls 1 of selectdifft1-7
// analyze_nbs          analyze neighbour states, calculating s and packed set of live indices
// disambiguate         disambiguate the cases where the canonical rotation does not uniquely identify a pattern start point : 1 of 8 methods
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
//......................................................  neighborhood processing  ......................................................................
// pack012neighbors     pack all up to 2nd neighbours in single word
// pack0123neighbors    pack all up to 3rd neighbours in single word
// pack49neighbors      fast routine to pack all up to 1st,2nd,3rd neighbours in single word : order of bits dictated by hierarchical assembly
// pack16neighbors      pack 4x4 blocks in single uint64_t word using 16 bits
// unpack16neighbors    unpack 16 bit word to 4x4 block at offset in full array of labels, marking with chosen label
// log2size             return log2 of linear size of pattern in integer power of 2 for small patterns 0-65535
// pack64neighbors      pack 8x8 blocks in single uint64_t (long) word with 64 bits
// unpack64neighbors    unpack 64 bit word to 8x8 block at offset in full array of labels, marking with chosen label
// compare_neighbors    compare packed pack neighbours with one given x,y shift of arbitrary size
// compare_all_neighbors compare packed pack neighbours with all nearest neighbour x,y shifts
// packandcompare       pack and compare either all 1-shifted 3-neighbourhoods with t=-1 or chosen (dx,dy,dt) 3-neighbourhoods
// golr_digest          digest information in golr (displacement record), extracting min and max mismatch, period, and x,y displacement of period
//......................................................  fast component labelling  ......................................................................
// lab_union            disjoint rank union of equivalence classes returning common root
// label_cell           label a cell (site) in the cellular automata with first pass label of connected component
// label_cell_genetic   label a cell (site) for connected component analysis taking differences in genes into account
// checklabels          check that the label tree consistently points to labels of lower values as we go via parents to root
// flattenlabels        flatten label tree so that each label points to its unique root
// label_components     do two-pass fast component labelling with 8-neighbour using Suzuki decision tree, rank union and periodic BCs, connect t-1 labels with t
// extract_components   extract labelled components to list of subimages embedded in square of side 2^n, each stored in a quadtree hash table
//...................................................... mapping of connected components .................................................................
// queue_push           implements one queue using an array, push
// queue_pop            implements one queue using an array, push
// bfs                  returns true if there is an augmenting path in graph represented in sparse matrix form using cclap,kklap and iilap
// hopcroftKarp         returns size of maximum matching using Hopcroft Karp algorithm
// maxmatch             preprocesses array to identify max vertex on right and then calls hopcroftKarp to get max match: only routine called from genelife
//........................................................  simulation update for different symmetries  ..................................................
// update_23            update gol, golg, golgstats for a single synchronous time step : for selection 0-7 with fixed GoL rule departures in repscheme
// update_lut_sum       update version for gene encoding look up table for totalistic survival and birth (disallowing 0 live neighbour entries) sel 8,9
// update_lut_dist      update version for gene encoding look up table for survival and birth based on corner & edge sums (2*19 states disallowing s=0,1,7,8): sel 10,11
// update_lut_canon_rot update version for gene encoding look up table for canonical rotation survival and birth (2*32 states, disallowing 0,1,7,8 entries) : sel 12,13
// update_lut_2D_sym    update version all different configurations under the standard 2D 4-rotation and 4-reflection symmetries are distinguished: sel 14,15
// genelife_update      master routine to call specific model update, collect statistics if required and rotate planes
//........................................................  initialization and file IO  ...................................................................
// initialize_planes    initialize periodic sequence of planes to record rolling time window of up to maxPlanes time points (≤8)
// readFile             read file of gol/golg array (32x32) data
// writeFile            write file of gol/golg array (32x32) data
// testmacros           test macros used to accelerate processing (usually not called): FIRST1INDEX, PATTERN4, PATTERN8
// initialize           initialize simulation parameters and arrays
//............................................................ spatial control ............................................................................
// vscroll              vertical scroll of array in addition to dynamics with wrapping block at top/bottom barrier
// random_influx        random influx of live genes in initialization containing square, in continuous or intermittent modes, with optiional boundary feathering
//......................................................... stash array state .............................................................................
// savegols             save current arrays to file
// retrievegols         retrieve arrays to current from file
// stash                stash current gol,golg, golb, golr, golgstats in stashgol, stshgolg, stashgolb, stashgolr, stashgolgstats
// unstash              retrieve current gol,golg, golb, golr, golgstats from stashed values
// label2stash          stash current gol,golg, golb, golr, golgstats from selected labelled component (either cumulatively or individually)
//.........................................................  set from python driver  ......................................................................
// set_colorfunction    set color function integer from GUI for use in patterning and coloring display
// setget_act_ymax      set activity ymax for scaling of gene activity plot
// setget_act_ymaxq     set activity ymax for scaling of quad activity plot
// set_selectedgene     set selected gene for highlighting from current mouse selection in graphics window
// set_offsets          set offsets for detection of glider structures in display for color function 8
// set_quadrant         set the pair of bits in repscheme (or survivalmask or overwritemask) used for quadrant variation 0-6
// set_randominflux     change the randominflux activation for rbackground if nonzero or continual updating of init field with random states and genes : 2,1,0
// set_rbackground      set the backround random live gene input rate per frame and site to rbackground/32768
// set_repscheme_bits   set the two of the repscheme (or survivalmask or overwritemask) bits corresponding to the selected quadrant
// set_repscheme        set repscheme from python
// set_rulemod          set rulemod from python
// set_surviveover      set the two masks for survival and overwrite from python (survivalmask, overwritemask)
// set_vscrolling       set vertical scrolling to track fronts of growth in vertical upwards direction
// set_noveltyfilter    set novelty filter for darkening already encountered components in connected component display (colorfunction 9)
// set_activity_size_colormode set colormode by size for colorfunction 10 : 0 by ID  1 log2 enclosing square size 2 use #pixels 3 use sqrt(#pixels)
// set_gcolors          set connected component colors as inherited colors from colliding connected components with random drift
// set_seed             set random number seed
// set_nbhist           set nbhist N-block of time points for trace from GUI for use in activity and population display traces
// set_genealogycoldepth set genealogycoldepth for colorfunction=11 display
// set_ancestortype     set ancestortype for genealogy display and return of first (0), clonal (1) or first & clonal in 2 windows (2)
// set_info_transfer_h  set information transfer histogram display value (0,1) from python
// set_activityfnlut    set collection of functional activity statistics corresponding to functional aggregate of genes by non-neutral bits
// set_colorupdate1     control update of colorgenes and regular print statements via flag colorupdate1
// set_colorfunction2   choice of colorfunction for window 2
//..........................................................  get to python driver  .....................................................................
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
// get_genealogies      // get current population genealogies (currently near end of code)
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
//..........................................................  comparison functions  ....................................................................
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
//..........................................................  gene and pattern analysis of dynamics .....................................................
// countconfigs         count the configs with python specified offsets in (x,y,t)
// tracestats           record the current stats in time trace
// countspecies1        count genes with gene array specified as input parameters
// countspecies         count different genes with genes specified at current time point
// countspecieshash     count different genes in current population from record of all species that have existed
// totalpoptrace        calculates and returns current population size and store in scrolling population trace array npopulation
// genefnindex          calculate index based on bits masked in from survival and birth masks only
//..........................................................  activity analysis of dynamics ..............................................................
// activitieshash       calculate array of current gene activities and update acttrace array of genes in activity plot format
// activitieshashquad   calculate array of current quad activities and update acttraceq array of patterns in activity plot format
//..........................................................  genealogy and clonal analysis of dynamics ..................................................
// get_genealogies      calculate and retrieve or display genealogies, depending on size of array passed (0: display only, >0: retrieve only)
// clonealogies         calculate and display clonealogies: genealogies of clones
//-------------------------------------------------------------------------------------------------------------------------------------------------------
#define OEX                             /* this is main file where global variables and constants are allocated, so no extern statement */
#include "genelife.h"
   
//-------------------------------------------------------------------------------------------------------------------------------------------------------
