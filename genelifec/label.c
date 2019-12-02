//
//  label.c
//  project genelife
//
//  connected component labelling and mapping between time steps
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
// lab_union            disjoint rank union of equivalence classes returning common root
// label_cell           label a cell (site) in the cellular automata with first pass label of connected component
// label_cell_genetic   label a cell (site) for connected component analysis taking differences in genes into account
// checklabels          check that the label tree consistently points to labels of lower values as we go via parents to root
// flattenlabels        flatten label tree so that each label points to its unique root
// label_components     do two-pass fast component labelling with 8-neighbour using Suzuki decision tree, rank union and periodic BCs, connect t-1 labels with t
// extract_components   extract labelled components to list of subimages embedded in square of side 2^n, each stored in a quadtree hash table
//----------------------------------------------------------- connected component labelling ------------------------------------------------------------
extern INLINE unsigned int lab_union(equivrec eqv[], unsigned int i, unsigned int j) {
// Combine two trees containing node i and j. Union by rank with halving - see https://en.wikipedia.org/wiki/Disjoint-set_data_structure
// Return the root of union tree in equivalence relation's disjoint union-set forest
    unsigned int root = i;
    while (eqv[root].pt<root) {
        eqv[root].pt = eqv[eqv[root].pt].pt;       // halving algorithm to compress path on the fly
        root = eqv[root].pt;
    }
    if (i != j) {
        unsigned int rooti = root;
        unsigned int rootj;
        root = j;
        while (eqv[root].pt<root) {
            eqv[root].pt = eqv[eqv[root].pt].pt;   // halving algorithm to compress paths on the fly
            root = eqv[root].pt;
        }
        if (root == rooti) return root;
        rootj=root;
        // if (eqv[rooti].rank < eqv[rootj].rank) {                 // swap rooti, rootj so that rooti has larger rank : disrupts ordering!
        if (rooti > rootj) {                        // swap rooti, rootj so that rooti is lower of two
            rootj=rooti;
            rooti=root;
        }
        eqv[rootj].pt = rooti;                      // merge rootj into rooti
        if (eqv[rooti].rank == eqv[rootj].rank) {
            eqv[rooti].rank++;
        }
    }
    return root;
}
//.......................................................................................................................................................
extern INLINE unsigned int label_cell(int ij,unsigned int *nlabel) {
    unsigned int clabel,clabel1,labelij,labelijold;
    labelijold=label[ij];
    clabel=label[DELTAXY(ij,0,-1)];                         // deltaxy takes periodic BCs into account (b in Suzuki Decision Tree, 2010)
    if(clabel) labelij=eqv[clabel].pt;         // copy b
    else {                                                  // N (=b) unlabelled
        clabel=label[DELTAXY(ij,-1,0)];                     // W (d in Suzuki DT)
        if(clabel) {
            clabel1=label[DELTAXY(ij,1,-1)];                // NE (c in Suzuki DT)
            if(clabel1) {
                labelij=lab_union(eqv,clabel1,clabel);      // resolve c,d
            }
            else {
                labelij=eqv[clabel].pt;        // copy d
            }
        }
        else {
            clabel=label[DELTAXY(ij,1,-1)];                 // NE (c in Suzuki DT)
            if(clabel) {
                clabel1=label[DELTAXY(ij,-1,-1)];           // NW (a in Suzuki DT)
                if(clabel1) {
                    labelij=lab_union(eqv,clabel,clabel1);  // resolve c,a
                }
                else {
                    labelij=eqv[clabel].pt;    // copy c
                }
            }
            else {
                clabel1=label[DELTAXY(ij,-1,-1)];           // NW (a in Suzuki DT)
                if(clabel1) {
                    labelij=eqv[clabel1].pt;   // copy a
                }
                else {
                    clabel=++*nlabel;                       // new label if a,b,c,d all without labels
                    eqv[clabel].pt=clabel;
                    labelij=clabel;
                    eqv[clabel].rank=0;
                }
            }
        }
    }
    if (labelijold) labelij= lab_union(eqv,labelij,labelijold); // resolve labelij,labelijold needed for periodic BCs wrap around
    return labelij;
}
//.......................................................................................................................................................
extern INLINE unsigned int label_cell_genetic(int ij,unsigned int *nlabel,uint64_t golg[]) {
    unsigned int clabel0,clabel1,labelij,labelijold;
    int ij0,ij1;
    uint64_t gene;
    
    labelijold=label[ij];
    gene = golg[ij];
    clabel0=label[ij0=DELTAXY(ij,0,-1)];                         // deltaxy takes periodic BCs into account (b in Suzuki Decision Tree, 2010)
    if(clabel0&&(golg[ij0]==gene))
        labelij=eqv[clabel0].pt;                // copy b
    else {                                                       // N (=b) unlabelled
        clabel0=label[ij0=DELTAXY(ij,-1,0)];                     // W (d in Suzuki DT)
        if(clabel0&&(golg[ij0]==gene)) {
            clabel1=label[ij1=DELTAXY(ij,1,-1)];                 // NE (c in Suzuki DT)
            if(clabel1&&(golg[ij1]==gene)) {
                labelij=lab_union(eqv,clabel1,clabel0);          // resolve c,d
            }
            else {
                labelij=eqv[clabel0].pt;        // copy d
            }
        }
        else {
            clabel0=label[ij0=DELTAXY(ij,1,-1)];                 // NE (c in Suzuki DT)
            if(clabel0&&(golg[ij0]==gene)) {
                clabel1=label[ij1=DELTAXY(ij,-1,-1)];            // NW (a in Suzuki DT)
                if(clabel1&&(golg[ij1]==gene)) {
                    labelij=lab_union(eqv,clabel0,clabel1);      // resolve c,a
                }
                else {
                    labelij=eqv[clabel0].pt;    // copy c
                }
            }
            else {
                clabel1=label[ij1=DELTAXY(ij,-1,-1)];            // NW (a in Suzuki DT)
                if(clabel1&&(golg[ij1]==gene)) {
                    labelij=eqv[clabel1].pt;   // copy a
                }
                else {
                    clabel0=++*nlabel;                           // new label if a,b,c,d all without labels
                    eqv[clabel0].pt=clabel0;
                    labelij=clabel0;
                    eqv[clabel0].rank=0;
                }
            }
        }
    }
    if (labelijold) labelij= lab_union(eqv,labelij,labelijold); // resolve labelij,labelijold needed for periodic BCs wrap around
    return labelij;
}
//.......................................................................................................................................................
void checklabels(equivrec eqv[],unsigned int *nlabel) {
    unsigned int xlabel = 0;
    unsigned int i;
    for (i = 1; i < *nlabel + 1; i++) {
        if (eqv[i].pt > i) {
            fprintf(stderr,"Error in label equivalencies at t=%d for i=%d, parent %d\n",totsteps,i,eqv[i].pt);
            xlabel++;
        }
    }
    if (xlabel) fprintf(stderr,"Error for %d cases\n",xlabel);
}
//.......................................................................................................................................................
void flattenlabels(equivrec eqv[],unsigned int *nlabel) {
// Flatten the Union-Find tree and relabel the components.
    unsigned int xlabel = 1;
    unsigned int i;
    for (i = 1; i <= *nlabel; i++) {
        if (eqv[i].pt < i) {
            eqv[i].pt = eqv[eqv[i].pt].pt;
        }
        else {
            eqv[i].pt = xlabel++;
        }
    }
    *nlabel = xlabel-1;
}
//.......................................................................................................................................................
unsigned int label_components(uint64_t gol[],uint64_t golg[]) {
// fast component labelling, updating global equivalence table equiv and placing labels in global array label, return number of labels
    int i,ij,sum,sumoverlap,maxoverlap,overallmaxoverlap;   //   abc   scan mask to calculate label at position e
    unsigned int nlabel = 0;                                //   de
    unsigned int oldnlabel;                                 //   with lapmod need also: ret
    unsigned int lab,conn,connprev,connf,maxlabel;
    float rsumoverlap;
    int nunique,nzconnect,nremconnect;
    int dx[9]={0,-1,0,1,1,1,0,-1,-1};
    int dy[9]={0,-1,-1,-1,0,1,1,1,0};
    static int connectout = 0;                              // whether to print lists of connected component mappings t-1 t
    void testflow(void);
    extern int maxmatch(int m, unsigned int kk[], unsigned int ii[], unsigned int xlap[], unsigned int ylap[], unsigned int dist[]);
    static int first = 1;
    if(!first) {
        oldnlabel = oldncomponents = ncomponents;
        for(ij=0;ij<N2;ij++) oldlabel[ij]=label[ij];
        for(i=0;i<NLM;i++) oldrelabel[i]=0;
        for(i=0;i<ncomponents;i++) oldrelabel[i]=relabel[i];
    }
    else {
        for(ij=0;ij<N2;ij++) oldlabel[ij]=0;
        oldnlabel = oldncomponents = 0;
        for(i=0;i<NLM;i++) oldrelabel[i]=0;
        first = 0;
    }
    for(ij=0;ij<N2;ij++) label[ij]=0;
    for(i=0;i<(NLM);i++) eqv[i].pt=0;
    for(i=0;i<(NLM);i++) eqv[i].rank=0;
    for(i=0;i<(NLM);i++) eqv[i].size=0;
    for(i=0;i<(NLM);i++) xlap[i]=ylap[i]=0;

    if(conn_genetic) {                                          // taking genetic differences into account : separate label if different gene
            for(ij=0;ij<N2;ij++) {                              // do first pass main array
            if (gol[ij]) {
                label[ij]=label_cell_genetic(ij,&nlabel,golg);
            }
        }
        for(ij=0;ij<N;ij++) {                                   // redo first row sites to take periodic wrap around into account
            if (gol[ij]) {
                label[ij]=label_cell_genetic(ij,&nlabel,golg);
            }
        }
        for(ij=N;ij<N2;ij+=N) {                                 // redo first column sites to take periodic wrap around into account
            if (gol[ij]) {
                label[ij]=label_cell_genetic(ij,&nlabel,golg);
            }
        }
    }
    else {                                                      // ignoring genetic differences,  looking only at live/empty state of gol array for connected components
        for(ij=0;ij<N2;ij++) {                                  // do first pass main array
            if (gol[ij]) {
                label[ij]=label_cell(ij,&nlabel);
            }
        }
        for(ij=0;ij<N;ij++) {                                   // redo first row sites to take periodic wrap around into account
            if (gol[ij]) {
                label[ij]=label_cell(ij,&nlabel);
            }
        }
        for(ij=N;ij<N2;ij+=N) {                                 // redo first column sites to take periodic wrap around into account
            if (gol[ij]) {
                label[ij]=label_cell(ij,&nlabel);
            }
        }
    }
    // checklabels(eqv,&nlabel);                            // check labels point to lower values in equivalence table
    flattenlabels(eqv,&nlabel);                             // single pass flatten of equivalence table

    for(ij=0;ij<N2;ij++) {                                  // do second pass of main array assigning final root labels and calculating component sizes in pixels
        if (gol[ij]) {
            label[ij]=eqv[label[ij]].pt;
            eqv[label[ij]].size++;
       }
    }

    for(i=0;i<NLM;i++) connlists[i]=0;                      // intialize all connection lists to zero
    for(i=0;i<NLM;i++) connlen[i]=0;
    for(i=0;i<NLM;i++) connlistsf[i]=0;
    for(i=0;i<NLM;i++) connlenf[i]=0;
    for(ij=0;ij<NLC;ij++) {                                  // these are the open memory reserve of connection elements to be linked, better to use memset perhaps
        connections[ij].next=connections[ij].nextf=0;
        connections[ij].oldlab=connections[ij].newlab=0;
        connections[ij].overlap=connections[ij].reserve=0;
    }
    connused=0;

    for(ij=0;ij<N2;ij++) {                                  // build up backwards connections for each label[ij] via nbs in t-1 frame, insert in increasing label nr
        if(label[ij]) {
            for(int k=0;k<9;k++) {
                if((lab=oldlabel[DELTAXY(ij,dx[k],dy[k])])) {
                    conn=connlists[label[ij]];
                    connprev=0;
                    while(conn && (connections[conn].oldlab<lab)) {
                        connprev = conn;
                        conn=connections[conn].next;
                    }
                    if (conn) {                             // connections[conn].oldlab>=lab : if oldlab>lab insert & increment overlap, else just increment overlap
                        if(connections[conn].oldlab>lab) {  // insert new node
                            if(connused>=NLC) fprintf(stderr,"Error, out of connection memory, need to increase NLC in subgenelife.c\n");
                            connections[connused].oldlab=lab;
                            connections[connused].newlab=label[ij];
                            connections[connused].next=conn;
                            connections[connused].overlap++;
                            if(connprev) connections[connprev].next = connused++;
                            else connlists[label[ij]]=connused++;
                        }
                        else connections[conn].overlap++;
                    }
                    else {                                  // insert as last node in connection list added
                        if(connused>=NLC) fprintf(stderr,"Error, out of connection memory, need to increase NLC in subgenelife.c\n");
                        connections[connused].oldlab=lab;
                        connections[connused].newlab=label[ij];
                        connections[connused].overlap++;
                        // connections[connused].next=0;
                        if(connprev) connections[connprev].next = connused++;
                        else connlists[label[ij]]=connused++;
                    }
                }
            }
        }
    }

    for(i=1;i<=nlabel;i++) {                                 // count number of connections to old components for each component
        sum = 0;
        conn=connlists[i];
        while(conn) {
            sum++;
            conn=connections[conn].next;
        }
        connlen[i]=sum;
    }

    for(i=1;i<=nlabel;i++) {                                 // count number of connections from each old component to new components
        conn=connlists[i];
        while(conn) {
            connlenf[connections[conn].oldlab]++;
            conn=connections[conn].next;
        }
    }
    connpref[0]=0;
    for(i=1;i<=nlabel;i++) {                                 // weave forward connections from each old component to new components & record preference
        conn=connlists[i];
        maxlabel=0; maxoverlap = 0;
        while(conn) {
            lab=connections[conn].oldlab;
            connf = connlistsf[lab];
            connprev=0;
            while(connf && connections[connf].newlab<i) {   // change to < for compatibility with lapmod sparse array order
                connprev = connf;
                connf = connections[connf].nextf;
            }
            if (connf) {                                    // connections[connf].newlab<=i : do nothing if newlab==i, otherwise insert
                if(connections[connf].newlab>i) {           // insert node conn
                    connections[conn].nextf=connf;
                    if(connprev) connections[connprev].nextf = conn;
                    else connlistsf[lab]=conn;
                }
            }
            else {                                          // insert as last node in connection list added
                if(connprev) connections[connprev].nextf = conn;
                else connlistsf[lab]=conn;
            }
            if(connections[conn].overlap>maxoverlap) {
                maxoverlap=connections[conn].overlap;
                maxlabel = connections[conn].oldlab;
            }
            conn=connections[conn].next;
        }
        connpref[i]=maxlabel;
    }

    overallmaxoverlap = 0;
    connpreff[0]=0;
    for(i=1;i<=oldnlabel;i++) {                             // calclulate weighted overlaps and max overlap connection and overall max overlap
        sumoverlap = 0;
        connf=connlistsf[i];
        maxlabel=0; maxoverlap = 0;
        while(connf) {
            sumoverlap+=connections[connf].overlap;
            if (connections[connf].overlap>maxoverlap) {
                maxoverlap=connections[connf].overlap;
                maxlabel = connections[connf].newlab;
            }
            connf=connections[connf].nextf;
        }
        if(maxoverlap>overallmaxoverlap) overallmaxoverlap = maxoverlap;
        connpreff[i]=maxlabel;
        
        rsumoverlap = 1.0/(float) sumoverlap;
        connf=connlistsf[i];                                // save weighted forward overlaps : .woverlap
        while(connf) {
            connections[connf].woverlap=rsumoverlap*connections[connf].overlap;
            connf=connections[connf].nextf;
        }
    }
    
    for(i=1;i<=nlabel;i++) {                                 // construct alternative overlap aoverlap from woverlap
        rsumoverlap = 0.;                                    // reuse rsumoverlap as floating point sum before doing reciprocal
        conn=connlists[i];
        while(conn) {
            rsumoverlap+=connections[conn].woverlap;
            conn=connections[conn].next;
        }
        rsumoverlap = 1.0/rsumoverlap;
        conn=connlists[i];                                  // save alternative weighted backward overlaps : .aoverlap
        while(conn) {
            connections[conn].aoverlap=rsumoverlap*connections[conn].woverlap;
            conn=connections[conn].next;
        }
    }

    for(i=1;i<=nlabel;i++) {                               // prune other connections for mutually preferred matchings
        if(connpreff[connpref[i]]==i) {
            conn=connlists[i];
            while(conn) {
                connlenf[connections[conn].oldlab]++;
                conn=connections[conn].next;
            }
        }
    }

    for(i=0;i<N2;i++) {                                    // reset arrays used in LAP
        cclap[i]=0;
        kklap[i]=0;
    }
    for(i=0;i<NLM;i++) iilap[i]=0;
    nlap = 0;                                              // initialize counters
    nclap = 0;
    iilap[0] = 0;
    for(i=1;i<=oldnlabel;i++) {                            // setup cost matrix as sparse array for lapmod
        if(connpref[connpreff[i]]==i) {                    // use only this connection, pruning all others if mutually preferred
            kklap[nclap]=connpreff[i];                     // for maxmatch, labels for sparse cost matrix column index kk are in range from 1 to nlabel
            cclap[nclap]=1;                                // minimal but non zero cost
            nclap++;
        }
        else {
            connf=connlistsf[i];
            while(connf) {
                lab = connections[connf].newlab;
                if(connpreff[connpref[lab]]!=lab) {           // only connect to labels that are not preassigned by mutuality to another label
                    kklap[nclap]=connections[connf].newlab;   // for maxmatch labels for sparse cost matrix column index kk are newlab : and run from 1 to nlabel
                    cclap[nclap]=connections[connf].overlap ?  (1+overallmaxoverlap-connections[connf].overlap) :  100 * overallmaxoverlap;   // use if minimizing cost
                    // if (cclap[nclap]<0) fprintf(stderr,"error in cost matrix from genelife, negative value at %d\n",nclap);
                    nclap++;
                }
                connf=connections[connf].nextf;
            }
        }
        nlap++;                                           // end of row, possibility of zero entries in row
        iilap[nlap]=nclap;
    }

    nmatched=maxmatch(nlap,kklap,iilap,xlap,ylap,dist);
                                                          // optionally, print connections, preferred, matched, and list of possible
    if (connectout) {
        fprintf(stderr,"BACKWARD\n");
        for(i=1;i<=nlabel;i++) {                          // print backward connections
            fprintf(stderr,"step %5d: %3d bwd conn's for newlabel %4d prefers %4d assigned y%4d:",totsteps,connlen[i],i,connpref[i],ylap[i]);
            conn = connlists[i];
            while(conn) {
                fprintf(stderr," %4d(%3d)",connections[conn].oldlab,connections[conn].overlap);
                conn=connections[conn].next;
            }
            fprintf(stderr,"\n");
        }
        fprintf(stderr,"FORWARD\n");
        for(i=1;i<=oldnlabel;i++) {                       // print forward connections
            fprintf(stderr,"step %5d: %3d fwd conn's for oldlabel %4d prefers %4d assigned x%4d:",totsteps,connlenf[i],i,connpreff[i],xlap[i]);
            conn = connlistsf[i];
            while(conn) {
                fprintf(stderr," %4d(%3d)",connections[conn].newlab,connections[conn].overlap);
                conn=connections[conn].nextf;
            }
            fprintf(stderr,"\n");
        }
    }
    
    for(i=0;i<NLM;i++) relabel[i]=0;                    // now using matching to implement relabelling of new components
    if(totsteps==0) {
        fprintf(stderr,"totsteps 0 nlabel %d\n",nlabel);
        for(i=1;i<=nlabel;i++) relabel[i]=i;
    }
    else {
        for(i=1;i<=NLM;i++) queue_array[i] = 0;
        for(i=1;i<=nlabel;i++) if(ylap[i]) {
            relabel[i]=oldrelabel[ylap[i]];             // keep old label for these matched components
            queue_array[relabel[i]]=1;                  // mark this label as taken
        }
        for(ij=i=1;i<=nlabel;i++)  {
            if(!ylap[i]) {                              // if unmatched component
                while(queue_array[ij]) ij++;            // find next free label ij with relabel[ij]==0, i.e. not yet assigned
                relabel[i]=ij;
                ij++;
            }
        }
    }
    //  for(i=1;i<=nlabel;i++) fprintf(stderr,"relabel[%4d]=%4d ylap[%4d]=%4d\n",i,relabel[i],i,ylap[i]);


    for(nunique=0,i=1;i<=oldnlabel;i++) if (connpref[connpreff[i]] == i) nunique++;
    for(nzconnect=0,i=1;i<=oldnlabel;i++) if (!connlistsf[i]) nzconnect++;
    for(nremconnect=0,i=1;i<=oldnlabel;i++) if (iilap[i]==iilap[i-1]) nremconnect++;       // this includes those components with all connections removed by nunique pairs
    if(!(totsteps % 10) && colorupdate1 && connectedprints) {
        fprintf(stderr,"connected cpts:  %d(%d) matched(unique) & %d(%d) with no-residual(no) connections i.e. %d out of %d(%d) old(new) components\n",
                        nmatched,nunique,nremconnect,nzconnect,nmatched+nremconnect,oldnlabel,nlabel);
    }

    return nlabel;
}
//.......................................................................................................................................................
unsigned int extract_components(uint64_t gol[],uint64_t golg[]) {
    int i,j,ij,ij1,log2n;
    short unsigned int k,nside,patt;
    unsigned int histside[log2N+1];
    unsigned int nlabel;
    int wpixels;
    uint64_t hashkey,rand;
    
    for(i=1;i<=ncomponents;i++) {
        oldcomplist[i]=complist[i];                                                                     // copy whole component struct to previous time step list
    }
    oldncomponents=ncomponents;
    ncomponents = label_components(gol,golg);                                                                // label connected components at this time step
    nlabel = ncomponents;

    for(i=1;i<=nlabel;i++) {                                                                            // initialize component structures for horizontal scan
        complist[i].lastrc=0;
        complist[i].label=0;
        complist[i].pixels=0;
    }
    for (i=0;i<N;i++) {                                                                                 // find lateral limits of each component in horizontal scan
      for (j=0; j<N; j++) {
        ij = j*N+i;
        if(label[ij]>nlabel) {
            fprintf(stderr,"in extract_components step %d label %d out of bounds (%d) lateral at %d\n",totsteps,label[ij],ncomponents,ij);
            exit(1);
        }
        if (label[ij]) {                                                                                // if site labelled
            if (!complist[label[ij]].label) {                                                           // if label encountered for first time
                complist[label[ij]].label=label[ij];
                complist[label[ij]].E=complist[label[ij]].W=i;                                          // set both horizontal bounds to first encountered x for label
                complist[label[ij]].lastrc=i;                                                           // row contains this label
            }
            else if (complist[label[ij]].lastrc != i) {                                                 // label reencountered for first time in row
                if (((complist[label[ij]].lastrc+1)&(N-1)) == i) {                                      // continuation of component from previous row
                    if(complist[label[ij]].E>=complist[label[ij]].W) complist[label[ij]].E=i;
                }
                else if (((complist[label[ij]].lastrc+1)&(N-1)) <  i) {                                 // component resumes after row gap
                    // fprintf(stderr,"HORIZ TRACK step %d setting W %d after gap at i %d j %d label %d lastrc %d\n",totsteps,i,i,j,label[ij],complist[label[ij]].lastrc);
                    complist[label[ij]].W=i;
                }
                complist[label[ij]].lastrc=i;                                                           // row contains this label
            }
            complist[label[ij]].pixels++;
        }
      }
    }
    for(i=1;i<=nlabel;i++) {                                                                            // initialize component structures for vertical scan
        complist[i].lastrc=0;
        complist[i].label=0;
    }
    for (j=0;j<N;j++) {                                                                                 // find vertical limits of each component
      for (i=0; i<N; i++) {
        ij = j*N+i;
        if(label[ij]>nlabel) {
            fprintf(stderr,"in extract_components step %d label %d out of bounds (%d) vertical at %d\n",totsteps,label[ij],ncomponents,ij);
            exit(1);
        }
        if (label[ij]) {                                                                                // if site labelled
            if (!complist[label[ij]].label) {                                                           // if label encountered for first time
                complist[label[ij]].label=label[ij];
                complist[label[ij]].N=complist[label[ij]].S=j;                                          // set both vertical bounds to first encountered x for label
                complist[label[ij]].lastrc=j;                                                           // col contains this label
            }
            else if (complist[label[ij]].lastrc != j) {                                                 // label reencountered for first time in col
                if (((complist[label[ij]].lastrc+1)&(N-1)) == j) {                                      // continuation of component from previous col
                    if(complist[label[ij]].S>=complist[label[ij]].N) complist[label[ij]].S=j;
                }
                else if (((complist[label[ij]].lastrc+1)&(N-1)) <  j) {                                 // component resumes after col gap
                    // fprintf(stderr,"VERT  TRACK step %d setting N %d after gap at i %d j %d label %d lastrc %d\n",totsteps,j,i,j,label[ij],complist[label[ij]].lastrc);
                    complist[label[ij]].N=j;
                }
                complist[label[ij]].lastrc=j;                                                           // col contains this label
            }
        }
      }
    }

    for(i=1;i<=nlabel;i++) {
        nside=((complist[i].E-complist[i].W)&(N-1));                                                    // modulo calculation required if component crosses periodic boundary
        k = ((complist[i].S-complist[i].N)&(N-1));
        nside = nside > k ? nside+1 : k+1;                                                              // side length of square one more than greater of vertical & horiz. difference
        for(log2n=0;(1<<log2n)<nside;log2n++);
        complist[i].log2n = log2n;                                                                      // log2n is smallest power of 2 for side of square enclosing component
        nside = 1<<log2n;                                                                               // convert nside to next power of 2 for quadtree analysis;
        wpixels = 0;
        for(ij=0;ij<nside*nside;ij++) {
            ij1 =   ((complist[i].W+(ij&(nside-1)))&(N-1)) +                                            // calculate coordinate ij1 in full array of ij index in the component square
                 N*((complist[i].N+(ij>>log2n))&(N-1));
            working[ij] = (label[ij1]==i) ? 0x1ull : 0ull;                                              // use working array to store binary component image
            if (working[ij]) wpixels++;
        }
        hashkey = quadimage(working,&patt,log2n);                                                       // quadtree hash code of component image: needs to work for all nside=2^n
        if(!hashkey) {                                                                                  // components either have (quad!=NULL and patt==0) or (quad==NULL and patt!=0)
            complist[i].patt=patt;
            complist[i].quad=0ull;
        }
        else {
            complist[i].patt=0;
            complist[i].quad=hashkey;
        }
    }
    
    if(totsteps==0) {                                                                                         // initialize colors to labels for t=0
        for(i=1;i<=nlabel;i++) {
            complist[i].gcolor=(float)(i-1)/ (float) nlabel;
        }
    }
    else {                                                                                              // for t>0 mixcolors from overlapping connected components
        for(i=1;i<=nlabel;i++) {
            RAND128P(rand);
            complist[i].gcolor=mixcolor(i,rand);
        }
    }
    
    for(i=0;i<=log2N;i++) {
        histside[i]=0;
    }
    for(i=1;i<=nlabel;i++) {
        histside[complist[i].log2n]++;
    }

    if(!(totsteps % 10) && colorupdate1 && connectedprints) {
        fprintf(stderr,"histogram log2n ");for(i=0;i<=log2N;i++) fprintf(stderr," %5d",i);fprintf(stderr,"\n");
        fprintf(stderr,"conn cmpt counts");for(i=0;i<=log2N;i++) fprintf(stderr," %5d",histside[i]);fprintf(stderr,"\n");
    }

    return nlabel;
}
//.......................................................................................................................................................
int novelcells() {
    int ij,d,popcount,n;
    quadnode *q;
    uint64_t quad;
    
    n=0;
    if (diagnostics & diag_hash_patterns) {
        for (ij=0; ij<N2; ij++) {
            if (label[ij]) {
                quad = complist[label[ij]].quad;
                if((d=complist[label[ij]].patt)) popcount=smallpatts[d].activity; // number of times small (<= 4x4) pattern encountered previously
                else if((q = (quadnode *) hashtable_find(&quadtable, quad)) != NULL) // if we reach here, quad should have been stored in hash table
                    popcount=q->activity;                       // number of times large pattern encountered previously (poss. also as part of larger patt)
                else popcount=1;                                // should never occur, but just in case, assume novel
                if (popcount==1) n++;
            }
        }
    }
    return(n);
}
