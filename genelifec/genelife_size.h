//
// genelife_size.h
// project genelife
// definition of model size via log2N and derived parameters N,N2,NLM,NLC,Nmask,N2mask: N is length of array side
// declaration of constants specifying array sizes using enum types, to facilitate run time checking of fixed array bounds as well as type checking
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
#ifndef genelife_size_h
#define genelife_size_h

#ifndef OEX
#define OEX extern
#endif  /* OEX */

enum {log2N = 9,                	// toroidal array of side length N = 2 to the power of log2N (minimum log2N is 6 i.e. 64x64)
  	  N = 0x1 << log2N,         	// only side lengths powers of 2 allowed to enable efficient implementation of periodic boundaries
 	  N2 = N*N,                 	// number of sites in square-toroidal array
 	  Nmask = N - 1,            	// bit mask for side length, used instead of modulo operation
 	  N2mask = N2 - 1};          	// bit mask for array, used instead of modulo operation
enum {NLM = N2,                     // maximum number of discrete components possible N*N
      NLC = N2<<2};                 // maximum number of connections N*N*4
#endif /* genelife_size_h */
