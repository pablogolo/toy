All the programs included here are small, minimalistic standalone programs, intended to illustrate general principles of doing Wagner trees and branch-swapping. 

These programs (in increasing order of size/efficiency) are MINITOY / TOY / INCTOY / DNATOY / PARATOY.

TOY reads a dataset, does a Wagner tree + TBR, saves final tree in parenthetical notation (TNT format) and [optionally] as a tree-diagram. 

Compile with 

    gcc -O3 -o toy toy.c

Usage: 

    toy data     
    
Data format is a simple TNT file (starting with "xread" on top). Results are written to standard output; hence, you can redirect results to a file with:

    toy data > outfile 

To run the example datasets (in a folder of their own), you can try:

    minitoy zilla.tnt > result.tre 
    toy zilla.tnt > result.tre 
    inctoy zilla.tnt > result.tre 
or 
    dnatoy dnazilla.tnt > result.tre 
    paratoy dnazilla.tnt > result.tre 

The resulting tree-file (result.tre) is ready to open with TNT.

Command-line options are: 
    N        a number following the dataset: random seed for the addition sequence (0 = as is)
    -        no branch-swapping (just a Wagner tree)
    +xxxx    read a tree in parenthetical notation from file "xxxx" (must be binary, 
             and with taxa numbered from 0 to T-1) and then start branch-swapping from that tree
    t        will draw the tree in ASCII characters (default=don't)
    /        use SPR branch-swapping instead of TBR.
  
You can find a modified version, INCTOY, in inctoy.c, with an extra 50 lines or so, in the fixupdn() function. This uses incremental reoptimization during the clipping/unclipping phases, which cuts the times needed for completing TBR to about half (and the time needed to complete SPR by much more than that). I am not aware of any other open-source program that uses incremental reoptimization, so this serves to illustrate both the advantage of using it (when compared to the more vanilla toy.c), and the subtleties involved. That code is compiled in the same way: 

     gcc -O3 -o inctoy inctoy.c 

DNATOY, a version identical to this (in dnatoy.c, also using incremental optimization) differs only in reading DNA data (in IUPAC format), which some users may prefer. Compile with, 

    gcc -O3 -o dnatoy dnatoy.c 

The most efficient version is PARATOY (in paratoy.c), which uses Farris' (1996) methods for parallel optimization, jointly with all the other algorithms of the other versions. This is the longest version (500 loc), but also the fastest, with a speed comparable to that of PAUP* (but still substantially slower than TNT, which uses other algorithms to speed up branch-swapping).

Finally, in the other end of the spectrum, MINITOY shows the simplest possible: just a Wagner tree plus TBR, no randomization of addition sequence, no error checking, no progress report, no tree-reading or drawing. This is what you find in minitoy.c (only 200 lines of code; never saw a program that small that could do an actual phylogenetic analysis!). As this does not use incremental optimization, it is slower than INCTOY and DNATOY, with about the same speed as TOY. 
