A minimalistic program: reads a dataset, does a Wagner tree + TBR, saves final tree in parenthetical notation (TNT format) and as a tree-diagram. 

Compile with gcc -O2 -o toy toy.c

Usage: 

    toy data     
    
Data format is a simple TNT file (starting with "xread" on top). Results are written to standard output; hence, you can redirect results to a file with:

    toy data > outfile 
  
A modified version (inctoy.c, with an extra 50 lines or so, in the fixupdn() function) uses incremental reoptimization during the clipping/unclipping phasaes. This cuts the times needed for completing TBR to about half. I am not aware of any other open-source program that uses incremental reoptimization, so this serves to illustrate both the advantage of using it (when compared to the more vanilla toy.c), and the subtleties involved. That code is compiled in the same way: 

     gcc -O2 -o inctoy inctoy.c 
     
