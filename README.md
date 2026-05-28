TOY - A minimalistic program: reads a dataset, does a Wagner tree + TBR, saves final tree in parenthetical notation (TNT format) and as a tree-diagram. 

Compile with gcc -O2 -o toy toy.c

Usage: 

    toy data     
    
Data format is a simple TNT file (starting with "xread" on top). Results are written to standard output; hence, you can redirect results to a file with:

    toy data > outfile 
  
You can find a modified version, INCTOY, in inctoy.c, with an extra 50 lines or so, in the fixupdn() function. This uses incremental reoptimization during the clipping/unclipping phasaes, which cuts the times needed for completing TBR to about half. I am not aware of any other open-source program that uses incremental reoptimization, so this serves to illustrate both the advantage of using it (when compared to the more vanilla toy.c), and the subtleties involved. That code is compiled in the same way: 

     gcc -O2 -o inctoy inctoy.c 
     
Finally, in the other end of the spectrum, MINITOY shows the simplest possible: just a Wagner tree plus TBR, no randomization of addition sequence, no error checking, no propgress report, no tree-reading or drawing. This is what you find in minitoy.c (only 200 lines of code; never saw a program that small that could do an actual phylogenetic analysis!). As this does not use incremental optimization, it has about the same speed as TOY. 
