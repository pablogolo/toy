A minimalistic program: reads a dataset, does a Wagner tree + TBR, saves final tree in parenthetical notation (TNT format) and as a tree-diagram. 

Compile with gcc -O2 -o toy toy.c

Usage: 

    toy data     
    
Results are written to standard output; hence, you can redirect results to a file with:

    toy data > outfile 
  
