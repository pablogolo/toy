/*********************************************************************\
|           A Toy program for phylogenetic analysis                   |
|                by Pablo A. Goloboff (2026)                          |
|         =============================================               |
|    A variaton of toy.c with incremental optimization and parallel   |
| optimization of characters (from Farris 1996; see Goloboff 2022).   |
| Free to use as is or modify (but please acknowledge authorship).    |
|         =============================================               |
|     The tree-printing algorithm is taken (simplified to handle      |
|     binary trees only) from Goloboff (2022 ISBN 9780367420277)      |
|   https://www.lillo.org.ar/phylogeny/eduscripts/treeplotting.pic    |
\*********************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <ctype.h>
FILE * inf , * tfile ; 
#define MAXT 5000
#define MAXC 500
#define MISSING 15
#define BIGINT 1000000000
#define TOP ( unsigned long long int ) 0x8888888888888888
#define BOT ( unsigned long long int ) 0x7777777777777777
void read_tree ( char * ) , preplot ( void ) , plot ( int ) ; 
int nt , nc , globest = BIGINT , dospr = 0 , wagstart = 1 , curnodin ; 
int anc[2*MAXT] , lef[2*MAXT] , rig[2*MAXT] , sis[2*MAXT] , oplist[2*MAXT] , trylist[2*MAXT] , rootlist[2*MAXT] , ladd[2*MAXT] ;
int atlin[2*MAXT] , gpsiz[2*MAXT] , btrack[2*MAXT] , rho = 218 , el = 192 , vert = 179 , hor = 196 , up = 193 , hook = 180 , linsdone = 0 ; 
char names[ MAXT ] [ 50 ] , inmask[256] , stepsin[ 65536 ] ;
unsigned long long int matrix[ 2*MAXT ][ MAXC ] , upmatrix[ 2*MAXT ][ MAXC ] , unmat[ 2*MAXT ][ MAXC ] , bakmat[ 2*MAXT ][ MAXC ] ;
unsigned long int rseed = 1 ;
void fixupdn( int , int , int ) ; 

void fillstepsin ( void ) {
    int a , b , val ;
    for ( a = 0 ; a < 65536 ; ++ a ) {
        val = 0 ; 
        for ( b = 0 ; b < 16 ; ++ b ) 
            if ( ( a & ( 1 << b ) ) ) ++ val ;
        stepsin[ a ] = val ; }}    
        
void errout ( char * txt ) {
    fprintf ( stderr , "Can't run:\n    %s\n" , txt ) ;
    exit ( 0 ) ; }

void save ( int i ) {
    if ( i < nt ) {
        printf ( "%i " , i ) ;
        return ; }
    if ( i == nt ) printf ( "\ntread\n" ) ;
    printf( "(" ) ;
    save( lef[ i ] ) ; 
    save( rig [ i ] ) ; 
    printf( ")" ) ;
    if ( i == nt ) printf ( ";\nproc/;\n" ) ; }

void tbreroot ( int old , int new ) {
   int i , j , k , m , n , nx , first ;
   if ( anc [ new ] == old || old < nt ) return ; 
   j = nx = first = anc [ i = new ] ;
   n = sis [ new ] ;
   while ( nx != old ) {
       k = i ; 
       i = j ;
       j = nx ;
       m = n ;
       n = sis [ j ] ; 
       if ( k == lef [ i ] ) { lef [ i ] = j ; sis [ sis [ rig [ i ] ] = j ] = rig [ i ] ; } 
       else { rig [ i ] = j ; sis [ sis [ lef [ i ] ] = j ] = lef [ i ] ; }
       nx = anc [ j ] ; 
       anc [ j ] = i ; }
   anc [ lef [ j ] = m ] = anc [ rig [ j ] = n ] = j ;
   sis [ sis [ n ] = m ] = n ;
   anc [ lef [ old ] = new ] = anc [ rig [ old ] = first ] = old ;
   sis [ sis [ new ] = first ] = new ; 
   return ; }

void inzert( int what , int root , int where ) {
    int dasis , danc , join ; 
    if ( root >= 0 ) 
       if ( anc[ root ] != what && root != what )
          tbreroot ( what , root ) ; 
    dasis = sis[ where ] ;
    danc = anc [ where ] ;
    join = anc [ what ] ; 
    lef [ anc [ where ] = join ] = what ;
    rig [ join ] = where ;
    sis [ sis [ where ] = what ] = where ;
    anc [ join ] = danc ;
    lef [ danc ] = dasis ;
    rig [ danc ] = join ;
    sis [ sis [ dasis ] = join ] = dasis ; }    

void setmask( void ) {
    int i ;
    for ( i = 0 ; i < 256 ; ++ i ) inmask[ i ] = 0 ;
    inmask [ 'a' ] = inmask [ 'A' ] = 1 ; 
    inmask [ 'g' ] = inmask [ 'G' ] = 2 ; 
    inmask [ 'c' ] = inmask [ 'C' ] = 4 ;
    inmask [ 'u' ] = inmask [ 'U' ] = inmask [ 't' ] = inmask [ 'T' ] = 8 ; 
    inmask [ 'R' ] = inmask [ 'r' ] = inmask [ 'a' ] | inmask [ 'g' ] ; 
    inmask [ 'Y' ] = inmask [ 'y' ] = inmask [ 't' ] | inmask [ 'c' ] ; 
    inmask [ 'W' ] = inmask [ 'w' ] = inmask [ 'a' ] | inmask [ 't' ] ; 
    inmask [ 'S' ] = inmask [ 's' ] = inmask [ 'c' ] | inmask [ 'g' ] ; 
    inmask [ 'M' ] = inmask [ 'm' ] = inmask [ 'a' ] | inmask [ 'c' ] ; 
    inmask [ 'K' ] = inmask [ 'k' ] = inmask [ 'g' ] | inmask [ 't' ] ; 
    inmask [ 'B' ] = inmask [ 'b' ] = inmask [ 'c' ] | inmask [ 'g' ] | inmask [ 't' ] ; 
    inmask [ 'D' ] = inmask [ 'd' ] = inmask [ 'a' ] | inmask [ 'g' ] | inmask [ 't' ] ; 
    inmask [ 'H' ] = inmask [ 'h' ] = inmask [ 'a' ] | inmask [ 'c' ] | inmask [ 't' ] ; 
    inmask [ 'V' ] = inmask [ 'v' ] = inmask [ 'a' ] | inmask [ 'c' ] | inmask [ 'g' ] ; 
    inmask [ 'N' ] = inmask [ 'n' ] = inmask [ '-' ] = inmask [ '?' ] = MISSING ; }
    
void readata ( void ) {
    int a , b , ispoly , usnc , inn ;
    unsigned long long int set , curshif ; 
    char ch , command[20] ; 
    setmask () ;
    fscanf( inf , " %s" , &command ) ;
    if ( strcmp( command, "xread" ) ) errout ( "No xread?" ) ; 
    fscanf ( inf , " %i %i" , &usnc , &nt ) ;
    if ( nt > MAXT ) errout ( "Too many taxa" ) ;
    if ( nt < 4 ) errout ( "Too few taxa" ) ; 
    if ( ( nc + 15 ) / 16 > MAXC ) errout ( "Too many characters" ) ; 
    if ( usnc < 1 ) errout ( "Too few characters " ) ;
    nc = ( usnc + 15 ) / 16 ; 
    for ( a = 0 ; a < nt ; ++ a )
       for ( b = 0 ; b < nc ; ++ b )
          matrix[ a ] [ b ] = 0 ; 
    for ( a = 0 ; a < nt ; ++ a ) {
       fscanf ( inf , " %s" , &names[ a ] ) ;
       for ( b = curshif = inn = 0 ; b < usnc ; ++ b ) {
           fscanf ( inf , " %c" , &ch ) ;
           ispoly = 0 ; 
           if ( ch == '[' ) { ispoly = 1 ; fscanf( inf , " %c" , &ch ) ; }
           set = 0 ; 
           while ( ch != ']' ) {
               if ( !inmask[ ( int ) ch ] ) errout ( "Bad symbol" ) ;
               set |= ( unsigned long long int ) inmask[ ( int ) ch ] ;
               if ( !ispoly ) ch = ']' ;
               else fscanf ( inf , " %c" , &ch ) ; }
           matrix[ a ] [ inn ] |= set << curshif ;
           upmatrix[ a ] [ inn ] = matrix [ a ] [ inn ] ; 
           if ( !a ) unmat[ a ] [ inn ] = matrix[ a ] [ inn ] ;
           curshif += 4 ;
           if ( curshif >= 64 ) { curshif = 0 ; ++ inn ; }}}
    if ( usnc % 16 ) { // Fill empty fields... 
        curshif = ( usnc % 16 ) * 4 ;
        set = 0 ;
        inn = nc - 1 ;
        while ( curshif < 64 ) {
            set |= ( unsigned long long int ) 15 << curshif ;
            curshif += 4 ; }
        for ( a = 0 ; a < nt ; ++ a ) {
           matrix[ a ] [ inn ] |= set ; 
           upmatrix[ a ] [ inn ] = matrix[ a ] [ inn ] ;
           if ( !a ) unmat[ a ] [ inn ] = matrix[ a ] [ inn ] ; }}}

int listabove( int * list , int from , int doterms ) {
    int * done , * todo , at ;
    * ( done = list ) = from ;
    todo = done + 1 ;
    while ( done < todo ) {
       at = * done ++ ;
       if ( at >= nt ) {
          if ( lef [ at ] >= nt || doterms ) * todo ++ = lef [ at ] ; 
          if ( rig [ at ] >= nt || doterms ) * todo ++ = rig [ at ] ; }}
    return done - list ; }

int fichop( int passes , int from ) {
    int i , n , b , at , steps = 0 , blow ;
    unsigned long long int x , y , z ; 
    for ( i = listabove( oplist , from , 0 ) ; i -- ; ) {
       at = oplist[ i ] ;
       for ( b = 0 ; b < nc ; ++ b ) {
          x = matrix[ lef [ at ] ] [ b ] & matrix [ rig [ at ] ] [ b ] ;
          y = matrix[ lef [ at ] ] [ b ] | matrix [ rig [ at ] ] [ b ] ;
          z = TOP & ~ ( x | ( ( x & BOT ) + BOT ) ) ;  // <-- From Farris (1996!)
          if ( z )
              steps += stepsin[ 65535 & ( z | ( z >> 17 ) | ( z >> 34 ) | ( z >> 51 ) ) ] ;
          z |= z - ( z >> 3 ) ;   // <-- From Farris (1996!)
          bakmat[ at ] [ b ] = y | z ; 
          matrix[ at ] [ b ] = x | ( y & z ) ; }}
    if ( passes == 2 ) {
        n = listabove( oplist , from , 1 ) ;
        for ( b = 0 ; b < nc ; ++ b ) {
           unmat [ at ] [ b ] = matrix[ at ] [ b ] ; 
           upmatrix [ at ] [ b ] = matrix[ at ] [ b ] ; }
        for ( i = 1 ; i < n ; ++ i ) 
          for ( b = 0 , blow = anc [ at = oplist[ i ] ] ; b < nc ; ++ b ) {
            if ( !at ) continue ;
            if ( at < nt ) {
               x = upmatrix[ blow ] [ b ] & ~ matrix[ at ] [ b ] ;
               x = TOP & ( x | ( ( x & BOT ) + BOT ) ) ;
               x |= x - ( x >> 3 ) ;
               unmat[ at ] [ b ] = upmatrix[ blow ] [ b ] | ( matrix[ at ] [ b ] & x ) ; } 
            else {
               x = upmatrix[ blow ] [ b ] & ~ matrix[ at ] [ b ] ;
               x = ( x | ( ( x & BOT ) + BOT ) ) & TOP ;
               x |= x - ( x >> 3 ) ;
               upmatrix[ at ] [ b ] = ( upmatrix[ blow ] [ b ] & ~ x ) | ( x & ( matrix[ at ] [ b ] | upmatrix[ blow ] [ b ] & bakmat[ at ] [ b ] ) ) ;
               unmat[ at ] [ b ] = upmatrix[ at ] [ b ] | upmatrix[ blow ] [ b ] ; }}}
    return steps ; }        

int testamove ( int what , int root , int where ) {
    int b , steps = 0 ;
    unsigned long long int x ; 
    for ( b = 0 ; b < nc && steps < globest ; ++ b ) {
        if ( what == root && what < nt ) 
             x = matrix[ root ] [ b ] & unmat[ where ] [ b ] ;
        else x = unmat[ root ] [ b ] & unmat[ where ] [ b ] ;
        x = TOP & ~ ( x | ( ( x & BOT ) + BOT ) ) ;
        if ( x ) steps += stepsin[ 65535 & ( x | ( x >> 17 ) | ( x >> 34 ) | ( x >> 51 ) ) ] ; }
    return steps ; }

void randomize( void ) {
   int i , n , q , * lp = ladd + 1 , left = nt - 1 ; 
   for ( i = 0 ; i < nt ; ++ i ) ladd [ i ] = i ;
   if ( !rseed ) return ; 
   srand ( rseed ) ;
   for ( i = 1 ; i < nt ; ++ i ) {
      q = * lp ; 
      * lp = lp [ n = rand () % left -- ] ;
      lp [ n ] = q ;
      ++ lp ; }}

void nitnet( void ) {
    randomize () ;
    anc[ lef [ nt ] = 0 ] = anc [ rig[ nt ] = nt + 1 ] = nt ;
    anc[ lef[ nt + 1 ] = ladd [ 1 ] ] = anc[ rig[ nt + 1 ] = ladd [ 2 ] ] = nt + 1 ;
    sis [ sis[ 0 ] = nt + 1 ] = 0 ;
    sis [ sis [ ladd [ 1 ] ] = ladd [ 2 ] ] = ladd [ 1 ] ; }

void wagner( void ) {
   int i , now , at , next = 3 , this , besloc , curnod = nt + 2 ; 
   nitnet () ;
   fichop ( 2 , nt ) ; 
   while ( next < nt ) {
       anc[ now = ladd [ next ++ ] ] = curnod ++ ;
       globest = BIGINT ;
       for ( i = listabove( trylist , sis[ 0 ] , 1 ) ; i -- ; ) {
          at = trylist[ i ] ;
          this = testamove( now , now , at ) ;
          if ( this < globest ) {
              besloc = at ;
              if ( !( globest = this ) ) break ; }}
       inzert( now , -1 , besloc ) ;
       fixupdn( anc[ now ] , -1 , 0 ) ; }
   globest = fichop ( 1 , nt ) ; 
   fprintf ( stderr , "Formed a Wagner tree of %i steps\n" , globest ) ; }    

void unzert( int cut ) {
    int danc = anc[ cut ] , ddanc , brother , cousin ;
    if ( danc == nt ) return ; 
    ddanc = anc [ danc ] ;
    brother = sis[ cut ] ;
    cousin = sis [ danc ] ;
    lef[ ddanc ] = brother ;
    rig [ ddanc ] = cousin ;
    sis [ sis [ cousin ] = brother ] = cousin ;
    anc [ brother ] = anc [ cousin ] = ddanc ; }    

void tbrswap( void ) {
    int i , j , rounds , cut = 0 , besroot , besloc , this , nitroot , nitloc , nroots , nlocs , lasreport = 0 , nmoves = 0 , nswitches = 0 ; 
    fichop ( 2 , nt ) ; 
    for ( rounds = 1 ; rounds < 2*nt - 1 ; ++ rounds ) {
        if ( ++ cut > 2*nt-2 ) cut = 1 ;
        if ( cut == nt ) continue ;
        i = ( rounds * 100 ) / ( 2*nt-1 ) ;
        if ( ( i - lasreport ) >= 5 ) {
            fprintf ( stderr , "\r" ) ; 
            for ( j = 0 ; j < i ; j += 5 ) fprintf ( stderr , "X" ) ;
            for ( j = i ; j < 100 ; j += 5 ) fprintf ( stderr , "=" ) ;
            fprintf ( stderr , " TBR" ) ;
            lasreport = i ; } 
        if ( cut >= nt ) nitroot = besroot = lef [ cut ] ;
        else nitroot = besroot = cut ;
        nitloc = besloc = sis[ cut ] ; 
        if ( !nitloc ) {
            j = listabove( trylist , cut , 1 ) ;
            fichop ( 2 , sis[ 0 ] ) ;
            globest = BIGINT ;
            globest = testamove( 0 , 0 , cut ) ;
            for ( i = 1 ; i < j ; ++ i ) {
               ++ nmoves ; 
               this = testamove( 0 , 0 , trylist[ i ] ) ;
               if ( this < globest ) {
                   rounds = 1 ; 
                   besroot = 0 ; 
                   besloc = trylist[ i ] ;
                   if ( !( globest = this ) ) break ; }}}
        else {
            unzert( cut ) ;
            nlocs = listabove( trylist , sis[ 0 ] , 1 ) ;
            if ( cut < nt ) { nroots = 3 ; rootlist[ 2 ] = cut ; }
            else nroots = listabove( rootlist , cut , 1 ) ;
            fixupdn( anc[ nitloc ] , cut , 0 ) ; 
            if ( dospr ) nroots = 3 ; 
            globest = BIGINT ; 
            globest = testamove ( cut , rootlist[ 2 ] , nitloc ) ;
            for ( i = 2 ; i < nroots ; ++ i )
               for ( j = 0 ; j < nlocs ; ++ j ) {
                   ++ nmoves ; 
                   this = testamove( cut , rootlist[ i ] , trylist[ j ] ) ;
                   if ( this < globest ) {
                       rounds = 1 ; 
                       besroot = rootlist[ i ] ; 
                       besloc = trylist[ j ] ;
                       if ( !( globest = this ) ) { i = nroots ; break ; }}}}
        if ( nitloc != besloc || besroot != nitroot ) ++ nswitches ;
        if ( !nitloc ) {
            if ( nitloc != besloc )
               tbreroot( cut , besloc ) ;
            fichop ( 2 , nt ) ; }
        else {
            inzert( cut , besroot , besloc ) ;
            if ( cut < nt ) fixupdn ( anc[ cut ] , -1 , 0 ) ; 
            else fixupdn ( anc[ cut ] , -1 , nitroot ) ; }}
  globest = fichop( 1 , nt ) ;
  if ( dospr ) fprintf ( stderr , "\rCompleted SPR (%i moves, %i accepted), best score %i\n" , nmoves , nswitches , globest ) ; 
  else fprintf ( stderr , "\rCompleted TBR (%i moves, %i accepted), best score %i\n" , nmoves , nswitches , globest ) ; }

void main ( int argc , char ** argv ) {
    int doswap = 1 , i , dotree = 0 ; 
    clock_t initime ;
    if ( sizeof ( unsigned long long int ) != 8 ) errout ( "Code requires compilation in a\nsystem where \"sizeof ( long long)\" is 8.\nThis is no such!!" ) ; 
    if ( argc == 1 ) errout ( "Give file name" ) ;
    if ( ( inf = fopen ( argv[ 1 ] , "rb" ) ) == NULL ) errout ( "Can't open file" ) ;
    if ( argc > 2 && argv[ 2 ] [ 0 ] != '-' && argv[ 2 ] [ 0 ] != '+' ) 
        rseed = atoi ( argv[ 2 ] ) ;
    fillstepsin () ; 
    readata () ;
    for ( i = 2 ; i < argc ; ++ i ) { 
        if ( argv[ i ][ 0 ] == '-' ) doswap = 0 ; 
        if ( argv[ i ][ 0 ] == '/' ) dospr = 1 ; 
        if ( argv[ i ][ 0 ] == '+' ) {
            wagstart = 0 ; 
            read_tree ( argv[ i ] + 1 ) ; }
        if ( argv[ i ][ 0 ] == 't' ) dotree = 1 ; }
    initime = clock () ; 
    if ( wagstart ) wagner () ;
    else {
        int len = fichop ( 1 , nt ) ;
        fprintf ( stderr , "Start swapping from user tree (%i steps)\n" , len ) ; }
    if ( doswap ) tbrswap () ; 
    fprintf ( stderr , "Time used: %.2f sec\n" , ( double ) ( clock() - initime ) / CLOCKS_PER_SEC ) ;
    fflush ( stderr ) ; 
    save( nt ) ;
    if ( dotree ) { 
       preplot () ;
       plot ( nt ) ; }}

/*** SOME EXTRA GOODIES : read trees in parenthetical notation, draw tree diagrams ***/

void add ( int what , int where ) {
    anc [ what ] = where ; 
    if ( lef [ where ] < 0 ) lef[ where ] = what ;
    else if ( rig[ where ] < 0 ) sis [ sis [ rig [ where ] = what ] = lef[ where ] ] = what ;
         else errout ( "Tree has polytomies" ) ; }

void treein ( int at ) {
    char c ;
    int i ;
    while ( 1 ) {
        fscanf ( tfile , " %c" , &c ) ;
        if ( c == ')' ) return ;
        if ( c == '(' ) {
            if ( ++ curnodin >= 2*nt-1 ) errout ( "Too many parentheses" ) ;
            add ( curnodin , at ) ;
            treein ( curnodin ) ; }
        else {
            ungetc ( c , tfile ) ;
            fscanf ( tfile , "%i" , &i ) ;
            if ( i >= nt ) errout ( "Taxon number" ) ;
            add ( i , at ) ; }}}             
    
void read_tree ( char * fnam ) {
    int i ;
    char command[20] ; 
    if ( ( tfile = fopen ( fnam , "rb" ) ) == NULL ) errout ( "Tree file" ) ;
    for ( i = 0 ; i < 2*nt-1 ; ++ i )
        anc[ i ] = lef[ i ] = rig[ i ] = sis[ i ] = -1 ;
    fscanf ( tfile , " %s" , &command ) ;
    if ( strcmp ( command , "tread" ) ) errout ( "No \"tread\" in file" ) ;
    i = ' ' ; while ( isspace ( i ) ) i = getc ( tfile ) ;
    if ( i == 39 ) {
        i = ' ' ; while ( i != 39 ) i = getc ( tfile ) ;
        i = ' ' ; while ( isspace ( i ) ) i = getc ( tfile ) ; }
    if ( i != '(' ) errout ( "Parentheses" ) ;
    treein ( curnodin = nt ) ; }

void initlines ( int which ) {  // Find out at which line every terminal is printed 
  if ( which >= nt ) {
     initlines ( lef [ which ] ) ;  
     initlines ( rig [ which ] ) ;  
     return ; }
  atlin [ which ] = ++ linsdone ; }     

void calculines ( int which ) {  // Find out at which line every internal node is printed 
  if ( which < nt ) return ; 
  calculines ( lef [ which ] ) ;  
  calculines ( rig [ which ] ) ; 
  atlin [ which ] = ( 1 + atlin [ lef [ which ] ] +  atlin [ rig [ which ] ] ) / 2 ; }     

void preplot ( void ) {
    int i , m , n ;
    for ( i = 0 ; i < nt ; ++ i ) gpsiz[ i ] = 1 ;
    for ( n = listabove( trylist , nt , 0 ) ; n -- ; ) { // rotate branches for aesthetics
        i = trylist[ n ] ; 
        if ( gpsiz[ lef[ i ]] > gpsiz[ rig [ i ]] ) { m = lef[ i ] ; lef[ i ] = rig [ i ] ; rig [ i ] = m ; }
        gpsiz[ i ] = gpsiz[ lef[ i ]] + gpsiz[ rig [ i ]] ; }
    initlines ( nt ) ;
    calculines ( nt ) ;
    linsdone = 0 ; }

void plot ( int which ) {
  int i = which , j , n , joint , * retrop = btrack ;
  if ( which >= nt ) {
      plot ( lef[ which ] ) ;
      plot ( rig[ which ] ) ;  
      return ; }
  while ( i != nt ) * retrop ++ = i = anc [ i ] ; 
  linsdone ++ ; 
  while ( retrop > btrack ) {
      j = * -- retrop ;
      if ( linsdone == atlin [ j ] ) printf( "%c%c" , hor , hor ) ; 
      else printf ( "  " ) ; 
      joint = ' ' ; 
      if ( linsdone <= atlin [ rig [ j ] ] ) {
          if ( linsdone >= atlin [ lef [ j ] ] ) joint = vert ; 
          if ( linsdone == atlin [ j ] ) joint = hook ; 
          for ( i = lef [ j ] , n = 0 ; n < 2 ; i = sis [ i ] , n ++ ) 
            if ( atlin [ i ] == linsdone ) {
              if ( i == lef [ j ] ) joint = rho ; 
              else if ( gpsiz[ j ] == 2 ) joint = up ; 
                   else joint = el ; }}
      printf ( "%c%" , joint ) ; }
  printf( "%c%c %s\n" , hor , hor , names[ which ] ) ; }  

void fixupdn( int from , int side , int nitroot ) {
    int b , n , nn , at , blow , was , canbreak , hasbreak = 1 , startat = from ;
    unsigned long long int w , x , y , z ; 
    if ( nitroot > 0 ) {
           hasbreak = 0 ;
           startat = nitroot ;
           if ( startat < nt ) startat = anc[ nitroot ] ; } 
    for ( b = 0 ; b < nc ; ++ b ) {
       upmatrix[ at = startat ] [ b ] = 0 ;
       canbreak = hasbreak ; 
       while ( 1 ) {
          x = matrix[ lef [ at ] ] [ b ] & matrix [ rig [ at ] ] [ b ] ;
          y = matrix[ lef [ at ] ] [ b ] | matrix [ rig [ at ] ] [ b ] ;
          z = TOP & ~ ( x | ( ( x & BOT ) + BOT ) ) ;
          z |= z - ( z >> 3 ) ;
          bakmat[ at ] [ b ] = y | z ; 
          w = x | ( y & z ) ; 
          upmatrix[ at ] [ b ] = 0 ;
          if ( w == matrix[ at ] [ b ] && at != from ) 
              if ( canbreak ) break ;
          if ( at == from ) canbreak = 1 ;     
          matrix[ at ] [ b ] = w ;
          if ( sis [ at ] >= 0 ) upmatrix[ sis [ at ] ] [ b ] = 0 ;
          if ( at == nt ) break ; 
          at = anc [ at ] ; }
       if ( at == nt ) { oplist[ n = 0 ] = sis[ 0 ] ; upmatrix[ nt ] [ b ] = matrix[ nt ] [ b ] ; }
       else oplist[ n = 0 ] = at ;
       nn = 1 ;
       if ( side > 0 ) 
           if ( side < nt ) unmat[ side ] [ b ] = matrix[ side ] [ b ] ;
           else
              if ( upmatrix[ side ] [ b ] != matrix[ side ] [ b ] ) {
                upmatrix[ side ] [ b ] = matrix[ side ] [ b ] ; 
                oplist[ nn ++ ] = lef [ side ] ; 
                oplist[ nn ++ ] = rig [ side ] ; }
       while ( n < nn ) {
           blow = anc [ at = oplist [ n ++ ] ] ;
           if ( !at ) continue ;
           was = upmatrix[ at ] [ b ] ; 
           if ( at < nt ) {
               x = upmatrix[ blow ] [ b ] & ~ matrix[ at ] [ b ] ;
               x = ( x | ( ( x & BOT ) + BOT ) ) & TOP ;
               x |= x - ( x >> 3 ) ;
               unmat[ at ] [ b ] = upmatrix[ blow ] [ b ] | ( matrix[ at ] [ b ] & x ) ; } 
           else {
               x = upmatrix[ blow ] [ b ] & ~ matrix[ at ] [ b ] ;
               x = ( x | ( ( x & BOT ) + BOT ) ) & TOP ;
               x |= x - ( x >> 3 ) ;
               upmatrix[ at ] [ b ] = ( upmatrix[ blow ] [ b ] & ~ x ) | ( x & ( matrix[ at ] [ b ] | upmatrix[ blow ] [ b ] & bakmat[ at ] [ b ] ) ) ;
               unmat[ at ] [ b ] = upmatrix[ at ] [ b ] | upmatrix[ blow ] [ b ] ; }
           if ( upmatrix[ at ] [ b ] != was ) 
             if ( at >= nt ) {
                 oplist[ nn ++ ] = lef [ at ] ; 
                 oplist[ nn ++ ] = rig [ at ] ; }}}}

