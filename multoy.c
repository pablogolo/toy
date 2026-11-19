/*********************************************************************\
|           A Toy program for phylogenetic analysis                   |
|                by Pablo A. Goloboff (2026)                          |
|         =============================================               |
|      This is a variation of toy.c with incremental optimization     |
|      during clipping/unclipping phases, and the ability to save     |
|      multiple trees. Feel free to use as is, or modify.             |
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
#define MAXT 3000
#define MAXC 2500
#define MISSING 1023
#define BIGINT 1000000000
int maxtrees = 100 ; 
void read_tree ( char * ) , preplot ( void ) , plot ( int ) ; 
int nt , nc , matrix[ 2*MAXT ][ MAXC ] , upmatrix[ 2*MAXT ][ MAXC ] , unmat[ 2*MAXT ][ MAXC ] , bakmat[ 2*MAXT ][ MAXC ] , inmask[256] , globest , gloglobest = BIGINT , dospr = 0 , wagstart = 1 , curnodin ; 
int anc[2*MAXT] , lef[2*MAXT] , rig[2*MAXT] , sis[2*MAXT] , oplist[2*MAXT] , trylist[2*MAXT] , rootlist[2*MAXT] , ladd[2*MAXT] , tmptbuf[2*MAXT] , supp[2*MAXT] ;
int atlin[2*MAXT] , gpsiz[2*MAXT] , btrack[2*MAXT] , rho = 218 , el = 192 , vert = 179 , hor = 196 , up = 193 , hook = 180 , linsdone ; 
char names[ MAXT ] [ 50 ] ;
unsigned long int rseed = 1 ;
void fixupdn( int , int , int ) ;
int ** stortree , * tkeys , nstored = 0 , nrepls = 10 , nmoves = 0 , ngoodmvs , multrees = 0 , nowswapping = 0 , oplimit , verbose = 0 ;
typedef struct { int dest , root ; } GoodMvTyp ;
GoodMvTyp goodmvs[ ((MAXT+1)/2)*((MAXT+1)/2) ] ; 
void compare_n_save ( GoodMvTyp * , int , int ) ; 
#define NEGBIT ( 1 << 31 )
#define GETANC ( ( 1 << 30 ) - 1 )

void errout ( char * txt ) {
    fprintf ( stderr , "Can't run:\n    %s\n" , txt ) ;
    exit ( 0 ) ; }

void save ( int i ) {
    if ( i < nt ) {
        printf ( "%i " , i ) ;
        return ; }
    printf( "(" ) ;
    save( lef[ i ] ) ; 
    save( rig [ i ] ) ; 
    printf( ")" ) ; }

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
    for ( i = 0 ; i < 10 ; ++ i ) inmask[ '0' + i ] = 1 << i ;
    inmask[ '?' ] = inmask[ '-' ] = MISSING  ;
    inmask[ '[' ] = 1 ; }
    
void readata ( void ) {
    int a , b , ispoly ;
    char ch , command[20] ; 
    setmask () ;
    fscanf( inf , " %s" , &command ) ;
    if ( strcmp( command, "xread" ) ) errout ( "No xread?" ) ; 
    fscanf ( inf , " %i %i" , &nc , &nt ) ;
    if ( nt > MAXT ) errout ( "Too many taxa" ) ;
    if ( nt < 4 ) errout ( "Too few taxa" ) ; 
    if ( nc > MAXC ) errout ( "Too many characters" ) ; 
    if ( nc < 1 ) errout ( "Too few characters " ) ; 
    for ( a = 0 ; a < nt ; ++ a ) {
       fscanf ( inf , " %s" , &names[ a ] ) ;
       for ( b = 0 ; b < nc ; ++ b ) {
           fscanf ( inf , " %c" , &ch ) ;
           ispoly = 0 ; 
           if ( ch == '[' ) { ispoly = 1 ; fscanf( inf , " %c" , &ch ) ; }
           matrix[ a ] [ b ] = 0 ; 
           while ( ch != ']' ) {
               if ( !inmask[ ( int ) ch ] ) errout ( "Bad symbol" ) ;
               if ( inmask[ ( int ) ch ] != MISSING ) 
                    matrix[ a ] [ b ] |= inmask[ ( int ) ch ] ;
               else matrix[ a ] [ b ] |= MISSING ;
               if ( !ispoly ) ch = ']' ;
               else fscanf ( inf , " %c" , &ch ) ; }
           upmatrix[ a ] [ b ] = matrix [ a ] [ b ] ; 
           if ( !a ) unmat[ a ] [ b ] = matrix[ a ] [ b ] ; }}}

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
    int i , n , b , at , w , x , steps = 0 , bak , blow ;
    for ( i = listabove( oplist , from , 0 ) ; i -- ; ) {
       at = oplist[ i ] ;
       for ( b = 0 ; b < nc ; ++ b ) {
          bak = matrix[ lef [ at ] ] [ b ] | matrix [ rig [ at ] ] [ b ] ; 
          if ( ( x = matrix[ lef [ at ] ] [ b ] & matrix [ rig [ at ] ] [ b ] ) == 0 ) {
              x = bak ; 
              bak = 1023 ;
              ++ steps ; }
          bakmat[ at ] [ b ] = bak ; 
          matrix[ at ] [ b ] = x ; }}
    if ( passes == 2 ) {
        n = listabove( oplist , from , 1 ) ;
        for ( b = 0 ; b < nc ; ++ b ) {
           unmat [ at ] [ b ] = matrix[ at ] [ b ] ; 
           upmatrix [ at ] [ b ] = matrix[ at ] [ b ] ; }
        for ( i = 1 ; i < n ; ++ i ) 
          for ( b = 0 , blow = anc [ at = oplist[ i ] ] ; b < nc ; ++ b ) {
            if ( !at ) continue ; 
            w = matrix[ at ] [ b ] & upmatrix[ blow ] [ b ] ;
            if ( w == upmatrix[ blow ] [ b ] ) {
               if ( at > nt ) upmatrix [ at ] [ b ] = upmatrix[ blow ] [ b ] ;
               unmat[ at ] [ b ] = upmatrix[ blow ] [ b ] ; }
            else 
               if ( at < nt ) unmat[ at ] [ b ] = upmatrix[ blow ] [ b ] | matrix [ at ] [ b ] ; 
               else {
                   upmatrix[ at ] [ b ] = matrix[ at ] [ b ] | ( upmatrix[ blow ] [ b ] & bakmat[ at ] [ b ] ) ;
                   unmat [ at ] [ b ] = upmatrix[ at ] [ b ] | upmatrix[ blow ] [ b ] ; }}}
    return steps ; }        

int testamove ( int what , int root , int where ) {
    int b , steps = 0 ;
    for ( b = 0 ; b < nc && steps < oplimit ; ++ b ) 
        if ( what == root && what < nt ) {
            if ( ! ( matrix[ root ] [ b ] & unmat[ where ] [ b ] ) ) ++ steps ; }
        else 
            if ( ! ( unmat[ root ] [ b ] & unmat[ where ] [ b ] ) ) ++ steps ;
    return steps ; }

void randomize( void ) {
   int i , n , q , * lp = ladd + 1 , left = nt - 1 ;
   static int firstcall = 1 ; 
   if ( firstcall ) {
      if ( !rseed ) { firstcall = 0 ; return ; }
      srand ( rseed ) ;
      firstcall = 0 ; }
   for ( i = 0 ; i < nt ; ++ i ) ladd [ i ] = i ;
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
       oplimit = globest = BIGINT ;
       for ( i = listabove( trylist , sis[ 0 ] , 1 ) ; i -- ; ) {
          at = trylist[ i ] ;
          this = testamove( now , now , at ) ;
          if ( this < globest ) {
              besloc = at ;
              if ( !( globest = this ) ) break ;
              oplimit = globest ; }}
       inzert( now , -1 , besloc ) ;
       fixupdn( anc[ now ] , -1 , 0 ) ; }
   globest = fichop ( 1 , nt ) ; }    

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

int tbrswap( int ismult ) {
    int i , j , rounds , cut = 0 , besroot , besloc , this , nitroot , nitloc , nroots , nlocs , nswitches = 0 , totbest , totscore , success = 0 , divsteps ;
    GoodMvTyp * goodmvpt ; 
    totbest = globest = fichop ( 2 , nt ) ;
    for ( rounds = 1 ; rounds < 2*nt - 1 ; ++ rounds ) {
        if ( ++ cut > 2*nt-2 ) cut = 1 ;
        if ( cut == nt ) continue ;
        if ( cut >= nt ) nitroot = besroot = lef [ cut ] ;
        else nitroot = besroot = cut ;
        nitloc = besloc = sis[ cut ] ;
        if ( !nitloc ) {
            j = listabove( trylist , cut , 1 ) ;
            fichop ( 2 , sis[ 0 ] ) ;
            oplimit = BIGINT ;
            globest = testamove( 0 , 0 , cut ) ;
            if ( ismult && nstored < maxtrees ) oplimit = globest + 1 ;
            else oplimit = globest ;
            divsteps = totbest - globest ; 
            goodmvpt = goodmvs ;
            if ( oplimit ) 
             for ( i = 1 ; i < j ; ++ i ) {
               ++ nmoves ; 
               this = testamove( 0 , 0 , trylist[ i ] ) ;
               totscore = divsteps + this ; 
               if ( this < globest ) {
                   if ( ismult ) {
                       gloglobest = totbest = totscore ; nstored = 0 ; goodmvpt = goodmvs ; success = 1 ; }
                   rounds = 1 ; 
                   besroot = 0 ; 
                   besloc = trylist[ i ] ;
                   globest = this ; 
                   if ( !this && !ismult ) break ;
                   oplimit = this ;
                   if ( ismult && nstored < maxtrees ) ++ oplimit ; }
               if ( totscore == gloglobest && ismult ) {
                   if ( anc [ trylist[ i ] ] == cut ) continue ; 
                   goodmvpt->dest = 0 ;
                   goodmvpt->root = trylist[ i ] ;
                   ++ goodmvpt ; }}
            if ( goodmvpt > goodmvs )
                compare_n_save ( goodmvpt , cut , nitroot ) ; }
        else {
            unzert( cut ) ;
            nlocs = listabove( trylist , sis[ 0 ] , 1 ) ;
            if ( cut < nt ) { nroots = 3 ; rootlist[ 2 ] = cut ; }
            else nroots = listabove( rootlist , cut , 1 ) ;
            fixupdn( anc[ nitloc ] , cut , 0 ) ; 
            if ( dospr ) nroots = 3 ; 
            oplimit = BIGINT ; 
            globest = testamove ( cut , rootlist[ 2 ] , nitloc ) ;
            if ( ismult && nstored < maxtrees ) oplimit = globest + 1 ;
            else oplimit = globest ;
            divsteps = totbest - globest ;
            if ( oplimit ) 
             for ( i = 2 ; i < nroots ; ++ i ) {
               goodmvpt = goodmvs ;  
               for ( j = 0 ; j < nlocs ; ++ j ) {
                   ++ nmoves ; 
                   this = testamove( cut , rootlist[ i ] , trylist[ j ] ) ;
                   totscore = divsteps + this ; 
                   if ( this < globest ) {
                       if ( ismult ) {
                           gloglobest = totbest = totscore ; nstored = 0 ; goodmvpt = goodmvs ; success = 1 ; }
                       rounds = 1 ; 
                       besroot = rootlist[ i ] ; 
                       besloc = trylist[ j ] ;
                       globest = this ; 
                       if ( !this && !ismult ) { i = nroots ; break ; }
                       oplimit = this ;
                       if ( ismult && nstored < maxtrees ) ++ oplimit ; }
                   if ( totscore == gloglobest && ismult ) {
                       if ( trylist[ i ] == nitloc && anc[ rootlist[ i ] ] == cut ) continue ; 
                       goodmvpt->dest = trylist[ j ] ;
                       goodmvpt->root = rootlist [ i ] ;
                       ++ goodmvpt ; }}
               if ( goodmvpt > goodmvs )
                   compare_n_save ( goodmvpt , cut , nitroot ) ; }}
        if ( nitloc != besloc || besroot != nitroot ) ++ nswitches ;
        if ( !nitloc ) {
            if ( nitloc != besloc )
               tbreroot( cut , besloc ) ;
            fichop ( 2 , nt ) ; }
        else {
            inzert( cut , besroot , besloc ) ;
            if ( cut < nt ) fixupdn ( anc[ cut ] , -1 , 0 ) ; 
            else fixupdn ( anc[ cut ] , -1 , nitroot ) ; }}
  if ( !ismult ) globest = fichop( 2 , nt ) ;
  if ( verbose ) {
     if ( !ismult ) {
         if ( dospr ) fprintf ( stderr , "\rCompleted SPR (%i moves, %i accepted), %i steps\n" , nmoves , nswitches , globest ) ; 
         else fprintf ( stderr , "\rCompleted TBR (%i moves, %i accepted), %i steps\n" , nmoves , nswitches , globest ) ; }
     else fprintf ( stderr ,"\rSwapped tree %3i of %3i (%i steps)..." , nowswapping , nstored , gloglobest ) ; }
  return success ; }

/** Functions to handle multiple trees in global branch swapping ************/

void mksis ( void ) {  // N.B. requires bin trees as input...
    int a , b , c , * lassis = oplist , * nodfork = ladd ; 
    for ( a = 2*nt-1 ; a -- ; ) nodfork [ a ] = 0 ;
    for ( a = nt ; a -- ; ) {
      b = anc[ c = a ] ;   
      while ( 1 ) {
        if ( !nodfork [ b ] ) lef [ b ] = lassis [ b ] = c ;
        else sis [ sis [ c ] = lassis [ b ] ] = rig [ b ] = c ;
        if ( ++ nodfork [ b ] > 1 ) break ;
        c = b ;
        if ( b == nt ) break ; 
        b = anc [ b ] ; }}}

void getstord ( int i ) {
    int j ; 
    for ( j = 0 ; j < 2*nt-1 ; ++ j ) anc[ j ] = stortree[ i ] [ j ] & GETANC ;
    mksis  () ;
    anc[ nt ] = 2*nt-1 ; }

void alloc_maxtrees ( void ) {
    int i ; 
    stortree = malloc ( maxtrees * sizeof ( int * ) ) ;
    tkeys = malloc ( maxtrees * sizeof ( int ) ) ;
    if ( stortree == NULL ) errout ( "not enough memory to save that many trees...\n" ) ;
    for ( i = 0 ; i < maxtrees ; ++ i ) {
       stortree [ i ] = malloc ( 2 * nt * sizeof ( int ) ) ;
       if ( stortree [ i ] == NULL ) errout ( "not enough memory to save that many trees...\n" ) ; }}

void fillsupports ( void ) {
    int i , j , isup ;
    for ( i = nt ; i < 2*nt-1 ; ++ i ) {
       if ( i == nt || anc[ i ] == nt ) { supp[ i ] = 1 ; continue ; }
       for ( j = 0 , isup = 0 ; j < nc && !isup ; ++ j ) 
          if ( !( upmatrix[ i ] [ j ] & upmatrix[ anc[ i ] ] [ j ] ) ) isup = 1 ;
       supp[ i ] = isup ; }}   

int resetnodz ( void ) {
  int a , b , x , i , previous , tkey , * ptr0, * ptr1 , * done = oplist , countdown ;
  for( i = nt ; i < 2*nt-1 ; ++ i ) done[ i ] = 0 ;
  ptr1 = tmptbuf ; 
  done[nt] = ptr1[0] = nt ;
  x = nt + 1 ;
  countdown = 2*nt - 2 ;
  for( a = 1 ; a < nt ; ++ a ) {
     previous = b = a ; 
     while ( b != nt ) {
        b = anc [b] ;
        if ( supp [ b ] ) {
          if ( !done[b]) { 
              done[b] = ptr1 [previous] = x ;
              previous = x ++ ; }
          else { 
               ptr1 [previous] = done [ b ] ;
               b = nt ; }} 
        else {
         if( ! done [ b ] ) {
             done [ b ] = ptr1 [previous] = countdown;
             previous = countdown -- ; }
         else { ptr1 [previous] = done [ b ] ;
                 b = nt ; }}}}
  ++ countdown ;
  for( ; countdown < 2*nt-1 ; ++ countdown ) ptr1 [countdown] |= NEGBIT;
  for( a = 0 , ptr0 = tmptbuf , ptr1 = stortree[ nstored ] ; a < 2*nt-1 ; ++ a ) * ptr1 ++ = * ptr0 ++ ;
  ptr1 = stortree[ nstored ] ;
  ptr1 [ nt ] = 2*nt-1 ;
  tkey = 0 ;
  for ( a = 1 ; a < x ; ++ a ) {
      if ( a == nt ) continue ; 
      b = ptr1 [ a ] & GETANC ;
      while ( ( ptr1 [ b ] & NEGBIT ) ) b = ptr1 [ b ] & GETANC ;
      tkey ^= b ; }
  return tkey ; }

int equal_trees( int i, int ii ) {
  int * one , * the_other , a , b , c;
  one = stortree[ i ] ;
  the_other = stortree[ ii ] ;
  for ( a = 1 ; a < 2*nt-1 ; ++ a ) {
      if ( a == nt ) continue ; 
      if( one [ a ] < 0 ) { if ( the_other [a] > 0 ) return 0 ; continue; } 
      if ( the_other [ a ] < 0 ) return 0 ;
      b = one [ a ] & GETANC ;
      c = the_other [ a ] & GETANC ;
      if ( b == nt ) { if ( c != nt ) return 0 ; continue; }
      if ( c == nt ) return 0 ; 
      while ( one [ b ] < 0 ) b = one [ b ] & GETANC ;
      while ( the_other [ c ] < 0 ) c = the_other [ c ] & GETANC ;
      if ( c != b ) return 0 ; }
  return 1 ; }

void storit ( void ) {
   int i , iseq = 0 , curky ;
   if ( nstored == maxtrees ) return ; 
   fillsupports () ; 
   for ( i = 0 , curky = resetnodz () ; i < nstored && !iseq ; ++ i )
     if ( curky == tkeys[ i ] && equal_trees( i , nstored ) ) iseq = 1 ;  
   if ( !iseq ) 
      tkeys[ nstored ++ ] = curky ; }

void slide_down( int cut , int start , int doend ) {
    int i , x , at = start , bak ;
    if ( at < nt ) at = anc[ at ] ; 
    while( 1 ) {
       for ( i = 0 ; i < nc ; ++ i ) {
          x = matrix[ lef[ at ] ][ i ] & matrix[ rig[ at ] ][ i ] ;
          bak = matrix[ lef[ at ] ][ i ] | matrix[ rig[ at ] ][ i ] ;
          if ( !x ) { x = bak ; bak = 1023 ; } 
          matrix[ at ] [ i ] = x ;
          bakmat[ at ][ i ] = bak ; 
          upmatrix[ at ][ i ] = 0 ; }
       if ( at == cut ) break ; 
       at = anc [ at ] ; }}   

void compare_n_save ( GoodMvTyp * goodmvpt , int cut , int nitroot ) {
    GoodMvTyp * atmv = goodmvs - 1 ;
    int dest , root , otroot = sis[ nitroot ] , nitloc , chgdroot ;
    if ( nstored == maxtrees ) return ;
    while ( ++ atmv < goodmvpt ) {
        dest = atmv -> dest ;
        root = atmv -> root ;
        chgdroot = 0 ;
        if ( dest && cut > nt && anc[ root ] != cut ) chgdroot = 1 ; 
        if ( dest ) {
            inzert( cut , root , dest ) ;
            if ( cut < nt ) fixupdn ( anc[ cut ] , -1 , 0 ) ; 
            else {
                if ( chgdroot ) slide_down( cut , nitroot , 0 ) ;  
                fixupdn ( anc[ cut ] , -1 , nitroot ) ; }}
        else {
            tbreroot( cut , root ) ;
            fichop ( 2 , nt ) ; }
        storit () ;
        if ( dest ) {
           nitloc = sis [ cut ] ; 
           unzert( cut ) ; }
        if ( cut > nt ) 
           if ( anc[ otroot ] == nitroot ) tbreroot ( cut , otroot ) ;
           else tbreroot ( cut , nitroot ) ;
        if ( dest ) {
             if ( cut > nt && chgdroot ) slide_down( cut , root , 1 ) ; 
             fixupdn( anc[ nitloc ] , cut , 0 ) ; }
        else fichop ( 2 , sis[ 0 ] ) ; }}
   
void main ( int argc , char ** argv ) {
    int doswap = 1 , i , dotree = 0 ; 
    clock_t initime ; 
    if ( argc == 1 )
        errout ( "\nMust give arguments...\n\n"
                 "     rN    set random seed to N"
                 "     nN    do N starting points\n"
                 "     mN    save multiple trees (no N=100)\n"
                 "     +xxx  swap first tree in file xxx\n"
                 "     t     display tree diagrm(s)\n"
                 "     v     verbose mode\n" 
                 "     /     use SPR swapping instead of TBR\n"
                 "     -     no branch-swapping (i.e. just Wagners)\n\n"
                 "Default is 10 starting points + TBR (nomulpars),\n"
                 "and no multiple trees.\n\n" ) ; 
    if ( ( inf = fopen ( argv[ 1 ] , "rb" ) ) == NULL ) errout ( "Can't open file" ) ;
    readata () ;
    for ( i = 2 ; i < argc ; ++ i ) {
        if ( !strcmp ( argv[ i ] , "v" ) ) verbose = 1 ; 
        if ( argv[ i ] [ 0 ] == 'r' ) rseed = atoi ( argv[ i ] + 1 ) ;
        if ( argv[ i ][ 0 ] == '-' ) doswap = 0 ; 
        if ( argv[ i ][ 0 ] == '/' ) dospr = 1 ;
        if ( argv[ i ][ 0 ] == 'n' ) nrepls = atoi ( argv[ i ]+1 ) ;
        if ( argv[ i ][ 0 ] == 'm' ) {
            multrees = 1 ;
            if ( isdigit ( argv[ i ][ 1 ] ) ) 
              if ( ( maxtrees = atoi ( argv[ i ]+1 ) ) < 1 ) errout ( "maxtrees must be 1 or more!\n" ) ; }
        if ( argv[ i ][ 0 ] == '+' ) { wagstart = 0 ; read_tree ( argv[ i ] + 1 ) ; }
        if ( argv[ i ][ 0 ] == 't' ) dotree = 1 ; }
    alloc_maxtrees () ; 
    initime = clock () ;
    if ( wagstart ) {
        for ( i = 0 ; i < nrepls ; ++ i ) {
           wagner () ;
           if ( doswap ) tbrswap ( 0 ) ;
           if ( globest < gloglobest ) { nstored = 0 ; gloglobest = globest ; }
           if ( globest == gloglobest ) storit () ; }}
    else {
        int len = fichop ( 1 , nt ) ;
        fprintf ( stderr , "Start swapping from user tree (%i steps)\n" , len ) ; 
        if ( doswap ) tbrswap ( 0 ) ;
        gloglobest = globest ; 
        storit () ; }
    fprintf ( stderr , "Overall best score: %i steps (%i trees)\n" , gloglobest , nstored ) ;
    if ( multrees ) {
        int abettertree ;
        fprintf ( stderr , "Doing global swapping...\n" ) ; 
        for ( nowswapping = 0 ; nowswapping < nstored ; ++ nowswapping ) {
           getstord( nowswapping ) ;
           abettertree = tbrswap ( 1 ) ;
           if ( abettertree ) nowswapping = 0 ; }
        fprintf ( stderr , "\rFinished swapping %i trees (%i steps)           \n" , nstored , gloglobest ) ; }
    fprintf ( stderr , "Time used: %.2f sec\n" , ( double ) ( clock() - initime ) / CLOCKS_PER_SEC ) ;
    fflush ( stderr ) ;
    printf ( "\ntread\n'%i trees, %i steps'\n" , nstored , gloglobest ) ; 
    for ( i = 0 ; i < nstored ; ++ i ) {
      getstord( i ) ;   
      if ( i ) printf ( " *\n" ) ; 
      save( nt ) ; }
    printf ( " ;\nproc/;\n" ) ; 
    if ( dotree ) 
      for ( i = 0 ; i < nstored ; ++ i ) {
         getstord( i ) ; 
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
    linsdone = 0 ; 
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
      printf ( "%c" , joint ) ; }
  printf( "%c%c %s\n" , hor , hor , names[ which ] ) ; }  

void fixupdn( int from , int side , int nitroot ) {
    int b , n , nn , x , bak , at , w, blow , was , canbreak , hasbreak = 1 , startat = from ;
    if ( nitroot > 0 ) {
           hasbreak = 0 ;
           startat = nitroot ;
           if ( startat < nt ) startat = anc[ nitroot ] ; } 
    for ( b = 0 ; b < nc ; ++ b ) {
       upmatrix[ at = startat ] [ b ] = 0 ;
       canbreak = hasbreak ; 
       while ( 1 ) {
          bak = matrix[ lef [ at ] ] [ b ] | matrix [ rig [ at ] ] [ b ] ; 
          if ( ( x = matrix[ lef [ at ] ] [ b ] & matrix [ rig [ at ] ] [ b ] ) == 0 ) {
              x = bak ;
              bak = 1023 ; }
          upmatrix[ at ] [ b ] = 0 ;
          bakmat[ at ] [ b ] = bak ;
          if ( x == matrix[ at ] [ b ] && at != from ) 
              if ( canbreak ) break ;
          if ( at == from ) canbreak = 1 ;     
          matrix[ at ] [ b ] = x ;
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
           w = matrix[ at ] [ b ] & upmatrix[ blow ] [ b ] ;
           if ( w == upmatrix[ blow ] [ b ] ) {
              if ( at > nt ) upmatrix [ at ] [ b ] = upmatrix[ blow ] [ b ] ;
              unmat[ at ] [ b ] = upmatrix[ blow ] [ b ] ; }
           else 
              if ( at < nt ) unmat[ at ] [ b ] = upmatrix[ blow ] [ b ] | matrix [ at ] [ b ] ; 
              else {
                  upmatrix[ at ] [ b ] = matrix[ at ] [ b ] | ( upmatrix[ blow ] [ b ] & bakmat[ at ] [ b ] ) ;
                  unmat [ at ] [ b ] = upmatrix[ at ] [ b ] | upmatrix[ blow ] [ b ] ; }
           if ( at >= nt ) 
              if ( upmatrix[ at ] [ b ] != was ) {
                 oplist[ nn ++ ] = lef [ at ] ; 
                 oplist[ nn ++ ] = rig [ at ] ; }}}}

