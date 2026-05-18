/***************************************************\
|       A Toy program for phylogenetic analysis     |
|            by Pablo A. Goloboff (2026)            |
|                                                   |
|           Intended for pedagogic use              |
|           Feel free to use or modify,             |
|        just acknowledge authorsip if you do.      |
\***************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <conio.h>
#include <unistd.h> 
FILE * inf ; 
#define MAXT 1000
#define MAXC 1500
#define MISSING 1023 
int nt , nc , matrix[ 2*MAXT ][ MAXC ] , unmat[ 2*MAXT ][ MAXC ] , bakmat[ 2*MAXT ][ MAXC ] , inmask[256] , globest = 10e9 ; 
int anc[2*MAXT] , lef[2*MAXT] , rig[2*MAXT] , sis[2*MAXT] , oplist[2*MAXT] , trylist[2*MAXT] , rootlist[2*MAXT] , ladd[2*MAXT] ;
char names[ MAXT ] [ 50 ] ;
unsigned long int rseed = 1 ;

void errout ( char * txt ) {
    fprintf ( stderr , "Can't run:\n    %s\n(press key to exit)\n" , txt ) ;
    getch() ; 
    exit ( 0 ) ; 
}

void save ( int i ) {
    if ( i < nt ) {
        printf ( "%i " , i ) ;
        return ; }
    if ( i == nt ) printf ( "\ntread\n" ) ;
    printf( "(" ) ;
    save( lef[ i ] ) ; 
    save( rig [ i ] ) ; 
    printf( ")" ) ;
    if ( i == nt ) printf ( ";\nproc/;" ) ;
}

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
   return ;
}

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
    sis [ sis [ dasis ] = join ] = dasis ; 
}    

void setmask( void ) {
    int i ;
    for ( i = 0 ; i < 256 ; ++ i ) inmask[ i ] = 0 ;
    for ( i = 0 ; i < 10 ; ++ i ) inmask[ '0' + i ] = 1 << i ;
    inmask[ '?' ] = inmask[ '-' ] = MISSING  ;
    inmask[ '[' ] = 1 ; 
}
    
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
          if ( !a ) unmat[ a ] [ b ] = matrix[ a ] [ b ] ; }}
}

int listabove( int * list , int from , int doterms ) {
    int * done , * todo , at ;
    * ( done = list ) = from ;
    todo = done + 1 ;
    while ( done < todo ) {
       at = * done ++ ;
       if ( at >= nt ) {
          if ( lef [ at ] >= nt || doterms ) * todo ++ = lef [ at ] ; 
          if ( rig [ at ] >= nt || doterms ) * todo ++ = rig [ at ] ; }}
    return done - list ;        
}

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
        for ( b = 0 ; b < nc ; ++ b )
           unmat [ at ] [ b ] = matrix[ at ] [ b ] ; 
        for ( i = 1 ; i < n ; ++ i ) 
          for ( b = 0 , blow = anc [ at = oplist[ i ] ] ; b < nc ; ++ b ) {
            if ( !at ) continue ;   
            w = matrix[ at ] [ b ] & matrix[ blow ] [ b ] ;
            if ( w == matrix[ blow ] [ b ] ) {
               if ( at > nt ) matrix [ at ] [ b ] = matrix[ blow ] [ b ] ;
               unmat[ at ] [ b ] = matrix[ blow ] [ b ] ; }
            else 
               if ( at < nt ) unmat[ at ] [ b ] = matrix[ blow ] [ b ] | matrix [ at ] [ b ] ; 
               else {
                   matrix[ at ] [ b ] |= matrix[ blow ] [ b ] & bakmat[ at ] [ b ] ;
                   unmat [ at ] [ b ] = matrix[ at ] [ b ] | matrix[ blow ] [ b ] ; }}}
    return steps ;    
}        

int testamove ( int what , int root , int where ) {
    int b , steps = 0 ;
    for ( b = 0 ; b < nc && steps < globest ; ++ b ) 
        if ( what == root && what < nt ) {
            if ( ! ( matrix[ root ] [ b ] & unmat[ where ] [ b ] ) ) ++ steps ; }
        else 
            if ( ! ( unmat[ root ] [ b ] & unmat[ where ] [ b ] ) ) ++ steps ;
    return steps ; 
}

void randomize( void ) {
   int i , n , q , * lp = ladd + 1 , left = nt - 1 ; 
   for ( i = 0 ; i < nt ; ++ i ) ladd [ i ] = i ;
   if ( !rseed ) return ; 
   srand ( rseed ) ;
   for ( i = 1 ; i < nt ; ++ i ) {
      q = * lp ; 
      * lp = lp [ n = rand () % left -- ] ;
      lp [ n ] = q ;
      ++ lp ; }
}

void nitnet( void ) {
    randomize () ;
    anc[ lef [ nt ] = 0 ] = anc [ rig[ nt ] = nt + 1 ] = nt ;
    anc[ lef[ nt + 1 ] = ladd [ 1 ] ] = anc[ rig[ nt + 1 ] = ladd [ 2 ] ] = nt + 1 ;
    sis [ sis[ 0 ] = nt + 1 ] = 0 ;
    sis [ sis [ ladd [ 1 ] ] = ladd [ 2 ] ] = ladd [ 1 ] ; 
}
void wagner( void ) {
   int i , now , at , next = 3 , this , besloc , curnod = nt + 2 ; 
   nitnet () ;
   while ( next < nt ) {
       anc[ now = ladd [ next ++ ] ] = curnod ++ ;
       globest = 10e9 ;
       fichop ( 2 , nt ) ; 
       for ( i = listabove( trylist , sis[ 0 ] , 1 ) ; i -- ; ) {
          at = trylist[ i ] ;
          this = testamove( now , now , at ) ;
          if ( this < globest ) {
              besloc = at ;
              if ( !( globest = this ) ) break ; }}
       inzert( now , -1 , besloc ) ; }
   globest = fichop ( 1 , nt ) ; 
   fprintf ( stderr , "Formed a Wagner tree of %i steps\n" , globest ) ; 
}    

void unzert( cut ) {
    int danc = anc[ cut ] , ddanc , brother , cousin ;
    if ( danc == nt ) return ; 
    ddanc = anc [ danc ] ;
    brother = sis[ cut ] ;
    cousin = sis [ danc ] ;
    lef[ ddanc ] = brother ;
    rig [ ddanc ] = cousin ;
    sis [ sis [ cousin ] = brother ] = cousin ;
    anc [ brother ] = anc [ cousin ] = ddanc ;
}    

void tbrswap( void ) {
    int i , j , rounds , cut = 0 , besroot , besloc , this , nitroot , nitloc , nroots , nlocs , lasreport = 0 , nmoves = 0 , nswitches = 0 ; 
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
            globest = 10e9 ;
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
            fichop ( 2 , nt ) ;
            if ( cut >= nt ) fichop ( 2 , cut ) ;
            globest = 10e9 ; 
            globest = testamove ( cut , cut , nitloc ) ;
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
               tbreroot( cut , besloc ) ; }
        else 
            inzert( cut , besroot , besloc ) ; }
  globest = fichop( 1 , nt ) ; 
  fprintf ( stderr , "\rCompleted TBR (%i moves, %i accepted), best score %i\n" , nmoves , nswitches , globest ) ; 
}

void main ( int argc , char ** argv ) {
    int doswap = 1 ; 
    clock_t initime ; 
    if ( argc == 1 ) errout ( "Give file name" ) ;
    if ( ( inf = fopen ( argv[ 1 ] , "rb" ) ) == NULL ) errout ( "Can't open file" ) ;
    if ( argc > 2 && argv[ 2 ] [ 0 ] != '-' ) 
        rseed = atoi ( argv[ 2 ] ) ;
    if ( argv[ argc - 1 ][ 0 ] == '-' ) doswap = 0 ;     
    readata () ;
    initime = clock () ; 
    wagner () ;
    if ( doswap ) tbrswap () ; 
    fprintf ( stderr , "Time used: %.2f sec\n" , ( double ) ( clock() - initime ) / CLK_TCK ) ;
    fflush ( stderr ) ; 
    save( nt ) ;
}    

