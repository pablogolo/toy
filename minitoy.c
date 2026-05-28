#include <stdio.h>
#include <stdlib.h>
FILE * inf ; 
#define MAXT 1000
#define MAXC 1500
#define MISSING 1023 
#define TBR_MOVE( x , y , z ) \
    { if ( ( this = testamove( x , y , z ) ) < globscore ) { \
         rounds = 1 ;                                        \
         besroot = y ;                                       \
         besloc = z ;                                        \
         globscore = this ; }}
int nt , nc , matrix[ 2*MAXT ][ MAXC ] , unmat[ 2*MAXT ][ MAXC ] , bakmat[ 2*MAXT ][ MAXC ] , inmask[256] , globest = 10e9 , dospr = 0 , wagstart = 1 , curnodin , rseed = 1 , globscore ; 
int anc[2*MAXT] , lef[2*MAXT] , rig[2*MAXT] , sis[2*MAXT] , oplist[2*MAXT] , trylist[2*MAXT] , rootlist[2*MAXT] , ladd[2*MAXT] ;
char names[ MAXT ] [ 50 ] ;

void save ( int i ) {
    if ( i < nt ) { printf ( "%i " , i ) ; return ; }
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
    for ( i = 0 ; i < 10 ; ++ i ) inmask[ '0' + i ] = 1 << i ;
    inmask[ '?' ] = inmask[ '-' ] = MISSING  ;
    inmask[ '[' ] = 0 ; }
    
void readata ( char * filename ) {
    int a , b , ispoly ;
    char ch ;
    if ( ( inf = fopen ( filename , "rb" ) ) == NULL ) { fprintf ( stderr , "Data file?" ) ; exit ( 0 ) ; } 
    setmask () ;
    fscanf ( inf , " %i %i" , &nt , &nc ) ;
    for ( a = 0 ; a < nt ; ++ a ) {
       fscanf ( inf , " %s" , &names[ a ] ) ;
       for ( b = 0 ; b < nc ; ++ b ) {
           fscanf ( inf , " %c" , &ch ) ;
           ispoly = 0 ; 
           if ( ch == '[' ) { ispoly = 1 ; fscanf( inf , " %c" , &ch ) ; }
           matrix[ a ] [ b ] = 0 ; 
           while ( ch != ']' ) {
               if ( !inmask[ ( int ) ch ] ) { fprintf ( stderr , "Bad symbol" ) ; exit ( 0 ) ; }
               matrix[ a ] [ b ] |= inmask[ ( int ) ch ] ;
               if ( !ispoly ) ch = ']' ;
               else fscanf ( inf , " %c" , &ch ) ; }
          if ( !a ) unmat[ a ] [ b ] = matrix[ a ] [ b ] ; }}}

int listabove( int * list , int from , int doterms ) {
    int * done = list , * todo = list + 1 , at ;
    * done = from ;
    while ( done < todo ) 
       if ( ( at = * done ++ ) >= nt ) {
          if ( lef [ at ] >= nt || doterms ) * todo ++ = lef [ at ] ; 
          if ( rig [ at ] >= nt || doterms ) * todo ++ = rig [ at ] ; }
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
        for ( b = 0 ; b < nc ; ++ b )
           unmat [ at ] [ b ] = matrix[ at ] [ b ] ; 
        for ( i = 1 ; i < n ; ++ i ) 
          for ( b = 0 , blow = anc [ at = oplist[ i ] ] ; b < nc ; ++ b ) {
            if ( !at ) continue ;   
            if ( ( w = matrix[ at ] [ b ] & matrix[ blow ] [ b ] ) == matrix[ blow ] [ b ] ) {
               if ( at > nt ) matrix [ at ] [ b ] = matrix[ blow ] [ b ] ;
               unmat[ at ] [ b ] = matrix[ blow ] [ b ] ; }
            else 
               if ( at < nt ) unmat[ at ] [ b ] = matrix[ blow ] [ b ] | matrix [ at ] [ b ] ; 
               else {
                   matrix[ at ] [ b ] |= matrix[ blow ] [ b ] & bakmat[ at ] [ b ] ;
                   unmat [ at ] [ b ] = matrix[ at ] [ b ] | matrix[ blow ] [ b ] ; }}}
    return steps ; }        

int testamove ( int what , int root , int where ) {
    int b , steps = 0 ;
    for ( b = 0 ; b < nc && steps < globscore ; ++ b ) 
        if ( what == root && what < nt ) { if ( ! ( matrix[ root ] [ b ] & unmat[ where ] [ b ] ) ) ++ steps ; }
        else if ( ! ( unmat[ root ] [ b ] & unmat[ where ] [ b ] ) ) ++ steps ;
    return steps ; }

void wagner( void ) {
   int i , n , q , now , next = 3 , this , besloc , curnod = nt + 2 , * lp = ladd + 1 , left = nt - 1 ; 
   for ( i = 0 ; i < nt ; ++ i ) ladd [ i ] = i ;
   srand ( rseed ) ;
   for ( i = 1 ; i < nt ; ++ i ) {
      q = * lp ; 
      * lp = lp [ n = rand () % left -- ] ;
      lp [ n ] = q ;
      ++ lp ; }
   anc[ lef [ nt ] = 0 ] = anc [ rig[ nt ] = nt + 1 ] = nt ;
   anc[ lef[ nt + 1 ] = ladd [ 1 ] ] = anc[ rig[ nt + 1 ] = ladd [ 2 ] ] = nt + 1 ;
   sis [ sis[ 0 ] = nt + 1 ] = 0 ;
   sis [ sis [ ladd [ 1 ] ] = ladd [ 2 ] ] = ladd [ 1 ] ; 
   while ( next < nt ) {
       anc[ now = ladd [ next ++ ] ] = curnod ++ ;
       globscore = 10e9 ;
       fichop ( 2 , nt ) ; 
       for ( i = listabove( trylist , sis[ 0 ] , 1 ) ; i -- && globscore ; ) 
          if ( ( this = testamove( now , now , trylist[ i ] ) ) < globscore ) {
              globscore = this ; 
              besloc = trylist[ i ] ; }
       inzert( now , -1 , besloc ) ; }}    

void unzert( cut ) {
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
    int i , j , rounds , cut = 0 , besroot , besloc , this , nitroot , nitloc , nroots , nlocs ; 
    for ( rounds = 1 ; rounds < 2*nt - 1 ; ++ rounds ) {
        if ( ++ cut > 2*nt-2 ) cut = 1 ;
        if ( cut == nt ) continue ;
        if ( cut > nt ) nitroot = besroot = lef [ cut ] ;
        else nitroot = besroot = cut ;
        globscore = 10e9 ; 
        if ( ! ( nitloc = besloc = sis[ cut ] ) ) {
            fichop ( 2 , sis[ 0 ] ) ;
            globscore = testamove( 0 , 0 , cut ) ;
            for ( i = listabove( trylist , cut , 1 ) ; -- i && globscore ; )
                 TBR_MOVE( 0 , 0 , trylist[ i ] ) ; }
        else {
            unzert( cut ) ;
            fichop ( 2 , nt ) ;
            nlocs = listabove( trylist , sis[ 0 ] , 1 ) ;
            if ( cut < nt ) { nroots = 3 ; rootlist[ 2 ] = cut ; }
            else { fichop ( 2 , cut ) ; nroots = listabove( rootlist , cut , 1 ) ; } 
            globscore = testamove ( cut , cut , nitloc ) ;
            for ( i = 2 ; i < nroots && globscore ; ++ i )
               for ( j = 0 ; j < nlocs && globscore ; ++ j )
                   TBR_MOVE( cut , rootlist[ i ] , trylist[ j ] ) ; }
        if ( !nitloc ) { if ( nitloc != besloc ) tbreroot( cut , besloc ) ; }
        else inzert( cut , besroot , besloc ) ; }
  fprintf ( stderr , "\rCompleted Wagner + TBR (score %i)\n" , fichop( 1 , nt ) ) ; }

void main ( int argc , char ** argv ) {
    if ( argc > 2 ) rseed = atoi ( argv[ 2 ] ) ; 
    readata ( argv[ 1 ] ) ;
    wagner () ; 
    tbrswap () ; 
    save( nt ) ; } 
