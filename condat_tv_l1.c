#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <limits.h>
#include <wchar.h>
#include <time.h>

#define EPSILON 1e-10
#define IS_ZERO(x) (x < EPSILON & x > -EPSILON) 
#define IS_POSITIVE(x) (x > EPSILON)
#define IS_NEGATIVE(x) (x < -EPSILON)

int tautString_TV1(double *y,double lambda,double *x,int n) {
    /* Minorant and minorant slopes */
    double mn, mx;
    /* Relative height of majorant and minorant slopes at the current points w.r.t. the tube center */
    double mnHeight, mxHeight;
    /* Last break points of minorant and majorant */
    int mnBreak, mxBreak;
    /* Last break point of taut string */
    int lastBreak;
    /* Auxiliary variables */
    int i, j;
    /* Helpful constants */
    const double minuslambda = -lambda;
    const double lambda2 = 2*lambda;
    const double minuslambda2 = 2*minuslambda;
        
    #define CANCEL(txt,info) \
        printf("tautString_TV1: %s\n",txt); \
        return 0;
    
    /* Starting point */
    mnHeight = mxHeight = 0;
    mn = minuslambda + y[0];
    mx = lambda + y[0];
    lastBreak = -1;
    mnBreak = mxBreak = 0;
    
    #ifdef DEBUG
        fprintf(DEBUG_FILE,"Starting taut-string with length=%d and penalty=%lf\n",n,lambda); fflush(DEBUG_FILE);
    #endif
        
    /* Proceed along string */
    i = 0;
    while ( i < n ) {
        /* Loop over all points except the last one, that needs special care */
        while ( i < n ) {
            #ifdef DEBUG
                fprintf(DEBUG_FILE,"i = %d, mx = %.3f, mn = %.3f, mxHeight = %.3f, mnHeight = %.3f, mxBreak = %d, mnBreak = %d, lastBreak = %d\n",i,mx,mn,mxHeight,mnHeight,mxBreak,mnBreak,lastBreak); fflush(DEBUG_FILE);
            #endif
            
            /* Update height of minorant slope w.r.t. tube center */
            /* This takes into account both the slope of the minorant and the change in the tube center */
            mnHeight += mn - y[i];
        
            /* Check for ceiling violation: tube ceiling at current point is below proyection of minorant slope */
            /* Majorant is r + lambda (except for last point), which is computed on the fly */   
            if ( lambda < mnHeight ) {
                #ifdef DEBUG
                    fprintf(DEBUG_FILE,"CEILING VIOLATION i = %d, mxVal = %f, mnHeight = %f, mxBreak = %d, mnBreak = %d, lastBreak = %d\n",i,lambda,mnHeight,mxBreak,mnBreak,lastBreak); fflush(DEBUG_FILE);
                #endif
                /* Break segment at last minorant breaking point */
                i = mnBreak + 1;
                /* Build valid segment up to this point using the minorant slope */
                for ( j = lastBreak+1 ; j <= mnBreak ; j++ )
                    x[j] = mn;
                /* Start new segment after the break */
                lastBreak = mnBreak;
                /* Build first point of new segment, which can be done in closed form */
                mn = y[i]; 
                mx = lambda2+y[i];
                mxHeight = lambda;
                mnHeight = minuslambda;
                mnBreak = mxBreak = i;
                i++;
                continue;
            }
            
            /* Update height of minorant slope w.r.t. tube center */
            /* This takes into account both the slope of the minorant and the change in the tube center */
            mxHeight += mx - y[i];
            
            /* Check for minorant violation: minorant at current point is above proyection of majorant slope */
            /* Minorant is r - lambda (except for last point), which is computed on the fly */
            if ( minuslambda > mxHeight ) {
                #ifdef DEBUG
                    fprintf(DEBUG_FILE,"FLOOR VIOLATION i = %d, mnVal = %f, mxHeight = %f, mxBreak = %d, mnBreak = %d, lastBreak = %d\n",i,minuslambda,mxHeight,mxBreak,mnBreak,lastBreak); fflush(DEBUG_FILE);
                #endif
                /* If violated, break segment at last majorant breaking point */
                i = mxBreak + 1;
                /* Build valid segment up to this point using the majorant slope */
                for ( j = lastBreak+1 ; j <= mxBreak ; j++ )
                    x[j] = mx;
                /* Start new segment after the break*/
                lastBreak = mxBreak;
                /* Build first point of new segment, which can be done in closed form */
                mx = y[i]; 
                mn = minuslambda2+y[i];
                mxHeight = lambda;
                mnHeight = minuslambda;
                mnBreak = mxBreak = i;
                i++;
                continue;
            }
            
            /* No violations at this point */

            /* Check if proyected majorant height is above ceiling */
            if ( mxHeight >= lambda ) {
                /* Update majorant slope */
                mx += ( lambda - mxHeight ) / ( i - lastBreak );
                /* Get correct majorant height (we are touching it!) */
                mxHeight = lambda;
                /* This is a possible majorant breaking point */
                mxBreak = i;
            }
            
            /* Check if proyected minorant height is under actual minorant */
            if ( mnHeight <= minuslambda ) {
                /* Update minorant slope */
                mn += ( minuslambda - mnHeight ) / ( i - lastBreak );
                /* Compute correct minorant height (we are touching it!) */
                mnHeight = minuslambda;
                /* This is a possible minorant breaking point */
                mnBreak = i;
            }

            #ifdef DEBUG
                fprintf(DEBUG_FILE,"i = %d, mx = %.3f, mn = %.3f, mxHeight = %.3f, mnHeight = %.3f, mxBreak = %d, mnBreak = %d, lastBreak = %d\n",i,mx,mn,mxHeight,mnHeight,mxBreak,mnBreak,lastBreak); fflush(DEBUG_FILE);
            #endif
            
            /* At this point: no violations, so keep up building current segment */
            i++;
        }
    }
    
    /* Build last valid segment */
    for ( i = lastBreak+1 ; i < n ; i++ )
        x[i] = mn;
        
    /* Return */
    return 1;
    
    #undef CANCEL
}