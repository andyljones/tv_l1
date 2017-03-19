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
    
    /* Starting point */
    mnHeight = mxHeight = 0;
    mn = minuslambda + y[0];
    mx = lambda + y[0];
    lastBreak = -1;
    mnBreak = mxBreak = 0;
    
    /* Proceed along string */
    i = 0;
    while ( i < n ) {

        /* Loop over all points except the last one, that needs special care */
        while ( i < n-1 ) {

            /* Update height of minorant slope w.r.t. tube center */
            /* This takes into account both the slope of the minorant and the change in the tube center */
            mnHeight += mn - y[i];
        
            /* Check for ceiling violation: tube ceiling at current point is below proyection of minorant slope */
            /* Majorant is r + lambda (except for last point), which is computed on the fly */   
            if ( lambda < mnHeight ) {
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
            
            /* At this point: no violations, so keep up building current segment */
            i++;
        }
        
        /* Special case i == n-1 (last point) */
        /* We try to validate the last segment, and if we can, we are finished */
        /* The code is essentially the same as the one for the general case, 
           the only different being that here the tube ceiling and floor are both 0 */
        
        /* Update height of minorant slope w.r.t. tube center */
        /* This takes into account both the slope of the minorant and the change in the tube center */
        mnHeight += mn - y[i];
    
        /* Check for ceiling violation: tube ceiling at current point is below proyection of minorant slope */
        /* Majorant is 0 at this point */   
//        if ( IS_POSITIVE(mnHeight) ) { // 0 < mnHeight
//            /* Break segment at last minorant breaking point */
//            i = mnBreak + 1;
//            /* Build valid segment up to this point using the minorant slope */
//            for ( j = lastBreak+1 ; j <= mnBreak ; j++ )
//                x[j] = mn;
//            /* Start new segment after the break */
//            lastBreak = mnBreak;
//            /* Go back to main loop, starting a new segment */
//            /* We do not precompute the first point of the new segment here, as it might be n-1 and this leads to issues */
//            mn = y[i]; 
//            mx = lambda2+y[i];
//            mxHeight = mnHeight = minuslambda;
//            mnBreak = mxBreak = i;
//            continue;
//        }
        
        /* Update height of minorant slope w.r.t. tube center */
        /* This takes into account both the slope of the minorant and the change in the tube center */
        mxHeight += mx - y[i];
        
        /* Check for minorant violation: minorant at current point is above proyection of majorant slope */
        /* Minorant is 0 at this point */
//        if ( IS_NEGATIVE(mxHeight) ) { // 0 > mxHeight
//            /* If violated, break segment at last majorant breaking point */
//            i = mxBreak + 1;
//            /* Build valid segment up to this point using the majorant slope */
//            for ( j = lastBreak+1 ; j <= mxBreak ; j++ )
//                x[j] = mx;
//            /* Start new segment after the break*/
//            lastBreak = mxBreak;
//            /* Go back to main loop, starting a new segment */
//            /* We do not precompute the first point of the new segment here, as it might be n-1 and this leads to issues */
//            mx = y[i]; 
//            mn = minuslambda2+y[i];
//            mxHeight = mnHeight = lambda;
//            mnBreak = mxBreak = i;
//            continue;
//        }
        
        /* No violations at this point */
        
        /* Check if proyected minorant height is under actual minorant */
//        if ( mnHeight <= 0 ) {
//            /* Update minorant slope */
//            mn += ( - mnHeight ) / ( i - lastBreak );
//        }

        
        /* At this point: we are finished validating last segment! */
        i++;
    }
    
    /* Build last valid segment */
    for ( i = lastBreak+1 ; i < n ; i++ )
        x[i] = mn;
        
    /* Return */
    return 1;
    
    #undef CANCEL
}


void prox_dp(double *y, double lam, double *beta, int n) {
  // Take care of a few trivial cases
  if (n==0) return;
  if (n==1 || lam==0) {
    for (int i=0; i<n; i++) beta[i] = y[i];
    return;
  }
  
  // These are used to store the derivative of the
  // piecewise quadratic function of interest
  double afirst, alast, bfirst, blast;
  double *x = (double*)malloc(2*n*sizeof(double));
  double *a = (double*)malloc(2*n*sizeof(double));
  double *b = (double*)malloc(2*n*sizeof(double));
  int l,r;

  // These are the knots of the back-pointers
  double *tm = (double*)malloc((n-1)*sizeof(double));
  double *tp = (double*)malloc((n-1)*sizeof(double));

  // We step through the first iteration manually
  tm[0] = -lam+y[0];
  tp[0] = lam+y[0];
  l = n-1;
  r = n;
  x[l] = tm[0];
  x[r] = tp[0];
  a[l] = 1;
  b[l] = -y[0]+lam;
  a[r] = -1;
  b[r] = y[0]+lam;
  afirst = 1;
  bfirst = -lam-y[1];
  alast = -1;
  blast = -lam+y[1];

  // Now iterations 2 through n-1
  int lo, hi;
  double alo, blo, ahi, bhi;
  for (int k=1; k<n-1; k++) {
    // Compute lo: step up from l until the
    // derivative is greater than -lam
    alo = afirst;
    blo = bfirst;
    for (lo=l; lo<=r; lo++) {
      printf("lo %d\n", lo);
      if (alo*x[lo]+blo > -lam) break;
      alo += a[lo];
      blo += b[lo];
    }

    // Compute the negative knot
    tm[k] = (-lam-blo)/alo;
    l = lo-1;
    printf("l %d\n", l);
    x[l] = tm[k];

    // Compute hi: step down from r until the
    // derivative is less than lam
    ahi = alast;
    bhi = blast;
    for (hi=r; hi>=l; hi--) {
      printf("hi %d\n", hi);
      if (-ahi*x[hi]-bhi < lam) break;
      ahi += a[hi];
      bhi += b[hi];
    }

    // Compute the positive knot
    tp[k] = (lam+bhi)/(-ahi);
    r = hi+1;
    printf("r %d\n", r);
    x[r] = tp[k];

    // Update a and b
    a[l] = alo;
    b[l] = blo+lam;
    a[r] = ahi;
    b[r] = bhi+lam;
    afirst = 1;
    bfirst = -lam-y[k+1];
    alast = -1;
    blast = -lam+y[k+1];
  }

  // Compute the last coefficient: this is where 
  // the function has zero derivative

  alo = afirst;
  blo = bfirst;
  for (lo=l; lo<=r; lo++) {
    if (alo*x[lo]+blo > 0) break;
    alo += a[lo];
    blo += b[lo];
  }
  beta[n-1] = -blo/alo;

  // Compute the rest of the coefficients, by the
  // back-pointers
  for (int k=n-2; k>=0; k--) {
    if (beta[k+1]>tp[k]) beta[k] = tp[k];
    else if (beta[k+1]<tm[k]) beta[k] = tm[k];
    else beta[k] = beta[k+1];
  }

  // Done! Free up memory
  free(x);
  free(a);
  free(b);
  free(tm);
  free(tp);
}
