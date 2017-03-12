#include <math.h>
#include <stdio.h>
#include <stdlib.h>

void TV1D_denoise(double* input, double* output, const int width, const double lambda) {
    if (width>0 && lambda>=0) {                /*to avoid invalid memory access to input[0] and invalid lambda values*/
        int k=0, k0=0;            /*k: current sample location, k0: beginning of current segment*/
        double umin=lambda, umax=-lambda;    /*u is the dual variable*/
        double vmin=input[0]-lambda, vmax=input[0]+lambda;    /*bounds for the segment's value*/
        int kplus=0, kminus=0;     /*last positions where umax=-lambda, umin=lambda, respectively*/
        const double twolambda=2.0*lambda;    /*auxiliary variable*/
        const double minlambda=-lambda;        /*auxiliary variable*/
        for (;;) {                /*simple loop, the exit test is inside*/
            while (k==width-1) {    /*we use the right boundary condition*/
                if (umin<0.0) {            /*vmin is too high -> negative jump necessary*/
                    do output[k0++]=vmin; while (k0<=kminus);
                    umax=(vmin=input[kminus=k=k0])+(umin=lambda)-vmax;
                } else if (umax>0.0) {    /*vmax is too low -> positive jump necessary*/
                    do output[k0++]=vmax; while (k0<=kplus);
                    umin=(vmax=input[kplus=k=k0])+(umax=minlambda)-vmin;
                } else {
                    vmin+=umin/(k-k0+1); 
                    do output[k0++]=vmin; while(k0<=k); 
                    return;
                }
            }
            if ((umin+=input[k+1]-vmin)<minlambda) {        /*negative jump necessary*/
                do output[k0++]=vmin; while (k0<=kminus);
                vmax=(vmin=input[kplus=kminus=k=k0])+twolambda;
                umin=lambda; umax=minlambda;
            } else if ((umax+=input[k+1]-vmax)>lambda) {    /*positive jump necessary*/
                do output[k0++]=vmax; while (k0<=kplus);
                vmin=(vmax=input[kplus=kminus=k=k0])-twolambda;
                umin=lambda; umax=minlambda;
            } else {     /*no jump necessary, we continue*/
                k++;
                if (umin>=lambda) {        /*update of vmin*/
                    vmin+=(umin-lambda)/((kminus=k)-k0+1);
                    umin=lambda;
                } 
                if (umax<=minlambda) {    /*update of vmax*/
                    vmax+=(umax+lambda)/((kplus=k)-k0+1);
                    umax=minlambda;
                }     
            }
        }
    }
}