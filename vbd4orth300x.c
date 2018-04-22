/* ==========================================================================
%  Diffusion on 3-sphere, with drift variability, pools within 16 orthants.
% Vdb4orth300 has bias on target-absent orthants.   9/7/17
%  Separate across-trial drift variability parameters.
%  Modified loop starting points to equate mass in orthants 7/7/17.
%  ntheta must be divisible by 4 for this to work.
%
%  [T, Gt, Ptheta, Mt] = ved4orth300(P, tmax, badix);
%   P = [v1, v2, v3, v4, eta1, eta2, eta3, eta4, sigma, a]
%
%  Building:
%            
%  [T, Gt, Ptheta, Mt] = vbd4orth300([1.0,1.0,0,0,1.0, 0.9], 1.5, 5);
%  
%  Ptheta scaled so mass sums to 1.0, Gt scaled so mass in angle sums to 1.0,
%  still density in time [sum(sum(sum(sum(Gt)))) * h = 1.0].
% ===========================================================================
*/
#include <mex.h>
#include <math.h>
#include <gsl/gsl_sf_bessel.h>

#define kmax 50             /* Maximum number of eigenvalues in dhamana */
#define ntheta 36           /* 0 - 2 * pi */    /**** coarse grained */
#define nphi (ntheta / 2)   /* 0 - pi */
#define npsi (ntheta / 2)   /* 0 - pi */
#define sz 300              /* Number of time steps */
#define north 16
#define nret (north * sz) 
#define nw (npsi * nphi * ntheta)
#define nall (npsi * nphi * ntheta * sz) 
#define NP 11     /* Number of input parameters */
#define BIAS 10   /* Index of bias parameter in parameter vector */
 
const double pi = 3.141592653589793;
double Pijk[nw], Mijk[nw], J1k[kmax], J1k_squared[kmax], J2k[kmax], Commonscale[sz], 
       Gt0[sz], Gtijk[nall], Gtlo[nret], Gthi[nret], Pthetalo[north], Pthetahi[north],
       Mtlo[north], Mthi[north], h;

double lambdac(double c) {
    /* bias angular offset on (-pi / 4, +pi / 4) from c on (-inf, +inf) */
    const double a = 2.0;  /* Steepness of logistic */
                 /* range = pi * (1.0 - 2.0 / nw); - Compress because of discretisation (not needed) */
    double lambda;

    /* Express bias in discrete steps */
    lambda = (1 - exp(-a * c)) / (1 + exp(-a * c)) / 4.0;
    lambda = lambda * ntheta;
    /* mexPrintf("c = %6.3f lambda = %6.3f\n", c, lambda);  */
    return lambda;
}


void d4hamana(double *T, double *Gt, double *P, double h, int badix) {
    /* zero-drift four-dimensional diffusion */ 
    double a, a2, sigma, sigma2, scaler; 
    int i, k;

    a = P[0];
    sigma = P[1];
    sigma2 = sigma * sigma;
    a2 = a * a;
    scaler = sigma2 / a2 / 2.0;
    
    T[0] = 0;
    Gt[0] = 0;
   
    /* Roots of J1k */
    for (k = 0; k < kmax; k++) { 
        J1k[k] = gsl_sf_bessel_zero_J1(k+1);
    }

    /* Evaluate Bessel function at the roots */
    for (k = 0; k < kmax; k++) {
        J1k_squared[k] = J1k[k] * J1k[k];
        /* J1k[k] = j1(J0k[k]);  */ /* Old Gnu C library besselj */
        J2k[k] = gsl_sf_bessel_Jn(2, J1k[k]);  /* Gnu scientific library besselj1 */

    }
    for (i = 1; i < sz; i++) {
        Gt0[i] = 0;    
        T[i] = i * h;
        for (k = 0; k < kmax; k++) {
            Gt[i] += J1k_squared[k] * exp(-J1k_squared[k] * sigma2 * T[i] /(2.0 * a2)) / J2k[k]; 
        }
        Gt[i] *= scaler;
        if (i <= badix || Gt[i] < 0) {
            Gt[i] = 0;
        }
    }
   /* mass = 0;
   for (i = 1; i < sz; i++) {
        mass += Gt[i];
   }
   mass *= h;
   mexPrintf("Mass = %6.3f\n", mass); */
}


void ved4cirle300(double *T, double *Gtijk, double *Pijk, double *Mijk, 
                             double *P, double h, int badix) {
    /* -----------------------------------------------------------------------------------------------
       Calculate first-passage-time density and response probabilities for circular diffusion process
       Surface area of 3-sphere is 2 * pi^2 * r^3 ~ 19.739
       Spherical volume element is sin^2(psi) * sin(phi) * r^3 
      ------------------------------------------------------------------------------------------------ */
    const double eps = 1e-4;
    const double S3area = 2.0 * pi * pi;
    double P0[2];
    double dw, h2, v1, v2, v3, v4, eta1, eta2, eta3, eta4, sigma, a, sigma2, 
           tscale, theta, phi, psi, vol_ijk, dw_cubed,
           eta1onsigma2, eta2onsigma2, eta3onsigma2, eta4onsigma2,
           G11, G21, G31, G41, G12, G22, G32, G42, Girs1, Girs2, Girs3, Girs4, gijkl;

    int i, j, k, l, ijk, ijkl, ijkl_minus_one;

    dw = 2.0 * pi / ntheta;  /* Angular increment, same for all dimensions */
    dw_cubed = dw * dw * dw;  
    /* Parameters */
    h2 = h / 2.0;
    v1 = P[0];
    v2 = P[1];
    v3 = P[2];
    v4 = P[3];
    eta1 = P[4];
    eta2 = P[5];
    eta3 = P[6];
    eta4 = P[7];
    sigma = P[8];
    a = P[9];
    /* mexPrintf("w = %6.4f h =%6.4f \n", w, h);*/
    /* mexPrintf("v1 = %6.4f v2 =%6.4f v3 = %6.4f v4 = %6.4f\n", v1, v2, v3, v4);  */
    if (eta1 < eps) {
        eta1 = eps;
    }
    if (eta2 < eps) {
        eta2 = eps;
    }
    if (eta3 < eps) {
        eta3 = eps;
    }
    if (eta4 < eps) {
        eta4 = eps;
    }

    sigma2 = sigma * sigma;
    eta1onsigma2 = eta1 * eta1 / sigma2;
    eta2onsigma2 = eta2 * eta2 / sigma2;
    eta3onsigma2 = eta3 * eta3 / sigma2;
    eta4onsigma2 = eta4 * eta4 / sigma2;
    P0[0] = a;
    P0[1] = sigma;

    /* Density of zero-drift process */
    d4hamana(T, Gt0, P0, h, badix); 


   /* Joint RT distribution on 3-sphere - outer sum is time, ranges are 0-pi, 0-pi, 0-2*pi  */
   for (l = 0; l < sz; l++) {
       tscale = sqrt((1 / (1 + eta1onsigma2 * T[l])) * (1 / (1 + eta2onsigma2 * T[l])) *
                     (1 / (1 + eta3onsigma2 * T[l])) * (1 / (1 + eta4onsigma2 * T[l])));
       G11 = 2 * eta1 * eta1 * (1 + eta1onsigma2 * T[l]);   
       G21 = 2 * eta2 * eta2 * (1 + eta2onsigma2 * T[l]);   
       G31 = 2 * eta3 * eta3 * (1 + eta3onsigma2 * T[l]);   
       G41 = 2 * eta4 * eta4 * (1 + eta4onsigma2 * T[l]);   
       psi= dw / 2;
       for (i = 0; i < npsi; i++) {
           phi =  dw / 2; 
           for (j = 0; j < nphi; j++) { 
                theta = dw;
                for (k = 0; k < ntheta; k++) {
                     /* Angle-dependent part of Girsanov transformation in hyperspherical coordinates */
                     G12 = v1 + a * eta1onsigma2 * cos(psi);
                     G22 = v2 + a * eta2onsigma2 * sin(psi) * cos(phi);
                     G32 = v3 + a * eta3onsigma2 * sin(psi) * sin(phi) * cos(theta);
                     G42 = v4 + a * eta4onsigma2 * sin(psi) * sin(phi) * sin(theta);
                     Girs1 = exp((G12 * G12) / G11 - (v1 * v1) / (eta1 * eta1) / 2.0);
                     Girs2 = exp((G22 * G22) / G21 - (v2 * v2) / (eta2 * eta2) / 2.0);
                     Girs3 = exp((G32 * G32) / G31 - (v3 * v3) / (eta3 * eta3) / 2.0);
                     Girs4 = exp((G42 * G42) / G41 - (v4 * v4) / (eta4 * eta4) / 2.0);
                     ijk = (nphi * npsi) * k + (npsi * j) + i;
                     ijkl = (ntheta * nphi * npsi) * l + (nphi * npsi) * k + npsi * j + i;
                     vol_ijk = sin(psi) * sin(psi) * sin(phi) * dw_cubed;
                     Gtijk[ijkl] = tscale * Girs1 * Girs2 * Girs3 * Girs4 * Gt0[l] * vol_ijk / S3area;
                     theta += dw;   
                }
                phi += dw;   
            }
            psi += dw;
        }
    }
    /* Integrate joint densities to get hitting probabilities and mean RTs */
    psi = dw / 2; 
    for (i = 0; i < npsi; i++) {
        phi = dw / 2;     
        for (j = 0; j < nphi; j++) {
            for (k = 0; k < ntheta; k++) {
                 ijk = (nphi * npsi) * k + (npsi * j) + i;
                 Pijk[ijk] = 0;
                 Mijk[ijk] = 0;
                 /* vol_ijk = sin(psi) * sin(psi) * sin(phi) * dw_cubed; */
                 for (l = 1; l < sz; l++) {
                      ijkl = (ntheta * nphi * npsi) * l + (nphi * npsi) * k + npsi * j + i;
                      ijkl_minus_one = (ntheta * nphi * npsi) * (l - 1) + (nphi * npsi) * k + npsi * j + i;
                      gijkl = (Gtijk[ijkl] + Gtijk[ijkl_minus_one]) * h / 2.0;
                      Pijk[ijk] += gijkl;
                      Mijk[ijk] += (T[l] - h2) * gijkl;
                 }
                 Mijk[ijk] /= (Pijk[ijk] + eps);
            }
            phi += dw;
        }
        psi += dw; 
    }
} /* ved4cirle300 */

void pool_orthants(int lambda, double *Gtb, double *Pthetab, double *Mtb, double *Gtijk, double *Pijk, double *Mijk) {
    /* Pool mass within biased orthants - first half <= pi/2, second half greater than pi/2 */

    int orthx[16] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    int orth1, orth2, orth3, orth4, o, i,j,k,l, ijk, ijkl;

    /* Initialize space for RT distributions and choice probabilities pooled within orthants */
    for (o = 0; o < north; o++) {
        Pthetab[o] = 0;
        Mtb[o] = 0;
        for (l = 0; l < sz; l++) {
              Gtb[north * l + o] = 0;
        }
    }
    for (i = 0; i < npsi; i++) {
        if (i < npsi / 2 - lambda) {
             orth1 = 0;
        } else {
             orth1 = 1;
        }
        for (j = 0; j < nphi; j++) {
            if (j < nphi / 2 - lambda) {
                 orth2 = 0;
            } else {
                 orth2 = 1;
            }
            /* This one on 0-2pi */
           for (k = 0; k < ntheta; k++) {
                if (k >= lambda && k < ntheta / 4 - lambda) {  
                    orth3 = 0;  /* (+,+) */
                    orth4 = 0;
                } else if (k < ntheta / 2 - lambda) {  
                    orth3 = 0;  /* (-,+) */
                    orth4 = 1;
                } else if (k < 3 * ntheta / 4 + lambda) {  
                    orth3 = 1;  /* (-,-) */
                    orth4 = 0;
                } else { 
                    orth3 = 1; /* (+,-) */
                    orth4 = 1;
                }
                o = 8 * orth1 + 4 * orth2 + 2 * orth3 + orth4;
                if (orthx[o] == 0) {
                   orthx[o] = 1;
                }
                ijk = (nphi * npsi) * k + npsi * j + i;
                Pthetab[o] += Pijk[ijk];
                Mtb[o] += Mijk[ijk] * Pijk[ijk];
                for (l = 0; l < sz; l++) { 
                    ijkl = (ntheta * nphi * npsi) * l + (nphi * npsi) * k + npsi * j + i;
                    Gtb[north * l + o] += Gtijk[ijkl];
                }
            }
        }
    }
} /* pool_orthants */
 
void vbd4orth300(double *T, double *Gt, double *Ft, double *Ptheta, double *Mt, double *P, double tmax, int badix) {
/* -----------------------------------------------------------------------------------------------
    Pool the choice probabilities and the RT distributions over the 16 orthants of S^3     
------------------------------------------------------------------------------------------------ */
    int l, o, lo, hi, big_index, big_index_minus;
    double c, lambda, frac, cfrac;
    /* Diffusion on S^3 */

    h = tmax / sz; 
    ved4cirle300(T, Gtijk, Pijk, Mijk, P, h, badix);
  
    /* Initialize space for cumulatives only - others done in pool_orthants */
    for (o = 0; o < north; o++) {
        for (l = 0; l < sz; l++) {
              Gt[north * l + o] = 0;
        }
    }

    /* Bias for nontarget orthants */
    c = P[BIAS];
    lambda = lambdac(c);
    lo = floor(lambda);
    hi = ceil(lambda);
    frac = lambda - lo;
    cfrac = 1.0 - frac;

/*
#ifdef SINGLEORTHANT
    mexPrintf("Single orthant \n");
#else
    mexPrintf("Double orthant \n");
#endif
*/
    pool_orthants(lo, Gtlo, Pthetalo, Mtlo, Gtijk, Pijk, Mijk);
    pool_orthants(hi, Gthi, Pthetahi, Mthi, Gtijk, Pijk, Mijk);
    /* mexPrintf("lambda = %6.3f, lo = %6d, hi = %6d, frac = %6.3f \n", lambda, lo, hi, frac); */
    for (o = 0; o < north; o++) {
        Ptheta[o] = cfrac * Pthetalo[o] + frac * Pthetahi[o]; 
        /* Ptheta[o] = Pthetalo[o]; */

        if (Ptheta[o] > 0) {
             Mt[o] = (frac * Mtlo[o] + cfrac * Mthi[o]) / Ptheta[o];
         } else {
             Mt[o] = 0;
        }
        for (l = 1; l < sz; l++) {
            big_index = north * l + o;
            big_index_minus = north * (l - 1) + o;
            /* Interpolate the density */
            Gt[big_index] = cfrac * Gtlo[big_index] + frac * Gthi[big_index]; 
            /* Gt[big_index] = Gtlo[big_index]; */

            /* Cumulative via Simpson's rule */
            Ft[north * l + o] += Ft[big_index_minus] + (Gt[big_index] + Gt[big_index_minus]) * h / 2.0;
        }
        /* mexPrintf("orthx[o] = %6d, Ptheta[o] = %8.5f Mt[o]= %8.5f\n", orthx[o], Ptheta[o], Mt[o]); */
    }         
} /* vbd4orth300 */

 
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
 /*
     =======================================================================
     Matlab gateway routine.
     =======================================================================
 */
 
    int badix; 
    double *T, *Gt, *Ft, *Ptheta, *Mt, *P;  /* Phase angle indexes point on 3-sphere */
    double tmax, badi;
    unsigned n, m;

    if (nrhs != 3) {
         mexErrMsgTxt("vbd4orth300: Requires 3 input args.");
    } else if (nlhs != 5) {
        mexErrMsgTxt("vbd4orth300: Requires 5 output args."); 
    }

    /*
      -----------------------------------------------------------------------
      Check all input argument dimensions.
      -----------------------------------------------------------------------   
    */

    /* P */
    m = mxGetM(prhs[0]);
    n = mxGetN(prhs[0]);
    if (!mxIsDouble(prhs[0]) || !(m * n == NP)) {
        mexPrintf("P is %4d x %4d \n", m, n);
        mexErrMsgTxt("vbd4orth3300: Wrong size P");
    }
    P = mxGetPr(prhs[0]);
    /* tmax */
    m = mxGetM(prhs[1]);
    n = mxGetN(prhs[1]);
    if (!mxIsDouble(prhs[1]) || !(m * n == 1)) {
        mexErrMsgTxt("vbd4orth3300: tmax must be a scalar");
    }
    tmax = mxGetScalar(prhs[1]);
    if (tmax <= 0.0) {
        mexPrintf("tmax =  %6.2f \n", tmax);
        mexErrMsgTxt("tmax must be positive");
    } 

    /* badi */
    m = mxGetM(prhs[2]);
    n = mxGetN(prhs[2]);
    if (!mxIsDouble(prhs[2]) || !(m * n == 1)) {
        mexErrMsgTxt("vbd4orth3300: badi must be a scalar");
    }
    badi = mxGetScalar(prhs[2]);
    badix = (int)(badi+0.5); 
 
    /*
      -----------------------------------------------------------------------
      Create output arrays.
      -----------------------------------------------------------------------
    */
 
    /* T */
    plhs[0] = mxCreateDoubleMatrix(1, sz, mxREAL);
    T = mxGetPr(plhs[0]);
    
    /* Gt */
    plhs[1] = mxCreateDoubleMatrix(north, sz, mxREAL);
    Gt = mxGetPr(plhs[1]);
    
   /* Ft */
    plhs[2] = mxCreateDoubleMatrix(north, sz, mxREAL);
    Ft = mxGetPr(plhs[2]);
    
    /* Ptheta */
    plhs[3] = mxCreateDoubleMatrix(north, 1, mxREAL);
    Ptheta = mxGetPr(plhs[3]);

   /* Mt */
    plhs[4] = mxCreateDoubleMatrix(north, 1, mxREAL);
    Mt = mxGetPr(plhs[4]);

    /*
      -----------------------------------------------------------------------
      Run the C-function.
      -----------------------------------------------------------------------
    */
    vbd4orth300(T, Gt, Ft, Ptheta, Mt, P, tmax, badix); 
}


