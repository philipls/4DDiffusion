Matlab routines to fit the 4D hyperspherical diffusion model to double-target detection data
To accompany Smith, P. L. and Corbett, E. A. "Speeded Multielement Decision Making as Diffusion in a Hypersphere: Theory and Application to Double-Target Detection"


I. Resources (Matlab)
---------------------

1. d4d.m - Fitting routine for double-target task
2. d4s.m - Fitting routine for single-target task
3. duqpfd.m - Quantile-probability plot for double-target task
4. duqpfs.m - Quantile-probability plot for single-target task
5. DoubleFit.mat - Group and individual data for double-target task; example fit parameters
6. SingleFit.mat - Group and individual data for single-target task; example fit parameters

II. Dependencies
-----------------
d4d  --> ved4sphere --> vedorth300.c

d4s  --> vbs4sphere --> vbd4orth300x.c

The C-routines were compiled under Linux using the GNU scientific library, gsl, and
the gslcblas library. These routines provide Bessel functions and their roots. Other libraries may provide these functions in other environments. The C-routines are compiled under Matlab with the commands:

   mex ved4orth300.c -lgsl -lgslcblas  -lm

   mex vbd4orth300x.c -lgsl -lgslcblas  -lm 

III. Data Structures
--------------------

The data structure is a 10 x 14 array. The rows are stimulus conditions; the columns are data for a condition. The first 7 columns are target-present responses; the next 7 columns are target-absent responses. The columns are P(R) (response probability); Ni (number of trials) and Q1...Q5, the .1, .3, .5, .7, and .9 RT quantiles in seconds. The odd-numbered rows are target-present stimuli; the even-numbered rows are target-absent stimuli. 

In DoubleFit.mat, DuncanD is the group data for the double-target task and DuncanDi is a 6-element cell array whose components are the data for the individual subjects. In SingleFit.mat, DuncanS is the group data for the single-target task and DuncanSi is a 6-element cell array whose components are the data for the individual subjects.


IV. Calling Conventions
-----------------------
Help calls "help d4d", "help d4s," etc., give the calling conventions. For d4d the call is:


	[G2, B, Pred] = d4d(Pvar, Pfix, Sel, Dmatrix, trace);

The function takes a vector of variable parameters, Pvar, which are estimated during fitting, a vector of fixed parameters, Pfix, which remain fixed, a selector vector, Sel, used to select and assemble parameter vectors, a data matrix, Dmatrix, and an optional trace switch, which gives information about parameter bounds. How the elements of Pfix are treated depends on the internal logic of d4d and d4s. These parameters may either be held fixed to their starting values or forced to be equal to other parameters. This logic can be adapted to the user's needs.

DoubleFit.mat contains a typical example, which fits the double-target task. The fit is controlled by the 33-element selector vector, Sel. This has 0's in positions 7-10 and 12-15. The resulting model allows drift rates for 2-target signals and 1-target noise to vary freely with contrast, but forces the drifts for 2-target noise and 1-target signals to be the same at all levels of contrast. Only one 2-target noise and 1-target signal parameter is estimated (controlled by the 1's in Sel). The others are forced to be the same by tests internal to d4d. The call is

	[g2,b,Pred]=d4d(P(Sel==1), P(Sel==0), Sel, DuncanD)

Here P is a 33-element parameter vector that is assembled from its fixed and variable components as described below. The returned values are g2, some model selection statistics, and a 10 x 14 matrix of predictions, Pred. Pred has the same form as the data matrix, except it returns means in columns 2 and 9. The elements of b are:

	b = [trueG2, adjG2, BIC, QAIC, dfr]; 

The the QAIC and the overdispersion-adjusted G2 depend on an OD parameter q, which must be calculated externally. It is set inside d4d and d4s to q = 1.0 by default.

SingleFit.mat contains another example, which fits Model 8 to the single-target task. The values of Sel force all of target-present noise drift rates to zero and estimate only two eta parameters, one for target-present and one for target-absent trials. It is called in the same way, as

	[g2,b,Pred]=d4s(P(Sel==1), P(Sel==0), Sel, DuncanS)


V. Fitting Sequence
-------------------
The typical fitting sequence in Matlab is as follows. For d4d you need to begin with a 33-element vector of starting parameters, P, and a 33-element selector vector, Sel. For d4s, P and Sel must both be 29 elements. Then, for Nelder-Mead Simplex fitting:

	settopt % Creat an options structure and turn on iteration trace

	pest = fminsearch(@d4d, P(Sel==1), options, P(Sel==0), Sel, DuncanD) % Do the fit

	P(Sel==1) = pest % Update the working parameter vector

	[g2,b,Pred]=d4s(P(Sel==1), P(Sel==0), Sel, DuncanS)  % Generate predictions with the new parameters

	duqpfd(DuncanD, Pred) % Generate quantile-probability plot of the fit

VI. Miscellaneous Notes and Cautions
------------------------------------

There are a number of context-sensitive flags and variables that are set internally in d4d and d4s. These include q, the overdispersion (which needs be calculated externally from the session-by-session or block-by-block data), N, the number of trials per condition, (currently set at N=168) and tmax, the maximum time index (currently set to 4.0). There is also a combination of hard and soft (quadratically penalized) constraints on the parameters to keep Simplex out of parts of the space where things may go bad numerically. These are very ad hoc, but work well enough in practice, so long as you're aware of what the settings are. Setting the trace flag to 1 in d4d and d4s will give some information about where the current parameter vector is in relation to the constraint set. The number of time steps on the range [0, tmax] is currently set to 300 and the number of hitting point bins on the range [0, 2pi] is currently set to 36. Both d4d and d4s save a copy of the working parameter vector to disk each time they are called. These calls add some disk-access overheads but are convenient if you want to exit a fit prematurely or to diagnose where something went bad. Because of the context-sensitive nature of these flags, these routines should not be treated as black-boxes, but instead as starting points that users should adapt to their own requirements. 

These routines are supplied as is, for nonprofit research purposes only, without any assumed, intended, or implied liability on the part of the author.  










