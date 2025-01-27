# RBiCG With Application to Model Reduction and Convection-Diffusion PDE

For questions contact:

Kapil Ahuja  
kahuja@vt.edu
http://www.math.vt.edu/people/kahuja/ 

Eric de Sturler   
sturler@vt.edu  
http://www.math.vt.edu/people/sturler/
Parent folder:
- rbicg.m: Code for recycling BiCG.

- rcgs.m: Code for recycling CGS.

- rbicgstab.m: Code for recycling BiCGSTAB.

- applyPrecond.m: Applies M1 & M2 precond via mat-vecs as inv(M1)*A*inv(M2). i.e. central or split precond.
				   
- binormalize.m: Normalizes columns of C and C_tilde with c1*c1_tilde etc.
				  
- getGenEigenvecs.m: Solves the generalized eigenvalue problem to build the recycle space.
					  
- orthogonalize.m: Builds bi-orthogonality between the C's.


Model_red folder:
- rail_models: These contain the matrices for four model sizes 1357, 5177, 20209, and 79841

- mmread.m: Used to read-in the above models.

- irkaRailPlots.m: IRKA that uses rbicg, and plots the convergence curves.

- irkaRailTime.m: IRKA that uses rbicg, and provides the time comparison.

Convection folder: 
- All *.m files except testConvDiff.m: Sets-up the convection-diffusion problem.

- testConvDiff.m: Call rbicg to solve the convection-diffusion problem.

-------------------------------------------------------------------------------------------------

To run model reduction code (assuming you are in parent folder):
>> cd model_red
>> irkaRailPlots(1)			% 1 is for size 1357 , 2 is for 5177, 3 for 20209, and 4 for 79841
or
>> irkaRailTime(1,0)		% First argument: 1 for 20209 and 2 for 79841
							% Second argument: 0 means no recycling, 1 means use recycling


To run convection-diffusion problem (assuming you are in parent folder):
>> cd convection
>> testConvDiff(1)			% 1 means use preconditiong, 0 means unpreconditioned solve

## *Citation*
If you use this code in your research, please cite:

Ahuja, Kapil, et al. "Recycling BiCG with an application to model reduction." SIAM Journal on Scientific Computing 34.4 (2012): A1925-A1949.

[https://doi.org/10.1137/100801500](https://doi.org/10.1137/100801500)
