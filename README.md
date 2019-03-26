# PLMP
This repository hosts code for checking minimality and computing algebraic degrees of relative pose estimation problems in computer vision involving complete correspondences of points, lines, and their incidence relations across some number of views.

To get started, navigate to the "CODE" directory and execute the lines of "example.m2" in Macaulay2.

## Files in CODE/

* Numerical degree computations: 
  + degree-5cameras.m2	
  + degree-4cameras.m2  
* Symbolic minimality check and degree computations:
  + dim-5cameras.m2	
  + dim-6cameras.m2		   
  + dim-degree-F4-3cameras.m2  
  + dim-degree-F4-4cameras.m2
  + dim-degree-F4-2cameras.m2  
* Methods and functions used by the scripts above:
  + numerical-problem-builder.m2  
  + partial-view-builder.m2	
  + problem-builder.m2	     
  + tester.m2
  + common.m2	   
  + partial-line-builder.m2       
  + problem-builder-core.m2	
  + problem-builder-matrices.m2  
  + vision-toolbox.m2

## Caveats
   
* Certain computations may take a long time (couple of hours) and may exceed the RAM capacity (one example employing F4 consumed ~16Gb).  
* At the moment the numerical degree computations need a version of Macaulay2 at [this branch of this fork](https://github.com/timduff35/M2/tree/monodromy). 
