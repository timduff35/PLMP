-- EXAMPLE of numerical degree computation for problem 2111_1 (minimal, degree=40) 
restart
setRandomSeed 0 -- ensures reproducibility
m = 3 -- number of cameras
D = (3,2,{{1,2},{1,3},{1,4}}) -- encodes problem 2111_1 
-* line 0 is free
   line 1 passes through points 0,1,2
   line 2 passes through point 0
   line 3 passes through point 1 (ghost line)
   line 4 passes through point 2 (ghost line)
*-    

-*
global variables for configuring monodromy run:
  - Jpivots indexes a square subsystem of the minors equations (hardcoded for this example)
  - RERUNMONODROMY=true tells the builder script to re-run monodromy
*-
Jpivots = {0, 1, 4, 5, 8, 9, 14, 29, 33, 44, 48, 57, 58, 59}
RERUNMONODROMY = true;
needs "numerical-problem-builder.m2"
y = matrix V.BasePoint; -- parameters for specialized system solved by monodromy
sols = for x in points V.PartialSols list transpose matrix x; -- solutions of specialized system
residuals = for x in sols list norm evaluate(F,x||p)
max residuals
