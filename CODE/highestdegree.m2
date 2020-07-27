m = 6; -- number of cameras
D = (3,1,{{2,3}})
COFACTOR = true; -- SLPs for determinants are computed using cofactor expansion
 -- instead of evaluate a 246 x 14 Jacobian to select a square subsystem, use hard-coded pivots
JACOBIAN = true;
--Jpivots = {0, 1, 4, 5, 16, 17, 20, 21, 41, 51, 71, 111, 121, 141, 181, 191, 211, 242, 243, 244, 245}
RERUNMONODROMY = true; -- runs monodromy
SATURATE = true
needsPackage "NumericalAlgebraicGeometry"
setDefault(tStepMin=>1e-7)
setDefault(maxCorrSteps=>2)
needs "numerical-problem-builder.m2"
end
restart
setRandomSeed 0
load "highestdegree.m2"--0 sols found, tim laptop, M2-master, MS 1.16, 7/27
