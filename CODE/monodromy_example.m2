restart
setRandomSeed 2019
m = 3; -- number of cameras
D=(1,5,{{0,1},{0,2},{0,3},{4,5}})--D=(2,3,{{1,2},{1,3},{1,4}}) --degree=32??
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

-- four blocks of size 8?
netList sort points V.PartialSols

-- post-process the output of monodromySolver


writePermutations(V, "example" | ".txt")
