restart
setRandomSeed 2019
m = 4; -- number of cameras
D=(2,3,{{1,2},{1,3},{1,4}}) --degree=32??
COFACTOR = true; -- SLPs for determinants are computed using cofactor expansion
 -- instead of evaluate a 246 x 14 Jacobian to select a square subsystem, use hard-coded pivots
JACOBIAN = false;
Jpivots = {0, 1, 4, 5, 16, 17, 20, 21, 41, 51, 71, 111, 121, 141, 181, 191, 211, 242, 243, 244, 245}
RERUNMONODROMY = true; -- runs monodromy
SATURATE = true
needsPackage "NumericalAlgebraicGeometry"
setDefault(tStepMin=>1e-7)
setDefault(maxCorrSteps=>2)
needs "numerical-problem-builder.m2"

-- four blocks of size 8?
netList sort points V.PartialSols

G=V.Graph
V1=first G.Vertices
V2=last G.Vertices
E1=toList V1.Edges
-- "petal loops" based at V1
e1=first E1
perms = apply(drop(E1,1),e->values pCompose(e1.Correspondence12,e.Correspondence21))
writePermutations(perms,"32.tst")
