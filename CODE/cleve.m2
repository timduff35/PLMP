-- run from PLMP/CODE/
-- needs M2 >= 1.14 (preferably 1.15)
restart
setRandomSeed 0
m=3 -- number of cameras 
D = (4,0,{{0,1},{0,2},{1,2}}) -- encoding for CLEVELAND, a problem of degree 216

-- global variables configuring "numerical-problem-builder.m2"
COFACTOR = true -- "true" means 2x2 and 3x3 determinants are evaluated by Laplace expansion
JACOBIAN = true -- 
RERUNMONODROMY = false -- we will run it ourselves a few lines below

needsPackage "NumericalAlgebraicGeometry"
--  path-tracker settings that are "safer" than the default
setDefault(tStepMin=>1e-7)--default 1e-6
setDefault(maxCorrSteps=>2)--default 3

-- running the problem builder takes a few seconds
needs "numerical-problem-builder.m2"
-- fabricate a starting pair (y,c) for monodromy
-- c: random camera parameters
-- y: image data consistent with c (plus charts on projective quantities)
(y, c) = fabricateyc CC
--running monodromy
errorDepth = 0
elapsedTime (V,np)=monodromySolve(GS, 
    y, {c},Verbose=>true,
    FilterCondition=>filterRank,
    Randomizer=>gammify,
    NumberOfNodes=>NNODES,
    NumberOfEdges=>NEDGES,
    );

-- initial parameters
y0 = matrix V.BasePoint
-- start solutions
c0s = matrix \ (points V.PartialSols)

-- write to file
f = openOut "cleve.txt"
f << "y1=" << toExternalString(y0) << endl
f << "c1s=" << toExternalString(c0s) << endl
close f

-- read from file
g=openIn "cleve.txt"
linesOfg = get g;
value linesOfg;
norm evaluate(GS,y1,first c1s)

-- write to file w/0 macaulay2 format
f = openOut "cleve_raw.txt"
f << "parameters" << endl
f << toString(flatten entries y0) << endl
f << "solutions" << endl
for c0 in c0s do f << toString(flatten entries c0) << endl
close f

-- code for generating various evaluators 
PH = parametricSegmentHomotopy GS

-- HxHt
h=cCode(transpose(PH.GateHomotopy#"Hx"|PH.GateHomotopy#"Ht"),gateMatrix{cameraVars|{PH.GateHomotopy#"T"}|flatten entries PH#Parameters})
-- HxH
h=cCode(transpose(PH.GateHomotopy#"Hx"|PH.GateHomotopy#"H"),gateMatrix{cameraVars|{PH.GateHomotopy#"T"}|flatten entries PH#Parameters})
--todo: pipe cCode output to text file
