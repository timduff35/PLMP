-- run from PLMP/CODE/
-- neds M2 >= 1.14 (preferably 1.15)
restart
setRandomSeed 0
m=3
D = (4,0,{{0,1},{0,2},{1,2}}) -- 216, CLEVELAND, monodromy works
COFACTOR = true
JACOBIAN = true
RERUNMONODROMY = false
needsPackage "NumericalAlgebraicGeometry"
setDefault(tStepMin=>1e-7)
setDefault(maxCorrSteps=>2)
needs "numerical-problem-builder.m2"
(y, c) = fabricateyc CC
elapsedTime (V,np)=monodromySolve(GS, 
    y, {c},Verbose=>true,
    FilterCondition=>filterRank,
    Randomizer=>gammify,
    NumberOfEdges=>NEDGES
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
value linesOf g;
norm evaluate(GS,y1,first c1s)

-- code for generating various evaluators 
PH = parametricSegmentHomotopy GS

-- HxHt
h=cCode(transpose(PH.GateHomotopy#"Hx"|PH.GateHomotopy#"Ht"),gateMatrix{cameraVars|{PH.GateHomotopy#"T"}|flatten entries PH#Parameters})
-- HxH
h=cCode(transpose(PH.GateHomotopy#"Hx"|PH.GateHomotopy#"H"),gateMatrix{cameraVars|{PH.GateHomotopy#"T"}|flatten entries PH#Parameters})

