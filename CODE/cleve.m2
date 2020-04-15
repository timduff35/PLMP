-- run from PLMP/CODE/
-- neds M2 >= 1.14 (preferably 1.15)
restart
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
(V,np)=monodromySolve(GS, 
    y, {c},Verbose=>true,
    FilterCondition=>filterRank,
    Randomizer=>gammify,
    NumberOfEdges=>NEDGES
    );

-- initial parameters
p0 = matrix V.BasePoint
-- start solutions
x0s = matrix \ (points V.PartialSols)

-- code for generating various evaluators 
PH = parametricSegmentHomotopy GS

-- HxHt
h=cCode(transpose(PH.GateHomotopy#"Hx"|PH.GateHomotopy#"Ht"),gateMatrix{cameraVars|{PH.GateHomotopy#"T"}|flatten entries PH#Parameters})
-- HxH
h=cCode(transpose(PH.GateHomotopy#"Hx"|PH.GateHomotopy#"H"),gateMatrix{cameraVars|{PH.GateHomotopy#"T"}|flatten entries PH#Parameters})
