restart
needsPackage "MonodromySolver"
needs "common.m2"
declareVariable \ {x}
declareVariable \ {p}
X = gateMatrix{{x}}
P = gateMatrix{{p0,p1}}
F=gateSystem(P,X,gateMatrix{{p0*x^3+p1}})
(p0,x0)=createSeedPair F
evaluate(F,x0,p0)
V= first monodromySolve(F,point{{p0}},{point{{x0}}},)
points V.PartialSols


restart
setRandomSeed 2020
m = 3; -- number of cameras
D=(1,5,{{0,1},{0,2},{0,3},{4,5}})----degree=64
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

-- needed for what comes next
evaluate (GateMatrix, GateMatrix, Point) := (GM, varMat, xVals) -> (
    Peval := gateSystem(varMat, transpose gateMatrix{flatten entries GM});
    result := evaluate(Peval, xVals);
    matrix(result, numrows GM, numcols GM)
    )
p = V.BasePoint
xs = points V.PartialSols

-- all of the trifocal tensors are DISTINCT
sort apply(xs, x->(
(A,B,C) := toSequence pMatrices/(r->transpose evaluate(r, vars GS, x));
T := toList apply((0,0,0)..(2,2,2), (i,j,k) -> (-1)^(i+1) * det (submatrix'(A,{i}) | submatrix(B, {j}) | submatrix(C, {k})));
(1/T#0) * drop(T,1)
)
)
netList oo

-- essential matrices in second and third view give a hint
netList sort apply(xs, x->(
(E1,E2,E3) := toSequence pMatrices/(r->(
        c := evaluate(r, vars GS, x);
        essential(c_{0,1,2},c_{3})
        )
    );
(1/E2_(0,0))*drop(flatten entries E2,1) | (1/E3_(0,0))*drop(flatten entries E3,1)
)
)
netList oo

-- four blocks of size 16?
netList sort points V.PartialSols

-- essential matrices in second and third view give a hint
netList sort apply(xs, x->(
ts := pMatrices/(r->(
        c := evaluate(r, vars GS, x);
        flatten entries c_{3}
        )
    );
ts' := ts#1 | ts#2;
1/ts'#5 * drop(reverse ts',1)
)
)
concatenate ts


-- post-process the output of monodromySolver


writePermutations(V, "example" | ".txt")
