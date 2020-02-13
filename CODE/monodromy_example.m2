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

writePermutations = method(Options=>{})
writePermutations (HomotopyNode, String) := o -> (V, filename) -> (
    -- assume for now that G has 2 nodes
    G := V.Graph;
    V1 := first G.Vertices;
    V2 := last G.Vertices;
    E1 := toList V1.Edges;
    -- "petal loops" based at V1
    e1 := first E1;
    perms := apply(drop(E1,1),e->values pCompose(e1.Correspondence12,e.Correspondence21));
    f := openOut filename;
    writePermutations(perms,filename);
    close f;
    )

L = points V.PartialSols
coordinates L#0
coordinates L#43

FF = frac(QQ[y]/(y^2-3))
R=FF[b,x]
f=x^2+b*x+b^2
(x-b*
discriminant(f,x)
primaryDecomposition ideal(f)

S=first flattenRing R
discriminant(sub(f,S),x_S)
x0 = first roots sub(sub(sub(f, S),{b=>b0_S}),CC[gens R])
V = first monodromySolve(
    polySystem{f},
    point{{b0}},
    {point{{x0}}},
    NumberOfNodes=>4,
    Verbose=>true
    )
points V.PartialSols
writePermutations(V, 
discriminant(f,x)
