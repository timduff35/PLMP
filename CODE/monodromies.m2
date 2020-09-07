DIR = "monodromies/monodromies" |toString(m) |"/";

LIST = if (m==2) then
{    
(5,0,{{0,1},{1,2},{2,3},{3,4},{4,0}}) => ("5000", 20), -- IMPRIMITIVE: C2 wr S10 ^ A20
(4,1,{{0,1},{1,2},{2,3},{3,0},{0,4}}) => ("4100_3", 16), -- IMPRIMITIVE: C2 wr S8 ^ A16
(2,4,{{0,1},{0,2},{0,3},{1,4},{1,5}}) => ("3200_3", 12) -- IMPRIMITIVE
-*
result of StructureDescription: (C2 x C2) : S4 (":" indicates semidirect product)
-- !!! S4 acts on pairs (i,j) w/ 1<=i<j<=4 !!!
-- ie S4 = C2 wr S3 ^ A6
--- better to write as C2 wr (C2 wr S3 ^ A6) ^ A12
*-
} else if (m==3) then 
--todo: annotate w/ degrees
--writePermutations should _fail_ if the base point doesn't have expected degree many solutions
{
    (4,2,{{4,5}}) => ("1040_0", 360), -- full symmetric
    (5,0,{{3,4}}) => ("1032_2", 552), -- full symmetric
    (6,0,{{2,3,4,5}}) => ("1024_4", 480), -- full symmetric
    (4,1,{{2,3},{2,4}}) => ("2021_1", 264), -- full symmetric
    (5,0,{{0,1,2},{0,3}}) => ("2013_2", 432), -- full symmetric
    (5,1,{{0,1,2,3},{0,5}}) => ("2013_3", 328),  -- full symmetric
    (6,0,{{0,1,3,4},{0,2,5}}) => ("2005_3", 480), -- full symmetric
    (6,0,{{0,1,2,3,4},{0,5}}) => ("2005_4", 240), -- full symmetric
    (6,1,{{0,6},{0,1,2,3,4,5}}) => ("2005_5", 64), -- IMPRIMITIVE
    -*
     Order = 2^27 * (2^4 * 4!)
     Centralizer in S64 = "C2 x C2 x C2"
     taking successive quotients by nontrivial action on blocks, there is always a nontrivial centralizer
     "Base Group" identified as C2 wr S4     
     "fiber groups" not obvious to identify
     gap> Order(G)/Order(G1);
     16
     gap> Order(G1)/Order(G2);
     128
     gap> Order(G2)/Order(G3);
     256
    so perhaps C2 wr (C2 wr (C2 wr (C2 wr S4)) ^ A32) ^ (Z2^4 -> S64)
    this one should probably be rerun with more loops!
    *-
    (4,0,{{1,2},{2,3},{1,3}})=> ("3010_0", 216), -- full symmetric
    (5,0,{{0,1,3},{0,2,4},{1,2}}) => ("3002_1", 312), -- GAP didnt finish quickly
    (5,0,{{0,1,3,4},{1,2},{0,2}}) => ("3002_2", 224), -- full symmetric
    (3,2,{{1,2},{1,3},{1,4}})=> ("2111_1", 40), -- full symmetric
    (4,0,{{0,1},{0,2},{0,3}})=> ("2103_1", 144), -- full symmetric
    (4,1,{{0,1,2},{0,3},{0,4}})=>("2103_2", 144), -- full symmetric
    (4,2,{{0,1,2,3},{0,4},{0,5}})=>("2103_3", 144), -- GAP didnt finish quickly
    (1,5,{{0,1},{0,2},{0,3},{4,5}})=>("3100_0", 64) -- IMPRIMITIVE: C2 wr (C2 wr S16 ^ A32) ^ A64
    } else if (m==4) then
{    
    -- degree 4
    (3,1,{{0,1},{0,2},{0,3}}) => ("2102_1", 544), -- full symmetric
    (3,2,{{0,1,2},{0,3},{0,4}}) => ("2102_2", 544), -- full symmetric
    (2,3,{{1,2},{1,3},{1,4}}) => ("2110_0", 32) -- IMPRIMITIVE: -*
    -*
    structure description | (C2 x C2) : ((C2 x C2) : ((C2 x C2 x C2 x C2) : C2))
    using 46 loops didn't increase group size, so previous 1024 may have come from path jumping?
    *-
    } else error "invalid m"

for P in LIST do (
    D = first P;
    (file, rc) := (first last P, last last P);
    FILENAME = DIR|first last P|".txt";
    if (not fileExists FILENAME) then (
	<< D << endl;
        if (rc == 32) then NNODES = 4 else NNODES = 3;
        ROOTCOUNT = rc;
    	compTime = first elapsedTiming load("numerical-problem-builder.m2");
	stdio << "computation time = " << compTime << " sec" << endl;
	stdio << "------------------------------------------" << endl;
        elapsedTime writePermutations(V, rc, FILENAME)
	);
        
	);
    
end
restart
m=4
NNODES=3
NEDGES=3
RUNMONODROMY=true
SATURATE=true
COFACTOR=true
needsPackage "NumericalAlgebraicGeometry"
setDefault(tStepMin=>1e-8)
setDefault(maxCorrSteps=>2)
needs "monodromies.m2"

end
restart
setRandomSeed 0
COFACTOR = true; -- strategy for how to evaluate determinants
RUNMONODROMY = true;
NEDGES=5;
NNODES=5
--betti for complete graph with multiple edges
betti (ZZ, ZZ) := o -> (n, m) -> m*(n^2-n)/2-n+1
betti(NNODES, NEDGES)
SATURATE=true;
D=(6,1,{{0,6},{0,1,2,3,4,5}})
m=4
needsPackage "NumericalAlgebraicGeometry"
setDefault(tStepMin=>1e-8)
setDefault(maxCorrSteps=>2)
FILENAME="test.txt"
D=(2,3,{{1,2},{1,3},{1,4}})
compTime = first elapsedTiming load("numerical-problem-builder.m2"); -- time to initialize graph can be wayy shorter here
elapsedTime writePermutations(V, 32, FILENAME)

needs "monodromies.m2"


-*
note: random seed 999999 for problem 2102_2 results in contamination without the filter, tstepMin=1e-16, maxCorr=>3
*-



restart
setRandomSeed 0
COFACTOR = true; -- strategy for how to evaluate determinants
RUNMONODROMY = true;
NEDGES=5;
NNODES=5
--betti for complete graph with multiple edges
betti (ZZ, ZZ) := o -> (n, m) -> m*(n^2-n)/2-n+1
betti(NNODES, NEDGES)
SATURATE=true;

m=3
needsPackage "NumericalAlgebraicGeometry"
setDefault(tStepMin=>1e-8)
setDefault(maxCorrSteps=>2)
FILENAME="test.txt"
D=(6,1,{{0,6},{0,1,2,3,4,5}})
compTime = first elapsedTiming load("numerical-problem-builder.m2"); -- time to initialize graph can be wayy shorter here
elapsedTime writePermutations(V, 64, FILENAME)

 => ("2005_5", 64), -- IMPRIMITIVE
