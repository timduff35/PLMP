DIR = "./monodromies" |toString(m) |"/";

LIST = if (m==2) then
{    
(5,0,{{0,1},{1,2},{2,3},{3,4},{4,0}}) => ("5000", 20),
(4,1,{{0,1},{1,2},{2,3},{3,0},{0,4}}) => ("4100_3", 16),
(2,4,{{0,1},{0,2},{0,3},{1,4},{1,5}}) => ("3200_3", 12),

} else if (m==3) then 
--todo: annotate w/ degrees
--writePermutations should _fail_ if the base point doesn't have expected degree many solutions
{
    (4,2,{{4,5}}) => ("1040_0", 360),
    (5,0,{{3,4}}) => ("1032_2", 552),
    (6,0,{{2,3,4,5}}) => ("1024_4", 480),
    (4,1,{{2,3},{2,4}}) => ("2021_1", 264),
    (5,0,{{0,1,2},{0,3}}) => ("2013_2", 432),
    (5,1,{{0,1,2,3},{0,5}}) => ("2013_3", 328),        
    (3,2,{{1,2},{1,3},{1,4}})=> ("2111_1", ),
    (4,0,{{0,1},{0,2},{0,3}})=> ("2103_1",),
    (4,1,{{0,1,2},{0,3},{0,4}})=>"2103_2",
    (4,2,{{0,1,2,3},{0,4},{0,5}})=>"2103_3",
    (1,5,{{0,1},{0,2},{0,3},{4,5}})=>"3100_0",
    (2,3,{{0,1},{1,2},{1,3},{1,4}})=>"2201_1",
    (4,0,{{1,2},{2,3},{1,3}})=> "3010",
    (5,0,{{0,1,3},{0,2,4},{1,2}}) => "3002_1",
    (6,0,{{0,1,2,3,4},{0,5}}) => "2005_4",
    (6,1,{{0,6},{0,1,2,3,4,5}}) => ("2005_5", 64)
    } else if (m==4) then
{    
    -- degree 4
--    (3,1,{{0,1},{0,2},{0,3}}) => "2102_1", -- 544?
--    (3,2,{{0,1,2},{0,3},{0,4}}) => "2102_2", -- 544?
    (2,3,{{1,2},{1,3},{1,4}}) => "2110_0" -- 32?
    } else error "invalid m"

for P in LIST do (
    D = first P;
    (file, rc) := (first last P, last last P);
    FILENAME = DIR|last P|".txt";
    if (not fileExists FILENAME) then (
	<< D << endl;
        if (rc == 32) then NNODES = 4 else NNODES = 3;
    	compTime = first elapsedTiming load("numerical-problem-builder.m2");
	stdio << "computation time = " << setupTime << " sec" << endl;
	stdio << "------------------------------------------" << endl;
        writePermutations(V, rc, FILENAME)
	);
        
	);
    
    )

setRandomSeed 0
COFACTOR = true; -- strategy for how to evaluate determinants
RUNMONODROMY = true;
NEDGES=4;
NNODES=3
SATURATE=true;
D=(6,1,{{0,6},{0,1,2,3,4,5}})
m=4
needsPackage "NumericalAlgebraicGeometry"
setDefault(tStepMin=>1e-8)
setDefault(maxCorrSteps=>2)
FILENAME="test.txt"
D=(2,3,{{1,2},{1,3},{1,4}})
compTime = first elapsedTiming load("numerical-problem-builder.m2");
writePermutations(V, 32, FILENAME)

needs "monodromies.m2"


-*
note: random seed 999999 for problem 2102_2 results in contamination without the filter, tstepMin=1e-16, maxCorr=>3
*-
