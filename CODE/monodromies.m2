DIR = "./numerical/monodromies" |toString(m) |"/";
COFACTOR = true;
JACOBIAN = true;
RERUNMONODROMY = true;

LIST = if (m==2) then
{    
(5,0,{{0,1},{1,2},{2,3},{3,4},{4,0}}) => "5000",
(4,1,{{0,1},{1,2},{2,3},{3,0},{0,4}}) => "4100_3",
(2,4,{{0,1},{0,2},{0,3},{1,4},{1,5}}) => "3200_3",

} else if (m==3) then 
{
    (6,1,{{0,6},{0,1,2,3,4,5}}) => "2005_5",--64
    (4,2,{{4,5}}) => "1040_0",--360
    (5,0,{{3,4}}) => "1032_2",
    (6,0,{{2,3,4,5}}) => "1024_4",
    (3,2,{{1,2},{1,3},{1,4}})=>"2111_1",
    (4,0,{{0,1},{0,2},{0,3}})=>"2103_1",
    (4,1,{{0,1,2},{0,3},{0,4}})=>"2103_2",
    (4,2,{{0,1,2,3},{0,4},{0,5}})=>"2103_3",
    (1,5,{{0,1},{0,2},{0,3},{4,5}})=>"3100_0",
    (2,3,{{0,1},{1,2},{1,3},{1,4}})=>"2201_1",
    (4,0,{{1,2},{2,3},{1,3}})=> "3010",
    (5,0,{{0,1,3},{0,2,4},{1,2}}) => "3002_1",
    (6,0,{{0,1,2,3,4},{0,5}}) => "2005_4",
    (4,1,{{2,3},{2,4}}) => "2021_1"
    } else if (m==4) then
{    
    -- degree 4
    (3,1,{{0,1},{0,2},{0,3}}) => "2102_1", -- 544?
    (3,2,{{0,1,2},{0,3},{0,4}}) => "2102_2", -- 544?
    (2,3,{{1,2},{1,3},{1,4}}) => "2110_0" -- 32?
    } else error "invalid m"

for P in LIST do (
    D = first P;
    FILENAME = DIR|last P|".txt";
    if (not fileExists FILENAME) then (
	<< D << endl;
    	compTime = first elapsedTiming load("numerical-problem-builder.m2");
	stdio << "computation time = " << setupTime << " sec" << endl;
	stdio << "------------------------------------------" << endl;
        writePermutations(V, last P | ".txt")
	);
        
	);
    
    )

end--
restart
setRandomSeed 0
needsPackage "NumericalAlgebraicGeometry"
setDefault(tStepMin=>1e-8)
setDefault(maxCorrSteps=>2)
m=3
needs "monodromies.m2"


-*
note: random seed 999999 for problem 2102_2 results in contamination without the filter, tstepMin=1e-16, maxCorr=>3
*-
