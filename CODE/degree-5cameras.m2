m = 5; -- number of cameras
DIR = "./numerical/degree-5cameras/";
COFACTOR = true;
JACOBIAN = true;
RERUNMONODROMY = true;

-- we should also stash working pivots to avoid recomputing Jacobian....
LIST = {    
--(3,1,{{1,2},{2,3}}) => "2011_1", -- deg = 11306 (not divisible by 8) -- 11296
(4,1,{{0,1,2,3},{0,4}}) => "2003_3",-- deg = 11008
(4,0,{{0,1,2},{0,3}}) => "2003_2" -- deg = 26240
    }

for P in LIST do (
    D = first P;
    FILENAME = DIR|last P|".txt";
    if (not fileExists FILENAME) then (
	<< D << endl;
	if (D == (4,1,{{0,1,2,3},{0,4}})) then Jpivots = {0,1,4,5,20,21,40,41,44,45,60,61,80,81,84,100,120,124,140,5100,5110,5120,5200,5215,5216,5217,5218,5219};
    	compTime = first elapsedTiming load("numerical-problem-builder.m2");
	stdio << "computation time = " << setupTime << " sec" << endl;
	stdio << "------------------------------------------" << endl;
	);
    )

end--
restart
needsPackage "NumericalAlgebraicGeometry"
setDefault(tStepMin=>1e-7)
setDefault(maxCorrSteps=>2)
needs "degree-5cameras.m2"
