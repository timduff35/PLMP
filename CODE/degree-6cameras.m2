needsPackage "MonodromySolver"
setDefault(tStepMin=>1e-7)
setDefault(maxCorrSteps=>2)
m = 6; -- number of cameras
DIR = "./numerical/degree-6cameras/";
COFACTOR = true;
JACOBIAN = true;
RERUNMONODROMY = true;

-- we should also stash working pivots to avoid recomputing Jacobian....
LIST = {    
    (3,1,{{2,3}}) => "1021_1" -- degree > 450k
    }

for P in LIST do (
    D = first P;
    FILENAME = DIR|last P|".txt";
    if (not fileExists FILENAME) then (
	<< D << endl;
	Jpivots = {0, 1, 4, 5, 16, 17, 40, 41, 80, 81, 84, 85, 96, 97, 120, 121, 160, 161, 164, 165, 176, 177, 200, 201, 295, 330, 386, 470, 590, 735, 736, 737, 738, 739, 740};
    	compTime = first elapsedTiming load("numerical-problem-builder.m2");
	stdio << "computation time = " << setupTime << " sec" << endl;
	stdio << "------------------------------------------" << endl;
	);
    )

end--
restart
needs "degree-6cameras.m2"
Jpivots
