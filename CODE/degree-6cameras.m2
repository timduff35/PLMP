m = 6; -- number of cameras
DIR = "./numerical/degree-6cameras/";
COFACTOR = true;
JACOBIAN = true;
RERUNMONODROMY = true;

-- we should also stash working pivots to avoid recomputing Jacobian....
LIST = {    
    (3,1,{{2,3}}) => "1021_1", -- degree > 450k
    }

for P in LIST do (
    D = first P;
    FILENAME = DIR|last P|".txt";
    if (not fileExists FILENAME) then (
	<< D << endl;
    	compTime = first elapsedTiming load("numerical-problem-builder.m2");
	stdio << "computation time = " << setupTime << " sec" << endl;
	stdio << "------------------------------------------" << endl;
	);
    )

end--
restart
needs "degree-6cameras.m2"
