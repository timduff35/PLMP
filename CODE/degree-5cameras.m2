m = 5; -- number of cameras
DIR = "./numerical/degree-5cameras/";
COFACTOR = true;
JACOBIAN = true;
RERUNMONODROMY = true;

-- we should also stash working pivots to avoid recomputing Jacobian....
LIST = {    
(4,1,{{0,1,2,3},{0,4}}) => "2003_3",-- deg = ??
(3,1,{{1,2},{2,3}}) => "2011_1", -- deg = 11306 ?
(4,0,{{0,1,2},{0,3}}) => "2003_2" -- deg = 26240
    }

for P in LIST do (
    D = first P;
    FILENAME = DIR|last P|".txt";
    if (not fileExists FILENAME) then (
	<< D << endl;
    	compTime = first elapsedTiming load("degree.m2");
	stdio << "computation time = " << setupTime << " sec" << endl;
	stdio << "------------------------------------------" << endl;
	);
    )

end--
restart
needs "degree-5cameras.m2"
