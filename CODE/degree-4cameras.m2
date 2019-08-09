m = 4; -- number of cameras
DIR = "./numerical/degree-4cameras/";
COFACTOR = true;
JACOBIAN = true;
RERUNMONODROMY = true;

LIST = {    
    (2,3,{{1,2},{1,3},{1,4}}) => "2110_0", -- 32?
    (3,2,{{3,4}}) => "1030_0", ---- degree = 3040??
    (4,0,{{2,3}}) => "1022_2", -- degree = 4524? (not divisible by 8)
--    (5,0,{{1,2,3,4}}) => "1014_4", -- hedgehog + free line, not minimal
--    (6,0,{{0,1,2,3,4,5}}) => "1006_6", -- hedgehog, not minimal
    (4,0,{{0,1,3},{1,2},{0,2}}) => "3001_1", -- degree 1728??    
    (3,1,{{0,1},{0,2},{0,3}}) => "2102_1", -- 544?
    (3,2,{{0,1,2},{0,3},{0,4}}) => "2102_2" -- 544?
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
needs "degree-4cameras.m2"
