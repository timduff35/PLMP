m = 4 -- number of cameras
needs "tester.m2"
DIR = "./dim-degree-F4-4cameras/"
LIST = {    
    (3,2,{{3,4}}) => "1030_0", ---- degree = 3067??
    (4,0,{{2,3}}) => "1022_2", -- degree = 4525?
    (5,0,{{1,2,3,4}}) => "1014_4", -- hedgehog + free line, not minimal
    (6,0,{{0,1,2,3,4,5}}) => "1006_6", -- hedgehog, not minimal
    (4,0,{{0,1,3},{1,2},{0,2}}) => "3001_1", -- degree 1728??    
    (2,3,{{1,2},{1,3},{1,4}}) => "2110_0",
    (3,1,{{0,1},{0,2},{0,3}}) => "2102_1",
    (3,2,{{0,1,2},{0,3},{0,4}}) => "2102_2"
    }
FF = ZZ/nextPrime 10000
for p in LIST do (
    isParametric = true;
    minimal := isMinimal first p;
    << "isMinimal = " << minimal << endl; 
    if member(last p, {}) then(
	isParametric = false;
	elapsedTime computeDimDegreeF4(first p, FileName=>DIR|last p|".txt");
	);
    print "-------------------------------------------------------------";
    )
end --------------------------------------------------------

for p in LIST do 
computeDimDegreeF4(first p, FileName=>DIR|last p|".txt") 
D=(2,3,{{1,2},{1,3},{1,4}})
isParametric=false
computeDimDegreeF4(D, FileName=>DIR|"2110_0"|".txt") 

end --------------------------------------------------------

restart
load "dim-degree-F4-4cameras.m2"
