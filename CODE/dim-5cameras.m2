m = 5 -- number of cameras
needs "tester.m2"
DIR = "./dim-degree-F4-5cameras/"
LIST = {    
    (3,1,{{1,2},{2,3}}) => "2011_1", -- deg = 11306 ?
    (4,0,{{0,1,2},{0,3}}) => "2003_2", -- deg = 26240
    (4,1,{{0,1,2,3},{0,4}}) => "2003_3" -- deg = ??
    }
FF = ZZ/nextPrime 10000
isParametric = true
for p in LIST do (
    minimal := isMinimal first p;
    << "isMinimal = " << minimal << endl; 
    print "-------------------------------------------------------------";
    )
end --------------------------------------------------------

for p in LIST do 
computeDimDegreeF4(first p, FileName=>DIR|last p|".txt") 

end --------------------------------------------------------

restart
load "dim-5cameras.m2"
