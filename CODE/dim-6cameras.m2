m = 6 -- number of cameras
needs "tester.m2"
DIR = "./dim-degree-F4-6cameras/"
LIST = {    
    (3,1,{{2,3}}) => "1021_1", -- degree > 450k
    (4,0,{{1,2,3}}) => "1013_3",
    (5,0,{{0,1,2,3,4}}) => "1005_5"
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
load "dim-6cameras.m2"
