m = 2 -- number of cameras
needs "tester.m2"
DIR = "./dim-degree-F4-2cameras/"
LIST = {    
(5,0,{{0,1},{1,2},{2,3},{3,4},{4,0}}) => "5000",
(4,1,{{0,1},{1,2},{2,3},{3,0},{0,4}}) => "4100_3",
(2,4,{{0,1},{0,2},{0,3},{1,4},{1,5}}) => "3200_3",
(5,0,{{0,1},{0,2},{0,3},{0,4},{1,2,3,4}}) => "3200_4",
(1,5,{{0,1},{0,2},{0,3},{0,4},{0,5}}) => "2300_5"
}
FF = QQ
for p in LIST do (
    isParametric = true;
    minimal := isMinimal first p;
    << "isMinimal = " << minimal << endl; 
    print "-------------------------------------------------------";
    isParametric = false;
    computeDimDegreeF4(first p, FileName=>DIR|last p|".txt", Saturate=>true) 
    )
end --------------------------------------------------------

restart
load "dim-degree-F4-2cameras.m2"

--m = 2
--D = (5,0,{{0,1},{1,2},{2,3},{3,4},{4,0}}) -- 20, monodromy works
--D = (4,1,{{0,1},{1,2},{2,3},{3,0},{0,4}}) -- monodromy gets 16 -- saturation?
--D=(2,4,{{0,1},{0,2},{0,3},{1,4},{1,5}}) -- monodromy gets 12 -- saturation?
--D=(3,2,{{0,1},{0,2},{2,0},{0,3},{0,4}}) -- fails rank check
--D=(1,5,{{0,1},{0,2},{0,3},{0,4},{0,5}}) -- fails rank check
