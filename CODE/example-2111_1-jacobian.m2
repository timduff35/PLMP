restart
setRandomSeed 0;
m = 3;
FF = ZZ/nextPrime 10000
isParametric = true;
D = (3,2,{{1,2},{1,3},{1,4}});
needs "problem-builder-matrices.m2"
matrices = pointMatrices | lineMatrices;
(numgens R, numgens coefficientRing R)
xy = fabricatedVectorFLPQ D;
gm = goodMinors(pointMatrices,lineMatrices,xy);
#gm
