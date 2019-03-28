restart
setRandomSeed 0;
m = 3;
FF = ZZ/nextPrime 10000
isParametric = false;
D = (3,2,{{1,2},{1,3},{1,4}});
needs "problem-builder.m2"
G = groebnerBasis(Jz, Strategy => "F4"); -- exploit the F4 algorithm 

-- now eliminate z to study the original ideal I
gbI = selectInSubring(1,G); -- Groebner basis of I, which equals Jz intersected with  FF[r_(i,j),t_i]
inI = ideal sub(leadTerm gbI,R) -- initial ideal of I
dim inI 
degree inI
