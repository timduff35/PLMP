restart
setRandomSeed 0
m = 3 -- number of cameras
FF = ZZ/nextPrime 10000
isParametric = false;
D = (1,5,{{0,1},{0,2},{0,3},{4,5}}) -- encodes table entry 3100_0
needs "problem-builder-core.m2"
I = sum(lineIncidences,LINES->commonPoint(C,L/(l->l^LINES)))+sum(nLines,i->seeLine(C,L/(l->l^{i})));

-- we need cayley matrix denominators nonzero
-- rather than "saturate I ...", the following code exploits F4
Rz = FF[z, gens R, MonomialOrder=>{Eliminate 1}]
Jz = sub(I,Rz) + ideal(1+z*product(skewSyms/(S->sub(det(1-S),Rz))));
G = groebnerBasis(Jz, Strategy => "F4");

-- now eliminate z to study I
GBI = selectInSubring(1,G); -- GB for elimination ideal, Jz intersect FF[r_(i,j),t_i]
inI = ideal sub(leadTerm GBI,R)
dim inI
degree inI
