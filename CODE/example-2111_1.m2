-- EXAMPLE of computations for problem 2111_1 (minimal, degree=40) 
restart
setRandomSeed 0
m = 3 -- number of cameras
FF = ZZ/nextPrime 10000 -- a finite field
isParametric = false; -- parameters are specialized to generic numbers

D = (3,2,{{1,2},{1,3},{1,4}}) -- encodes problem 2111_1 
-* line 0 is free
   line 1 passes through points 0,1,2
   line 2 passes through point 0
   line 3 passes through point 1 (ghost line)
   line 4 passes through point 2 (ghost line)
*-    
needs "problem-builder-core.m2"
I = sum(lineIncidences,
    LINES->commonPoint(C,L/(l->l^LINES)) -- CP constraint for each point
    ) + sum(nLines,
    i->seeLine(C,L/(l->l^{i})) -- LC constraint for each (non-ghost) line
    );

-* Cayley matrix denominators do not vanish:
   auxiliary z ensures "invertibility" of these.
*-
Rz = FF[z, gens R, MonomialOrder=>{Eliminate 1}]
Jz = sub(I,Rz) + ideal(1+z*product(skewSyms/(S->sub(det(1-S),Rz)))); -- ensure Cayley denominators are nonzero
G = groebnerBasis(Jz, Strategy => "F4"); -- exploit the F4 algorithm 

-- now eliminate z to study the original ideal I
gbI = selectInSubring(1,G); -- Groebner basis of I, which equals Jz intersected with  FF[r_(i,j),t_i]
inI = ideal sub(leadTerm gbI,R) -- initial ideal of I
dim inI 
degree inI
