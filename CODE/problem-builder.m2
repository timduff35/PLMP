load "problem-builder-core.m2"
I = sum(lineIncidences, LINES->commonPoint(C,L/(l->l^LINES)))+ sum(nLines, i->seeLine(C,L/(l->l^{i})));

Rz = FF[z, gens R, MonomialOrder=>{Eliminate 1}];
Jz = sub(I,Rz) + ideal(1+z * product(skewSyms/(S -> 
	    sub(det(1-S),Rz) )));

end