-- builds a NONsymmetric problem from diagram D and list of non-visible points NV   

load "vision-toolbox.m2"
-- (lines, ghosts, list: per visible point, the indices of intersecting lines)
(nLines,nGhosts,lineIncidences) = D
linesNotVisibleSomewhere = flatten lineIncidences_(flatten NV)
linesVisibleEverywhere = toList(set(0..nLines-1) - set linesNotVisibleSomewhere)
nVisiblePoints = #lineIncidences
controlPoints = select(nLines,i->#select(lineIncidences,l->member(i,l))==1) 
setup3cameras(); -- R (ring), C (cameras); H2, H3 (skew symmetric matrices for Cayley)  
P = apply(#C,i->newPoints nVisiblePoints) 
Q = apply(#C,i->newPoints(#controlPoints)) 

L = apply(#C,i -> matrix( 
	apply(nLines+nGhosts, l->{
		randomLineThroughPoints(
		    if member(l,controlPoints) then 
		    P#i_{position(lineIncidences,x->member(l,x))} | Q#i_{position(controlPoints,x->x==l)}
		    else P#i_(positions(lineIncidences,x->member(l,x))) 
		    )
		}) 
	))
I = sum(#lineIncidences, 
    p->commonPoint(C,apply(#L,i->(L#i)^(if member(p,NV#i) then {} else lineIncidences#p)))
    ) + sum(linesVisibleEverywhere, i->seeLine(C,L/(l->l^{i})));

Rz = FF[z, gens R, MonomialOrder=>Eliminate 1];
Jz = sub(I,Rz) + ideal(1+z * product(skewSyms/(S -> 
	    sub(det(1-S),Rz) )));
end
