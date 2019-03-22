-- builds a NONsymmetric problem in _3 cameras_ 
-- from diagram D and list of non-visible points NV   

load "vision-toolbox.m2"
(nLines,nGhosts,lineIncidences) = D -- see problem-builder-core.m2 for description
-*
 NV is a list of lists of points hidden in corresponding cameras 

E.g: NV = {{0},{},{1}}
means point 0 does not show up in the first image and 1 in occluded from the last camera.
*-

linesNotVisibleSomewhere = flatten lineIncidences_(flatten NV)
linesVisibleEverywhere = toList(set(0..nLines-1) - set linesNotVisibleSomewhere)
nVisiblePoints = #lineIncidences
nPointsOnLines = apply(nLines+nGhosts,i->#select(lineIncidences,l->member(i,l)))
controlPoints = select(nLines,i->nPointsOnLines#i==1) 
nFreeLines = #select(nLines,i->nPointsOnLines#i==0)
<< "--building: " << nVisiblePoints << " points, " << #controlPoints << " pins, " << nFreeLines << " free lines." << endl 

setup3cameras(); -- R (ring), C (cameras); H2, H3 (skew symmetric matrices for Cayley)  
FL = apply(#C,i->newLines(nFreeLines))
P = apply(#C,i->newPoints nVisiblePoints) 
Q = apply(#C,i->newPoints(#controlPoints))

L = apply(#C,i -> matrix( 
	apply(nLines+nGhosts, l->{
		if nPointsOnLines#l==0 then FL#i^{l} 
		else randomLineThroughPoints(
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
