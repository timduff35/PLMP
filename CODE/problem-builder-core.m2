load "vision-toolbox.m2"
parseD D;

<< "--building: " 
<< nVisiblePoints << " points, " 
<< #controlPoints << " pins, " 
<< #freeLines << " free lines." << endl 
<< "--          dependent points:" << endl
print peek dependentPoints 

if isParametric then nParametersPerCamera = 2*(nVisiblePoints+#controlPoints+#freeLines)- #dependentPoints; --correct the number of parameters --- need to rewrite this!
setupcameras m; -- R (ring), C (cameras)

FL := apply(#C, i->newLines(#freeLines)); 
P = apply(#C,i->newPoints(nVisiblePoints, dependentPoints)) 
Q = apply(#C,i->newPoints(#controlPoints)) 
L = apply(#C,i -> matrix( 
	apply(nLines+nGhosts, l->{
		if member(l,freeLines) then FL#i^{l}
		else randomLineThroughPoints(
		    if member(l,controlPoints) then 
		    P#i_{position(lineIncidences,x->member(l,x))} | Q#i_{position(controlPoints,x->x==l)}
		    else P#i_(positions(lineIncidences,x->member(l,x))) 
		    )
		}) 
	))

