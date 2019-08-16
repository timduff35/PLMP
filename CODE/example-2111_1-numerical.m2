-- EXAMPLE of numerical degree computation for problem 2111_1 (minimal, degree=40) 
restart
setRandomSeed 0
m = 3
D = (3,2,{{1,2},{1,3},{1,4}})
Jpivots = {0, 1, 4, 5, 8, 9, 14, 29, 33, 44, 48, 57, 58, 59}
RERUNMONODROMY = false;
JACOBIAN = false;
COFACTOR = true;
needs "numerical-problem-builder.m2"

(yStart, cStarts) = readStartSys "2111_1-start";
max(cStarts/(c->norm evaluate(F,c||yStart)))
netList(cStarts/transpose)

-- fabrication
(cTargetList, PLTarget, CTarget) =  fabricatePair(D,RR,nvars);
rotations23 = {take(cTargetList,{0,3}),
              take(cTargetList,{4,7})}/(c->
	      Q2R((1/norm(2,c))*c,Normalized=>true));
translations23 = {{take(cTargetList,{8,10})},
                  {take(cTargetList,{11,13})}}/matrix/transpose;
L=PLTarget/last;
netList rotations23
netList translations23
netList L

-- encode fabricated data as (parameter, solution) pair
(yTarget, cTarget) = encodeyc(cTargetList, PLTarget, CTarget, RR)

--homotopy continuation
PH = parametricSegmentHomotopy GS
PH' = specialize(PH,yStart||transpose matrix yTarget);
cTargets = trackHomotopy(PH', cStarts);

--filtering
cTargets2 = cTargets/(c-> (
	R2params = take(c.Coordinates,{0,3});
	Q2R(R2params,Normalized=>true)
	)
    );	
i=minPosition(
    cTargets2/(R2->
	norm(R2-rotations23#0)));
y=(cTargets#i).Coordinates;
Q2R(take(y,{0,3}),Normalized=>true)
Q2R(take(y,{4,7}),Normalized=>true)
transpose matrix {take(y,{8,13})}
