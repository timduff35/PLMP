-- A point-line diagram D is 
-- (#lines, #ghosts, {{a,b},{c,d},...}) where a,b,c,d are line numbers

-- Global variables:
-- 
-- FF, a field
-- R, the "master" ring
--    coefficientRing of R is either FF or FF[a_1..]
--    variables: r_(1,1).., c_1.., d_1..
-- (created by setup3cameras)
-- H2, H3 -- skew symmetric matrices produced by Cayley parametrization 
-- C = list of cameras
-- nParametersPerCamera -- controls how many a_s are in coefficientRing R
--                         and creates R.NextParameter 

MAX'NUMERATOR'DENOMINATOR=50
old'random'Type = lookup(random,Type)
random Type := o -> R -> (
    if R === QQ then (
	p := random(-MAX'NUMERATOR'DENOMINATOR,MAX'NUMERATOR'DENOMINATOR);
	q := random(1,MAX'NUMERATOR'DENOMINATOR);
	p/q
	) 
    else (old'random'Type o) R
    ) 


adjugate = method()
adjugate Matrix := M -> (
    n := numcols M;
    assert(n == numrows M);
    matrix table(n,n,(i,j)->((-1)^(i+j))*det submatrix'(M,{j},{i}))
    )

cayley = method()
cayley Matrix := H -> (
    assert(transpose H + H == 0);
    (1+H)*adjugate(1-H)
    )

seeLineMatrix = method() -- line-...-line correspondence
seeLineMatrix (List,List) := (C,L) -> (
    assert(#C==#L);
    assert all(C,c->(numrows c, numcols c) == (3,4));
    assert all(L,l->(numrows l, numcols l) == (1,3));
    matrix apply(C,L,(c,l)->{l*c})
    )

seeLine = method() -- line-...-line correspondence
seeLine (List,List) := (C,L) -> (
    ret := minors(3,seeLineMatrix(C,L),Strategy=>Cofactor);
    << "-- seeLine: #equations = " << numgens ret << endl;
    ret
    )

commonPointMatrix = method() -- L are images in C of several lines passing through the same point 
commonPointMatrix (List,List) := (C,L) -> (
    assert(#C>=2 and #C==#L);
    assert all(C,c->(numrows c, numcols c) == (3,4));
    assert all(L,l->numcols l == 3);
    assert(sum(L,numrows) >= 4);  
    matrix apply(C,L,(c,l)->{l*c}) 
    )

commonPoint = method() -- L are images in C of several lines passing through the same point 
commonPoint (List,List) := (C,L) -> (
    ret := minors(4,commonPointMatrix(C,L),Strategy=>Cofactor);
    << "-- commonPoint: #equations = " << numgens ret << endl;
    ret
    )

-- methods for saturation by lower rank minors
commonPointSat = method()
commonPointSat (List,List) := (C,L) -> (
    ret := minors(3,commonPointMatrix(C,L),Strategy=>Cofactor);--,Limit=>1); -- is Limit=>1 enough
    << "-- commonPoint saturation by " << numgens ret << " equations" << endl;
    ret
    )
seeLineSat = method() -- line-...-line correspondence
seeLineSat (List,List) := (C,L) -> (
    ret := minors(2,seeLineMatrix(C,L),Strategy=>Cofactor);--,Limit=>1); -- is Limit=>1 enough?
    << "-- commonLine saturation by " << numgens ret << " equations" << endl;
    ret
    )


setupcameras = m -> (
    Coeff := if not isParametric then FF 
    else FF[a_1..a_(m*nParametersPerCamera)];
    R = Coeff[r_(1,1)..r_(m-1,3),t_(1,1)..t_(m-1,1), t_(1,2)..t_(m-1,2), t_(2,3)..t_(m-1,3)];
    if isParametric then R.NextParameter = 0;
    skewSyms = for i from 1 to m-1 list genericSkewMatrix(R,r_(i,1),3);
    Rs = {map(R^3)} | skewSyms/cayley;
    ts = {
	-* random(R^3,R^1), -- this complicates some computations *-
	transpose matrix{{0_R,0,0}},
	-* random(R^3,R^3) * -- this does not... but does not seem to be necessary *-  
	transpose matrix{{t_(1,1),t_(1,2),1}}} | 
    for i from 2 to m-1 list transpose matrix {{t_(i,1),t_(i,2),t_(i,3)}};
    C = apply(Rs,ts,(R,t)->R|t);
    )

setup3cameras = () -> setupcameras 3

-* 
D = (lines, ghosts, incidence-list)

items of incidence-list correspond to (visible) points,
each item is a list of incident lines.

ASSUMPTIONS: 

(1) free lines occur BEFORE any other lines 
E.g: 
D = (4,1,{{2,3},{3,4}}) 
has 2 points, 1 pin, 2 free lines;
free lines 0,1 in precede 2,3,4.

(2) dependent points are points on the lines that are listed _third_ or later in incidence-list     
E.g: 
D = (2,4,{{0,1},{0,2},{0,3},{1,4},{1,5}}) 
has 5 points with 012 and 034 collinear.
Points 0,1,3 are independent and 2,4 are dependent.
Note that    
D = (2,4,{{1,5},{0,2},{0,3},{1,4},{0,1}}) 
would be not admissable diagram, since point 4 belongs to both visible lines and can't be dependent on 03 and 12.
*-
parseD = D -> (
    -- (lines, ghosts, list: per visible point, the indices of intersecting lines)
    (nLines,nGhosts,lineIncidences) = D;
    nVisiblePoints = #lineIncidences;
    dependentPoints = new MutableHashTable;
    scan(nLines,i->(
	    local A, local B;
	    c := 0;	
    	    scan(nVisiblePoints,j->if member(i,lineIncidences#j) then (
		    if c == 0 then A = j
		    else if c == 1 then B = j
		    else dependentPoints#j = {A,B};
		    c = c + 1;
		    )) 
	    ));
    controlPoints = select(nLines,i->#select(lineIncidences,l->member(i,l))==1); 
    freeLines = select(nLines,i->#select(lineIncidences,l->member(i,l))==0); 
    assert(toList(0..#freeLines-1) == freeLines) -- check assumption (1)
    )

newPoints = method()
newPoints ZZ := n -> ( 
    if (n>0 and R.?NextParameter) then (
	C := coefficientRing R;
	ret := sub(genericMatrix(C,C_(R.NextParameter),2,n) || matrix{toList(n:1)},R);
	R.NextParameter = R.NextParameter + 2*n;
    	ret
	)
    else random(FF^3,FF^n) 
    )

getNextParameterOrRandom = () -> if R.?NextParameter then (
    ret := (coefficientRing R)_(R.NextParameter);
    R.NextParameter = R.NextParameter + 1;
    return ret
    ) else random FF

newPoints (ZZ,MutableHashTable) := (n,dependentPoints) -> ( 
    ret := map(R^3,R^0,0);
    scan(n,i-> ret = ret | if dependentPoints#?i then (
	    a := getNextParameterOrRandom();
	    ret_(dependentPoints#i) * matrix{{1-a},{a}} 
	    ) else (
	    a1 :=getNextParameterOrRandom();
	    a2 :=getNextParameterOrRandom();
	    matrix{
	    {a1},
	    {a2},
	    {1}
	    }
	)
	); 	   
    ret
    )
    
newLines = n -> (
    if R.?NextParameter then (
	C := coefficientRing R;
    	ret := sub(genericMatrix(C,C_(R.NextParameter),2,n) || matrix{toList(n:1)},R);
	R.NextParameter = R.NextParameter + 2*n;
    	transpose ret
	)
    else random(FF^n,FF^3) 
    )

randomLineThroughPoints = P -> ( 
    m := numrows P; -- m = dim + 1
    n := numcols P; -- n = number of points
    assert(m>=3 and m<=4);
    K := gens ker transpose P;
    if n<3 then assert(numcols K == m-n) -- true if points are distinct
    else assert(numcols K == 1); -- points have to be collinear
    transpose(K * random(FF^(m-min(2,n)),FF^(m-2)))
    )

-- returns (camera parameters, fabricated configuration) where
-- "fabricated configuration" is the list of triples 
--    (free lines, visible points, control points) 
-- which are camera projections in PP^2 of a configuration in PP^3 
fabricatePairFLPQ = D -> (
    (nLines,nGhosts,lineIncidences) := D;
    parseD D;
    wFLa := random(FF^4,FF^(#freeLines)); 
    wFLb := random(FF^4,FF^(#freeLines)); 
    wP := map(FF^4,FF^0,0);
    scan(nVisiblePoints,i-> wP = wP | if dependentPoints#?i then (
	    a := random FF;
	    wP_(dependentPoints#i) * matrix{{1-a},{a}} 
	    ) else random(FF^4,FF^1)
	);
    wQ := random(FF^4,FF^#controlPoints); 
    sampleCameraParameters := random(FF^1,FF^(numgens R));
    sampleC := apply(C,cam->sub(cam,sampleCameraParameters));
    (sampleCameraParameters,apply(sampleC, cam->(
		if #freeLines == 0 then map(FF^0,FF^3,0) 
		else matrix apply(#freeLines, l->
		    {randomLineThroughPoints(cam*wFLa_{l}|cam*wFLb_{l})}
		    ),
		cam*wP, 
		cam*wQ
		))) 
    )

-- determines t such that C = (1-t) A + t B 
dependentPointCoordinate = (C,A,B) -> (
    t := (C_(0,0)-A_(0,0))/(B_(0,0)-A_(0,0)); -- (B-A) t = C-A
    if C != (1-t)*A + t*B then error "C in not on the line AB";
    t         
    )

-- dehomogenize w.r.t last coordinate
toAffine = method()
toAffine Matrix := M -> transpose matrix ( 
    entries transpose M / toAffine -- do this to columns   
    )
toAffine List := abc -> (
    c = last abc;
    drop(abc,-1) / (x->x/c)
    );
    
fabricatedVectorFLPQ = D -> (
    (nLines,nGhosts,lineIncidences) := D;
    parseD D;
    wFLa := random(FF^4,FF^(#freeLines)); 
    wFLb := random(FF^4,FF^(#freeLines)); 
    wP := map(FF^4,FF^0,0);
    scan(nVisiblePoints,i-> wP = wP | if dependentPoints#?i then (
	    a := random FF;
	    wP_(dependentPoints#i) * matrix{{1-a},{a}} 
	    ) else random(FF^4,FF^1)
	);
    wQ := random(FF^4,FF^#controlPoints); 
    sampleCameraParameters := random(FF^1,FF^(numgens R));
    sampleC := apply(C,cam->sub(cam,sampleCameraParameters));
    conf := apply(sampleC, cam->(
		if #freeLines == 0 then map(FF^0,FF^3,0) 
		else matrix apply(#freeLines, l->
		    {randomLineThroughPoints(cam*wFLa_{l}|cam*wFLb_{l})}
		    ),
		cam*wP, 
		cam*wQ
		)); 
    freeLines := conf / first;
    visible := conf / (c->c#1);
    control := conf / last;
    sampleCameraParameters | matrix { 
	entries transpose matrix {freeLines/transpose} / toAffine // flatten |
	flatten apply(m, c->flatten apply(nVisiblePoints, i->if dependentPoints#?i 
		then (
		    (A,B) := toSequence dependentPoints#i;
		    P := (conf#c#1);
		    {dependentPointCoordinate(P_{i}//toAffine, P_{A}//toAffine, P_{B}//toAffine)}
		    ) 
		else toAffine flatten entries (conf#c#1)_{i}
		)) | 
	entries transpose matrix {control} / toAffine // flatten  
    	}
    )

-- returns (camera parameters, fabricated configuration) where
-- "fabricated configuration" is the list of pairs (visible points, control points) 
-- which are camera projections in PP^2 of a configuration in PP^3 
fabricatePairControlPoints = D -> (
    (nLines,nGhosts,lineIncidences) := D;
    nVisiblePoints := #lineIncidences;
    controlPoints := select(nLines,i->#select(lineIncidences,l->member(i,l))==1); 
    wP = random(FF^4,FF^nVisiblePoints); 
    wQ = random(FF^4,FF^#controlPoints); 
    sampleCameraParameters := random(FF^1,FF^(numgens R));
    sampleC := apply(C,cam->sub(cam,sampleCameraParameters));
    (sampleCameraParameters,apply(sampleC, cam->(cam*wP, cam*wQ))) 
    )

fabricatedVectorControlPoints = D -> (
    (cams, conf) := fabricatePairControlPoints D;
    visible := conf / first;
    control := conf / last;
    cams | matrix{ flatten apply(
	    entries transpose matrix {visible|control}, 
	    abc->(
	    	(a,b,c) := toSequence abc;
	    	{a/c,b/c}
	    	))
    	}
    )

rfold = L -> if (#L ==0) then random(FF^3,FF^0) else fold(L,(a,b)->a||b)

-- returns (camera parameters, fabricated configuration) where
-- "fabricated configuration" is the list of triples (visible points, control points, free lines) 
-- which area camera projections in PP^2 of a configuration in PP^3 
fabricateAll = D -> (
    (nLines,nGhosts,lineIncidences) := D;
    nVisiblePoints := #lineIncidences;
    controlPoints := select(nLines,i->#select(lineIncidences,l->member(i,l))==1); 
    freeLines := select(nLines,i->#select(lineIncidences,l->member(i,l))==0); 
    wP = random(FF^4,FF^nVisiblePoints); 
    wQ = random(FF^4,FF^#controlPoints); 
    helperPoints := freeLines/(fl->random(FF^4,FF^(2)));
    sampleCameraParameters := random(FF^1,FF^(numgens R));
    sampleC := apply(C,cam->sub(cam,sampleCameraParameters));
    (sampleCameraParameters,apply(sampleC, cam->(cam*wP, cam*wQ, rfold(helperPoints/(h->randomLineThroughPoints(cam*h)))))) 
    )

fabricatedVectorAll = D -> (
    (cams, conf) := fabricateAll D;
    visible := conf / first;
    control := conf/(c->c#1);
    free := conf / last;
    cams | matrix{ delete(null, flatten apply(
	    (entries transpose matrix {visible|control})|(free/entries/flatten),
	    abc-> if (#abc>0) then (
	    	(a,b,c) := toSequence abc;
	    	{a/c,b/c}
	    	)))
    	}
    )

 
-- returns (camera parameters, fabricated configuration) where
-- "fabricated configuration" is list of elements in the form (points,lines)
-- which are camera projections in PP^2 of a configuration in PP^3 
fabricatePair = D -> (
    (nLines,nGhosts,intersections) := D;
    worldPoints := random(FF^4,FF^(#intersections));
    pointsOnLineIndices := apply(nLines+nGhosts, l->positions(intersections,i->member(l,i)));
    helperPoints := apply(pointsOnLineIndices, pp->random(FF^4,FF^(2-#pp)));
    sampleCameraParameters := random(FF^1,FF^(numgens R));
    sampleC := apply(C,cam->sub(cam,sampleCameraParameters));
    (
	sampleCameraParameters,
	apply(sampleC, cam->(
	    	P := cam * worldPoints;
	    	L := matrix apply(nLines+nGhosts, l->(
	    	       line := randomLineThroughPoints(P_(pointsOnLineIndices#l)|(cam*helperPoints#l));
		       {(1/line_(0,2))*line} -- make last coordinate = 1 
		       ));
	       (P,L)
	       ))
	) 
    )
    
-- returns a list of maximal minors whose jacobian does not vanish at x 
-- OLD
goodMinors = (matrices,x) -> (
    labels := {};
    J := null;
    scan(#matrices, i->(
	    (cols,jA) := jacobianMaxMinors(matrices#i,x);
	    labels = labels | apply(cols, c->i=>c);
	    if J===null then J = jA else J = J||jA;
	    ));
    labels_(first \ pivots gens gb J)
    )

-- returns a list of maximal minors whose jacobian does not vanish at x 
goodMinors = (pointMatrices,lineMatrices,x) -> (
    labels := {};
    J := null;
    scan(#pointMatrices, i->(
	    (cols,jA) := jacobianMaxMinors(pointMatrices#i,x);
	    labels = labels | apply(cols, c->i=>c);
	    if J===null then J = jA else J = J||jA;
	    ));
    scan(#lineMatrices, i->(
	    (cols,jA) := jacobianMinors(lineMatrices#i,x,3);
	    labels = labels | apply(cols, c->i=>c);
	    if J===null then J = jA else J = J||jA;
	    ));
    labels_(first \ pivots gens gb J)
    )


-- returns (labels for rows, jacobian matrix of size-k minors of A evaluated at x)
jacobianMinors = (A,x,k) -> (
    rowMinorSet := subsets(numrows A, k);
    colMinorSet := subsets(numcols A, k);
    minorSet := flatten(
	rowMinorSet/(r -> 
	    colMinorSet/(c -> (r,c)
	    	)
	    )
	);
    vA := sub(A,x);
    vAx := apply(gens ring A, X->sub(diff(X,A),x));
    minorSet, transpose matrix apply(numgens ring A, i->apply(minorSet, minor->(
		(rows, cols) := (minor#0, minor#1);
	    	M := entries vA_cols^rows;
	    	Mx := entries (vAx#i)_cols^rows;
		sum(#M, j->det matrix replace(j,Mx#j,M)) 
		)) 
	)
    )   


-- returns (labels for rows, jacobian matrix of max minors of A evaluated at x)
jacobianMaxMinors = (A,x) -> (
    if numrows A > numcols A then A = transpose A;
    jacobianMinors(A, x, numrows A)
    )   



end

-- SOME TESTS --
restart
load "vision-toolbox.m2"
FF = ZZ/14253954677
R = FF[a,b,c]
H = genericSkewMatrix(R,3)
adjugate (1-H) * (1-H) == det(1-H)
cayley H

C = apply(3,i->random(FF^3,FF^4))
L = apply(3,i->random(FF^6,FF^3))
seeLine(C,L/(l->l^{0}))

