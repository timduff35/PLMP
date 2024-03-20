-- IMPORTS 
debug needsPackage "NumericalAlgebraicGeometry"
debug needsPackage "MonodromySolver"
debug SLPexpressions
debug needsPackage "Core"

-- FUNCTIONS

size GateMatrix := M -> (numrows M, numcols M)
size Matrix := M -> (numrows M, numcols M)

-*
-- evaluate a gateMatrix G at a matrix x
-- don't use this a lot...
evaluate (GateMatrix, Matrix) := (G, x) -> (
    M := mutableMatrix(FF,numrows G,numcols G);
    E' := makeEvaluator(G,matrix{cameraVars|dataParams});
    evaluate(E',mutableMatrix(x),M);
    matrix M
    )
*-

--random diagonal matrix
randDiag = n -> diagonalMatrix for i from 1 to n list random FF

dehomogenize = method(Options=>{})
dehomogenize (Matrix, ZZ) := o -> (v, n) -> (
    --assumes column vector
    (1/v_(n, 0))*v^(toList splice(0..n-1,n+1..numrows v-1))
    )
dehomogenize Matrix := o -> v -> dehomogenize(v, numrows v -1)

summary = L -> (
    n := #L;
    H := sort L;
    Q1 := (1/2) * (H#(floor((n-1)/3))+H#(ceiling((n-1)/3)));
    med := (1/2) * (H#(floor((n-1)/2))+H#(ceiling((n-1)/2)));
    Q3 := (1/2) * (H#(floor(2*(n-1)/3))+H#(ceiling(2*(n-1)/3)));
    mean := (sum L)/n;
    var := sum(L/(x-> (x - mean)^2))/(n-1);
    << "Min: " << toString(min L) << endl;
    << "1Q: " << toString(Q1) << endl;
    << "Med: " << toString(med) << endl;
    << "Avg: " << toString(sub(mean,RR)) << endl;
    << "3Q: " << toString(Q3) << endl;
    << "Max: " << toString(max L) << endl;
    << "Std Dev: " << toString(sqrt(var)) << endl;
    )    

-- random element in the kernel of M
randKernel = method(Options=>{Tolerance=>1e-4})
randKernel (Matrix, InexactFieldFamily) := o -> (M, FF) -> (
    K := numericalKernel(M, Tolerance=>o.Tolerance);
    K*random(FF^(numcols K), FF^1)
    )
randKernel Matrix := o -> M -> randKernel(M, FF, oo)

reshapeCol = p -> if (numrows p == 1) then transpose p else p
reshapeRow = p -> if (numcols p == 1) then transpose p else p

-- RANDOMIZATION FOR PARAMETER POINT p
-- assumes random CC has unit modulus!
gammify = method()
gammify Point := p -> gammify reshapeCol matrix p
gammify Matrix := p -> (
    gammas := for i from 1 to m*#indepLines list random FF;
    indDiag := flatten(gammas/(g->{g,g,g}));
    -- abstract to arbitrary diagram, number of cameras
    depWInd := depLines/(l->first select(1,last D,i->member(l,i)));
    mfoldIntersections := flatten(depWInd/(l->for i from 0 to m-1 list l/(x->m*x+i)));
    -- next line assumes indepenedent lines come first!
    depDiag := flatten(
	mfoldIntersections/(ind -> {
	    conjugate gammas#(ind#0), 
	    conjugate gammas#(ind#1)
	    }
	)
    );
    tChartdiag := toList((3*(m-1)+1):random(FF)); -- t chart gamma
    qChartDiags := flatten for i from 0 to m-2 list toList(5:random(FF));
    p' := diagonalMatrix(indDiag|depDiag|tChartdiag|qChartDiags)*p;
    p'
    )


-- fold along rows
rfold = L -> if (#L ==0) then random(FF^0,FF^0) else fold(L, (a,b) -> a||b)

-- fold along cols
cfold = L -> fold(L, (a,b) -> a|b)


-- write starting parameters and solutions to file
writeStartSys = method(Options=>{Filename=>"startSys"})
writeStartSys (Matrix, List) := o -> (M, sols) -> writeStartSys(point M, sols, o)
writeStartSys (Point, List) := o -> (p, sols) -> (
   assert(instance(o.Filename,String));
   f := openOut o.Filename;
   f << "Parameter values: " << endl;
   f << toExternalString p << endl;
   f << "Solutions : " << endl;
   for s in sols do f << toExternalString s << endl;
   close f;
   )

readStartSys = filename -> (
    l := separate("\n", get filename);
    p0 := value l#1;
    sols := for i from 3 to #l-2 list value l#i;
    (transpose matrix p0, sols/matrix/transpose)
    )

-- for testing the contents of a start system file
startSysTester = (p,sols) -> (
    p0 := (transpose matrix V.BasePoint);
    p1 := random(FF^(#dataParams),FF^1);
    P01 = p0||p1;
    Pspec01 := specialize(PH,P0);
    target01 := trackHomotopy(Pspec01, sols);
    Pspec10 := (gammify p1)|(gammify p0);
    trackHomotopy(Pspec10, target01)
    )
    

adjugate = method()
adjugate Thing := M -> (
    n := numcols M;
    assert(n == numrows M);
    matrix table(n,n,(i,j)->((-1)^(i+j))*det submatrix'(M,{j},{i}))
    )

-- not printing to high precision -- deprecated?
sol2String = p -> replace("\\{|\\}","",toString p.Coordinates)

-- produces gates for "small" determinants"
det2 = M -> M_(0,0)*M_(1,1)-M_(1,0)*M_(0,1)
det3 = M -> M_(0,0)*det2(M_{1,2}^{1,2})-M_(0,1)*det2(M_{0,2}^{1,2})+M_(0,2)*det2(M_{0,1}^{1,2})
det4 = M -> M_(0,0)*det3(M_{1,2,3}^{1,2,3})-M_(0,1)*det3(M_{0,2,3}^{1,2,3})+M_(0,2)*det3(M_{0,1,3}^{1,2,3})-M_(0,3)*det3(M_{0,1,2}^{1,2,3})

laplaceDet = M -> (
    (m, n) := size M;
    if (m=!=n) then error("not square matrix")
    else if (m>5) then error("no Laplace for matrices larger than 4x4")
    else if (m==2) then det2 M
    else if (m==3) then det3 M
    else -* m==4 *- det4 M
    )

-- jacobian of GateMatrix wrt. a list of inputGates
jacobian (GateMatrix, List) := (F,inGates) -> fold(apply(inGates,g->diff(g,F)),(a,b)->a|b)

-- get rotation matrix from cayley parameters
cay2R = method(Options=>{Normalized=>false})
cay2R (Thing,Thing,Thing) := o -> (X,Y,Z) -> (
    if instance(X, RR) then x := X_FF else x = X;
    if instance(Y, RR) then y := Y_FF else y = Y;
    if instance(Z, RR) then z := Z_FF else z = Z;
    M := matrix{
    {1+x*x-(y*y+z*z), 2*(x*y-z), 2*(x*z+y)},
    {2*(x*y+z), 1+y^2-(x*x+z*z), 2*(y*z-x)},
    {2*(x*z-y), 2*(y*z+x), 1 +z*z -(x*x+y*y)}
	};
    if o.Normalized then (1/(1+x^2+y^2+z^2)) * M else M
    )
cay2R List := o -> L -> cay2R(L#0, L#1, L#2, o)

-- get Cayley parameters from rotation matrix
R2Cay = method(Options=>{UnNormalize=>false})
R2Cay Matrix := o -> R -> (
    assert(numcols R == 3);
    assert(numrows R == 3);
    S := (R-id_(FF^3))*(R+id_(FF^3))^-1;
    (S_(2,1), S_(0,2), S_(1,0))
    )

-*/// TEST
restart
needs "common.m2"
(x, y, z) = (random RR, random RR, random RR)
R = cay2R(x, y, z)
(x',y',z') = R2Cay R
R = cay2R(x', y', z')
R2Cay R
///*-

-- get rotation matrix from quaternion parameters
Q2R = method(Options=>{Normalized=>false, FF=>FF})
Q2R (Thing,Thing,Thing, Thing) := o -> (W, X,Y,Z) -> (
    if instance(W, RR) then w := W_FF else w = W;
    if instance(X, RR) then x := X_FF else x = X;
    if instance(Y, RR) then y := Y_FF else y = Y;
    if instance(Z, RR) then z := Z_FF else z = Z;
    M := matrix{
    {w*w+x*x-(y*y+z*z), 2*(x*y-w*z), 2*(x*z+w*y)},
    {2*(x*y+w*z), w^2+y^2-(x*x+z*z), 2*(y*z-w*x)},
    {2*(x*z-w*y), 2*(y*z+w*x), w^2 +z*z -(x*x+y*y)}
	};
    if o.Normalized then (1/(w^2+x^2+y^2+z^2)) * M else M
    )
Q2R List := o -> L -> Q2R(L#0, L#1, L#2, L#3, o)

-- get Cayley parameters from rotation matrix
R2Q = method(Options=>{UnNormalize=>false,FF=>FF})
R2Q Matrix := o -> R -> (
    assert(numcols R == 3);
    assert(numrows R == 3);
    c := (R_(2,1) - R_(1,2));
    b := (R_(0,2) - R_(2,0));
    a := (R_(1,0) - R_(0,1));
    w := (1/2)*sqrt(R_(0,0)+R_(1,1)+R_(2,2)+1);
    x := 1/(4*w) * c;
    y := 1/(4*w) * b;
    z := 1/(4*w) * a;
--    << w^2+x^2+y^2+z^2 << endl;
    (w, x, y, z)
    )

-*/// TEST
R=FF[W]
netList solveSystem {W^4-W^2+1/16}

clean T
T=QQ[a..d]
R=Q2R gens T
S = (R-id_(((QQ)^3)))*adjugate(R+id_((QQ)^3));
S
((first x)/(first L))*L
1/sqrt(sum(x/(y->y^2)))*x
L
///*-


-- cross product of col vectors -- takes Matrice or GateMatrix pair
crossProduct = (y,q) -> matrix{{y_(1,0)*q_(2,0)-y_(2,0)*q_(1,0)},{y_(2,0)*q_(0,0)-y_(0,0)*q_(2,0)},{y_(0,0)*q_(1,0)-y_(1,0)*q_(0,0)}}

cmat = method(Options=>{Normalized=>true})
cmat List := o -> L -> cmat(L#0,L#1,L#2,o)
cmat (Thing,Thing,Thing) := o -> (a,b,c) -> (
    tx := matrix{{0,-c,b},{c,0,-a},{-b,a,0}};
    if o.Normalized and not areEqual(L#2,0.0) then tx = (-1/L#2) * tx;
    tx
    )

--
randomLineThroughPoints = (P, FF) -> ( 
    m := numrows P; -- m = dim + 1
    n := numcols P; -- n = number of points
    assert(m>=3 and m<=4);
    K := numericalKernel(transpose P,Tolerance=>1e-6);
    --assert(numcols K == m-n); -- true if points are distinct
    transpose(K * random(FF^(numcols K),FF^(m-2)))
    )

-- constructs problem data given a PL diagram D (complete visibility)
-- returns (camera parameters, lines, camera matrices
-- ASSUMES: "intersections" are sorted
fabricatePair = (D, FF, nvars) -> (    
    (nLines,nGhosts,intersections) := D;
    depPoints := set {};
    scan(nLines,l->(
	    ptsOnl := positions(last D,i->member(l,i));
	    depPoints = depPoints + set drop(ptsOnl,2);
	    )
	);    
    pointsOnLineIndices := apply(nLines+nGhosts, l->positions(intersections,i->member(l,i)));
    worldPointsFF := random(FF^4,FF^0);
    scan(#intersections,i->(
	    pointi := if member(i,depPoints) then (
		li := first select(1,pointsOnLineIndices,l->member(i,l));
		<< li << endl;
		a := random FF;
		a*worldPointsFF_{li#0}+(1-a)*worldPointsFF_{li#1}
		) else random(FF^4,FF^1);
	    worldPointsFF = worldPointsFF | pointi
	    )
	);
    worldPoints := sub(worldPointsFF,FF);    
    helperPoints := apply(pointsOnLineIndices, pp->sub(random(FF^4,FF^(2-min(2,#pp))),FF));
    -- future (line below): may be interesting to sample space of variables differently depending on the field we fabricate data over
    sampleCameraParameters := for i from 1 to nvars list sub(random FF,FF);
    subTable := apply(sampleCameraParameters, cameraVars, (a,b) -> b=>inputGate a);
    sampleC := apply(C,cam -> (
	    M := mutableMatrix(FF, 3, 4);
	    evaluate(cam, mutableMatrix{sampleCameraParameters}, M);
	    matrix M
	    )
	    );
    (
	sampleCameraParameters,
	apply(sampleC, cam->(
	    	P := cam * worldPoints;
	    	L := matrix apply(nLines+nGhosts, l->(
		      Hl := helperPoints#l;
		      Pl := P_(pointsOnLineIndices#l);
		      Pl = Pl|cam*Hl;
	    	      line := randomLineThroughPoints(Pl, FF);
		      {(1/norm(2,line))*line} -- normalization that seemed to benefit the chicago solver
		       ));
	       (P,L)
	       )),
       sampleC
	) 
    )

-- functions which fabricate input of the form (parameter, solution)
-- notation for supp materials is (y,c)---previous notation is (p,x)
encodey = (P,L,projs,FF) -> (
    c := transpose matrix{P};
    allLines := L/last;
    yIndLineMatrix := cfold(
	    allLines/(m->m^(toList indepLines))
	    );
    yInd := matrix(yIndLineMatrix,product toList size yIndLineMatrix,1);
    yDep := if (#depLines == 0) then random(FF^0,FF^1) else 
    rfold(
	depLines/(l-> (
	    	lineInds := first select(1,last D,i->member(l,i));
		triplet := take(lineInds, 2) | {l};
	    	rfold(allLines/(m -> (
	    		n :=numericalKernel(transpose m^triplet, Tolerance=>kTol);
	    		(1/n_(2,0))*n^{0,1}
			))
	    	    )
		)
    	    )
	);
    ytChart := sub(randKernel(transpose(c^{4*(m-1)..(4*(m-1)+3*(m-1)-1)}||matrix{{1_FF}}), FF),FF);
    yqChart := sub(rfold(
	for i from 0 to m-2 list randKernel(transpose(c^{4*i..4*i+3}||matrix{{1_FF}}), FF)
	),FF);
    yInd||yDep||ytChart||yqChart
    )

encodeyc = (P, L, projs,FF) -> (
    c := transpose matrix{P};
    y := encodey(P, L, projs,FF);
    (point y, point c)    
    )    

fabricateyc = FF -> (
    (P, L, projs) := fabricatePair(D, FF, nvars); -- these variable names are confusing
    encodeyc(P, L, projs,FF)
    )

-- convenience functions for minors
minors (GateMatrix, ZZ, Sequence, Boolean) := o ->  (M, k, S, laplace) -> (
    (Sm, Sn) := (first S, last S);
    (m,n) := (numrows M, numcols M);
    assert(k<=min(m,n));
    assert(all(Sm,s->#s==k));
    assert(all(Sn,s->#s==k));
    flatten apply(Sm,sm->apply(Sn, sn -> 
	    if (laplace) then laplaceDet submatrix(M,sm,sn)
	    else det submatrix(M,sm,sn)
	    ))
    )

allMinors = method(Options=>{Laplace=>false})
allMinors (GateMatrix, ZZ) := o -> (M, k) -> (
    (m, n ) := (numrows M, numcols M);
    s := (subsets(0..m-1,k),subsets(0..n-1,k));
    minors(M, k, s, o.Laplace)
    )

maxMinors = method(Options=>{Laplace=>false})
maxMinors GateMatrix := o -> M -> allMinors(M,min(numrows M, numcols M), Laplace=>o.Laplace)

-- this seems to work
complexQR = M -> (
    A := mutableMatrix M;
    k := ring A;
    Q := mutableMatrix(k,0,0,Dense=>true);
    R := mutableMatrix(k,0,0,Dense=>true);
    rawQR(raw A, raw Q, raw R, true);
    assert(areEqual(Q*R,A)); -- idk if it will work every time!
    (matrix Q,matrix R)
    )

leverageScores = M -> (
    Q = first complexQR M;
    rsort apply(numrows Q,i->(norm(2,Q^{i}),i))
    )

leverageScoreRowSelector = J0 -> (
    sortedRows := (leverageScores J0)/last;
    r := rowSelector J0^(sortedRows);
    sort(r/(i->sortedRows#i))
    )

log10 = x -> log(x)/log(10)

argFF = z -> atan((imaginaryPart z)/(realPart z))

-- complex number whose real and imag parts are standard normal
gaussFF = () -> (
    (u1,u2):=(random RR,random RR);
    sqrt(-2*log(u1))*cos(2*pi*u2)+ii*sqrt(-2*log(u1))*sin(2*pi*u2)
    )

-- random sample drawn from normal distriution N(mu, var^2)
rNorm = (mu,var) -> mu+var*(realPart gaussFF())_FF

-- random sample from (n-1)-sphere with radius r
sphere = (n,r) -> (
    l:=apply(n,i->rNorm(0,1));
    matrix{r/norm(2,l)*l}
    )

-- assumes "u" of unit length
householder=method()
householder (Matrix,ZZ) := (u,n) -> (
    if (numrows u > 1) then error("householder takes a row vector");
    R:=ring u;
    k:=numcols u;
    id_(R^(n-k))++(id_(R^k)-2*(transpose u)*u)
    )
householder (List,ZZ) := (u,n) -> householder(matrix {u},n)

randomOn = n -> diagonalMatrix(toList((n-1):1_RR)|{(-1)^(random 2)}) * fold(reverse apply(2..n,i->householder(sphere(i,1),n)),(a,b)->a*b)

randomCameraNormalized = () -> (
    R := randomOn 3;
    t := matrix{{random FF},{random FF},{random FF}};
--    t := transpose matrix{sphere(3,1)};
    tnorm := (1 / t_(2,0))* t;
    (R|tnorm)
    )

randomCameraNormalizedCayley = () -> (
    R := cay2R(random FF, random FF, random FF,Normalized=>true);
    t := matrix{{random FF},{random FF},{random FF}};
--    t := transpose matrix{sphere(3,1)};
    tnorm := (1 / t_(2,0))* t;
    (R|tnorm)
    )


randomCamera = () -> (
    R := randomOn 3;
    t := transpose matrix{sphere(3,1)};
    (R|t)
    )

ranks = method(Options=>{})
ranks (Matrix, Matrix) := o -> (x, p) -> (
    if (numcols x > 1) then x = matrix(x, numcols x,1);
    if (numcols p > 1) then p = matrix(p, numcols p,1);
    xpp := mutableMatrix sub(x||p,CC);
    a := PE/( m -> (
	    evaluate(first m, xpp, last m);
	    numericalRank matrix last m
	    )
	    );
    b := LE/( m -> (
	    evaluate(first m, xpp, last m);
	    numericalRank matrix last m
	    )
	    );
   (a, b)
   )
ranks (Point,Point) := o -> (x,p) -> ranks(matrix x, matrix p)

rankCheck = method(Options=>{Hard=>true})
rankCheck (Matrix, Matrix) := o -> (x, p) -> (
   (a, b) := ranks(x, p);
   if (o.Hard) then (all(a,x->x==3) and all(b,x->x==2))
     else (all(a,x->x<=3) and all(b,x->x<=2))
   )

cpMatrix = t -> matrix{{0,-t_(2,0),t_(1,0)},{t_(2,0),0,-t_(0,0)},{-t_(1,0),t_(0,0),0}}

essential = (R,t) -> R * cpMatrix t
