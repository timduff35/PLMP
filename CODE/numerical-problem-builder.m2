-- IMPORTS
needs "common.m2"
needsPackage "MonodromySolver"

-- GLOBALS
FF = CC
kTol=1e-4
COFACTOR = true

nLines = D#0 + D#1
depLines = flatten((last D)/(i -> drop(i,2)))
indepLines = sort toList(set(0..nLines-1)-set(depLines))


-*
-- SET TRACKER OPTIONS HERE
-- null indicates default value 
scan({CorrectorTolerance=>null,
	EndZoneFactor=>null,
	InfinityThreshold => null, 
	maxCorrSteps => null, 
	NoOutput => null,
	numberSuccessesBeforeIncrease => null,
	Precision => null,
	Predictor => null,
	stepIncreaseFactor => null,
	tStep => null,
	tStepMin => null
	}, 
    opt -> setDefault(opt))
*-

-- setting up input gates
-- "l" indices go (view, feature, coordinate)
varInSymbs = (flatten (for i from 2 to m list {w_i,x_i,y_i,z_i}))|
    (flatten (for i from 2 to m list {t_(i,1),t_(i,2),t_(i,3)}))

nvars = #varInSymbs -- should change variable name
inSymbs=varInSymbs|
    (flatten for j in indepLines list flatten for i from 1 to m list for k from 1 to 3 list l_(i,j,k))|
	(flatten for j in depLines list flatten for i from 1 to m list for k from 1 to 2 list a_(i,j,k))|
	(flatten (for i from 2 to m list(
		     {ct_(i,1), ct_(i,2), ct_(i,3)}
		     )))|{ct_0}|
	 (flatten (for i from 2 to m list toList(cq_(i,0)..cq_(i,4))))
scan(inSymbs,s->value(toString(s)|" = inputGate " | toString(s)))
inGates = toList(inSymbs)/value

-- variable groups
cameraVars = take(inGates,nvars)
centerVars = drop(cameraVars,4*(m-1))
qVars = take(cameraVars,4*(m-1))
varMatrix = matrix{cameraVars}

-- parameter groups
dataParams = drop(inGates,nvars)
chartParams = drop(dataParams,
    #dataParams-(3*(m-1)+1+5*(m-1)))
tchartParams = take(chartParams,3*(m-1)+1)
qchartParams = drop(chartParams,3*(m-1)+1)
paramMatrix = matrix{dataParams}
inputMatrix = varMatrix | paramMatrix

-- chart equations
tChart = matrix{tchartParams}*transpose (matrix{centerVars}|matrix{{1_CC}})
qCharts = rfold for i from 0 to m-2 list (
	matrix{take(qchartParams,{5*i,5*i+4})}*
    transpose (matrix{take(qVars,{4*i,4*i+3})}|matrix{{1_CC}})
    )
charts = tChart||qCharts

-- rotation and projection matrices
R = {gateMatrix(id_(FF^3))} | for i from 2 to m list Q2R(w_i,x_i,y_i,z_i)
T = {gateMatrix{{0},{0},{0}}} | for i from 2 to m list (w_i^2+x_i^2+y_i^2+z_i^2)*transpose matrix{{t_(i,1),t_(i,2),t_(i,3)}}
pMatrices = apply(R,T,(r,t)->r|t)

-- camera evaluator used for fabrication routine
C = pMatrices/(P->makeEvaluator(P, varMatrix))

--row vectors giving implicit equations of 2d visibleLines
-- better variable name?
visibleLines = for i from 1 to m list for j from 0 to nLines-1 list (
    if member(j,indepLines) then matrix{{l_(i,j,1),l_(i,j,2),l_(i,j,3)}} 
    else (
	whichIncs = select(last D, i -> member(j,i));
	assert(#whichIncs == 1);
	(k1,k2) = (whichIncs#0#0, whichIncs#0#1);
	a_(i,j,1)*matrix{{l_(i,k1,1),l_(i,k1,2),l_(i,k1,3)}} + a_(i,j,2)*matrix{{l_(i,k2,1),l_(i,k2,2),l_(i,k2,3)}}
	)
    )

-- backprojected planes of visiblelines for each view
planes = for i from 0 to m-1 list apply(visibleLines#i,l->l*pMatrices#i)

-- matrices whose rank deficicency encodes a common point on world lines
CoPmatrices = (last D)/(i->
    rfold(i/(l->
	rfold(planes/(p->
	    p#l
	    ))
	))
    )

-- matrices whose rank deficicency encodes a common world line
-- not needed for m==2
CoLmatrices = if (m<3) then {} else apply(D#0,l->
    rfold(planes/(p->
	    p#l
	    ))
    )

-- master system w/ all max minors from CoL, CoP matrices
elapsedTime F=transpose gateMatrix{
    flatten(
	CoLmatrices/(M -> allMinors(M, 3, Laplace=>COFACTOR)) |
	CoPmatrices/(M -> maxMinors(M, Laplace=>COFACTOR))
	)
    } || charts;
<< " number of polynomials is " << numrows F << endl
-- filter path jumps during monodromy
filterEval = (p,x) -> (
    -- false iff residual small
    resid := norm evaluate(F,x||p);
--    << "residual: " << resid << endl;
    (resid > 4e-4)
    )

-- matrix evaluators: used for rank filter
PE=apply(#CoPmatrices,i->(
	G := CoPmatrices#i;
	(makeEvaluator(G, inputMatrix),
	mutableMatrix(FF,numrows G,numcols G)
	)
    )
)
LE=apply(#CoLmatrices,i->(
	G := CoLmatrices#i;
	(makeEvaluator(G, inputMatrix),
	mutableMatrix(FF,numrows G,numcols G)
	)
    )
)


filterRank = (p,x) -> (
    -- false iff residual small
    (cop,col) := ranks(x,p);
    not(all(cop,x->x==3) and all(col,x->x==2))
    )

filterRankCoP = (p,x) -> (
    -- false iff residual small
    a := PE/( m -> (
	    evaluate(first m, mutableMatrix(x||p), last m);
	    numericalRank matrix last m
	    )
	    );
    not all(a,x->x==3)
    )



--setRandomSeed 31452345342
(p, x) = fabricatepx CC
norm evaluate(F,x||p) -- ~0?
if (instance(Jpivots, Symbol) and JACOBIAN) then (
    -- better to have this precomputed
    << "differentiating" << endl;
    elapsedTime J = diff(varMatrix,F);
    elapsedTime J0 = matrix evaluate(J,x||p);
    elapsedTime Jpivots = rowSelector(J0,Threshold=>5e-5);
    elapsedTime S = first SVD J0^Jpivots;
    )

elapsedTime F'= F^Jpivots;
<< " preparing homotopy " << endl;
elapsedTime PH = parametricSegmentHomotopy(F', cameraVars, dataParams);


if RERUNMONODROMY then elapsedTime (V,np)= monodromySolve(PH, 
    point p, {point x},Verbose=>true,
    FilterCondition=>filterRank,
    Randomizer=>gammify);
if (not instance(FILENAME,Symbol)) then writeStartSys(V.BasePoint, points V.PartialSols, Filename => FILENAME);
stdio << #(points V.PartialSols) << " solutions found!" << endl;

-- clear symbols for next run
w=symbol w
x=symbol x
y=symbol y
z=symbol z
t=symbol t
l=symbol l
a=symbol a
ct=symbol ct
cq=symbol cq
end--

restart
--m=6
--D = (3,1,{{2,3}})

m=5
--D = (3,1,{{1,2},{2,3}}) -- deg = 11306 ?
--D = (4,0,{{0,1,2},{0,3}}) -- deg = 26240
--D = (4,1,{{0,1,2,3},{0,4}}) -- deg = 11008 -- 24077.1 seconds elapsed (AL desktop) 

-- m = 4
-- D = (4,0,{{0,1,3},{1,2},{0,2}}) -- degree 1728??
-- D = (5,0,{{0,1,2,3}}) -- not minimal??
--D = (4,0,{{0,1}}) -- degree = 4525?
--D = (3,2,{{3,4}}) -- degree = 3067??
--D=(2,3,{{1,2},{1,3},{1,4}}) --degree=32??
--D=(3,1,{{0,1},{0,2},{0,3}}) -- degree = 544??
--D =(3,2,{{0,1,2},{0,3},{0,4}}) -- degree = 544??

--m=3
--D = (4,0,{{0,1},{0,2},{1,2}}) -- 216, CLEVELAND, monodromy works
--D = (6,0,{{0,1,2,3,4},{0,5}}) -- 240 -- D.C.-ish, monodromy works
--D = (4,1,{{0,1},{0,4}}) -- 264 -- ANGUILLA, monodromy works
--D = (5,0,{{0,1,3},{0,2,4},{1,2}}) -- 312 -- CHICAGO, monodromy works
--D = (5,1,{{0,1,2,3},{0,5}}) -- 328 -- LOS ANGELES, monodromy works
--D = (4,2,{{4,5}}) -- 360 -- SEATTLE-ish, monodromy works
--D = (5,0,{{0,1,2},{0,3}}) -- 432 -- ODESSA1, monodromy works
--D = (6,0,{{2,3,4,5}}) -- ODESSA2, monodromy works
--D = (6,0,{{0,1,3,4},{0,2,5}}) -- 480 -- PHOENIX, monodromy works
--D = (5,0,{{3,4}}) -- 552 -- DES MOINES, monodromy works
--D = (6,1,{{0,1,2,3,4,5},{0,6}}) -- monodromy gets 64 -- saturation?

-- D = (5,0,{{0,1,3,4},{1,2},{0,2}}) -- 224, monodromy works
-- cases with dependent points
-- D = (4,0,{{0,1},{0,2},{0,3}}) -- 144, monodromy works
-- D=(2,3,{{0,1},{1,2},{1,3},{1,4}}) -- fails rank check
-- D=(3,2,{{1,2},{1,3},{1,4}}) -- 40, monodromy works
-- D=(4,0,{{0,1},{0,2},{0,3}}) --144, monodromy works
-- D=(4,1,{{0,1,2},{0,3},{0,4}}) -- 144, monodromy works
-- D=(4,2,{{0,1,2,3},{0,4},{0,5}})-- 144, monodromy works
--D=(1,5,{{0,1},{0,2},{0,3},{4,5}})--64, monodromy works

--m = 2
--D = (5,0,{{0,1},{1,2},{2,3},{3,4},{4,0}}) -- 20, monodromy works
--D = (4,1,{{0,1},{1,2},{2,3},{3,0},{0,4}}) -- monodromy gets 16 -- saturation?
--D=(2,4,{{0,1},{0,2},{0,3},{1,4},{1,5}}) -- monodromy gets 12 -- saturation?
--D=(3,2,{{0,1},{0,2},{2,0},{0,3},{0,4}}) -- fails rank check
--D=(1,5,{{0,1},{0,2},{0,3},{0,4},{0,5}}) -- fails rank check

COFACTOR = true
JACOBIAN = true
needs "degree.m2"


Jpivots
#Jpivots, nvars
netList apply(CoLmatrices,n->first SVD evaluate(n,x||p))
netList apply(CoPmatrices,n->first SVD evaluate(n,x||p))

-*
--domninance check
Jdom = diff(paramMatrix, F);
elapsedTime Jdom0 = evaluate(Jdom, x||p);
SS=first SVD Jdom0
(max SS)/(min SS)
*-

--elapsedTime filterEval(gammify p,x) -- false?
elapsedTime filterRankCoP(gammify p,x) -- false?
(max S)/(min S) -- small-ish?

elapsedTime (V,np)= monodromySolve(PH, 
    point p, {point x},Verbose=>true,
    FilterCondition=>filterRank,
    Randomizer=>gammify)

writeStartSys(V.BasePoint, points V.PartialSols, Filename => "./numerical/degree-4cameras/3001_1.txt")
-- ,NumberOfNodes=>4,NumberOfRepeats=>15) -- insanity check. for small examples!


G = V.Graph
keys G

-- robustness check
setDefault(tStepMin=>1e-8)
setDefault(tStepMin=>maxCorrSteps=>2)
P0=(transpose matrix V.BasePoint)||random(CC^(#dataParams),CC^1)
Pspec = specialize(PH,P0)
elapsedTime targetSols=trackHomotopy(Pspec,points V.PartialSols);
#(solutionsWithMultiplicity select(targetSols,x->status x == Regular))

-- repeat for other node?