-*
--symbolic: degrees of abs pose problems
restart
(npoints, nlines) = (2,1)
FF=QQ
R=FF[r_(1,1)..r_(3,3),t_1..t_3]
Rot=genericMatrix(R,3,3)
Cam=Rot|matrix{{t_1},{t_2},{t_3}}
pts3d = for i from 1 to npoints list random(FF^4,FF^1)
pts2d = for i from 1 to npoints list random(FF^3,FF^1)
lns3d = for i from 1 to nlines list random(FF^2,FF^4)
lns2d = for i from 1 to nlines list random(FF^1,FF^3)
I=ideal(Rot*transpose Rot-id_(R^3), det Rot -1) + 
    ideal apply(pts3d,pts2d,(p,q) -> minors(2,Cam*p|q)) + 
    ideal apply(lns3d,lns2d,(l,m) -> minors(3,l||m*Cam))
dim I
degree I
*-

--numerical
needs "common.m2"
--needs "galois.m2"
--monodromy options
msOptions={Verbose=>true,NumberOfNodes=>3}

-- setting up absolute pose equations
Rot = matrix for i from 1 to 3 list for j from 1 to 3 list declareVariable r_(i,j)
Cam=Rot|transpose matrix{for i from 1 to 3 list declareVariable t_i}
pts3d = for i from 1 to npoints list transpose matrix{for j from 1 to 4 list declareVariable P_(i,j)}
pts2d = for i from 1 to npoints list transpose matrix{for j from 1 to 3 list declareVariable Q_(i,j)}
lns3d = for i from 1 to nlines list matrix{for j from 1 to 4 list declareVariable L1_(i,j)} || matrix{for j from 1 to 4 list declareVariable L2_(i,j)}
lns2d = for i from 1 to nlines list matrix{for j from 1 to 3 list declareVariable M_(i,j)}
FF=CC
ptParam = transpose(rfold pts3d || rfold pts2d)
lnParam = if instance(lns3d, Symbol) then random(FF^0, FF^0) else cfold apply(lns3d,l->matrix{flatten entries l}) | cfold lns2d
Param = if npoints == 0 then lnParam else if nlines == 0 then ptParam else ptParam | lnParam
X = gateMatrix{flatten entries Cam}
eqs=flatten entries(Rot*transpose Rot - id_(CC^3))|
   flatten apply(pts3d,pts2d,(p,q) -> allMinors(Cam*p|q,2,Laplace=>true)) |
   flatten apply(lns3d,lns2d,(l,m) -> allMinors(l||m*Cam,3,Laplace=>true))
G=gateSystem(Param,X,transpose gateMatrix{eqs})

--sampling the incidence variety
R0=(-1)*randomOn 3
t0 = random(CC^3,CC^1)
Cam0 = R0|t0
pts3d0 = for i from 1 to npoints list random(CC^4,CC^1)
pts2d0 = for i from 1 to npoints list Cam0* pts3d0#(i-1)
lns3d0 = for i from 1 to nlines list random(CC^2,CC^4)
lns2d0 = apply(lns3d0, l -> (
        p1p2 := numericalKernel(l,1e-5);
        transpose numericalKernel(transpose(Cam0*p1p2),1e-5)
        )
    )
lnParam0 = if instance(lns3d0, Symbol) then random(FF^0, FF^0) else cfold apply(lns3d0,l->matrix{flatten entries l}) | cfold lns2d0
ptParam0=transpose(rfold pts3d0 || rfold pts2d0)
Param0 = point if npoints == 0 then lnParam0 else if nlines == 0 then ptParam0 else ptParam0 | lnParam0
X0= point matrix{flatten entries Cam0}
evalX0=evaluate(G,Param0,X0)
resid=norm evalX0
assert areEqual(0,resid)
Gsq = squareDown(Param0,X0,numVariables G,G)
end--
restart
needsPackage "MonodromySolver"
for i from 0 to 3 do (
    npoints = i;
    nlines = 3-npoints;
    load "abs-pose.m2";
    monodromyGroup(Gsq, Param0, {X0}, "msOptions"=>msOptions, FileName=>"npts-"|toString(npoints)|"nlns"|toString(nlines));
    )

-- let's have a closer look at the one w/ 4 solutions
npoints = 2;
nlines = 3-npoints;
setRandomSeed 0
load "abs-pose.m2";
netList entries gateMatrix G
V=first monodromySolve(Gsq, Param0, {X0})
sols = points V.PartialSols
Rots = apply(coordinates \ sols, x -> matrix{{x#0,x#1,x#2},{x#4,x#5,x#6},{x#8,x#9,x#10}})
translations = apply(coordinates \ sols, x -> matrix{{x#3},{x#7},{x#11}})
v = transpose Rots#0 * translations#0 - transpose Rots#1 * translations#1
SVD(v|A_{0})
eigenvalues(A-id_(CC^3))
(eval,evec)=eigenvectors(A-id_(CC^3))
SVD(evec_{0}|v)
A=transpose Rots#0 *Rots#1 + id_(CC^3)

-*
M=x y^t = [y1 x | y2 x | y3 x]
      = 
      
M^t M = (x^t x) y y^t -> eig-vals = 1, 0, 0
                      -> eigenvector for e-val 1: a y for some scalar a
                          M^t M * a y = a (x^t x) y = 1 * y
                                      => a = 1/(x^t x)
I diag(*, 0, 0) [x ]^t

*-

first SVD A -- rots 0 and rots 2 differ by a reflection
-- rk A = 1
-- R0^t R1 + I = A(R0,t0) <-> R1 = R0 * ( A(R0, t0) - I) = R0 (2 R0 * R0^t
-- d: (R, t) -> R 
-*
lets sample many problem instances and try to interpolate the deck transformation
we just need to see how the rank-1 matrix A depends on the input data: m parameters and n unknowns, w/ (m, n) = (25, 12)
*-
nsamples = (m, n, d) -> 1+2*binomial(d+m+n,m+n)
for d from 0 to 5 list nsamples(25,12,d)--with parameters
for d from 0 to 5 list nsamples(0,12,d)--without parameters
--a=p(x)/q(x), p(x) = sum(ci x^i) q(x)=sum(-di x^i) => ci x^i + a dix^i = (ci+a) x^i => M(x) C = a(x) N(x) D => C|D in ker(M(x) ++ a(x) N(x))
p0=transpose matrix V.BasePoint
x0=sols#0
x1=sols#1
RNG = CC[rr_(1,1)..rr_(3,3),tt_1..tt_3,vv]

dnum=2
Bnum = matrix{for i from 0 to dnum list basis(i,RNG,Variables=>drop(gens RNG,-1))};
inputGate v
ddenom=2
Bdenom = vv*matrix{for i from 0 to ddenom list basis(i,RNG,Variables=>drop(gens RNG,-1))};
B=Bnum|Bdenom;
GInterp=gateSystem(
    gateMatrix{getVarGates RNG}
    , transpose gateMatrix{gatePolynomial \ (flatten entries(B))})
elapsedTime M = matrix for i from 0 to numcols B list (
    p1 := random(CC^25,CC^1);
    H := specialize(parametricSegmentHomotopy Gsq, p0||p1);
    sols := trackHomotopy(H,{x0,x1});
    Rots := apply(coordinates \ sols, x -> matrix{{x#0,x#1,x#2},{x#4,x#5,x#6},{x#8,x#9,x#10}});
    translations := apply(coordinates \ sols, x -> matrix{{x#3},{x#7},{x#11}});
    A := transpose Rots#0 *Rots#1 + id_(CC^3);
    {evaluate(GInterp, matrix first sols | matrix{{A_(0,0)}})}
    );
elapsedTime K = numericalKernel(M,1e-5);
NRMS = sort for i from 0 to numcols K -1 list (norm K_{i}^{0..numcols Bnum-1}, norm K_{i}^{numcols Bnum..numcols Bnum + numcols Bdenom-1});
last NRMS -- = max
netList NRMS
