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
lnParam = cfold apply(lns3d,l->matrix{flatten entries l}) | cfold lns2d
Param = if npoints == 0 then lnParam else if nlines == 0 then ptParam else ptParam | lnParam
X = gateMatrix{flatten entries Cam}
eqs=flatten entries(Rot*transpose Rot - id_(CC^3))|
   flatten apply(pts3d,pts2d,(p,q) -> allMinors(Cam*p|q,2,Laplace=>true)) |
   flatten apply(lns3d,lns2d,(l,m) -> allMinors(l||m*Cam,3,Laplace=>true))
G=gateSystem(Param,X,transpose gateMatrix{eqs})

--sampling the incidence variety
R0=randomOn 3
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
lnParam0 = cfold apply(lns3d0,l->matrix{flatten entries l}) | cfold lns2d0
ptParam0=transpose(rfold pts3d0 || rfold pts2d0)
Param0 = point if npoints == 0 then lnParam0 else if nlines == 0 then ptParam0 else ptParam0 | lnParam0
X0= point matrix{flatten entries Cam0}
assert areEqual(0,norm evaluate(G,Param0,X0))
Gsq = squareDown(Param0,X0,numVariables G,G)
end--
restart
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
V=first monodromySolve(Gsq, Param0, {X0})
sols = points V.PartialSols
Rots = apply(coordinates \ sols, x -> matrix{{x#0,x#1,x#2},{x#4,x#5,x#6},{x#8,x#9,x#10}})
translations = apply(coordinates \ sols, x -> matrix{{x#3},{x#7},{x#11}})
first SVD(transpose Rots#0 *Rots#1 + id_(CC^3)) -- rots 0 and rots 2 differ by a reflection
(S,U,Vt)=SVD(transpose Rots#0 *Rots#1 + id_(CC^3))
n=transpose (transpose (U*diagonalMatrix(S)))_{0}|matrix{{0}}
clean_(1e-3) n
--other pair?
(S,U,Vt)=SVD(transpose Rots#2 *Rots#3 + id_(CC^3))
n=transpose (transpose (U*diagonalMatrix(S)))_{0}|matrix{{0}}
clean_(1e-3) n
-- deck transformation is given by reflection thru hyperplane perp to (1 0 0 0) in P^3
