needs "common.m2"
X=vars(x,y,z)
P=vars(a_1..a_9)
f1 = a_1*(y^2 + z^2) - a_2 *y*z - a_3
f2 = a_4*(z^2 + x^2) - a_5 *x*z - a_6
f3 = a_7*(y^2 + x^2) - a_8 *x*y - a_9

GS = gateSystem(
P, X,transpose gateMatrix{{f1,f2,f3}}
)

V = first monodromySolve(GS, NumberOfNodes=>5, NumberOfEdges=>5)
length V.PartialSols
writePermutations(V,8,"p3p.txt") -- we get G = S_2 wr S_4 ^ A_8


restart
FF = ZZ/101
FFF = frac(FF[a_1..a_9])
R=FFF[x,y,z,MonomialOrder=>Lex]
f1 = a_1*(y^2 + z^2) - a_2 *y*z - a_3
f2 = a_4*(z^2 + x^2) - a_5 *x*z - a_6
f3 = a_7*(y^2 + x^2) - a_8 *x*y - a_9
G = gb ideal(f1,f2,f3)
netList flatten entries gens G


restart
FF = QQ
R=FF[x,y,z,a..f,MonomialOrder=>Eliminate 3]
f1 = (y^2 + z^2) - a *y*z - b
f2 = (z^2 + x^2) - c *x*z - d
f3 = (y^2 + x^2) - e *x*y - f
I=ideal(f1,f2,f3)
Icrit = I+ideal det (jacobian I)^{0..2}
elapsedTime G = groebnerBasis(Icrit, Strategy => "F4");
D1 = (selectInSubring(1, G))_(0,0);
D1

-- extra chec
elapsedTime J = saturate(Icrit, ideal(a..f));
elapsedTime G2 = groebnerBasis(J, Strategy => "F4");
D2 = (selectInSubring(1, G2))_(0,0)

D1-D2
# factor D1
D1-D2

netList apply(flatten entries first coefficients D1, c -> degree c)

needsPackage "SumsOfSquares"

solveSOS D1



