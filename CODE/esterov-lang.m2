-- esterov-lang
restart
needsPackage "MonodromySolver"
X = gateMatrix{toList vars(x,y)}
P = gateMatrix{toList vars(a..k)}
u = inputGate 1
f1 = a + b *x*y + c *x*u/y + d*y*u/x + e * u/(x*y)
f2 = f + g *x*y + h *x*u/y + i*y*u/x + j * u/(x*y)
GS = gateSystem(P,X,transpose gateMatrix{{f1,f2}})
V = first monodromySolve(GS, NumberOfNodes=>3)
length V.PartialSols
writePermutations(V, 8, "esterov-lang.txt")

A = matrix{{0,0},{1,1},{1,-1},{-1,-1},{-1,1}}
smithNormalForm transpose A

-- related experiment: support set with lattice index = 4
restart
needsPackage "MonodromySolver"
X = gateMatrix{toList vars(x,y)}
P = gateMatrix{toList vars(a..k)}
u = inputGate 1
f1 = a + b *x^2 + c *u/y^2 + d*u/x^2 + e * y^2
f2 = f + g *x^2 + h *u/y^2 + i*u/x^2 + j * y^2
GS = gateSystem(P,X,transpose gateMatrix{{f1,f2}})
V = first monodromySolve(GS, NumberOfNodes=>3)
length V.PartialSols
needs "common.m2"
writePermutations(V, 16, "16.txt")



-- related experiment: support set with lattice index = 4
restart
needsPackage "MonodromySolver"
X = gateMatrix{toList vars(x,y)}
P = gateMatrix{toList vars(a..k)}
u = inputGate 1
k=4
f1 = a + b *x^k + c *u/y^k + d*u/x^k + e * y^k
f2 = f + g *x^k + h *u/y^k + i*u/x^k + j * y^k
GS = gateSystem(P,X,transpose gateMatrix{{f1,f2}})
V = first monodromySolve(GS, NumberOfNodes=>k+1)
assert(4*k^2 == length V.PartialSols)
needs "common.m2"
writePermutations(V, 4*k^2, "nonreduced " | toString(4*k^2) | ".txt")
