-- EXAMPLE of computations for problem 2111_1 (minimal, degree=40) 
restart
setRandomSeed 0
FF = ZZ/10007
R = FF[r_(1,1), r_(1,2), r_(1,3), r_(2,1), r_(2,2), r_(2,3), t_(1,1), t_(2,1),t_(1,2), t_(2,2), t_(2,3)]
Rs = {matrix{{1, 0, 0}, {0, 1, 0}, {0, 0, 1}}, 
matrix{{-r_(1,1)^2-r_(1,2)^2+r_(1,3)^2+1,-2*r_(1,2)*r_(1,3)+2*r_(1,1),2*r_(1,1)*r_(1,3)+2*r_(1,2)},
{-2*r_(1,2)*r_(1,3)-2*r_(1,1),-r_(1,1)^2+r_(1,2)^2-r_(1,3)^2+1, -2*r_(1,1)*r_(1,2)+2*r_(1,3)},
{2*r_(1,1)*r_(1,3)-2*r_(1,2), -2*r_(1,1)*r_(1,2)-2*r_(1,3),r_(1,1)^2-r_(1,2)^2-r_(1,3)^2+1}},
matrix{{-r_(2,1)^2-r_(2,2)^2+r_(2,3)^2+1, -2*r_(2,2)*r_(2,3)+2*r_(2,1),
      	2*r_(2,1)*r_(2,3)+2*r_(2,2)}, {-2*r_(2,2)*r_(2,3)-2*r_(2,1),
      -r_(2,1)^2+r_(2,2)^2-r_(2,3)^2+1, -2*r_(2,1)*r_(2,2)+2*r_(2,3)},
      {2*r_(2,1)*r_(2,3)-2*r_(2,2), -2*r_(2,1)*r_(2,2)-2*r_(2,3),
      r_(2,1)^2-r_(2,2)^2-r_(2,3)^2+1}}};
netList Rs
ts={matrix{{0},{0},{0}},matrix{{t_(1,1)},{t_(1,2)},{1}},matrix{{t_(2,1)},{t_(2,2)},{t_(2,3)}}};
C = apply(Rs,ts,(R,t)->R|t);

P = {matrix{{-2639, -4936, 1789}, {2653, -591, -643}, {1, 1, 1}}, matrix{{-3868, -1776, 3174}, {3669, -4143, -1982}, {1, 1, 1}}, matrix{{-1889,
     -1604, 4629}, {-473, -4513, -4210}, {1, 1, 1}}}
Q = {matrix{{-3777}, {-974}, {-4900}}, matrix{{-609}, {4098}, {-4609}},matrix {{-4458}, {-1449}, {2627}}}

L = {matrix{{107, 4376, 3187}, {3897, -4638, 2998}, {-1009, 2785, -4328},
      {2883, -1540, 1031}, {-1517, -713, 3879}}, matrix{{3783, -1437, 275},
      {2926, -687, -1310}, {-2582, -2466, 1236}, {-2721, -1085, -1127}, {3072,
      4259, 1727}}, matrix{{-1563, -4868, -1852}, {-968, -960, -1036}, {-3676,
      382, 1454}, {1795, 2177, -4909}, {2696, -1409, 1206}}}

cl = (C,L') -> minors(3,matrix apply(C,L',(c,l)->{l*c}),Strategy=>Cofactor);
cp = (C,L') -> minors(4,matrix apply(C,L',(c,l)->{l*c}),Strategy=>Cofactor);


D = (3,2,{{1,2},{1,3},{1,4}}) -- encodes problem 2111_1 
-* line 0 is free
   line 1 passes through points 0,1,2
   line 2 passes through point 0
   line 3 passes through point 1 (ghost line)
   line 4 passes through point 2 (ghost line)
*-    
nLines = D#0 + D#1
lineIncidences = D#2

Icl = sum(nLines,i->cl(C,L/(l->l^{i}))); -- LC constraint for each (non-ghost) line
Icp = sum(lineIncidences,LINES->cp(C,L/(l->l^LINES))); -- CP constraint for each point
I=Icl+Icp;

--inequations
Rz = FF[z, gens R, MonomialOrder=>{Eliminate 1}]
Jz = sub(I,Rz) + ideal(1+z*(r_(2,1)^2+r_(2,2)^2+r_(2,3)^2+1)*(r_(1,1)^2+r_(1,2)^2+r_(1,3)^2+1));
G = groebnerBasis(Jz, Strategy => "F4"); -- exploit the F4 algorithm 

-- now eliminate z to study the original ideal I
gbI = selectInSubring(1,G); -- Groebner basis of I, which equals Jz intersected with  FF[r_(i,j),t_i]
inI = ideal sub(leadTerm gbI,R) -- initial ideal of I
dim inI 
degree inI
