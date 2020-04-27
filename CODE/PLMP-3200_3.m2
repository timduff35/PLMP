restart
Rng = QQ[r_(1,1)..r_(3,3),t_1..t_3,n_1..n_3,d]
R = transpose(genericMatrix(Rng,r_(1,1),3,3))
T = genericMatrix(Rng,t_1,3,1)
N = genericMatrix(Rng,n_1,3,1)
A = random(QQ^3,QQ^3)

-- formulation with ||N||=1
I = ideal(A - d*R-T*transpose(N)) +
    ideal(transpose(R)*R-id_(QQ^3)) + 
    ideal(det(R)-1) +ideal(transpose(N)*N-1)
dim I, degree I

-- formulation with N_3 =1
I = ideal(A - d*R-T*transpose(N)) +ideal(transpose(R)*R-id_(QQ^3)) +ideal(det(R)-1) +ideal(N_(2,0)-1)
dim I, degree I
