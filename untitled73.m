clc;
A = input('Enter the adjacency matrix of the graph: ');
degree = sum(A);

final_polynomial = 0;

for i = 1:size(A,1)
 for j = i+1:size(A,2)
 if A(i,j) == 1

 degree_i = degree(i);
 degree_j = degree(j);
 
 if degree_i ~= degree_j
 x_power = min(degree_i, degree_j);
 y_power = max(degree_i, degree_j);
 else
 x_power = degree_i;
 y_power = degree_i;
 end

  poly_x = sym('x')^ x_power ;
            poly_y = sym('y')^ y_power;
            polynomial = poly_x * poly_y;
 final_polynomial = final_polynomial + polynomial;
           
 end
 end
end

degree = sum(A);
syms x y k;
final_polynomial = 0;
Sfinal_polynomial = 0;
Rfinal_polynomial = 0;

for i = 1:size(A,1)
 for j = i+1:size(A,2)

 if A(i,j) == 1

 degree_i = degree(i);
 degree_j = degree(j);

  if degree_i ~= degree_j
 x_power = min(degree_i, degree_j);
 y_power = max(degree_i, degree_j);
 else
 x_power = degree_i;
 y_power = degree_i;
  end

  poly_x = sym('x')^ x_power ;
            poly_y = sym('y')^ y_power;
            polynomial = poly_x * poly_y;
 
 final_polynomial = final_polynomial + polynomial;

  D = x_power * y_power;
  Dk = D^k;
  
 poly_x = sym('x')^ x_power ;
            poly_y = sym('y')^ y_power;
            Rpolynomial = Dk * poly_x * poly_y;
            Spolynomial = 1/Dk * poly_x * poly_y;

 Rfinal_polynomial = Rfinal_polynomial + Rpolynomial;
 Sfinal_polynomial = Sfinal_polynomial + Spolynomial;
 end
 end
end

Rk = subs(Rfinal_polynomial, [x,y], [1,1]);
RRk = subs(Sfinal_polynomial, [x,y], [1,1]);

disp('The M-polynomial for the graph is M(x,y) = :');
disp(final_polynomial);

Dy = simplify(y*diff(final_polynomial,y));
D2y = simplify(y*diff(Dy,y));
D3y = simplify(y*diff(D2y,y));
Dx = simplify(x*diff(D3y,x));
D2x = simplify(x*diff(Dx,x));
D3x = simplify(x*diff(D2x,x));
 Apolynomial = subs(D3x, y, x);
 Apolynomial = simplify(x^(-2)*Apolynomial);
 Apolynomial = int(Apolynomial/x,x);
 Apolynomial = int(Apolynomial/x,x);
 Apolynomial = int(Apolynomial/x,x);
A = subs(Apolynomial, x,1);
disp(['The Augmented Zagreb index is A(G) = :' char(expand(Apolynomial))]);
disp([' A(G) = ' char(A)]);

syms x y k;
Dx = simplify(x*diff(final_polynomial,x));
Dy = simplify(y*diff(final_polynomial,y));
M1 = simplify(Dx + Dy);
Z1 = subs(M1, [x,y], [1,1]);
disp(['The first Zagreb index is M1 = :', char(expand(M1))]);
disp([' M1(G) = ' char(Z1)]);

Dy = simplify(y*diff(final_polynomial,y));
DxDy = simplify(x*diff(Dy,x));
Z2 = subs(DxDy, [x,y], [1,1]);
disp(['The second Zagreb index is M2 = : ', char(expand(DxDy))]);
disp([' M2(G) = ' char(Z2)]);

Sy = simplify(int(final_polynomial/y, y));
SxSy = simplify(int(Sy/x, x));
mM2 = subs(SxSy, [x,y], [1,1]);
disp(['The second modified Zagreb index is mM2 = : : ', char(expand(SxSy))]);
disp([' mM2(G) = ' char(mM2)]);


final_polynomial = subs(final_polynomial, y, x);
Q = int(final_polynomial/x, x);
Harmonic = simplify(2*Q);
H = subs(Harmonic, [x,y], [1,1]);
disp(['The Harmonic index is H = : ', char(expand(Harmonic))]);
disp([' H(G) = ' char(H)]);

DyM1 = simplify(y*diff(M1,y));
RezG = simplify(x*diff(DyM1,x));
RezG3 = subs(RezG, [x,y], [1,1]);
disp(['The RezG3 index is: ', char(expand(RezG))]);
disp([' RezG3 = ' char(RezG3)]);

DxSy = simplify(x*diff(Sy,x));
SxDy = simplify(int(Dy/x, x));
SDD1= simplify(DxSy + SxDy);
SDD = subs(SDD1, [x,y], [1,1]);
disp(['The SDD index is: ', char(expand(SDD1))]);
disp([' SDD(G) = ' char(SDD)]);

Q = subs(DxDy, y, x);
IG = simplify(int(Q/x, x));
I = subs(IG, [x,y], [1,1]);
disp(['The Inverse sum indeg index is: ', char(expand(IG))]);
disp([' I(G) = ' char(I)]);

D2x = simplify(x*diff(Dx,x));
D2y = simplify(y*diff(Dy,y));
Forgotten = simplify(D2x + D2y);
F = subs(Forgotten, [x,y], [1,1]);
disp(['The Forgotten index is: ', char(expand(Forgotten))]);
disp([' F(G) = ' char(F)]);

disp('The Randi´ c index is R_k(G) = :');
disp(Rfinal_polynomial);
disp([' R_k(G) = ' char(Rk)]);

disp('The Inverse Randi´ c index is RR_k(G) = :');
disp(Sfinal_polynomial);
disp([' RR_k(G) = ' char(RRk)]);



