clc
clear all

s=tf('s');
M = 0.5;
m = 0.2;
b = 0.1;
I = 0.006;
g = 9.8;
l = 0.3;
q = (M+m)*(I+m*l^2)-(m*l)^2;

sys = (m*l*s/q)/((s^3 + (b*(I + m*l^2))*s^2/q - ((M + m)*m*g*l)*s/q - b*m*g*l/q));

sys1 = (4.45*s)/((s+5.7)*(s-5.6)*(s+0.143));
sys2 = ((4.45*s)*(s^2+8*s+26.66))/(s*(s+5.7)*(s-5.6)*(s+0.143)*(s+50));
rlocus(sys2);
%bode(sys);

% A,B,C,D = tf2ss(b,a);
A = [0,1,0,0 ; 0,-0.1818,2.6727,0; 0,0,0,1;0,-0.4545,31.1818,0;];
B = [0;1.8182;0;4.5455];
C = [0,0,1,0];
D = 0;

% syms k1;
% syms k2;
% syms k3;
% syms k4;
% K=[k1,k2,k3,k4];
p =[-10-1.5i,-2.5+ 1.5i,-2.5-1.5i,-10+1.5i];
K = place(A,B,p);
%  -1.0000   -1.6567   18.6854    3.4594
% K=[-70.7107 , -37.8345 , 105.5298  , 20.9238];
Ac = (A-B*K);
Bc = B;
Cc = C;
Dc = D;
Co = ctrb(A,B);
rank(Co)
Delta_c = eig(Ac);
[b1,a1] = ss2tf(Ac,Bc,Cc,Dc);
% g1=(s*eye(size(Ac,1))-Ac)^(-1);
% g1=(Cc*g1);
% g1=g1*Bc;
GC = tf(b1,a1);
figure(2)
rlocus(GC);
% bode(GC);
figure(3)
GC_s = GC * ((s^2+8*s+26.66));
rlocus(GC_s);

a1=A(1:4:16);
a2=A(2:4:16);
a3=A(3:4:16);
a4=A(4:4:16);

b1=B(1);
b2=B(2);
b3=B(3);
b4=B(4);

c1=C(1);
c2=C(2);
c3=C(3);
c4=C(4);

k1=K(1);
k2=K(2);
k3=K(3);
k4=K(4);
