
clear; xdel(winsid()); clc;

Ts = 0.1; 
tfinal = 25;
N=round( (tfinal+Ts)/Ts ); t = 0:Ts:N*Ts-Ts;

// Discrete-time state space model
Ad = [0.9359 0.2428 0.0132 -.0076;
     -.2136 0.8093 0.0440 -.0252;
     0.0018 -.0015 0.971 0.1965; 
     0.0060 -.0049 -.0960 0.6551];
Bd = [.0628 .0002;
     .2093 .0007; 
     .0020 .0274; 
     .0067 .0913];
Cd = [1 0 0 0
     0 0 1 0];
Dd = [0 0;
     0 0;]; 

// matrizes aumentadas
Aa = [eye(2,2) Cd*Ad; 
      zeros(4,2) Ad];
Ba = [Cd*Bd;
        Bd ];       
Ca = [1 0 0 0 0 0
      0 1 0 0 0 0]; 
Da = [0 0]; 

// matrizes de peso
R = 10*eye(2,2);
Q = diag([1 10 10 100 10 1]);

// Using the recursive Riccati Difference Equation (RDE)
//P é a matriz de solução única da equação algébrica de Riccati
P = 1000*eye(6,6);
for i=1:1000
    P = Aa'*P*Aa -Aa'*P*Ba*inv(Ba'*P*Ba+R)*Ba'*P*Aa+Q;
end

//ganho do controlador
K = ( Aa'*P*Ba*inv(Ba'*P*Ba+R) )';

disp("P obtained by the RDE: "); disp(P);
disp('LQR solution for K:'); disp(K);

// Kalman filter (KF)

// matrizes de peso
Qkf = diag([1 10 10 10 10 1]); 
Rkf=100*eye(2,2);

S = 1000*eye(6,6);
for i=1:1000
    S = Aa*S*Aa' -Aa*S*Ca'*inv(Ca*S*Ca'+Rkf)*Ca*S*Aa'+Qkf;
end
disp("S obtained by the RDE: "); disp(S);

L = Aa*S*Ca'*inv(Ca*S*Ca'+Rkf);
disp('Kalman filter solution for L:'); disp(L);

// pertubações
var_w1 = 0.01; var_w2 = 0.001; 
w1 = grand(N, "mn", 0, var_w1);
w2 = grand(N, "mn", 0, var_w2);
var_v = 0.1;
v = grand(N, "mn", 0, var_v);

// initial conditions
r1(1:10)=0; r1(11:N)=1; 
r2(1:10)=0; r2(11:N)=1;   
x(:,:,1) = [0;0;0;0]+[w1(1);w2(1); w1(1); w2(1)];
xest(:,:,1)=[0;0;0;0];  
y(:,1)= Cd*x(:,:,1)+ v(1);
xa(:,:,1)=[0;0;0;0;0;0];
ya(:,1) = Ca*xa(:,:,1);
u = zeros(2,1);
du = zeros(2,1);

//Controlador LQG modelo

for k=2:N
    // Plant simulation using the state-space model
    x(:,:,k) = Ad*x(:,:,k-1) +Bd*u(:,k-1)+[w1(k-1);w2(k-1);w1(k-1);w2(k-1)];
    y(:,k) = Cd*x(:,:,k)+1*v(k); 
    // Kalman filter estimator
    xa(:,:,k) = Aa*xa(:,:,k-1) +Ba*du(:,k-1) +L*( y(:,k-1)-ya(:,k-1) );
    ya(:,k) = Ca*xa(:,:,k);
    // Regulator control law
    du(:,k) = K(1)*r1(k) + K(2)*r2(k)-K*xa(:,:,k);
    u(:,k) = u(:,k-1) +du(:,k);
end

    subplot(211);
    plot(t,r1,':k',t,y(1,:),'b',t,ya(1,:), t, ya(2,:),'r'); title('Discrete-time case');
      legend('Reference','y(t), sensor','y_a(t), estimated',4);
      ylabel('Output (V)'); xlabel('Time (s)');
   subplot(212);
    plot(t,u,'b');
      ylabel('Control (V)'); xlabel('Time (s)');
