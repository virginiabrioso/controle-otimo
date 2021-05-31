clear all
close all
clc

%parametros iniciais
K=0.9; wn=2; qsi=0.3;
%tempo de amostragem
tal = 1/(qsi*wn); xmin=tal/10; xmax=tal/5; ts=xmin; %0.1667

%função de transferencia
num = K*wn^2; den =[1 2*qsi*wn wn^2];
Gs=tf(num,den); [z,p,k] = tf2zp(num,den);
%função de transferencia em Z
Gz = c2d(Gs,ts);

%retirando os coeficientes
[numz denz] = tfdata(Gz, 'v'); b0 = numz(1,2); b1 = numz(1,3); a1 = denz(1,2); a2 = denz(1,3);

%parametros
kp = 1; ki = 0.5; kd = 0.3;

%discretizacao 
s=tf('s'); cs = kp + ki/s + kd*s;
z=tf('z'); cz = kp + ki/((1-z^(-1))/(ts)) + kd*((1-z^(-1))/(ts)); [numS denS] = tfdata(cz, 'v');
s0 = numS(1,1)/ts; s1 = numS(1,2)/ts; s2 = numS(1,3)/ts; %s0=(kp+ki*ts+kd/ts); s1=-(kp+2*kd/ts); s2=kd/ts;
t0=ki*ts;
%Parametros de simulacao
tfinal = 20;  N = round( tfinal/ts ); t = 0:ts:N*ts-ts;
%Degrau Unitário
yd(1:10) = 0; yd(11:N) = 1;

%Perturbações no sistema
pi(1:N) = 0.1;
po(1:round(N/2)) = 0; po(round(N/2)+1 : N) = -0.1;
    
%servo
ys(1:2)=0; us(1:2)=0; es(1:2)=0; ISEs = 0; ISUs = 0;

%regulatorio
yr(1:2)=0; ur(1:2)=0; urp(1:2)=0; er(1:2)=0; ISEr = 0; ISUr = 0;

for k = 3:N
  % Caso Servo
        ys(k) = -a1*ys(k-1) -a2*ys(k-2) +b0*us(k-1) +b1*us(k-2);
        es(k) = yd(k) -ys(k);
        us(k) = us(k-1) + t0*yd(k) -s0*ys(k) -s1*ys(k-1) -s2*ys(k-2);
        ISEs = ISEs +es(k)^2;
        ISUs = ISUs +us(k)^2;
        
   % Caso Regulatório
        yr(k) = -a1*yr(k-1) -a2*yr(k-2) +b0*urp(k-1) +b1*urp(k-2) + po(k) + a1*po(k-1) + a2*po(k-2);
        er(k) = yd(k) -yr(k);
        ur(k) = ur(k-1) + t0*yd(k) -s0*yr(k) -s1*yr(k-1) -s2*yr(k-2);
        urp(k) = ur(k) + pi(k);
        ISEr = ISEr +er(k)^2;
        ISUr = ISUr +urp(k)^2;
end

disp('Servo:');
disp('ISE = '); disp(ISEs);
disp('ISU = '); disp(ISUs);
disp('ISE+ISU = '); disp(ISEs+ISUs);

disp('Regulatorio:');
disp('ISE = '); disp(ISEr);
disp('ISU = '); disp(ISUr);
disp('ISE+ISU = '); disp(ISEr+ISUr);

subplot(3,1,1);
plot(t, ys, 'b', t, yd, 'black')
ylabel('y(t) [Unid.]'); xlabel('Tempo [s]');
legend({'Saída','Referência'}, 'Location','northwestoutside')  % inserir legenda, na ordem
title('Caso servo')

subplot(3,1,2); 
plot(t, yr, 'r', t, yd, 'black')
ylabel('y(t) [Unid.]'); xlabel('Tempo [s]');
legend({'Saída','Referência'}, 'Location','northwestoutside') % inserir legenda, na ordem
title('Caso regulatório')

subplot(3,1,3); 
plot(t, us, 'b', t, ur, 'r')
ylabel('u(t) [Unid.]'); xlabel('Tempo [s]');
legend({'Caso servo','Caso regulatório'},'Location','northwestoutside') 
title('Sinais de controle')