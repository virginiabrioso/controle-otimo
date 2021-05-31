// Exemplo: regulador de variância mínima
clear; xdel(winsid()); clc;

// Condições Iniciais 
z = %z;
Ts = 0.1;
tfinal = 50; 
N = round( tfinal/Ts );
t = 0:Ts:N*Ts-Ts; 
yr(1:50)=0; yr(51:N+1)=1;
y(1:3)=0; u(1:3)=0; e(1:3)=0;

// Parametros fornecidos
Bz = [0.0234015 0.5489666 0.0751607]
b0 = Bz(1); b1 = Bz(2); b2 = Bz(3);
Az = [1 -0.4641429 0.077525 0.0429111];
a1 = Az(2); a2 = Az(3); a3 = Az(4);

// Ruido branco
Variancia = 0.0007362;
xi = grand(N, "mn", 0, Variancia); 

// Projeto do GMV Posicional
f0 = -a1;
f1 = -a2;
f2 = -a3;
q0 = 1; 

for k = 4:N
    
    // Simula modelo do processo
    y(k) = -a1*y(k-1) -a2*y(k-2) -a3*y(k-3) +b0*u(k-1) +b1*u(k-2) +b2*u(k-3)+xi(k);
    
    // Parcela do regulador de variância mínima
    u(k) = (1/(b0+q0))*( -b1*u(k-1) -b2*u(k-2) +yr(k+1) -f0*y(k) -f1*y(k-1) -f2*y(k-2));
    e(k) = yr(k) - y(k);
end

subplot(3,1,1);
plot(t,yr(1:N),'k')
ylabel('yr(t) [Unid.]'); xlabel('Tempo [s]');
legend({'Referência'}, 'Location','northwestoutside') 

subplot(3,1,2); 
plot(t,u,'b')
ylabel('u(t) [Unid.]'); xlabel('Tempo [s]');
legend({'Controle'},'Location','northwestoutside') 

subplot(3,1,3); 
plot(t, y, 'r')
ylabel('y(t) [Unid.]'); xlabel('Tempo [s]');
legend({'Saída'}, 'Location','northwestoutside')



