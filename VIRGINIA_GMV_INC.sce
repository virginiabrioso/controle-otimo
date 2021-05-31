// Patient models
    s=%s;
     // Patient 1
    theta = 3; // time delay for the drug to "kick in"
    num = (s+0.055)*(s+0.0033)*1.15;
    den = 2.0164*(s+0.0792)*(s+0.114)*(s+0.0419)*(s+1.15);
 
   /*     // Patient 2
    theta = 3; // time delay for the drug to "kick in"
    num = (s+0.055)*(s+0.0033)*0.46;
    den = 2.0164*(s+0.0792)*(s+0.114)*(s+0.0419)*(s+0.46);    

   
  // Patient 3
    theta = 8; // time delay for the drug to "kick in"
    num = (s+0.055)*(s+0.0033)*0.88;
    den = 2.0164*(s+0.0792)*(s+0.114)*(s+0.0419)*(s+0.88);
 */ 
  
      Gs = syslin('c', num, den); // no delay at this time!        
      // Patient 1 ZOH equivalent (without the delay yet)      
      Gz = ss2tf( dscr( tf2ss(Gs), Ts) ); disp(Gz);
       // Extracting polynomial parameters
       Bz = coeff(Gz.num); Bz = flipdim(Bz,2);
        b0=Bz(1); b1=Bz(2); b2=Bz(3); b3=Bz(4);
        
       Az = coeff(Gz.den); Az = flipdim(Az,2);
        a1=Az(2); a2=Az(3); a3=Az(4);
      
        // Adding the discrete time delay d
        d = theta/Ts; // discrete time delay
        z=%z;
        Gdz = Gz*z^(-d); // ZOH equivalent with time delay
           
// Projeto do GMV Incremental
   // Eq. Diofantina do MVR:  1 = Delta*A(z) + (z^-1) * F(z)
   // no caso incremental, nf = na
   Delta = [1 -1]; // Delta = 1-z^-1
   Abarz = conv(Delta,Az); // Delta*Az
      abar1 = Abarz(2); 
      abar2 = Abarz(3); 
      abar3 = Abarz(4); 
      abar4 = Abarz(5); 
      abar5 = Abarz(6);
      
f0 = -abar1;
f1 = -abar2;
f2 = -abar3;
f3 = -abar4;
f4 = -abar5;
q0 = 55.5; // Ponderacao do esforço de controle
        
   // Reference sequence for the Degree of Hypnosis (DOH)
   yr(1:4)=90; yr(5:N+1)=50; // 50% of DOH
   
y(1:4+d)=90; // Note d participating in the 
u(1:4+d)=200; // initial conditions
du(1:4+d)=0; // do controlador: Delta*u(k)
ISE = 0; ISU = 0;

for k = 5+d:N
    // without noise
    y(k) = -a1*y(k-1) -a2*y(k-2) -a3*y(k-3) -a4*y(k-4) ...
      +b0*u(k-1-d) +b1*u(k-2-d) +b2*u(k-3-d) +b3*u(k-4-d);
      
        //with noise
   /* y(k) = -a1*y(k-1) -a2*y(k-2) -a3*y(k-3) -a4*y(k-4) ...
      +b0*u(k-1-d) +b1*u(k-2-d) +b2*u(k-3-d) +b3*u(k-4-d)+xi(k);
      */   
    
    // Parcela do regulador de variância mínima
    du(k) = (1/(b0+q0))*( -b1*du(k-1) -b2*du(k-2) -b3*du(k-3) ...
             +yr(k+1) -f0*y(k) -f1*y(k-1) -f2*y(k-2) -f3*y(k-3) -f4*y(k-4));
             
    u(k) = u(k-1) +du(k);
    e(k)=yr(k)-y(k);
   
    //parametros de desempenho
        ISE = ISE +e(k)^2;
        ISU = ISU +u(k)^2;
end

//parametros

disp(ISE)
disp(ISU)
J= ISE + ISU; disp(J)
disp(variance(e))
disp(variance(u))

// Plots dos resultados
fig1 = scf();
   subplot(211);
      plot(t,yr(1:N),':k');
      plot(t,60,'k--');
      plot(t,40,'k--');
      plot(t,y); ylabel('Saída [Unid.]'); xlabel('Tempo [min]');
    title(['Variância de e(k) = ' string( variance(e) ) ...
      ', Variância de u(k) = ' string( variance(u))]);
      
   subplot(212);
      plot(t,u); ylabel('Controle [Unid.]'); xlabel('Tempo [min]');
