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
    Gz = ss2tf( dscr( tf2ss(Gs), Ts) ); disp(Gz);
       // Extracting polynomial parameters
       Bz = coeff(Gz.num); Bz = flipdim(Bz,2);
        b0=Bz(1); b1=Bz(2); b2=Bz(3); b3=Bz(4);
       Az = coeff(Gz.den); Az = flipdim(Az,2);
        a1=Az(2); a2=Az(3); a3=Az(4); a4=Az(5);
      
        // Adding the discrete time delay d
        d = theta/Ts; // discrete time delay
        z=%z;
        Gdz = Gz*z^(-d); // ZOH equivalent with time delay

   // Reference sequence for the Degree of Hypnosis (DOH)
   yr(1:4)=90; yr(5:N)=50; // 50% of DOH
   ISE = 0; ISU = 0;
   
y(1:4+d)=90; // Note d participating in the 
u(1:4+d)=200; // initial conditions
for k=5+d:N  
    // without noise
    y(k) = -a1*y(k-1) -a2*y(k-2) -a3*y(k-3) -a4*y(k-4) ...
      +b0*u(k-1-d) +b1*u(k-2-d) +b2*u(k-3-d) +b3*u(k-4-d);
    
    //with noise
   /* y(k) = -a1*y(k-1) -a2*y(k-2) -a3*y(k-3) -a4*y(k-4) ...
      +b0*u(k-1-d) +b1*u(k-2-d) +b2*u(k-3-d) +b3*u(k-4-d)+xi(k);
      */      
    kp=0.25; ki=0.03; kd=0.6;
    s0 = kp+ki*Ts+kd/Ts; s1=-kp-2*kd/Ts; s2=kd/Ts;
    e(k)=yr(k)-y(k);
    u(k)=u(k-1)+s0*e(k)+s1*e(k-1)+s2*e(k-2);
    
    ISE = ISE +e(k)^2;
    ISU = ISU +u(k)^2;
end

disp(ISE)
disp(ISU)
J= ISE + ISU; disp(J)
disp(variance(e))
disp(variance(u))

// Ploting results
figure1=scf();
 subplot(211);
  plot(t,60,'k--');
  plot(t,40,'k--');
  plot(t,yr,'k');
  plot(t,y,'b');
title(['Variância de e(k) = ' string( variance(e) ) ...
      ', Variância de u(k) = ' string( variance(u))]);
   ylabel('DOH [WAVcns]'); xlabel('Time [min]');
   legend('Setpoint','DOH (%)');
 subplot(212);
  plot(t,u,'b'); 
   ylabel('Infusion rate [mcg/kg/min]')
   xlabel('Time [min]'); 
