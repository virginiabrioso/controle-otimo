clear; xdel(winsid()); clc;
        
      //frquência de interesse pacientes 1, 2 e 3
      f1=1.15; f2 = 0.46; f3=0.88;
      
      //frequencias de amostragem mínimas (rad/s)
      ws1=2*f1; //2.3
      ws2 = 2*f2;  // 0.92
      ws3=2*f3;   //1.76
      //frequencias de amostragem mínimas (Hz)
      f1min = ws1/(2*%pi);  //0.3660564
      f2min = ws2/(2*%pi); //0.1464225
      f3min = ws3/(2*%pi); //0.2801127
      
      //maiores períodos de amostragem
      Ts1max=1/f1min; //2.7318197
      Ts2max=1/f2min; //6.8295492
      Ts3max=1/f3min; //3.5699917
      
      //valores utilizados no projeto
      Ts = 1;
      tfinal = 120; // tempo total de simulação em minutos
      N = round( tfinal/Ts ); // total de amostras
      t = 0:Ts:N*Ts-Ts; // vetor de tempo para os plots
      
      // Simulacao de um processo Gaussiano
      noisegen(Ts, tfinal, sqrt(1e-5)); // Matlab:  wgn
      xi = 1*feval(t, Noise); // ruido de processo 
