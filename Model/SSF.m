% ===============================================================
%  FUNCTION: Split-step Fourier waveguide propagation
%  Release year: 2020
%  Author: Stijn Cuyvers
%  Affiliation: Photonics Research Group, Ghent University - imec, Belgium
% ===============================================================

function Prop=SSF(Signal_in,T_span,lambda0, L,betas,loss,n2,Aeff,MODEL_TPA_and_FCA,MODEL_RAMAN)
% REMARK: normally, trace should have a length that is a power of 2 for the
% FFT, MATLAB 2011 and newer can however handle non power of 2 signal
% lengths
% It is also possible to interpolate the data to match a power of 2 length,
% but care should be taken not to distort phase and or amplitude
%tic

showevolution=0; % DEBUGGING, show pulse evolution

global Nc
global Nc_avg

%global Nc_monitor
%global counter

% -- ARGUMENTS --
% Signal_in: input signal
% T_span: corresponding time grid
% lambda0: wavelength [nm]
% L: length of propagation
% betas: dispersion
% loss: linear losses
% n2: nonlinear Kerr nonlinearity [m^2/W]
% Aeff: effective mode area [m^2]

% --- PARAMETERS ---
c       = 3e8;      % speed of light [m/s]
lambdap = lambda0;  % wavelength [m]
npas    = 100;     % number of steps == NEEDS OPTIMIZATION ==
w0 = (2.0*pi*c)/lambdap;    % ref frequency
btpa = 6e-12; %0.8e-12;             % beta_TPA [m/W] (two-photon absorption parameter)
k0 = 2*pi/(lambdap);        
gamma = k0*n2/Aeff + 1i*btpa/(2*Aeff); % Nonlinearity, TPA

if MODEL_TPA_and_FCA==0
    gamma = k0*n2/Aeff;
end

sigma = 1.45e-21;           % Free carrier absorption [m^2]

kc = 1.35e-27;              % Free carrier dispersion [m^3]
mu = 2*kc*2*pi/lambdap/sigma;

% NEEDS TO BE ADDRESSED IN FUTURE: NO 'RESET' BETWEEN PULSES (RECOVERY TIME
% SLOW COMPARED TO PULSE REP RATE)
tau = 1e-9;                 % Free carrier lifetime [s]
planck = 6.626e-34;
D = 2*pi*btpa/(2*planck*w0*Aeff^2);


original_length=size(Signal_in,2);
npt=original_length;
t = linspace(-T_span/2,T_span/2,original_length);                 % time grid [ps]
dt = t(2)-t(1);                                                             % time grid spacing [ps]
f = (-original_length/2:original_length/2-1)/(original_length*dt);     % angular frequency deviation from w0 [10^12 rad/s]


lambda = c./(c/lambdap+f);

alpha = log(10.^(loss/10));     % attenuation coefficient

h       = L/npas;               % propagation step

lambdadisp = [1200 3300];       % wavelength range for plot
lambdaidx  = find((lambda*1e9>=lambdadisp(1)) & (lambda*1e9<=lambdadisp(2)));





%% ----------------------------------------------------------------------

om = fftshift(2*pi*f.');

% input signal
E = Signal_in;

% Nc: TPA-generated free-carrier density
%Nc = E.*0; %%%%% initialisation Nc;


% %%%%%%%%%%%%% Dispersion %%%%%%%%%%%%%%%%
B = 0;
for i = 1:size(betas,2)        % Taylor expansion of betas
B = B + betas(i)/factorial(i+1).*om.^(i+1);
end
opdisp = exp(1i*h*B);

%%%%%%%%%%%%% Raman response %%%%%%%%%%%%%%%%
wr = 15.6e12*2*pi;
dw = 105e9*pi;
gRmax = 1.5e-10; %m/W
fr = gRmax*dw/(wr*k0*n2); % see Nonlinear optical phenomena in silicon waveguides: modeling and applications
if MODEL_RAMAN==0
    fr = 0;
end

RW = wr^2./(wr^2-om.^2-2*1i.*om*dw);
RT = real(fft(RW));
RT = RT/(trapz(t,RT));
RW = npt*ifft(RT);


idxp = linspace(2,npt+1,npt);
idxp(end) = 1;
idxm = linspace(0,npt-1,npt);
idxm(1)=npt;

for k=1:npas,
   
  %%%%%%%%%%%%% nonlinear step (écrit par Stef (voir blow and wood))  %%%%%%%%%  
    
    U0 = E.';
    TI = ifft(abs(U0).^2);
    int0 = (1-fr)*abs(U0).^2 + dt*fr*fft(RW.*TI);
    Nd = U0.*int0;
    U1 = U0-h/2*gamma/w0*(Nd(idxp)-Nd(idxm))/(2*dt);   % On est au 1/2 pas de RK
    TI = ifft(abs(U1).^2); int1 = (1-fr)*abs(U1).^2+dt*fr*fft(RW.*TI);
    Nd = U1.*int1;
    U1 = U0+ h*1i*gamma*U1.*(int1-int0)-h*gamma/w0*(Nd(idxp)-Nd(idxm))/(2*dt);
    %U1 = U0;
    E  = U1.*exp(1i*gamma*h*int0);
    E = E.';
    
  
  %%%%%%%%%%%%%%%% carrier density %%%%%%%%%%%%%%%%%%%%%%      
    
   %Nc(2:npt)=(Nc(1:npt-1)).*exp(-dt/tau)+(1-exp(-dt/tau)).*tau*D*abs(E(2:npt)).^4;
    
  %%%%%%%%%%%%%%%%%%% losses %%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
  if MODEL_TPA_and_FCA~=0  
    E = E.*exp(h*(-alpha/2 - sigma/2*(1+1i*mu)*Nc_avg));
  else
    E = E.*exp(h*(-alpha/2));  
  end
   %%%%%%%%%%%%%%%%% dispersion %%%%%%%%%%%%%%%%%%%%%%%
   
   E=fft(ifft(E).*opdisp.');
   
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   if showevolution~=0 && (mod(k,npas/100)==0),
            
       disp([num2str(k*h/L*100) '% completed' '      ' 'beta2=' num2str(betas(1)*1e24)])
 
       figure(9);
       subplot(211); 
       plot(t*1e12,abs(E).^2); set(gca,'ylim',[0 max(abs(E).^2)]);
       set(gca,'xlim',[-1 1]);
       xlabel('t [ps]'); ylabel('P (intracavity) [W]');

       TEsh=fftshift(ifft(E));
       Sp=abs(TEsh).^2;
       mSp = max(Sp);
       PSp = angle(TEsh);
       
       subplot(212); 
       plot(lambda(lambdaidx)*1e9,10*log10(Sp(lambdaidx)./mSp));
       xlabel('Wavelength [nm]'); ylabel('Spectrum');axis ([lambdadisp(1) lambdadisp(2) -30 0])
       %set(gca,'ylim',[-300 0],'ytick',[-300:20:0]);
       drawnow;

   end;
   
 
end;

% update Nc
AVG_N=floor((tau)/dt); % averaging time in n.o.samples
for i=1:length(E)
   Nc=Nc+dt*(D*abs(E(i))^4-Nc/tau);
   Nc_avg=Nc_avg*(AVG_N-1)/AVG_N+Nc/AVG_N;
   %Nc_monitor(counter)=Nc_avg;
   %counter=counter+1;
end


Prop=E;
%toc
end


