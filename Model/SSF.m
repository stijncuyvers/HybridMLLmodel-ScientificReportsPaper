% ===============================================================
%  FUNCTION: Split-step Fourier waveguide propagation
%  Release year: 2020
%  Affiliation: Photonics Research Group, Ghent University - imec, Belgium
% ===============================================================

function Prop=SSF(Signal_in,T_span,lambda0, L,betas,loss,n2,Aeff,MODEL_TPA_and_FCA,MODEL_RAMAN)
% REMARK: normally, trace should have a length that is a power of 2 for the
% FFT, MATLAB 2011 and newer can however handle non power of 2 signal
% lengths
% It is also possible to interpolate the data to match a power of 2 length,
% but care should be taken not to distort phase and or amplitude


%disp(['Split-step window size: ',num2str(T_span*1e12), ' ps']);

% -- ARGUMENTS --
% Signal_in: input signal
% T_span: corresponding time grid
% lambda0: wavelength [nm]
% L: length of propagation
% betas: dispersion
% loss: linear losses
% n2: nonlinear Kerr nonlinearity [m^2/W]
% Aeff: effective mode area [m^2]
% MODEL_TPA_and_FCA: model TPA and FCA
% MODEL_RAMAN: include Raman effect

global Nc
global Nc_avg

%global Nc_monitor
%global counter

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

tau = 1e-9;                 % Free carrier lifetime [s]
planck = 6.626e-34;
D = 2*pi*btpa/(2*planck*w0*Aeff^2);

original_length=size(Signal_in,2);
nTime=original_length;
alpha = log(10.^(loss/10));     % attenuation coefficient
step      = L/npas;               % propagation step
T = linspace(-T_span/2,T_span/2,nTime);     % time grid
dT = T_span/nTime;                          % time grid spacing [ps]
V = 2*pi*(-nTime/2:nTime/2-1)/T_span;       % angular frequency deviation from w0 [10^12 rad/s]%


if step>200e-6
   % require minimum stepsize of 200 microns
   step=200e-6;
   npas=round(L/step);
end
  

%% ----------------------------------------------------------------------


% input signal
A = Signal_in;


%%%%%%%%%%%%% Raman response %%%%%%%%%%%%%%%%
wr = 15.6e12*2*pi; % Raman shift (wr/2pi=15.6 THz), denoted as OMEGA_R
dw = 105e9*pi; % Raman gain bandwidth (more precisely: dw/pi=FWHM of the Raman-gain spectrum), denoted as GAMMA_R
gRmax = 1.5e-10; %m/W % Raman gain coefficient
fR = gRmax*dw/(wr*k0*n2); % see Nonlinear optical phenomena in silicon waveguides: modeling and applications
if MODEL_RAMAN==0
    fR = 0;
end

% see Nonlinear optical phenomena in silicon waveguides: modeling and applications
% Define Raman response
tau1 = 1/sqrt(wr^2-dw^2); 
tau2 = 1/dw; % Raman response function fitting parameters
hR = (tau1^2+tau2^2)/(tau1*tau2^2) * exp(-T/tau2) .* sin(T/tau1);    % Raman time response fit
hR(T<0) = 0;    % heaviside function for causality
hR = hR/trapz(T,hR);    % normalize response function


 % Define propagation constant
B = 0;
for i=1:length(betas)
    B = B + betas(i)/factorial(i+1)*V.^(i+1);
end
% Define propagation operator in frequency domain
Disp = fftshift(1i*B - alpha/2- sigma/2*(1+1i*mu)*Nc_avg);  % V=0 NOT central --> compatible for multiplication with fft(A)

% Iterate over fiber segments
for i = 1:npas
    % First dispersion + losses...
    A = ifft(exp(step*Disp/2).*fft(A));
    
    % ... then nonlinearity: Kerr effect operator + Raman contribution
     K = 1i*gamma*((1-fR)*abs(A).^2 + fR*nTime*dT*fft(ifft(fftshift(hR)).*ifft(abs(A).^2)));
     A = exp(step*K).*A;
      
    % ... and again dispersion + losses.
    A = ifft(exp(step*Disp/2).*fft(A));   
end

% Note: for now an 'average' Nc is used, for more accurate results,
% ...one could include a fitted model to more accurately capture the 
% ...free-carrier density at each location of the waveguide

% update Nc 
AVG_N=floor((tau)/dT); % averaging time in n.o.samples
for i=1:length(A)
   Nc=Nc+dT*(D*abs(A(i))^4-Nc/tau);
   Nc_avg=Nc_avg*(AVG_N-1)/AVG_N+Nc/AVG_N;
   %Nc_monitor(counter)=Nc_avg;
   %counter=counter+1;
end


Prop=A;

end


