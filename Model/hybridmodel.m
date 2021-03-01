% ===============================================================
%% Semiconductor-based mode-locked laser modeling
% Release year: 2020
% Author: Stijn Cuyvers
% Affiliation: Photonics Research Group, Ghent University - imec, Belgium

% Hybrid simulation model:
% - Traveling-wave model for active semiconductor sections
%    (similar to paper "Numerical Study of Dynamical Regime 
% ... in a Monolithic Passively Mode-Locked Semiconductor Laser" 
% ... from A. G. Vladimirov + papers from Javaloyes et al. and Agrawal et al.)
% - Split-step Fourier method for passive laser cavity

% ===============================================================

clear all
close all


% ---------------------------------------------------------------
%% === Parameters ===
% global parameters used for Free-carrier absorption modeling
global Nc
global Nc_avg

Nc=0;
Nc_avg=0;

% -- Physical parameters --
lambda0=1.57*1e-6;              % wavelength [m]
c=3*1e8;                        % speed of light [m/s]
omega0=2*pi*c/lambda0;          % angular frequency [rad s^-1]
n_g=3.85;                       % group index
v_g=c/n_g;                      % group velocity [m/s], assumed constant 
h=6.626070040*1e-34;            % planck's constant [J*s]
hv=h*c/lambda0;                 % photon energy [J] 
e=1.60217662*1e-19;             % electron charge [C]
neff=3;                         % effective mode index in SOA/SA
conf=0.075;                     % confinement factor (dimensionless)
Nt=8.7*1e23;                    % transparency carrier density
Aqw=0.54e-12;                   % Mode area in semiconductor section [m^2]
kappa_1=0.5;                    % reflectivity left mirror
kappa_2=0.99;                   % reflectivity right mirror

% -- Gain parameters --
gamma_g=10^9;                   % carrier relax. rate gain section, [s^-1] 
chi0g=0.070;                    % material gain constant
beta_g=2.56*10^11;              % linear internal losses gain section, [s^-1]
lg=850*1e-6;                    % length of gain section [m]
%**********
Ig=45*1e-3;                     % injection current gain section [A]
%**********
Vg=7.56*1e-17;                  % active volume of gain medium [m^3]
omegag_g=0;                     % bandgap offset (rad s^-1)
omegat_g=90*1e12;               % span of the band (rad s^-1)
intrabandrelax_g=4*1e12;        % intraband relaxation rate (s^-1)

% -- Absorber parameters --
gamma_q=10^11;                  % carrier relax. rate absorber section, [s^-1]
beta_q=2.56*10^11;              % linear internal loss absorber section, [s^-1]
chi0q=0.48;                     % material absorption constant 0.35
lq=60*1e-6;                     % length of absorber section [m]
alpha_q=0;                      % linewidth enhancement factor absorber, not used 
Vq=Vg*(lq/lg);                  % active volume of absorber medium [m^3]
omegag_q=5*1e12;                % bandgap offset (rad s^-1)
omegat_q=90*1e12;               % span of the band (rad s^-1)
intrabandrelax_q=8*1e12;        % intraband relaxation rate (s^-1)

% -- Isolation between SA and SOA --
gamma_iso=10^9;                 % carrier relax. rate isolation section, [s^-1] 
chi0_iso=chi0g;                 % material constant  
beta_iso=2.56*1e11;             % linear internal loss isolation section, [s^-1]
l_iso=30*1e-6;                  % length of isolation section [m]
Viso=Vg*(l_iso/lg);             % active volume of isolation section
omegag_iso=0;                   % bandgap offset (rad s^-1)
omegat_iso=90*1e12;             % span of the band (rad s^-1)
intrabandrelax_iso=8*1e12;      % intraband relaxation rate (s^-1)

% -- Passive waveguide parameters --
lp=14*1e-3;                                                             % length of passive section [m]
loss_passive= 0.7*1e2;                                                  % time-independent loss [dB/m]

beta2=1.3e-24;                                                          % [s^2/m], positive means normal dispersion, negative anomalous
beta3=0.0042e-36;                                                       % [s^3/m], third  order dispersion
betas_passive=[beta2 beta3];
n2=5*1e-18;                                                             % nonlinear coefficient n2 [m^2/W] (according to PhD Bart:  n2=6e-18 m^2/W)
Aeff=0.29*1e-12;                                                        % effective mode area [m^2]
gamma_passive=n2*(2*pi/lambda0)/Aeff;                                   % Kerr nonlinearity parameter gamma [1/ (W*m)]

MODEL_TPA_and_FCA=1;
MODEL_RAMAN=1;

% -- Discretization --
dt=20*1e-15;                                % time step [s] (Recommended stepsize: <40fs for convergence)
dz=dt;                                      % space step [s] (Courant–Friedrichs–Lewy condition)
sim_time=100*1e-9;                          % total simulation time [s]
L=lq+l_iso+lg+lp;                           % total length [m]
L_TWM=lq+l_iso+lg;                          % length of TWM model (gain+absorber region)
RT=2*L/v_g;                                 % roundtriptime [s]
n_o_timesteps=round(sim_time/dt);           % number of time steps
n_o_spacesteps=round((L_TWM/v_g)/dz);       % number of spatial steps
Z=linspace(0,L_TWM,n_o_spacesteps).*1e6;    % space positions [um]
bq=round((lq/v_g)/dz);                      % boundary of absorber
biso=round(((lq+l_iso)/v_g)/dz);            % boundary of isolation section
bg=round(((lg+l_iso+lq)/v_g)/dz);           % boundary of gain

% -- Spectral filter --
delta=1*1e12;                              % bandwidth parameter of spectral filter [Hz] (FWHM=2*delta)
A=log(2)/(delta^2);                        % bandwidth-dependent paramter
filter_memory=10/delta;                    % determines how many samples from the past are considered [s] (this approximation is actually related to the intraband relaxation time)
n_o_filtersamples=round(filter_memory/dt); % number of samples used for filter 
GAUSSIAN=0;                                % choose appropriate filter, gaussian or Lorentzian
if GAUSSIAN==1
    % Gaussian filter:
    % filter samples in steady state (precalculated to boost simulation speed)
    filter_steadystate=(1/(sqrt(2*pi*A))).*exp(-([0:dt:(n_o_filtersamples-1)*dt]).^2./(2*A)); 
    filter_nonsteadystate=@(timearray)(1/(sqrt(2*pi*A))).*exp(-(timearray).^2./(2*A));
else
    % Lorentzian filter
    filter_steadystate=delta*exp(-delta.*[0:dt:(n_o_filtersamples-1)*dt]);     
    filter_nonsteadystate=@(timearray)delta*exp(-delta.*timearray);
end

% -- Simulation settings and parameters --
DEBUG_MODE=0;                               % debug purposes, show additional information during simulation
WRITEDATATOFILE=0;                          % Write monitor data to .csv file     
SAVEFIGURES=1;                              % Save plots
cavitymonitorpos=bg;                        % monitor position for forward-propagating field in cavity
amplitude_threshold=(1e-6)/sqrt(Aqw);       % threshold for field amplitude
nearzero_amplitude=1e-60;                   % if all field samples are below this threshold, stop simulation
FLAG_RETRY=0;                               % if numerical error results, retry simulation with smaller timestep; NOT YET IMPLEMENTED
% ---------------------------------------------------------------

% ---------------------------------------------------------------
%% === Initialization ===
display(['================================================================================='])
display(['Hybrid TWM/Split-step Fourier simulation of semiconductor mode-locked laser.'])
display(['Number of roundtrips simulating: ',num2str(floor(sim_time/RT))])
display(['SOA injection current [mA]: ' num2str(Ig*1e3)])
timertot=0; % timer

figurename=['MLLsim_','Current_',num2str(Ig*1000),'X0g_',num2str(chi0g),'X0q_',num2str(chi0q),'beta2_',num2str(beta2*1e24),'beta3_',num2str(beta3*1e36),'gamma_',num2str(gamma_passive),'TPAandFCA_',num2str(MODEL_TPA_and_FCA),'RAMAN_',num2str(MODEL_RAMAN),'dt_',num2str(dt)];
filename=[figurename,'.txt'];
fid=fopen(filename,'w');
fprintf(fid, [ '-----------------------------------------------------------------------------------', '\n']);
fprintf(fid, [ 'Injection current: ', num2str(Ig), ' A', '\n']);
fprintf(fid, [ 'Number of roundtrips simulated: ', num2str(floor(sim_time/RT)), '\n']);
fprintf(fid, [ 'Timestep [s]: ', num2str(dt), '\n']);
fprintf(fid, [ 'Dispersion beta_2 [ps^2/m]: ', num2str(beta2*1e24), '\n']);
fprintf(fid, [ 'Dispersion beta_3 [ps^3/m]: ', num2str(beta3*1e36), '\n']);
fprintf(fid, [ 'Nonlinearity gamma [1/W/m]: ', num2str(gamma_passive), '\n']);
fprintf(fid, [ 'TPA and FCA included: ', num2str(MODEL_TPA_and_FCA), '\n']);
fprintf(fid, [ 'Raman effect included: ', num2str(MODEL_RAMAN), '\n']);

% -- Variables --
const1=dz*1i*omega0*v_g*conf/(2*neff*c);    % constant used for fields
const2=conf/(Nt*h*neff*c/(2*pi));           % constant used for gain
Reservoir_counter=0;                        % track how full reservoir is
Queue_counter=1;                            % tracks current queue element 

% Define fields, gain, parameter-vectors
Af=(zeros(2,n_o_spacesteps));                       % forward traveling waves
Ab=(zeros(2,n_o_spacesteps));                       % backward traveling waves 
Abpos1=(zeros(1,n_o_filtersamples));                % backward traveling wave at position 1
OUTPUTmonitor_real=(zeros(1,n_o_timesteps));        % output wave monitor
OUTPUTmonitor_imag=(zeros(1,n_o_timesteps));        % output wave monitor
cavitymonitor_real=(zeros(1,n_o_timesteps));        % cavity wave monitor
cavitymonitor_imag=(zeros(1,n_o_timesteps));        % cavity wave monitor
N=zeros(2,n_o_spacesteps);                          % carrier density vector
N_SOA=zeros(1,n_o_timesteps);                       % carrier density monitor of SOA
N_ISO=zeros(1,n_o_timesteps);                       % carrier density monitor of isolation section
N_SA=zeros(1,n_o_timesteps);                        % carrier density monitor of SA section

Is=zeros(1,n_o_spacesteps);                         % current vector
Vs=zeros(1,n_o_spacesteps);                         % active volume vector
Pump=zeros(1,n_o_spacesteps);                       % semiconductor driving term 
Chis=zeros(1,n_o_spacesteps);                       % chi0's vector
betas=zeros(1,n_o_spacesteps);                      % beta vector
gammas=zeros(1,n_o_spacesteps);                     % gamma vector
intrabandrelax=zeros(1,n_o_spacesteps);             % intraband relaxation rate vector
omegag=zeros(1,n_o_spacesteps);                     % omega_g vector (NOT USED)
omegat=zeros(1,n_o_spacesteps);                     % omega_t vector

% Define currents
Is(1,1:bq)=zeros(1,bq);
Is(1,bq+1:biso)=zeros(1,biso-bq);
Is(1,biso+1:bg)=Ig*ones(1,bg-biso);

% Define active volume
Vs(1,1:bq)=Vq*ones(1,bq);
Vs(1,bq+1:biso)=Viso*ones(1,biso-bq);
Vs(1,biso+1:bg)=Vg*ones(1,bg-biso);

% Define driving term
Pump=Is./(Nt*e.*Vs);

% Define material chi0 values
Chis(1,1:bq)=chi0q*ones(1,bq);
Chis(1,bq+1:biso)=chi0_iso*ones(1,biso-bq);
Chis(1,biso+1:bg)=chi0g*ones(1,bg-biso);

% Define beta vector
betas(1:bq)=exp(-0.5*dz*beta_q.*ones(1,bq));
betas(bq+1:biso)=exp(-0.5*dz*beta_iso.*ones(1,biso-bq));
betas(biso+1:bg)=exp(-0.5*dz*beta_g.*ones(1,bg-biso));

% Define gamma vector
gammas(1:bq)=gamma_q.*ones(1,bq);
gammas(bq+1:biso)=gamma_iso.*ones(1,biso-bq);
gammas(biso+1:bg)=gamma_g.*ones(1,bg-biso);

% Define intraband relaxation times
intrabandrelax(1:bq)=intrabandrelax_q.*ones(1,bq);
intrabandrelax(bq+1:biso)=intrabandrelax_iso.*ones(1,biso-bq);
intrabandrelax(biso+1:bg)=intrabandrelax_g.*ones(1,bg-biso);

% Define omega_g's 
omegag(1:bq)=0-omegag_q.*ones(1,bq);
omegag(bq+1:biso)=0-omegag_iso.*ones(1,biso-bq);
omegag(biso+1:bg)=0-omegag_g.*ones(1,bg-biso);

% Define omeg_t's
omegat(1:bq)=omegat_q.*ones(1,bq);
omegat(bq+1:biso)=omegat_iso.*ones(1,biso-bq);
omegat(biso+1:bg)=omegat_g.*ones(1,bg-biso);

% Initialize Reservoir and Queue for hybrid model
n_o_spacesteps_passive=round(2*lp/(v_g*dt)); % twice length of passive cavity because roundtrip!
Reservoir_capacity=n_o_spacesteps_passive;
Reservoir=zeros(1,Reservoir_capacity);
Queue=zeros(1,n_o_spacesteps_passive);

% ~~~ choose initial gain values ~~~
sigma=0.1;
N(1,1:bq)=0.5*ones(1,bq)+sigma*( randn(1,bq) + 1i*randn(1,bq)   ) ;
N(1,bq+1:biso)=1*ones(1,biso-bq)+sigma*( randn(1,biso-bq) + 1i*randn(1,biso-bq)   ) ;
N(1,biso+1:bg)=2*ones(1,bg-biso)+sigma*( randn(1,bg-biso) + 1i*randn(1,bg-biso)   ) ;

sigma=0.1;
Af(1,:)=sigma*( randn(1,n_o_spacesteps) + 1i*randn(1,n_o_spacesteps)   ) ;
Ab(1,:)=sigma*( randn(1,n_o_spacesteps) + 1i*randn(1,n_o_spacesteps)   ) ;
Queue(1,:)=sigma*( randn(1,n_o_spacesteps_passive) + 1i*randn(1,n_o_spacesteps_passive)   ) ;

%{
% Useful in case you want to use specific initial conditions
Af(1,:)=table2array(readtable('Af.csv')); %Afinitialization;
Ab(1,:)=table2array(readtable('Ab.csv'));%Abinitialization;
Queue(1,:)=table2array(readtable('Queue.csv'));%Queueinitialization;
N(1,:)=table2array(readtable('N.csv'));%Ninitialization;
%}

% Initialize values at position 1 of grid
Abpos1(1)=Ab(1,1);
OUTPUTmonitor_real(1)=sqrt(1-kappa_1)*real(Ab(1,1));
OUTPUTmonitor_imag(1)=sqrt(1-kappa_1)*imag(Ab(1,1));
cavitymonitor_real(1)=real(Ab(1,cavitymonitorpos));
cavitymonitor_imag(1)=imag(Ab(1,cavitymonitorpos));
N_SA(1)=N(1,round(bq/2));
N_ISO(1)=N(1,round(bq+(biso-bq)/2));
N_SOA(1)=N(1,round(biso+(bg-biso)/2));
% ---------------------------------------------------------------

% ---------------------------------------------------------------
%% === Actual simulation ===

% LOOP
tic;                                        % start timer
timer_cycles=zeros(1,100);
percentage=1;                               % tracks simulation completion
percentage_cycle=floor(n_o_timesteps/100);
disp(['Simulation started...']);

for i=2:n_o_timesteps
    if mod(i,percentage_cycle)==0
        timer_cycles(percentage)=toc;
        timertot=timertot+timer_cycles(percentage);
        disp(['...  ',num2str(percentage), '% complete   [', num2str(timer_cycles(percentage)),'s]']);
        percentage=percentage+1;
        tic;
    end
    
    % check validity of data
    if any(isnan(Af(1,:)))
        disp('Error: field contains NaN elements.')
        break;
    elseif any(isnan(Ab(1,:)))
        disp('Error: field contains NaN elements.')
        break;
    elseif any(isnan(N(1,:)))
        disp('Error: gain contains NaN elements.')
        break;
    end
    
    % check if signals are converging to zero 
    if mod(i,percentage_cycle)==0 && max(abs(Af(1,:)))<nearzero_amplitude && max(abs(Ab(1,:)))<nearzero_amplitude && max(abs(Queue))<nearzero_amplitude
        disp('Signal converges to zero, no mode-locking regime obtained.')
        break;
    end
    
    % Check if Reservoir is full
    % ---- SPLIT-STEP FOURIER /  ---
    if Reservoir_counter==Reservoir_capacity
        % reservoir full - select segment for split-step propagation    
        if any(isnan(Reservoir))
            % check if Reservoir has no NaN values
            disp('Error: Reservoir contains NaN elements.')
            break;
        end
        
        % algorithm:
        % generate grid of 100 cells over the entire reservoir, calculate
        % normalized energy for each cell
        % find peaks of grid
        % CASES:
        %   1) if first peak is in cell > 50 (further than second half of grid)
        %       -> take segment end before the first peak, as such that the 'center
        %       of gravity' (in terms of energy) shifts to the center of the grid
        %   2) else, if only 1 peak
        %       -> take segment end after the peak, as such that energy is as much
        %       centralized in grid as possible
        %   3) else, if multiple peaks
        %       -> iterate points in between peaks, and after last peak
        %        for each calculate energy centralization, take best one


        % construct energy grid
        cellgrid=zeros(1,100);
        tot_energy=sum(Aqw.*abs(Reservoir).^2);
        step=round(Reservoir_capacity/100);
        if 99*step>Reservoir_capacity
            step=floor(Reservoir_capacity/100);
        end
        for m=1:99
            cellgrid(m)=sum(Aqw.*abs(Reservoir(1+(m-1)*step:m*step)).^2);
        end
        cellgrid(100)=sum(Aqw.*abs(Reservoir(1+(99-1)*step:Reservoir_capacity)).^2);

        % find peaks of cellgrid
        minpeakheight=max(cellgrid)/100;
        [peaks,peakpos]=findpeaks(cellgrid,'MinPeakHeight',minpeakheight);

        % CASES
        % case 1: there are no peaks -> no spit-step fourier
        % case 2: if first peak is in second half of reservoir:
        %   -> take segment end before peak so that peak is in center
        % case 3: there is 1 peak
        %   -> if peak position is in first half of reservoir, compare energy
        %   centralization for case 1) samen distance after as before peak 2)
        %   entire reservoir
        % case 4: multiple peaks:
        %   -> iterate all minima and end of reservoir, take the one with best
        %   energy centralization

        if isempty(peaks)==1
            % case: no peaks are found
            max_pos=find(abs(Reservoir)==max(abs(Reservoir)));
            if max_pos==Reservoir_capacity
                min_pos=round(Reservoir_capacity/2);
            else
                min_pos=Reservoir_capacity;
            end
            FLAG_PULSE=0;
        else
            peakpostolerance=5; % tolerance for peak centralization
            if peakpos(1)>50+peakpostolerance % if first peak is in second half of grid
                energycentrum=round(sum([1:1:100].*cellgrid)/tot_energy);
                cell_pos=energycentrum-50;
                if cell_pos<=0
                    % case: first peak is in second half, but energy
                    % center of mass is not
                    cell_pos=min(2*energycentrum,100);
                    if cell_pos==100
                        minimum=min(abs(Reservoir(1+(cell_pos-1)*step:Reservoir_capacity)));
                    else
                        minimum=min(abs(Reservoir(1+(cell_pos-1)*step:cell_pos*step)));
                    end
                    min_pos=find(abs(Reservoir)==minimum);

                    % check energy centralization to determine if spit-step
                    % fourier can be used
                    etot = sum(Aqw.*abs(Reservoir(1:min_pos)).^2);
                    interval=round(min_pos*0.25);
                    ecenter= sum(Aqw.*abs(Reservoir(interval:min_pos-interval)).^2);
                    if ecenter/etot>=0.999
                        FLAG_PULSE=1;
                    else
                        FLAG_PULSE=0;
                    end
                else
                    minimum=min(abs(Reservoir(1+(cell_pos-1)*step:cell_pos*step)));
                    min_pos=find(abs(Reservoir)==minimum);
                    FLAG_PULSE=0;
                end 
            elseif length(peaks)==1
                % case: only 1 peak, but this peak is in already in first
                % half of grid or already at center of grid
                cell_pos=min(2*peakpos(1),100);
                if cell_pos==100
                    minimum=min(abs(Reservoir(1+(cell_pos-1)*step:Reservoir_capacity)));
                else
                    minimum=min(abs(Reservoir(1+(cell_pos-1)*step:cell_pos*step)));
                end
                min_pos=find(abs(Reservoir)==minimum);

                % check energy centralization to determine if spit-step
                % fourier can be used
                etot = sum(Aqw.*abs(Reservoir(1:min_pos)).^2);
                interval=round(min_pos*0.25);
                ecenter= sum(Aqw.*abs(Reservoir(interval:min_pos-interval)).^2);
                if ecenter/etot>=0.999
                    FLAG_PULSE=1;
                else
                    FLAG_PULSE=0;
                end

            else
                % case: multiple peaks
                % start looking from cell in between first two peaks
                startcell=floor((peakpos(1)+peakpos(2))/2);
                centralenergyfraction=0;
                optcellposition=0;
                % iterate to find out best segmentation
                for m=startcell:100
                    temp_totenergy=sum(cellgrid(1:m));
                    interval=round(m*0.25);
                    temp_centralenergyfraction=sum(cellgrid(interval:m-interval))/temp_totenergy;
                    % critertion: highest energy fraction in center 
                    if temp_centralenergyfraction>centralenergyfraction
                        optcellposition=m;
                    end
                end
                if optcellposition==100
                    minimum=min(abs(Reservoir(1+(optcellposition-1)*step:Reservoir_capacity)));
                else
                    minimum=min(abs(Reservoir(1+(optcellposition-1)*step:optcellposition*step)));
                end
                min_pos=find(abs(Reservoir)==minimum);
                
                % check energy centralization to determine if spit-step
                % fourier can be used
                etot = sum(Aqw.*abs(Reservoir(1:min_pos)).^2);
                interval=round(min_pos*0.25);
                ecenter= sum(Aqw.*abs(Reservoir(interval:min_pos-interval)).^2);
                if ecenter/etot>=0.999
                    FLAG_PULSE=1;
                else
                    FLAG_PULSE=0;
                end
            end
        end
       
        if FLAG_PULSE==0 
            % if segment is considered 'noise' rather than 'pulse',
            % ... simply model loss, (phase shift if desired)(do not use split-step)
            Propagated_field=sqrt(kappa_2).*Reservoir(1:min_pos).*exp(-2*lp*loss_passive/4.343/2).*exp(1i*2*pi*2*lp/lambda0/n_g);
            
            if DEBUG_MODE 
                disp(['No split-step Fourier; ']);                   
                plot([1:1:Reservoir_capacity],abs(Reservoir).^2);
                hold on;
                plot([1:1:min_pos],abs(Propagated_field).^2);
                v=get(gca, 'ylim');
                patch([1 min_pos min_pos 1],[v(1) v(1) v(2) v(2)],'k','FaceAlpha',.09,'EdgeColor','none')
                hold off;
            end
            
        else
            % if segment is considered as pulse: use split step fourier
            Tspan=min_pos*dz; % time span of signal trace [s]
            if sqrt(kappa_2)>0.98
                Propagated_field=sqrt(kappa_2).*SSF(Reservoir(1:min_pos).*sqrt(Aqw),Tspan,lambda0,2*lp,betas_passive,loss_passive,n2,Aeff,MODEL_TPA_and_FCA,MODEL_RAMAN);
            else
                Propagated_field=sqrt(kappa_2).*SSF(Reservoir(1:min_pos).*sqrt(Aqw),Tspan,lambda0,lp,betas_passive,loss_passive,n2,Aeff,MODEL_TPA_and_FCA,MODEL_RAMAN);
                Propagated_field=SSF(Propagated_field,Tspan,lambda0,lp,betas_passive,loss_passive,n2,Aeff,MODEL_TPA_and_FCA,MODEL_RAMAN,SSFstep); 
            end
            Propagated_field=Propagated_field./sqrt(Aqw); % convert to Watt/m^2 again
            
            if DEBUG_MODE 
                disp(['Split-step Fourier; ']);                    
                plot([1:1:Reservoir_capacity],abs(Reservoir).^2);
                hold on;
                plot([1:1:min_pos],abs(Propagated_field).^2);
                v=get(gca, 'ylim');
                patch([1 min_pos min_pos 1],[v(1) v(1) v(2) v(2)],'k','FaceAlpha',.09,'EdgeColor','none')
                hold off;
            end
  
        end

        % Clean reservoir
        Reservoir=[Reservoir(min_pos+1:Reservoir_capacity),zeros(1,min_pos)];
        Reservoir_counter=Reservoir_capacity-min_pos;

        % add propagated field to queue of backward prop. field
        Queue(1:min_pos)=Propagated_field;
        Queue_counter=1;
    end    
    

    % prepare filter samples
    if i==n_o_filtersamples
        filter=filter_steadystate;
    elseif i<n_o_filtersamples
        filter=filter_nonsteadystate([0:dt:(i-1)*dt]);
    end
    
    % TRAVELING-WAVE MODEL
    % update forward- and backward propagating fields in TWM section
    % to save memory, not all samples are saved 
    % only samples at time t (index i) and previous time (index i-1)
    % are saved
    % i   -> 2-mod(i,2)
    % i-1 -> mod(i,2)+1


    % update forward propagating field
    % Original: temp_factor=(log(1+1i*(N(mod(i,2)+1,1:n_o_spacesteps-1)))+log(1+1i*(N(mod(i,2)+1,2:n_o_spacesteps)))-log((intrabandrelax(1:n_o_spacesteps-1)+1i*omegat(1:n_o_spacesteps-1))./intrabandrelax(1:n_o_spacesteps-1)));
    
    temp_factor=2*log(1+1i*(N(mod(i,2)+1,1:n_o_spacesteps-1))./(1-1i*omegag(1:n_o_spacesteps-1)./intrabandrelax(1:n_o_spacesteps-1)))-log(1-omegat(1:n_o_spacesteps-1)./(omegag(1:n_o_spacesteps-1)+1i*intrabandrelax(1:n_o_spacesteps-1)));
    Af(2-mod(i,2),2:n_o_spacesteps)=Af(mod(i,2)+1,1:n_o_spacesteps-1).*betas(1:n_o_spacesteps-1).*exp(-const1.*Chis(1:n_o_spacesteps-1).*temp_factor);
    % update backward propagating field
    Ab(2-mod(i,2),1:n_o_spacesteps-1)=Ab(mod(i,2)+1,2:n_o_spacesteps).*betas(2:n_o_spacesteps).*exp(-const1.*Chis(2:n_o_spacesteps).*temp_factor);

    % for backward propagating field at position 1, save more past samples because this is needed for the convolution with filter
    Abpos1=circshift(Abpos1,[0,1]);
    Abpos1(1)=Ab(2-mod(i,2),1);

    % boundary conditions
    if i>n_o_filtersamples
        Af(2-mod(i,2),1)=sqrt(kappa_1)*dt*sum(filter.*(Abpos1));
        
        % TEST (can be ignored)
        % plan to use distributed filter instead of single filter (future work)
        %x1=[filter zeros(1,length(Abpos1)-1)];
        %y1=[flip(Abpos1) zeros(1,length(filter)-1)];
        
        %tempvalue=sqrt(kappa_1)*dt*ifft(fft(x1).*fft(y1));
        %Af(2-mod(i,2),1)=tempvalue(length(filter));
        
        %{
        if abs(sqrt(kappa_1)*dt*sum(filter.*(Abpos1))-tempvalue(length(filter)))>1e-12 
            disp('FILTER ERROR')
            convolution=sqrt(kappa_1)*dt*sum(filter.*(Abpos1))
            Fourierdomain=tempvalue(length(filter))
            difference=convolution-Fourierdomain
            
        end
        %}
    else
        Af(2-mod(i,2),1)=sqrt(kappa_1)*dt*sum(filter.*(Abpos1(1:i)));
        
        % TEST (can be ignored)
        % plan to use distributed filter instead of single filter (future work)
        %x1=[filter zeros(1,length(Abpos1(1:i))-1)];
        %y1=[flip(Abpos1(1:i)) zeros(1,length(filter)-1)];
        
        %tempvalue=sqrt(kappa_1)*dt*ifft(fft(x1).*fft(y1));
        %Af(2-mod(i,2),1)=tempvalue(length(filter));

        %{
        if abs(sqrt(kappa_1)*dt*sum(filter.*(Abpos1(1:i)))-tempvalue(length(filter)))>1e-12 
            disp('FILTER ERROR')
            convolution=sqrt(kappa_1)*dt*sum(filter.*(Abpos1(1:i)))
            Fourierdomain=tempvalue(length(filter))
            difference=convolution-Fourierdomain
        end
        %}
    end


    % update output and cavity monitor
    OUTPUTmonitor_real(i)=real(Af(2-mod(i,2),1)*sqrt(1-kappa_1)/sqrt(kappa_1));
    OUTPUTmonitor_imag(i)=imag(Af(2-mod(i,2),1)*sqrt(1-kappa_1)/sqrt(kappa_1));

    cavitymonitor_real(i)=real(Ab(2-mod(i,2),cavitymonitorpos));
    cavitymonitor_imag(i)=imag(Ab(2-mod(i,2),cavitymonitorpos));

    % concat signal samples in queue to backward prop. field
    Ab(2-mod(i,2),n_o_spacesteps)=Queue(Queue_counter);
    Queue_counter=Queue_counter+1;

    % update gain
    N(2-mod(i,2),:)=N(mod(i,2)+1,:)+dt.*(Pump-gammas.*N(mod(i,2)+1,:)+const2.*(abs((Af(2-mod(i,2),:))+(Ab(2-mod(i,2),:))).^2).*(-Chis).*(2*atan(N(mod(i,2)+1,:))-atan(omegat./intrabandrelax)));
    
    % update carrier monitors
    N_SA(i)=N(2-mod(i,2),round(bq/2));
    N_ISO(i)=N(2-mod(i,2),round(bq+(biso-bq)/2));
    N_SOA(i)=N(2-mod(i,2),round(biso+(bg-biso)/2));

    % fill reservoir
    Reservoir_counter=Reservoir_counter+1;
    Reservoir(Reservoir_counter)=Af(mod(i,2)+1,n_o_spacesteps);

end
% ---------------------------------------------------------------

% ---------------------------------------------------------------
%% == Process simulation results: Data visualization ==
% Plot settings (colors, fontsize, etc.)
default_blue=[0 0.4470 0.7410];
default_red=[178/256 24/256 43/256];
default_darkgrey=[100/256,100/256,100/256];
Fontsize=14;

% Define complex output and cavity monitor
OUTPUTmonitor=(OUTPUTmonitor_real+1i*OUTPUTmonitor_imag)*sqrt(Aqw);
cavitymonitor=(cavitymonitor_real+1i*cavitymonitor_imag)*sqrt(Aqw);

% write data to file
if WRITEDATATOFILE==1
    csvfilename=[figurename,'.csv'];
    csvwrite(['OUTPUT_',csvfilename],OUTPUTmonitor);
    csvwrite(['CAVITY_',csvfilename],cavitymonitor);
    csvwrite(['N_SOA',csvfilename],N_SOA);
    csvwrite(['N_SA',csvfilename],N_SA);
    csvwrite(['N_ISO',csvfilename],N_ISO);
end

% Time domain trace over entire simulation time
figure;
OUTPUTmonitor_selection=OUTPUTmonitor(1:n_o_timesteps);
taxis=linspace(0,n_o_timesteps*dt*1e12,length(OUTPUTmonitor_selection));
plot(taxis./1e3,abs(OUTPUTmonitor_selection).^2,'Color',default_blue);
xlabel('Time (ns)', 'FontSize', Fontsize);
ylabel('Envelope (W)', 'FontSize', Fontsize);
set(gca,'FontSize',Fontsize);
title('Output pulse train', 'FontSize', Fontsize);
xlim([0 max(taxis)./1e3]);
grid on;
set(gca,'box','off');
if SAVEFIGURES==1 saveas(gcf,[figurename,'_1','.png']); end;

% Time domain trace over last 2ns (to show a few pulses)
time_span=2000*1e-12;
figure;
OUTPUTmonitor_selection=OUTPUTmonitor(n_o_timesteps-round(time_span/dt):n_o_timesteps);
taxis=linspace((n_o_timesteps-round(time_span/dt))*dt*1e12,n_o_timesteps*dt*1e12,length(OUTPUTmonitor_selection));
plot(taxis./1e3,abs(OUTPUTmonitor_selection).^2,'Color',default_blue);
xlabel('Time (ns)', 'FontSize', Fontsize);
ylabel('Envelope (W)', 'FontSize', Fontsize);
set(gca,'FontSize',Fontsize);
title('Output pulse train', 'FontSize', Fontsize);
xlim([min(taxis)./1e3 max(taxis)./1e3]);
grid on;
set(gca,'box','off');
if SAVEFIGURES==1 saveas(gcf,[figurename,'_2','.png']); end;

% Time domain trace of normalized carrier densities
figure;
taxis=linspace(0,n_o_timesteps*dt*1e12,n_o_timesteps);
plot(taxis./1e3, N_SOA,'Color',default_darkgrey);
hold on;
plot(taxis./1e3, N_SA);
plot(taxis./1e3, N_ISO);
xlabel('Time (ns)', 'FontSize', Fontsize);
ylabel('Envelope (W)', 'FontSize', Fontsize);
set(gca,'FontSize',Fontsize);
xlim([min(taxis)./1e3 max(taxis)./1e3]);
grid on;
set(gca,'box','off');
legend('SOA','SA','ISO');
title('Normalized carrier densities at the center of each section', 'FontSize', Fontsize);
if SAVEFIGURES==1 saveas(gcf,[figurename,'_6','.png']); end;

% Plot spectrum
%{
% spectrum
[ validdata, lambdas, spectrum, lambdas_envelope, envelope, lambdas_phase, phases ] = getspectrum( OUTPUTmonitor, lambda0, dt );
if validdata==1
    freq_center_pos=find(spectrum==max(spectrum));
    lambda_center=lambdas(freq_center_pos);
    figure;
    plot(lambdas,spectrum,'Color',[178/256 24/256 43/256]);
    xlabel('Wavelength (nm)');
    ylabel('Power spectral density (dBm/Hz)');
    set(gca,'FontSize',Fontsize);
    set(gca,'box','off');
    xlim([lambda_center-5 lambda_center+5]);
    set(gca,'box','off');
    ax=gca();
    tick_size=0.010;
    ax.TickLength = [tick_size, tick_size]; % Make tick marks longer.
    ax.LineWidth = 100*tick_size; % Make tick marks thicker.
    set(gca, 'Layer', 'top')
    if SAVEFIGURES==1 saveas(gcf,[figurename,'_7','.png']); end;
end
%}

%{
% In case phase of spectrum is to be visualized
figure;
yyaxis left;
plot(lambdas,spectrum);
ylabel('Power spectral density (dBm/Hz)');
yyaxis right;
plot(lambdas_phase,(phases)./pi,'.');
xlim([lambda_center-5 lambda_center+5]);
xlabel('Wavelength (nm)');
ylabel('Modal phase (\pi)');
set(gca,'FontSize',Fontsize);
set(gca,'box','off');
set(gca,'box','off');
ax=gca();
tick_size=0.010;
ax.TickLength = [tick_size, tick_size]; % Make tick marks longer.
ax.LineWidth = 100*tick_size; % Make tick marks thicker.
set(gca, 'Layer', 'top')
if SAVEFIGURES==1 saveas(gcf,[figurename,'_8','.png']); end;
%}

time_span=2000*1e-12;
OUTPUTmonitor_selection=OUTPUTmonitor(n_o_timesteps-round(time_span/dt):n_o_timesteps);
[peaks,peakpositions]=findpeaks(abs(OUTPUTmonitor_selection).^2, 'MinPeakHeight',max(abs(OUTPUTmonitor_selection).^2)/10);
MLL_regime=1; %0: Q-switching/noisy/chaotic state, 1: fundamental mll, 2: harmonic mll
if isempty(peaks)==0 && length(peaks)>2 && max(abs(OUTPUTmonitor_selection).^2)>1e-3
    pos_1=0.5*(peakpositions(length(peaks)-1)+peakpositions(length(peaks)-2));   
    pos_2=0.5*(peakpositions(length(peaks)-1)+peakpositions(length(peaks)));   

    steadystatepulsenergy=sum(abs(OUTPUTmonitor_selection(pos_1:pos_2)).^2)*dt*1e12;
    pulsereprate=1e-9/((pos_2-pos_1)*dt);
    display(['steady-state pulse energy: ' num2str(steadystatepulsenergy) ' pJ']);
    display(['pulse repetition rate: ' num2str(pulsereprate) ' [GHz]']);

    if pulsereprate>1.1*(1/RT)*1e-9
        MLL_regime=2; % harmonic mode-locking
    end

    % determine pulsewidth of output pulses
    pulsewidth=-1;
    startpeakpos=peakpositions(length(peaks)-1);
    peakpower=abs(OUTPUTmonitor_selection(startpeakpos)).^2;
    for i=startpeakpos:length(OUTPUTmonitor_selection)
        if abs(OUTPUTmonitor_selection(i))^2<=peakpower/2
            pulsewidth=(2*(i-startpeakpos)*dt)*1e12;
            break;
        end
    end
   
    % plot single output pulse
    nr_samples_selection=ceil((10*pulsewidth*1e-12)/dt);
    lower_bound=max(1,round(startpeakpos-nr_samples_selection/2));
    upper_bound=min(length(OUTPUTmonitor_selection),round(startpeakpos+nr_samples_selection/2));
    toffset=(length(OUTPUTmonitor)-length(OUTPUTmonitor_selection))*dt*1e12;
    taxis=linspace(toffset+lower_bound*dt*1e12,toffset+upper_bound*dt*1e12,nr_samples_selection+1);

    % Fit sech to output pulse
    xdata=taxis./1e3;
    ydata=abs(OUTPUTmonitor_selection(lower_bound:upper_bound)).^2;
    P0=max(ydata);
    fun=@(x,xdata) P0.*sech((x(1)-xdata)/x(2)).^2;
    x0=[xdata(ydata==max(ydata)) 1e-3*pulsewidth/1.763];
    if length(xdata)==length(ydata)
        x=lsqcurvefit(fun,x0,xdata,ydata,[],[]);
        fit_sech=fun(x,xdata);

        pulsewidth=x(2)*1.763*1e3;
        display(['Output pulse width: ' num2str(pulsewidth) ' [ps]']);

        figure;
        plot(xdata,ydata,'Color',default_blue);
        hold on;
        plot(xdata,fit_sech,'r--');
        xlabel('Time (ns)','FontSize',Fontsize);
        ylabel('Envelope (W)','FontSize',Fontsize);
        title('Output pulse','FontSize',Fontsize);
        grid on;
        legend('Simulation','Sech^2 fit');
        set(gca,'FontSize',Fontsize);
        set(gca,'box','off');
        if SAVEFIGURES==1 saveas(gcf,[figurename,'_3','.png']); end;
    end
    
    % determine pulse width intracavity pulse 
    intracavpulsewidth=-1;
    cavitymonitor_selection=cavitymonitor(n_o_timesteps-round(time_span/dt):n_o_timesteps);
    [peaks_intracavity,peakpositions_intracavity]=findpeaks(abs(cavitymonitor_selection).^2, 'MinPeakHeight',max(abs(cavitymonitor_selection).^2)/10);
    startpeakpos=peakpositions_intracavity(length(peaks_intracavity)-1);
    intracavitypeakpower=abs(cavitymonitor_selection(startpeakpos)).^2;
    for i=startpeakpos:length(cavitymonitor_selection)
        if abs(cavitymonitor_selection(i))^2<=intracavitypeakpower/2
            intracavpulsewidth=(2*(i-startpeakpos)*dt)*1e12;
            break;
        end
    end

    % plot single intracavity pulse
    nr_samples_selection=ceil((10*intracavpulsewidth*1e-12)/dt);
    lower_bound=max(1,round(startpeakpos-nr_samples_selection/2));
    upper_bound=min(length(cavitymonitor_selection),round(startpeakpos+nr_samples_selection/2));
    toffset=(length(cavitymonitor)-length(cavitymonitor_selection))*dt*1e12;
    taxis=linspace(toffset+lower_bound*dt*1e12,toffset+upper_bound*dt*1e12,nr_samples_selection+1);

    % Fit sech to intracavity pulse
    xdata=taxis./1e3;
    ydata=abs(cavitymonitor_selection(lower_bound:upper_bound)).^2;
    xdata=taxis./1e3;
    ydata=abs(cavitymonitor_selection(lower_bound:upper_bound)).^2;
    P0=max(ydata);
    fun=@(x,xdata) P0.*sech((x(1)-xdata)/x(2)).^2;
    x0=[xdata(ydata==max(ydata)) 1e-3*intracavpulsewidth/1.763];
    if length(xdata)==length(ydata)
        x=lsqcurvefit(fun,x0,xdata,ydata,[],[]);
        fit_sech=fun(x,xdata);

        intracavpulsewidth=x(2)*1.763*1e3;
        display(['Intracavity pulse width: ' num2str(intracavpulsewidth) ' [ps]']);

        figure;
        plot(taxis./1e3,abs(cavitymonitor_selection(lower_bound:upper_bound)).^2,'Color',default_blue);
        hold on;
        plot(xdata,fit_sech,'r--');
        xlabel('Time (ns)','FontSize',Fontsize);
        ylabel('Envelope (W)','FontSize',Fontsize);
        title('Intracavity pulse','FontSize',Fontsize);
        grid on;
        set(gca,'FontSize',Fontsize);
        set(gca,'box','off');
        legend('Simulation','Sech^2 fit');
        if SAVEFIGURES==1 saveas(gcf,[figurename,'_5','.png']); end;
    end
else
    MLL_regime=0;
end    

% Write results to .txt file
if MLL_regime==0
    fprintf(fid, [ 'MLL regime: Q-switching', '\n']);
elseif MLL_regime==2
    fprintf(fid, [ 'MLL regime: Harmonic mode-locking', '\n']);
    fprintf(fid, [ 'Pulse repetition rate: ', num2str(pulsereprate), ' GHz', '\n']);
else
    fprintf(fid, [ 'Output pulse steady-state energy: ', num2str(steadystatepulsenergy), ' pJ', '\n']);
    fprintf(fid, [ 'Output pulse steady-state peak power: ', num2str(peakpower), ' W', '\n']);
    fprintf(fid, [ 'Intracavity pulse steady-state peak power: ', num2str(intracavitypeakpower), ' W', '\n']);
    fprintf(fid, [ 'Pulse repetition rate: ', num2str(pulsereprate), ' GHz', '\n']);
    fprintf(fid, [ 'Output Pulse width: ', num2str(pulsewidth), ' ps', '\n']);
    fprintf(fid, [ 'Intracavity pulse width: ', num2str(intracavpulsewidth), ' ps', '\n']);
end
   

display(['total simulation time: ', num2str(timertot), ' [s]']);  % stop timer
if DEBUG_MODE
    % plot evolution of simulation time/cycle
    figure;
    title('Simulation timing analysis');
    plot([1:1:100],timer_cycles);
    xlabel('Simulation cycle');
    ylabel('Simulation time/cycle [s]');
    grid on;
end
display(['================================================================================='])
% ---------------------------------------------------------------


