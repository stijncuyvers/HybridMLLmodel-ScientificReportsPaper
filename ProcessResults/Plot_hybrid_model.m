% =====================
% Plot simulated data
% =====================
close all
clear all

% Plot settings (colors, fontsize, etc.)
default_blue=[0 0.4470 0.7410];
default_red=[1 0 0];
default_darkgrey=[100/256,100/256,100/256];
default_orange=[231/256 111/256 81/256];
default_lightorange=[244/256 162/256 97/256];
default_yellow=[233/256 196/256 106/256];
color1=[29/256 46/256 126/256];
color2=[57/256 32/256 182/256];
color3=[133/256 203/256 207/256];

color4=[178/256 24/256 43/256];
color5=[213/256 96/256 76/256];
color6=[244/256 165/256 129/256];
Fontsize=18;

% read datafiles
%outputmonitor=csvread('OUTPUT_Current0.075chi0g0.048chi0q0.2D-11000gamma55.2004.csv');
%outputmonitor=csvread('OUTPUT_Current0.075chi0g0.048chi0q0.2D0gamma55.2004.csv');
%outputmonitor=csvread('OUTPUT_Current0.075chi0g0.048chi0q0.2D-800gamma0.csv');

%NSOA=csvread('N_SOACurrent0.075chi0g0.048chi0q0.2D-11000gamma55.2004.csv');
%NSA=csvread('N_SACurrent0.075chi0g0.048chi0q0.2D-11000gamma55.2004.csv');
%NISO=csvread('N_ISOCurrent0.075chi0g0.048chi0q0.2D-11000gamma55.2004.csv');
folder='C:\Users\Stijn\Documents\MATLAB\MLL_modeling\DATASWEEP_CONFIRM\';
outputmonitor=csvread([folder,'OUTPUT_MLLsim_Current_45X0g_0.07X0q_0.48beta2_1.3beta3_0.0042gamma_69.0005TPAandFCA_1RAMAN_1.csv']); 
NSOA=csvread([folder,'N_SOAMLLsim_Current_45X0g_0.07X0q_0.48beta2_1.3beta3_0.0042gamma_69.0005TPAandFCA_1RAMAN_1.csv']);
NSA=csvread([folder,'N_SAMLLsim_Current_45X0g_0.07X0q_0.48beta2_1.3beta3_0.0042gamma_69.0005TPAandFCA_1RAMAN_1.csv']);
NISO=csvread([folder,'N_ISOMLLsim_Current_45X0g_0.07X0q_0.48beta2_1.3beta3_0.0042gamma_69.0005TPAandFCA_1RAMAN_1.csv']);

dt=30*1e-15;
Nt=8.7*1e23;                    % transparency carrier density

% OUTPUT PULSE TRAIN
OUTPUTmonitor_selection=outputmonitor;
n_o_samples=max(size(OUTPUTmonitor_selection));
taxis=linspace(0,n_o_samples*dt,n_o_samples);
taxis=taxis';

figure;
plot(taxis.*1e9,abs(OUTPUTmonitor_selection).^2,'Color',default_blue);
xlabel('Time (ns)', 'FontSize', Fontsize);
ylabel('Envelope (W)', 'FontSize', Fontsize);
set(gca,'FontSize',Fontsize);
%xlim([0 max(taxis).*1e9]);
xlim([0 100]);
%ylim([0 0.5]);
set(gca,'box','off');
ax=gca();
tick_size=0.010;
ax.TickLength = [tick_size, tick_size]; % Make tick marks longer.
ax.LineWidth = 100*tick_size; % Make tick marks thicker.
set(gca, 'Layer', 'top')

% OUTPUT PULSE

t_capture=79.9*1e-9;
t_interval=10*1e-9;
OUTPUTmonitor_selection=outputmonitor(round(t_capture/dt):round((t_capture+t_interval)/dt));
[peaks,peakpositions]=findpeaks(abs(OUTPUTmonitor_selection), 'MinPeakHeight',max(abs(OUTPUTmonitor_selection))/10);
t_pulseinterval=20*1e-12;
t_peakposition=peakpositions(2);
OUTPUTmonitor_selection=OUTPUTmonitor_selection(t_peakposition-round(0.5*t_pulseinterval/dt):t_peakposition+round(0.5*t_pulseinterval/dt));
n_o_samples=max(size(OUTPUTmonitor_selection));
taxis=linspace(-0.5*n_o_samples*dt,0.5*n_o_samples*dt,n_o_samples);
taxis=taxis';

figure;
plot(taxis.*1e12,abs(OUTPUTmonitor_selection).^2,'Color',default_blue,'LineWidth',1.5);
xlabel('Eigentime (ps)', 'FontSize', Fontsize);
ylabel('Envelope (W)', 'FontSize', Fontsize);
set(gca,'FontSize',Fontsize);
ylim([0 0.6]);
xlim([-10 10]);
set(gca,'box','off');
ax=gca();
tick_size=0.010;
ax.TickLength = [tick_size, tick_size]; % Make tick marks longer.
ax.LineWidth = 100*tick_size; % Make tick marks thicker.
set(gca, 'Layer', 'top')


% CARRIER DENSITY

NSOA_selection=NSOA(round(t_capture/dt):round((t_capture+t_interval)/dt));
NSA_selection=NSA(round(t_capture/dt):round((t_capture+t_interval)/dt));
NISO_selection=NISO(round(t_capture/dt):round((t_capture+t_interval)/dt));

t_pulseinterval=100*1e-12;
NSOA_selection=NSOA_selection(t_peakposition-round(0.5*t_pulseinterval/dt):t_peakposition+round(0.5*t_pulseinterval/dt));
NSA_selection=NSA_selection(t_peakposition-round(0.5*t_pulseinterval/dt):t_peakposition+round(0.5*t_pulseinterval/dt));
NISO_selection=NISO_selection(t_peakposition-round(0.5*t_pulseinterval/dt):t_peakposition+round(0.5*t_pulseinterval/dt));

n_o_samples=max(size(NSOA_selection));
taxis=linspace(-0.5*n_o_samples*dt,0.5*n_o_samples*dt,n_o_samples);
taxis=taxis';

% Time domain trace of normalized carrier densities
figure;
plot(taxis.*1e12, NSOA_selection,'k','LineWidth',1.5);
hold on;
plot(taxis.*1e12, NSA_selection,'k--','LineWidth',1.5);
plot(taxis.*1e12, NISO_selection,'k-.','LineWidth',1.5);
xlabel('Eigentime (ps)', 'FontSize', Fontsize);
xlim([-0.5*t_pulseinterval*1e12 0.5*t_pulseinterval*1e12])
ylabel('Norm. Carrier Density', 'FontSize', Fontsize);
set(gca,'FontSize',Fontsize);
%grid on;
set(gca,'box','off');
legend('SOA','SA','ISO');
legend boxoff
set(gca,'box','off');
ax=gca();
tick_size=0.010;
ax.TickLength = [tick_size, tick_size]; % Make tick marks longer.
ax.LineWidth = 100*tick_size; % Make tick marks thicker.
set(gca, 'Layer', 'top')

% SPECTRUM
lambda0=1.57*1e-6;
[ validdata,lambdas, spectrum, lambdas_envelope, envelope, lambdas_phase, phases  ] = getspectrum( outputmonitor, lambda0, dt );

freq_center_pos=find(spectrum==max(spectrum));
lambda_center=lambdas(freq_center_pos);
figure;
plot(lambdas,spectrum,'Color',default_darkgrey);
%plot(lambdas,spectrum,'Color',color4);;
xlabel('Wavelength (nm)');
ylabel('Power (10 dB per div)');
ylim([-125 -85])
set(gca, 'YTick', [-120 -110 -100 -90]);
set(gca,'Yticklabel',[]) 
set(gca,'FontSize',Fontsize);
set(gca,'box','off');
xlim([lambda_center-5 lambda_center+5]);
set(gca,'box','off');
ax=gca();
tick_size=0.010;
ax.TickLength = [tick_size, tick_size]; % Make tick marks longer.
ax.LineWidth = 100*tick_size; % Make tick marks thicker.
set(gca, 'Layer', 'top')
xlim([1565 1572.2])
set(gca, 'XTick', [1566 1568 1570 1572]);

Fontsize=30;
figure;
plot(lambdas,spectrum,'Color',default_darkgrey);
%plot(lambdas,spectrum,'Color',color4);;
%xlabel('Wavelength (nm)');
%ylabel('3 dB per div');
ylim([-120 -90])
set(gca, 'YTick', [-120 -110 -100 -90]);
set(gca,'Yticklabel',[]) 
set(gca,'FontSize',Fontsize);
set(gca,'box','off');
xlim([lambda_center-5 lambda_center+5]);
set(gca,'box','off');
ax=gca();
tick_size=0.010;
ax.TickLength = [tick_size, tick_size]; % Make tick marks longer.
ax.LineWidth = 100*tick_size; % Make tick marks thicker.
set(gca, 'Layer', 'top')
xlim([1567-0.25 1567+0.25])
set(gca, 'XTick', [1567]);
set(gca,'FontSize',Fontsize);




%{
OUTPUTmonitor_selection=outputmonitor(round(150.6*1e-9/dt):max(size(outputmonitor)));
%OUTPUTmonitor_selection=outputmonitor(round(100*1e-9/dt):max(size(outputmonitor)));
OUTPUTmonitor_selection=OUTPUTmonitor_selection+0.01.*( randn(length(OUTPUTmonitor_selection),1) + 1i*randn(length(OUTPUTmonitor_selection),1)   );
freq_window=2;
c=3*1e8;

V=linspace(-0.5*2*pi/dt,0.5*2*pi/dt,length(OUTPUTmonitor_selection));
W=V+2*pi*c/lambda0;
f=(W./(2*pi)).*1e-12;
spectrum=(dt/length(V))*abs(fftshift(fft(OUTPUTmonitor_selection))).^2;

envelope_n_o_points=1000;
spectrum_envelope=zeros(1,envelope_n_o_points);
f_envelope=zeros(1,envelope_n_o_points);
n_o_samp=floor(length(spectrum)/envelope_n_o_points);
for i=1:envelope_n_o_points
    if i==length(spectrum)
        [temp_max,pos_max]=max(10*log10(spectrum(1+(i-1)*n_o_samp:length(spectrum))/1e-3));
    else
        [temp_max,pos_max]=max(10*log10(spectrum(1+(i-1)*n_o_samp:i*n_o_samp)/1e-3));
    end
    spectrum_envelope(i)=temp_max;
    f_envelope(i)=f(pos_max+(i-1)*n_o_samp);
end
freq_center_pos=find(spectrum==max(spectrum));
freq_center=f(freq_center_pos);
figure;
plot((c./f).*1e-3,10*log10(spectrum/1e-3),'Color',default_darkgrey);
hold on;
plot((c./f_envelope).*1e-3,spectrum_envelope,'r');
xlabel('Wavelength (nm)');
ylabel('Power (10 dB per div)');
ylim([-120 -85])
set(gca, 'YTick', [-120 -110 -100 -90]);
set(gca,'Yticklabel',[]) 
set(gca,'FontSize',Fontsize);
set(gca,'box','off');
xlim([1e-3*c/(freq_center+freq_window/2) 1e-3*c/(freq_center-freq_window/2)]);
set(gca,'box','off');
ax=gca();
tick_size=0.010;
ax.TickLength = [tick_size, tick_size]; % Make tick marks longer.
ax.LineWidth = 100*tick_size; % Make tick marks thicker.
set(gca, 'Layer', 'top')
%}
