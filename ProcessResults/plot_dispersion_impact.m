%clear all
close all

color1=[29/256 46/256 126/256];
color2=[57/256 32/256 182/256];
color3=[133/256 203/256 207/256];

color4=[178/256 24/256 43/256];
color5=[213/256 96/256 76/256];
color6=[244/256 165/256 129/256];
Fontsize=18;

folder='C:\Users\Stijn\Documents\MATLAB\MLL_modeling\DATASWEEP_CONFIRM\';
disp1=dlmread([folder,'OUTPUT_MLLsim_Current_45X0g_0.07X0q_0.48beta2_1.3beta3_0.0042gamma_69.0005TPAandFCA_1RAMAN_1.csv']);
disp2=dlmread([folder,'OUTPUT_MLLsim_Current_45X0g_0.07X0q_0.48beta2_-4beta3_0.0042gamma_69.0005TPAandFCA_1RAMAN_1.csv']);
disp3=dlmread([folder,'OUTPUT_MLLsim_Current_45X0g_0.07X0q_0.48beta2_4beta3_0.0042gamma_69.0005TPAandFCA_1RAMAN_1.csv']);

[ selection1, fit, t ]= extractpulse( disp1 );
[ selection2, fit2, t2 ]= extractpulse( disp2 );
[ selection3, fit3, t3 ]= extractpulse( disp3 );

% manual location selection of pulse trace 3
taxiss=linspace(-10,10,(20*1e-12)/dt);
posss=round(195.633*1e-9/dt);
posssmin=round(posss-length(taxiss)/2+1);
posssmax=round(posss+length(taxiss)/2);

timespan=max(t)-min(t);
timeaxis=linspace(-timespan/2,timespan/2,length(t));
figure;
plot(timeaxis,selection1,'Color',color3,'LineWidth',1.5);
hold on;
%plot(t,fit,'r--');
plot(timeaxis,selection2,'-.','Color',color2,'LineWidth',1.5);
%plot(t,selection50,'.');
plot(timeaxis,selection3,'Color',color1,'LineWidth',1.5);
%plot(taxiss',abs(disp3(posssmin:posssmax)).^2,'--','Color',color1,'LineWidth',1.5);
xlabel('Time (ps)','Fontsize',Fontsize);
xlim([-5 6]);
set(gca, 'XTick', [-5 0 5]);
set(gca,'Fontsize',Fontsize);
ylabel('Envelope (W)','Fontsize',Fontsize);
%title('Output pulse');
%grid on;
legend({'\beta_{2}=1.3 ps^2/m','\beta_{2}=-4 ps^2/m','\beta_{2}=4 ps^2/m'},'Fontsize',Fontsize-2);
legend boxoff
%title('Impact of dispersion D (ps/nm/km)');  
set(gca,'box','off');
set(gca,'FontSize',Fontsize);
ax=gca();
tick_size=0.010;
ax.TickLength = [tick_size, tick_size]; % Make tick marks longer.
ax.LineWidth = 100*tick_size; % Make tick marks thicker.
set(gca, 'Layer', 'top')

lambda0=1.57e-6;
dt=30*1e-15;
[validdata,lambdas, spectrum, lambdas_envelope, envelope, lambdas_phase, phases] = getspectrum( disp1, lambda0, dt );
[validdata,lambdas2, spectrum2, lambdas_envelope2, envelope2, lambdas_phase2, phases2] = getspectrum( disp2, lambda0, dt );
[validdata3,lambdas3, spectrum3, lambdas_envelope3, envelope3, lambdas_phase3, phases3] = getspectrum( disp3, lambda0, dt );

freq_center_pos=find(spectrum==max(spectrum));
lambda_center=lambdas(freq_center_pos);
figure;
plot(lambdas_envelope,envelope,'Color',color4,'LineWidth',1.5);
hold on;
plot(lambdas_envelope2,envelope2,'Color',color5,'LineWidth',1.5);
plot(lambdas_envelope3,envelope3,'Color',color6,'LineWidth',1.5);
xlabel('Wavelength (nm)');
ylabel('Power (10 dB per div)');
ylim([-120 -85])
set(gca, 'YTick', [-120 -110 -100 -90]);
legend('\beta_{2}=1.3 ps^2/m','\beta_{2}=-4 ps^2/m','\beta_{2}=4 ps^2/m');
legend boxoff
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

freq_center_pos=find(spectrum==max(spectrum));
lambda_center=lambdas(freq_center_pos);
figure;
plot(lambdas,spectrum,'Color',color4,'LineWidth',1.5);
xlabel('Wavelength (nm)');
ylabel('Power (10 dB per div)');
ylim([-120 -85])
set(gca, 'YTick', [-120 -110 -100 -90]);
legend('\beta_{2}=1.3 ps^2/m','\beta_{2}=-4 ps^2/m','\beta_{2}=4 ps^2/m');
legend boxoff
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

freq_center_pos=find(spectrum==max(spectrum));
lambda_center=lambdas(freq_center_pos);
figure;
plot(lambdas2,spectrum2,'Color',color4,'LineWidth',1.5);
xlabel('Wavelength (nm)');
ylabel('Power (10 dB per div)');
ylim([-120 -85])
set(gca, 'YTick', [-120 -110 -100 -90]);
legend('\beta_{2}=1.3 ps^2/m','\beta_{2}=-4 ps^2/m','\beta_{2}=4 ps^2/m');
legend boxoff
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

freq_center_pos=find(spectrum==max(spectrum));
lambda_center=lambdas(freq_center_pos);
figure;
plot(lambdas3,spectrum3,'Color',color4,'LineWidth',1.5);
xlabel('Wavelength (nm)');
ylabel('Power (10 dB per div)');
ylim([-120 -85])
set(gca, 'YTick', [-120 -110 -100 -90]);
legend('\beta_{2}=1.3 ps^2/m','\beta_{2}=-4 ps^2/m','\beta_{2}=4 ps^2/m');
legend boxoff
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

