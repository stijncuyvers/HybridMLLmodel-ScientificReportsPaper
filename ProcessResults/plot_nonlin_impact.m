clear all
close all

%color1=[255/255 221/255 0];
%color2=[30/255 100/255 200/255];
%color3=[0/255 0/255 0];

color1=[29/256 46/256 126/256];
color2=[57/256 32/256 182/256];
color3=[133/256 203/256 207/256];

color4=[178/256 24/256 43/256];
color5=[213/256 96/256 76/256];
color6=[244/256 165/256 129/256];
Fontsize=18;

folder='C:\Users\Stijn\Documents\MATLAB\MLL_modeling\DATASWEEP_CONFIRM\';
disp1=dlmread([folder,'OUTPUT_MLLsim_Current_45X0g_0.07X0q_0.48beta2_1.3beta3_0.0042gamma_0TPAandFCA_1RAMAN_0.csv']);
disp2=dlmread([folder,'OUTPUT_MLLsim_Current_45X0g_0.07X0q_0.48beta2_1.3beta3_0.0042gamma_69.0005TPAandFCA_1RAMAN_1.csv']);
disp3=dlmread([folder,'OUTPUT_MLLsim_Current_45X0g_0.07X0q_0.48beta2_1.3beta3_0.0042gamma_150TPAandFCA_1RAMAN_1.csv']);

[ selection1, fit, t ]= extractpulse( disp1 );
[ selection2, fit2, t2 ]= extractpulse( disp2 );
[ selection3, fit3, t3 ]= extractpulse( disp3 );
timespan=max(t)-min(t);
timeaxis=linspace(-timespan/2,timespan/2,length(t));
figure;
plot(timeaxis,selection1,'Color',color3,'LineWidth',1.5);
hold on;
plot(timeaxis,selection2,'-.','Color',color2,'LineWidth',1.5);
plot(timeaxis,selection3,'--','Color',color1,'LineWidth',1.5);
xlabel('Time (ps)','Fontsize',Fontsize);
xlim([-5 6]);
set(gca, 'XTick', [-5 0 5]);
set(gca,'Fontsize',Fontsize);
ylabel('Envelope (W)','Fontsize',Fontsize);
%grid on;
legend({'\gamma_{NL}=0 m^{-1}W^{-1}','\gamma_{NL}=69 m^{-1}W^{-1}','\gamma_{NL}=150 m^{-1}W^{-1}'},'Fontsize',Fontsize-2);
legend boxoff  
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
ylabel('Power (20 dB per div)');
legend('\gamma_{NL}=0 m^{-1}W^{-1}','\gamma_{NL}=69 m^{-1}W^{-1}','\gamma_{NL}=150 m^{-1}W^{-1}');
legend boxoff  
ylim([-120 -80])
set(gca, 'YTick', [-120 -100 -80]);
set(gca,'Yticklabel',[]) 
set(gca,'FontSize',Fontsize);
set(gca,'box','off');
xlim([lambda_center-7 lambda_center+7]);
set(gca,'box','off');
ax=gca();
tick_size=0.010;
ax.TickLength = [tick_size, tick_size]; % Make tick marks longer.
ax.LineWidth = 100*tick_size; % Make tick marks thicker.
set(gca, 'Layer', 'top')


