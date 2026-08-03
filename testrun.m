clear                                   
clc                                    
close all                        

addpath(genpath('EPGX_functions'))      


theta = deg2rad([90 repmat(180,1,32)]); % Create RF pulse train:
                                        % 90° excitation
                                        % followed by thirty-two 180° refocusing pulses
ESP = 10;                          
T1_WM = 2933;                         
T2_WM = 119;                          
T1_GM = 2653;                         
T2_GM = 99;                            
T1_CSF = 4300;                      
T2_CSF = 1200;                         

TI = 3000;  % TInull=T1ln2 = 4300*0.693 = 2980 ms                           

[F0_WM,Fn_WM,Zn_WM,F_WM] = EPG_FLAIR(theta,ESP,T1_WM,T2_WM,TI);
[F0_GM,Fn_GM,Zn_GM,F_GM] = EPG_FLAIR(theta,ESP,T1_GM,T2_GM,TI);
[F0_CSF,Fn_CSF,Zn_CSF,F_CSF] = EPG_FLAIR(theta,ESP,T1_CSF,T2_CSF,TI);

ETL = length(F0_WM);       %number of echoes in the echo train
echo_number = 1:ETL;       %echo indices
TE = echo_number * ESP;    %echo time of each echo (ms)
effectiveEcho = 16;        %echo filling centre of k-space
effectiveTE = TE(effectiveEcho);

fprintf('Effective TE = %.1f ms\n',effectiveTE)
%plot echo trains
% expect WM highest, GM middle, CSF lowest, after every echo decay
figure
plot(abs(F0_WM),'b','LineWidth',2)      
hold on
plot(abs(F0_GM),'r','LineWidth',2)      
plot(abs(F0_CSF),'k','LineWidth',2)     
xlabel('Echo Number')  % 32 refocuing pulses so 32 echoes                 
ylabel('|F_0|')      % F0 is transverse magnetisation of zero gradient dephasing, physically this is exactly the magnetisation that induces a voltage in the MRI receiver coil so every point in the graph is how much MRI signal wouls scanner receive at this echo
title('EPG FLAIR Echo Train at 7T')
legend('White Matter','Grey Matter','CSF')
grid on


%how much longitudinal magnetisation has recovered by the time excitation
%pulse ocurs? this is beforer imaging begins Mxy
t = 0:1:(TI+500); %added +500 because couldnt see green line
Mz_WM = 1 - 2*exp(-t/T1_WM);
Mz_GM = 1 - 2*exp(-t/T1_GM);
Mz_CSF = 1 - 2*exp(-t/T1_CSF);

figure
plot(t,Mz_WM,'b','LineWidth',2)
hold on
plot(t,Mz_GM,'r','LineWidth',2)
plot(t,Mz_CSF,'k','LineWidth',2)
xline(TI,'--g', 'LineWidth', 2)
xlabel('Time after inversion (ms)')
ylabel('M_z')
title('FLAIR Inversion Recovery')
legend('White Matter','Grey Matter','CSF')
grid on



%transverse decay after the 90 deg excitation
t2 = 0:1:250;
Mxy_WM = abs(Mz_WM(end))*exp(-t2/T2_WM);
Mxy_GM = abs(Mz_GM(end))*exp(-t2/T2_GM);
Mxy_CSF = abs(Mz_CSF(end))*exp(-t2/T2_CSF);

figure
plot(t2,Mxy_WM,'b','LineWidth',2)
hold on
plot(t2,Mxy_GM,'r','LineWidth',2)
plot(t2,Mxy_CSF,'k','LineWidth',2)
xlabel('Time after excitation (ms)')
ylabel('M_{xy}')
title('Transverse Magnetisation after Excitation')
legend('White Matter','Grey Matter','CSF')
grid on


%TSE k-space trajectory and signal weighting
% In a linear TSE acquisition:
% - Each echo fills one phase-encoding line of k-space
% - The centre of k-space controls image contrast
% - The effective echo determines which echo fills the centre
ETL = length(F0_WM); %number of phase encoding lines,assume one echo fills one k-space line
% Negative values = outer k-space
% Zero = centre of k-space
% Positive values = opposite outer k-space

kspace_position = -floor(ETL/2):ceil(ETL/2)-1;
echo_number = 1:ETL;
effectiveEcho = 16; %fills centre of k-soace 
% x-axis = position in k-space
% y-axis = echo that fills that position

figure
stem(kspace_position,effectiveEcho*ones(size(kspace_position)),'filled')
hold on
stem(kspace_position,effectiveEcho + ...
 (echo_number-effectiveEcho),...
    'filled')
xlabel('k-space position')
ylabel('Echo number')
title('Linear TSE k-space trajectory')
grid on

%assign echo signals into k-space, each echo amplitude becomes one k-space line value

kspace_WM = zeros(1,ETL);
kspace_GM = zeros(1,ETL);
kspace_CSF = zeros(1,ETL); %preallocate k-space arrays



for ii = 1:ETL

    line = ii; %current k space location
    kspace_WM(line)=abs(F0_WM(ii)); %filling k space with corresponding echo signal
    kspace_GM(line)=abs(F0_GM(ii));
    kspace_CSF(line)=abs(F0_CSF(ii));
end

%shift k-space so that centre is displayed in the middle
kspace_WM_shifted = fftshift(kspace_WM);
kspace_GM_shifted = fftshift(kspace_GM);
kspace_CSF_shifted = fftshift(kspace_CSF);

figure
plot(kspace_position,kspace_WM_shifted,...
    'b','LineWidth',2)
hold on
plot(kspace_position,kspace_GM_shifted,...
    'r','LineWidth',2)
plot(kspace_position,kspace_CSF_shifted,...
    'k','LineWidth',2)
xlabel('k-space position')
ylabel('Signal magnitude')
title('TSE k-space signal weighting')
legend('White Matter','Grey Matter','CSF')
grid on

%point spread function measures exactly how blurred a point becomes: late 
%echoes are weak, outer k space weak hence blurrier image

K_WM=kspace_WM/max(kspace_WM);
K_GM=kspace_GM/max(kspace_GM);
K_CSF=kspace_CSF/max(kspace_CSF);
PSF_WM=abs(fftshift(ifft(K_WM)));
PSF_GM=abs(fftshift(ifft(K_GM)));
PSF_CSF=abs(fftshift(ifft(K_CSF)));
%normalise
PSF_WM=PSF_WM/max(PSF_WM);
PSF_GM=PSF_GM/max(PSF_GM);
PSF_CSF=PSF_CSF/max(PSF_CSF);


figure
plot(PSF_WM,'b','LineWidth',2)
hold on
plot(PSF_GM,'r','LineWidth',2)
plot(PSF_CSF,'k','LineWidth',2)
xlabel('Spatial position')
ylabel('Normalised PSF')
title('TSE FLAIR Point Spread Function')
legend('White Matter','Grey Matter','CSF')
grid on

%TI optimisation

TI_values = 1000:50:3500;
WM = zeros(size(TI_values));
GM = zeros(size(TI_values));
CSF = zeros(size(TI_values));
contrast = zeros(size(TI_values));
effectiveEcho=16;

for ii = 1:length(TI_values)
TItest = TI_values(ii);
[F0,~,~,~] = EPG_FLAIR(theta,ESP,T1_WM,T2_WM,TItest);
WM(ii) = abs(F0(effectiveEcho));
[F0,~,~,~] = EPG_FLAIR(theta,ESP,T1_GM,T2_GM,TItest);
GM(ii) = abs(F0(effectiveEcho));
[F0,~,~,~] = EPG_FLAIR(theta,ESP,T1_CSF,T2_CSF,TItest);
CSF(ii) = abs(F0(effectiveEcho));
contrast(ii)=abs(WM(ii)-GM(ii))/(WM(ii)+GM(ii));
%cant call abd(EPG_FLAIR) directly as it retirns multiple outputs and abs
%would only work on first output vector, instead take first output and
%extract first echo
end

figure
yyaxis left
h1 = plot(TI_values,contrast,'m','LineWidth',2);
ylabel('WM-GM Contrast')
yyaxis right
h2 = plot(TI_values,CSF,'k','LineWidth',2);
ylabel('CSF Signal')
xlabel('Inversion Time (ms)')
title('FLAIR TI Optimisation')
legend([h1 h2],{'WM-GM Contrast','CSF Signal'})
grid on
