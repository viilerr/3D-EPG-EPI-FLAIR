clear
clc
close all

% This is the file to run for the project.
% EPG_FLAIR.m is the function used by this test script.

flipAngle = deg2rad(90);
ESP = 0.8; % EPI echo spacing in ms
nPE = 64; % ky lines in the EPI readout
nPartitions = 32; % kz partitions in the 3D volume
effectiveEcho = floor(nPE/2) + 1;
TI = 3000; % close to CSF null: 4300*log(2) = 2980 ms

tissue(1) = struct('name','White matter','short','WM','T1',2933,'T2',119, ...
'T2star',32,'color',[0.000 0.278 0.671]);
tissue(2) = struct('name','Grey matter','short','GM','T1',2653,'T2',99, ...
'T2star',42,'color',[0.800 0.157 0.157]);
tissue(3) = struct('name','CSF','short','CSF','T1',4300,'T2',1200, ...
'T2star',450,'color',[0.050 0.050 0.050]);

signals = struct();
for ii = 1:numel(tissue)
[signals(ii).F0,signals(ii).K3D,signals(ii).Zn,signals(ii).F,signals(ii).info] = ...
EPG_FLAIR(flipAngle,ESP,tissue(ii).T1,tissue(ii).T2,TI, ...
'sequence','epi3d', ...
'T2star',tissue(ii).T2star, ...
'nPE',nPE, ...
'nPartitions',nPartitions, ...
'effectiveEcho',effectiveEcho);
end

echoNumber = signals(1).info.echoNumber;
effectiveTE = signals(1).info.effectiveTE;
ky = signals(1).info.kyPositions;
kz = signals(1).info.kzPositions;
[KY,KZ] = ndgrid(ky/max(abs(ky)),kz/max(abs(kz)));
centralKspaceWeight = exp(-3.5*(KY.^2 + KZ.^2));

wmEff = abs(signals(1).F0(effectiveEcho));
gmEff = abs(signals(2).F0(effectiveEcho));
csfEff = abs(signals(3).F0(effectiveEcho));
contrastWMGM = abs(wmEff-gmEff)/(wmEff+gmEff);
csfRatio = csfEff/max([wmEff gmEff]);

fig = figure('Name','3D FLAIR EPI test','Color','w','Units','normalized', ...
'Position',[0.08 0.08 0.82 0.82]);
tabs = uitabgroup(fig);

% 1. FLAIR inversion recovery
tab = uitab(tabs,'Title','1 FLAIR recovery');
axes('Parent',tab)
t = 0:1:(TI+800);
hold on
for ii = 1:numel(tissue)
Mz = 1 - 2*exp(-t/tissue(ii).T1);
plot(t,Mz,'LineWidth',2,'Color',tissue(ii).color)
end
xline(TI,'--','TI','LineWidth',1.4,'Color',[0.05 0.45 0.20])
yline(0,':','Color',[0.25 0.25 0.25])
xlabel('Time after inversion (ms)')
ylabel('M_z')
title('1. FLAIR recovery')
legend({tissue.short},'Location','southeast')
grid on
box on

% 2. EPI echo train
tab = uitab(tabs,'Title','2 EPI echo train');
axes('Parent',tab)
hold on
for ii = 1:numel(tissue)
plot(echoNumber,abs(signals(ii).F0),'LineWidth',2,'Color',tissue(ii).color)
end
xline(effectiveEcho,'--','k_0','LineWidth',1.4,'Color',[0.10 0.10 0.10])
xlabel('EPI echo / ky line')
ylabel('|F_0|')
title(sprintf('2. EPI echo train, TE = %.1f ms',effectiveTE))
legend({tissue.short},'Location','northeast')
grid on
box on

% 3. Tissue signal at the k-space centre
tab = uitab(tabs,'Title','3 Effective echo');
axes('Parent',tab)
b = bar([wmEff gmEff csfEff],'FaceColor','flat');
for ii = 1:numel(tissue)
b.CData(ii,:) = tissue(ii).color;
end
set(gca,'XTickLabel',{tissue.short})
ylabel('|F_0| at k-space centre')
title(sprintf('3. Contrast %.3f, CSF ratio %.3f',contrastWMGM,csfRatio))
grid on
box on

% 4. 3D ky-kz weighting
tab = uitab(tabs,'Title','4 3D k-space');
axes('Parent',tab)
imagesc(kz,ky,centralKspaceWeight)
axis image
colorbar
colormap(gca,'parula')
xlabel('k_z partition')
ylabel('k_y line')
title('4. 3D k-space weighting: centre bright, outer k-space darker')
grid on
box on

% 5. Point-spread function from EPI T2* decay
tab = uitab(tabs,'Title','5 PSF');
axes('Parent',tab)
hold on
for ii = 1:numel(tissue)
K = abs(signals(ii).F0);
K = K/max(K);
psf = abs(fftshift(ifft(ifftshift(K))));
psf = psf/max(psf);
plot(linspace(-0.5,0.5,numel(psf)),psf,'LineWidth',2, ...
'Color',tissue(ii).color)
end
xlabel('Normalised phase position')
ylabel('Normalised PSF')
title('5. EPI PSF')
legend({tissue.short},'Location','northeast')
grid on
box on

% 6. TI optimisation
tab = uitab(tabs,'Title','6 TI optimisation');
axes('Parent',tab)
TI_values = 1800:25:4200;
WM = zeros(size(TI_values));
GM = zeros(size(TI_values));
CSF = zeros(size(TI_values));
contrast = zeros(size(TI_values));

for ii = 1:numel(TI_values)
TItest = TI_values(ii);
[F0,~,~,~] = EPG_FLAIR(flipAngle,ESP,tissue(1).T1,tissue(1).T2,TItest, ...
'sequence','epi3d','T2star',tissue(1).T2star, ...
'nPE',nPE,'nPartitions',nPartitions,'effectiveEcho',effectiveEcho);
WM(ii) = abs(F0(effectiveEcho));

[F0,~,~,~] = EPG_FLAIR(flipAngle,ESP,tissue(2).T1,tissue(2).T2,TItest, ...
'sequence','epi3d','T2star',tissue(2).T2star, ...
'nPE',nPE,'nPartitions',nPartitions,'effectiveEcho',effectiveEcho);
GM(ii) = abs(F0(effectiveEcho));

[F0,~,~,~] = EPG_FLAIR(flipAngle,ESP,tissue(3).T1,tissue(3).T2,TItest, ...
'sequence','epi3d','T2star',tissue(3).T2star, ...
'nPE',nPE,'nPartitions',nPartitions,'effectiveEcho',effectiveEcho);
CSF(ii) = abs(F0(effectiveEcho));
contrast(ii) = abs(WM(ii)-GM(ii))/(WM(ii)+GM(ii));
end

[~,csfMinIdx] = min(CSF);
bestNullTI = TI_values(csfMinIdx);

purple = [0.45 0.10 0.60];
black = [0.05 0.05 0.05];
green = [0.05 0.45 0.20];

yyaxis left
hContrast = plot(TI_values,contrast,'LineWidth',2.2,'Color',purple);
ylabel('WM-GM contrast')
ax = gca;
ax.YAxis(1).Color = purple;

yyaxis right
hCSF = plot(TI_values,CSF,'LineWidth',2.2,'Color',black);
ylabel('CSF signal')
ax.YAxis(2).Color = black;

hNull = xline(bestNullTI,'--','minimum CSF','LineWidth',1.4,'Color',green);
xlabel('TI (ms)')
title(sprintf('6. TI optimisation, best TI = %.0f ms',bestNullTI))
legend([hContrast hCSF hNull],{'Purple: WM-GM contrast', ...
'Black: CSF signal','Green: best CSF null'},'Location','best')
grid on
box on

