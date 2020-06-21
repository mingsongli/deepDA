%% This code is for PETM ensemble
%% PETM014: 
% test sensitivity, 2 myr run from cold start with
% with variable pCO2 and delta F 2x (2, 4 default, 6 some, 8, and 12)
%% PETM015: 
% test sensitivity, 2 myr run from cold start with
% with variable pCO2 and delta F 2x (2, 4 default, 6 some, 8, 10, and 12)
%% Steps
% 1. read folder name, get pCO2, delf2x, values
% 2. read interested states: such as pCO2, ALK, SLT, SST, Global temperature, pH, CaCO3
% 3. Save as a matrix
% 4. plot
%%
% ensemble directory
ens_dir = 'D:\cGENIE\ML.petm\ML.petm015\';

% CESM CAM5 ECS
cesm =[1,3.5; 3,6.6; 6,9.7];
% xPAL to ppm
cesm(:,1) = cesm(:,1) * 285;
% working directory
wrk_dir = pwd;
%
int_dir = 'biogem';
int_file = {'biogem_series_atm_pCO2.res',...    
    'biogem_series_ocn_ALK.res',...
    'biogem_series_atm_temp.res',...
    'biogem_series_ocn_temp.res',...
    'biogem_series_misc_surpH.res',...
    'biogem_series_sed_CaCO3.res'};
% cd ens. dir and read list
cd(ens_dir);
foldnames = dir;
% number of folders within
foldn = size(foldnames);
int_filen = length(int_file);
outmat = [];

ids  = 1009:1343;  % last 500 kyr

for i = 3 : foldn
    fname = foldnames(i).name;
    % id
    outmat(i-2, 1) = i-2;
    % outgas
    outmat(i-2, 2) = str2double(fname(30)) + str2double(fname(32))/10;
    delf2x_raw = fname(end-1:end);
    % delf2x * log(2) = Wm-2 radiative forcing
    if strcmp(delf2x_raw(1),'x')
        outmat(i-2, 3) = str2double(fname(end));
    else
        outmat(i-2, 3) = str2double(fname(end-1:end));
    end
    
    for j = 1:6
        int_file_j = int_file{j};
        fulldir = fullfile(ens_dir, fname,int_dir, int_file_j);
        int_var = load(fulldir);
        try
            % die exp will be skipped
            if strcmp(int_file_j, int_file{1})
                % pCO2
                outmat(i-2, 4) = max(int_var(ids,end)) * 1E6;
            end
            if strcmp(int_file_j, int_file{2})
                % ALK
                outmat(i-2, 5) = max(int_var(ids,end))*1000;
            end
            if strcmp(int_file_j, int_file{3})
                % SAT
                outmat(i-2, 6) = max(int_var(ids,end));
            end
            if strcmp(int_file_j, int_file{4})
                % SST
                outmat(i-2, 7) = max(int_var(ids,end));
                % global
                outmat(i-2, 8) = max(int_var(ids,2));
                % benthic
                outmat(i-2, 9) = max(int_var(ids,4));
            end
            if strcmp(int_file_j, int_file{5})
                % surface pH
                outmat(i-2, 10) = min(int_var(ids,end));
            end
            if strcmp(int_file_j, int_file{6})
                % CaCO3
                outmat(i-2, 11) = max(int_var(ids,end));
            end
        catch
        end
    end
end

%% plot log2 fit
figure('Renderer', 'painters', 'Position', [50 50 1000 800])
hold on
option = 1; % 1 = from 0.6; 2 = 0.8; 3 = 1.0
%idsplus = 3;  % 0 = x10;  1 = x12; 2 = x2; 3 = x4; 4 = x6; 5 = x8
slopei =[];idspi =[10,12,2,4,6,8];
for idsplus = 0:5
    idid = {'x10', 'x12', 'x2', 'x4', 'x6', 'x8'};
    if option == 1
        if idsplus == 1
            ids = 7:6:60;
        else
            ids = 1:6:60;
        end
    elseif option == 2
        ids = 13:6:60;
        idslen = length(ids);
    elseif option == 3
        ids = 19:6:60;
        idslen = length(ids);
    elseif option == 4
        ids = 25:6:60;
        idslen = length(ids);
    elseif option == 5
        ids = 31:6:60;
        idslen = length(ids);
    end
    
    ids = idsplus + ids;
    sf1 = fit(log2(outmat(ids,4)),outmat(ids,7),'poly1');
    plot(sf1,log2(outmat(ids,4)),outmat(ids,7));
    text(max(log2(outmat(ids,4)))+0.1,max(outmat(ids,7))+0.2,num2str(sf1.p1,'%2.2f'))
    slopei(idsplus+1,1) = idspi(idsplus+1);
    slopei(idsplus+1,2) = sf1.p1;
    slopei(idsplus+1,3) = sf1.p2;
end
hold off
xlabel('log_2(pCO_2) (ppm)');ylabel('SAT (degree C)');set(gcf,'color','white')
legend off
title(['log_2(pCO_2) vs. SAT & delf2x '])

% PLOT 
sf2 = fit(slopei(:,1),slopei(:,2),'poly1');
figure('Renderer', 'painters', 'Position', [350 50 1000 800])
y = slopei(:,2);
yCalc1 = slopei(:,1) * sf2.p1 + sf2.p2;
Rsq = 1 - sum((y - yCalc1).^2)/sum((y - mean(y)).^2);
plot(sf2,slopei(:,1),slopei(:,2));
legend off
title(['Climate sensitivity vs. ECS'])
xlabel('delf2x * ln(2) (W/m^2)');ylabel(['ECS (', char(176),'C)']);set(gcf,'color','white')
xlim([0, 14]); ylim([0, 10])
eq1 = ['ECS = ', num2str(sf2.p1),' * x + ', num2str(sf2.p2)];
text(6,8,eq1)
eq2 = ['R^2 = ',num2str(Rsq)];
text(6,7.5,eq2)

% fit CESM
figure('Renderer', 'painters', 'Position', [350 50 1000 800])
sf3 = fit(log2(cesm(:,1)), cesm(:,2),'poly1');
plot(sf3,log2(cesm(:,1)), cesm(:,2),'ro');
y = cesm(:,2);
yCalc1 = log2(cesm(:,1)) * sf3.p1 + sf3.p2;
Rsq = 1 - sum((y - yCalc1).^2)/sum((y - mean(y)).^2);
eq1 = ['ECS = ', num2str(sf3.p1),' * x + ', num2str(sf3.p2)]; text(9,8,eq1)
eq2 = ['R^2 = ',num2str(Rsq)]; text(9,7.5,eq2)
xlabel('log_2CO_2 (x1 PAL)'); ylabel(['ECS (', char(176),'C)']);title('iCESM + CAM5');set(gcf,'color','white');legend off
hold on; plot(log2(278),3,'bo');

% fit to 278 , 3 degree C
CESMi = cesm;
CESMi(1,:) = [278, 3.15];
%CESMi(1,:) = [278, 3.13];
%figure('Renderer', 'painters', 'Position', [350 50 1000 800])
sf4 = fit(log2(CESMi(:,1)), CESMi(:,2),'poly1');
plot(sf4,log2(CESMi(:,1)), CESMi(:,2));
y = CESMi(:,2);
yCalc1 = log2(CESMi(:,1)) * sf4.p1 + sf4.p2;
Rsq = 1 - sum((y - yCalc1).^2)/sum((y - mean(y)).^2);
eq1 = ['ECS = ', num2str(sf4.p1),' * x + ', num2str(sf4.p2)]; text(9,4,eq1)
eq2 = ['R^2 = ',num2str(Rsq)]; text(9,3.5,eq2)
xlabel('log_2CO_2 (x1 PAL)'); ylabel(['ECS (', char(176),'C)']);title('iCESM + CAM5 fit 278,3C');set(gcf,'color','white');legend off
hold on; plot(log2(278),3,'ro');hold off



%wm2 = delf2x*log(2);
co2 = 278; wm2 = (2.463661 * log2(co2) - 17.002343 + 0.201808) / 0.737505;
co2 = 278: 10: 5000; co2 = co2';
delf2x = 4.81937 * log2(co2) -32.86485;
figure; plot(log2(co2), delf2x);xlabel('log2 pco2'),ylabel('delf2x');title('pCO2 vs. delf2x');legend off; set(gcf,'color','white')

%   2     3       4   5   6    7     8      9      10  11
% outgas delf2x pco2 ALK SAT  SST  global benthic  pH  CaCO3

figure; scatter3(outmat(:,2),outmat(:,3),outmat(:,6),[],outmat(:,4),'filled'); title('SAT vs. outgas & Radiative forcing | color=pco2')
xlabel('outgas (x 3pal outgas)');ylabel('Climate forcing (W/m^2)');zlabel('SAT (degree C)');set(gcf,'color','white')
view(0,0),zlim([15,32])

figure; 
scatter3(outmat(:,2),outmat(:,6),outmat(:,11),[],outmat(:,3),'filled'); title('Outgas vs SAT vs. CaCO3 | color= W/m^2')
xlabel('Outgas (x 3pal outgas)'), ylabel('SAT (degree C)');zlabel('CaCO3');set(gcf,'color','white')
view(0,0)

figure; 
scatter3(log2(outmat(:,4)),outmat(:,6),outmat(:,5),[],outmat(:,3),'filled'); title('log2pco2 vs SAT vs. ALK | color= W/m^2')
xlabel('log2(pco2)'), ylabel('SAT (degree C)');zlabel('ALK');set(gcf,'color','white');%xlim([8 12.5]);zlim([1 4])
view(0,0)