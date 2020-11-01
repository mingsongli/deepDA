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
%ens_dir = 'D:\cGENIE\ML.petm\ML.petm014\';
ens_dir = 'D:\cGENIE\ML.petm\ML.petm015\';

% working directory
wrk_dir = pwd;
%
int_dir = 'biogem';
int_file = {'biogem_series_atm_pCO2.res',...    
    'biogem_series_ocn_ALK.res',...
    'biogem_series_misc_SLT.res',...
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
        if strcmp(int_file_j, int_file{1})
            % pCO2
            outmat(i-2, 4) = int_var(end,end) * 1E6;
        end
        if strcmp(int_file_j, int_file{2})
            % ALK
            outmat(i-2, 5) = int_var(end,3)*1000;
        end
        if strcmp(int_file_j, int_file{3})
            % SLT
            outmat(i-2, 6) = int_var(end,end);
        end
        if strcmp(int_file_j, int_file{4})
            % SST
            outmat(i-2, 7) = int_var(end,end);
            % global
            outmat(i-2, 8) = int_var(end,2);
            % benthic
            outmat(i-2, 9) = int_var(end,4);
        end
        if strcmp(int_file_j, int_file{5})
            % surface pH
            outmat(i-2, 10) = int_var(end,end);
        end
        if strcmp(int_file_j, int_file{6})
            % CaCO3
            outmat(i-2, 11) = int_var(end,end);
        end
    end
end

%   2     3       4   5   6    7     8      9      10  11
% outgas delf2x pco2 ALK SLT  SST  global benthic  pH  CaCO3
% plot pco2
figure; scatter3(outmat(:,2),outmat(:,3),log2(outmat(:,4)),sqrt(outmat(:,4)),log2(outmat(:,4)),'filled'); title('log2pCO2 vs. outgas & Wm-2 | color=pco2')
xlabel('Outgas (x 3pal outgas)');ylabel('Radiative forcing (W/m^2)');zlabel('log2 pCO2 (ppm)');set(gcf,'color','white')

% plot SST
figure; scatter3(outmat(:,3),outmat(:,4),outmat(:,7),[],outmat(:,7),'filled'); title('SST vs. pCO2 & Wm-2 | color=SST')
xlabel('Radiative forcing (W/m^2)');ylabel('pCO2 (ppm)');zlabel('SST (degree C)')
ylim([0 5000])
set(gcf,'color','white')

figure; 
scatter3(log2(outmat(:,4)),outmat(:,7),outmat(:,5),[],outmat(:,3),'filled'); title('log2pco2 vs SST vs. ALK | color= W/m^2')
xlabel('log2(pco2)'), ylabel('SST (degree C)');zlabel('ALK');set(gcf,'color','white');xlim([8 12.5]);zlim([1 4])
view(0,0)

% plot SST vs. outgas & Wm-2
figure; scatter3(outmat(:,2),outmat(:,3),outmat(:,7),[],outmat(:,4),'filled'); title('SST vs. outgas & Radiative forcing | color=pco2')
xlabel('outgas (x 3pal outgas)');ylabel('Radiative forcing (W/m^2)');zlabel('SST (degree C)');set(gcf,'color','white')

% CaCO3 SST pCO2
figure; 
scatter3(log2(outmat(:,4)),outmat(:,7),outmat(:,11),[],outmat(:,3),'filled'); title('pco2 vs SST vs. CaCO3 | color= W/m^2')
xlabel('log2(pCO2)');ylabel('SST (degree C)');zlabel('CaCO3');set(gcf,'color','white')
view(90,0)

sf2 = fit([log2(outmat(:,4)),outmat(:,7)],outmat(:,11),'poly23');
figure;plot(sf2,[log2(outmat(:,4)),outmat(:,7)],outmat(:,11));
title('pco2 vs SST vs. CaCO3 | color= W/m^2')
xlabel('log2(pCO2)');ylabel('SST (degree C)');zlabel('CaCO3');set(gcf,'color','white')
zlim([0 100])

figure; 
scatter3(outmat(:,2),outmat(:,7),outmat(:,11),[],outmat(:,3),'filled'); title('Outgas vs SST vs. CaCO3 | color= W/m^2')
xlabel('Outgas (x 3pal outgas)'), ylabel('SST (degree C)');zlabel('CaCO3');set(gcf,'color','white')

figure; 
scatter3(outmat(:,2),outmat(:,7),outmat(:,10),[],outmat(:,3),'filled'); title('Outgas vs SST vs. pH | color= W/m^2')
xlabel('Outgas (x 3pal outgas)'), ylabel('SST (degree C)');zlabel('pH');set(gcf,'color','white')

figure; 
scatter3(log2(outmat(:,4)),outmat(:,7),outmat(:,10),[],outmat(:,3),'filled'); title('log2pco2 vs SST vs. pH | color= W/m^2')
xlabel('log2(pco2)'), ylabel('SST (degree C)');zlabel('pH');set(gcf,'color','white')
view(90,0)
view(0,0)
%% plot linear fit
figure;
hold on
for idsplus = 0:5
%idsplus = 3;  % 0 = x10;  1 = x12; 2 = x2; 3 = x4; 4 = x6; 5 = x8
ids = 1:6:60;
idid = {'x10', 'x12', 'x2', 'x4', 'x6', 'x8'};
ids = idsplus + ids;
sf1 = fit(log(outmat(ids,4)),outmat(ids,7),'poly1');
plot(sf1,log(outmat(ids,4)),outmat(ids,7));
end
hold off
title(['log(pCO2) vs. log(pCO2)']);set(gcf,'color','white')
xlabel('log2(pCO2)');ylabel('SST (degree C)');
%title(['log(pCO2) vs. log(pCO2) delf ', idid{idsplus+1}]);set(gcf,'color','white')

%% plot log2 fit
figure;
hold on
slopei =[];
for idsplus = 0:5
    %idsplus = 3;  % 0 = x10;  1 = x12; 2 = x2; 3 = x4; 4 = x6; 5 = x8
    ids = 1:6:60;
    idid = {'x10', 'x12', 'x2', 'x4', 'x6', 'x8'};
    ids = idsplus + ids;
    sf1 = fit(log2(outmat(ids,4)),outmat(ids,7),'poly1');
    plot(sf1,log2(outmat(ids,4)),outmat(ids,7));
    text(max(log2(outmat(ids,4)))+.2,max(outmat(ids,7))+.2,num2str(sf1.p1))
    slopei(idsplus+1,1) = sf1.p1;
    slopei(idsplus+1,2) = sf1.p2;
end
hold off
title(['log2(pCO2) vs. log2(pCO2) delf ', idid{idsplus+1}])
xlabel('log2(pCO2) ppm');ylabel('SST (degree C)');set(gcf,'color','white')


%% unused
% sflogco2 = fit([outmat(:,3),log(outmat(:,4))],outmat(:,7),'poly23');
% figure;plot(sflogco2,[outmat(:,3),log(outmat(:,4))],outmat(:,7));title('SST vs. ln(pCO2) & Wm-2')
% xlabel('Wm-2');ylabel('ln(pCO2)');zlabel('SST (degree C)')
% zlim([20 40])
% 
% sflog2co2 = fit([outmat(:,3),log2(outmat(:,4))],outmat(:,7),'poly23');
% figure;plot(sflog2co2,[outmat(:,3),log2(outmat(:,4))],outmat(:,7));title('SST vs. log2(pCO2) & Wm-2')
% xlabel('Wm-2');ylabel('log2(pCO2)');zlabel('SST (degree C)')
% zlim([20 40])