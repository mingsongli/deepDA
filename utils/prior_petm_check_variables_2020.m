% check ensenmbe PETM015

% ensemble directory
%ens_dir = 'D:\cGENIE\ML.petm\ML.petm014\';
ens_dir = 'D:\cGENIE\ML.petm\ML.petm015\';
% working directory
wrk_dir = pwd;
%
int_dir = 'biogem';

var = 'atm_pCO2';
unit = ' (ppm)';

%var = 'atm_temp';
%unit = ' (degree C)';


biogemseries = ['biogem_series_',var,'.res'];

% cd ens. dir and read list
cd(ens_dir);
foldnames = dir;

% number of folders within
foldn = size(foldnames);
outmat = [];
y=[];
for i = 3 : foldn
    fname = foldnames(i).name;
    % id
    for j = 1:6
        fulldir = fullfile(ens_dir, fname,int_dir, biogemseries);
        int_var = load(fulldir);
        if strcmp(var, 'atm_pCO2')
            y(1:length(int_var(:,end)),i-1) = int_var(:,end) * 1E6;
        elseif strcmp(var, 'atm_temp')
            y(1:length(int_var(:,end)),i-1) = int_var(:,end);
        end
    end
end
y(:,1) = int_var(:,1);

for i = 2:6:61
    figure('Renderer', 'painters', 'Position', [50 50 1000 800])
    for j = 0:5
        if j > 1
            subplot(3,2,j-1)
        else
            subplot(3,2,j+5)
        end
        plot(y(:,1), y(:,i+j))
        title(foldnames(i+j+1).name)
        xlabel('Time (yr)')
        ylabel([var, unit], 'Interpreter', 'none')
    end
    %print(gcf,[foldnames(i+j+1).name(1:32),'.png'],'-dpng','-r200');
    fnamen = [foldnames(i+j+1).name(1:32),'.',var,'.fig'];
    set(gcf,'Name',fnamen,'NumberTitle','off');
    %saveas(gcf,fnamen)
end