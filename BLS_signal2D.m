clearvars;
clc;

% constants
hbar=6.626e-34/(2*pi); %kg m s^-1
ee=1.602e-19; %C
me=9.109e-31; %kg
c=2.998e8; %m s^-1
KB=1.38e-23;
Eps0=8.85e-12;
mu0=4*pi*1e-7;

Color=[18/255,103/255,221/255;247/255,66/255,66/255;117/255,117/255,117/255];

%% Loading data from FDTD - Data1: illumination path, Data2: detection path

FileName='MieStripes_X'; % polarization in x
[Data,f]=DataLoader(FileName);

FileName='MieStripes_RefX'; % reference (without structures), adding the reference to Data
tmp=DataLoader(FileName);
Data=cat(1,Data,tmp);

FileName='MieStripes_Y'; % polarization in y
Data2=DataLoader(FileName);

FileName='MieStripes_RefY'; % reference (without structures), adding the reference to Data2
tmp=DataLoader(FileName);
Data2=cat(1,Data2,tmp);

M=length(Data); % getting the length of Data, equals length of Data2

%% Calculating magnetization data

lambda=c/f;
k0=2*pi/lambda;

fm = 401; % number of frequency points

kmax = 300; % maximal wavenumber in both x and y
dk = 2; % wavenumber step in both x and y
    kx = linspace(-kmax,kmax,2*kmax/dk+1)*1e6;
    ky = linspace(-kmax,kmax,2*kmax/dk+1)*1e6;

mx = zeros(fm,length(kx),length(ky));
my = zeros(fm,length(kx),length(ky));
mz = zeros(fm,length(kx),length(ky));


% function [ff, f00fkx] = SpinWaveGreen(kx, Bext, d, N, Nf, n, mu),
% material parameters have to be adjusted directly in the function
Bext = 0.05; % external field
d = 25.99885e-9; % NiFe thickness
mu = -1e12; % chemical potencial of magnons (in frequency)
% mu = 0; % chemical potencial of magnons (in frequency)

[ffn0, f00fkxn0] = SpinWaveGreen(kx+0.0001, Bext, d, length(kx), fm, 0, mu); 
[ffn1, f00fkxn1] = SpinWaveGreen(kx+0.0001, Bext, d, length(kx), fm, 1, mu);

ffE = linspace(min(ffn0(:)), max(ffn1(:)), fm); % calculated freq range (preparing for interpolation)

f00fkxn0Interp = interp1(ffn0, f00fkxn0, ffE, 'linear', 0);
f00fkxn1Interp = interp1(ffn1, f00fkxn1, ffE, 'linear', 0);

f00kxSumInterp = f00fkxn0Interp + f00fkxn1Interp;
f00kxSumInterp(isnan(f00kxSumInterp))=0;

ii = 0;
for ii=1:length(ffE)

    %Calculation of polarization
    mx(ii,:,:) = mx(ii,:,:) + (f00kxSumInterp(ii,:,:));
    my(ii,:,:) = my(ii,:,:) + zeros(size(my(ii,:,:)));
    mz(ii,:,:) = mz(ii,:,:) - 1j*(f00kxSumInterp(ii,:,:));

end

mx = flip(permute(mx,[3,2,1]),1); % flip not necessary in this instance
my = flip(permute(my,[3,2,1]),1);
mz = flip(permute(mz,[3,2,1]),1);


Nfm=length(ffE);

chiNiFe={zeros(size(mz)),1i*mz,-1i*my;-1i*mz,zeros(size(mz)),1i*mx;1i*my,-1i*mx,zeros(size(mz))}; % dynamic magnetic susceptibility tensor

%% Plot magnetization

Fig=figure(1); % plot of the magnon dispersion for a chosen magnon frequency
Fig.OuterPosition=[200 100 1500 450];

Label={{'$m_{x} (k_{x},k_{y})$'};{'$m_{y} (k_{x},k_{y})$'};{'$m_{z} (k_{x},k_{y})$'}};
fmChosen=27; % chosen frequency
[~,fmIndex]=min(abs(ffE(:)-fmChosen),[],1);
PlotData={mx(:,:,fmIndex);my(:,:,fmIndex);mz(:,:,fmIndex)};

for kk=1:3
    subplot(1,3,kk);
    
    tmp=abs(PlotData{kk});
    imagesc(kx*1e-6,ky*1e-6,tmp.');

    set(gca,'FontSize',14,'FontName','Times New Roman');
    set(gca,'YDir','normal');
    set(gca,'XGrid','on');
    set(gca,'YGrid','on');
    set(gca,'XMinorTick','on');
    set(gca,'YMinorTick','on');
    set(gca,'TickDir','both');
    xlabel('$k_{x}\,\mathrm{(1/\mu{m})}$','Interpreter','latex','FontSize',14,'FontName','Times New Roman');
    ylabel('$k_{y}\,\mathrm{(1/\mu{m})}$','Interpreter','latex','FontSize',14,'FontName','Times New Roman','VerticalAlignment','bottom');
    colormap(inferno(2048));

    hh=colorbar;
    set(hh,'Box','off','FontSize',14,'FontName','Times New Roman','TickDirection','out');
    hh.Label.String = Label{kk};
    hh.Label.Interpreter='latex';
    hh.Label.VerticalAlignment='top';
end
%% Overlap integral stemming from reciprocity theorem

WindowSize=4e-6; % size of the window in which the overlap integral is evaluated
Ei = cell(M,1);
Ej = cell(M,1);

sigmaSW=zeros(Nfm,M);
for mm=1:M

    xi=Data{mm}.xi; % size of simulation area in x
    yi=Data{mm}.yi; % size of simulation area in y

    x=xi(abs(xi)<WindowSize/2); % crop of simulation area in x
    y=yi(abs(yi)<WindowSize/2); % crop of simulation area in y
    Nx=length(x); % number of steps in cropped simulation area in x
    Ny=length(y); % number of steps in cropped simulation area in y
    [X,Y]=ndgrid(x,y); % grid 

    % calculating each step
    dx=x(2:end)-x(1:end-1); %dx not equidistant, from second to second from last
    dx=cat(1,0,dx,0); % zero in the beginning and end
    dx=(dx(1:end-1)+dx(2:end))/2; % (without last zero + without first zero)/2

    dy=y(2:end)-y(1:end-1);
    dy=cat(1,0,dy,0);
    dy=(dy(1:end-1)+dy(2:end))/2;

    dS=dx*dy'; %areas of the surface differentials
    
    Ei{mm}=cell(3,1);
    Ej{mm}=cell(3,1);
    for ii=1:3 % loads E-field data from Lumerical
        Ei{mm}{ii}=squeeze(Data{mm}.Ei(abs(xi)<WindowSize/2,abs(yi)<WindowSize/2,ii)); % polarization x
        Ej{mm}{ii}=squeeze(Data2{mm}.Ei(abs(xi)<WindowSize/2,abs(yi)<WindowSize/2,ii)); % polarization y
    end

    EjEi=cell(3,3);
    for uu=1:3
        for vv=1:3
            EjEi{uu,vv}=Ej{mm}{uu}.*Ei{mm}{vv}; % calculating matrix of different polarization combinations (Ej, Ei outer product calculation)
        end
    end

    Nkm=length(kx)^2; %total number of evaluated magnon wavevectors
    [Kxm,Kym]=ndgrid(kx,ky);
    qmEiEj{mm}=repmat({zeros(size(Kxm))},[3,3]);
    % multiWaitbar( 'Reciprocity theorem', 0, 'Color', 'g', 'CanCancel', 'on' );
    fprintf('Reciprocity theorem: 0%% complete\n');
    for ii=1:Nkm
        % multiWaitbar( 'Reciprocity theorem', ii / Nkm );
        if mod(ii, round(Nkm/20)) == 0
        fprintf('Reciprocity theorem: %d%% complete\n', round(ii/Nkm*100));
        end

        WaveFactor=dS.*exp(1i*(Kxm(ii)*X+Kym(ii)*Y)); %magnon wavefunction
        for uu=1:3
            for vv=1:3
                if uu~=vv
                     qmEiEj{mm}{uu,vv}(ii)=sum(EjEi{uu,vv}.*WaveFactor,[1,2]); % Fourier transformation of non-equidistant system
                end
            end
        end
    end
    % multiWaitbar( 'CloseAll' );
    % toc(Clock1);

    for ii=1:Nfm
        tmp=0;
        for uu=1:3
            for vv=1:3
                tmp=tmp+qmEiEj{mm}{uu,vv}.*chiNiFe{uu,vv}(:,:,ii); %calculation of the BLS spectrum
            end
        end
        sigmaSW(ii,mm)=sum(abs(tmp).^2,[1,2]); % summing the intensity on the detector
    end
    % sigmaSW(:,mm)=sigmaSW(:,mm)/max(sigmaSW(:,mm)); % normalization of the signal, for comparison between this and experimental data just multiply with the highest experimental BLS signal value
end

%% Plot of qmEiEj matrix, last step before multiplication with chi (before taking magnetization into account)
% great way to visualize wavenumbers accessible by incident light (does not take magnetization into account -> no way of knowing if there are magnons there to be accessed)

Fig=figure(2);
    Fig.OuterPosition=[0 0 1550 1550];

    Label={{'$G_{xx} (k_{x},k_{y})$'},{'$G_{xy} (k_{x},k_{y})$'},{'$G_{xz} (k_{x},k_{y})$'};{'$G_{yx} (k_{x},k_{y})$'},{'$G_{yy} (k_{x},k_{y})$'},...
        {'$G_{yz} (k_{x},k_{y})$'};{'$G_{zx} (k_{x},k_{y})$'},{'$G_{zy} (k_{x},k_{y})$'},{'$G_{zz} (k_{x},k_{y})$'}};

    Ordering=[1,4,7,2,5,8,3,6,9];

    for ii=1:9
        subplot(3,3,ii);

        tmp=abs(qmEiEj{1}{Ordering(ii)});
        imagesc(kx*1e-6,kx*1e-6,tmp.');

        set(gca,'FontSize',14,'FontName','Times New Roman');
        set(gca,'YDir','normal');
        set(gca,'XGrid','on');
        set(gca,'YGrid','on');
        set(gca,'XMinorTick','on');
        set(gca,'YMinorTick','on');
        set(gca,'XLim',[min(kx) max(kx)]*1e-6);
        set(gca,'YLim',[min(kx) max(kx)]*1e-6);
        set(gca,'TickDir','both');
        xlabel('$k_{x}\,\mathrm{(1/\mu{m})}$','Interpreter','latex','FontSize',14,'FontName','Times New Roman');
        ylabel('$k_{y}\,\mathrm{(1/\mu{m})}$','Interpreter','latex','FontSize',14,'FontName','Times New Roman','VerticalAlignment','bottom');
        colormap(inferno(2048));
        axis equal

        hh=colorbar;
        set(hh,'Box','off','FontSize',14,'FontName','Times New Roman','TickDirection','out');
        hh.Label.String = Label{Ordering(ii)};
        hh.Label.Interpreter='latex';
        hh.Label.VerticalAlignment='top';
    end
%% ChiNiFe plot, overlap of this matrix and qmEiEj makes the final BLS signal

Fig=figure(3); %plot of the maps that highlight the ability of the system to access magnons at different kms
    Fig.OuterPosition=[0 0 1550 1550];

    Label={{'$\chi_{xx} (k_{x},k_{y})$'},{'$\chi_{xy} (k_{x},k_{y})$'},{'$\chi_{xz} (k_{x},k_{y})$'};{'$\chi_{yx} (k_{x},k_{y})$'},{'$\chi_{yy} (k_{x},k_{y})$'},...
        {'$\chi_{yz} (k_{x},k_{y})$'};{'$\chi_{zx} (k_{x},k_{y})$'},{'$\chi_{zy} (k_{x},k_{y})$'},{'$\chi_{zz} (k_{x},k_{y})$'}};

    Ordering=[1,4,7,2,5,8,3,6,9];

    for ii=1:9
        subplot(3,3,ii);

        tmp=abs(chiNiFe{Ordering(ii)}(:,:,20));
        imagesc(kx*1e-6,kx*1e-6,tmp.');

        set(gca,'FontSize',14,'FontName','Times New Roman');
        set(gca,'YDir','normal');
        set(gca,'XGrid','on');
        set(gca,'YGrid','on');
        set(gca,'XMinorTick','on');
        set(gca,'YMinorTick','on');
        set(gca,'XLim',[min(kx) max(kx)]*1e-6);
        set(gca,'YLim',[min(kx) max(kx)]*1e-6);
        set(gca,'TickDir','both');
        xlabel('$k_{x}\,\mathrm{(1/\mu{m})}$','Interpreter','latex','FontSize',14,'FontName','Times New Roman');
        ylabel('$k_{y}\,\mathrm{(1/\mu{m})}$','Interpreter','latex','FontSize',14,'FontName','Times New Roman','VerticalAlignment','bottom');
        colormap(inferno(2048));
        axis equal

        hh=colorbar;
        set(hh,'Box','off','FontSize',14,'FontName','Times New Roman','TickDirection','out');
        hh.Label.String = Label{Ordering(ii)};
        hh.Label.Interpreter='latex';
        hh.Label.VerticalAlignment='top';
    end

%% Plot of BLS signal, comparison between structures and reference thin film

Limits = [0,50]; % what frequency you are interested in
Fig=figure(4);
Fig.OuterPosition=[600 100 500 425];

Legendary=cell(M,1);
for ii=1:M
    plot(ffE,sigmaSW(:,ii),'-','Color',Color(ii,:),'LineWidth',1.5);
    hold on;
    Legendary{ii}=sprintf('$%.0f$\\,nm',Data{ii}.a*1e9);
end
Legendary{end}='bare film';
hold off;

set(gca,'FontSize',14,'FontName','Times New Roman');
set(gca,'XGrid','on');
set(gca,'YGrid','on');
set(gca,'XMinorTick','on');
set(gca,'YMinorTick','on');
set(gca,'XLim',Limits);

set(gca,'TickDir','both');
xlabel('frequency (GHz)','Interpreter','latex','FontSize',14,'FontName','Times New Roman');
ylabel('BLS signal','Interpreter','latex','FontSize',14,'FontName','Times New Roman','VerticalAlignment','bottom');
hh=legend(Legendary);
set(hh,'Box','off','Interpreter','latex','FontSize',14,'FontName','Times New Roman','Location','northeastoutside');   

%%

% Define folder name based on the current date and time
DateTimeString = datestr(now, 'yyyy-mm-dd_HH-MM-SS');
FolderName = ['Results_', DateTimeString];

% Create the folder if it doesn't exist
if ~exist(FolderName, 'dir')
    mkdir(FolderName);
end

% Save all figures to the folder
FigHandles = findall(0, 'Type', 'figure');
for i = 1:length(FigHandles)
    figHandle = FigHandles(i);
    
    % Set figure to invisible if running in a non-graphical environment
    % set(figHandle, 'Visible', 'off');
    
    % Define the file name for saving
    FileName = sprintf('Figure_%d.png', i);
    FilePath = fullfile(FolderName, FileName);
    
    % Save the figure as a .png file
    print(figHandle, FilePath, '-dpng', '-r300'); % Save as PNG with 300 DPI resolution
end

% Save variables ffE, qmEiEj, sigmaSW
save(fullfile(FolderName, 'Results.mat'), 'ffE', 'qmEiEj', 'sigmaSW');

% Save a copy of the current script
CurrentScriptName = mfilename('fullpath'); % get the full path of the current script
ScriptCopyName = fullfile(FolderName, [mfilename, '.m']);
copyfile([CurrentScriptName, '.m'], ScriptCopyName); % copy the script

% Save a copy of SpinWaveGreen.m
SpinWaveGreenPath = which('SpinWaveGreen'); % find SpinWaveGreen.m path
if ~isempty(SpinWaveGreenPath)
    copyfile(SpinWaveGreenPath, fullfile(FolderName, 'SpinWaveGreen.m')); % copy SpinWaveGreen.m to the folder
else
    warning('SpinWaveGreen.m not found. Please make sure it is in the MATLAB path.');
end

disp(['All figures, variables, and code copies saved in folder: ', FolderName]);



