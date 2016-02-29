function plotBeamCool()
rehash path 
rehash toolbox
set(0, 'DefaultUICOntrolFontSize', 16)                        % Let the font in listdlg bigger;
disp('Thank you for using frankliuao''s programs!')
nowPwd = pwd;                                                 % Get the file path where the m file is called;
[beamFileName, beamFilePath] = uigetfile({'*.beam','G4BL ASCII beam files (*.beam)';'*.*',...
                                                                     'All files (*.*)'},'Choose the beam file;', ...
                                                                     nowPwd);
% If the user pressed "Cancel":
if ~ischar(beamFileName) || ~ischar(beamFilePath)
    disp('Cancelling the program...')
    return;
end
%                                                               
beamFile = strcat(beamFilePath,beamFileName);
try
    disp(['Reading beam file',beamFile,', please wait......'])
    beamFile = dlmread(beamFile,' ',3,0);               % Read the beamFile to a matrix;
catch ME
    disp(['The file',beamFile,'can not be opened, please check again.'])
    return;
end
disp('Reading complete!')
beamFile = sortrows(beamFile,[8,9]);                                       % Sort the beam first by its particle name and then eventID;

plotList = {'X (mm)','Y (mm)', 'X'' (rad)', 'Y'' (rad)','Momentum (MeV/c)','t (ns)',...
                'Number of particles'};                                    % A list of available plot
[temp_s, temp_v] = listdlg('PromptString','Select what to be plotted from the list:',...
                                            'SelectionMode','Multiple',...
                                            'ListString',plotList,...
                                            'ListSize',[300,450]);
                                                                           % Let user choose which plots they want to make;
flagList = zeros(1,7);
if isempty(temp_s)
    disp({'Warning: No plot is chosen!!';'Exit.'})
    return
end
flagList(temp_s) = 1;                                         % Flag the indexes of plots to be plotted;
xFlag = flagList(1);
yFlag = flagList(2);
xPhaseFlag = flagList(3);
yPhaseFlag = flagList(4);
momFlag = flagList(5);
tFlag = flagList(6);
numberFlag = flagList(7);
%%
%
total2bplot = 0;
if xFlag == 1
    xFlagList = {'Y-Y'' phase space','t-dp/p longitudinal phase space','Just the histogram'};
    while (1)
        [xFlag_temp_s, xFlag_temp_v] = listdlg('PromptString','Plot X on which base?',...
                                                                      'SelectionMode','Multiple',...
                                                                      'ListString',xFlagList,...
                                                                      'ListSize',[300,450]);
        if ~isempty(xFlag_temp_s)
            break
        end
    end
    xFlagList = zeros(1,3);
    xFlagList(xFlag_temp_s) = 1;
    total2bplot = total2bplot + length(xFlag_temp_s);
end
if yFlag == 1
    yFlagList = {'X-X'' phase space','t-dp/p longitudinal phase space','Just the histogram'};
    while (1)
        [yFlag_temp_s, yFlag_temp_v] = listdlg('PromptString','Plot Y on which base?',...
                                                                      'SelectionMode','Multiple',...
                                                                      'ListString',yFlagList,...
                                                                      'ListSize',[300,450]);
        if ~isempty(yFlag_temp_s)
            break
        end
    end
    yFlagList = zeros(1,3);
    yFlagList(yFlag_temp_s) = 1;
    total2bplot = total2bplot + length(yFlag_temp_s);
end
if xPhaseFlag == 1
    xPhaseFlagList = {'Y-Y'' phase space','t-dp/p longitudinal phase space','Just the histogram'};
    while (1)
        [xPhaseFlag_temp_s, xPhaseFlag_temp_v] = listdlg('PromptString','Plot x'' on which base?',...
                                                                      'SelectionMode','Multiple',...
                                                                      'ListString',xPhaseFlagList,...
                                                                      'ListSize',[300,450]);
        if ~isempty(xPhaseFlag_temp_s)
            break
        end
    end
    xPhaseFlagList = zeros(1,3);
    xPhaseFlagList(xPhaseFlag_temp_s) = 1;
    total2bplot = total2bplot + length(xPhaseFlag_temp_s);
end
if yPhaseFlag == 1
    yPhaseFlagList = {'X-X'' phase space','t-dp/p longitudinal phase space','Just the histogram'};
    while (1)
        [yPhaseFlag_temp_s, yPhaseFlag_temp_v] = listdlg('PromptString','Plot y'' on which base?',...
                                                                      'SelectionMode','Multiple',...
                                                                      'ListString',yPhaseFlagList,...
                                                                      'ListSize',[300,450]);
        if ~isempty(yPhaseFlag_temp_s)
            break
        end
    end
    yPhaseFlagList = zeros(1,3);
    yPhaseFlagList(yPhaseFlag_temp_s) = 1;
    total2bplot = total2bplot + length(yPhaseFlag_temp_s);
end
if momFlag == 1
    momFlagList = {'X-X'' phase space','Y-Y'' phase space','X-Y real space','Just the histogram'};
    while (1)
        [momFlag_temp_s, momFlag_temp_v] = listdlg('PromptString','Plot momentum on which base?',...
                                                                      'SelectionMode','Multiple',...
                                                                      'ListString',momFlagList,...
                                                                      'ListSize',[300,450]);
        if ~isempty(momFlag_temp_s)
            break
        end
    end
    momFlagList = zeros(1,4);
    momFlagList(momFlag_temp_s) = 1;
    total2bplot = total2bplot + length(momFlag_temp_s);
end
if tFlag == 1
    tFlagList = {'X-X'' phase space','Y-Y'' phase space','X-Y real space','Just the histogram'};
    while (1)
        [tFlag_temp_s, tFlag_temp_v] = listdlg('PromptString','Plot time of arrival on which base?',...
                                                                      'SelectionMode','Multiple',...
                                                                      'ListString',tFlagList,...
                                                                      'ListSize',[300,450]);
        if ~isempty(tFlag_temp_s)
            break
        end
    end
    tFlagList = zeros(1,4);
    tFlagList(tFlag_temp_s) = 1;
    total2bplot = total2bplot + length(tFlag_temp_s);
end
if numberFlag == 1
    numberFlagList = {'X-X'' phase space','Y-Y'' phase space','X-Y real space','t-dp/p longitudinal phase space'};
    while (1)
        [numberFlag_temp_s,numberFlag_temp_v] = listdlg('PromptString','Plot number of particles on which base?',...
                                                                      'SelectionMode','Multiple',...
                                                                      'ListString',numberFlagList,...
                                                                      'ListSize',[300,450]);
        if ~isempty(numberFlag_temp_s)
            break
        end
    end
    numberFlagList = zeros(1,4);
    numberFlagList(numberFlag_temp_s) = 1;
    total2bplot = total2bplot + length(numberFlag_temp_s);
end

%%
% number of bins:
numBin = inputdlg({'The 2D histogram has n^2 bins.'},'Enter the number of bins of the 1D or 2D histogram',1,{'100'},'on');
numBin = str2double(numBin);
while (1)
    pCut_Bin = inputdlg({'In the form of [p_min, p_max]   Default [0,0] means no cut'},...
        'Enter the momentum histogram upper and lower cuts',1,{'[0,0]'},'on');
    % If histogram of momentum needed,
    % Ask user if a momentum cut is also needed;
    % Default: no cut.
    if ~isempty(pCut_Bin)
        break
    end
end
pCut_Bin = str2num(pCut_Bin{1}); %#ok<*ST2NM>
%%
particleKind = unique(beamFile(:,8));               % See how many kinds of particles are in the file;
kindCount = length(particleKind);
kindList = cell(1,kindCount);
for ii = 1:kindCount
    kindList{ii} = num2str(particleKind(ii));       % Store the particleKind in a cell array, to 
                                                                            % be used in listdlg below;
end
[temp_s, temp_v] = listdlg('PromptString','Select the particle ID you want to plot:',...
    'SelectionMode','Multiple',...
    'ListString',kindList,...
    'ListSize',[300,450]);
% Let user choose which particles to plot;
if ~temp_v
    return
end

beamGroup = cell(1,length(temp_s));
for ii = 1:length(temp_s)
    beamGroup{ii} = beamFile(beamFile(:,8)==particleKind(temp_s(ii)),:);
                                                                            % Divide the beam into groups;
                                                                            % Each group contains one kind of particles;
end

global figureHandle axesHandle
figureHandle = zeros(length(temp_s),total2bplot);          % figure handle recorder;
axesHandle = zeros(length(temp_s),total2bplot);            % axes handle recorder;

screenSize = get(0,'ScreenSize');                       % Get the screen size;
global figureLeft figureBottom screenWidth screenHeight
screenWidth = screenSize(1,3);                         % Width and height;
screenHeight = screenSize(1,4);
figureLeft = screenWidth/2-0.8*screenHeight/2;  % Set up the left edge position for figures;
figureBottom = screenHeight/2-0.8*screenHeight/2; % Bottom edge

for ii = 1:length(temp_s)
    % PDGID of the current kind of particles to be plotted:
    pdgID = kindList{temp_s(ii)};
    tot_mom = sqrt(beamGroup{ii}(:,4).^2+beamGroup{ii}(:,5).^2+...
                       beamGroup{ii}(:,6).^2);            % calculate the list of total momentum for
                                                                           % the current particle;
    if pCut_Bin(1)~=0 && pCut_Bin(2)~=0
        tot_mom = tot_mom(tot_mom<=pCut_Bin(2)&tot_mom>=pCut_Bin(1),:);
        beamGroup{ii} = beamGroup{ii}(tot_mom<=pCut_Bin(2)&tot_mom>=pCut_Bin(1),:);
    end
    xAngle = beamGroup{ii}(:,4)./beamGroup{ii}(:,6);    % x', angle, in rad;
    yAngle = beamGroup{ii}(:,5)./beamGroup{ii}(:,6);    % y', angle, in rad;
    beamX = beamGroup{ii}(:,1);
    beamY = beamGroup{ii}(:,2);
    beamTime = beamGroup{ii}(:,7);
    beamWeight = beamGroup{ii}(:,12);
    [xEdges,xCenter] = fr_findedge(beamGroup{ii}(:,1),numBin);
    [yEdges,yCenter] = fr_findedge(beamGroup{ii}(:,2),numBin);
    [xAngleEdges,xAngleCenter] = fr_findedge(xAngle(:,1),numBin);
    [yAngleEdges,yAngleCenter] = fr_findedge(yAngle(:,1),numBin);
    [momEdges,momCenter] = fr_findedge(tot_mom(:,1),numBin);
    [tEdges,tCenter] = fr_findedge(beamTime(:,1),numBin);
    toBplot = 1;
    if exist('xFlagList','var') 
        if length(unique(beamX))>1
            if xFlagList(1) ~= 0
                if length(unique(beamY))==1 || length(unique(yAngle))==1
                    display('Either the y or the y'' of the beam is the same for all the particles, y-y'' can not be plotted')
                else
                    x_yyHist = fr_histogram({yEdges,yAngleEdges,[beamY,yAngle]},beamX.*beamWeight,beamWeight);
                    figureHandle(ii,toBplot) = figure;
                    axesHandle(ii,toBplot) = axes('Parent',figureHandle(ii,toBplot),'FontSize',40);
                    plot1 = fr_plotScatter(x_yyHist,yCenter,yAngleCenter,figureHandle(ii,toBplot),axesHandle(ii,toBplot),pdgID,'x (mm)','yy');
                    toBplot = toBplot + 1;
                end
            end
            if xFlagList(2) ~= 0
                if length(unique(beamTime))==1 || length(unique(tot_mom))==1
                    display('Either the time or the momentum of the beam is the same for all the particles, t-p can not be plotted')
                else
                    x_tpHist = fr_histogram({tEdges,momEdges,[beamTime,tot_mom]},beamX.*beamWeight, beamWeight);
                    figureHandle(ii,toBplot) = figure;
                    axesHandle(ii,toBplot) = axes('Parent',figureHandle(ii,toBplot),'FontSize',40);
                    plot1 = fr_plotScatter(x_tpHist,tCenter,momCenter,figureHandle(ii,toBplot),axesHandle(ii,toBplot),pdgID,'x (mm)','tp');
                    toBplot = toBplot + 1;
                end
            end
            if xFlagList(3) ~= 0
                x_numberHist = fr_histogram2({xEdges,beamX},beamWeight);
                figureHandle(ii,toBplot) = figure;
                axesHandle(ii,toBplot) = axes('Parent',figureHandle(ii,toBplot),'FontSize',40);
                plot1 = fr_plotScatter2(x_numberHist,xCenter,figureHandle(ii,toBplot),axesHandle(ii,toBplot),pdgID,'x (mm)');
                toBplot = toBplot + 1;
            end
        end
    end
    if exist('xPhaseFlagList','var')
        if length(unique(xAngle))>1
            if xPhaseFlagList(1) ~= 0
                if length(unique(beamX))==1 || length(unique(xAngle))==1
                    display('Either the x or the x'' of the beam is the same for all the particles, x-x'' can not be plotted')
                else
                    xPhase_yyHist = fr_histogram({yEdges,yAngleEdges,[beamY,yAngle]},xAngle.*beamWeight, beamWeight);
                    figureHandle(ii,toBplot) = figure;
                    axesHandle(ii,toBplot) = axes('Parent',figureHandle(ii,toBplot),'FontSize',40);
                    plot1 = fr_plotScatter(xPhase_yyHist,yCenter,yAngleCenter,figureHandle(ii,toBplot),axesHandle(ii,toBplot),pdgID,'x'' (rad)','yy');
                    toBplot = toBplot + 1;
                end
            end
            if xPhaseFlagList(2) ~= 0
                if length(unique(beamTime))==1 || length(unique(tot_mom))==1
                    display('Either the time or the momentum of the beam is the same for all the particles, t-p can not be plotted')
                else
                    xPhase_tpHist = fr_histogram({tEdges,momEdges,[beamTime,tot_mom]},xAngle.*beamWeight, beamWeight);
                    figureHandle(ii,toBplot) = figure;
                    axesHandle(ii,toBplot) = axes('Parent',figureHandle(ii,toBplot),'FontSize',40);
                    plot1 = fr_plotScatter(xPhase_tpHist,tCenter,momCenter,figureHandle(ii,toBplot),axesHandle(ii,toBplot),pdgID,'x'' (rad)','tp');
                    toBplot = toBplot + 1;
                end
            end
            if xPhaseFlagList(3) ~= 0
                xPhase_numberHist = fr_histogram2({xAngleEdges,xAngle},beamWeight);
                figureHandle(ii,toBplot) = figure;
                axesHandle(ii,toBplot) = axes('Parent',figureHandle(ii,toBplot),'FontSize',40);
                plot1 = fr_plotScatter2(xPhase_numberHist,xAngleCenter,figureHandle(ii,toBplot),axesHandle(ii,toBplot),pdgID,'x'' (rad)');
                toBplot = toBplot + 1;
            end
        end
    end
    if exist('yFlagList','var')
        if length(unique(beamY))>1
            if yFlagList(1) ~= 0
                if length(unique(beamX))==1 || length(unique(xAngle))==1
                    display('Either the x or the x'' of the beam is the same for all the particles, x-x'' can not be plotted')
                else
                    y_xxHist = fr_histogram({xEdges,xAngleEdges,[beamX,xAngle]},beamY.*beamWeight, beamWeight);
                    figureHandle(ii,toBplot) = figure;
                    axesHandle(ii,toBplot) = axes('Parent',figureHandle(ii,toBplot),'FontSize',40);
                    plot1 = fr_plotScatter(y_xxHist,xCenter,xAngleCenter,figureHandle(ii,toBplot),axesHandle(ii,toBplot),pdgID,'y (mm)','xx');
                    toBplot = toBplot + 1;
                end
            end
            if yFlagList(2) ~= 0
                if length(unique(beamTime))==1 || length(unique(tot_mom))==1
                    display('Either the time or the momentum of the beam is the same for all the particles, t-p can not be plotted')
                else
                    y_tpHist = fr_histogram({tEdges,momEdges,[beamTime,tot_mom]},beamY.*beamWeight, beamWeight);
                    figureHandle(ii,toBplot) = figure;
                    axesHandle(ii,toBplot) = axes('Parent',figureHandle(ii,toBplot),'FontSize',40);
                    plot1 = fr_plotScatter(y_tpHist,tCenter,momCenter,figureHandle(ii,toBplot),axesHandle(ii,toBplot),pdgID,'y (mm)','tp');
                    toBplot = toBplot + 1;
                end
            end
            if yFlagList(3) ~= 0
                y_numberHist = fr_histogram2({yEdges,beamY},beamWeight);
                figureHandle(ii,toBplot) = figure;
                axesHandle(ii,toBplot) = axes('Parent',figureHandle(ii,toBplot),'FontSize',40);
                plot1 = fr_plotScatter2(y_numberHist,yCenter,figureHandle(ii,toBplot),axesHandle(ii,toBplot),pdgID,'y (mm)');
                toBplot = toBplot + 1;
            end
        end
    end
    if exist('yPhaseFlagList','var')
        if length(unique(yAngle))>1
            if yPhaseFlagList(1) ~= 0
                if length(unique(beamX))==1 || length(unique(xAngle))==1
                    display('Either the x or the x'' of the beam is the same for all the particles, x-x'' can not be plotted')
                else
                    yPhase_xxHist = fr_histogram({xEdges,xAngleEdges,[beamX,xAngle]},yAngle.*beamWeight, beamWeight);
                    figureHandle(ii,toBplot) = figure;
                    axesHandle(ii,toBplot) = axes('Parent',figureHandle(ii,toBplot),'FontSize',40); %#ok<*LAXES>
                    plot1 = fr_plotScatter(yPhase_xxHist,xCenter,xAngleCenter,figureHandle(ii,toBplot),axesHandle(ii,toBplot),pdgID,'y'' (rad)','xx'); %#ok<*NASGU>
                    toBplot = toBplot + 1;
                end
            end
            if yPhaseFlagList(2) ~= 0
                if length(unique(beamTime))==1 || length(unique(tot_mom))==1
                    display('Either the time or the momentum of the beam is the same for all the particles, t-p can not be plotted')
                else
                    yPhase_tpHist = fr_histogram({tEdges,momEdges,[beamTime,tot_mom]},yAngle.*beamWeight, beamWeight);
                    figureHandle(ii,toBplot) = figure;
                    axesHandle(ii,toBplot) = axes('Parent',figureHandle(ii,toBplot),'FontSize',40);
                    plot1 = fr_plotScatter(yPhase_tpHist,tCenter,momCenter,figureHandle(ii,toBplot),axesHandle(ii,toBplot),pdgID,'y'' (rad)','tp');
                    toBplot = toBplot + 1;
                end
            end
            if yPhaseFlagList(3) ~= 0
                yPhase_numberHist = fr_histogram2({yAngleEdges,yAngle},beamWeight);
                figureHandle(ii,toBplot) = figure;
                axesHandle(ii,toBplot) = axes('Parent',figureHandle(ii,toBplot),'FontSize',40);
                plot1 = fr_plotScatter2(yPhase_numberHist,yAngleCenter,figureHandle(ii,toBplot),axesHandle(ii,toBplot),pdgID,'y'' (rad)');
                toBplot = toBplot + 1;
            end
        end
    end
    if exist('momFlagList','var')
        if length(unique(tot_mom))>1
            if momFlagList(1) ~= 0
                if length(unique(beamX))==1 || length(unique(xAngle))==1
                    display('Either the x or the x'' of the beam is the same for all the particles, x-x'' can not be plotted')
                else
                    mom_xxHist = fr_histogram({xEdges,xAngleEdges,[beamX,xAngle]},tot_mom.*beamWeight, beamWeight);
                    figureHandle(ii,toBplot) = figure;
                    axesHandle(ii,toBplot) = axes('Parent',figureHandle(ii,toBplot),'FontSize',40);
                    plot1 = fr_plotScatter(mom_xxHist,xCenter,xAngleCenter,figureHandle(ii,toBplot),axesHandle(ii,toBplot),pdgID,'Momentum P (MeV/c)','xx');
                    toBplot = toBplot + 1;
                end
            end
            if momFlagList(2) ~= 0
                if length(unique(beamY))==1 || length(unique(yAngle))==1
                    display('Either the y or the y'' of the beam is the same for all the particles, y-y'' can not be plotted')
                else
                    mom_yyHist = fr_histogram({yEdges,yAngleEdges,[beamY,yAngle]},tot_mom.*beamWeight, beamWeight);
                    figureHandle(ii,toBplot) = figure;
                    axesHandle(ii,toBplot) = axes('Parent',figureHandle(ii,toBplot),'FontSize',40);
                    plot1 = fr_plotScatter(mom_yyHist,yCenter,yAngleCenter,figureHandle(ii,toBplot),axesHandle(ii,toBplot),pdgID,'Momentum P (MeV/c)','yy');
                    toBplot = toBplot + 1;
                end
            end
            if momFlagList(3) ~= 0
                if length(unique(beamX))==1 || length(unique(beamY))==1
                    display('Either the x or the y of the beam is the same for all the particles, x-y can not be plotted')
                else
                    mom_xyHist = fr_histogram({xEdges,yEdges,[beamX,beamY]},tot_mom.*beamWeight, beamWeight);
                    figureHandle(ii,toBplot) = figure;
                    axesHandle(ii,toBplot) = axes('Parent',figureHandle(ii,toBplot),'FontSize',40);
                    plot1 = fr_plotScatter(mom_xyHist,xCenter,yCenter,figureHandle(ii,toBplot),axesHandle(ii,toBplot),pdgID,'Momentum P (MeV/c)','xy');
                    toBplot = toBplot + 1;
                end
            end
            if momFlagList(4) ~= 0
                mom_numberHist = fr_histogram2({momEdges,tot_mom},beamWeight);
                figureHandle(ii,toBplot) = figure;
                axesHandle(ii,toBplot) = axes('Parent',figureHandle(ii,toBplot),'FontSize',40);
                plot1 = fr_plotScatter2(mom_numberHist,momCenter,figureHandle(ii,toBplot),axesHandle(ii,toBplot),pdgID,'Momentum P (MeV/c)');
                toBplot = toBplot + 1;
            end
        end
    end
    %
    if exist('tFlagList','var')
        if length(unique(beamTime))>1
            if tFlagList(1) ~= 0
                if length(unique(beamX))==1 || length(unique(xAngle))==1
                    display('Either the x or the x'' of the beam is the same for all the particles, x-x'' can not be plotted')
                else
                    t_xxHist = fr_histogram({xEdges,xAngleEdges,[beamX,xAngle]},beamTime.*beamWeight, beamWeight);
                    figureHandle(ii,toBplot) = figure;
                    axesHandle(ii,toBplot) = axes('Parent',figureHandle(ii,toBplot),'FontSize',40);
                    plot1 = fr_plotScatter(t_xxHist,xCenter,xAngleCenter,figureHandle(ii,toBplot),axesHandle(ii,toBplot),pdgID,'Time of Arrival (ns)','xx');
                    toBplot = toBplot + 1;
                end
            end
            if tFlagList(2) ~= 0
                if length(unique(beamY))==1 || length(unique(yAngle))==1
                    display('Either the y or the y'' of the beam is the same for all the particles, y-y'' can not be plotted')
                else
                    t_yyHist = fr_histogram({yEdges,yAngleEdges,[beamY,yAngle]},beamTime.*beamWeight, beamWeight);
                    figureHandle(ii,toBplot) = figure;
                    axesHandle(ii,toBplot) = axes('Parent',figureHandle(ii,toBplot),'FontSize',40);
                    plot1 = fr_plotScatter(t_yyHist,yCenter,yAngleCenter,figureHandle(ii,toBplot),axesHandle(ii,toBplot),pdgID,'Time of Arrival (ns)','yy');
                    toBplot = toBplot + 1;
                end
            end
            if tFlagList(3) ~= 0
                if length(unique(beamX))==1 || length(unique(beamY))==1
                    display('Either the x or the y of the beam is the same for all the particles, x-y can not be plotted')
                else
                    t_xyHist = fr_histogram({xEdges,yEdges,[beamX,beamY]},beamTime.*beamWeight, beamWeight);
                    figureHandle(ii,toBplot) = figure;
                    axesHandle(ii,toBplot) = axes('Parent',figureHandle(ii,toBplot),'FontSize',40);
                    plot1 = fr_plotScatter(t_xyHist,xCenter,yCenter,figureHandle(ii,toBplot),axesHandle(ii,toBplot),pdgID,'Time of Arrival (ns)','xy');
                    toBplot = toBplot + 1;
                end
            end
            if tFlagList(4) ~= 0
                t_numberHist = fr_histogram2({tEdges,beamTime},beamWeight);
                figureHandle(ii,toBplot) = figure;
                axesHandle(ii,toBplot) = axes('Parent',figureHandle(ii,toBplot),'FontSize',40);
                plot1 = fr_plotScatter2(t_numberHist,tCenter,figureHandle(ii,toBplot),axesHandle(ii,toBplot),pdgID,'Time of Arrival (ns)');
                toBplot = toBplot + 1;
            end
        end
    end
    %
    if exist('numberFlagList','var')
        if numberFlagList(1) ~= 0
            if length(unique(beamX))==1 || length(unique(xAngle))==1
                display('Either the x or the x'' of the beam is the same for all the particles, x-x'' can not be plotted')
            else
                number_xxHist = fr_histogram({xEdges,xAngleEdges,[beamX,xAngle]},beamWeight,zeros(length(beamWeight),1));
                figureHandle(ii,toBplot) = figure;
                axesHandle(ii,toBplot) = axes('Parent',figureHandle(ii,toBplot),'FontSize',40);
                plot1 = fr_plotScatter(number_xxHist,xCenter,xAngleCenter,figureHandle(ii,toBplot),axesHandle(ii,toBplot),pdgID,'Number of particles','xx');
                toBplot = toBplot + 1;
            end
        end
        if numberFlagList(2) ~= 0
            if length(unique(beamY))==1 || length(unique(yAngle))==1
                display('Either the y or the y'' of the beam is the same for all the particles, y-y'' can not be plotted')
            else
                number_yyHist = fr_histogram({yEdges,yAngleEdges,[beamY,yAngle]},beamWeight,zeros(length(beamWeight),1));
                figureHandle(ii,toBplot) = figure;
                axesHandle(ii,toBplot) = axes('Parent',figureHandle(ii,toBplot),'FontSize',40);
                plot1 = fr_plotScatter(number_yyHist,yCenter,yAngleCenter,figureHandle(ii,toBplot),axesHandle(ii,toBplot),pdgID,'Number of particles','yy');
                toBplot = toBplot + 1;
            end
        end
        if numberFlagList(3) ~= 0
            if length(unique(beamX))==1 || length(unique(beamY))==1
                display('Either the x or the y of the beam is the same for all the particles, x-y can not be plotted')
            else
                number_xyHist = fr_histogram({xEdges,yEdges,[beamX,beamY]},beamWeight,zeros(length(beamWeight),1));
                figureHandle(ii,toBplot) = figure;
                axesHandle(ii,toBplot) = axes('Parent',figureHandle(ii,toBplot),'FontSize',40);
                plot1 = fr_plotScatter(number_xyHist,xCenter,yCenter,figureHandle(ii,toBplot),axesHandle(ii,toBplot),pdgID,'Number of particles','xy');
                toBplot = toBplot + 1;
            end
        end
        if numberFlagList(4) ~= 0
            if length(unique(beamTime))==1 || length(unique(tot_mom))==1
                display('Either the time or the momentum of the beam is the same for all the particles, t-p can not be plotted')
            else
                number_tpHist = fr_histogram({tEdges,momEdges,[beamTime,tot_mom]},beamWeight,zeros(length(beamWeight),1));
                figureHandle(ii,toBplot) = figure;
                axesHandle(ii,toBplot) = axes('Parent',figureHandle(ii,toBplot),'FontSize',40);
                plot1 = fr_plotScatter(number_tpHist,tCenter,momCenter,figureHandle(ii,toBplot),axesHandle(ii,toBplot),pdgID,'Number of particles','tp');
                toBplot = toBplot + 1;
            end
        end
    end
end

end

%%
function [fr_Edge,fr_Center] = fr_findedge(x,fr_Bins)
% First create fr_Bins bins in [min(x),max(x)], e.g. [0,1,2,3] where fr_Bins=1
% Then shift the array by binsize/2, e.g. [-0.5,0.5,1.5,2.5];
% add the last number in the array. e.g. [-0.5,0.5,1.5,2.5,3.5]
fr_Center = linspace(min(x),max(x),fr_Bins);
fr_binSize = fr_Center(2)-fr_Center(1);
fr_Edge = [fr_Center-(fr_binSize)/2,fr_Center(end)+(fr_binSize)/2];
end

function [fr_hist] = fr_histogram(x,z_weight,weight)
% x is a cell, which contains 3 data sets, the first two are the bin edges of x and y, the 3rd one is the [x,y] data (??-by-2)
% z_weight is the value of z times the weights of the particles, weight_x is the weights of the particles. both (??-by-1)
% e.g. If there're two particles, we want to get the histogram of their momentum. 
%      z_weight = [P_1*W_1;P_2*W_2]
%      weight = [W_1;W_2]
    fr_xEdge = x{1};
    fr_yEdge = x{2};
    fr_data = x{3};
    fr_xEdgeMin = min(fr_xEdge);
    fr_yEdgeMin = min(fr_yEdge);
    fr_xBinWidth = fr_xEdge(2)-fr_xEdge(1);
    fr_yBinWidth = fr_yEdge(2)-fr_yEdge(1);
    fr_hist = zeros(length(fr_xEdge)-1,length(fr_yEdge)-1);
    fr_weight = zeros(length(fr_xEdge)-1,length(fr_yEdge)-1);
    for fr_ii = 1:1:size(fr_data,1)
        % Go through the particles and put them into the bins.
        temp1 = floor((fr_data(fr_ii,1)-fr_xEdgeMin)/fr_xBinWidth)+1;
        temp2 = floor((fr_data(fr_ii,2)-fr_yEdgeMin)/fr_yBinWidth)+1;
        fr_hist(temp1,temp2) = fr_hist(temp1,temp2) + z_weight(fr_ii);
        fr_weight(temp1,temp2) = fr_weight(temp1,temp2) + weight(fr_ii);
    end
    if all(~weight)==1
        fr_hist(fr_hist==0) = NaN;
    else
        fr_hist = fr_hist./fr_weight;
    end
end

function [fr_hist] = fr_histogram2(x,weight_x)
    fr_xEdge = x{1};
    fr_data = x{2};
    fr_xEdgeMin = min(fr_xEdge);
    fr_xBinWidth = fr_xEdge(2)-fr_xEdge(1);
    fr_hist = zeros(1,length(fr_xEdge)-1);
    for fr_ii = 1:length(fr_data)
       temp = floor((fr_data(fr_ii)-fr_xEdgeMin)/fr_xBinWidth)+1;
       fr_hist(temp) = fr_hist(temp) + weight_x(fr_ii);
    end
end

function [plot1] = fr_plotScatter(fr_histogram_data,fr_xBin,fr_yBin,fr_figureHandle,fr_axesHandle,fr_kind,fr_what,fr_base)
    global figureLeft figureBottom screenHeight
    fr_xBins = size(fr_histogram_data,1);
    fr_yBins = size(fr_histogram_data,2);
    scatterZ = zeros(1,fr_xBins*fr_yBins);
    scatterX = zeros(1,fr_xBins*fr_yBins);
    scatterY = zeros(1,fr_xBins*fr_yBins);
    vv = 1;
    for tt = 1:fr_xBins
        for uu = 1:fr_yBins
            scatterX(vv) = fr_xBin(tt);
            scatterY(vv) = fr_yBin(uu);
            scatterZ(vv) = fr_histogram_data(tt,uu);
            vv = vv+1;
        end
    end
    set(fr_figureHandle,'Color','w',...
              'Units','pixel',...
              'OuterPosition',[figureLeft,figureBottom,0.8*screenHeight,0.8*screenHeight])
   
    axes(fr_axesHandle)
    validPoints = ~isnan(scatterZ);
    scatterX = scatterX(validPoints);
    scatterY = scatterY(validPoints);
    scatterZ = scatterZ(validPoints);
    plot1 = scatter(scatterX,scatterY,90,scatterZ,'fill');
    switch fr_base
        case 'xx'
            xlabel('X (mm)','FontSize',40);
            ylabel('X''(rad)','FontSize',40);
            title({['The Horizontal Phase Space Distribution Plot for PDGid: ', fr_kind, '  '],['of ',fr_what,' on x-x''']},'FontSize',35,'FontName','Times');
            
        case 'yy'
            xlabel('Y (mm)','FontSize',40);
            ylabel('Y''(rad)','FontSize',40);
            title({['The Vertical Phase Space Distribution Plot for PDGid: ', fr_kind,'  '],['of ',fr_what,' on y-y''']},'FontSize',35,'FontName','Times');
        
        case 'tp'
            xlabel('t (ns)','FontSize',40);
            ylabel('p (MeV/c)','FontSize',40);
            title({['The Longitudinal Space Distribution Plot for PDGid: ', fr_kind,'  '],['of ',fr_what,' on t-p']},'FontSize',35,'FontName','Times');
            
        case 'xy'
            xlabel('X (mm)','FontSize',40);
            ylabel('Y (mm)','FontSize',40);
            title({['The Real Space Distribution Plot for PDGid: ', fr_kind,'  '],['of ',fr_what, 'on x-y']},'FontSize',35,'FontName','Times');
        
        otherwise
            error('Error Happened')
    end
    set(fr_axesHandle,'XGrid','off');
    set(fr_axesHandle,'YGrid','off');
    set(fr_axesHandle,'XMinorGrid','off');
    set(fr_axesHandle,'YMinorGrid','off');
    set(fr_axesHandle,'XMinorTick','on');
    set(fr_axesHandle,'YMinorTick','on');
    set(fr_axesHandle,'FontSize',40);
    colormap jet
    colorbar('FontSize',40)
    axis auto
    box(fr_axesHandle,'on');
end

function [plot1] = fr_plotScatter2(fr_histogram_data,fr_xBin,fr_figureHandle,fr_axesHandle,fr_kind,fr_what)
    global figureLeft figureBottom screenHeight
    set(fr_figureHandle,'Color','w',...
              'Units','pixel',...
              'OuterPosition',[figureLeft,figureBottom,0.8*screenHeight,0.8*screenHeight])
    axes(fr_axesHandle)
    plot1 = bar(fr_xBin,fr_histogram_data);
    histHandle = findobj(fr_axesHandle,'Type','patch');
    set(histHandle,'FaceColor',[0.4,0.4,0.4],'EdgeColor','w')
    set(plot1,'BarWidth',1)
    xlabel(fr_what,'FontSize',40);
    ylabel('Number of Particles','FontSize',40);
    title({['The Histogram of ',fr_what], [' for PDGid: ', fr_kind,'  ']},'FontSize',30,'FontName','Times');
    set(fr_axesHandle,'XGrid','off');
    set(fr_axesHandle,'YGrid','off');
    set(fr_axesHandle,'XMinorGrid','off');
    set(fr_axesHandle,'YMinorGrid','off');
    set(fr_axesHandle,'XMinorTick','on');
    set(fr_axesHandle,'YMinorTick','on');
    set(fr_axesHandle, 'FontSize', 40);
    box(fr_axesHandle,'on');
    
end
