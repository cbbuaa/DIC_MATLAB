% function []=DIC_main()
% x should along vertical, and y along horizontal
clear; clc; close all;

%% --------------input parameters-------------------%
% set parameters, including subset size, iteration algorithm, strain
% calculation window size, initializing method, etc.
Params               = paramset();
subplot(1,2,1);
ImRef                = imageRead('ref.bmp',Params,1);
imshow(uint8(ImRef))

Params.sizeX         = size(ImRef,1);
Params.sizeY         = size(ImRef,2);


% calculation points on the reference image
[calPtX,calPtY,Params] = calcuPt(Params);
comptPoints          = Params.comptPoints;

Params.NumCalPt      = Params.Lx * Params.Ly;
pCoord               = Params.InitP;
% Jacobian
Params.gridW2P       = [Params.localSub(1,:)',Params.localSub(2,:)'];
% Gradient of reference image, needs noly be calculated once=
[gradxImR, gradyImR] = gradImg(ImRef);


%% Here is a loop. Reading deformed image and registrating
% Implement of DIC-2D
subplot(1,2,2);
ImDef                = imageRead('def_3.bmp',Params,2);
imshow(uint8(ImDef)),drawnow;
% Initialization
if Params.fixed_seedPts
    InitMatP         = Params.InitP(1:2,1);
else
    InitMatP         = findInitP(ImRef,ImDef,Params,Params.InitP(:,1),Params.InitP(:,1));
end
Params.InitDispP     = InitMatP-Params.InitP(1:2,1);

% Registration
switch Params.IterMethod
    case 'IC-GN'
        outP                 = zeros(Params.NumCalPt,6);
        
        p                    = [Params.InitDispP(1),0,0,Params.InitDispP(2),0,0]';
        [Disp,strain,ZNCC,outIter] = DIC2D(ImRef,ImDef,pCoord,p,Params,...
                                   gradxImR,gradyImR,comptPoints,outP);
    case 'ILS'
        [gradxImD, gradyImD] = gradImg(ImDef);
        % InitP  
        outP                 = zeros(Params.NumCalPt,8);
        p                    = [Params.InitDispP(1),0,0,Params.InitDispP(2),0,0,1,0]';
        [Disp,strain,ZNCC,outIter] = DIC2D(ImRef,ImDef,pCoord,p,Params,...
                                gradxImD,gradyImD,comptPoints,outP);
    case 'N-R'
        [gradxImD, gradyImD] = gradImg(ImDef);
        % InitP  
        outP                 = zeros(Params.NumCalPt,6);
        p                    = [Params.InitDispP(1),0,0,Params.InitDispP(2),0,0]';
        [Disp,strain,ZNCC,outIter] = DIC2D(ImRef,ImDef,pCoord,p,Params,...
                                gradxImD,gradyImD,comptPoints,outP);
    case 'Wznssd NR'
        [gradxImD, gradyImD] = gradImg(ImDef);
        % InitP  
        outP                 = zeros(Params.NumCalPt,8);
        p                    = [Params.InitDispP(1),0,0,Params.InitDispP(2),0,0,floor(Params.subset/2)*sqrt(2)]';
        [Disp,strain,ZNCC,outIter] = DIC2D(ImRef,ImDef,pCoord,p,Params,...
                                gradxImD,gradyImD,comptPoints,outP);
    case 'IC-GN2'
        % InitP  
        outP                 = zeros(Params.NumCalPt,12);
        p                    = [Params.InitDispP(1),0,0,0,0,0,Params.InitDispP(2),0,0,0,0,0]';
        [Disp,strain,ZNCC,outIter] = DIC2D(ImRef,ImDef,pCoord,p,Params,...
                                gradxImR,gradyImR,comptPoints,outP);  
end

%% plot the results
plotResults(Params, ImRef, calPtX, calPtY, Disp, strain,'exx');
load('Param.mat');
strain_real = Param.strain;
ux_real     = reshape(strain_real(:,1),[Params.sizeX,Params.sizeY]);

strain_err  = ux_real(calPtX,calPtY)*1e6 - reshape(strain(:,1),Params.Lx, Params.Ly);
mean(strain_err(:))
std(strain_err(:))
max(strain_err(:))
min(strain_err(:))

disp_real  = Param.disp;
u_real     = reshape(disp_real(:,1),[Params.sizeX,Params.sizeY]);
u_real_mat = u_real(calPtX,calPtY);
u_meas_mat = reshape(Disp(:,1),Params.Lx, Params.Ly);
disp_err   = u_real_mat - u_meas_mat;

mean(disp_err(:))
std(disp_err(:))
max(disp_err(:))
min(disp_err(:))

numel(find(ZNCC < 0.9))
mean(outIter(find(ZNCC >= 0.9)))
mean(ZNCC(find(ZNCC >= 0.9)))
