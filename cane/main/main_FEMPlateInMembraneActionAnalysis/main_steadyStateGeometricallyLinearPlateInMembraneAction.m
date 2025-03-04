%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Script documentation
%
% Task : Plane stress analysis for a rectangular plate subject to uniform
%        pressure on its top edge
%
% Date : 19.02.2014
%
%% Preamble
clear;
clc;
close all;

%% Includes

% Add general math functions
addpath('../../generalMath/');

% Add all functions related to parsing
addpath('../../parsers/');

% Add all functions related to the low order basis functions
addpath('../../basisFunctions/');

% Add all equation system solvers
addpath('../../equationSystemSolvers/');

% Add all the efficient computation functions
addpath('../../efficientComputation/');

% Add all functions related to plate in membrane action analysis
addpath('../../FEMPlateInMembraneActionAnalysis/solvers/',...
        '../../FEMPlateInMembraneActionAnalysis/solutionMatricesAndVectors/',...
        '../../FEMPlateInMembraneActionAnalysis/loads/',...
        '../../FEMPlateInMembraneActionAnalysis/graphics/',...
        '../../FEMPlateInMembraneActionAnalysis/output/',...
        '../../FEMPlateInMembraneActionAnalysis/postprocessing/',...
        '../../FEMPlateInMembraneActionAnalysis/errorComputation/');
    
% Include performance optimzed functions
addpath('../../efficientComputation/');

%% Parse data from GiD input file

% Define the path to the case
pathToCase = '../../inputGiD/FEMPlateInMembraneActionAnalysis/';
% caseName = 'infinitePlateWithHoleQuadrilaterals';
% caseName = 'cantileverBeamPlaneStress';
% caseName = 'PlateWithAHolePlaneStress';
% caseName = 'PlateWithMultipleHolesPlaneStress';
% caseName = 'InfinitePlateWithAHolePlaneStress';
% caseName = 'unitTest_curvedPlateTipShearPlaneStress';
% caseName = 'gammaStructureMixedElementsPlaneStress';
% caseName = 'NACA2412_AoA5_CSD';
caseName = 'trial_slab_with_openings';

% Parse the data from the GiD input file
[strMsh, homDOFs, inhomDOFs, valuesInhomDOFs, propNBC, propAnalysis, ...
    parameters, propNLinearAnalysis, ~, propGaussInt] = ...
    parse_StructuralModelFromGid(pathToCase, caseName, 'outputEnabled');

%% GUI

% On the body forces
computeBodyForces = @computeConstantVerticalStructureBodyForceVct;

% Choose equation system solver
solve_LinearSystem = @solve_LinearSystemMatlabBackslashSolver;
% solve_LinearSystem = solve_LinearSystemGMResWithIncompleteLUPreconditioning;

% Output properties
propOutput.isOutput = true;
propOutput.writeOutputToFile = @writeOutputFEMPlateInMembraneActionToVTK;
propOutput.VTKResultFile = 'undefined';

% Choose computation of the stiffness matrix
computeStiffMtxLoadVct = @computeStiffMtxAndLoadVctFEMPlateInMembraneActionCST;
% computeStiffMtxLoadVct = @computeStiffMtxAndLoadVctFEMPlateInMembraneActionMixed;

% Quadrature for the stiffness matrix and the load vector of the problem
% 'default', 'user'
propIntDomain.type = 'default';
propIntDomain.noGP = 1;

% Quadrature for the L2-norm of the error
intError.type = 'user';
intError.noGP = 4;

% Linear analysis
propStrDynamics = 'undefined';

% Initialize graphics index
graph.index = 1;

%% Output data to a VTK format
pathToOutput = '../../outputVTK/FEMPlateInMembraneActionAnalysis/';

%% Compute the load vector
t = 0;
F = computeLoadVctFEMPlateInMembraneAction ...
    (strMsh, propAnalysis, propNBC, t, propGaussInt, 'outputEnabled');

%% Visualization of the configuration
graph.index = plot_referenceConfigurationFEMPlateInMembraneAction ...
    (strMsh, propAnalysis, F, homDOFs, [], graph, 'outputEnabled');

%% Initialize solution
numNodes = length(strMsh.nodes(:,1));
numDOFs = 2*numNodes;
dHat = zeros(numDOFs,1);

%% Solve the plate in membrane action problem
[dHat, FComplete, minElSize] = solve_FEMPlateInMembraneAction ...
    (propAnalysis, strMsh, dHat, homDOFs, inhomDOFs, valuesInhomDOFs, ...
    propNBC, computeBodyForces, parameters, computeStiffMtxLoadVct, ...
    solve_LinearSystem, propNLinearAnalysis, propIntDomain, propOutput, ...
    caseName, pathToOutput, 'outputEnabled');

% Stiffness matrix calculation 
loadFactor = 1;
numNodes = length(strMsh.nodes(:,1));
numDOFs = 2*numNodes;
u = zeros(numDOFs, 1);
DOFNumbering = 1:numDOFs;
precomResVec = 'undefined';
precompStiffMtx = 'undefined';
uMeshALE = 'undefined';
uDot = 'undefined';
uDotSaved = 'undefined';
uSaved = 'undefined';
[K, F, minElEdgeSize] = computeStiffMtxLoadVct(propAnalysis, u, uSaved, uDot, uDotSaved, uMeshALE, ...
    precompStiffMtx, precomResVec, DOFNumbering, strMsh, F, ...
    loadFactor, computeBodyForces, propStrDynamics, t, ...
    parameters, propGaussInt);

%% Postprocessing
graph.visualization.geometry = 'reference_and_current';
resultant = 'stress';
component = '1Principal';
nodeIDs_active = 'undefined';
contactSegments = 'undefined';
graph.index = plot_currentConfigurationAndResultants ...
    (propAnalysis, strMsh, homDOFs, dHat, nodeIDs_active, contactSegments, ...
    parameters, resultant, component, graph);

% Compute the error in the L2-norm for the case of the plane stress
% analysis over a quarter annulus plate subject to tip shear force
if strcmp(caseName,'unitTest_curvedPlateTipShearPlaneStress')
    nodeNeumann = strMsh.nodes(propNBC.nodes(1, 1), 2:end);
    funHandle = str2func(propNBC.fctHandle(1, :));
    forceAmplitude = norm(funHandle(nodeNeumann(1, 1),nodeNeumann(1, 2),nodeNeumann(1, 3), 0));
    internalRadius = 4;
    externalRadius = 5;
    propError.resultant = 'stress';
    propError.component = 'tensor';
    errorL2 = computeRelErrorL2CurvedBeamTipShearFEMPlateInMembraneAction ...
        (strMsh, dHat, parameters, internalRadius, externalRadius, ...
        forceAmplitude, propError, intError, 'outputEnabled');
end

%% Load Path Calculation
% Calculation of Strain Energy, U, of Total System 

    LoadNodes = propNBC.nodes;
    DOFsLoadNodes = zeros(2*length(LoadNodes),1);
    DOFsLoadNodes(2:2:2*length(LoadNodes)) = 2*LoadNodes;
    DOFsLoadNodes(1:2:end) = 2*LoadNodes-1;

    SupportNodes = zeros(length(homDOFs),1);
    for i = 1:length(homDOFs)
        if floor(homDOFs(i)/2) == homDOFs(i)/2
            SupportNodes(i) = homDOFs(i)/2;
        else 
            SupportNodes(i) = (homDOFs(i)+1)/2;
        end
    end
    SupportNodes = unique(SupportNodes);
    RightEdgeSupportNodes = SupportNodes(find((strMsh.nodes(SupportNodes,2) == max(strMsh.nodes(:,2)))));
    DOFsSupportNodes = homDOFs';
    DOFsRightEdgeSupportNodes = zeros(length(RightEdgeSupportNodes),1);
    DOFsRightEdgeSupportNodes(1:end) = 2*RightEdgeSupportNodes-1;

    % SupportAndLoadNodes_Ordered = sort([LoadNodes; SupportNodes]);
    SupportAndLoadNodes_Ordered = sort([LoadNodes; RightEdgeSupportNodes]);
    DOFsSupportAndLoadNodes_Ordered = sort([DOFsLoadNodes; DOFsRightEdgeSupportNodes]);
    FreeNodes = (1:numNodes)';
    FreeNodes(SupportAndLoadNodes_Ordered) = [];
    DOFsFreeNodes = 1:2*numNodes;
    DOFsFreeNodes(DOFsSupportAndLoadNodes_Ordered) = [];
    DOFsFreeNodes = DOFsFreeNodes';
    

% Calculating U for the total system simply by using the calculated displacements:
    U_index_total = (1/2)*transpose(K*dHat)*dHat;

%% Calculation of Strain Energy, U, For Each Node
    U_free_nodes = zeros(numNodes,1);
    U_index = zeros(numNodes,1);


% Using an updated displacement vector by solving for each system:
    for kk = 1:numNodes
        fixedNodeDOF1 = (kk)*2-1;
        
        if ~ismember(fixedNodeDOF1,DOFsSupportAndLoadNodes_Ordered) && ~ismember(fixedNodeDOF1+1,DOFsSupportAndLoadNodes_Ordered)
            DOFsKnown_Ordered = sort([DOFsSupportAndLoadNodes_Ordered; fixedNodeDOF1; fixedNodeDOF1+1]);
            dHatFixedNode = dHat;
            dHatFixedNode(fixedNodeDOF1:fixedNodeDOF1+1,1) = 0;
            dHatKnown = dHatFixedNode(DOFsKnown_Ordered);

            DOFsUnknown = (1:numDOFs)';
            DOFsUnknown(DOFsKnown_Ordered) = [];

            K_fixed_node = K;
            K_fixed_node(DOFsKnown_Ordered,:) = [];

            K_known_displ = K_fixed_node(:, DOFsKnown_Ordered);
            F_known_displ = -1*(K_known_displ * dHatKnown);


            K_fixed_node(:,DOFsKnown_Ordered) = [];


            dHatUnknown = K_fixed_node\F_known_displ;

            dHat_fixed_node = zeros(numDOFs,1);
            dHat_fixed_node(DOFsKnown_Ordered) = dHatKnown;
            dHat_fixed_node(DOFsUnknown) = dHatUnknown;

            U_free_nodes(kk,1) = (1/2)*transpose(K*dHat_fixed_node)*dHat_fixed_node;
            U_index(kk,1) = 1- U_index_total/U_free_nodes(kk);

        end

    end

U_index(LoadNodes) = 1;

%% Visualization of load path - updated

XNodesPlot = strMsh.nodes(:,2);
YNodesPlot = strMsh.nodes(:,3);

[X,Y]=meshgrid(linspace(min(XNodesPlot),max(XNodesPlot),50),linspace(min(YNodesPlot),max(YNodesPlot),50));
zq = griddata(XNodesPlot, YNodesPlot, U_index, X,Y,'cubic');

figure;
surf(X, Y, zq);
xlabel('X');
ylabel('Y');
zlabel('U* Index');
title('U* Index');
%% Streamlines - updated
UX_zq = gradient(reshape(zq,1,[]), reshape(X,1,[]));
UY_zq = gradient(reshape(zq,1,[]), reshape(Y,1,[]));

UX_zq = reshape(UX_zq,size(X));
UY_zq = reshape(UY_zq,size(X));
UX_zq(isinf(UX_zq)) = 1;
UX_zq(isnan(UX_zq)) = 0;


RightEdgeSupportNodes = SupportNodes(find((strMsh.nodes(SupportNodes,2) == max(strMsh.nodes(:,2)))));

% XSupportNodes = strMsh.nodes(SupportNodes,2);
% YSupportNodes = strMsh.nodes(SupportNodes,3);
XSupportNodes = strMsh.nodes(RightEdgeSupportNodes,2)-0.1;
YSupportNodes = strMsh.nodes(RightEdgeSupportNodes,3);
[startX,startY] = meshgrid(XSupportNodes,YSupportNodes);
figure;
lineobj = stream2(X, Y,UX_zq,UY_zq,startX,startY);
lineobj = streamline(lineobj);
%% Visualization of load path
XNodesPlot = strMsh.nodes(:,2);
YNodesPlot = strMsh.nodes(:,3);

[X,Y]=meshgrid(XNodesPlot,YNodesPlot);
Z = griddata(XNodesPlot, YNodesPlot, U_index, X, Y);
figure;
[C, hCont] = contourf(X,Y,Z,'linecolor','none','LineWidth',0.01,"ShowText","on","LabelSpacing",600);
daspect([2 2 2]);
[hLeg, levMat] = contourLegend(hCont);
hLeg.Location = 'southeastoutside'; 
hLeg.FontSize = 7;
hold on

% hold off
U_x_grad = gradient(U_index,XNodesPlot);
U_x_grad(isinf(U_x_grad)) = 1;
U_x_grad(isnan(U_x_grad)) = 0;
U_y_grad = gradient(U_index,YNodesPlot);



% Using streamlines
XSupportNodes = strMsh.nodes(SupportNodes,2);
YSupportNodes = strMsh.nodes(SupportNodes,3);
[startX,startY] = meshgrid(XSupportNodes,YSupportNodes);
X_unique = unique(strMsh.nodes(:,2));
Y_unique = unique(strMsh.nodes(:,3));

UX  = zeros(length(Y_unique),length(X_unique))+1e-16;
UY  = zeros(length(Y_unique),length(X_unique))+1e-16;



for kk = 1:numNodes
    XCoordinate = strMsh.nodes(kk,2);
    YCoordinate = strMsh.nodes(kk,3);
    
    posX_unique = find(X_unique == XCoordinate);
    posY_unique = find(Y_unique == YCoordinate);

    UX(posY_unique, posX_unique) = U_x_grad(kk);
    UY(posY_unique, posX_unique) = U_y_grad(kk);
end

[XGrid, YGrid] = meshgrid(X_unique,Y_unique);

% figure;
lineobj = stream2(XGrid, YGrid,UX,UY,startX,startY);

lineobj = streamline(lineobj);
% lineobj(1).Color = "k";
% lineobj(2).Color = "m";
for i = 1:length(lineobj)
    lineobj(i).Color = "k";
    lineobj(i).LineWidth = 3;
    % legend('lineobj(i)','')
end
% legend('lineobj','')
XLoadNodes = strMsh.nodes(LoadNodes,2);
YLoadNodes = strMsh.nodes(LoadNodes,3);
% [startX_inv,startY_inv] = meshgrid(propStr.XLx-0.01,strMsh.nodes(LoadNodes,2));
[startX_inv,startY_inv] = meshgrid(XLoadNodes,YLoadNodes);
verts_inv = stream2(XGrid, YGrid,-UX,-UY,startX_inv,startY_inv);
lineobj_inv = streamline(verts_inv);
% for i = 1:length(lineobj_inv)
%     lineobj_inv(i).Color = "#7E2F8E";
%     lineobj_inv(i).LineWidth = 2;
% end

% hold off


% % Streamlines properties:
% linobj1 = lineobj(1);
% linobj2 = lineobj(2);
% 
% linobj1X = linobj1.XData;
% linobj1Y = linobj1.YData;
% 
% linobj2X = linobj2.XData;
% linobj2Y = linobj2.YData;
% % Streamline lengths, found using Pythagoran theorem:
% n1 = numel(linobj1X);
% length1 = 0;
% for i = 1:n1-1
%     length1 = length1 + sqrt( (linobj1X(i+1)-linobj1X(i))^2 + (linobj1Y(i+1)-linobj1Y(i))^2 );
% end
% n2 = numel(linobj2X);
% length2 = 0;
% for i = 1:n2-1
%     length2 = length2 + sqrt( (linobj2X(i+1)-linobj2X(i))^2 + (linobj2Y(i+1)-linobj2Y(i))^2 );
% end

% Streamline curvatures:
% C1 = LineCurvature2D([transpose(linobj1X),transpose(linobj1Y)]);
% C2 = LineCurvature2D([transpose(linobj2X),transpose(linobj2Y)]);

% plot(linspace(0,1,length(C1)), C1)
% plot(linspace(0,1,length(C2)), C2)

%% Using streamslice
% Flow from supports to load nodes:
figure;
[sliceverts, sliceaverts] = streamslice(XGrid, YGrid,UX,UY);
slicelinobj = streamslice(XGrid, YGrid,UX,UY);
% xlim([propStr.X0,propStr.XLx])
% ylim([propStr.Y0,propStr.YLy])
hold off
%% Flow from load nodes to supports:
figure;
streamslice(XGrid, YGrid,-UX,-UY)
% xlim([propStr.X0,propStr.XLx])
% ylim([propStr.Y0,propStr.YLy])
hold off
%% END OF THE SCRIPT
