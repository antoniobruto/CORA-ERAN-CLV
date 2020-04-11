function [Rcont,options,completed] = example_nonlinear_reach_TORA(u_min,u_max,lb,ub,R_zono)
%  example_nonlinear_reach_01_TORA - example of nonlinear reachability
% analysis; this example is also a unit test function.
%
% This example can be found at http://www.cs.colorado.edu/~srirams/papers/neural-reachability-hscc2019.pdf of
%
%
% Syntax:  
%    example_nonlinear_reach_01_TORA()
%
% Inputs:
%    no
%
% Outputs:
%    res - boolean
%
% Example:9
%
%
% Author:       Nikolaos Kekatos
% Written:      3-Oct-2019
% Last update:  ---
% Last revision:---


%------------- BEGIN CODE --------------


%set options---------------------------------------------------------------
options.tStart=0; %start time
%options.tFinal=0.01; %final time
options.tFinal=1; %final time
mptopt('lpsolver', 'CDD', 'qpsolver', 'QUADPROG')
%options.linprog = optimset('linprog');
% cplex gets broken for t_final=1.5s

if nargin<1
    options.x0=[0.6; -0.7;-0.4;0.59]; %initial state for simulation
    Z0_center=[ 0.605;-0.695;-0.395;0.595];
    %Z0=zonotope([Z0_center, 0.005*diag([1;1;1;1])]);
    %Z0_center=[ -0.886730353745207;  -0.228282763115134; 0.262973979334848;  -0.249234458087221];  
   % Z0_center=[ -0.274137231365497;  -0.891059033891944;   0.141295604189229;   0.492591208378458];

    Z0=zonotope([Z0_center, 0.0001*diag([1;1;1;1])]);
elseif nargin==5
    Z0=R_zono;
elseif nargin==4
    Xt_min=[lb];
    Xt_max=[ub];
    Xt_average=(Xt_min+Xt_max)/2;
    try        
        Z0_center=Xt_average;
        Xt_deviation=Xt_max-Xt_average;
        Z0=zonotope([Z0_center,diag(Xt_deviation)]);
        
    catch
        warning('Check again for different dimensions');
    end
else
    error('incorrect inputs');
end
% verify ranges by interval(Z0);
options.R0=Z0;
%options.timeStep=0.005; %time step size for reachable set computation
options.timeStep=0.01; %time step size for reachable set computation
options.taylorTerms=20; %number of taylor terms for reachable sets
options.zonotopeOrder=30; %zonotope order
options.polytopeOrder=4; %polytope order

options.advancedLinErrorComp = 0;
options.tensorOrder=1;
%options.tensorOrder=3;
%options.errorOrder=3;
options.reductionTechnique='girard';
%options.filterLength = [10, 5];
options.maxError=0.01*[0.1; 1;0.1;1];
options.reductionInterval=2;
options.verbose = 1;

%uncertain inputs
if nargin<1
  %  options.uTrans = 9.89;%9.87643<=u<=9.91252
      options.uTrans=-0.097408791621541;

    options.U = zonotope([options.uTrans,0.02]); %input for reachability analysis
 %   options.uTrans=-0.741825666465679;
 %   options.U = zonotope([options.uTrans,eps]); 
else
    NN_min=[u_min];
    NN_max=[u_max];
    NN_average=(NN_min+NN_max)/2;
    try        
        options.uTrans=NN_average;
        NN_deviation=NN_max-NN_average;
        options.U=zonotope([options.uTrans,diag(NN_deviation)]);
        
    catch
        warning('Check again for different dimensions');
    end 
end
%uTrans=[1; 0; 0; 0.5; -0.5];
%options.uTrans=uTrans; %center of input set
%options.U=0.5*zonotope([zeros(5,1),diag([0.2, 0.5, 0.2, 0.5, 0.5])]); %input for reachability analysis
%---------------------------
%--------------------------------------------------------------------------

%specify continuous dynamics-----------------------------------------------
TORA=nonlinearSys(4,1,@TORA_ode,options); %initialize TORA model
%--------------------------------------------------------------------------
      
tic
%compute reachable set
Rcont = reach(TORA, options);
tComp = toc;
disp(['computation time of reachable set: ',num2str(tComp)]);

% 
% %create random simulations; RRTs would provide better results, but are
% %computationally more demanding
% runs = 60;
% fractionVertices = 0.5;
% fractionInputVertices = 0.5;
% inputChanges = 6;
% simRes = simulate_random(TORA, options, runs, fractionVertices, fractionInputVertices, inputChanges);

%plotting message
%disp('Start plotting; takes a while since plotting acceleration for zonotope bundles not yet implemented');

%plot results-------------------741825666465679-------------------------------------------
projectedDimensions=[1 3];
plotOrder = 20;

%figure
%figure;
hold on

%plot reachable sets 
for i=1:length(Rcont)
    for j=1:length(Rcont{i})
        Zproj = reduce(Rcont{i}{j},'girard',plotOrder);
        plotFilled(Zproj,projectedDimensions,[.3 .3 .3],'EdgeColor','none');drawnow;
        %j = j+2;
    end
end

%plot initial set
plotFilled(options.R0,projectedDimensions,'w','EdgeColor','k');
disp(options)

% %plot simulation results      
% for i=1:length(simRes.t)
%     plot(simRes.x{i}(:,projectedDimensions(1)),simRes.x{i}(:,projectedDimensions(2)),'Color',0*[1 1 1]);
% end

%öabel plot
xlabel(['x_{',num2str(projectedDimensions(1)),'}']);
ylabel(['x_{',num2str(projectedDimensions(2)),'}']);
%--------------------------------------------------------------------------


%example completed
completed = 1;

%------------- END OF CODE --------------
%{
figure
%plot results-------------------741825666465679-------------------------------------------
projectedDimensions=[2 4];
plotOrder = 20;
    
%figure;
hold on

%plot reachable sets 
for i=1:length(Rcont)
    for j=1:length(Rcont{i})
        Zproj = reduce(Rcont{i}{j},'girard',plotOrder);
        plotFilled(Zproj,projectedDimensions,[.2 .3 .3],'EdgeColor','none')%;drawnow;
        %j = j+2;
    end
end

%plot initial set
plotFilled(options.R0,projectedDimensions,'w','EdgeColor','k');
disp(options)

% %plot simulation results      
% for i=1:length(simRes.t)
%     plot(simRes.x{i}(:,projectedDimensions(1)),simRes.x{i}(:,projectedDimensions(2)),'Color',0*[1 1 1]);
% end

%öabel plot
xlabel(['x_{',num2str(projectedDimensions(1)),'}']);
ylabel(['x_{',num2str(projectedDimensions(2)),'}']);
%--------------------------------------------------------------------------
%}