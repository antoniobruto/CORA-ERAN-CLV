%% Calling ERAN and CORA
clear all;
close all;
clc;
%!python3 . --netname ../nets/inv_pendulum.tf
%ERAN_PATH='/Users/kekatos/Files/Projects/Gitlab/Tools/eran/our_tf_verify/';
timer_main = tic;

ERAN_PATH='/home/antonio/Workspace/Work_Verimag/eran/tf_verify';

CORA_PATH='/home/antonio/Workspace/Work_Verimag/eran/CORA';
addpath(genpath(CORA_PATH));
current_PATH=pwd;
current_PATH='/home/antonio/Workspace/Work_Verimag/linking';

% Adding Python3 to the PATH
PATH=getenv('PATH');
setenv('PATH', [PATH ':/usr/bin']);
setenv('PATH', [PATH ':/usr/local/bin']);

% Model
addpath(genpath('CORA_models'));

%%% Only Part that needs to be changed (corresponds to the model)

model_file='example_nonlinear_reach_TORA.m';
model_short_name='TORA';
NN_file=' /home/antonio/Workspace/Work_Verimag/eran/nets/neural_network_tora.tf';

%%%
ERAN_input_lb='/home/antonio/Workspace/Work_Verimag/eran/nets/lower_bounds.txt';
ERAN_input_ub='/home/antonio/Workspace/Work_Verimag/eran/nets/upper_bounds.txt';
ERAN_input_PATH='/home/antonio/Workspace/Work_Verimag/eran/nets/';

format long
zonotope_output=0;
%global iterations model_shoDroprt_name
%% Algorithm

% Prepare for ERAN execution

%!python3 . --netname ../nets/inv_pendulum.tf6

% Give IC for states or read them from CORA file
logFP = fopen("timeLog.txt","w");

iterations=2;
switch model_short_name
    case 'TORA'
        initial_ranges=[0.6 0.61; -0.7 -0.69; -0.4 -0.39; 0.5 0.51];
end
ranges{1}=initial_ranges;
timerStart = tic;
global timerCORA CORA_time timerERAN ERAN_time total_time; 
for iteration=1:iterations
    

    disp('')
    disp('')
    disp(['iteration: ' num2str(iteration)])
    disp('')
    disp('')
    %First Iteration - No ERAN
    
        %Other Iterations - With ERAN
        %Execute ERAN with Zonotope/Interval inputs
        timerERAN=tic;
        if ~zonotope_output
            % Interval Inputs:
            % We need lower and upper bounds for each variable that is
            % computed by CORA and then used as an input to ERAN.

            lb=ranges{iteration}(:,1);
            ub=ranges{iteration}(:,2);
           
            % we create a txt file with the lower bounds and 
            % one with the upper bounds
            fid = fopen( 'lower_bounds.txt', 'wt' );
            for i=1:length(lb)
                fprintf( fid, '%.8f \n', lb(i));
            end
            fclose(fid);

            fid = fopen( 'upper_bounds.txt', 'wt' );
            for i=1:length(lb)
                fprintf( fid, '%.8f \n', ub(i));
            end
            fclose(fid);

            % Move files to ERAN folder
            copyfile lower_bounds.txt '/home/antonio/Workspace/Work_Verimag/eran/nets/'
            copyfile upper_bounds.txt '/home/antonio/Workspace/Work_Verimag/eran/nets/'
            
            %Execute ERAN
            cd(ERAN_PATH)
            
            switch model_short_name
                case 'TORA'
                    !python3 . --netname ../nets/new_nets/bigger_controller.tf --domain deepzono --input text
            end
            copyfile 'NN_output_lower_bounds.txt' '/home/antonio/Workspace/Work_Verimag/eran/linking'
            copyfile 'NN_output_upper_bounds.txt' '/home/antonio/Workspace/Work_Verimag/eran/linking'

            cd(current_PATH)

            fileID = fopen('NN_output_lower_bounds.txt','r');
            formatSpec = '%f';
            NN_lb{iteration} = fscanf(fileID,formatSpec);
            fclose(fileID);

            fileID = fopen('NN_output_upper_bounds.txt','r');
            formatSpec = '%f';
            NN_ub{iteration} = fscanf(fileID,formatSpec);
            fclose(fileID);
            
        else
            %Zonotope Inputs
            
            %Create Zonotope for ERAN
            R_zono=zonotope(R_final_zonotope);
            R_values=get(R_zono,'Z');
            disp('The final zonotope is described below.')
            disp('')
            R_zono
            % need to write to a function to transform zonotopes from CORA to
            % ERAN format. CORA stores the centers as the 1st column (1 column *
            % number of dimensions). Next d columns describe generators (d *
            % number of dimensions).
            % ERAN is a text file described by one column (separated by a
            % delimiter). First row is number of inputs, second is the number of
            % generators d+1 (corresponding to the center). Then it goes dimension
            % by dimension. For dimension i, first is the center, then the first
            % generator, then the second, the third,..., the dth generator.
            
            cd(ERAN_input_PATH)
            
            % fix the numbering- TO DO
            
            zono_output_file=strcat('ERAN_inputs_',model_short_name,'_',num2str(iterations),'_iter','.txt');
            zono_output_file=strcat('ERAN_inputs_',model_short_name,'.txt');
            
            %                !python3 . --netname ../nets/neural_network_tora.tf --domain deepzono --input text
            !python3 . --netname ../nets/new_nets/bigger_controller.tf --domain deepzono --input text 
            
            fid = fopen(zono_output_file, 'wt' );
            no_outputs=length(ub);
            fprintf(fid,'%i\n',no_outputs);
            no_generators=size(R_values,2)-1;
            fprintf(fid,'%i\n',no_generators);
            matrix_zono=R_values';
            vector_zono=matrix_zono(:);
            for ii=1:length(vector_zono)
                fprintf(fid,'%i\n',vector_zono(ii));
            end
            cd(current_PATH)
            
            
            cd(ERAN_PATH)
            switch model_short_name
                case 'TORA'
                    !python3 . --netname ../nets/new_nets/bigger_controller.tf --domain deepzono --zonotope ../nets/ERAN_inputs_TORA.txt --complete 0 --input text
                    %!python3 . --netname ../nets/eran_wt/neural_network_tora.tf --domain deepzono --zonotope ../nets/ERAN_inputs_TORA.txt --complete 0 --input text
            end
            copyfile 'NN_output_lower_bounds.txt' '/home/antonio/Workspace/Work_Verimag/eran/linking'
            copyfile 'NN_output_upper_bounds.txt' '/home/antonio/Workspace/Work_Verimag/eran/linking'

            cd(current_PATH)


            fileID = fopen('NN_output_lower_bounds.txt','r');
            formatSpec = '%f';
            NN_lb{iteration} = fscanf(fileID,formatSpec);
            fclose(fileID);

            fileID = fopen('NN_output_upper_bounds.txt','r');
            formatSpec = '%f';
            NN_ub{iteration} = fscanf(fileID,formatSpec);
            fclose(fileID);
            
        
        end
        ERAN_time=toc(timerERAN);
        
        %Execute Plant
        %Take ERAN outputs and feed to CORA
        timerCORA=tic;
        if ~zonotope_output
            %Run CORA
            cd(current_PATH)
            cd CORA_models
            switch model_short_name
                case 'TORA'
                    %if it==1
                    %     [R,opt,~]=example_nonlinear_reach_TORA();
                    % else
                    [R,opt,~]=example_nonlinear_reach_TORA(NN_lb{iteration},NN_ub{iteration},lb,ub);
                    % end
            end
            % keep reachable set at the last time instance
            R_final_zonotope=R{end}{1};
            R_final_interval=interval(R_final_zonotope);

            % get(interval) does not exist
            % temporary solution - trick
            R_temp_zono=zonotope(R_final_interval);
            R_temp_values=get(R_temp_zono,'Z');
            j=2;
            for i=1:size(R_temp_values,1)
                final_intervals=[R_temp_values(i,1) - R_temp_values(i,j), R_temp_values(i,1) + R_temp_values(i,j)];
                j=j+1;
                intervals(i,:)=final_intervals;
            end
            ranges{iteration+1}=intervals;
            cd(current_PATH)
        else
            %Run CORA
            
            cd(current_PATH)
            cd CORA_models
            switch model_short_name
                case 'TORA'
                    %if it==1
                    %     [R,opt,~]=example_nonlinear_reach_TORA();
                    % else
                    [R,opt,~]=example_nonlinear_reach_TORA(NN_lb{iteration},NN_ub{iteration},[],[],R_zono);
                    % end
            end
            % keep reachable set at the last time instance
            R_final_zonotope=R{end}{1};
            R_final_interval=interval(R_final_zonotope);

            cd(current_PATH)
        end
        CORA_time = toc(timerCORA);
    
    fprintf(logFP,'Iteration %d \n',iteration);
    if exist('ERAN_time','var')
        fprintf(logFP,'\t ERAN: %f\n',ERAN_time);
    end
    if exist('CORA_time','var')
        fprintf(logFP,'\t CORA: %f\n',CORA_time);
    end
    fprintf(logFP,'\n');
    
end
fprintf(logFP,'\nTotal Time = %f\n',toc(timer_main));
total_time = toc(timer_main);

%% PRINT outputs to file

cd outputs
if ~zonotope_output
    output_file=strcat('output_',model_short_name,'_',num2str(iterations),'_iter_',num2str(opt.tFinal),'_w_dt_',num2str(opt.timeStep),'_sec.txt');
    output_figure=strcat('output_',model_short_name,'_',num2str(iterations),'_iter_',num2str(opt.tFinal),'_w_dt_',num2str(opt.timeStep),'_sec');
    
    fid = fopen(output_file, 'wt' );
    fprintf(fid,'\n');
    fprintf(fid,'The model under examination is %s.\n',model_short_name);
    fprintf(fid,'\n');
    
    fprintf(fid,'X_0: %s\n',mat2str(ranges{1}));
    fprintf(fid,'\n');
    
    fprintf(fid,'The time step is: %s\n',mat2str(opt.tFinal));
    fprintf(fid,'\n');
    for iteration=1:iterations
        fprintf( fid, 'Iteration: %d \n', iteration);
        fprintf(fid,'\n');
        NN_combined_ranges=[cell2mat(NN_lb);cell2mat(NN_ub)]';
        NN_combined_ranges_string=mat2str(NN_combined_ranges(iteration,:));
        fprintf(fid, 'ERAN outputs the controller ranges: %s\n',NN_combined_ranges_string);
        fprintf(fid,'\n');
        fprintf(fid, 'CORA outputs the state ranges: %s\n',mat2str(ranges{iteration+1}));
        fprintf(fid,'\n');
        
    end
    
    fprintf(fid,'The total computation time is %s.\n',total_time);
    fprintf(fid,'\n');
    fclose(fid);
    
    savefig(gcf,strcat(output_figure,'.fig'))
    saveas(gcf,strcat(output_figure,'.png'))
else
    disp("Need to fix print options for zonotopes")
    
    output_file=strcat('output_',model_short_name,'_',num2str(iterations),'_iter_',num2str(opt.tFinal),'_w_dt_',num2str(opt.timeStep),'_sec_zono.txt');
    output_figure=strcat('output_',model_short_name,'_',num2str(iterations),'_iter_',num2str(opt.tFinal),'_w_dt_',num2str(opt.timeStep),'_sec_zono');
    
    savefig(gcf,strcat(output_figure,'.fig'))
    saveas(gcf,strcat(output_figure,'.png'))
end
fclose(logFP);
