%% Calling ERAN and CORA

clear all;
close all;
clc;
%!python3 . --netname ../nets/inv_pendulum.tf

nnv='eran';
% nnv='sherlock';
ERAN_PATH='/Users/kekatos/Files/Projects/Gitlab/Tools/eran/tf_verify_no_specs/';
ERAN_PATH='/Users/kekatos/Files/Projects/Gitlab/Tools/eran/our_tf_verify/';

SHERLOCK_PATH='/Users/kekatos/Files/Projects/Gitlab/Tools/sherlock';

FLOWSTAR_PATH='/Users/kekatos/Files/Projects/Gitlab/Tools/flowstar-2.1.0';

CORA_PATH='/Users/kekatos/Files/Projects/Gitlab/Tools/CORA_2018';
addpath(genpath(CORA_PATH));
current_PATH=pwd;
current_PATH='/Users/kekatos/Files/Projects/Gitlab/Matlab_Python_Interfacing';

% Adding Python3 to the PATH
PATH=getenv('PATH');
setenv('PATH', [PATH ':/usr/local/bin']);

% Model
addpath(genpath('CORA_models'));

%%% Only Part that needs to be changed (corresponds to the model)

model_file='example_nonlinear_reach_TORA.m';
model_short_name='TORA';
NN_file='/Users/kekatos/Files/Projects/Gitlab/Tools/eran/nets/eran_wt/neural_network_tora.tf';
NN_file='/Users/kekatos/Files/Projects/Gitlab/Tools/eran/nets/eran_wt/bigger_controller.tf';

%%%
ERAN_input_lb='/Users/kekatos/Files/Projects/Gitlab/Tools/eran/nets/lower_bounds.txt';
ERAN_input_ub='/Users/kekatos/Files/Projects/Gitlab/Tools/eran/nets/upper_bounds.txt';
ERAN_input_PATH='/Users/kekatos/Files/Projects/Gitlab/Tools/eran/nets/';

format long
zonotope_output=1;
%global iterations model_short_name

if zonotope_output==1 && strcmp(nnv,'sherlock')
    error('Sherlock only outputs intervals. It cannot use zonotopes.')
end
%
%% Algorithm

% Prepare for ERAN execution

%!python3 . --netname ../nets/inv_pendulum.tf

% Give IC for states or read them from CORA file
total_time=tic;

iterations=10;
switch model_short_name
    case 'TORA'
        initial_ranges=[0.6 0.61; -0.7 -0.69; -0.4 -0.39; 0.59 0.6];
%         initial_ranges=[0.6 0.7; -0.7 -0.6; -0.4 -0.3; 0.5 0.6];
        
end
ranges{1}=initial_ranges;
for it=1:iterations
    
    disp('')
    disp('')
    disp(['iteration: ' num2str(it)])
    disp('')
    disp('')
    if ~zonotope_output || it==1
        % in case of intervals we need lower and upper bounds for each
        % variable that is computed by CORA and then used as an input to
        % ERAN.
        
        lb=ranges{it}(:,1);
        ub=ranges{it}(:,2);
       
        % we create a txt file with the lower bounds and one with the upper
        % bounds
        % name_lb=strcat(model_short_name,'_lower_bounds.txt');
        % fid = fopen( name_lb, 'wt' );
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
        
        % move files to ERAN folder
        
        %movefile upper_bounds.txt ERAN_input_PATH
        %movefile lower_bounds.txt ERAN_input_PATH
        copyfile lower_bounds.txt '/Users/kekatos/Files/Projects/Gitlab/Tools/eran/nets/'
        copyfile upper_bounds.txt '/Users/kekatos/Files/Projects/Gitlab/Tools/eran/nets/'
        
        
        % call ERAN
        if strcmp(nnv,'eran')
            cd(ERAN_PATH)
            switch model_short_name
                case 'TORA'
                    !python3 . --netname ../nets/eran_wt/neural_network_tora.tf --domain deepzono --input text
            end
            copyfile 'NN_output_lower_bounds.txt' '/Users/kekatos/Files/Projects/Gitlab/Matlab_Python_Interfacing'
            copyfile 'NN_output_upper_bounds.txt' '/Users/kekatos/Files/Projects/Gitlab/Matlab_Python_Interfacing'
            
            cd(current_PATH)
        elseif strcmp(nnv,'sherlock')
            sherlock_benchmark=99;
            cd(SHERLOCK_PATH)
            switch model_short_name
                case 'TORA'
                    command_temp=strcat("./run_file ", num2str(sherlock_benchmark))
                    %lb -> lower bounds, ub -> upper bounds
                    inputs_sherlock=[lb,ub]';
                    inputs_sherlock_str=num2str(inputs_sherlock(:)');
                    inputs_sherlock_str_space=strcat({' '},inputs_sherlock_str);
                    %matlab ignores trailing white space but keeps leading
                    command_sherlock=strcat(command_temp,inputs_sherlock_str_space);
                    system(command_sherlock)
            end
            copyfile 'NN_output_lower_bounds.txt' '/Users/kekatos/Files/Projects/Gitlab/Matlab_Python_Interfacing'
            copyfile 'NN_output_upper_bounds.txt' '/Users/kekatos/Files/Projects/Gitlab/Matlab_Python_Interfacing'
            cd(current_PATH)
    
        end
        
        fileID = fopen('NN_output_lower_bounds.txt','r');
        formatSpec = '%f';
        NN_lb{it} = fscanf(fileID,formatSpec);
        fclose(fileID);
        
        fileID = fopen('NN_output_upper_bounds.txt','r');
        formatSpec = '%f';
        NN_ub{it} = fscanf(fileID,formatSpec);
        fclose(fileID);
        
        % saved NN_lb for minimum values of controller
        % saved NN_ub for maximum values
        
        cd CORA_models
        switch model_short_name
            case 'TORA'
                %if it==1
                %     [R,opt,~]=example_nonlinear_reach_TORA();
                % else
                [R,opt,~]=example_nonlinear_reach_TORA(NN_lb{it}-10,NN_ub{it}-10,lb,ub);
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
        ranges{it+1}=intervals;
        cd(current_PATH)
        %total_time=toc;
        
        % for the first iteration we need to have intervals and construct
        % the zonotope file in CORA and translate it to boxes.
        if it==1
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
            
            
            fid = fopen(zono_output_file, 'wt' );
            no_outputs=length(ub);
            fprintf(fid,'%i\n',no_outputs);
            no_generators=size(R_values,2)-1;
   % ERAN requires the number of generators +1 for the center
        fprintf(fid,'%i\n',no_generators+1);
            matrix_zono=R_values';
            vector_zono=matrix_zono(:);
            for ii=1:length(vector_zono)
                fprintf(fid,'%i\n',vector_zono(ii));
            end
            cd(current_PATH)
        end
    else
        %in case of zonotopes we can no longer use lower_bounds.txt and
        %upper_bounds.txt. We have to create a new txt file and call eran
        %with the flag --zonotope.
        
        % call ERAN
        %lb=ranges{it}(:,1);
        %ub=ranges{it}(:,2);
        cd(ERAN_PATH)
        switch model_short_name
            case 'TORA'
                !python3 . --netname ../nets/eran_wt/neural_network_tora.tf --domain deepzono --zonotope ../nets/ERAN_inputs_TORA.txt --complete 0 --input text
        end
        copyfile 'NN_output_lower_bounds.txt' '/Users/kekatos/Files/Projects/Gitlab/Matlab_Python_Interfacing'
        copyfile 'NN_output_upper_bounds.txt' '/Users/kekatos/Files/Projects/Gitlab/Matlab_Python_Interfacing'
        
        cd(current_PATH)
        
        
        fileID = fopen('NN_output_lower_bounds.txt','r');
        formatSpec = '%f';
        NN_lb{it} = fscanf(fileID,formatSpec);
        fclose(fileID);
        
        fileID = fopen('NN_output_upper_bounds.txt','r');
        formatSpec = '%f';
        NN_ub{it} = fscanf(fileID,formatSpec);
        fclose(fileID);
        
        cd CORA_models
        switch model_short_name
            case 'TORA'
                %if it==1
                %     [R,opt,~]=example_nonlinear_reach_TORA();
                % else
                [R,opt,~]=example_nonlinear_reach_TORA(NN_lb{it}-10,NN_ub{it}-10,[],[],R_zono);
                % end
        end
        % keep reachable set at the last time instance
        R_final_zonotope=R{end}{1};
        R_final_interval=interval(R_final_zonotope);
        
        warning("Zonotope outputs require modifications");
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
        zono_output_file=strcat('ERAN_inputs_',model_short_name,'_',num2str(iterations),'_iter','.txt');
        
        fid = fopen(zono_output_file, 'wt' );
        no_outputs=length(ub);
        fprintf(fid,'%i\n',no_outputs);
        no_generators=size(R_values,2)-1;
        % ERAN requires the number of generators +1 for the center
        fprintf(fid,'%i\n',no_generators+1);
        matrix_zono=R_values';
        vector_zono=matrix_zono(:);
        for ii=1:length(vector_zono)
            fprintf(fid,'%i\n',vector_zono(ii));
        end
        cd(current_PATH)
    end
    
    
end
toc(total_time);

fprintf(" The total elapsed time for %i iterations is %f seconds.\n",iterations,toc(total_time));

%% PRINT outputs to file

% stop
cd outputs
if ~zonotope_output
    output_file=strcat(nnv,'_output_',model_short_name,'_',num2str(iterations),'_iter_',num2str(opt.tFinal),'_w_dt_',num2str(opt.timeStep),'_sec.txt');
    output_figure=strcat(nnv,'_output_',model_short_name,'_',num2str(iterations),'_iter_',num2str(opt.tFinal),'_w_dt_',num2str(opt.timeStep),'_sec');
    
    fid = fopen(output_file, 'wt' );
    fprintf(fid,'\n');
    fprintf(fid,'The model under examination is %s.\n',model_short_name);
    fprintf(fid,'\n');
    
    fprintf(fid,'X_0: %s\n',mat2str(ranges{1}));
    fprintf(fid,'\n');
    
    fprintf(fid,'The time step is: %s\n',mat2str(opt.tFinal));
    fprintf(fid,'\n');
    for it=1:iterations
        fprintf( fid, 'Iteration: %d \n', it);
        fprintf(fid,'\n');
        NN_combined_ranges=[cell2mat(NN_lb);cell2mat(NN_ub)]';
        NN_combined_ranges_string=mat2str(NN_combined_ranges(it,:));
        fprintf(fid, '%s outputs the controller ranges: %s\n',nnv,NN_combined_ranges_string);
        fprintf(fid,'\n');
        fprintf(fid, 'CORA outputs the state ranges: %s\n',mat2str(ranges{it+1}));
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
