load('nn_results_100_100_100.mat')

% y = 10.98;
% fileID = fopen('exptable.txt','w');
% fprintf(fileID, 'Exponential Function\n\n');
% fprintf(fileID,'%f %f\n',y);
% fclose(fileID);

fileID = fopen('bigger_controller','w');
data = 4;
fprintf(fileID,'%d \n',data);
data = 1;
fprintf(fileID,'%d \n',data);
data = 3;
fprintf(fileID,'%d \n',data);
data = 100;
fprintf(fileID,'%d \n',data);
data = 100;
fprintf(fileID,'%d \n',data);
data = 100;
fprintf(fileID,'%d \n',data);

%  Layer 1 data : 
for i=1:100
    for j = 1:4
        data = W1(i,j);
        fprintf(fileID,'%f \n',data);
    end
    data = b1(i);
    fprintf(fileID,'%f \n',data);
end

%  Layer 2 data : 
for i=1:100
    for j = 1:100
        data = W2(i,j);
        fprintf(fileID,'%f \n',data);
    end
    data = b2(i);
    fprintf(fileID,'%f \n',data);
end

%  Layer 3 data : 
for i=1:100
    for j = 1:100
        data = W3(i,j);
        fprintf(fileID,'%f \n',data);
    end
    data = b3(i);
    fprintf(fileID,'%f \n',data);
end

%  Layer 3 data : 
for i=1:1
    for j = 1:100
        data = W4(i,j);
        fprintf(fileID,'%f \n',data);
    end
    data = b4(i);
    fprintf(fileID,'%f \n',data);
end

fclose(fileID);

