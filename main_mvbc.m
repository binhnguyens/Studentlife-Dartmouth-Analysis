% Evaluate on simulation data
% Modified by Binh Nguyen 10/1/2021

clear;
clc;

%% Data

addpath ('/Users/binhnguyen/Documents/MATLAB/Data Preparation');
main_DP;

% First cluster
sz1 = 9;
sv1 = [5; 5; 5];
ini_v1 = 1;

% Second cluster
sz2 = 7;
sv2 = [5;5;5];
ini_v2 = 1;

%% Load data

avg_data = avg_view;
loc_data =  loc_view;
trend_data = trend_view;

M_avg_raw = avg_data;
M_loc_raw = loc_data;
M_trend_raw = trend_data;

n = size(avg_data, 1);
ad = size(avg_data, 2);
ld = size(loc_data, 2);
td = size(trend_data,2);

%% Label data
idx=find(ismember(phq_values{1}',phq_values{2}','rows'));
post = phq_values{4};
pre = phq_values{3}(idx); 
% Have to take idx because pre has more samples than post

% lbl = pre;
% lbl = post;
lbl = (pre+post)./2;

%% Normalize
% The function normc retains the relative magnitudes of the data in 
% each column with respect to the length of the vector they describe. 
% Normalising between min and max destroys those relationships.

M_avg_norm = normc(M_avg_raw);
M_trend_norm = normc(M_trend_raw);
M_loc_norm = normc(M_loc_raw);


M_avg = M_avg_norm;
M_trend = M_trend_norm;
M_loc = M_loc_norm;

M = cell(1, 1);
M{1} = M_avg;
M{2} = M_trend;
M{3} = M_loc;

%% Proposed Method
% Inputs:
%   M - a cell array of data matrixes from multiple views. Rows in each
%       matrix represent samples, columns represent features. It is assumed
%       samples are aligned well among all the matrix. In other words, rows with
%       same row index represent exactly the same sample but charaterized
%       from different views. All columns in matrixes are considered as
%       features in the analysis, so do not include sample id in the
%       matrixes.
%   sz - a scalar, hyperparameter, it controls the maximum non-zeros in vector z 
%   sv - a vector, hyperprameter, it constrols the maximum non-zeros in
%       each vector v_i, so it should have a length exact same as input M 
%   iSeedV1 - the index of the seed feature in view 1. This feature will be 
%       the one most likely remained active in view 1 (i.e., with a non-zero
%       component in v_1). This parameter helps the initialization and kick
%       off the iterative searching. Basically v1 is initialized with all zeros 
%       except this feature. 

% Outputs:
%   z - a vector, exactly the z in the formulation
%   U - a matrix, the i-th column is u_i in the formulation
%   V - a cell array, the i-th cell is v_i in the formulation
%   obj - the value of the objective function calculated with returned z,
%       U and V (singular value)

%% Cluster 1
[z1, U1, V1, obj1] = mvlrrl0(M, sz1, sv1, ini_v1);
fprintf('\nThe distribution of the true labels of the first identified cluster.\n');
holder = (lbl(z1 ~= 0));
tabulate (holder)

%%%% Info
% variable(z~=0) gets the indices of where the points are in this cluster
% avg_data(z1~=0,:) is the avg view in cluster 1
% z~=0 receives the indices of that cluster
% z==0 recieves the indices not in the cluster

%%%%% START TESTING
% Data storage for cluster 1 
label_c1 = lbl(z1 ~= 0);
data_c1 = cell (1,1);
for i =1:3
    data_c1{i} = M{i}(z1~=0,:);
end

% Visualization of data
cluster = 1; % Doesn't change as it is used for title
choice = 3; % Determines what view you have for all three clusters
rmv(M{choice},M{choice},lbl,lbl,z1,cluster,choice);
%%%%% END TESTING

M2_dv = cell(1, 1);
M2_dv{1} = M{1}(z1 == 0, :);
M2_dv{2} = M{2}(z1 == 0, :);
M2_dv{3} = M{3}(z1==0, :);
lbl2 = lbl(z1 == 0);

%% Cluster 2
[z2, U2, V2, obj2] = mvlrrl0(M2_dv, sz2, sv2, ini_v2);
fprintf('\nThe distribution of the true labels of the second identified cluster.\n');
% I think it should be lbl2 because it's labels that don't include the
% samples from the first cluster
holder = (lbl2(z2 ~= 0)); 
tabulate (holder)


%%%%% START TESTING
% Data storage for cluster 2
label_c2 = lbl2(z2 ~= 0);
data_c2 = cell (1,1);
for i =1:3
    data_c2{i} = M2_dv{i}(z2~=0,:);
end

% Visualization 
cluster = 2; % Doesn't change
rmv(M2_dv{choice},M{choice},lbl2,lbl,z2,cluster,choice);

M3_dv = cell(1, 1);
M3_dv{1} = M2_dv{1}(z2 == 0, :);
M3_dv{2} = M2_dv{2}(z2 == 0, :);
M3_dv{3} = M2_dv{3}(z2==0, :);
lbl3 = lbl2(z2 == 0);
%%%%% END TESTING

%% Cluster 3
%%%%% START TESTING
cluster = 3; % Title purposes 

% Z3 is everything leftover in M
z3 = ones (length(lbl3),1); 
rmv(M3_dv{choice},M{choice},lbl3,lbl,z3,cluster,choice);

% Data storage for cluster 3
label_c3 = lbl3(z3 ~= 0);
data_c3 = cell (1,1);
for i =1:3
    data_c3{i} = M3_dv{i}(z3~=0,:);
end
%%%%% END TESTING


%% Label density curve
% lbl(z1~=0)
% lbl2(z2~=0)
% lbl3(z3~=0)

figure
hold on;
[f,xi] = ksdensity(lbl(z1~=0)); 
plot(xi,f);
[f,xi] = ksdensity(lbl2(z2~=0)); 
plot(xi,f);
[f,xi] = ksdensity(lbl3(z3~=0)); 
plot(xi,f);
title ('Density curve of PHQ-9 -Pre');
xlabel ('PHQ score');
ylabel ('Frequency');
legend ('C1','C2','C3');

figure;
subplot (3,1,1);
histfit(lbl(z1~=0));
title('C1');
subplot (3,1,2);
histfit(lbl2(z2~=0));
title('C2');
subplot (3,1,3);
histfit(lbl3(z3~=0));
title('C3');

%% Key features

% M is the matrix of all the views
% V is the feature matrix from the algorithm
% Must transpose V since it's 12x1... need to be 1x12

% Key features from lit:  Convd, Darkd, Darkc, Audiov, Audioq and PhoneLockd

key_feature_cluster1 = cell(1, 1);
names = {'Act_s','Act_w','Act_r','Convo_d','Convo_c','Dark_d',...
            'Dark_c','Aud_q','Aud_v','Aud_n','Lock_d','Lock_c'};
names = names (:,[1 4 6 7 5 2 3 8 10 9 12 11]);

for i =1:3
    col_kf = transpose(V1{i});
    key_feature_cluster1{i} = M{i}(:,col_kf~=0);
    disp (names(:,col_kf~=0));
end

%% Comparison of label derived and label given
lbl_dv = double(z1 ~= 0);
lbl_dv(z1 == 0) = double(z2 ~= 0);
lbl_dv(z1 == 0 & lbl_dv == 1) = 2;


%% Save the results
rs = struct();
rs.sz1 = sz1;
rs.sv1 = sv1;
rs.z1 = z1;
rs.U1 = U1;
rs.V1 = V1;
rs.sz2 = sz2;
rs.sv2 = sv2;
rs.z2 = z2;
rs.U2 = U2;
rs.V2 = V2;
rs.lbl_dv = lbl_dv;



