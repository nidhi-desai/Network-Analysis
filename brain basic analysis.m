%% Brain Functional Connectivity graph analysis
% by Nidhi Desai
%% Loading data
load("FC.mat");
FC_thres=FC;
FC_thres(FC_thres<=0.2)=0;
figure, image(FC-FC');colorbar; % The differences in values is in the range of 10^-16
                                %so making the matrix symmetric makes sense 
A = (FC_thres + FC_thres.')/2; % To make the graph undirected 
A(logical(eye(size(A)))) = 0; % Removing self loaping
imagesc(A);

%% General network properties
% Is this a connected graph?
SP = batch_paths(A);
nnz(isinf(SP.distance));% since greater than zero so not a connected graph
[comps,comp_sizes] = get_components(A); 
length(comp_sizes)% since this is not 1, the graph is not connected

% Finding the indices of nodes which are not connected to any other nodes
rows_has_all_zeros = ~any(A,2);
indices_zeros = find(rows_has_all_zeros);
% Removing nodes which are not connected to any other nodes
A(indices_zeros,:)=[];
A(:,indices_zeros)=[];

% Now check if the graph is connected
[comps,comp_sizes] = get_components(A); 
length(comp_sizes)

%% Question 2 - Network visualization
convert2cytoscape(A,'cyto_input.csv');%Converting Adjacency matrix to cytoscape input
plot(graph(A));% Looking at the graph in matlab

% Import region names
[num, region_names] = xlsread('Book1.xlsx');

% Effect of thresholding shown through adjacency matrix visualization
FC1=FC;
FC1(FC1<=0)=0;
FC1(logical(eye(size(FC1)))) = 0; 
FC1(indices_zeros,:)=[];
FC1(:,indices_zeros)=[];
figure;
imagesc(FC1);colorbar;
figure;
imagesc(A);colorbar;

%% Question 3 - Assesment of network integration
SP = batch_paths(A);

% Measures of network integration
[SP.distance,SP.length,SP.B] = get_shortest_path_lengths(1./A); %distance can be expressed as -log(A), as 1./A and in other ways...
mask_ut = triu(true(size(A)),1);
Erout_A = mean(1./SP.distance(mask_ut)); %efficency of routing (usually called global efficiency)
[SP] = batch_paths(A);
char_path_length_A = mean2(SP.distance);
P = f_markov_chain(A);
mfpt_A = mean2(f_mfpt(P));

%% Question 4 - Assesment of network segregation thorugh modularity and network organization into communities 
gamma = [0;0.25;0.5;0.75;1;2;3;4;5;6;7;8;9;10]% gamma has to be a column vectors
[Q,S] = run_louvain(A,gamma); % Q is Q-modularity score - higher the better
figure;
plot(gamma,Q); title('Modularity score vs. penalization factor'); xlabel('Penalization term (gamma)'); ylabel('Modularity score');
figure;
plot(gamma,S); title('Number of modules vs. penalization factor'); xlabel('Penalization term (gamma)'); ylabel('Number of modules'); % Number of modules
% Choose gamma and give reason (platau region)

% Measures of network segregation
Clus_Coef_A = mean(clustering_coef_wu(A));
Q_A = run_louvain(A,1); 
    
%% Question 5 - Randomization procedures with respect to your network 
[diss, Q_score, path_len] = track_with_swaps(A,14,1);
num_of_changes = 8000;
% Creating ensemble of randomized networks
rand_network = zeros(100,159,159);
for i =1:100
    [rug,e,d] = rand_xswap_wu(A,num_of_changes,num_of_changes*50,i);
    rand_network(i,:,:) = full(rug);
end
rand1 = rand_network(i,:,:)           
%Visualizing 3 random neetworks
rand_vec_3 = round(1+100*rand(1,3)); 
for i = 1:3
    figure;
    imagesc(squeeze(rand_network(rand_vec_3(i),:,:)));colorbar;
end

%% Question 6 - Is network different than the ensemble of random networks
% Measures of segregation
Clus_Coef = zeros(100,1); %Clustering coefficiency
Q_ensemble = zeros(100,1); % Q is Q-modularity score 

% Measures of integration
Erout = zeros(100,1); %efficiency of routing - global efficiency
char_path_length = zeros(100,1); % Characteristic path length
mfpt = zeros(100,1);

for i = 1:100
    B = squeeze(rand_network(i,:,:));
    [SP.distance,SP.length,SP.B] = get_shortest_path_lengths(1./B); %distance can be expressed as -log(A), as 1./A and in other ways...
    mask_ut = triu(true(size(B)),1);
    Erout(i) = mean(1./SP.distance(mask_ut)); %efficency of routing (usually called global efficiency)
    Clus_Coef(i) = mean(clustering_coef_wu(B));
    Q_ensemble(i) = run_louvain(B,5); 
    [SP] = batch_paths(B);
    char_path_length(i) = mean2(SP.distance);
    P = f_markov_chain(B);
    mfpt(i) = mean2(f_mfpt(P));
end

% Z-score
Clus_Coef_z = (0.225-0.0686)/0.0015;
Q_z = (0.0617-0.011)/0.0016;
Erout_z = (0.1729-0.2414)/(6.1887*10^-4);
char_path_length_z = (7.1456-4.526)/0.0243;
mfpt_z = (357.049-276.0624)/13.4117;

% Histogram distribution 
x_pos = mfpt_A; y_pos = 0;
histfit(mfpt); hold on; plot(x_pos,y_pos,'r*');
xlabel('mean first passage time'); ylabel('frequency'); 
title('histogram of mean first passage time of ensemble');

%% Question 7 - Power law fit of network degree distribution
nnodes = length(A);
degree_node = zeros(nnodes,1);
for n=1:nnodes
    degree_node(n) = nnz(A(n,:)); %Degree of each node
end
degree = unique(degree_node);
% How many nodes have that degree
freq = zeros(length(degree),1);
for n=1:length(degree)
    freq(n) = length(find(degree(n)==degree_node));
end
[slope, intercept,MSE, R2] = logfit(degree,freq,'loglog');
title('log-log plot of frequency of degree vs. degree');
xlabel('log(frequency of degree)');
ylabel('log(degree)');
y_log_estimate = (10^intercept)*(degree.^slope);
poly = polyfit(degree,freq,1);
y_line_estimate = poly(1)*degree + poly(2);
figure;
histogram(degree_node); hold on;
plot(degree, y_log_estimate); hold on;
plot(degree, y_line_estimate);
title('Power law fit for network degree-distribution');
xlabel('Degree of a node');
ylabel('Frequency of a degree');

%% Question 8 - Small World Property

Clus_Coef_A = mean(clustering_coef_wu(A));
[SP_A]=batch_paths(A);
char_path_length_A=mean2(SP_A.distance);

Clus_Coef_rand = zeros(100,1);
char_path_length_rand = zeros(100,1);
small_world = zeros(100,1);
for i=1:100
    Clus_Coef_rand(i) = mean(clustering_coef_wu(squeeze(rand_network(i,:,:))));
    [SP_rand] = batch_paths(squeeze(rand_network(i,:,:)));
    char_path_length_rand(i) = mean2(SP_rand.distance);
    small_world(i) = (Clus_Coef_A/Clus_Coef_rand(i))/(char_path_length_A/char_path_length_rand(i));
end

boxplot(small_world);
xlabel('Random Networks'); ylabel('Measure of small-worldness');
title('Box plot of measure of small-worldness for random networks');

%% Question 9 - Network measure not covered in class
% Leverage centrality 
lev_vec_A = leverage_centrality(A);
lev_vec_rand = zeros(100,159);
for i = 1:100
    lev_vec_rand(i,:) = leverage_centrality(squeeze(rand_network(i,:,:)));
end

% Export leverage centrality to cytospace
nodeNumbers = 1:159; nodeNumbers = nodeNumbers';
positive_lev = 10*abs(lev_vec_A);
csvwrite('leverage.txt',[nodeNumbers,positive_lev]);

lev_sign = sign(lev_vec_A);
csvwrite('lev_sign.txt',[nodeNumbers,lev_sign]);

%% Example of exporting data to cytoscape
% Exporting things to cytoscape
degree = sum(A)';
csvwrite('degree.txt',[nodeNumbers,degree]);
% Import table with this data, change node styles to visualize
