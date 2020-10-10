% Leverage centrality for all nodes of the network
% Input Adjacency matrix
% by Nidhi Desai
function [lev_vec] = leverage_centrality(A)
G = graph(A);
degree_vec = degree(G);
lev_vec = zeros(size(A,1),1);
for i=1:size(A,1)
    neigh = neighbors(G,i);
    for j=1:size(neigh,1)
    lev_vec(i) = (degree_vec(i)-degree_vec(neigh(j)))/(degree_vec(i)+degree_vec(neigh(j)));
    end
    lev_vec(i) = lev_vec(i)/degree_vec(i);
end

