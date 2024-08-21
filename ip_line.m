function [x,w] = ip_line(N)
%  http://www.comlab.ox.ac.uk/people/nick.trefethen/gauss.m
  beta = .5./sqrt(1-(2*(1:N-1)).^(-2));
  T = diag(beta,1) + diag(beta,-1);
  [V,D] = eig(T);
  x = diag(D); [x,i] = sort(x);
  w = 2*V(1,i).^2;
