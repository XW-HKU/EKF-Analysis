function [ out ] = A( u )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
out=eye(3)+(1-cos(norm(u)))/(norm(u)^2)*cross_mat(u)+(1-sin(norm(u))/norm(u))*(cross_mat(u)^2)/(norm(u)^2);

end

