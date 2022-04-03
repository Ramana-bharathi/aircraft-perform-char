function D=drag(rho,v,S,C_Do,K,W)
% function to calculate drag force.
    D=0.5*rho.*v.^2*S*C_Do+(K*W^2./(0.5*rho.*v.^2*S));
end 
