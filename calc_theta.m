function theta=calc_theta(lambda1,X,tau)
theta=4*sin(lambda1)*cos(lambda1*X)*exp(-(lambda1*lambda1)*tau)/((2*lambda1)+sin(2*lambda1));
end