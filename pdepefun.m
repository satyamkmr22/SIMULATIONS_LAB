function [c, f, s] = pdepefun(x, t, u, dudx)
    c = 1;
    f = dudx;
    s = 0;
end
% u represents theta and t represents tau