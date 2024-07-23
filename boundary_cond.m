function [pl, ql, pr, qr] = boundary_cond(xl, ul, xr, ur, t)
    pl = 0;
    ql = 1;
    pr = 10 * ur;
    qr = 1;
end