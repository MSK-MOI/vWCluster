%% Install CVX: http://cvxr.com/cvx/doc/install.html prior to running this function.
function distanceM= dist_cvx(drho,D,F,gamma,m,n,c,d)
% Solve individual vector OMT problem via CVX

% Jiening Zhu, 06/19/2021

cvx_solver sdpt3
cvx_begin quiet
    variables u(c*m,1) p(d*n,1)
    minimize (norm(u,1)+gamma*norm(p,1))
    subject to
        drho+D*u+F*p == 0;
cvx_end

    distanceM=cvx_optval;