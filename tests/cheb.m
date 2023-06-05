function [D,x] = cheb(N)

if N == 0
    D=0;
    x=1;
    warning('You''re kidding, right?')
    return
end

a = (0:N)';

x = cos(pi*a/N);
c = [2 ; ones(N-1,1) ; 2] .*(-1).^a;
X = repmat(x,1,N+1);
dX = X-X';
D  = (c*(1./c)')./(dX+(eye(N+1)));                                          % off-diagonal entries
D  = D - diag(sum(D,2));                                                    % diagonal entries

end