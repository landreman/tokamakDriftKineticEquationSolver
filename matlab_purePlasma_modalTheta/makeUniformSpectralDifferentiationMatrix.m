function [x, D, NFinal] = makeUniformSpectralDifferentiationMatrix(N, xMin, xMax)

if mod(N,2)==1
    fprintf('Warning: N must be even. Adding 1.\n')
    N=N+1;
end
NFinal=N;

Delta = xMax-xMin;

h=2*pi/N;
x=((1:N)*h)/(2*pi)*Delta+xMin;
col=[0, 0.5*(-1).^(1:N-1).*cot((1:N-1)*h/2)]';  
D = toeplitz(col, col([1, N:-1:2]))*2*pi/Delta;

end