% matlab_res = [];
% for n = 15:25
%     U = randU(n);
%     f = @() ryser(U);
%     matlab_res = [matlab_res, timeit(f)]
% end

fileName = 'data_matlab.txt';
fid = fopen(fileName, 'w');
for i = 1:length(matlab_res)
    fprintf(fid, '%f ', matlab_res(i));
end
fclose(fid);

function U = randU(n)
X = (randn(n) + i*randn(n))/sqrt(2);
[Q,R] = qr(X);
R = diag(diag(R)./abs(diag(R)));
U = Q*R;
end
