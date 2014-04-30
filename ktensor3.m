function x = ktensor3(vec,mat)
% Chi Wang's own implementation of tensor(ktensor(vec,mat,mat,mat))
%any(any(~isreal(vec)))
%any(any(~isreal(mat)))
n = size(mat,2);
m = size(mat,1);
x = zeros(m*m*m,1);
for i=1:n
    mi = mat(:,i);
    p = mi * mi';
    p = p(:) * mi';
    %oldx = x;
    x = x + vec(i)*p(:);
    %if any(isnan(x))
    %   fprintf('i=%d, n=%d\n', i, n);
    %   x
    %   oldx
    %   vec(i)
    %   p(:)
    %   error('sss');
    %end
end
%any(isnan(x))
x = reshape(x,[m m m]);

