function x = ktensor2(vec, mat)
n = size(mat, 2);
m = size(mat, 1);
x = zeros(m, m);

for i = 1:n
   mi = mat(:, i);
   p = mi * mi';
   x = x + vec(i) * p;
end
