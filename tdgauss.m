function mumap = tdgauss(dwmat,T)
global L N lr powereps

eye3 = tensor(sptensor([1 1 1], 1));
doclen = sum(dwmat,2);
ind = doclen == 0;
dwmat(ind,:)=[];   % filter documents with length<3
doclen(ind,:)=[];
[D,W] = size(dwmat)
doclen3 = 1./doclen./doclen./doclen;
doclen = sparse(1:D, 1:D, 1./doclen);
wdmat = dwmat';

M1 = full(sum(doclen*dwmat)')/D; 
punorm2 = doclen * dwmat;

E2 = punorm2' * punorm2 / D;
% find eigen values of E[xx^T] - E[x]xE[x]^T
[U00, D00] = eigs(E2, T, 'la');
x = U00' * M1;
DD = D00 - x * x';
[U01, D01] = eig(DD);
D0 = D01;
%[U0, D0] = eigs(E2 - M1 * M1', T, 'la');
sigma = min(diag(D0));

M2 = E2 - sigma * sparse(1:W, 1:W, ones(W,1));
[U1, D1] = eigs(M2, T, 'la');
D1 = diag(D1);

W1 = diag(1./sqrt(abs(D1))) * U1';
tdmat = W1 * wdmat;
WM1 = W1 * M1;
x2M1 = squeeze(ttt(tensor(W1 * W1'), tensor(WM1)));

E3 = tensor(ktensor3(doclen3, tdmat));
%tenmatrix = tenmat(E3, [1 2], 3);
%any(any(isnan(tenmatrix.data)))

x11 = x2M1 + permute(x2M1, [1 3 2]) + permute(x2M1, [3 1 2]);

M3 = E3/D - sigma * x11;
ws = zeros(1, T);
thetas = zeros(T, T);
for kk = 1:T
   kk
   maxlambda = -1000;
   thetastar = 0;
   tao = 0;
   while true
      tao = tao + 1; 
      if tao > L
         break;
      end
      theta = W1 * wdmat(:, ceil(rand()*D));
      theta = theta/norm(theta);
      tmplambda = -inf;
      for t = 1:N
        oldlambda = tmplambda; 
        update = ttsv(M3, theta, -1);
        tmplambda = dot(theta, update);
        delta = abs(tmplambda - oldlambda);
        converge = delta < powereps * oldlambda;
        theta = theta + lr * update;
        theta = theta ./norm(theta);
      end
      if tmplambda > maxlambda
         maxlambda = tmplambda;
         thetastar = theta;
      end
   end
   theta = thetastar;
   lambda = maxlambda;
   ws(kk) = (1/lambda)^2; 
   thetas(:, kk) = theta;
   M3 = M3 - lambda * ttm(eye3, {theta, theta, theta});
end

mumap = normalize_cols(U1 * (diag(sqrt(D1))*thetas));
mumap = bsxfun(@rdivide, mumap, sum(mumap, 1));

function M = normalize_cols(M)
sp = sum(max(M,0),1);
sn = sum(max(-M,0),1);
M = max(M*diag(2*double(sp>sn)-1),0);
M = bsxfun(@rdivide, M, sum(M, 1));
