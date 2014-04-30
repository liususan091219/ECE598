function mumap = tdgauss_img(dwmat,T)
global L N lr powereps

dwmat = bsxfun(@rdivide, dwmat, sum(dwmat, 2));
eye3 = tensor(sptensor([1 1 1], 1));
doclen = sum(dwmat, 2);
[D,W] = size(dwmat)
doclen = round(doclen * 1000)/1000;
wdmat = dwmat';
doclen2 = 1./doclen./doclen;
doclen3 = doclen2 ./doclen;
doclen = diag(1./doclen);

E1 = sum(doclen * dwmat)'/D; 

punorm2 = doclen * dwmat;
E2 = ktensor2(doclen2, wdmat) / D;
disp 'computing first eigen decomposition...'
[U01, D01] = eig(E2 - E1 * E1');
disp 'finished computing first eigen decomposition.'
[sigma, idx] = min(diag(D01));
vector = U01(:, idx);
vtxmex = vector' * bsxfun(@minus, wdmat, E1);
M1 = sum(wdmat * diag(vtxmex .* vtxmex), 2)/D;

M2 = E2 - sigma * eye(W, W);
[U1, D1] = eigs(M2, T, 'la');
D1 = diag(D1);

if any(any(D1 <= 0))
   error('D1 is negative');
end

W1 = diag(1./sqrt(D1)) * U1';
tdmat = W1 * wdmat;
WM1 = W1 * M1;
x2M1 = squeeze(ttt(tensor(W1 * W1'), tensor(WM1)));

E3 = tensor(ktensor3(doclen3, tdmat));
%tenmatrix = tenmat(E3, [1 2], 3);
%any(any(isnan(tenmatrix.data)))

x11 = x2M1 + permute(x2M1, [1 3 2]) + permute(x2M1, [3 1 2]);

M3 = E3/D - x11;
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

function M = normalize_cols(M)
sp = sum(max(M,0),1);
sn = sum(max(-M,0),1);
M = max(M*diag(2*double(sp>sn)-1),0);
M = bsxfun(@rdivide, M, sum(M, 1));
