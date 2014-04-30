function main_img_tensor(dataidx)
clearvars -except dataidx
global L N lr powereps name

if dataidx == 1
   filepath = '~/data/yale_img/Yale_32x32.mat';
   classnum = 15;
else
   filepath = '~/data/pie_img/PIE_32x32.mat';
   classnum = 68;
end
path('./tensor_toolbox_2.5/', path);
load(filepath);
L = 30.0;
N = 30.0;
lr = 1.0;
powereps = 1e-5;

[nSmp,nFea] = size(fea);
for i = 1:nSmp
     fea(i,:) = fea(i,:) ./ max(1e-12,norm(fea(i,:)));
end

maxValue = max(max(fea));
fea = fea/maxValue;

T = classnum;
mumap = tdgauss_img(fea, T);
fea2 = fea * mumap * diag(1./sqrt(sum(mumap.*mumap, 1)));
[C, result] = max(fea2, [], 2);
score = nmi(gnd, result)
