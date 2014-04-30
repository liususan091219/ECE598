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

T = classnum;
mumap = tdgauss_img(fea, T);
fea2 = fea * mumap;
[C, result] = max(fea2, [], 2);
score = nmi(gnd, result)
