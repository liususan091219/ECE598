function main_img_em(dataidx)
clearvars -except dataidx

if dataidx == 1
   filepath = '~/data/yale_img/Yale_32x32.mat';
   classnum = 15;
else
   filepath = '~/data/pie_img/PIE_32x32.mat';
   classnum = 68;
end

load(filepath);
path('./emgm/', path);
[result, model, llh] = emgm(fea', 68);
score = nmi(gnd, result)
