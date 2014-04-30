function main(folder, dataname)

global L N lr powereps name

path('./tensor_toolbox_2.5/', path);
path('./emgm/', path);
stem = 0;
filter = 0;
prepare_dblp;
emptyidx = find(sum(pumat, 2) == 0);
pumat(emptyidx, :) = [];
label = load([folder 'all.label']);
label(emptyidx, :) = [];
[docnum, wordnum] = size(pumat);
fprintf('docnum %d, wordnum %d', docnum, wordnum);
%testnum = round(docnum / 5);
%y = randomsample(docnum, testnum); 
%testmat = pumat(y, :);
%trainmat = pumat;
%trainmat(y, :) = [];

name = values(uname);
L = 30;
N = 30;
lr = 1.0;
powereps = 1e-5;

T = 20;
%mumap = tdgauss(pumat, T);
%mumap = mumap';

%writeflattopics('tree', mumap, 20);
%error(' ');
%pumat2 = pumat * mumap;
%[C, result] = max(pumat2, [], 2);
%
%label = load([folder 'all.label']);
%score = nmi(label, result);
%
%fprintf('score=%f\n', score); 
randommat = rand(wordnum, 40);
randommat = bsxfun(@rdivide, randommat, sqrt(sum(randommat.*randommat,1)));
sm_pumat = pumat * randommat;

sm_pumat = bsxfun(@rdivide, sm_pumat, sqrt(sum(sm_pumat .*sm_pumat, 2)));
[result, model, llh] = emgm(sm_pumat', T);
score = nmi(label, result)

