function writeflattopics(filename, mumap, topN)

global name

fid = fopen(filename, 'w');

childnum = size(mumap, 1);

for i = 1:childnum
    [~,ind] = sort(mumap(i, :), 'descend');
    ind = reshape(ind, max(size(ind)), 1);
    childtopics = name(ind(1:min(topN,size(ind,1)))');
    topicsize = max(size(childtopics));
    childtopics = reshape(childtopics, topicsize, 1);
    for j = 1:topicsize
        fprintf(fid, '%s ', childtopics{j, 1});
    end
    fprintf(fid, '\n');
end

fclose(fid);
