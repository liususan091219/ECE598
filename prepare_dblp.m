% prepare data
% folder is given
% stem - whether to stem the plural forms; only valid when filter = 1
% filter - whether to filter the network
global folder

[pu, pumat] = ReadEdge([folder 'pt.txt']);

%[pa, pamat] = ReadEdge([folder '/PA.txt']); % paper-author
%[pc, pcmat] = ReadEdge([folder '/PC.txt']); % paper-conference

uname = ReadName([folder 'term.txt']);

%aname = ReadName([folder '/author.txt']);
%cname = ReadName([folder '/conf.txt']);
uname = containers.Map(uname{1},uname{2});
%aname = containers.Map(aname{1},aname{2});
%cname = containers.Map(cname{1},cname{2});
round = 0;
while filter
    round = round + 1
    stable = true;
    deg = sum(pumat); %degree of term in PT network
    uremain = find(deg>1);  % remove terms with only one occurrence
    stable = stable & length(uremain)==length(deg);
    name = values(uname);
    uname = containers.Map(1:length(uremain),name(uremain));
    pumat = pumat(:,uremain);
    uold2new = zeros(length(deg),1);
    uold2new(uremain) = 1:length(uremain);
    % merge single and plural forms into one word
    if stem
        pp = PluralPair(uname);
        name = values(uname);
        deg = sum(pumat);
        [~,ind]=sort(deg(pp),2);
        for i=1:size(pp,1)
            pp(i,:) = pp(i,ind(i,:));
        end
        pumat(:,pp(:,2)) = pumat(:,pp(:,2)) | pumat(:,pp(:,1));
        toreplace=zeros(size(uold2new));
        toreplace(uremain(pp(:,1)))=uremain(pp(:,2));
        ind = toreplace(pu(:,2));
        pu(ind>0,2) = nonzeros(ind);
        remain = setdiff(1:length(uremain),pp(:,1));
        uname = containers.Map(1:length(remain),name(remain));
        pumat = pumat(:,remain);
        uold2new(uremain(remain))=1:length(remain);
        uold2new(uremain(pp(:,1)))=uold2new(uremain(pp(:,2)));
        uremain = uremain(remain);
        stem = false;
    end

    deg = sum(pamat); %degree of author in PA network
    aremain = find(deg>1);
    stable = stable & length(aremain)==length(deg);
    name = values(aname);
    aname = containers.Map(1:length(aremain),name(aremain));
    pamat = pamat(:,aremain);
    aold2new = zeros(length(deg),1);
    aold2new(aremain) = 1:length(aremain);

    deg = sum(pcmat);
    cremain = find(deg>100); % conferences with more than 100 papers
    stable = stable & length(cremain)==length(deg);
    name = values(cname);
    cname = containers.Map(1:length(cremain),name(cremain));
    pcmat = pcmat(:,cremain);
    cold2new = zeros(length(deg),1);
    cold2new(cremain) = 1:length(cremain);

    deg1 = sum(pumat,2);
    deg2 = sum(pamat,2);
    deg3 = sum(pcmat,2);
    premain = deg1>2 & deg2>0 & deg3>0; 
    % more than 2 terms and at least one author, one conference
    remain = find(premain);
    stable = stable & length(remain)==length(premain);
    pumat=pumat(remain,:); % remove incomplete papers
    pamat=pamat(remain,:);
    pcmat=pcmat(remain,:);
    pold2new = zeros(length(premain),1);
    pold2new(remain) = 1:length(remain);

    pu(:,1)=pold2new(pu(:,1));
    pu(:,2)=uold2new(pu(:,2));
    pu = pu(pu(:,1) & pu(:,2),:);
    pa(:,1)=pold2new(pa(:,1));
    pa(:,2)=aold2new(pa(:,2));
    pa = pa(pa(:,1) & pa(:,2),:);
    pc(:,1)=pold2new(pc(:,1));
    pc(:,2)=cold2new(pc(:,2));
    pc = pc(pc(:,1) & pc(:,2),:);

    if stable 
        break;
    end
end
size(pumat)
if filter    
    f = fopen([folder '/stemmed/pt.txt'],'w');
    fprintf(f,'%d\t%d\n',pu(:,1:2)');
    fclose(f);
    f = fopen([folder '/stemmed/pa.txt'],'w');
    fprintf(f,'%d\t%d\n',pa(:,1:2)');
    fclose(f);
    f = fopen([folder '/stemmed/pc.txt'],'w');
    fprintf(f,'%d\t%d\n',pc(:,1:2)');
    fclose(f);
    WriteDict([folder '/stemmed/term.txt'],[keys(uname);values(uname)]);
    WriteDict([folder '/stemmed/author.txt'],[keys(aname);values(aname)]);
    WriteDict([folder '/stemmed/conf.txt'],[keys(cname);values(cname)]);
end
