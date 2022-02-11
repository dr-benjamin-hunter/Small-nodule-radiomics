
ind0=find(map_b==0);
ind00=find(map_m==0);
map_b(map_b==0)=nan;
map_m(map_m==0)=nan;


ind1 = find(map_b==inf);
if isempty(ind1)
    [x1, y1, z1] = ind2sub(size(map_b), ind1);
    s1 = length(x1);
    for p = 1:s1
        map_b(x1(p),y1(p),z1(p))=nan;
    end
end

ind2 = find(map_b==-inf);
if isempty(ind2)
    [x1, y1, z1] = ind2sub(size(map_b), ind2);
    s1 = length(x1);
    for p = 1:s1
     map_b(x1(p),y1(p),z1(p))=nan;
    end
end

ind3 = find(map_m==inf);
if isempty(ind3)
    [x2, y2, z2] = ind2sub(size(map_m), ind3);
    s2 = length(x2);
    for pp = 1:s2
     map_m(x2(pp),y2(pp),z2(pp))=nan;
    end
end

ind4 = find(map_m==-inf);
if isempty(ind4)
    [x2, y2, z2] = ind2sub(size(map_m), ind4);
    s2 = length(x2);
    for pp = 1:s2
        map_m(x2(pp),y2(pp),z2(pp))=nan;
    end
end



map = cat(3,map_b,map_m);


starmap_b = nan(size(map_b));
starmap_m = nan(size(map_m));

mean_v = nan(size(map_b,1),size(map_b,2));
std_v = mean_v;
mean_i = nan(size(map_b,1),size(map_b,2));
std_i= mean_i;

for i=1:size(map_b,1)
    for j=1:size(map_b,2)
        data(:)=map(i,j,:);
        mean_v(i,j)=mean(data(:));     
        std_v(i,j)=std(data(:));
    end
end

for l=1:size(map_b,3)
starmap_b(:,:,l) = (map_b(:,:,l) - mean_v)./(std_v);
end

for k=1:size(map_m,3)
starmap_m(:,:,k) = (map_m(:,:,k) - mean_v)./(std_v);
end





n = size(map_b,3);

TEST_f = cat(3, starmap_b, starmap_m);
temp=TEST_f(3,666,:);
% Wilcoxon Rank sum test (CN vs AD)
TEST_f(isnan(TEST_f))=0;
for j = 1:size(TEST_f,1)
    for i = 1:size(TEST_f,2)
        A(:) = TEST_f(j, i, :);
        if ~isnan(A(1,1))
            bn=A(1:n);
            mn=A((n+1):end);            
            [p,~]=ranksum(bn, mn);
            P(j,i)=p;
        else
            P(j,i)=nan;
        end
    end
end


FDR64tot = nan(size(P));
for i=1:size(P,1)
    prow(:) = P(i,:);
    if isempty(prow(1,1))==0
        [~, ~, ~, FDR64tot(i,:)]=fdr_bh(prow,0.05,'pdep','no');
    else
        FDR64tot(i,:) = nan;
    end
end

fdr=FDR64tot;
fdr(fdr>0.05)=nan;
fdr=fdr';
          
indf=find(~isnan(fdr));
[x, y]=ind2sub(size(fdr),indf);

FM2=nan(size(TEST_f,3),length(y));

mask = nan(1,length(y));

for j = 1:size(TEST_f,3)
    subj = TEST_f(:,:,j);
    subj2=subj';
    for i=1:length(x)
        mask(1,i) =subj2(x(i), y(i));
    end
    FM2(j,:) = mask;
end
fdr2=fdr';
FM2_f=nan(size(TEST_f));
for i =1:3
    for j=1:666
        if ~isnan(fdr2(i,j))
            for k=1:422
                FM2_f(i,j,k)=TEST_f(i,j,k);
            end
        end
    end
end

%LASSO
score=ones(size(TEST_f,3),1);
score(1:n)=0;
FM2(isnan(FM2))=0;

% [B, FitInfo] = lassoglm(FM2, score, 'binomial','CV',10,'Standardize',false);

load('FitInfo_stad_split.mat');
load('B_stad_split.mat');
RPVa=nan(1,size(TEST_f,3));
for i=1:size(TEST_f,3)
    testt = FM2(i,:);
    index = FitInfo.IndexMinDeviance;
    notredindex = find(B(:,index));
    b = B(:,index);
    b(b==0)=[];
    feat=nan(length(notredindex),1);
    for k=1:length(notredindex)
        feat(k) = testt(notredindex(k));
    end
    feat=feat(:);
    b = b(:);
    prod = feat.*b;
    a = find(isnan(prod));
    if ~isempty(a)
        prod(isnan(prod))=[];
    end
    som = sum(prod);
   
    RPVa(1,i)=som;
end
rpv=RPVa';