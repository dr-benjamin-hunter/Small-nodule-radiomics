
T = readtable('/Users/ca1117/Documents/Lung_docs/train/results.csv');
name_t = T(:,1:4);
num_t = T(:,5:end);
t=table2array(num_t);
t_t=t';
x_m = zeros(3,666,422);
x_m(1,:,:)=t_t(:,1:3:end);
x_m(2,:,:)=t_t(:,2:3:end);
x_m(3,:,:)=t_t(:,3:3:end);


z = readtable('/Users/ca1117/Documents/Lung_docs/train/Patientss_train.csv');

ben = 0;
mal = 0;
ben_i=zeros(1,422);
mal_i=zeros(1,422);

for i=1:422
    if z{i,6}==0
        ben=ben+1;
        ben_i(i)=i;

    elseif z{i,6}==1
        mal=mal+1;
        mal_i(i)=i;

    end
end
ben_i(ben_i==0)=[];
mal_i(mal_i==0)=[];

map_b=x_m(:,:,ben_i(:));
map_m=x_m(:,:,mal_i(:));
temp(:,:)= map_m(1,:,:);
