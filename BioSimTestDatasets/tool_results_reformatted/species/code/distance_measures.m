cd .
h=dir('../temp/*_table');
outname=('../temp/l1_averageabundance_results.txt');
outname2=('../temp/l1_truthabundance_results.txt');
fid2=fopen(outname,'wt');
fprintf(fid2,'L1 Distance\t Sample Name\t Tool\n');
fid3=fopen(outname2,'wt');
fprintf(fid3,'L1 Distance\t Sample Name\t Tool\n');


truth_set_vector=[];
average_distances=[];
truth_distances=[];
average_distance_table={};
truth_distance_table={};
for i=1:numel(h)
    truth_index=0;
    fid = fopen(['../temp/' h(i).name]);
    datasetname=regexprep(h(i).name,'_table','');
    average_distance_table{1,(i-1)*2+1}=h(i).name;
    
    header=fgetl(fid);
    header_fields=regexp(header,'\t','split');
    numformat=[];
    for j=1:length(header_fields)-1
        numformat=[numformat '%f '];
    end
    C = textscan(fid, ['%f' numformat]);
    fclose(fid);
    thetable_old=cell2mat(C);
    
    %Check fidelity of data
    sum_to_one=sum(thetable_old);
    if sum(isnan(sum_to_one))>0
        disp(['Found NaN''s in: ' h(i).name ' -- overwrite with 0!'])
        thetable_old(isnan(thetable_old))=0;
    end
    disp('The Classifiers that do not match the sum-to-one criteria')
    binvector=round(sum_to_one*1000)~=1000;
    thetable=thetable_old(:,1);
    for j=2:length(binvector)
        if binvector(j)
            sprintf('%s %s %s %s%s','In', header_fields{j}, 'only', num2str(sum_to_one(j)*100), '% of the reads found')
        end
        %correct to be real distributions (normalize)
        %thetable(:,j)=thetable_old(:,j)/sum_to_one(j);    
    end
    % if not normalize
    thetable=thetable_old;
    
    for j=1:numel(header_fields)
       if j>1
            average_distance_table{j,(i-1)*2+1}=header_fields{j};
            truth_distance_table{j,(i-1)*2+1}=header_fields{j};
       end
       if(~isempty(regexp(header_fields{j},'TRUTH','match')))
           truth_index=j;
           truth_set_vector=[truth_set_vector;i];
           truth_abundance=thetable(:,j);
       end
    end
    
    if truth_index==0
        average_abundance=mean(thetable(:,2:end),2);
        %avg_temp=KLDiv(thetable(:,2:end)',average_abundance');
        avg_temp=pdist2_again(thetable(:,2:end)',average_abundance','L1');
        average_distances=[average_distances;avg_temp];
        for j=1:length(avg_temp)
            average_distance_table{j+1,i*2}=avg_temp(j);
             temp1=regexprep(header_fields{j+1},[datasetname '_'],'');
             toolname=regexprep(temp1,'.txt','');
             fprintf(fid2,'%f\t %s\t %s\n',avg_temp(j),datasetname,toolname);
        end
        
    else
        truth_distance_table{1,(i-1)*2+1}=h(i).name;
        if (~isempty(thetable(:,2:truth_index-1)))
            tempavg1=mean(thetable(:,2:truth_index-1),2);
        else
            tempavg1=zeros(size(thetable(:,2:truth_index-1)));
        end
        if (~isempty(thetable(:,truth_index+1:end)))
            tempavg2=mean(thetable(:,truth_index+1:end),2);
        else
            tempavg2=zeros(size(thetable(:,truth_index+1:end)));
        end
        
        average_abundance=mean([tempavg1 tempavg2],2);
        %avg_temp=KLDiv([thetable(:,2:truth_index-1) thetable(:,truth_index+1:end)]',average_abundance');
        avg_temp=pdist2_again([thetable(:,2:truth_index-1) thetable(:,truth_index+1:end)]',average_abundance','L1');
        average_distances=[average_distances;avg_temp];
        for j=1:length(avg_temp)
            average_distance_table{j+1,i*2}=avg_temp(j);
            if ~((j+1)==truth_index)
                temp1=regexprep(header_fields{j+1},[datasetname '_'],'');
                toolname=regexprep(temp1,'.txt','');
                fprintf(fid2,'%f\t %s\t %s\n',avg_temp(j),datasetname,toolname);
            end
        end
        %truth_temp=KLDiv([thetable(:,2:truth_index-1) thetable(:,truth_index+1:end)]',truth_abundance');
        truth_temp=pdist2_again([thetable(:,2:truth_index-1) thetable(:,truth_index+1:end)]',truth_abundance','L1');
        truth_distances=[truth_distances;truth_temp];
        for j=1:length(truth_temp)
            truth_distance_table{j+1,i*2}=truth_temp(j);
            temp1=regexprep(header_fields{j+1},[datasetname '_'],'');
            toolname=regexprep(temp1,'.txt','');
            fprintf(fid3,'%f\t %s\t %s\n',truth_temp(j),datasetname,toolname);
        end
    end
        
end

fclose(fid2);
fclose(fid3);
