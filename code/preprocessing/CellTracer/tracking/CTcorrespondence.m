function result = CTcorrespondence(method,label1,label2,nb1,nb2,d1,d2,minratio,nneighbor)
%create initial working variables
[m n] = size(label1);
n1 = max(label2(:));
m1 = max(label1(:));
denom1 = zeros(m1,n1); denom2 = zeros(n1,m1);
for i = 1:m1
    [r c] =find(label1 == i);
    denom1(i,:) = length(r);
end
for i = 1:n1
    [r c] =find(label2 == i);
    denom2(i,:) = length(r);
end

d = d1 * d2;
container = zeros(m+2* d, n+ 2 * d);
t1 = container; t1(d+1:d+m, d+1:d+n) = label1; 

[c1 dummy] = size(nb1);
[c2 dummy] = size(nb2);
result = zeros(c1,c2,4);
best = zeros(c1,c2);
bestk = zeros(c1,1);
bestp = zeros(c2,1);

for k = 1:c1
    indexk{k} = find(nb1(k,:)>0);
end 
for p = 1:c2
    indexp{p}= find(nb2(p,:) > 0);
end
best = zeros(c1,c2);
for i = 1:2 * d1+1
    row = (i-1) * d2 + 1;
    for j = 1:2* d1+1
        col = (j - 1 ) * d2 + 1;
        t2 = container; t2(row:row + m-1, col:col+n-1) = label2;
        bw = t1 & t2;
        num1 = zeros(m1,n1);
        r = t1(bw);c= t2(bw);
        for k=1:length(r)
            num1(r(k),c(k)) = num1(r(k),c(k))+1; 
        end
        mat1 = num1 ./ denom1;
        mat2 = num1' ./ denom2;
        mat1(mat1<minratio) = 0;
        mat2(mat2<minratio) = 0;
        
        if method == -3
            score = (mat1 + mat2') .* (mat1 + mat2');
            
            for k = 1:c1
                for p = 1:c2
                    if mat1(k,p) >= minratio && mat2(p,k) >= minratio
                        submat = score(indexk{k},indexp{p});
                        temp = sum(submat(:));
                        if temp > best(k,p)
                            best(k,p) = temp;
                            result(k,p,1) = mat1(k,p);
                            result(k,p,2) = mat2(p,k);
                            result(k,p,3) = row;
                            result(k,p,4) = col;
                        end
                    end
                end
            end
        elseif method == -2
            for k = 1:c1
                for p = 1:c2
                    submat1 = mat1(indexk{k},indexp{p});
                    [v iv] = sort(submat1,2,'descend');
                    submat1(submat1 < repmat(v(:,2),1,length(indexp{p}))) = 0;
                    record1 = result(indexk{k},indexp{p},1);
                    newf = max(record1,submat1); 
                    cf = newf > record1;
                    
                    submat2 = mat2(indexp{p},indexk{k});                    
                    [v iv] = sort(submat2,2,'descend');
                    submat2(submat2 < repmat(v(:,1),1,length(indexk{k}))) = 0;
                    record2 = result(indexk{k},indexp{p},2)';
                    newb = max(record2,submat2);
                    cb = newb > record2;
                    if length(find(max(cb')>0)) >= min(length(indexp{p}),nneighbor) & length(find(max(cf')>0)) >= min(length(indexk{k}),nneighbor)
                        result(indexk{k},indexp{p},1) = newf;
                        result(indexk{k},indexp{p},2) = newb';
                    end
                end
            end
        elseif method == -1
            [v iv] = sort(mat1,2,'descend');
            mat1(mat1 < repmat(v(:,2),1,n1)) = 0;
            [v iv] = sort(mat2,2,'descend');
            mat2(mat2 < repmat(v(:,1),1,m1)) = 0;
            newf = max(result(:,:,1),mat1); 
            newb = max(result(:,:,2)',mat2);
            cf = newf > result(:,:,1);
            cb = newb > result(:,:,2)';
            
            if (length(find(max(cf')>0)) >= min(m1,nneighbor) && length(find(max(cb')>0)) >= min(n1,nneighbor)) 
                result(:,:,1) = newf;
                result(:,:,2) = newb';
            end 
        elseif method == 0
            score = (mat1 + mat2') .* (mat1 + mat2');
            for k = 1:c1
                for p = 1:c2
                    if mat1(k,p) >= minratio && mat2(p,k) >= minratio
                        submat = score(indexk{k},indexp{p});
                        submat = submat(:);
                        submat = sort(submat,'descend');
                        lastindex = max(length(indexp{p}),length(indexk{k}));
                        if score(k,p) >=submat(lastindex)
                            temp = sqrt(sum(submat(1:lastindex)) / lastindex)/2; 
                            if temp > best(k,p) 
                                best(k,p) = temp;
                                result(k,p,3) = row;
                                result(k,p,4) = col;
                                 result(k,p,1) = (temp + mat1(k,p)) / 2;
                                 result(k,p,2) = (temp + mat2(p,k)) / 2;
                            end
                        end
                    end
                end
            end
        elseif method == 1
            score = mat2';
            for p = 1:c2
                submat = score(:,indexp{p}); 
                submat = sort(submat,1,'descend');
                submat = submat(1,:); 
                temp = sum(submat(:));
                if temp > best(1,p)
                    best(1,p) = temp;
                    result(:,p,3) = row;
                    result(:,p,4) = col;
                    result(:,p,1) = mat1(:,p);
                    result(:,p,2) = mat2(p,:);
                end
            end
        elseif method == 2
            for p = 1:c2
                submat = mat2(indexp{p},:);
                temp = sum(submat(:));
                if temp > bestp(p)
                    bestp(p) = temp;
                    result(:,p,2) = mat2(p,:)';
                end
            end
            for k = 1:c1
                submat = mat1(indexk{k},:);
                temp = sum(submat(:));
                if temp > bestk(k)
                    bestk(k) = temp;
                    result(k,:,1) = mat1(k,:);
                end
            end
        elseif method == 3
            score = (mat1 + mat2') .* (mat1 + mat2');
            for k = 1:c1
                for p = 1:c2
                    if mat1(k,p) >= minratio && mat2(p,k) >= minratio
                        submat = score(indexk{k},indexp{p});
                        submat = submat(:);
                        submat = sort(submat,'descend');
                        lastindex = max(length(indexp{p}),length(indexk{k}));
                        if score(k,p) >=submat(lastindex)
                            temp = sqrt(sum(submat(1:lastindex)) / lastindex)/2; 
                            if temp > best(k,p) 
                                best(k,p) = temp;
                                result(k,p,3) = row;
                                result(k,p,4) = col;
                                result(k,p,1) = mat1(k,p);
                                result(k,p,2) = mat2(p,k);
                             end
                        end
                    end
                end
            end
        end
    end
end
