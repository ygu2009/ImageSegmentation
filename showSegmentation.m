function [seg, state]=showSegmentation(r_state, K, imdata, low_bkg, high_bkg, mycolor)

[rows, cols]=size(imdata);


n=0;
for k=1:K
    for j=1:cols
        for i=1:rows
            n=n+1;
            state(i,j,k)=r_state(n);
        end
    end
end

seg=zeros(rows,cols);
for i=1:rows
    for j=1:cols
        ind=find(state(i,j,:)==max(state(i,j,:)));
        seg(i,j)=ind(1);
        if imdata(i,j)<low_bkg | imdata(i,j)>high_bkg
            seg(i,j)=0;
        end
    end
end
