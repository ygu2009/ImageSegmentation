clc
clear
close all

%% load a 2D CT image
load CTbag.mat
figure, imshow(imdata,[])

%% choose proper intensity ranges, removing background for finding the number of clusters/objects
y=imdata(:);
low_bkg=500;
high_bkg=4000;
y(find(y<low_bkg))=[];
y(find(y>high_bkg))=[];
bkg_num=length(find(imdata<low_bkg | imdata>high_bkg));

mycolor=lines;


%% Using kmeans clustering method
K=10;

[IDX,C] = kmeans(y,K,'emptyaction','singleton','MaxIter',500);
for i=1:K
    mu(i)=mean(y(find(IDX(:)==i),:));
    sigma(i)=(cov(y(find(IDX(:)==i),:)));
    pipi(i)=length(y(find(IDX(:)==i),:));
end
pipi=pipi/sum(pipi);

figure
bins=low_bkg:5:high_bkg;
[histxx, xout]=hist(y,bins);
hold on
plot(bins,histxx/trapz(xout,histxx),'b')

hold on
plot(mu(:), 0,'k*')
for i=1:K
    xx = bins;
    yy = pipi(i)*normpdf(xx, mu(i), (sigma(i))^(1/2));
    plot(xx, yy, 'r')
end
title(['KMeans: ', num2str(K), ' classes'])
box on
drawnow;

[rows, cols]=size(imdata);

%% Image segmentation based on results from Kmeans clustering method
mex Kmeans_2Dseg.c
[r_state]=Kmeans_2Dseg(imdata, mu, sigma, pipi, low_bkg, high_bkg);

[seg_Kmeans]=showSegmentation(r_state, K, imdata, low_bkg, high_bkg, mycolor);
figure, imshow(label2rgb(seg_Kmeans, mycolor)),title('Kmeans clustering');
drawnow;

%% using MAP EM
% initial priors
[N,D]=size(y);
alpha=0.001*ones(1,K);
beta_0=5*ones(1,K);
m0=mean(y)*ones(K,D);
S0=1*repmat(eye(D),[1 1 K]);
v0=3*ones(1,K);
%rubber sheet; bulk rubber; saline; clay
targets_mean=[1147 1253 1122 1544];
targets_var=[80590 41251 28678 71742]/5;

vb1400=length(find(y<1400));
va1400=length(find(y>=1400));

% adding known priors
for j=1:K
    for n=1:4
        dist(n)=mean((mu(j)-targets_mean(n)).^2);
    end
    [c ind]=min(dist);
    m0(j)=targets_mean(ind);
    if m0(j)<1400
        v0(j)=vb1400/25;
        beta_0(j)=5;
        S0(j)=targets_var(ind)*v0(j);
    else
        v0(j)=va1400/5;
        beta_0(j)=20;
        S0(j)=targets_var(ind)*v0(j);
    end
end

maxiter=1000;
[mu_mapEM, sigma_mapEM, pipi_mapEM]=MAPEM_1D(y, mu, sigma, pipi, maxiter, m0,  S0, v0, alpha, beta_0, bins);
mu=mu_mapEM;
sigma=sigma_mapEM;
pipi=pipi_mapEM;
ind=find(pipi<0.01);
mu(ind)=[];
sigma(ind)=[];
pipi(ind)=[];
K=length(mu);

%% MAP EM clustering process
[histxx, xout]=hist(y,bins);
figure, plot(bins,histxx/trapz(xout,histxx),'b')
hold on
for i=1:K
    xx = bins;
    yy = pipi(i)*normpdf(xx, mu(i), (sigma(i))^(1/2));
    plot(xx, yy, 'LineWidth',4,'Color',mycolor(i,:))
end
title(['MAPEM: ', num2str(K), ' classes'])
box on
drawnow;

%% Image segmentation based on results from MAPEM clustering
[r_state]=Kmeans_2Dseg(imdata, mu, sigma, pipi, low_bkg, high_bkg);

[seg_MAPEM, state]=showSegmentation(r_state, K, imdata, low_bkg, high_bkg, mycolor);
figure,imshow(label2rgb(seg_MAPEM, mycolor)),title('MAPEM clustering'),drawnow;

% initials segmentation from MAPEM clustering
fig=figure,
imshow(label2rgb(seg_MAPEM, mycolor)),title('MAPEM clustering');

winsize = get(fig,'Position');
winsize(1:2) = [0 0];
mm=1;
mov(mm) = getframe(fig, winsize);


% MAP Priors
D=2;
MAP_alpha=0.001*ones(1,K);
MAP_beta_0=10*ones(1,K);
MAP_m0=mean(y)*ones(K,D);
MAP_S0=10000*800*repmat(eye(D),[1 1 K]);
MAP_v0=800*ones(1,K);

% add known priors
% post-processing the kmeans results (find the initial cluster that is
% targets mode and variance:
%rubber sheet; bulk rubber; saline; clay
targets_mean=[1150 1250 1120 1540];
targets_var=[90000 50000 30000 70000];

vtotal=length(y);
vb1400=length(find(y<1500));
va1400=length(find(y>=1400));

for n=1:4
    for j=1:K
        dist(j)=mean((mu(j)-targets_mean(n)).^2);
    end
    [c ind]=min(dist);
    MAP_m0(ind)=targets_mean(n);
    if MAP_m0(ind)<1500
        MAP_v0(ind)=vb1400/3;
        MAP_beta_0(ind)=10;
        MAP_S0(ind)=targets_var(n)*MAP_v0(ind);
    else
        MAP_v0(ind)=va1400/5;
        MAP_beta_0(ind)=20;
        MAP_S0(ind)=targets_var(n)*MAP_v0(ind);
    end
end


% % Regular EM initials 
% MAP_alpha=1*ones(1,K);
% MAP_beta_0=0*ones(1,K);
% MAP_m0=0*ones(K,D);
% MAP_S0=0*repmat(eye(D),[1 1 K]);
% MAP_v0=(-D-2)*ones(1,K);


%% MAP EM MRF segmentation process
mex MAPEM_MRF_2Dseg.c
flag=(-1)*ones(1,K);
iter=0;
diff_Q=1;
MRF_beta=1;
while iter<100 %& diff_Q>1e-10
    
    iter=iter+1
    
    %%% EM step
    [r_state, mu_hat, sigma_hat, pipi_hat, nc]=MAPEM_MRF_2Dseg(imdata, mu, sigma, pipi, state, flag, MRF_beta, MAP_alpha, MAP_beta_0, MAP_m0, MAP_S0, MAP_v0, low_bkg, high_bkg, bkg_num); 

    %deal with extinguished clusters
    for j=1:K
        if nc(j)<=10
            flag(j)=1;
        end
    end

    mu_hat=mu_hat';
    sigma_hat=sigma_hat';
    pipi_hat=pipi_hat';
    
    % convergence
    Q=0;
    for j=1:K
        if flag(j)==-1
            Q=Q+abs(mu(j)-mu_hat(j))/abs(mu(j))+abs(sigma(j)-sigma_hat(j))/abs(sigma(j))+abs(pipi(j)-pipi_hat(j))/abs(pipi(j));
        end
    end
    Q
    
    
    % update mu, sigma, pi
    mu=mu_hat;
    sigma=sigma_hat;
    pipi=pipi_hat;
    
    [seg, state]=showSegmentation(r_state, K, imdata, low_bkg, high_bkg, mycolor);
    
    hold off
    imshow(label2rgb(seg, mycolor))
    title(strcat('iter:',num2str(iter)))
    drawnow;
    
    if iter<10 || mod(iter, 5)==0
        mm=mm+1;
        mov(mm) = getframe(fig, winsize);
    end

    
end

movie2gif(mov(1:mm), 'MAPEMMRF.gif')
