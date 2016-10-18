function [mu, sigma, PI]=MAP_EM_1D(y, mu, sigma, PI, maxiter, m0, S0, v0, alpha, beta_0, bins)

% initials from the results of kmeans clustering
figure
[histxx, xout]=hist(y,bins);
hold on
plot(bins,histxx/trapz(xout,histxx),'b')

K=length(mu);
plot(mu(:), 0,'k*')
for i=1:K
    xx = bins;
    yy = PI(i)*normpdf(xx, mu(i), (sigma(i))^(1/2));
    plot(xx, yy, 'r')
end
title(['KMeans: ', num2str(K), ' classes'])
box on
drawnow;

K=length(mu);
[N,D]=size(y);

% initial likelihood
% Q0=0;
for j=1:K
    t=bsxfun(@minus, y, mu(j));
    p_y(:,j)=(1/((2*pi)*sqrt(det(sigma(j)))))*exp(-0.5*(bsxfun(@times, t*inv(sigma(j)), t)));
end

Q0=sum(log(p_y*PI'));

diff_Q=1;
iter=0;

flag=zeros(1,K);


% figure
while iter<maxiter & diff_Q>1e-15
    
    iter=iter+1;
    
    %%% E step
    
    PI_x=PI;
    ind=find(flag==1);
    PI_x(ind)=0;
    
    r=(p_y.*repmat(PI_x,N,1))./repmat(p_y*PI_x',1,K);  %unchanged E step (page. 356)
    
    
    %%% M step
    
    % estimate mu
    for j=1:K
        if flag(j)==0
            y_bar(j,:)=sum(bsxfun(@times, y, r(:,j)))/sum(r(:,j)); %(11.45) from Kevin Murphy
            mu_hat(j)=(sum(r(:,j))*y_bar(j,:)+beta_0(j)*m0(j))/(sum(r(:,j))+beta_0(j)); %(11.43) from Kevin Murphy
        end
    end
    
    
    % estimate sigma
    for j=1:K
        if flag(j)==0
            y_t=bsxfun(@minus, y, y_bar(j,:));
            S(:,:,j)=bsxfun(@times, y_t, r(:,j))'*y_t; %(11.47)  from Kevin Murphy
            Stemp(:,:,j)=(sum(r(:,j))*beta_0(j)/(sum(r(:,j))+beta_0(j)))*((y_bar(j,:)-m0(j,:))'*(y_bar(j,:)-m0(j,:)));
            sigma_hat(j)=(S0(:,:,j)+S(:,:,j)+Stemp(:,:,j))/(v0(j)+sum(r(:,j))+D+2); %(11.46) from Kevin Murphy
        end
    end
    
    % estimate pi
    K_x=K-sum(flag);
    pi_hat=(sum(r)+alpha-1)/(N+sum(alpha)-K_x);  %(11.41) from kevin murphy
    
    %deal with extinguished clusters
    for j=1:K
        if sum(r(:,j))<=0.001 
            flag(j)=1;
        end
    end
    
    % convergence ?
    
    Q=sum(abs(mu-mu_hat)./abs(mu))+sum(abs(PI-pi_hat)./abs(PI))+sum(abs(sigma-sigma_hat)./abs(sigma));
    diff_Q=Q/(3*K_x);
    
    % update mu, sigma, pi
    mu=mu_hat;
    PI=pi_hat;
    sigma=sigma_hat;
    
    for j=1:K
        if flag(j)==0
            t=bsxfun(@minus, y, mu(j));
            p_y(:,j)=(1/((2*pi)*sqrt(det(sigma(j)))))*exp(-0.5*(bsxfun(@times, t*inv(sigma(j)), t)));
        end
    end
    
    % %     Q=sum(log(p_y*PI'));
    % %
    % %     diff_Q=abs(Q0-Q);
    
    % %     Q0=Q;
    
    
        hold off
        plot(bins,histxx/trapz(xout,histxx),'b')
        hold on
        for i=1:K
            xx = 500:1:3000;
            yy = PI(i)*normpdf(xx, mu(i), (sigma(i))^(1/2));
            plot(xx, yy, 'r')
            plot(mu(i), 0,'k*')
        end
        title(strcat('MAP EM iteration#', num2str(iter)),'fontsize',15)
        drawnow;
    
end


% figure,
% bar(PI)