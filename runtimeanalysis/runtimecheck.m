%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Records the runtime of gradient descent with Armijo Rule(newdivv_estimator) and the hybrid  %%%
%%% method(gradient descent and newton (ddiv_estimator3)) with non robust                       %%%
%%% warmstarts                                                                                  %%%
%%%             (C) Andrews T. Anum and Michael Pokojovy (202                                   %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear

y = [30 50 100 500 1000 5000 10000 20000 50000];
alpha = 0.5;
simsize = 50000;
nbdp = ceil(0.2.*y);

for j=1:length(y)
    avgt_gd = 0;
    avgt2_gd = 0;
                                                          
    for i=1:simsize
        x = randn(y(j), 1);
	    x(randperm(y(j),nbdp(j))) = 300;
        mu0 = mean(x);
        sigma0 = std(x);
        tic
        Hn = newddiv_estimator(mu0, sigma0, x, alpha);
        t = toc;
        avgt_gd = avgt_gd + t/simsize;
        avgt2_gd = avgt2_gd + (t^2)/(simsize-1);
    end
                                                          
    var_gd = max(0, avgt2_gd - simsize/(simsize -1)*avgt_gd^2);
    h(j) = avgt_gd;
    s(j) = sqrt(var_gd/simsize);
end

for l=1:length(y)
    avgt = 0;
    avgt2 = 0;
  
    for m=1:simsize
        x = randn(y(l), 1);
        x(randperm(y(l),nbdp(l))) = 300;
	    mu0 = mean(x);
        sigma0 = std(x);
        tic
        Hn = ddiv_estimator3(mu0, sigma0, x, alpha);
        t = toc;
        avgt = avgt + t/simsize;
        avgt2 = avgt2 + (t^2)/(simsize-1);
    end
                                                          
    var_gdm = max(0, avgt2 - simsize/(simsize -1)*avgt^2);
    h1(l) = avgt;
    s1(l) = sqrt(var_gdm/simsize);
end

A = [y; h; s; h1; s1];
fileID = fopen('runtimetimeContaminated.txt','w');
fprintf(fileID,'%6s %12s %12s\n','n','GD', 'NM');
fprintf(fileID,'%6.0f %6.5f(%6.5f) %6.5f(%6.5f)\n',A);
fclose(fileID);
gd_mean_rt = h;
gd_ste_rt = s;
nm_mean_rt = h1;
nm_ste_rt = s1;

save('runtimetableContaminated.mat','gd_mean_rt','gd_ste_rt','nm_mean_rt', 'nm_ste_rt', 'y');