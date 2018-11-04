%%
% fit_MRI_test_MCMC
if 0
    % Test with simple linear model
    X = rand(10,3);
    % tht = rand(3,4,3);
    tht = repmat(rand(3,1),[1 ,4 , 3] );
    thtinit = ones(size(tht));
    noiselevel = 1;

    predfun = @(par) verysimpletestfitfun( par, X ,2);
else
    % Test with T2 decay model, large enough A to avoid sign flips
    % (hopefully)
    TE = [10 20 40 60 80 100]'/1000;
    tht = repmat([10; 1000/50],[1 3 3]);
    thtinit = 10*ones(size(tht));
    noiselevel = .7;
    predfun = @(par) predict_MRI_T2( par, TE ,false);
end;

ynoisefree = predfun(tht(:,:)); 
if 0 
    ynoise = ynoisefree+noiselevel*randn(size(ynoisefree)) 
else
    ynoise = AddRiceNoise( ynoisefree, noiselevel );
end;
y= reshape( ynoise, [size(ynoisefree,1), size(tht,2) size(tht,3)] );
%%
% profile on;
clear fitopt;
fitopt.optimizer = 4; % MCMC
fitopt.noiseLevel = noiselevel;
fitopt.imageType = 'magnitude';
fitopt.blockSize = [1 1 ];
tstpar = cell(1,10000);
tstpar{1} = thtinit;
progressbar('start',numel(tstpar));
for k=2:numel(tstpar)
    tstpar{k} = fit_MRI( predfun, y, tstpar{k-1}, fitopt);
    if mod(k,100)==0
        progressbar(k);
    end;
end;
progressbar('ready');
tmp = cat(4,tstpar{2:end});
tmp1 =permute(tmp(:,1,1,:),[4 1 2 3]);
c=cov(tmp1)
tmp2 = tmp(:,:)';

rejectratio = (mean(all(diff(tmp,[],4)==0,1),4))
% profile viewer

%% 'correct' sign:
sgnadj = sum(X*tmp(:,:),1)<0;
tmp(:,sgnadj) = -tmp(:,sgnadj);
%% 
mn = min(tmp2);
mx = max(tmp2);
ndim = numel(mn);
if ndim==2
    nsamp = 200;
else
    nsamp = 50;
end;
x_in = cell(1,ndim);
for k=1:ndim
    x_in{k} = linspace(mn(k),mx(k),nsamp);
end;
x = cell(1,ndim);
[x{:}] = ndgrid(x_in{:});
szx = size(x{1});
x = permute( cat(ndim+1, x{:}),[ ndim+1 1:ndim]);
szp = size(tstpar{2});
numvox = prod(szp(2:end));
LL = zeros([szx numvox]);

kst = 1;
for voxid = 1: numvox
    y_tst = y(:,voxid);
    for k=1:prod(szx);
        LL(kst) = sum( logricepdf( y_tst, predfun(x(:,k)), noiselevel) );
        kst = kst+1;
    end;
end;
%%
pnts = max(1,min(nsamp,bsxfun(@times, bsxfun(@minus, tmp2, mn ) , (nsamp-1)./(mx-mn))+1));

LLsmp = zeros([szx numvox]);
imgsel = repmat({':'},1,ndim);
for k=1:numvox
    LLsmp(imgsel{:},k)  = genAddGaussBumps( szx, pnts(k:numvox:end,:), [] , .6 );
end;

%% display log likelihood true and estimated:
imagebrowse(cat(ndim+2,LL-max(LL(:)),log(LLsmp./max(LLsmp(:)))),[-10 0],'xrange',[mn 1 1;mx numvox 2] )
%% display log ratio:
% sumLLsmp = sum(LLsmp,4);
imagebrowse(LL-max(LL(:))-log(LLsmp./max(LLsmp(:))),[-5 5],'xrange',[mn 1;mx numvox] )