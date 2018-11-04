function test_suite = test_fit_MRI
% This function provides the tests for fit_MRI
% Uses the xUnit test suite: 
% http://www.mathworks.com/matlabcentral/fileexchange/22846-matlab-xunit-test-framework
%
% Created by Dirk Poot, Erasmus MC, 2-2-2012

initTestSuite;
% test_project_2D
% test_simple_2Dn
% test_simple_1D_sigmest
% test_project_2D

function test_simple_scalar1D
%% Simple test of Rice data, 1 measurement per point, 1D 'image'
% Effectively, fit_MRI 'corrects' for bias due to rice noise.
% we have to provide the noise level, since it can't be estimated from 1 measurment.
npnttrace = 1;
fun = @(par) verysimpletestfitfun(par, npnttrace);
data = linspace(1,10,10);
sigma = 1;
init = 1;  % scalar initial value, let fit_MRI determine correct size; which it obtains from data.
tht = fit_MRI(fun, data, init,'noiseLevel',sigma);
% Compare with true optimum; fit_MRI will not reach this exactly:
tht_true = [0.000000000000000   1.662924049508514   2.816438372461607   3.868512411768108   4.896804442808278 5.914854010161516   6.927443261715339   7.936749866960533   8.943920278071543   9.949619262232044];
% Obtain true optimimum by newton iterations (and manually check that they converge):
% [L,dL, hL] = logricepdf(data, tht, sigma,[false true false]);tht = tht - (dL./hL)';
assertElementsAlmostEqual(tht, tht_true ,'absolute',1e-3);

function test_simple_scalar2D
%% Simple test of Rice data, 1 measurement per point, 2D 'image'
% Effectively, fit_MRI 'corrects' for bias due to rice noise.
% we have to provide the noise level, since it can't be estimated from 1 measurment.
npnttrace = 1;
fun = @(par) verysimpletestfitfun(par, npnttrace);
data = reshape(linspace(1,10,25),[1 5 5]); % 2D data, scalar per voxel.
sigma = 1;
init = 1;  % scalar initial value, let fit_MRI determine correct size; which it obtains from data.
tht = fit_MRI(fun, data, init,'noiseLevel',sigma);
% Compare with true optimum; fit_MRI will not reach this exactly:
tht_true = reshape([ -0.000000000000000   0.000000000000000   1.286147903556077   1.824591934380471   2.267174336388405   2.681475884689648   3.083305878131754   3.478134357979410   3.868512411768106   4.255844271640894   4.640982632120651   5.024481110503180   5.406717358025005   5.787959057572419   6.168401873162328   6.548192494958466   6.927443261715339   7.306241781379571  7.684657456819911   8.062746029193750   8.440552813369003   8.818115047923291   9.195463632103765   9.572624429801195   9.949619262232044],[1 5 5]);
% Obtain true optimimum by newton iterations (and manually check that they converge):
% [L,dL, hL] = logricepdf(data(:), tht(:),sigma,[false true false]);tht(:) = tht(:) - (dL./hL);
assertElementsAlmostEqual(tht, tht_true ,'absolute',1e-3);

function test_simple_2D
%% Simple test of Rice data, 2 measurement per point, 2D 'image'
% Effectively, fit_MRI 'corrects' for bias due to rice noise.
% we provide the noise level.
npnttrace = 2;
fun = @(par) verysimpletestfitfun(par, npnttrace);
[x,y]=ndgrid([1:5],[1:4]);
data = [permute(x,[3 1 2]);permute(y,[3 1 2])]; % 2D data, 2 data elements per voxel.
sigma = 1;
init = 1;  % scalar initial value, let fit_MRI determine correct size; which it obtains from data.
tht = fit_MRI(fun, data, init,'noiseLevel',sigma);
% Compare with true optimum; fit_MRI will not reach this exactly:
tht_true = reshape([ 0.000000000000000   0.815502151236574   1.654920615540983   2.253179304388592   2.806504347228211   0.815502151236574   1.662924049508514   2.266079522975673  2.815226099638210   3.346017219424885   1.654920615540983   2.266079522975673   2.816438372461607   3.347017793270996   3.868330681965017   2.253179304388592   2.815226099638210   3.347017793270996   3.868512411768106   4.384413517876916],[1 5 4]);
% Obtain true optimimum by newton iterations (and manually check that they converge):
% [L,dL, hL] = logricepdf(data, repmat(tht,[2 1 1]),sigma,[false true false]);dL = dL(1:2:end)+dL(2:2:end);hL = hL(1:2:end)+hL(2:2:end);tht(:) = tht(:) - (dL./hL);
assertElementsAlmostEqual(tht, tht_true ,'absolute',1e-3);

function test_simple_2Dn
%% Simple test of Rice data, 2 measurement per point, 2D 'image'
% Effectively, fit_MRI 'corrects' for bias due to rice noise.
% we estimate the noise level.
npnttrace = 3;
fun = @(par) verysimpletestfitfun(par, npnttrace);
[x,y]=ndgrid([1:5],[1:4]);
data = [permute(1+x,[3 1 2]);permute(x+y/2,[3 1 2]);permute((x+y/3+1),[3 1 2])]; % 2D data, 3 data elements per voxel.
init = 4;  % scalar initial value, let fit_MRI determine correct size; which it obtains from data.
% fit tht with automatically estimated global noise level:
[tht0 , resid, opts]= fit_MRI(fun, data, init); 
[tht1 , resid, opts]= fit_MRI(fun, data, tht0); 
[tht2 , resid, opts]= fit_MRI(fun, data, tht1); 
[tht , resid, opts]= fit_MRI(fun, data, tht2); 
if any(tht0(:)-tht(:)>2e-2)
    warning('fit_MRI should be improved to iterate sufficiently to reach a true optimum');
end;
% Compare with true optimum; fit_MRI will not reach this exactly:
tht_true = reshape([ 1.897123660579510   2.913871600406928   3.921783877297499   4.926424434276874   5.929481351985710   2.181203922897392   3.194359727884693   4.201076251979444   5.205170534570016   6.207931889940691   2.463759820558987   3.474400358405967   4.480177306125979   5.483817345224460   6.486324314045040   2.745298319804716   3.754097344080219   4.759121440593551   5.762379480791726   6.764665874094002],[1 5 4]);
sigma_true = 0.420977548173955;
% Obtain true optimimum by newton iterations (and manually check that they converge):
if 0
    sigma = .5;
    hLsigm = 2*-2.674653260149436e+002; % appoximate hessian w.r.t. sigma
%%  
    [L,dLsigm] = logricepdf(data, repmat(tht,[3 1 1]),sigma,[false false true]);
    sigma = sigma - sum(dLsigm)/hLsigm
    sum(dLsigm)
%%
    [L,dL, hL] = logricepdf(data, repmat(tht,[3 1 1]),sigma,[false true false]);
    dL = dL(1:3:end)+dL(2:3:end)+dL(3:3:end);hL = hL(1:3:end)+hL(2:3:end)+hL(3:3:end);
    tht(:) = tht(:) - (dL./hL);
    plot(dL)
end;
assertElementsAlmostEqual(tht, tht_true ,'absolute',2e-2);
assertElementsAlmostEqual(opts.noiseLevel, sigma_true,'absolute',.045);

function test_simple_1D_sigmest
%% Simple test of Rice data, 3 measurements per point, 1D 'image'
%  1 intensity parameter per point & estimate noise level per point
npnttrace = 3;
fun = @(par) verysimpletestfitfun(par, npnttrace);
data = [1:10;2:11;linspace(3,9,10)]; % 1D data, 3 data elements per voxel.
init = [1:10;ones(1,10)];  % full initial value. (this is the first test that uses that)

[tht1 , resid, opts]= fit_MRI(fun, data, init,'numPDFoptpar',1); 
[tht , resid, opts]= fit_MRI(fun, data, tht1,'numPDFoptpar',1); 
% Note: sign of sigma in tht is ambiguous (for rice distribution); set to positive:
tht(2,:) = abs(tht(2,:));

%% Compare validated known-good optimum:
tht_true = [   1.732492716658367   2.800770202055111   3.734013263376019   4.642543994760954   5.539889300692178   6.430963869559172   7.318102609944162   8.202585238150029   9.085180943393423   9.966383376083606
               0.912451691381977   0.696599271495924   0.569972659227605   0.472644293686604   0.416331362024573   0.416177493461763   0.471897335413947   0.567239516292729   0.685918620231256   0.817884174246390];
% Obtain true optimimum by newton iterations (and manually check that they converge):
assertElementsAlmostEqual(tht, tht_true ,'absolute',1e-8);
if any(tht(:)-tht1(:)>1e-7)
    warning('fit_MRI did not converge the first time, improve iteration/convergence criteria');
end;

function test_simple_1D_real
%% Simple test of Rice data, 3 measurements per point, 1D 'image'
%  1 intensity parameter per point & estimate noise level per point
npnttrace = 3;
fun = @(par) verysimpletestfitfun(par, npnttrace);
data = [1:10;2:11;linspace(3,9,10)]; % 1D data, 3 data elements per voxel.
init = [1:10;];  % full initial value. (this is the first test that uses that)

warning('fit_MRI should update noise final returned noise level and should not require an extra call for that.');
[tht , resid, opts]= fit_MRI(fun, data, init,'imageType','real'); 
[thtb , residb, optsb]= fit_MRI(fun, data, tht,'imageType','real'); 

%% Compare validated known-good optimum:
tht_true = mean(data,1);
% Obtain true optimimum by newton iterations (and manually check that they converge):
assertElementsAlmostEqual(tht, tht_true ,'absolute',1e-8);
% std(reshape(bsxfun(@minus, data, mean(data,1)),[],1)) * sqrt(numel(data)/(numel(data)-numel(tht)))
sigm_orig = 0.693888666488711;
assertElementsAlmostEqual(optsb.noiseLevel, sigm_orig,'absolute',1e-4);

function test_laplacian_1D
%% test of Rice data, 3 measurements per point, 1D 'image'
% estimate noise level per point, laplacian regularization
npnttrace = 3;
fun = @(par) verysimpletestfitfun(par, npnttrace);
data = [1:10;2:11;linspace(3,9,10)]; % 1D data, 3 data elements per voxel.
init = [1:10;ones(1,10)];  % full initial value. (this is the first test that uses that)
data = data(:,[1:2:end 2:2:end]); % permute data;

warning('test_laplacian_1D: Need to redo fitting for accurate results; should improve fit_MRI to do this in 1 iteration');
tht0=fit_MRI(fun, data, init, 'numPDFoptpar',1, 'spatialRegularizer','laplacian','spatialRegularizerWeights',.1*[1;1]);
tht1=fit_MRI(fun, data, tht0, 'numPDFoptpar',1, 'spatialRegularizer','laplacian','spatialRegularizerWeights',.1*[1;1]);
tht=fit_MRI(fun, data, tht1, 'numPDFoptpar',1, 'spatialRegularizer','laplacian','spatialRegularizerWeights',.1*[1;1],'tolFun',1e-9);

% value from reference evaluation:
% Old laplacian implementation:
% tht_true = [1.744060409456432   3.730042221401462   5.544193673905747   7.359851374574939   8.610212470068737   3.171067010147854   4.595872939648574   6.427975184754629   8.202964054346023   9.966726256338362
%             0.904672480059735   0.571082043691280   0.416041363400607   0.475592170862546   0.837506019440313   0.746345904281026   0.478780289155977   0.416635695740093   0.567148826163887   0.816765191844148];
% new laplacian implementation (/definition)
tht_true = [1.818530733129475   3.711108803956268   5.543347680496735   7.359965464572171   8.610281278386939   3.171046451832239   4.595762505503915   6.428552998163638   8.219681110134424   9.930181672804078
            0.874915860424290   0.574279313103472   0.416186030193200   0.475588111279389   0.837466716378263   0.746337500602165   0.478792652601683   0.416626141809292   0.568134224660504   0.816345423944215];
assertElementsAlmostEqual(tht, tht_true ,'absolute',1e-4);

function test_TV_1D
%% test of Rice data, 3 measurements per point, 1D 'image'
% estimate noise level per point, laplacian regularization
npnttrace = 3;
fun = @(par) verysimpletestfitfun(par, npnttrace);
data = [1:10;2:11;linspace(3,9,10)]; % 1D data, 3 data elements per voxel.
init = [1:10;ones(1,10)];  % full initial value. (this is the first test that uses that)
data = data(:,[1:2:end 2:2:end]); % permute data;
tht=fit_MRI(fun, data, init, 'numPDFoptpar',1, 'spatialRegularizer','TV','spatialRegularizerWeights',.1*[1;1]);
tht(2,:) = abs(tht(2,:));
% currently need to redo to accurately compute:
warning('test_TV_1D: Need to redo fitting for accurate results; should improve fit_MRI to do this in 1 iteration');
tht0=fit_MRI(fun, data, tht, 'numPDFoptpar',1, 'spatialRegularizer','TV','spatialRegularizerWeights',.1*[1;1]);
tht=fit_MRI(fun, data, tht0, 'numPDFoptpar',1, 'spatialRegularizer','TV','spatialRegularizerWeights',.1*[1;1]);

% value from reference evaluation:
tht_true = [1.923951213110853   3.733845304271746   5.535106658285712   7.322558532849872   8.916459357182035   2.966634477812451   4.637126144061640   6.435360225235911   8.200486637038875   9.848558145952412
            0.847266445515007   0.571933718825787   0.418285744481692   0.474079539977682   0.707722150435091   0.690911339668232   0.474960291337821   0.417804074073988   0.569141748225943   0.822190185585068];
assertElementsAlmostEqual(tht, tht_true ,'absolute',1.5e-4);

function test_project_2D
%% test of normal data, 3 measurements per point, 3D 'image'
% 
npnttrace = 3;
fun = @(par) verysimpletestfitfun(par, npnttrace);
[x,y]=ndgrid([1:5],[1:4]);
trueimg = [permute(1+x,[3 1 2]);permute(x+y/2,[3 1 2]);permute((x+y/3+1),[3 1 2])]; % 2D data, 3 data elements per voxel.
init = ones(1,5,4);  % full initial value. 
project = build_project(npnttrace, [5 4],[0 2]);
data = cell(size(project,1),1);
for k=1:numel(data)
    data{k} = project{k,1}(squeeze(trueimg(k,:,:)),[]);
end;

warning('test_project_2D: Noise level initialisation should not be needed in this call.') 
% using fminunc:
tht=fit_MRI(fun, data, init, 'imageType','real', 'spatialRegularizer','laplacian','spatialRegularizerWeights',.1*[1],'project',project,'noiseLevel',.4, 'projectOptimization.skip',true,'optimizeBlocks',false,'projectParameters',cell(npnttrace,1),'tolFun',1e-7);
% using my own optimizer:
tht2=fit_MRI(fun, data, init, 'imageType','real', 'spatialRegularizer','laplacian','spatialRegularizerWeights',.1*[1],'project',project,'noiseLevel',.4, 'projectOptimization.skip',true,'optimizeBlocks',false,'projectParameters',cell(npnttrace,1),'tolFun',1e-9,'optimizer',3,'maxIter',23);
%% value from reference evaluation:
% Reference with old Laplace implementation:
% tht_true = reshape([1.828766025499089   1.978419684210352   2.161022958153525   2.106777315116103
%                     2.529742767313631   2.948994635257290   3.397827134165548   3.443304723444858
%                     3.511899281657698   3.962283229882416   4.362543307047687   4.437943547745730
%                     4.590055384360963   4.985716975726203   5.444205330152977   5.735105847350146
%                     5.457692349635800   6.172767281570230   6.694968541483679   7.195921946405083],[1 5 4]);
% Reference with new Laplace implementation, modified noise level parameter and imageType set to 'real'                
tht_true = reshape([1.966031213122669   2.412911364141812   3.072141944541731   3.300134223224392
   2.663259837214939   3.015849464664815   3.895230090577886   3.578433029951802
   3.785201396066925   4.002055771458704   4.669062649204865   4.399248205914655
   4.967644371193494   5.087281607189643   5.736446597228306   5.972261861279280
   5.535656109726159   6.220586449396960   6.780539585482639   7.242337288905348],[1 5 4]);
assertElementsAlmostEqual(tht, tht_true ,'absolute',1.5e-3);
assertElementsAlmostEqual(tht2, tht_true ,'absolute',1.5e-3);


function project = build_project( nimg, spatialsz, opt )
% nimg = number of images
% spatialsz = spatial size of each image
% opt(1) = fill in threshold; as threshold for a normal distributed value.
% opt(2) : 0 => multiplication with 1
%          1 => multiplication with random magnitude
%          2 => multiplication with random magnitude + random offset

project = cell(nimg, 5);
oldstate = randn('state');
randn('state',[10 1000]'); % just some fixed state.
mask = cell(nimg,1);
offsets = cell(nimg,1);
for k=1:numel(mask);
    mask{k} = randn(spatialsz)>opt(1); 
    offsets{k} =0;
    if opt(2)>=1
        mask{k} = mask{k} .* (1+randn(spatialsz));
    end;
    if opt(2)>=2
        offsets{k} = randn(nnz(mask{k}),1);
    end;
end;
for k=1:nimg
    project{k,1} = @(x,par) projectfunction(x, mask{k}, offsets{k});
    project{k,2} = @(x, g ) projectfunction(x, g, 0);
    project{k,3} = @(x,par) projectfunctionAdj(x, mask{k});
    % don't specify columns 4 and 5 since we don't have project adjust parameters
end;
randn('state',oldstate);

function [y , gradifo] = projectfunction( x, mask, offset )
tmp = x.*mask;
y = tmp(logical(mask))+offset;
gradifo = mask;

function [x , gradifo] = projectfunctionAdj( y, mask )
tmp = y.*mask(logical(mask));
x = zeros(size(mask));
x(logical(mask)) = tmp;
gradifo = mask;


function test_extraInitialValues
%%
opts =struct;
opts.initialValueSpecifierVect=[1 0 0 0 1000];
opts.initialValueRange = [0.3 4;10 50;20 100];
TI = (100:100:2000)';
fun = @(tht) predict_IRT1(tht, TI);
thtGT = repmat([.7 .98 1.4 .7;10 23 43 15;15 30 79 30],[1 1]);
data = fun( thtGT );
rngstate = rng; rng(5488,'twister');
datan = AddRiceNoise( repmat(data,[1 1 4]),1);
rng(rngstate)
opts.noiseLevel = .5;
% tic
thtest = fit_MRI( fun, datan, zeros(3,size(datan,2),size(datan,3)), opts);
opts.blockSize = [1 1]; 
thtest2 = fit_MRI( fun, datan, zeros(3,size(datan,2),size(datan,3)), opts);
% toc
tht_verified = cat(3,[ 0.328889877577874   1.019364016991038   1.411935079822267   0.548034793626360
  19.354500944948445  22.363042928357796  43.789654101365706  18.190465398624625
  24.151669214611914  29.182904343554849  83.201128050995734  31.469201086779467],...
[  0.311698949695022   0.965620001906880   1.338703270009090   0.475771451053402
  20.393490427051731  23.119127191840899  43.383567759417268  22.381498414793803
  25.345014326352118  30.158850100370092  77.487403820503076  36.372711342444134],...
[  0.904911524467240   0.914496395266664   1.375776803032240   0.718467810872458
   8.699382565285774  24.055201713498132  43.614277513224600  13.565729719750339
  14.869670228150307  29.880164195769677  79.923238567767186  28.209700202552590],...
[  0.901855080815093   0.993121235361927   1.361718639330074   0.217402641226252
   8.158910893361092  22.844015146051937  43.108945133325847  49.570293708298294
  13.696072873253147  29.964909679737168  78.707807289031905  63.037738843503874]);
assertElementsAlmostEqual(thtest, tht_verified ,'absolute',.01052);
assertElementsAlmostEqual(thtest2(:,:,1:3), tht_verified(:,:,1:3) ,'absolute',.0113);

doplotfitfigure = false;
if doplotfitfigure
    %% compare to fit with just a single initial value:
    tic
    thtest_noexinit = fit_MRI( fun, datan, [1.5;20;40],'noiseLevel',.5);
    % thtest_noexinit = fit_MRI( fun, datan, thtest_noexinit);
    toc
    %%
    rp = 1;
    TIpred = (100:10:2000)';
    fun_fit = @(tht) predict_IRT1(tht,TIpred);
    GTpred = fun_fit(thtGT);
    fitpred = fun_fit(thtest(:,:,rp));
    fitpred2 = fun_fit(thtest_noexinit(:,:,rp));
    lh = plot(TI,datan(:,:,1),'*',TIpred, GTpred(:,:), '-',TIpred, fitpred, '--',TIpred, fitpred2, ':');
    for k=1:size(datan,2)
        set(lh(k+[1 2 3]*size(datan,2)),'color',get(lh(k),'color'))
    end;
end;


function test_LinearParameterInitialValues
%% test option 'linearParameters'
opts =struct;
opts.initialValueSpecifierVect=[1 0 0 0 50];
opts.initialValueRange = [0.3 4;1 50;0 5];
opts.linearParameters = [1 3];
TE = ([0:10:100 300])'/1000;
fun = @(tht) predict_MRI_T2(tht, TE, false);
thtGT = repmat([20  21 22 53;2 5 8.5 20;0 1 2 3],[1 1]);
data = fun( thtGT );
rngstate = rng; rng(5488,'twister');
datan = AddRiceNoise( repmat(data,[1 1 4]),1);
rng(rngstate)
opts.noiseLevel = .5;
% tic
thtest = fit_MRI( fun, datan, zeros(3,size(datan,2),size(datan,3)), opts);
%%
opts.blockSize = [1 1]; 
thtest2 = fit_MRI( fun, datan, zeros(3,size(datan,2),size(datan,3)), opts);
% toc
tht_verified = cat(3, [15.206584200058549  27.505286852348405  21.998056347448802  53.691518673639997
   2.883597675582596   2.267993143516489   8.243735375742009  18.860897829741464
   4.682448542688193  -7.099769421785617   1.781405433979219   2.412236372558072],[
   8.807368882426585  24.880807414226823  19.716482562504726  51.331232391926690
   4.459596597319917   3.095333164410389   9.345015050453430  20.963395456879240
  11.103006689343935  -4.189768601595298   3.798711691944131   3.768123284073045],[
   8.567738628765616  18.854425560725922  22.032002341590214  53.586845055285679
   6.257306250936035   5.439731861506655   8.217777733591641  20.830668302556756
  11.710851218919920   2.752310888721676   1.151350915432404   3.119582385667560],[
  19.638105760664693  20.076447046156240  22.851750288177264  50.283679383742715
   2.278889398883035   6.640851179070810   7.704282910336846  20.119914969624705
   0.607966118116301   3.018108667920771   0.861550205946661   4.102146375812036]);
assertElementsAlmostEqual(thtest, tht_verified ,'absolute',.0012);
assertElementsAlmostEqual(thtest2, tht_verified ,'absolute',.003);

%%
function test_ConstrainedOptimization
%%
data= 10* randn(10,20,20);
X=[ones(10,1) (0:9)'];
fun = @(x) verysimpletestfitfun(x,X,2);
clear fop
% fop.constraints.lb=[0;0]
% fop.constraints.ub=[1;1]
% fop.constraints.Aineq = [eye(2);-eye(2)]
% fop.constraints.bineq = [1;1;0;0]
fop.imageType = 'real';
tht0 = fit_MRI( fun, data, [0;1],fop);

fop1 = fop;
fop1.constraints.lb=[-.1;-.2];
fop1.constraints.ub=[1.1;1.4];
tht1 = fit_MRI( fun, data, [0;1],fop1);

fop2 = fop;
fop2.constraints.Aineq = [eye(2);-eye(2)];
fop2.constraints.bineq = [1.1;1.4;.1;.2];
tht2 = fit_MRI( fun, data, [0;1],fop2);
% assertElementsAlmostEqual(tht1, tht2,'absolute',.00001);
if any(any(bsxfun(@gt, fop2.constraints.Aineq*tht2(:,:),fop2.constraints.bineq)))
    error('constraints not satisfied');
end;
if any( median(abs(tht1(:,:)'-tht2(:,:)')) > .0015 )
    error('Difference between implicit and explicity specification of upper and lower bounds.');
end;