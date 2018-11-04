%% computeStatisticsFromTensor

function [Dvector,Dmatrix,EV,EW,FA,notposdef] = computeStatisticsFromTensor(D)

% notposdef = 0 (D is positive definite), notposdef = 1 (D is not positive definite)
notposdef = 0;

% Set inf and nan entries to zero
contains_inf = sum(isinf(D(:)))>0;
contains_nan = sum(isnan(D(:)))>0;
if contains_inf || contains_nan
    notposdef = 1;
    D(isinf(D)) = 0;
    D(isnan(D)) = 0;
end

% Turn the vector into a matrix, or the matrix into a vector
if numel(D)==6
    Dvector = D(:);
    Dmatrix = [D(1) D(4) D(5);D(4) D(2) D(6);D(5) D(6) D(3)];
elseif numel(D)==9
    Dvector = [D(1,1);D(2,2);D(3,3);D(1,2);D(1,3);D(2,3)];
    Dmatrix = D;
else
    error('Input D does not contain 6 or 9 elements');
end

% Compute and sort the eigenvalues EW and eigenvectors EV
[EV,EW] = eig(Dmatrix);
[EW,ind] = sort(diag(EW),'descend');
EV = EV(:,ind);

% Set negative eigenvalues to zero
contains_negativeEW = sum(EW<0)>0;
if contains_negativeEW
    notposdef = 1;
    EW(EW<0) = 0;
    
    % Recompute Dvector and Dmatrix
    Dmatrix = EV*diag(EW)*EV';
    Dvector = [Dmatrix(1,1);Dmatrix(2,2);Dmatrix(3,3);Dmatrix(1,2);Dmatrix(1,3);Dmatrix(2,3)];
end

% Compute the fractional anisotropy
if EW(1)^2+EW(2)^2+EW(3)^2>0 
    FA = sqrt((EW(1)-EW(2))^2+(EW(2)-EW(3))^2+(EW(1)-EW(3))^2)/(sqrt(2)*sqrt(EW(1)^2+EW(2)^2+EW(3)^2));
else
    FA = 0;
end

end