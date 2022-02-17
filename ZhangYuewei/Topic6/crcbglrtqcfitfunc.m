function [fitVal,varargout] = crcbglrtqcfitfunc(xVec,params)
%Fitness function for quadratic chirp regression
%F = CRCBQCFITFUNC(X,P)
%Compute the fitness function (sum of squared residuals function after
%maximimzation over the amplitude parameter) for data containing the
%quadratic chirp signal. X.  The fitness values are returned in F. X is
%standardized, that is 0<=X(i,j)<=1. 
% The fields of P are:
% 'dataY': The data vector (a uniformly sampled time series).
% 'dataX': The time stamps of the data samples.
% 'dataXSq': dataX.^2
% 'dataXCb': dataX.^3
% 'sampFreq': sampling frequency
% 'psdPosFreq': power spectral density
% 'rmin', 'rmax': The minimum and maximum values of the three parameters
%                 a1, a2, a3 in the candidate signal:
%                 a1*dataX+a2*dataXSq+a3*dataXCb
% rmin and rmax are used
%to convert X(i,j) internally before computing the fitness:
%X(:,j) -> X(:,j)*(rmax(j)-rmin(j))+rmin(j).
%The fields P.dataY and P.dataX are used to transport the data and its
%time stamps. The fields P.dataXSq and P.dataXCb contain the timestamps
%squared and cubed respectively.
%
%[F,R] = CRCBQCFITFUNC(X,P)
%returns the quadratic chirp coefficients corresponding to the rows of X in R. 
%
%[F,R,S] = CRCBQCFITFUNC(X,P)
%Returns the quadratic chirp signals corresponding to the rows of X in S.

%Soumya D. Mohanty
%June, 2011
%April 2012: Modified to switch between standardized and real coordinates.

%Shihan Weerathunga
%April 2012: Modified to add the function rastrigin.

%Soumya D. Mohanty
%May 2018: Adapted from rastrigin function.

%Soumya D. Mohanty
%Adapted from QUADCHIRPFITFUNC
%==========================================================================

% xVec: popsize*D; params:1*1 structure containing rmin and rmax of D
% parameters
%rows: points
%columns: coordinates of a point
[nVecs,~]=size(xVec);

%storage for fitness values
fitVal = zeros(nVecs,1);

%Check for out of bound coordinates and flag them
validPts = crcbchkstdsrchrng(xVec);
%Set fitness for invalid points to infty
fitVal(~validPts)=inf;
xVec(validPts,:) = s2rv(xVec(validPts,:),params);

for lpc = 1:nVecs
    if validPts(lpc)
    % Only the body of this block should be replaced for different fitness
    % functions
        x = xVec(lpc,:);                        % its size shoule be 1*3 as [a1, a2, a3]
        fitVal(lpc) = glrtqcsig_Mod(x, params); % My Modification 
    end
end

%Return real coordinates if requested
if nargout > 1
    varargout{1}=xVec;
end
