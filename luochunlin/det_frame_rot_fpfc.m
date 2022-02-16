function [fPlus,fCross]=det_frame_rot_fpfc(polAngleTheta,polAnglePhi,Psi)

addpath '../GWSIG'

%Number of locations requested
nLocs = length(polAngleTheta);
if length(polAnglePhi) ~= nLocs
    error('Number of theta and phi values must be the same');
end

%Obtain the components of the unit vector pointing to the source location
sinTheta = sin(polAngleTheta(:));
vec2Src = [sinTheta.*cos(polAnglePhi(:)),...
           sinTheta.*sin(polAnglePhi(:)),...
           cos(polAngleTheta(:))];
       
%Get the wave frame vector components (for multiple sky locations if needed)
xVec = vcrossprod(repmat([0,0,1],nLocs,1),vec2Src);
yVec = vcrossprod(xVec,vec2Src);
%Normalize wave frame vectors
for lpl = 1:nLocs
    xVec(lpl,:) = xVec(lpl,:)/norm(xVec(lpl,:));
    yVec(lpl,:) = yVec(lpl,:)/norm(yVec(lpl,:));
end
%rotation: multiply a rotation matrix
for lpl = 1:nLocs
    rot_mat = [cos(Psi(lpl)),sin(Psi(lpl));
                -sin(Psi(lpl)),cos(Psi(lpl))];
    tmp = rot_mat * [xVec(lpl,:);yVec(lpl,:)];
    xVec(lpl,:) = tmp(1,:);
    yVec(lpl,:) = tmp(2,:);
end


%Detector tensor of a perpendicular arm interferometer 
detTensor = [1,0,0]'*[1,0,0]-[0,1,0]'*[0,1,0];
fPlus = zeros(1,nLocs);
fCross = zeros(1,nLocs);
%For each location ...
for lpl = 1:nLocs
    %ePlus contraction with detector tensor
    waveTensor=xVec(lpl,:)'*xVec(lpl,:)-yVec(lpl,:)'*yVec(lpl,:);
    fPlus(lpl) = sum(waveTensor(:).*detTensor(:));
    %eCross contraction with detector tensor
    waveTensor = xVec(lpl,:)'*yVec(lpl,:)+yVec(lpl,:)'*xVec(lpl,:);
    fCross(lpl) = sum(waveTensor(:).*detTensor(:));
end


end
