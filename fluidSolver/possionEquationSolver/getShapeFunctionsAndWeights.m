function [shapeFunction, GW, cord]= getShapeFunctionsAndWeights(numGP)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This function returns the guassian shape function matrix for a given number of
% gaussian points.
% Input: numGP - Number of Guassian Points
% Output: shapeFunction - shape function values at the guasstion points.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if (numGP == 0)
              shapeFunction(:,:,1) = [0] ;
              GW = [0];
              cord(:,:,1) = [0];
elseif(numGP == 1)
              shapeFunction(:,:,1) = [1/3, 0, 1/3, 0, 1/3, 0;
                                      0, 1/3, 0, 1/3, 0, 1/3];
              GW = [1];
              cord(:,:,1) = [1/3, 1/3, 1/3];
elseif(numGP == 3)
              shapeFunction(:,:,1) = [1/2, 0, 1/2, 0, 0, 0;
                                      0, 1/2, 0, 1/2, 0, 0];
                                  
              shapeFunction(:,:,2) = [0, 0, 1/2, 0, 1/2, 0;
                                      0, 0, 0, 1/2, 0, 1/2];
                                  
              shapeFunction(:,:,3) = [1/2, 0, 0, 0, 1/2, 0;
                                      0, 1/2, 0, 0, 0, 1/2];
              GW = [1/3, 1/3, 1/3];
              cord(:,:,1) = [1/2, 1/2, 0];
              cord(:,:,2) = [0, 1/2, 1/2];
              cord(:,:,3) = [1/2, 0, 1/2];
              
elseif(numGP == 4) 
              shapeFunction(:,:,1) = [1/3, 0, 1/3, 0, 1/3, 0;
                                      0, 1/3, 0, 1/3, 0, 1/3];
                                  
              shapeFunction(:,:,2) = [0.6, 0, 0.2, 0, 0.2, 0;
                                      0, 0.6, 0, 0.2, 0, 0.2];
                                  
              shapeFunction(:,:,3) = [0.2, 0, 0.6, 0, 0.2, 0;
                                      0, 0.2, 0, 0.6, 0, 0.2];
                                  
              shapeFunction(:,:,4) = [0.2, 0, 0.2, 0, 0.6, 0;
                                      0, 0.2, 0, 0.2, 0, 0.6];
                                  
              GW = [-27/48, 25/48, 25/48, 25/48];
              cord(:,:,1) = [1/3, 1/3, 1/3];
              cord(:,:,2) = [0.6, 0.2, 0.2];
              cord(:,:,3) = [0.2, 0.6, 0.2];
              cord(:,:,4) = [0.2, 0.2, 0.6];
              
elseif(numGP == 7)
              shapeFunction(:,:,1) = [0, 0, 0, 0, 1, 0;
                                      0, 0, 0, 0, 0, 1];
                                  
              shapeFunction(:,:,2) = [1/2, 0, 0, 0, 1/2, 0;
                                      0, 1/2, 0, 0, 0, 1/2];
                                  
              shapeFunction(:,:,3) = [1, 0, 0, 0, 0, 0;
                                      0, 1, 0, 0, 0, 0];
                                  
              shapeFunction(:,:,4) = [1/2, 0, 1/2, 0, 0, 0;
                                      0, 1/2, 0, 1/2, 0, 0];
                                  
              shapeFunction(:,:,5) = [0, 0, 1, 0, 0, 0;
                                      0, 0, 0, 1, 0, 0];
                                  
              shapeFunction(:,:,6) = [0, 0, 1/2, 0, 1/2, 0;
                                      0, 0, 0, 1/2, 0, 1/2];
                                  
              shapeFunction(:,:,7) = [1/3, 0, 1/3, 0, 1/3, 0;
                                      0, 1/3, 0, 1/3, 0, 1/3];
                                  
              GW = [1/40, 1/15, 1/40, 1/15, 1/40, 1/15, 9/40];
              cord(:,:,1) = [0, 0, 1];
              cord(:,:,2) = [0.5, 0, 0.5];
              cord(:,:,3) = [1, 0, 0];
              cord(:,:,4) = [0.5, 0.5, 0];
              cord(:,:,5) = [0, 1, 0];
              cord(:,:,6) = [0, 0.5, 0.5];
              cord(:,:,7) = [1/3, 1/3, 1/3];
end

end