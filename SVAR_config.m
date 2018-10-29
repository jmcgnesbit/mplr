classdef SVAR_config
    
   properties 
       % Data
       data
       names
       constant     % Include constant
       nlag        % Number of lags
       nhorizon    % Number of horizons
       
       restrictionMatrix
       
       epsilon
       delta
       
       drawMatrix
       
       Y
       X
       
       
       nobs
       nvar
       T
       
       A
       Sigma
       uhat
       
       cholsig
   end
   % Constructor
   methods
       function obj = SVAR_config(data, names, constant, nlag, nhorizon, restrictionMatrix, epsilon, delta, drawMatrix)
           % Data
           obj.data = data;
           obj.names = names;
           
           obj.constant = constant;
           obj.nlag = nlag;
           obj.nhorizon = nhorizon;
           
           obj.restrictionMatrix = restrictionMatrix;
           
           obj.epsilon = epsilon;
           obj.delta = delta; 
          
           obj.drawMatrix = drawMatrix;
           
           [obj.Y, obj.X] = VARmakexy(data, nlag, constant);
           obj.T = size(data,1);
           
           [obj.nobs, obj.nvar] = size(obj.Y);
           
           [obj.A, obj.Sigma, obj.uhat] = VARestimate(obj.Y, obj.X);
           
           obj.cholsig = chol(obj.Sigma, 'lower');
       end
       
    end
end