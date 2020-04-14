function T=stepT(Tseg,noseg,noNL)
%""""
    % Creating step temperature profile w/ temperature as inputs
        % Note: The function is defined to process even number of segments. 
        % Should be extended to uneven if necessary.
    % Args:
        % Tseg: constant temperature in each segment
        % noseg: no of segments
        % noNL: resolution in length steps
    % Return:    
        % T: Step temperature profile
%""""

for n=1:1:noseg-1
    startTseg=(Tseg(n)+Tseg(n+1))/2;
    diffTseg=(Tseg(n+1)-Tseg(n))/2;
    Two=startTseg+diffTseg*tanh(0:1:(2*noNL-1)-noNL);
        % Two temperature segments
        
    if n==1 % beginning
        T=Two';
    else
        T=vertcat(T,Two'); % combine two segments
    end
end