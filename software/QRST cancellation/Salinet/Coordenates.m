

% X_f = X_f (1:64, 2:33);
% Y_f = Y_f (1:64, 2:33);
% Z_f = Z_f (1:64, 2:33);
% 
% X_f = reshape(X_f,[1 2048]);
% Y_f = reshape(Y_f,[1 2048]);
% Z_f = reshape(Z_f,[1 2048]);


a = -20.06;
b = 27.05;
c = 12.05;


for coordx= 1:2048
      if X_f (coordx) == a && Y_f (coordx) == b && Z_f (coordx) == c
      x=coordx
          break
      end
end

% for coordx= 1:2048
%       if X_f (coordx) == a 
%       x=coordx
%           break
%       end
% end    
% % 
% for coordx= 1:2048
%       if  Y_f (coordx) == b 
%       x=coordx
%           break
%       end
% end
% % 
% for coordx= 1:2048
%       if Z_f (coordx) == c
%       x=coordx
%           break
%       end
% end