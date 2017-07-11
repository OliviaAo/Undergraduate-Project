function [A]=OMPerr(D,X,errorGoal,vecOfMeans); 
%=============================================
% Sparse coding of a group of signals based on a given 
% dictionary and specified number of atoms to use. 
% input arguments: D - the dictionary
%                  X - the signals to represent
%                  errorGoal - the maximal allowed representation error for
%                  each siganl.
% output arguments: A - sparse coefficient matrix.
%=============================================

%????????????????????????????????
%????????????????????
template{1} = uint8(255 - double(rgb2gray( imread('templatex.bmp','bmp'))));
for i = 1:7
    template{i+1} = imrotate(template{1}, 22.5*i,'crop');
end

[n,P]=size(X);
[n,K]=size(D);
% E2 = errorGoal^2*n;
% maxNumCoef = n/2;
E2 = errorGoal;
maxNumCoef = n;
% maxNumCoef=K;
A = sparse(size(D,2),size(X,2));
for k=1:1:P,
    a=[];
    x=X(:,k);
    residual=x;
	indx = [];
	a = [];
    
    [isWide,indexpos, indexneg] = decideVesselType ( x+vecOfMeans(k), template);
    
    if isWide  % wide vessel
        currResNorm2 = sum(residual.^2);
        j = 0;
        while currResNorm2>E2 & j < maxNumCoef,
            j = j+1;
            proj=D'*residual;
            pos=find(abs(proj)==max(abs(proj)));
            pos=pos(1);
            indx(j)=pos;
            a=pinv(D(:,indx(1:j)))*x;
            residual=x-D(:,indx(1:j))*a;
            currResNorm2 = sum(residual.^2);
        end;
        if (length(indx)>0)
             A(indx,k)=a;
        end
    else  % small vessel
        coef = 1/max(x);
        residual(indexpos) = residual(indexpos) * coef;
        currResNorm2 = sum(residual.^2);
        j = 0;

    %     figure(5); 
    %     subplot(141); imshow(reshape(x,8,8),[]); colorbar;  temp1=caxis;  title('orginal patch');

        while currResNorm2>E2/3 & j < maxNumCoef,
            j = j+1;
            proj=D'*residual;
            pos=find(abs(proj)==max(abs(proj)));
            pos=pos(1);
            indx(j)=pos;
            a=pinv(D(:,indx(1:j)))*x;
            residual=x-D(:,indx(1:j))*a;

            residual(indexpos) = residual(indexpos) * coef;
             currResNorm2 = sum(residual.^2);
    %         subplot(142); imshow(reshape(residual,8,8),[]); colorbar;  title('residual');
    %         subplot(143); imshow(reshape(D(:,indx(j)),8,8),[]);  colorbar;  title('current dictionary');
        end;
        if (length(indx)>0)
             A(indx,k)=a;
        end
        
    end

end;
return;
