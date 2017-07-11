% x is a patch
% decide if x is a small vessel, if yes, strength it

function [isWide, indexpos, indexneg] = decideVesselType ( x, template )

% width = sqrt ( size(x,1) );
patch = reshape( x, 8,8);

for i = 1: size(template,2)
    indexposx = find(template{i});
    indexnegx = find ( ~template{i} );
    posE = sum(patch (indexposx).^2);
    negE = sum(patch (indexnegx).^2);
    allNorm(i) = (posE - negE) ;
    
end
for i=1:size(template,2)
    if(allNorm(i)==max(allNorm))
        idxx = i;
        break;
    end
end
% idxx = find ( allNorm == max ( allNorm ) );
if ( allNorm(idxx) > 0.001)
    isWide = 0;
else
    isWide = 1;
end

% if(find(idxx)>1)
%     idxx=idxx(1);
% end
% disp(['idxx=   ',num2str(idxx)]);
indexposx = find(template{idxx});
indexnegx = find ( ~template{idxx} );
indexpos = indexposx;
indexneg = indexnegx;

% figure(6);  
% subplot(121); imshow( reshape(x,8,8),[]); colorbar;
% subplot(122); imshow( template{idx},[]);
% title ( num2str(max ( allNorm )));
% pause;

