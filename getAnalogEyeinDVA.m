function [eye_dva,bhv,eye_ana] = getAnalogEyeinDVA(nsfilename,bhvfilename)

bhv     = concatBHV(bhvfilename); 
tform   = bhv.EyeTransform; 

eye_ana      = getBNCData({'ainp2'},nsfilename,1); % in microV 
eye_ana(:,2) = getBNCData({'ainp3'},nsfilename,1);
eye_ana      = eye_ana./10E2; %convert to 10-4V

eye_dva = convertANAtoDVA(tform,eye_ana); 

function eye_dva = convertANAtoDVA(tform,eye_ana)
    eye_dva = [eye_ana ones(size(eye_ana,1),1)] * tform.tdata.T; 
    eye_dva = eye_dva(:,1:2) ./ repmat(eye_dva(:,3),1,2);
end

% function eye_ana = convertDVAtoANA(tform,eye_dva)
%     eye_ana = [eye_dva ones(size(eye_dva,1),1)] * tform.tdata.invT;
%     eye_ana = eye_ana(:,1:2) ./ repmat(eye_ana(:,3),1,2);
% end

end

%obj.tform{2}.origin = obj.tform{2}.origin + offset*obj.tform{2}.rotation_rev_t./obj.tform{2}.gain;
 
 
%  function q = fwd_projective(tform,p)
%     q = [p ones(size(p,1),1)] * tform.tdata.T;
%     q = q(:,1:2) ./ repmat(q(:,3),1,2);
% end
% function p = inv_projective(tform,q)
%     p = [q ones(size(q,1),1)] * tform.tdata.Tinv;
%     p = p(:,1:2) ./ repmat(p(:,3),1,2);
% end
