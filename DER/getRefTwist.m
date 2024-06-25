function refTwist = getRefTwist(d1, tangent, refTwist)

ne = length(refTwist);

for c=2:ne
    u0 = d1(c-1,:);
    u1 = d1(c,:);
    t = tangent(c,:);
    t0 = tangent(c-1,:);
    
    refTwist(c) = computeReferenceTwist(u0, u1, t0, t);
    
%     ut = parallel_transport(u0, t0, t);
%     % ut = rotateAxisAngle(ut, t, refTwist(c) );
%     
%     % get signed angle
%     sgnAngle = signedAngle(ut, u1, t);
%     % refTwist(c) = refTwist(c) + sgnAngle;
%     refTwist(c) = sgnAngle;
end

end
