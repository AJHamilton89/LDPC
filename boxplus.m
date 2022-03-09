function c=boxplus(a,b)
        c=sign(a)*sign(b)*min(abs(a),abs(b))+log(1+exp(-abs(a+b)))-log(1+exp(-abs(a-b)));  % SPA
%         c=sign(a)*sign(b)*min(abs(a),abs(b));                                                % MSA
end