function [kf_r, kb_r] = rateCoefficients(lnkf_r,lnkb_r)
    kf_r = exp(lnkf_r);
    kb_r = exp(lnkb_r);
end