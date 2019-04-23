function [x] = signedlog(x)
    epspow = 15;
    x(abs(x)<10^(-epspow))=0;
    x = sign(x).*(log10(abs(x))+epspow+1);
    x(~isfinite(x))=0;

end