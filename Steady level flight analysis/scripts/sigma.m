function sig=sigma(h_)
% function to calculate relative density as a variation of altitude.
    sig=ones(size(h_));
    for i=1:1:length(h_)
        h=h_(i);
        if h<11000
            e=1;
            ha=0;
            B=9296;
        else
            e=exp(-11000/9296);
            ha=11000;
            B=6216;
        end
        sig(i)=e.*exp(-(h-ha)/B);
    end
end