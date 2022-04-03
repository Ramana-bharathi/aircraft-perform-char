function h=siginv(sig)
% function to calculate altitude as a variation of relative density.

    assert(sig<=1 & sig>=0);
    if sig>exp(-11000/9296)
        h=-9296*(log(sig));
    else
        h=11000-6216*(log(sig)-log(0.3063));
    end
end