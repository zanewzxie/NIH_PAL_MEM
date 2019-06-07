function zr=fisherz(r)

r(r==1)=0.99;
zr = .5.*[log(1+r)-log(1-r)];

end