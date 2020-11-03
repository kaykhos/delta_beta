function ham = reduce_ham(Hin)
%REDUCE_HAM Summary of this function goes here
%   Detailed explanation goes here
    Hin = Hin - 11660*eye(size(Hin));
    ham = Hin(1:7,1:7);
end

