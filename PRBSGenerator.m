
classdef PRBSGenerator < handle
    properties
        G = [];
        Seed = 1;
        Reg = [];
        BitMask = 1;
        BitShifts = [];
    end
    methods
        function obj = PRBSGenerator(g, seed)
            obj.G = g;
            obj.Seed = seed;
            
            obj.BitMask = 2^(length(g) - 1) - 1;
            obj.BitShifts = length(obj.G) - find(obj.G==1);
            
            obj.Reset();
        end
        
        function nextBits = Generate(obj, n)
            nextBits = zeros(1, n);
            for k = 1:n
                nextBit = mod(sum(bitand(bitshift(obj.Reg, -(obj.BitShifts - 1)), 1)), 2);
                obj.Reg = bitand(bitor(bitshift(obj.Reg, 1), nextBit), obj.BitMask);
                nextBits(k) = nextBit;
            end
        end
        
        function Reset(obj)
            obj.Reg = bitand(obj.Seed, obj.BitMask);
        end
    end
end