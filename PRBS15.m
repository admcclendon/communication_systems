
classdef PRBS15 < PRBSGenerator
   methods
       function obj = PRBS15()
           seed = 1;
           g = [1 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0];
           
           obj@PRBSGenerator(g, seed);
       end
   end
end