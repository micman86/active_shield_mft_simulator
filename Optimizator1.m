classdef Optimizator1
   properties
      new_factor;
      old_factor;
      optimization_factor;
      optimization_factor_last;
      num_loops;
      variable_to_optimize;
      counter_change_sign;
      range_optimization_search_sigma; %percent value of sigma for search of optimum value
       num_of_variable_to_optimize; %modulus and angle of SIGMA
      last_evaluation_index;
      evaluation_index;
      
   end
   
   methods
       function obj = Optimizator1(numLoops)
          obj.num_loops = numLoops;
          
          obj.new_factor = 1;
          obj.old_factor=1;
          obj.optimization_factor(1:2*obj.num_loops) = 100*ones(2*obj.num_loops,1); %percent of modulus and angle
           obj.optimization_factor_last(1:2*obj.num_loops) = 100*ones(2*obj.num_loops,1); %percent of modulus and angle
          obj.variable_to_optimize=1;
          obj.counter_change_sign=1;
           obj.range_optimization_search_sigma = 90; %percent value of sigma for search of optimum value
       obj.num_of_variable_to_optimize = 2*obj.num_loops; %modulus and angle of SIGMA
       obj.last_evaluation_index = 0;
       

       end
       
       function obj = step(obj, eval_index)
            
            obj.evaluation_index = eval_index;
           
            if (obj.counter_change_sign > 4)
                    obj.counter_change_sign = 0;
                    obj.new_factor = 1;
                    obj.old_factor = 0;
                    obj.variable_to_optimize = obj.variable_to_optimize + 1;
                    if(obj.variable_to_optimize > obj.num_of_variable_to_optimize) 
                        obj.variable_to_optimize = 1;
                    end
            end
                    
            
             if (obj.old_factor == 0)
                    obj.old_factor = 1;
                    obj.new_factor = 1;
             else
                    if obj.evaluation_index < obj.last_evaluation_index
                        obj.new_factor = 1;
                    else
                        obj.new_factor = -1;
                    end
                    obj.last_evaluation_index = obj.evaluation_index;
                   
             end
                
             
              if(obj.new_factor == -1) 
                    obj.counter_change_sign = obj.counter_change_sign +1;
                end
                
                
                for u=1 : obj.num_of_variable_to_optimize
                    if u==obj.variable_to_optimize

                obj.optimization_factor(u) = obj.optimization_factor_last(u) + obj.old_factor*obj.new_factor*2;
                obj.old_factor=obj.old_factor*obj.new_factor;
                    else
                obj.optimization_factor(u) = obj.optimization_factor_last(u);        
                    end
                end
                
                
                 if obj.optimization_factor(obj.variable_to_optimize) > 100+obj.range_optimization_search_sigma 
                        obj.optimization_factor(obj.variable_to_optimize) = 100+obj.range_optimization_search_sigma;
                 elseif obj.optimization_factor(obj.variable_to_optimize) < 100-obj.range_optimization_search_sigma 
                        obj.optimization_factor(obj.variable_to_optimize) = 100-obj.range_optimization_search_sigma;
                 end
                    
                 if obj.optimization_factor(obj.variable_to_optimize) < 0 
                        obj.optimization_factor(obj.variable_to_optimize) = 0;
                 end
                     obj.optimization_factor_last = obj.optimization_factor;
                     
                  
           
       end
      
   end
end