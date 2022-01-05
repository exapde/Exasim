#****************************************************************************                               
#                     Molecular Dynamics Potentials (MDP)
#                            CESMIX-MIT Project  
# 
#  Contributing authors: Ngoc-Cuong Nguyen (cuongng@mit.edu, exapde@gmail.com)
#****************************************************************************

function contains(ustr, s)

itsin = false;    
if occursin(s, string(ustr))
    itsin = true;
end

return itsin;

end

