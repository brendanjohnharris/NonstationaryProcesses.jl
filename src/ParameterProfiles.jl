

# ------------------------------------------------------------------------------------------------ #
#                           Functions for constructing parameter profiles                          #
# ------------------------------------------------------------------------------------------------ #

function heaviside(x::Real, zeroval::Real=1.0)
    if x < 0
        y = 0
    elseif x > 0
        y = 1
    else
        y = convert(typeof(x), zeroval)
    end
end


