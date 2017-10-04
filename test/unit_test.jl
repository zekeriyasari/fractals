"""
UnitTest module for unit tests.
"""
module UnitTest

export EPSILON, pass, fail, equal, almost_equal, raises

EPSILON = 1e-6
# EPSILON = 100 * eps()


"""
    pass()

Indicated that the test has passed.
"""
function pass()
end  # End of pass function


"""
    fail(message::AbstractString="")

Indicates that the test has failed.
"""
function fail(message::AbstractString="")
    throw(ErrorException("Test failed! $message"))
end  # End of fail function


"""
    equal(a::Real, b::Real)

Checks if `a` is equal to `b`
"""
function equal(a::Any, b::Any)
    return a == b
end


"""
    almost_equal(a, b; tol=EPSILON)

Checks that if ``
"""
function almost_equal(a::Real, b::Real; tol::Real=EPSILON)
    if tol <= 0
        throw(ArgumentError("tol must be positive."))
    end
    abs(a - b) <= tol
end


"""
    raises(func::Function, exc::DataType, args...; kwargs...)

Checks if `exc` is raised when `func(args...; kwargs...)` is called.
"""
function raises(func::Function, exc::DataType, args...; kwargs...)
    try
        func(args...; kwargs...)
    catch call_exc
        isa(call_exc, exc)
        return true
    end
    return false
end  # End of raises function

end  # End of UnitTest module
