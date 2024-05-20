
#= ==========================================================================================
=============================================================================================
functions
=============================================================================================
========================================================================================== =#

#---------------------------------------- kinematics ----------------------------------------

g = 9.81; # m/s²
potentialEnergy(pos::Tuple) = g * pos[3]; # method when poisition is provided
potentialEnergy(z::Real) = g * z; # method when height is provided 
speedFromKineticEnergy(K::Real) = (K < 0) ? (error("Kinetic energy became negative!")) : (√(2K))

#---------------------------------- cart auxilary functions --------------------------------#

"A cart is initialized with the appropriate type (`cartStruct`) taking the curve of its roller coaster."
function cartInit(curve::curveStruct; excessEnergy = 0.1)
    # the initial position of the cart is equal to the start of the track
    s = curve.start;
    pos = curve.r(s);
    # the total energy is defined
    max_height = [curve[3] for curve ∈ curve.r.(range(curve.start, stop = curve.stop, length = curve.length))] |> maximum;
    totalEnergy = excessEnergy + potentialEnergy(max_height);
    # the cart is initialized
    return cartStruct(s, pos, totalEnergy, curve)
end

"This function copies a cart onto a new cart (memory management)"
function copy(cart::cartStruct)
    newCart = cartInit(cart.curve)
    newCart.s = cart.s
    newCart.pos = cart.pos
    newCart.totalEnergy = cart.totalEnergy
    newCart.curve = cart.curve

    return newCart
end

"A cart is moved with a provided displacement along the track."
function cartMove!(cart::cartStruct, Δs::Float64)
    # s and r(s) are mutated for the cart
    cart.s += Δs;
    cart.pos = cart.curve.r(cart.s);
end

"The cart is moved with a provided time step using a 4-th order Runge Kutta integration scheme."
function cartRK4Step!(cart::cartStruct, Δt::Float64)
    # ds/dt(s) = F(s), F is defined
    F(s) = cart.totalEnergy - potentialEnergy(cart.curve.r(s)) |> speedFromKineticEnergy |> x -> x / dr_ds(cart.curve, s; Δs = 0.001)

    # 4-th order Runge Kutta is used to determine Δs
    k1 = Δt .* F(cart.s)
    k2 = Δt .* F(cart.s + k1/2)
    k3 = Δt .* F(cart.s + k2/2)
    k4 = Δt .* F(cart.s + k3)

    Δs = (k1 + 2k2 + 2k3 + k4) / 6

    # the cart is moved according to Δs
    cartMove!(cart, Δs);
end

"The cart is moved for a large amount of time; if stop is a number, then r(t) is found for t ∈ [0, stop], otherwise r(t) is found until the end of the track is reached"
function cartPropagate(cart::cartStruct; stop = true, Δt = 0.01)
    # the state of the cart will be recorded for each time stamp
    carts = Array{cartStruct}(undef, 0)
    if stop isa Bool
        # the time it takes to traverse the track is not known a priori; it is found adding Δt each step
        # it is initialized at -Δt since the last step corresponds to t = 0
        stop = -Δt;
        # the cart is moved for Δt until the end of the track is reached 
        while cart.s < cart.curve.stop
            # a new cart is appended to the list of carts
            append!(carts,[cartInit(cart.curve)]);
            # the cart copied in the recordings
            carts[end] = copy(cart);
            # the cart is moved
            cartRK4Step!(cart, Δt)
            # the stop time is updated
            stop += Δt
        end
        # the time is defined
        t = range(0, stop = stop, step = Δt)
    else
        # the time is defined
        t = range(0, stop = stop, step = Δt)
        # for each timestep, the cart is moved
        for _ ∈ eachindex(t)[1:end]
            # a new cart is appended to the list of carts
            append!(carts,[cartInit(cart.curve)]);
            # the cart copied in the recordings
            carts[end] = copy(cart);
            # the cart is moved
            cartRK4Step!(cart, Δt)
        end
    end
    
    # the time and the recordings are returned
    return t, carts
end

#------------------------------------ numerical analysis -----------------------------------#

"Calculates |dr/ds(s)| either analytically (`dr_dsCuve(s::Real)` must be defined for this to happen) or numerically (6th order finite differences are used)."
function dr_ds(curve::curveStruct, s::Real; Δs = 0.01)
    # the function attempts to use the analytical derivative; if it is not defined, it will
    # perform a numerical derivative
    try
        return dr_dsCurve(s) |> normTuple
    catch
        coeff = [45; -45; -9; 9; 1; -1];
        dr = [[curve.r(s + Δs), curve.r(s - Δs)] for Δs ∈ Δs:Δs:3Δs] |> v -> vcat(v...)
        dr = [coeff[i] .* dr[i] for i ∈ eachindex(dr)];
        
        return sumVectorOfTuples(dr) ./ (60Δs) |> normTuple
    end
end

"The matrix representation of the second order derivative is found using `t`; optionally, periodic boundary conditions may be used"
function Diff2(t; pbc = true)
    # f''(tᵢ) = (f(tᵢ₋₁) - 2f(tᵢ) + f(tᵢ₊₁))/Δt^2 is used

    # the matrix is found before dividing by Δt^2; it will be of the form
        #⌈-2   1            ⌉
        #| 1  -2   1        |
        #|     1  -2   1    |
        #|         1  -2   1|
        #⌊             1  -2⌋
    D2 = [-2 * (i == j) + (i == j+1) + (i+1 == j) for i ∈ eachindex(t), j ∈ eachindex(t)]
    
    # if periodic boundary conditions are not needed, the matrix is returned straight away
    (!pbc) ? (return D2 / step(t)^2) : nothing
    
    # the matrix is adjusted to consider periodic boundary conditions; it will be of the form
        #⌈-2   1           1⌉
        #| 1  -2   1        |
        #|     1  -2   1    |
        #|         1  -2   1|
        #⌊ 1           1  -2⌋
    D2[1,end] = D2[end,1] = 1;

    # the matrix is returned
    return D2 / step(t)^2
end

"An auxilary function was defined to sum the elemetns of arrays of tuples; they are of the form:
[(3,1,4), (1,5,9), (2,6,5)]
In this case the expected result is
(6, 12, 18)"
function sumVectorOfTuples(V)
    C = V[1];

    for v ∈ V[2:end]
        C = C .+ v
    end

    return C
end

"The norm of a tuple is found."
normTuple(T) = √(sum(T^2 for T ∈ T))

"A tuple is made unitary."
unitaryTuple(T) = T ./ normTuple(T)

import Base.*
"Multiplication of a matrix times a vector of 3-tuples is defined."
*(M::Matrix{Float64}, V::Vector{Tuple{Float64, Float64, Float64}}) = [[M[i,j] .* V[j] for j ∈ eachindex(V)] |> sumVectorOfTuples for i ∈ eachindex(V)]

#---------------------------------- functions for graphics ---------------------------------#

"The curve is plotted; this method simply plots the curve."
function plotCurve(curve::curveStruct; linewidth = 4, alternateView = false)
    r = curve.r.(range(curve.start, stop = curve.stop, length = curve.length));

    return lines(r,
        axis=(type=Axis3, 
            aspect = :data, 
            azimuth = (alternateView) ? (2.9π/4) : (5.1π/4),
        ),
        linewidth = linewidth,
    )
end

"The curve is plotted; this method colors the track according to the acceleration of the cart."
function plotCurve(carts::Vector{cartStruct}; linewidth = 4, pbcError = 0.01, colormap = :bamako, alternateView = false, profile = "abs")
    r = [cart.pos for cart ∈ carts]

    # it is checked whether the track is periodic or not to define the differential operator accordingly
    (r[1] .- r[end] |> normTuple <= pbcError) ? (pbc = true) : (pbc = false)
    D2 = Diff2(t; pbc = pbc)

    # the acceleration is found in units of g
    acceleration = D2 * [cart.pos for cart ∈ carts] .|> T -> T ./ g

    # for graphic stability at the borders, the initial and final accelerations are not used.
    # The track will be designed such it is smooth at all points, so this shouldn't cause problems.
    acceleration[1] = acceleration[2]
    acceleration[end] = acceleration[end-1]

    # according to the profile, the data to be plotted is selected
    if profile == "abs"
        acceleration = acceleration .|> normTuple
    elseif profile == "x"
        acceleration = acceleration .|> T -> T[1]
    elseif profile == "y"
        acceleration = acceleration .|> T -> T[2]
    elseif profile == "z"
        acceleration = acceleration .|> T -> T[3]
    end

    # the track is plotted 
    fig, ax, hm = lines(r,
        axis=(type=Axis3, 
        aspect = :data, 
        azimuth = (alternateView) ? (2.9π/4) : (5.1π/4),
        ),
        linewidth = linewidth,
        # colorrange = (0, 6), highclip = :red, # truncate the colormap in case the acceleration becomes too large 
        color = acceleration,
        colormap = colormap
    )

    title = Dict(
        "abs"   =>  "acceleration magnitude, in units of g",
        "x"     =>  "acceleration towards x, in units of g",
        "y"     =>  "acceleration towards y, in units of g",
        "z"     =>  "acceleration towards z, in units of g"
    )
    Colorbar(fig[1, 2], hm, label = title[profile])

    # the plot is returned
    return fig, ax
end

"The roller coaster is plotted with a red dot at the initial poisition of the cart."
function previewRollerCoaster(curve::curveStruct; linewidth = 4, alternateView = false)
    fig, ax = plotCurve(curve; linewidth = linewidth, alternateView = alternateView)
    scatter!(ax, cartInit(curve; excessEnergy = 0).pos, color = :red, markersize = 15)

    return fig, ax
end

"The animation of the cart moving along the roller coaster is created."
function anim8rollerCoaster(t, carts::Vector{cartStruct}; desired_fps = 30, linewidth = 4, plotAcceleration = true, pbcError = 0.01, alternateView = false)
    # for the animation, only the position of the cart will vary from frame to frame. This Observable will deal with that
    currentCartPos = Observable((0., 0., 0.)) 

    # the step and frame rate of the animation that can be achieved with the provided data are determined from a requested frame rate
    animStep = step(t) * desired_fps |> inv |> round |> Int64
    achieved_fps = t[1:animStep:end] |> step |> inv |> round |> Int64

    # it may or may not be needed to plot the acceleration throughout the plot
    if plotAcceleration
        animationFig, animationAx = plotCurve(carts; linewidth = linewidth, pbcError = pbcError, alternateView = alternateView)
    else
        animationFig, animationAx = plotCurve(carts[1].curve; linewidth = linewidth, alternateView = alternateView)
    end

    # a dot is used to represent the cart
    scatter!(animationAx, currentCartPos, color = :red, markersize = 15)

    # the animation is created; this function is essentially a for loop of the form:
    # for i ∈ 1:animStep:length(carts)
        # the animation is mutated for each frame
    # end
    record(animationFig, 
        "animations/$(today())/roller coaster dynamics ($(achieved_fps) fps).mp4", 
        1:animStep:length(carts); framerate = achieved_fps) do i
        # mutating the Observable declared above, the cart's dot is moved around the scene
        currentCartPos[] = carts[i].pos
    end
end