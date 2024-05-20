include("preamble.jl")
include("structs.jl")
include("functions.jl")
include("geometryEg.jl")


# the roller coaster is defined
rollerCoaster = curveStruct(s_lb, s_ub, 1000, s -> curve(s))

# preview
begin
    fig, ax = previewRollerCoaster(rollerCoaster, alternateView = false)
    save("figures/$(today())/roller coaster preview.png", fig) 
    fig
end

# a cart is initialized with an excess energy of 1 Joule, and its dynamics are found for the whole track
t, carts = cartInit(rollerCoaster; excessEnergy = 1) |> cartPropagate


# acceleration plot - magnitude
begin
    fig, ax = plotCurve(carts, alternateView = false, profile = "abs")
    save("figures/$(today())/acceleration magnitude.png", fig) 
    fig, ax = plotCurve(carts, alternateView = true, profile = "abs")
    save("figures/$(today())/acceleration magnitude 2.png", fig) 
end
# acceleration plot - x
begin
    fig, ax = plotCurve(carts, alternateView = false, profile = "x")
    save("figures/$(today())/acceleration towards x.png", fig) 
    fig, ax = plotCurve(carts, alternateView = true, profile = "x")
    save("figures/$(today())/acceleration towards x 2.png", fig) 
end
# acceleration plot - y
begin
    fig, ax = plotCurve(carts, alternateView = false, profile = "y")
    save("figures/$(today())/acceleration towards y.png", fig) 
    fig, ax = plotCurve(carts, alternateView = true, profile = "y")
    save("figures/$(today())/acceleration towards y 2.png", fig) 
end
# acceleration plot - z
begin
    fig, ax = plotCurve(carts, alternateView = false, profile = "z")
    save("figures/$(today())/acceleration towards z.png", fig) 
    fig, ax = plotCurve(carts, alternateView = true, profile = "z")
    save("figures/$(today())/acceleration towards z 2.png", fig) 
end

# animation
anim8rollerCoaster(t, carts; desired_fps = 30, alternateView = false)
anim8rollerCoaster(t, carts; desired_fps = 30, alternateView = true)
