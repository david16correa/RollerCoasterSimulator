
"The limits of the sections of the track are defined."
s_lims = [  0, # 1
            2.5, # 2
            10, # 3
            30π + 10, # 4
            30π + 17.5, # 5
            37.5π + 17.5, # 6
            37.5π + 18.5, # 7
            37.5π + 37.5, # 8
            45π + 37.5, # 9
            45π + 40 
]
"The initial and final s of the track are found."
s_lb, s_ub = s_lims[1], s_lims[end]

"The curve r(s) of the track is defined as a piecewise function."
function curve(s::Real)
    if s < s_lims[1] # 0
        # Due to numerical differentiation, r(-Δs) might be needed.
        # This is okay for small values of Δs, but a warning is thrown 
        # to alert the user anyways.
        @warn "r(-Δs) - $(s) is not in range! r ∈ ($(s_lims[1]), $(s_lims[end]))"
        (
            22.5 - s, 
            15, 
            0.5cos(2π/5 * s) + 10.5
        )
    elseif s_lims[1] <= s <= s_lims[2] # 1
        (
            22.5 - s, 
            15, 
            0.5cos(2π/5 * s) + 10.5
        )
    elseif s_lims[2] < s <= s_lims[3] # 2
        (
            22.5 - s, 
            15, 
            10
        )
    elseif s_lims[3] < s <= s_lims[4] # 3
        (
            -5*sin((s - 10)/5) + 12.5, 
            -5*cos((s - 10)/5) + 20, 
            5*cos((s - 10)/30)+5
        )
    elseif s_lims[4] < s <= s_lims[5] # 4
        (
            22.5 + 30π - s,
            15,
            0
        )
    elseif s_lims[5] < s <= s_lims[6] # 5
        (
            -7.5*sin(s/7.5 - 30π/7.5 - 7/3) + 5,
            7.5*cos(s/7.5 - 30π/7.5 - 7/3) + 7.5,
            -5*cos(s/7.5 - 30π/7.5 - 7/3) + 5
        )
    elseif s_lims[6] < s <= s_lims[7] # 6
        (
            s - 37.5π - 12.5,
            0,
            10
        )
    elseif s_lims[7] < s <= s_lims[8] # 7
        (
            s - 37.5π - 12.5,
            0,
            3.5*cos(2π/19*(s - 37.5π - 18.5)) + 6.5
        )
    elseif s_lims[8] < s <= s_lims[9] # 8
        (
            7.5*sin((s - 37.5π - 37.5)/7.5)+25,
            -7.5*cos((s - 37.5π - 37.5)/7.5)+7.5,
            10
        )
    elseif s_lims[9] < s <= s_lims[10] # 9
        (
            -s + 45π + 62.5,
            15,
            -0.5cos(2π/5 * (s - 45π - 37.5)) + 10.5
        )
    else
        # Due to numerical differentiation, r(s+Δs) might be needed.
        # Moreover, when determining the time needed to traverse the track,
        # the cart may exceed s_lb by a small amount.
        # This is okay for small values of Δs, but a warning is thrown 
        # to alert the user anyways.
        @warn "r(S + Δs) - $(s) is not in range! r ∈ ($(s_lims[1]), $(s_lims[end]))"
        (
            -s + 45π + 62.5,
            15,
            -0.5cos(2π/5 * (s - 45π - 37.5)) + 10.5
        )
    end
end

"dr/ds is defined as a piecewise function."
function dr_dsCurve(s::Real)
    if s < s_lims[1] # 0
        # There's no reason for dr_ds(-Δs) to be needed, but this scenario is considered just in case. 
        @warn "dr/ds(-Δs) - $(s) is not in range! r ∈ ($(s_lims[1]), $(s_lims[end])). This shouldn't happen."
        (
            -1, 
            0, 
            - 2π/10 * sin(2π/5 * s)
        )
    elseif s_lims[1] <= s <= s_lims[2] # 1
        (
            -1, 
            0, 
            - 2π/10 * sin(2π/5 * s)
        )
    elseif s_lims[2] < s <= s_lims[3] # 2
        (
            -1, 
            0, 
            0
        )
    elseif s_lims[3] < s <= s_lims[4] # 3
        (
            -cos((s - 10)/5), 
            sin((s - 10)/5),
            -5/30*sin((s - 10)/30)
        )
    elseif s_lims[4] < s <= s_lims[5] # 4
        (
            -1,
            0,
            0
        )
    elseif s_lims[5] < s <= s_lims[6] # 5
        (
            -cos(s/7.5 - 30π/7.5 - 7/3),
            -sin(s/7.5 - 30π/7.5 - 7/3),
            5/7.5*sin(s/7.5 - 30π/7.5 - 7/3)
        )
    elseif s_lims[6] < s <= s_lims[7] # 6
        (
            1,
            0,
            0
        )
    elseif s_lims[7] < s <= s_lims[8] # 7
        (
            1,
            0,
            -7π/19*sin(2π/19*(s - 37.5π - 18.5))
        )
    elseif s_lims[8] < s <= s_lims[9] # 8
        (
            cos((s - 37.5π - 37.5)/7.5),
            sin((s - 37.5π - 37.5)/7.5),
            0
        )
    elseif s_lims[9] < s <= s_lims[10] # 9
        (
            -1,
            0,
            2π/10*sin(2π/5 * (s - 45π - 37.5))
        )
    else # 9
        # When determining the time needed to traverse the track,
        # the cart may exceed s_lb by a small amount.
        # This is okay for small values of Δs, but a warning is thrown 
        # to alert the user anyways.
        @warn "dr/ds(S + Δs) - $(s) is not in range! r ∈ ($(s_lims[1]), $(s_lims[end]))"
        (
            -1,
            0,
            2π/10*sin(2π/5 * (s - 45π - 37.5))
        )
    end
end