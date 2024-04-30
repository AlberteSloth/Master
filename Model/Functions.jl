using ModelingToolkit, DifferentialEquations, Random, Distributions, Plots, CSV, DataFrames

function SimData(odeSys, u0, Cons, L, V, K, tspan)
    p = [Cons 
        L 
        V 
        K]
    prob = ODEProblem(sys, u0, tspan, p, jac = true)
    sol = solve(prob)

    flux_ss = flux_vector(sol[end], p)
    return sol[end], flux_ss
end


function flux_vector(u, p_flux)
    if typeof(p_flux) == Vector{Pair{Num, Float64}}
        flux_01 = zeros(length(p_flux))
        for i in 1:length(p_flux)
            flux_01[i] = p_flux[i][2]
        end
    else
        flux_01 = p_flux
    end

    if typeof(u) == Vector{Pair{Num, Float64}}
        u_01 = zeros(length(u))
        for i in 1:length(u)
            u_01[i] = u[i][2]
        end
    else
        u_01 = u
    end 
    
    C_ATP, C_ADP, C_cit, C_AMP, C_PFKM, C_AMPK, C_Pi, C_G6P, C_GAP, C_DHAP, L_PFKM, L_PFK2, Vf_GPI, Vr_GPI, Vf_PFK2, Vr_PFK2, Vf_FBP, Vr_FBP, Vf_ALD, Vr_ALD, Kf_GPI, Kr_GPI, Ki_ATP, Ki_cit, Ka_f6p, Ka_f16bp, Ka_AMP, Ka_f26bp, Kcf_PFKM, Kcr_PFKM, K_f6p, K_ATP, K_f16bp, K_ADP, K_cat, v_P, K2_ATP, K2_f6p, K2_f26bp, K2_ADP, K_FBP_f26bp, K_FBP_f6p, K_Pi, KAld_f16bp, K_GAP, K_DHAP = flux_01
    f6p, f16bp, f26bp = u_01

    r_GPI = (Vf_GPI*C_G6P/Kf_GPI - Vr_GPI*f6p/Kr_GPI)/(1 + C_G6P/Kf_GPI + f6p/Kr_GPI)

    N_PFKM = 1 + L_PFKM*(1 + C_ATP/Ki_ATP)^4*(1 + C_cit/Ki_cit)^4/
        ((1 + f6p/Ka_f6p + f16bp/Ka_f16bp)^4 * (1 + C_AMP/Ka_AMP)^4 * (1 + f26bp/Ka_f26bp)^4)

    r_PFKM = C_PFKM*((Kcf_PFKM*C_ATP*f6p/(K_f6p*K_ATP) - Kcr_PFKM*C_ADP*f16bp/(K_f16bp*K_ADP))/
        ((1 + f6p/K_f6p)*(1 + C_ATP/K_ATP) + (1 + f16bp/K_f16bp)*(1 + C_ADP/K_ADP) -1))*
        (1/N_PFKM)

    psi = (K_cat*C_AMPK)/((K_cat*C_AMPK)+v_P)
    N_PFK2 = 1 + L_PFK2*(psi/(1-psi))^2

    r_PFK2 = ((Vf_PFK2*C_ATP*f6p)/(K2_ATP*K2_f6p) - (Vr_PFK2*C_ADP*f26bp)/(K2_f26bp*K2_ADP)) / 
        ((1+f6p/K2_f6p)*(1+C_ATP/K2_ATP) + (1+f26bp/K2_f26bp)*(1+C_ADP/K2_ADP) -1) *
        (1-(1/N_PFK2))

    r_FBP = ((Vf_FBP*f26bp)/(K_FBP_f26bp) - (Vr_FBP*C_Pi*f6p)/(K_Pi*K_FBP_f6p))/
        ((1+f26bp/K_FBP_f26bp) + (1+f6p/K_FBP_f6p)*(1+C_Pi/K_Pi)-1) * 
        (1/N_PFK2)

    r_ALD = ((Vf_ALD*f16bp)/(KAld_f16bp) - (Vr_ALD*C_GAP*C_DHAP)/(K_GAP*K_DHAP))/ 
        ((1 + f6p/KAld_f16bp) + (1 + C_GAP/K_GAP)*(1 + C_DHAP/K_DHAP) -1)
    
    f6p_flux = r_GPI - r_PFKM - r_PFK2 + r_FBP
    f16bp_flux = r_PFKM - r_ALD
    f26bp_flux = r_PFK2 - r_FBP
    
    df = DataFrame(a = ["r_GPI", "r_PFKM", "r_PFK2", "r_FBP", "r_ALD", "f6p_flux", "f16bp_flux", "f26bp_flux"], values = [r_GPI, r_PFKM, r_PFK2, r_FBP, r_ALD, f6p_flux, f16bp_flux, f26bp_flux])

    return[df
    ]
end 


function flux_vec_eq(u, p_flux)
    if typeof(p_flux) == Vector{Pair{Num, Float64}}
        flux_01 = zeros(length(p_flux))
        for i in 1:length(p_flux)
            flux_01[i] = p_flux[i][2]
        end
    else
        flux_01 = p_flux
    end

    if typeof(u) == Vector{Pair{Num, Float64}}
        u_01 = zeros(length(u))
        for i in 1:length(u)
            u_01[i] = u[i][2]
        end
    else
        u_01 = u
    end 
    
    C_ATP, C_ADP, C_cit, C_AMP, C_PFKM, C_AMPK, C_Pi, C_G6P, C_GAP, C_DHAP, L_PFKM, L_PFK2, Vm_GPI, Vm_PFK2, Vm_FBP, Vm_ALD, Keq_GPI, Keq_PFKM, Keq_PFK2, Keq_FBP, Keq_ALD, Kf_GPI, Kr_GPI, Ki_ATP, Ki_cit, Ka_f6p, Ka_f16bp, Ka_AMP, Ka_f26bp, Kc_PFKM, K_f6p, K_ATP, K_f16bp, K_ADP, Kc_AMPK, v_P, K2_ATP, K2_f6p, K2_f26bp, K2_ADP, K_FBP_f26bp, K_FBP_f6p, K_Pi, KAld_f16bp, K_GAP, K_DHAP = flux_01
    f6p, f16bp, f26bp = u_01

    r_GPI = (Vm_GPI/Kf_GPI)*(C_G6P-(f6p/Keq_GPI))/
        (1 + C_G6P/Kf_GPI + f6p/Kr_GPI)
    
    N_PFKM = 1 + L_PFKM*(1 + C_ATP/Ki_ATP)^4*(1 + C_cit/Ki_cit)^4/
        ((1 + f6p/Ka_f6p + f16bp/Ka_f16bp)^4 * (1 + C_AMP/Ka_AMP)^4 * (1 + f26bp/Ka_f26bp)^4)

    r_PFKM = ((C_PFKM*(Kc_PFKM/(K_f6p*K_ATP))*(C_ATP*f6p - C_ADP*f16bp/Keq_PFKM))/
        ((1 + f6p/K_f6p)*(1 + C_ATP/K_ATP) + (1 + f16bp/K_f16bp)*(1 + C_ADP/K_ADP) -1))*
        (1/N_PFKM)

    psi = (Kc_AMPK*C_AMPK)/((Kc_AMPK*C_AMPK)+v_P)
    N_PFK2 = 1 + L_PFK2*(psi/(1-psi))^2

    r_PFK2 = (Vm_PFK2/(K2_ATP*K2_f6p))*((C_ATP*f6p)-(C_ADP*f26bp)/Keq_PFK2)/
        ((1+f6p/K2_f6p)*(1+C_ATP/K2_ATP) + (1+f26bp/K2_f26bp)*(1+C_ADP/K2_ADP) -1) *
        (1-(1/N_PFK2))

    #(Vm_FBP/(K_FBP_f26bp))*(f26bp-(f6p*C_Pi)/Keq_FBP)/

    r_FBP = (Vm_FBP/(K_FBP_f26bp))*(f26bp-(f6p*C_Pi)/Keq_FBP)/
        ((1+f26bp/K_FBP_f26bp) + (1+f6p/K_FBP_f6p)*(1+C_Pi/K_Pi)-1) * 
        (1/N_PFK2)


    r_ALD = (Vm_ALD/KAld_f16bp)*(f16bp-(C_GAP*C_DHAP)/Keq_ALD)/
        ((1 + f6p/KAld_f16bp) + (1 + C_GAP/K_GAP)*(1 + C_DHAP/K_DHAP) -1)
    
    f6p_flux = r_GPI - r_PFKM - r_PFK2 + r_FBP
    f16bp_flux = r_PFKM - r_ALD
    f26bp_flux = r_PFK2 - r_FBP
    
    df = DataFrame(a = ["r_GPI", "r_PFKM", "r_PFK2", "r_FBP", "r_ALD", "f6p_flux", "f16bp_flux", "f26bp_flux", "psi", "N_PFK2", "N_PFKM"], values = [r_GPI, r_PFKM, r_PFK2, r_FBP, r_ALD, f6p_flux, f16bp_flux, f26bp_flux, psi, N_PFK2, N_PFKM])
    return[df
    ]

end 