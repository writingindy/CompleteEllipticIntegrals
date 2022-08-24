function AGM(a₁ , g₁ , maxIter = 100)
    aₙ = 0.5*(a₁ + g₁)
    gₙ = sqrt(a₁*g₁)

    error = abs(aₙ - gₙ)

    i = 1
    while i < maxIter && error > 1e-15

        aₙ₊₁ = 0.5*(aₙ + gₙ)
        gₙ₊₁ = sqrt(aₙ*gₙ)

        aₙ = aₙ₊₁ 
        gₙ = gₙ₊₁

        error = abs(aₙ - gₙ)

        i += 1
    end

    if (i == maxIter)
        print("AGM function hit maximum iteration limit", maxIter, ", check printed error to ensure convergence")
    end
    
    #print(error)

    return(aₙ)
end


function elliptic_K(m, maxIter = 100)
    return(π / (2*AGM(1, sqrt(1-m), maxIter)))
end

function elliptic_E(m, maxIter = 50)
    
    K = elliptic_K(m, maxIter)


    a₀ = 1
    g₀ = sqrt(1- m)

    a₁ = 0.5*(a₀ + g₀)
    g₁ = sqrt(a₀*g₀)

    aₙ = 0.5*(a₁ + g₁)
    gₙ = sqrt(a₁*g₁)

    cₙ = sqrt(aₙ^2 - gₙ^2)
    
    partialSum = 2*cₙ^2

    error = abs(aₙ - gₙ)

    i = 2
    while i < maxIter && error > 1e-15

        aₙ₊₁ = 0.5*(aₙ + gₙ)
        gₙ₊₁ = sqrt(aₙ*gₙ)

        cₙ₊₁ = cₙ^2/(4*aₙ₊₁)

        aₙ = aₙ₊₁ 
        gₙ = gₙ₊₁

        cₙ = cₙ₊₁

        partialSum += (2^i)*cₙ^2

        error = abs(aₙ - gₙ)

        i += 1
    end

    return(K*(a₁^2 - partialSum))
end

function elliptic_Π(α², m, maxIter = 50)
    M = AGM(1, sqrt(1- m), maxIter)


    a₀  = 1
    g₀  = sqrt(1 - m)
    p₀² = 1 - α²
    ε₀  = (p₀² - a₀*g₀)/(p₀² + a₀*g₀)
    Q₀  = 1

    aₙ = 0.5*(a₀ + g₀)
    gₙ = sqrt(a₀*g₀)
    pₙ = (p₀² + a₀*g₀)/(2*sqrt(p₀²))
    Qₙ = (1/2)*Q₀*ε₀

    partialSum = Q₀ + Qₙ

    error = abs(aₙ - gₙ)

    i = 1
    while i < maxIter && error > 1e-15

        aₙ₊₁ = 0.5*(aₙ + gₙ)
        gₙ₊₁ = sqrt(aₙ*gₙ)
        pₙ₊₁ = (pₙ^2 + aₙ*gₙ)/(2*pₙ)
        εₙ   = (pₙ^2 - aₙ*gₙ)/(pₙ^2 + aₙ*gₙ)
        Qₙ₊₁ = (1/2)*Qₙ*εₙ

        aₙ = aₙ₊₁ 
        gₙ = gₙ₊₁
        pₙ = pₙ₊₁
        Qₙ = Qₙ₊₁

        partialSum += Qₙ

        error = abs(aₙ - gₙ)

        i += 1
    end

    return((π/(4*M))*(2 + (α²/(1 - α²))*partialSum))
end

