module Cournot

using ForwardDiff
using NLsolve

function simulateCournot(playerCount, variableCount, scoreWeight, score, costWeight, cost)
    q′ = []
    for i in 1:playerCount
        function f!(F, x)
            for j in 1:variableCount
                F[j] = scoreWeight[i][j] * ForwardDiff.derivative(score[i][j], x[j]) - costWeight[i][j] * ForwardDiff.derivative(cost[i][j], x[j])
            end
        end
        r = nlsolve(f!, zeros(Float32, variableCount), autodiff = :forward)
        @assert converged(r) "Error: nl solve not converging"
        push!(q′, r.zero)
    end
    return q′
end

function tScoreNormalizer(x)
    if x <= 0
        return 0
    elseif 0 < x && x < 1
        return x
    else
        return 1
    end
end

function tCostNormalizer(x)
    e = 2.7182818284590
    return 1 / (2 * (1 + e ^ (-(x-10))))
end

println(simulateCournot(
    2, 2,
    [[1, 2], [3, 4]],
    [[tScoreNormalizer, tScoreNormalizer], [tScoreNormalizer, tScoreNormalizer]],
    [[5, 6], [7, 8]],
    [[tCostNormalizer, tCostNormalizer], [tCostNormalizer, tCostNormalizer]]
))

end  # module GenderNeutralLavatories
