function gradient(f, h = 1.)
    n = length(x)
    df = zeros(n)

    # First point
    df[1] = (f[2] - f[1])/h
    df[end] = (f[end] - f[end-1])/h

    # Last point
    return df
end

g[1] = (F[2] - F[1]) / (h[2] - h[1])
        g[n] = (F[n] - F[n-1]) / (h[end] - h[end-1])
        if n > 2
            h = h[3:n] - h[1:n-2]
            g[2:n-1] = (F[3:n] - F[1:n-2]) ./ h
        end
    end