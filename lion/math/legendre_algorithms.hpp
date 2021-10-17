#ifndef __LEGENDRE_ALGORITHMS_HPP__
#define __LEGENDRE_ALGORITHMS_HPP__

//
//     ----------------------------------------------------------------------
//     Compute the  Legendre Polynomial by the three point recursion
//
//           L (x) = 2k-1 xL     - k-1 L
//            k      ----   k-1    ---  k-2
//                     k            k
//
//     Compute the Legendre Polynomial of degree k and its derivative
//     ------------------------------------------------------------------------
//
inline constexpr scalar legendre_polynomial(size_t N, scalar x )
{      
    if ( N == 0 ) return 1.0;
    else if ( N == 1 ) return x;
    else
    {
        scalar Lnm2 = 1.0;
        scalar Lnm1 = x;
        scalar Ln   = 0.0;

        for (size_t k = 2; k <= N; ++k )
        {
            Ln = ((2.0*k-1.0)*x*Lnm1 - (k-1.0)*Lnm2)/k;
            Lnm2 = Lnm1;
            Lnm1 = Ln;    
        } 

        return Ln;
    }
}

inline std::vector<scalar> legendre_polynomials(size_t N, scalar x )
{      
    if ( N == 0 ) return {1.0};
    else if ( N == 1 ) return {1.0,x};
    else
    {
        std::vector<scalar> result(N+1);
        result[0] = 1.0;
        result[1] = x;

        for (size_t k = 2; k <= N; ++k )
        {
            result[k] = ((2.0*k-1.0)*x*result[k-1] - (k-1.0)*result[k-2])/k;
        } 

        return result;
    }
}

//     Compute the Legendre Polynomial of degree k and its derivative
inline constexpr std::pair<scalar,scalar> legendre_poly_and_derivative(size_t N, scalar x)
{
    if ( N == 0 )  return {1.0, 0.0};
    else if ( N == 1 ) return {x,1.0};
    else
    {
        scalar Lnm2 = 1.0;   scalar dLnm2 = 0.0;
        scalar Lnm1 = x ;    scalar dLnm1 = 1.0;
        scalar Ln   = 0.0;   scalar dLn   = 0.0;
        for (size_t k = 2; k <= N; ++k)
        {
            Ln = ((2.0*k-1.0)*x*Lnm1 - (k-1.0)*Lnm2)/k;
            dLn = dLnm2 + (2.0*k-1.0)*Lnm1;

            Lnm2 = Lnm1;  dLnm2 = dLnm1;
            Lnm1 = Ln  ;  dLnm1 = dLn;
        }

        return {Ln,dLn};
    }
}


inline std::pair<std::vector<scalar>,std::vector<scalar>> gauss_legendre_nodes_and_weights(size_t N)
{
    const size_t n_newton_iterations = 10;
    const scalar tolerance_factor = 4.0;

    const scalar tolerance = tolerance_factor*eps;

    if ( N == 0 )
        return  {std::vector<scalar>{0.0}, std::vector<scalar>{2.0}} ;

    else if ( N == 1 )
        return {std::vector<scalar>{-std::sqrt(1.0/3.0),std::sqrt(1.0/3.0)}, std::vector<scalar>{1.0, 1.0}};

    else
    {
      
        const size_t n_div_2 = (N+1)/2;

        std::vector<scalar> nodes(N+1);
        std::vector<scalar> weights(N+1);

//      Iterate on half the interior nodes
//      ----------------------------------
        for (size_t j = 0; j < n_div_2; ++j)
        {
            scalar xj = -cos( (2.0*j+1.0)*pi/(2.0*N+2.0) );
            std::pair<scalar,scalar> Ln({0.0,0.0});

            for (size_t k = 0; k < n_newton_iterations; ++k)
            {
                Ln =  legendre_poly_and_derivative( N+1, xj);
                const scalar delta = -Ln.first/Ln.second;
                xj += delta;

                if ( std::abs(delta) <= tolerance*std::abs(xj) ) break;
            }
    
            nodes[j] = xj;
            weights[j] = 2.0/( (1.0-xj*xj)*Ln.second*Ln.second);
            nodes[N-j] = -nodes[j]; 
            weights[N-j] = weights[j];
        }
//
//      ---------------------------
//      Fill in middle if necessary
//      ---------------------------
//
        if ( (N % 2) == 0 )
        {
            std::pair<scalar,scalar> Ln = legendre_poly_and_derivative(N+1,0.0);
            nodes[N/2]   = 0.0;
            weights[N/2] = 2.0/(Ln.second*Ln.second);
        }

        return {nodes,weights};
    }
}


inline scalar lagrange_polynomial(scalar x, size_t j, size_t N, const std::vector<scalar>& xj)
{
    if ( j == 0 )
    {
        scalar p = (x - xj[1])/(xj[0] - xj[1]);
  
        for (size_t k = 2; k <= N; ++k)
            p *= (x - xj[k])/(xj[0] - xj[k]);
  
        return p;
    }
    else
    {
        scalar p = (x - xj[0])/(xj[j] - xj[0]);

        for ( size_t k = 1; k < j; ++k )
            p *= (x - xj[k])/(xj[j] - xj[k]);

        for ( size_t k = j+1; k <= N; ++k)
            p *= (x - xj[k])/(xj[j] - xj[k]);

        return p;
    }
}


inline std::pair<std::vector<scalar>,std::vector<scalar>> gauss_legendre_lobatto_nodes_and_weights(const size_t N)
{
    const size_t n_newton_iterations = 10;
    const scalar tolerance_factor = 4.0;
      
    const scalar tolerance = tolerance_factor*eps;

    if ( N == 0 )
        throw std::runtime_error("Order must be N>0 for GL points");

    else if ( N == 1 ) 
        return std::pair(std::vector<scalar>{-1.0,1.0}, std::vector<scalar>{1.0,1.0});

    else
    {
        std::vector<scalar> x(N+1), w(N+1);
    
        x[0] = -1.0;
        w[0] = 2.0/(N*(N+1.0));

        x[N] = 1.0;
        w[N] = w[0];

        const size_t n_div_2 = (N+1.0)/2.0;

        for (size_t j = 1; j < n_div_2; ++j)
        {
            scalar xj = -cos( (j+0.25)*pi/N - 3.0/(8.0*N*pi*(j+0.25)));

            for (size_t k = 0; k < n_newton_iterations; ++k)
            {
                auto q = q_and_L_evaluation(N,xj);
                const scalar delta = -std::get<0>(q)/std::get<1>(q);
                xj += delta;
    
                if ( std::abs(delta) <= tolerance*std::abs(xj) )
                    break;
            }

            auto l_n = legendre_polynomial(N, xj);

            x[j] = xj;
            w[j] = 2.0/(N*(N+1)*l_n*l_n);

            x[N-j] = -xj;
            w[N-j] = w[j];
        }
        // Fill in the middle if necessary
        if ( (N % 2) == 0 )
        {
            scalar l_n = legendre_polynomial(N, 0.0);
            x[N/2] = 0.0;
            w[N/2] = 2.0/(N*(N+1)*l_n*l_n);
        }

        return {x,w};
    }
}


inline std::tuple<scalar,scalar,scalar> q_and_L_evaluation(const size_t N, const scalar x)
{
    if ( N == 0 )
        throw std::runtime_error("q_and_L_evaluation cannot be called for N=0");

    else if ( N == 1 )
        throw std::runtime_error("q_and_L_evaluation cannot be called for N=1");

    else
    {
        scalar L_kM2 = 1.0;
        scalar L_prime_kM2 = 0.0;

        scalar L_kM1 = x;
        scalar L_prime_kM1 = 1.0;

        scalar L_k;
        scalar L_prime_k;

        for (size_t k = 2; k <= N; ++k)
        {
            L_k = ((2.0*k-1.0)*x*L_kM1 - (k-1.0)*L_kM2)/k;    
            L_prime_k = L_prime_kM2 + (2.0*k-1.0)*L_kM1;
            L_kM2 = L_kM1;
            L_kM1 = L_k;
            L_prime_kM2 = L_prime_kM1;
            L_prime_kM1 = L_prime_k;
        }

        const size_t k = N+1;
        L_k = ((2.0*k-1.0)*x*L_kM1 - (k-1.0)*L_kM2)/k;
        L_prime_k = L_prime_kM2 + (2.0*k-1.0)*L_kM1;

               /*   Q              Qprime             L_N */
        return { L_k-L_kM2, L_prime_k - L_prime_kM2, L_kM1 };
    }
}


#endif
