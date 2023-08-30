#' @useDynLib seir
#' @importFrom deSolve ode
NULL

#' @export
seir_c <- function(
    alpha,
    beta,
    gamma,
    contact_matrix,
    demography_vector,
    initial_conditions,
    time_end,
    vac_start_times = NULL,
    vac_end_times = NULL,
    nu = NULL,
    increment = 1
) {
    # set the max number of parameters:
    #     - if increased, remember to increase in C src as well!
    MAX <- 500L

    # scale the contact matrix
    mat <- contact_matrix / max(Re(eigen(contact_matrix)$values))
    mat <- mat / demography_vector

    # pull out the number of age categories
    n <- nrow(contact_matrix)

    # scale the initial conditions to the population level and convert to vector
    init <- c(initial_conditions * demography_vector)

    # handle no vaccination
    if (is.null(vac_start_times) || is.null(vac_end_times) || is.null(nu))
        vac_start_times <- vac_end_times <- nu <- rep.int(0, nrow(contact_matrix))

    # bundle the parameters in to a vector
    p <- c(n, alpha, beta, gamma, vac_start_times, vac_end_times, nu, mat)

    # check the length of parameters
    if (length(p) > MAX)
        stop("Flattened parameters of length %d (only space for %d preallocated)", length(p), MAX)
    length(p) <- MAX

    ode(
        y = init,
        times = seq(from = 0, to = time_end, by = increment),
        func = "derivs",
        parms = p,
        dllname = "seir",
        initfunc = "initmod"
    )
}

#' @export
seir_r <- function(
    alpha,
    beta,
    gamma,
    contact_matrix,
    demography_vector,
    initial_conditions,
    time_end,
    vac_start_times = NULL,
    vac_end_times = NULL,
    nu = NULL,
    increment = 1
) {

    # scale the contact matrix
    mat <- contact_matrix / max(Re(eigen(contact_matrix)$values))
    mat <- mat / demography_vector

    # pull out the number of age categories
    n <- nrow(contact_matrix)

    # scale the initial conditions to the population level and convert to vector
    init <- c(initial_conditions * demography_vector)

    # handle no vaccination
    if (is.null(vac_start_times) || is.null(vac_end_times) || is.null(nu))
        vac_start_times <- vac_end_times <- nu <- rep.int(0, nrow(contact_matrix))

    # pull out the indices that will be used for the odes
    s_index <- 1:n
    e_index <- s_index + n
    i_index <- e_index + n
    r_index <- i_index + n

    # pull out the
    .ode <- function(t, y, params, ...) {

        S <- y[s_index]
        E <- y[e_index]
        I <- y[i_index]

        idx <- vac_start_times <= t & t < vac_end_times
        StoV <- nu * idx * S
        StoE <- beta * S * mat %*% I
        EtoI <- alpha * E
        ItoR <- gamma * I

        dS <- -StoE - StoV
        dE <- StoE - EtoI
        dI <- EtoI - ItoR
        dR <- ItoR

        list(c(dS, dE, dI, dR))
    }

    ode(
        y = init,
        times = seq(from = 0, to = time_end, by = increment),
        func = .ode,
        params = NULL
    )
}
