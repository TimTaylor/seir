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
    intervention_start_times = NULL,
    intervention_end_times = NULL,
    intervention = NULL,
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

    # calculate the cumulative intervention effect
    if (is.null(intervention_start_times) || is.null(intervention_end_times) || is.null(intervention)) {
        nint <- 0
        int_start <- NULL
        int_end <- NULL
        int <- NULL

    } else {
        tmp <- .cumulate_interventions(intervention_start_times, intervention_end_times, intervention)
        int_start <- tmp$start
        int_end <- tmp$end
        int <- unlist(tmp$out)
        nint <- length(int_start)
    }

    # bundle the parameters in to a vector
    p <- c(n, alpha, beta, gamma, vac_start_times, vac_end_times, nu, mat, nint, int_start, int_end, int)

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
    intervention_start_times = NULL,
    intervention_end_times = NULL,
    intervention = NULL,
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

    # calculate the cumulative intervention effect
    if (is.null(intervention_start_times) || is.null(intervention_end_times) || is.null(intervention)) {
        nint <- 0
    } else {
        tmp <- .cumulate_interventions(intervention_start_times, intervention_end_times, intervention)
        int_start <- tmp$start
        int_end <- tmp$end
        int <- tmp$out
        nint <- length(int_start)
    }

    # ode solver
    .ode <- function(t, y, params, ...) {

        S <- y[s_index]
        E <- y[e_index]
        I <- y[i_index]

        cm <- mat
        cm <- .intervention_on_cm(t,cm,int_start,int_end,int,nint)

        idx <- vac_start_times <= t & t < vac_end_times
        StoV <- nu * idx * S
        StoE <- beta * S * cm %*% I
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


.cumulate_interventions <- function(start, end, values) {
    if (!is.matrix(values))
        values <- matrix(values, nrow = length(values))
    bounds <- sort(unique(c(start, end)))
    out <- matrix(0, nrow = nrow(values), ncol = length(bounds) - 1L)
    for (j in seq_len(length(bounds) - 1L)) {
        lb <- bounds[j]
        ub <- bounds[j+1L]
        for (k in seq_along(start)) {
            s <- start[[k]]
            e <- end[[k]]
            if (lb >= s && ub <= e) {
                out[,j] <- out[,j] + values[,k]
            }
        }
    }
    out[out>1] <- 1
    out <- lapply(seq_len(ncol(out)), function(j) out[,j,drop=TRUE])
    list(start = bounds[-length(bounds)], end = bounds[-1], out = out)
}

# create function to adjust the contact matrix
.intervention_on_cm <- function(t, cm, int_start, int_end, int_mat, nint) {
    for (i in seq_len(nint)) {
        if (int_start[i] <= t && t < int_end[i]) {
            cm <- cm * (1 - int_mat[[i]])
            break
        }
    }
    cm
}

