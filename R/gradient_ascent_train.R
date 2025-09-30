#' @title Gradient Ascent Optimization for Model Training
#'
#' @description
#' This function performs gradient ascent optimization to maximize the objective function for training the gene regulatory network prediction model.
#' It iteratively updates model parameters using gradient information and adaptive learning rate selection to find optimal parameter values that maximize the model's predictive performance.
#'
#' @param iterations Iterations Integer. Maximum number of gradient ascent iterations to perform.
#' @param alpha_lower Numeric. Lower bound for learning rate search range. Default is 1e-5.
#' @param alpha_upper Numeric. Upper bound for learning rate search range. Default is 1e-3.
#' @param alpha_guess Numeric. Initial guess for learning rate in the optimization from alpha_lower to alpha_upper. Default is 1e-4.
#' @param init_theta
#' Numeric vector. Initial parameter values for optimization.
#' If NULL, initializes with zeros based on gradient function structure. Default is NULL.
#' @param h_max Numeric. Maximum step size for numerical gradient computation. Default is 0.01.
#' @param h_min_percent Numeric. Minimum step size as percentage of parameter range for numerical gradient. Default is 0.01 (1%).
#' @param setseed See parameter set_seed in ?spilt_dataset.
#' @param ncores
#' See in ?get_interest_cell_type_data.
#' While this code can be executed on Windows systems, it is still strongly recommended to run it on Linux systems, as the execution time without parallel computing would be astronomically long.
#' The actual runtime depends on the total available RAM, the number of CPU cores, the previously determined number of cell groups, and the number of transcription factors (TFs) and target genes (TGs) obtained.
#' Larger core counts, fewer cell groups, and fewer genes correspond to shorter runtime, while the total available RAM dictates the maximum number of cores that can participate in parallel computing.
#'
#' @returns
#' A list containing two elements:
#' \itemize{
#' \item \code{parameters}: Optimized parameter vector that maximizes the objective function
#' \item \code{objective_value}: The maximized objective function value
#' }
#' The optimization process stops early if no improvement is detected in consecutive iterations.
#'
#' @details
#' The function implements the following optimization strategy:
#' \enumerate{
#' \item Initialization: Sets initial parameters (zeros if not provided)
#' \item Gradient Computation: Uses numerical gradient approximation with parallel processing
#' \item Learning Rate Optimization: Finds optimal step size using Brent's method for each iteration
#' \item Parameter Update: Updates parameters with gradient ascent, respecting parameter bounds
#' \item Convergence Check: Stops early if objective function fails to improve
#' }
#' The algorithm uses adaptive learning rate selection and parameter bounding to ensure stable convergence.
#' The optimization maximizes the objective function defined in \code{f_train_cal}.
#'
#' @note
#' Important considerations:
#' \itemize{
#' \item Requires global variables: \code{branch_params_vec_lower} and \code{branch_params_vec_upper} for parameter bounding
#' \item Uses numerical gradient approximation which may be less accurate than analytical gradients
#' \item The optimization process is computationally intensive, especially with many parameters
#' \item Early stopping helps prevent overfitting and saves computation time
#' \item Learning rate bounds may need adjustment for different problem domains
#' \item Parallel processing significantly speeds up gradient computation
#' }
#'
#' @seealso \code{\link{f_train_cal}}, \code{\link{f_train_cal_gr}}, \code{\link[stats]{optim}}
#'
#' @export
#'
#' @examples
#' \dontrun{
#' iterations_grad = 3
#' setseed = 123 # Reproducible
#' # setseed = stats::runif(1) # Not Reproducible
#' ncores = parallel::detectCores() - 1 # in Linux
#' # ncores = 1 # in Windows
#' t = base::Sys.time(); result_gradient = gradient_ascent_train(
#'   iterations = iterations_grad,
#'   alpha_lower = 1e-5,
#'   alpha_upper = 1e-3,
#'   alpha_guess = 1e-4,
#'   init_theta = as.numeric(result_ga@solution[1,]),
#'   h_max = 1e-9,
#'   h_min_percent = 0.01,
#'   setseed = setseed,
#'   ncores = ncores
#' ); base::print(base::Sys.time(); - t)
#' }
gradient_ascent_train = function(iterations,
                                 alpha_lower = 1e-5,
                                 alpha_upper = 1e-3,
                                 alpha_guess = 1e-4,
                                 init_theta = NULL,
                                 h_max = 0.01,
                                 h_min_percent = 0.01,
                                 setseed = stats::runif(1),
                                 ncores = 1)
{
  # Initialing.
  if (base::is.null(init_theta)) {
    temp_grad = scDEDS::f_train_cal_gr(branch_params_vec = numeric(0),
                                       h_max = h_max,
                                       h_min_percent = h_min_percent,
                                       setseed = setseed,
                                       ncores = ncores)
    theta = base::rep(0, base::length(temp_grad))
  } else {theta = base::as.numeric(init_theta)}

  # Gradient ascent iteration.
  f0 = scDEDS::f_train_cal(theta)
  for (i in 1:iterations) {
    grad = scDEDS::f_train_cal_gr(branch_params_vec = theta,
                                  h_max = h_max,
                                  h_min_percent = h_min_percent,
                                  setseed = setseed,
                                  ncores = ncores)
    grad_need_cal = grad != 0

    # Selecting the optimal learning rate.
    best_alpha = stats::optim(
      par = alpha_guess,
      fn = function(al) {
        theta_new = theta
        theta_new[grad_need_cal] = base::pmin(
          base::pmax(theta[grad_need_cal] + al * grad[grad_need_cal], branch_params_vec_lower[grad_need_cal]),
          branch_params_vec_upper[grad_need_cal])
        -scDEDS::f_train_cal(theta_new)
      },
      lower = alpha_lower,
      upper = alpha_upper,
      method = "Brent"
    )

    # Parameter iteration.
    theta[grad_need_cal] = base::pmin(
      base::pmax(theta[grad_need_cal] + best_alpha$par * grad[grad_need_cal], branch_params_vec_lower[grad_need_cal]),
      branch_params_vec_upper[grad_need_cal]
    )

    # Printing current optimal value.
    current_f = -best_alpha$value
    if(f0 >= current_f) {
      message("Grad | iter = ", i, " | Best = ", round(f0, 6), "| learning_rate = ", round(best_alpha$par, 6), " | ncores = ", ncores)
      return(base::list(parameters = theta, objective_value = f0))
    } else {
      message("Grad | iter = ", i, " | Best = ", round(current_f, 6), "| learning_rate = ", round(best_alpha$par, 6), " | ncores = ", ncores)
      f0 = current_f
    }
  }

  # Returning result.
  base::list(parameters = theta, objective_value = f0)
}
