#' @title This function is the minimized objective function in the model trained using the training set data.
#'
#' @param branch_params_vec Vector. A flattened vector containing all model parameters for a specific branch, with parameter names defined before.
#'
#' @returns
#' A numeric value representing the negative weighted objective function (multiplied by -100 for optimization purposes).
#' Lower values indicate better model fit. The objective function combines:
#' \itemize{
#' \item Mean squared prediction errors for TFE, TGA, and TGE values across time points
#' \item Variance regularization terms for beta parameters and expression states
#' \item Weighting factors that emphasize different components of the model
#' }
#'
#' @details
#' The function performs the following complex operations for each TG-TF pair in the training data:
#' \enumerate{
#' \item Parameter Organization: Assigns names to the parameter vector
#' \item TG-TF Pair Processing: Extracts and processes each regulator-target pair from training data
#' \item Objective Computation: Calculates weighted squared errors between predicted and observed values The mathematical formulation incorporates, where the error terms include squared differences for TFE, TGA, and TGE predictions, and regularization terms include variances of beta parameters and expression states.
#' }
#'
#' @note
#' Important considerations:
#' \itemize{
#' \item This function assumes the existence of several global variables:
#' \code{branch_params_names}, \code{train_data_pos}, \code{tao}, and \code{Tslot_K}
#' \item The function returns negative values for compatibility with minimization algorithms
#' \item Parameter names must follow the specific convention defined in \code{branch_params_names}
#' \item The function is computationally intensive due to multiple nested operations
#' \item Regularization terms help prevent overfitting by penalizing parameter variations
#' \item The weighting scheme (100 times and -100 times) is designed for optimization stability
#' }
#'
#' @seealso \code{\link{R_cal}}, \code{\link{Hill_cal}}, \code{\link{S_cal}}
#'
#' @export
#'
#' @examples
#' \dontrun{
#' f0 = f_train_cal(theta)
#' }
f_train_cal = function(branch_params_vec = branch_params_vec)
{
  base::names(branch_params_vec) = branch_params_names
  -100 * base::sum(base::sapply(base::rownames(train_data_pos), function(TG_TF) {
    TG = base::sub(".*_([^~]+)~.*", "\\1", TG_TF)
    TF = base::sub(".*~([^_]+)_.*", "\\1", TG_TF)

    TFE = train_data_pos[TG_TF, grep("TFE_", base::colnames(train_data_pos), value = TRUE)]
    TGA = train_data_pos[TG_TF, grep("TGA_", base::colnames(train_data_pos), value = TRUE)]
    TGE = train_data_pos[TG_TF, grep("TGE_", base::colnames(train_data_pos), value = TRUE)]
    Length = base::length(TFE)
    seqLength = seq_len(Length - 1)

    R = scDEDS::R_cal(
      theta_TF_TG = train_data_pos[TG_TF, "theta_s"],
      tao = tao,
      r1 = 200,
      r2 = (branch_params_vec[base::paste0("r.r2_TG.",TG)] +
              branch_params_vec[base::paste0("r.r2_TF.",TF)]) / 2,
      r3 = 1,
      r4 = 1000,
      r5 = 1
    )

    alpha1 = (branch_params_vec[base::paste0("alpha.alpha1_TG_K-1.", TG, ".", seqLength)] +
                branch_params_vec[base::paste0("alpha.alpha1_TF_K-1.", TF, ".", seqLength)]) / 2
    alpha2 = (branch_params_vec[base::paste0("alpha.alpha2_TG_K-1.", TG, ".", seqLength)] +
                branch_params_vec[base::paste0("alpha.alpha2_TF_K-1.", TF, ".", seqLength)]) / 2

    beta1 = (branch_params_vec[base::paste0("beta.beta1_TG_K-1.", TG, ".", seqLength)] +
               branch_params_vec[base::paste0("beta.beta1_TF_K-1.", TF, ".", seqLength)]) / 2
    beta2 = (branch_params_vec[base::paste0("beta.beta2_TG_K-1.", TG, ".", seqLength)] +
               branch_params_vec[base::paste0("beta.beta2_TF_K-1.", TF, ".", seqLength)]) / 2
    beta3 = branch_params_vec[base::paste0("beta.beta3_TG_K-1.", TG, ".", seqLength)]

    s1 = (branch_params_vec[base::paste0("s.s1_TG.", TG)] +
            branch_params_vec[base::paste0("s.s1_TF.", TF)]) / 2
    s2 = (branch_params_vec[base::paste0("s.s2_TG.", TG)] +
            branch_params_vec[base::paste0("s.s2_TF.", TF)]) / 2

    u11 = branch_params_vec[base::paste0("u.u11_TG.", TG)]
    u21 = branch_params_vec[base::paste0("u.u21_TG.", TG)]
    u311 = branch_params_vec[base::paste0("u.u311_TG.", TG)]
    u321 = branch_params_vec[base::paste0("u.u321_TG.", TG)]
    u331 = branch_params_vec[base::paste0("u.u331_TG.", TG)]
    u12 = branch_params_vec[base::paste0("u.u12_TG.", TG)]
    u22 = branch_params_vec[base::paste0("u.u22_TG.", TG)]
    u312 = branch_params_vec[base::paste0("u.u312_TG.", TG)]
    u322 = branch_params_vec[base::paste0("u.u322_TG.", TG)]
    u332 = branch_params_vec[base::paste0("u.u332_TG.", TG)]

    TGE_tilde = branch_params_vec[base::paste0("E~.E~_TG_K-1.", TG, ".", seqLength)]
    TFE_tilde = branch_params_vec[base::paste0("E~.E~_TF_K-1.", TF, ".", seqLength)]

    U1_tilde = scDEDS::Hill_cal(x = TGE_tilde, Dissociation_Constant = u11, Hill_Coefficient = u12)
    U2_tilde = scDEDS::Hill_cal(x = TFE_tilde, Dissociation_Constant = u21, Hill_Coefficient = u22)
    U1 = scDEDS::Hill_cal(x = TGE[-Length], Dissociation_Constant = u11, Hill_Coefficient = u12)
    U2 = scDEDS::Hill_cal(x = TFE[-Length], Dissociation_Constant = u21, Hill_Coefficient = u22)
    U31 = scDEDS::Hill_cal(x = TGA[-Length], Dissociation_Constant = u311, Hill_Coefficient = u312)
    U32 = scDEDS::Hill_cal(x = TGE[-Length], Dissociation_Constant = u321, Hill_Coefficient = u322)
    U33 = scDEDS::Hill_cal(x = TGE_tilde[-Length], Dissociation_Constant = u331, Hill_Coefficient = u332)

    v1 = branch_params_vec[base::paste0("v.v1_TG.", TG)]
    v2 = branch_params_vec[base::paste0("v.v2_TG.", TG)]
    v3 = branch_params_vec[base::paste0("v.v3_TG.", TG)]

    # Computing the objective function for this sample.
    100 * base::mean(base::as.numeric(base::sapply(2:Length, function(K) {
      prev_idx = K - 1
      2 * (TFE[prev_idx] + alpha1[prev_idx] * R * S_cal(s = s1, U = U1[prev_idx], U_tilde = U1_tilde[prev_idx]) + beta1[prev_idx] - TFE[K])^2 +
        2 * (TGA[prev_idx] + alpha2[prev_idx] * R * S_cal(s = s2, U = U2[prev_idx], U_tilde = U2_tilde[prev_idx]) + beta2[prev_idx] - TGA[K])^2 +
        (TGE[prev_idx] + R * (v1 * U31[prev_idx] - v2 * U32[prev_idx] - v3 * U33[prev_idx]) * Tslot_K[K] + beta3[prev_idx] - TGE[K])^2
    }))) +
      stats::var(beta1) + stats::var(beta2) + stats::var(beta3) + stats::var(TGE_tilde) + stats::var(TFE_tilde)
  }))
}
