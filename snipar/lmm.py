import numpy as np
from scipy.optimize import fmin_l_bfgs_b

class model(object):
    """Define a linear model with within-class correlations.

    Args:
        y : :class:`~numpy:numpy.array`
            1D array of phenotype observations
        X : :class:`~numpy:numpy.array`
            Design matrix for the fixed mean effects.
        labels : :class:`~numpy:numpy.array`
            1D array of sample labels

    Returns:
        model : :class:`snipar.model`

    """

    def __init__(self, y, X, labels, add_intercept=False):
        if y.shape[0] == X.shape[0] and X.shape[0] == labels.shape[0]:
            pass
        else:
            raise (ValueError('inconsistent sample sizes of response, covariates, and labels'))
        # Get sample size
        self.n = X.shape[0]
        if X.ndim == 1:
            X = X.reshape((self.n, 1))
        if add_intercept:
            X = np.hstack((np.ones((self.n, 1), dtype=X.dtype), X))
        self.X = X
        # Label mapping
        self.label_counts = dict()
        self.label_indices = dict()
        for l in range(0, labels.shape[0]):
            if labels[l] not in self.label_counts:
                self.label_counts[labels[l]] = 1
                self.label_indices[labels[l]] = [l]
            else:
                self.label_counts[labels[l]] += 1
                self.label_indices[labels[l]].append(l)
        self.y_lab = dict()
        self.X_lab = dict()
        for label in self.label_indices.keys():
            self.y_lab[label] = y[self.label_indices[label]]
            self.X_lab[label] = X[self.label_indices[label], :]
        self.n_labels = len(self.y_lab.keys())
        # response
        self.y = y
        self.labels = labels

    def alpha_mle(self, tau, sigma2, compute_cov=False, xtx_out=False):
        """
        Compute the MLE of alpha given variance parameters

        Args:
            sigma2 : :class:`float`
                variance of model residuals
            tau : :class:`float`
                ratio of variance of model residuals to variance explained by mean differences between classes

        Returns:
            alpha : :class:`~numpy:numpy.array`
                MLE of alpha

        """
        X_T_X = np.zeros((self.X.shape[1], self.X.shape[1]), dtype=np.float64)
        X_T_y = np.zeros((self.X.shape[1]), dtype=np.float64)

        for label in self.y_lab.keys():
            sigma_u = sigma2 / tau
            Sigma_lab = sigma_u * np.ones((self.label_counts[label], self.label_counts[label]))
            np.fill_diagonal(Sigma_lab, sigma_u + sigma2)
            Sigma_lab_inv = np.linalg.inv(Sigma_lab)
            X_T_X = X_T_X + np.dot(self.X_lab[label].T, Sigma_lab_inv.dot(self.X_lab[label]))
            X_T_y = X_T_y + np.dot(self.X_lab[label].T, Sigma_lab_inv.dot(self.y_lab[label]))

        if xtx_out:
            return [X_T_X, X_T_y.reshape((self.X.shape[1]))]
        else:
            alpha = np.linalg.solve(X_T_X, X_T_y)
            alpha = alpha.reshape((alpha.shape[0],))

            if compute_cov:
                alpha_cov = np.linalg.inv(X_T_X)
                return [alpha, alpha_cov]
            else:
                return alpha

    # Compute likelihood of data given beta, alpha
    def likelihood_and_gradient(self, sigma2, tau):
        """
        Compute the loss function, which is -2 times the likelihood along with its gradient

        Args:
            sigma2 : :class:`float`
                variance of model residuals
            tau : :class:`float`
                ratio of variance of model residuals to variance explained by mean differences between classes

        Returns:
            L, grad : :class:`float`
                loss function and gradient, divided by sample size

        """
        ## Likelihood
        alpha = self.alpha_mle(tau, sigma2)
        resid = self.y - self.X.dot(alpha)
        RSS = np.sum(np.square(resid))

        L = self.n * np.log(sigma2) + RSS / sigma2

        ## Gradient with respect to sigma2
        grad_sigma2 = self.n / sigma2 - RSS / np.square(sigma2)

        ## Gradient with respect to tau
        grad_tau = 0

        for label in self.y_lab.keys():
            resid_label = resid[self.label_indices[label]]
            resid_sum = np.sum(resid_label)
            resid_square_sum = np.square(resid_sum)
            # Add to likelihood
            L = L - resid_square_sum / (sigma2 * (tau + self.label_counts[label])) + np.log(
                1 + self.label_counts[label] / tau)
            # Add to grad sigma2
            grad_sigma2 += resid_square_sum / (np.square(sigma2) * (tau + self.label_counts[label]))
            # Add to grad tau
            grad_tau += (resid_square_sum / sigma2 - self.label_counts[label] * (
                        1 + self.label_counts[label] / tau)) / np.square(tau + self.label_counts[label])

        # Overall gradient vector
        grad = np.hstack((grad_sigma2, grad_tau))

        return L / self.n, grad / self.n

    def optimize_model(self, init_params):
        """
        Find the parameters that minimise the loss function for a given regularisation parameter

        Args:
            init_param : :class:`array`
                initial values for residual variance (sigma^2_epsilon) followed by ratio
                of residual variance to within-class variance (tau)

        Returns:
            optim : :class:`dict`
                dictionary with keys: 'success', whether optimisation was successful (bool);
                'warnflag', output of L-BFGS-B algorithm giving warnings; 'sigma2', MLE of
                residual variance; 'tau', MLE of ratio of residual variance to within-class variance;
                'likelihood', maximum of likelihood.

        """
        # Paramtere boundaries
        parbounds = [(0.00001, None), (0.00001, None)]
        # Optimize
        optimized = fmin_l_bfgs_b(func=lik_and_grad, x0=init_params,
                                  args=(self.y, self.X, self.labels),
                                  bounds=parbounds)

        # Get MLE
        optim = {}
        optim['success'] = True
        optim['warnflag'] = optimized[2]['warnflag']
        if optim['warnflag'] != 0:
            print('Optimization unsuccessful.')
            optim['success'] = False
        optim['sigma2'] = optimized[0][0]
        optim['tau'] = optimized[0][1]
        # Get parameter covariance
        optim['likelihood'] = -0.5 * np.float64(self.n) * (optimized[1] + np.log(2 * np.pi))

        return optim

    def sigma_inv_root(self, tau, sigma2):
        sigma_u = sigma2 / tau
        sigma2_nsqrt = dict()
        famsizes = np.unique(list(self.label_counts.values()))
        sigma2_nsqrt[1] = np.power(sigma_u + sigma2, -0.5)
        famsizes = famsizes[famsizes > 1]
        for famsize in famsizes:
            Sigma_lab = sigma_u * np.ones((famsize, famsize))
            np.fill_diagonal(Sigma_lab, sigma_u + sigma2)
            vals, vectors = np.linalg.eigh(Sigma_lab)
            vals = np.power(vals, 0.25)
            vectors = vectors / vals
            sigma2_nsqrt[famsize] = vectors.dot(vectors.T)
        return sigma2_nsqrt

    def predict(self, X):
        """
        Predict new observations based on model regression coefficients

        Args:
            X : :class:`array`
                matrix of covariates to predict from

        Returns:
            y : :class:`array`
                predicted values

        """
        if hasattr(self, 'alpha'):
            return X.dot(self.alpha)
        else:
            raise (AttributeError('Model does not have known regression coefficients. Try optimizing model first'))

    def set_alpha(self, alpha):
        self.alpha = alpha


def lik_and_grad(pars, *args):
    # Wrapper for function to pass to L-BFGS-B
    y, X, labels = args
    mod = model(y, X, labels)
    return mod.likelihood_and_gradient(pars[0], pars[1])

def simulate(n, alpha, sigma2, tau):
    """Simulate from a linear model with correlated observations within-class. The mean for each class
     is drawn from a normal distribution.

    Args:
        n : :class:`int`
            sample size
        alpha : :class:`~numpy:numpy.array`
            value of regression coefficeints
        sigma2 : :class:`float`
            variance of residuals
        tau : :class:`float`
            ratio of variance of residuals to variance of distribution of between individual means

    Returns:
        model : :class:`regrnd.model`
            linear model with repeated observations

    """
    c = alpha.shape[0]
    # X = np.random.randn((n * c)).reshape((n, c))
    X_cov = np.ones((c, c))
    np.fill_diagonal(X_cov, 1.2)
    X = np.random.multivariate_normal(np.zeros((c)), X_cov, n).reshape((n, c))
    labels = np.random.choice(n // 10, n)
    random_effects = np.sqrt(sigma2 // tau) * np.random.randn(n)
    y = X.dot(alpha) + random_effects[labels - 1] + np.random.randn(n) * np.sqrt(sigma2)
    return model(y, X, labels)

def fit_model(y, X, fam_labels, add_intercept=False, tau_init=1, return_model=True, return_vcomps=True, return_fixed=True):
    """Compute the MLE for the fixed effects in a family-based linear mixed model.

    Args:
        y : :class:`~numpy:numpy.array`
            vector of phenotype values
        X: :class:`~numpy:numpy.array`
            regression design matrix for fixed effects
        fam_labels : :class:`~numpy:numpy.array`
            vector of family labels: residual correlations in y are modelled between family members (that share a fam_label)
        add_intercept : :class:`bool`
            whether to add an intercept to the fixed effect design matrix

    Returns:
       model : :class:`snipar.model`
            the snipar model object, if return_model=True
       vcomps: :class:`float`
            the MLEs for the variance parameters: sigma2 (residual variance) and tau (ratio between sigma2 and family variance), if return_vcomps=True
       alpha : :class:`~numpy:numpy.array`
            MLE of fixed effects, if return_fixed=True
       alpha_cov : :class:`~numpy:numpy.array`
            sampling variance-covariance matrix for MLE of fixed effects, if return_fixed=True

    """
    # Optimize  model
    sigma_2_init = np.var(y)*tau_init/(1+tau_init)
    null_model = model(y, X, fam_labels, add_intercept=add_intercept)
    null_optim = null_model.optimize_model(np.array([sigma_2_init,tau_init]))
    # Create return list
    return_list = []
    if return_model:
        return_list.append(null_model)
    if return_vcomps:
        return_list += [null_optim['sigma2'], null_optim['tau']]
    if return_fixed:
        return_list += null_model.alpha_mle(null_optim['tau'], null_optim['sigma2'], compute_cov=True)
    return return_list