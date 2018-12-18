import numpy as np
from scipy.optimize import fmin_l_bfgs_b

class model(object):
    """Define a linear model with repeated observations.

    Parameters
    ----------
    y : :class:`~numpy:numpy.array`
        1D array of phenotype observations
    X : :class:`~numpy:numpy.array`
        Design matrix for the fixed mean effects.
    labels : :class:`~numpy:numpy.array`
        1D array of sample labels

    Returns
    -------
    model : :class:`regrnd.model`

    """
    def __init__(self,y,X,labels):
        # Get sample size
        self.n = X.shape[0]
        self.X=X
        # Label mapping
        self.label_counts = dict()
        self.label_indices = dict()
        for l in xrange(0,labels.shape[0]):
            if labels[l] not in self.label_counts:
                self.label_counts[labels[l]]=1
                self.label_indices[labels[l]] = [l]
            else:
                self.label_counts[labels[l]]+=1
                self.label_indices[labels[l]].append(l)
        self.y_lab = dict()
        self.X_lab = dict()
        for label in self.label_indices.iterkeys():
            self.y_lab[label]=y[self.label_indices[label]]
            self.X_lab[label]=X[self.label_indices[label],:]
        self.n_labels = len(self.y_lab.keys())
        # response
        self.y=y
        self.labels=labels


    # Compute likelihood of data given beta, alpha
    def likelihood_and_gradient(self, alpha, sigma2, tau, l):
        """
        Compute the loss function, which is -2 times the likelihood plus a L2 regularisation term,
        along with its gradient

        Parameters
        ----------
        alpha : :class:`~numpy:numpy.array`
            value of the regression coefficients
        sigma2 : :class:`float`
            variance of model residuals
        tau : :class:`float`
            ratio of variance of model residuals to variance explained by mean differences between individuals

        Returns
        -------
        L, grad : :class:`float`
            loss function and gradient, divided by sample size

        """
        l = np.array(l)
        if l.shape == () or len(l) == self.X.shape[1]:
            pass
        else:
            raise(ValueError('Incorrect length of regularisation vector'))
        alpha_T_scaled = alpha.T * l
        ## Likelihood
        resid = self.y - self.X.dot(alpha)
        RSS = np.sum(np.square(resid))

        L = self.n * np.log(sigma2)+RSS/sigma2+alpha_T_scaled.dot(alpha)

        ## Gradient with respect to alpha
        grad_alpha = -2 * resid.T.dot(self.X)/sigma2+2*alpha_T_scaled

        ## Gradient with respect to sigma2
        grad_sigma2 = self.n/sigma2-RSS/np.square(sigma2)

        ## Gradient with respect to tau
        grad_tau = 0

        for label in self.y_lab.iterkeys():
            resid_label=resid[self.label_indices[label]]
            resid_sum = np.sum(resid_label)
            resid_square_sum = np.square(resid_sum)
            # Add to likelihood
            L = L - resid_square_sum/(sigma2*(tau+self.label_counts[label]))+np.log(1+self.label_counts[label]/tau)
            # Add to alpha gradient
            grad_alpha = grad_alpha + (2/sigma2)*resid_sum*np.sum(self.X_lab[label],axis=0)/(tau+self.label_counts[label])
            # Add to grad sigma2
            grad_sigma2+=resid_square_sum/(np.square(sigma2)*(tau+self.label_counts[label]))
            # Add to grad tau
            grad_tau+=(resid_square_sum/sigma2-self.label_counts[label]*(1+self.label_counts[label]/tau))/np.square(tau+self.label_counts[label])

        # Overall gradient vector
        grad = np.hstack((grad_alpha,grad_sigma2,grad_tau))

        return L/self.n, grad/self.n

    def optimize_model(self,l,init_params = None):
        """
        Find the parameters that minimise the loss function for a given regularisation parameter

        Parameters
        ----------
        l : :class:`array`
            array of regularisation parameters for regression coefficients

        Returns
        -------
        optim : :class:`list`
            the output of the scipy.fmin_l_bfgs_b function, first element has optimised parameters
        """
        # Initialise parameters
        if init_params is None:
            init_params=np.zeros((self.X.shape[1]+2))
        # Optimize
        optimized = fmin_l_bfgs_b(func=lik_and_grad,x0=init_params,
                                args=(self.y, self.X, self.labels, l))
        self.alpha = optimized[0][0:self.X.shape[1]]

        return optimized

    def predict(self,X):
        """
        Predict new observations based on model regression coefficients

        Parameters
        ----------
        X : :class:`array`
            matrix of feature observations to predict from

        Returns
        -------
        y : :class:`array`
            predicted values
        """
        if hasattr(self,'alpha'):
            return X.dot(self.alpha)
        else:
            raise(AttributeError('Model does not have known regression coefficients. Try optimizing model first'))

    def set_alpha(self,alpha):
        self.alpha = alpha

def lik_and_grad(pars,*args):
    # Wrapper for function to pass to L-BFGS-B
    y, X, labels, l = args
    mod = model(y,X,labels)
    return mod.likelihood_and_gradient(pars[0:X.shape[1]],np.exp(pars[X.shape[1]]),np.exp(pars[X.shape[1]+1]), l)

def simulate(n,alpha,sigma2,tau):
    """Simulate from a linear model with repeated observations from individuals. The mean for each individual
     is drawn from a normal distribution.

    Parameters
    ----------

    n : :class:`int`
        sample size
    alpha : :class:`~numpy:numpy.array`
        value of regression coefficeints
    sigma2 : :class:`float`
        variance of residuals
    tau : :class:`float`
        ratio of variance of residuals to variance of distribution of between individual means

    Returns
    -------
    model : :class:`regrnd.model`
        linear model with repeated observations
    """
    c = alpha.shape[0]
    X = np.random.randn((n * c)).reshape((n, c))
    labels = np.random.choice(n,n)
    random_effects = np.sqrt(sigma2/tau)*np.random.randn(n)
    y = X.dot(alpha)+random_effects[labels-1]+np.random.randn(n)*np.sqrt(sigma2)
    return model(y,X,labels)