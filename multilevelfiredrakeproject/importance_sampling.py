"""

Importance Sampling

Alastair Gregory 2016

"""


def UpdateWeights(Ensemble_Hierarchy, ProposalDistribution):
    """ This updates the weights of the :class:`EnsembleHierarchy' by using an importance sampling update with the proposal distribution (this is a function taking a value and giving a probability) given

    """

    # Weights
    Weights = Ensemble_Hierarchy.Weights
    # Now use a weighting inverse to proposal distribution
    NewWeights = []
    for i in range(len(Weights)):
        for j in range(len(Weights[i])):
            for k in range(len(Weights[i][j])):
                NewWeights.append([])


class ImportanceSampling(object):

    def __init__(
            self,
            EnsembleStructure,
            OriginalWeights,
            ProposalDistribution):

        # Proposal Distribution must take Ensemble Structure In as arg. Whether
        # this is prior or posterior
        self.ProposalDistribution = ProposalDistribution
        # Array of all samples sizes in each ensemble
        self.N_Hierarchy = np.zeros(len(OriginalWeights))
        for i in range(len(OriginalWeights)):
            self.N_Hierarchy[i] = len(OriginalWeights[i][0])
        #
        self.L = len(self.N_Hierarchy)

        super(ImportanceSampling, self).__init__()

    def UpdateWeights(self):
        self.Weights = []
        for i in range(self.L):
            self.Weights.append([self.ProposalDistribution(EnsembleStructure[i][0]), self.ProposalDistribution(EnsembleStructure[
                                i][1])])  # This distribution might take all dimensions, or just state space dimensions into account
            # Rad-Nikson Deriv.
            self.Weights[i][0] = np.multiply(
                OriginalWeights[i][0], self.Weights[i][0])
            self.Weights[i][1] = np.multiply(
                OriginalWeights[i][1], self.Weights[i][1])
            # Normalize
            self.Weights[i][0] = (
                self.Weights[i][0] /
                np.sum(
                    self.Weights[i][0]))
            self.Weights[i][1] = (
                self.Weights[i][1] /
                np.sum(
                    self.Weights[i][1]))
        return self.Weights
    #### USE THE BELOW FUNCTIONALITY TO REPLACE THE CURRENT ONLINE MEAN IN ERR

    def WeightedOnlineMean(self, PosteriorEnsemble):
        El = []
        for i in range(self.L):
            N = np.shape(PosteriorEnsemble[i][0])[-1]
            if i == 0:
                for j in range(N):
                    if j == 0:
                        El.append(PosteriorEnsemble[i][
                                  1][..., j] * self.Weights[i][1][j])
                    else:
                        El[i] += PosteriorEnsemble[i][1][..., j] * \
                            self.Weights[i][1][j]
            else:
                for j in range(N):
                    if j == 0:
                        El.append((PosteriorEnsemble[i][1][..., j] * self.Weights[i][1][j]) - (
                            PosteriorEnsemble[i][0][..., j] * self.Weights[i][0][j]))
                    else:
                        El[i] += (PosteriorEnsemble[i][1][..., j] * self.Weights[i][1][j]) - \
                            (PosteriorEnsemble[i][0][..., j] * self.Weights[i][0][j])
        return
