from __future__ import print_function, division
import numpy as np
import pandas as pd
import sys
from collections import OrderedDict


def sample_normal_dis(loc, scale, repeat):
    # data = np.zeros((loc.size, repeat), dtype=np.float)
    data = list()
    for i in range(repeat):
        # data[:, i] = np.random.normal(loc=loc, scale=scale)
        data.append(np.random.normal(loc=loc, scale=scale))

    return(data)


if __name__ == "__main__":
    n_features = 5
    n_genes = 1000000
    n_PROseq_data = 2
    n_RNAseq_data = 2

    # step 1: generate feature matrix
    design_matrix = np.random.normal(size=(n_genes, n_features))
    one = np.ones((n_genes, 1))
    design_matrix = np.concatenate((one, design_matrix), axis=1)

    trans_weights = np.random.normal(size=n_features+1)
    halflife_weights = np.random.normal(size=n_features+1)

    # step 2: generate unobserved transcription rates
    # and halflives
    trans = (np.dot(design_matrix, trans_weights)
             + np.random.normal(scale=3., size=n_genes))
    halflife = (np.dot(design_matrix, halflife_weights)
                + np.random.normal(scale=2., size=n_genes))

    # step 3: sample measurements
    PROseq_data = sample_normal_dis(trans, 0.2, n_PROseq_data)
    RNAseq_data = sample_normal_dis(trans + halflife, 0.5, n_RNAseq_data)

    # put everything to a pandas data frame
    data_dict = OrderedDict()
    data_dict["gene"] = range(n_genes)
    for i in range(n_PROseq_data):
        data_dict["sample_" + str(i) + "_PROseq"] = PROseq_data[i]

    for i in range(n_RNAseq_data):
        data_dict["sample_" + str(i) + "_Exon"] = RNAseq_data[i]

    for i in range(1, design_matrix.shape[1]):
        data_dict["covariate_" + str(i - 1)] = design_matrix[:, i]

    final_data = pd.DataFrame.from_dict(data_dict)
    final_data.to_csv(sys.argv[1], index=False)

    # print out weights
    print("transcription weight", trans_weights[1:])
    print("halflife weight", halflife_weights[1:])
