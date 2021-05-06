from warnings import simplefilter
# ignore all future warnings
simplefilter(action='ignore', category=FutureWarning)
import pandas as pd
import numpy as np

import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import train_test_split
from sklearn.feature_selection import RFE
from sklearn.linear_model import LogisticRegression
from sklearn.decomposition import PCA
from sklearn.metrics import roc_auc_score
from sklearn.metrics import roc_curve

import statsmodels.api as sm
from imblearn.over_sampling import SMOTE
from sklearn import svm
import gseapy

###################################################
# This script uses logistic regression for the incidence of ORN
# (yes or no) and considers no grade or Dur_to_ORN in the analysis.
###################################################

def do_gene_enrichment(gene_df, dataname):
    '''
    This function consumes a data frame of gene names and values and uses the
    `gseapy` package to assess any function in the genes selected
    '''
    gl = gene_df.sort_values(by='value', ascending=False)['gene']
    res = gseapy.enrichr(gene_list=gl, description='pathway', gene_sets='KEGG_2016', \
                         outdir=dataname)
    print(res)
    return res

def prepareData(datapath_prot):
    # Reading the proteomic dataset
    ###################################################
    datapath_prot = "data/S054_HNSCC_imputed_0920.tsv"

    df_prot = pd.read_csv(datapath_prot, header=None, sep='\t')
    data_top = df_prot.head()
    #pIDs =df.index        # Patient IDs
    df_prot = df_prot.transpose()
    new_header = df_prot.iloc[0]
    df_prot = df_prot[1:]
    df_prot.columns = new_header
    df_prot = df_prot.rename(columns={df_prot.columns[0]: "pID"})
    df_prot.reindex

    # Cleaning the data and adding case_id
    ###################################################
    #--Removing meaningless pIDs
    pIDsToBeRemoved = ['Channel 128C', 'QC2', "QC3", "QC4",
                   "C3L-02617-T-duplicate", "Channel 129N",
                   "LungTumor1", "Pooled-sample14", "LungTumor2",
                   "QC6", "LungTumor3", "Pooled-sample17", "QC7",
                   "Pooled-sample19"]
    pIDsToBeRemoved_indices = [4, 34, 62, 81, 91, 115, 121, 133, 142, 147, 166,
                           170, 178, 187]

    df_prot.drop(index=pIDsToBeRemoved_indices, inplace = True)

    #--Removing spaces from pID
    pID = pd.Series(data = df_prot['pID'])
    pID = pID.str.strip()
    pID = pID.str.replace(' ', '', regex=False)
    df_prot.drop(columns=['pID'], inplace=True)
    df_prot['pID'] = pID.to_list()


    #--Adding case_id to df_prot: this means cutting -T or -N from pID
    case_id = [str[:-2].strip() for str in df_prot['pID'].tolist()]
    #i = 0
    #for str in case_id:
    #    str = str[:-2]
    #    case_id[i] = str.strip()
    #    i = i+1
    df_prot['case_id'] = case_id

    #df_prot.columns = df_prot.columns.str.replace(' ', '')

    #--Reindexing columns to have pID, case_id, ...
    cols = df_prot.columns
    cols =  cols[:-2].insert(0, 'case_id')
    cols =  cols.insert(0, 'pID')
    df_prot = df_prot.reindex(columns = cols)

    # Reading the target variables
    ###################################################
    datapath_target = "data/HNSCC_Feb2021_clinical_data.csv"
    df_target = pd.read_csv(datapath_target)
    target_cols = ['case_id','baseline/lymph_nodes_extranodal_extension',\
                 'baseline/lymph_vascular_invasion',\
                 'baseline/number_of_lymph_nodes_positive_for_tumor_by_he']

    df_target = df_target[target_cols]
    # Merging the datasets
    ###################################################
    df_main = pd.merge(df_prot, df_target, on=["case_id"])

    #--Reindexing columns to have pID, case_id,
    #      baseline_lymph_nodes_extranodal_extension,
    #      baseline_lymph_vascular_invasion,
    #      baseline_number_of_lymph_nodes_positive_for_tumor_by_he

   # cols = df_prot.columns
   # cols = cols[:-3].insert(2, 'baseline_lymph_nodes_extranodal_extension')
   # cols = cols.insert(3, 'baseline_lymph_vascular_invasion')
   # cols = cols.insert(4, 'baseline_number_of_lymph_nodes_positive_for_tumor_by_he')
   # df_main = df_main.reindex(columns = cols)

    #--Renaming the target columns
    df_main.rename(columns={"baseline/lymph_nodes_extranodal_extension": "ExtraNodalExtension",\
                            "baseline/lymph_vascular_invasion": "VasInv",\
                            "baseline/number_of_lymph_nodes_positive_for_tumor_by_he": "numLymphPositiveForTumor"},\
                   errors="raise", inplace = True)

    ###################################################
    # At this time, df is ready for ML. It contains three target variables,
    # we make sure that we work with only one of them.
    ###################################################
    return df_main

def buildBinClassifier(df_main,
                        ExtraNodalExtension_switch,
                        VasInv_switch,
                        useSmoth_switch=True,
                        test_size_ratio = 1/5.0,
                        PCA_switch=False,
                        LinReg_switch=True,
                        SVM_switch=False):
    if (ExtraNodalExtension_switch == VasInv_switch):
        print("Choose exactly one target variable.\n")
        return -1

    all_targs = ['pID','case_id','VasInv','ExtraNodalExtension','numLymphPositiveForTumor']

    if (ExtraNodalExtension_switch == True):
        #target_col_to_be_dropped = 'VasInv'
        target_col = 'ExtraNodalExtension'
    else:
        #target_col_to_be_dropped = 'ExtraNodalExtension'
        target_col ='VasInv'

    #drop everything BUT the target column
    df = df_main.drop(columns=[a for a in all_targs if a != target_col])

    # Next, we get rid of "Indeterminate" labels
    df = df[df[target_col] != 'Indeterminate']

    # Setting up labels
    df[target_col].replace('Not identified', '0', inplace=True)
    df[target_col].replace('Present', '1', inplace=True)
   # print(df[target_col])
    df = df.astype({target_col: 'int32'})

    # setting up y and X
    features = df.columns[3:]

    y = df.loc[:, df.columns == target_col]
    X = df.loc[:, df.columns != target_col]


    # Using SMOTE
    ###################################################
    if (useSmoth_switch == True):
        os  = SMOTE(random_state = 0)
        #X_train, y_train = os.fit_resample(X_train, y_train)
        X, y = os.fit_resample(X, y)


        # Train/test split
    ###################################################
    X_train, X_test, y_train, y_test = train_test_split(X, y,
                                                        test_size = test_size_ratio,
                                                        random_state=0)


    # Normalization of the numeric data
    ###################################################
    stdscaler = StandardScaler()
    stdscaler.fit(X_train)
    X_train = stdscaler.transform(X_train)
    X_test = stdscaler.transform(X_test)

    # Applying PCA
    ###################################################
    if (PCA_switch==True):
        varianceToBeCovered = 0.80
        pca = PCA(varianceToBeCovered)
        pca.fit(X_train)
        print("To cover {} of the variance, we need {} principle components."
              .format(varianceToBeCovered, pca.n_components_))

        X_train = pca.transform(X_train)
        X_test = pca.transform(X_test)


    # Training using logistic regression and creating performance report
    ###################################################
    if (LinReg_switch==True):
        logreg = LogisticRegression(max_iter = 1000)
        logreg.fit(X_train, y_train.values.ravel())

        # Performance study
        ###################################################
        y_pred = logreg.predict(X_test)
        print('Accuracy of logistic regression classifier on test set: {:.2f}'.format(logreg.score(X_test, y_test)))

        #Confusion Matrix
        from sklearn.metrics import confusion_matrix
        confusion_matrix = confusion_matrix(y_test, y_pred)
        print(confusion_matrix)

        #Compute precision, recall, F-measure and support
        from sklearn.metrics import classification_report
        print(classification_report(y_test, y_pred))

        from sklearn import metrics
        metrics.plot_roc_curve(logreg, X_test, y_test)
        plt.plot([0, 1], [0, 1],'r--')
        plt.xlabel('False Positive Rate')
        plt.ylabel('True Positive Rate')
        if (ExtraNodalExtension_switch==True):
            if (PCA_switch==True):
                titleStr = '(Target variable: ENE, PCA: Active)'
            else:
                titleStr = '(Target variable: ENE)'
        else:
            if (PCA_switch==True):
                titleStr = '(Target variable: VasInv)'
            else:
                titleStr = '(Target variable: VasInv, PCA: Active)'

        plt.title('Log Regression ROC '+titleStr)
        plt.legend(loc="lower right")
        plt.show()

        #now get gene sets
        genevals = pd.DataFrame({'value':logreg.coef_[0],'gene':X.columns})
        do_gene_enrichment(genevals,titleStr)

    if (SVM_switch==True):
        weights = {0: 1.00, 1: 1.00}
        if (useSmoth_switch==False):
            # Cost-sensitive SVM
            numOneLabel = np.sum(y_train, axis = 0)
            numZeroLabel = np.shape(y_train) - numOneLabel
            weights = {0: (numOneLabel + numZeroLabel)/numZeroLabel, 1: (numOneLabel + numZeroLabel)/numOneLabel}

        svm_model = svm.SVC(kernel='linear', degree=2, class_weight=weights)
        svm_model.fit(X_train, y_train.values.ravel())
        ##I'm not sure how to get the weights of the genes from this model

        y_pred = svm_model.predict(X_test)
        print('Accuracy of SVM on test set: {:.2f}'.format(svm_model.score(X_test, y_test)))

        #Confusion Matrix
        from sklearn.metrics import confusion_matrix
        confusion_matrix = confusion_matrix(y_test, y_pred)
        print(confusion_matrix)

        #Compute precision, recall, F-measure and support
        from sklearn.metrics import classification_report
        print(classification_report(y_test, y_pred))

        #ROC Curve
        from sklearn import metrics
        metrics.plot_roc_curve(svm_model, X_test, y_test)
        plt.plot([0, 1], [0, 1],'r--')
        plt.xlabel('False Positive Rate')
        plt.ylabel('True Positive Rate')
        if (ExtraNodalExtension_switch==True):
            if (PCA_switch==True):
                titleStr = '(Target variable: ENE, PCA: Active)'
            else:
                titleStr = '(Target variable: ENE)'

        else:
            if (PCA_switch==True):
                titleStr = '(Target variable: VasInv)'
            else:
                titleStr = '(Target variable: VasInv, PCA: Active)'

        plt.title('SVM ROC '+titleStr)
        plt.legend(loc="lower right")
        plt.show()





if __name__ == "__main__":
    df_main = prepareData('')
    buildBinClassifier(df_main,
                        ExtraNodalExtension_switch=True,
                        VasInv_switch= False,
                        useSmoth_switch=True,
                        test_size_ratio = 1/5.0,
                        PCA_switch=False,
                        LinReg_switch=True,
                        SVM_switch=True)
    buildBinClassifier(df_main,
                       ExtraNodalExtension_switch=True,
                       VasInv_switch= False,
                       useSmoth_switch=True,
                       test_size_ratio = 1/5.0,
                       PCA_switch=False,
                       LinReg_switch=True,
                       SVM_switch=True)
