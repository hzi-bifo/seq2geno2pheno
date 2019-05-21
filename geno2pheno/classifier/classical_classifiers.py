__author__ = "Ehsaneddin Asgari"
__license__ = "Apache 2"
__version__ = "1.0.0"
__maintainer__ = "Ehsaneddin Asgari"
__email__ = "asgari@berkeley.edu"
__project__ = "GENO2PHENO of SEQ2GENO2PHENO"
__website__ = ""

import sys

sys.path.append('../')
from sklearn.svm import LinearSVC, SVC
from classifier.cross_validation import KFoldCrossVal, PredefinedFoldCrossVal
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.neighbors import KNeighborsClassifier
from utility.file_utility import FileUtility
import numpy as np
import codecs
import math
import operator
from utility.list_set_util import argsort

class SVM:
    '''
        Support vector machine classifier
    '''

    def __init__(self, X, Y, clf_model='LSVM'):
        if clf_model == 'LSVM':
            self.model = LinearSVC(C=1.0)  # , multi_class='ovr'
            self.type = 'linear'
        else:
            self.model = SVC(C=1.0, kernel='rbf')
            self.type = 'rbf'
        self.X = X
        self.Y = Y

        SVM.params_tuning = [{'C': [1000, 500, 200, 100, 50, 20, 10, 5, 2, 1, 0.2, 0.5, 0.01, 0.02, 0.05, 0.001],
                               'penalty': ['l1'], "tol": [1e-06, 1e-04], 'dual': [False, True], "fit_intercept": [True],
                               'loss': ['l2'], 'class_weight': ['balanced', None]}]

    def tune_and_eval(self, results_file,
                      params=None, njobs=50, kfold=10, feature_names=None, optimized_for='f1_macro'):
        '''
        K-fold cross-validation
        :param results_file: file to save the results
        :param params: parameters to be tuned
        :param njobs: number of cores
        :param kfold: number of folds
        :return:
        '''
        if params==None:
            params=SVM.params_tuning
        CV = KFoldCrossVal(self.X, self.Y, folds=kfold)
        CV.tune_and_evaluate(self.model, parameters=params, score=optimized_for, file_name=results_file + '_SVM',
                             n_jobs=njobs)
        if feature_names is not None:
            [label_set, conf, label_set, best_score_, best_estimator_,
                              cv_results_, best_params_,  (cv_predictions_pred,cv_predictions_trues,isolates ), (Y_test_pred, Y_test) ] = FileUtility.load_obj(results_file + '_SVM.pickle')
            self.generate_SVM_important_features(best_estimator_, feature_names, results_file)

    def tune_and_eval_predefined(self, results_file, isolates, folds_file, test_file, params=None, njobs=50,feature_names=None, optimized_for='f1_macro'):
        '''
        :param results_file:
        :param isolates:
        :param folds:
        :param params:
        :param njobs:
        :return:
        '''
        if params==None:
            params=SVM.params_tuning
        self.CV = PredefinedFoldCrossVal(self.X, self.Y, isolates, folds_file, test_file)
        self.CV.tune_and_evaluate(self.model, parameters=params, score=optimized_for, file_name=results_file + '_SVM',
                                  n_jobs=njobs)
        if feature_names is not None:
            [label_set, conf, label_set, best_score_, best_estimator_,
                              cv_results_, best_params_,  (cv_predictions_pred,cv_predictions_trues ,isolates), (Y_test_pred, Y_test) ] = FileUtility.load_obj(results_file + '_SVM.pickle')
            self.generate_SVM_important_features(best_estimator_, feature_names, results_file)

    def generate_SVM_important_features(self, clf_SVM, feature_names, results_file, N=1000):
        '''
        :param clf_SVM:
        :param feature_names:
        :param results_file:
        :param N:
        :return:
        '''

        results_file=results_file.replace('/classifications/','/feature_selection/classifications/')
        FileUtility.ensure_dir(results_file)
        file_name = results_file + '_SVM'

        idxs=argsort(np.abs(clf_SVM.coef_.tolist()[0]).tolist(),rev=True)[0:N]

        f = codecs.open(file_name, 'w')
        f.write('\t'.join(['feature', 'score']) + '\n')
        for idx in idxs:
            f.write('\t'.join([feature_names[idx], str(clf_SVM.coef_.tolist()[0][idx])]) + '\n')
        f.close()



class LogRegression:
    '''
        LR classifier
    '''

    def __init__(self, X, Y):
        '''

        :param X:
        :param Y:
        '''
        self.model = LogisticRegression(C=1.0)
        self.X = X
        self.Y = Y
        LogRegression.params_tuning=[{'C': [1000, 500, 200, 100, 50, 20, 10, 5, 2, 1, 0.2, 0.5, 0.01, 0.02, 0.05, 0.001],
                               'penalty': ['l1'], "tol": [1e-06, 1e-04], 'dual': [False, True], "fit_intercept": [True],
                               'class_weight': ['balanced', None], 'solver': ['liblinear']}]

    def tune_and_eval(self, results_file,
                      params=None, njobs=50, kfold=10, feature_names=None, optimized_for='f1_macro'):
        '''
        :param results_file:
        :param params:
        :param njobs:
        :param kfold:
        :return:
        '''
        if params==None:
            params=LogRegression.params_tuning
        CV = KFoldCrossVal(self.X, self.Y, folds=kfold)
        CV.tune_and_evaluate(self.model, parameters=params, score=optimized_for, file_name=results_file + '_LR',
                             n_jobs=njobs)
        if feature_names is not None:
            [label_set, conf, label_set, best_score_, best_estimator_,
                              cv_results_, best_params_,  (cv_predictions_pred,cv_predictions_trues ,isolates), (Y_test_pred, Y_test) ] = FileUtility.load_obj(results_file + '_LR.pickle')
            self.generate_LR_important_features(best_estimator_, feature_names, results_file)

    def tune_and_eval_predefined(self, results_file, isolates, folds_file, test_file, params=None, njobs=50, feature_names=None, optimized_for='f1_macro'):
        '''
        :param results_file:
        :param isolates:
        :param folds:
        :param params:
        :param njobs:
        :return:
        '''
        if params==None:
            params=LogRegression.params_tuning
        self.CV = PredefinedFoldCrossVal(self.X, self.Y, isolates, folds_file, test_file)
        self.CV.tune_and_evaluate(self.model, parameters=params, score=optimized_for, file_name=results_file + '_LR',
                                  n_jobs=njobs)

        if feature_names is not None:
            [label_set, conf, label_set, best_score_, best_estimator_,
                              cv_results_, best_params_,  (cv_predictions_pred,cv_predictions_trues,isolates ), (Y_test_pred, Y_test) ] = FileUtility.load_obj(results_file + '_LR.pickle')
            self.generate_LR_important_features(best_estimator_, feature_names, results_file)

    def generate_LR_important_features(self, clf_LR, feature_names, results_file, N=1000):
        '''
        :param clf_logistic_regression:
        :param feature_names:
        :param results_file:
        :param N:
        :return:
        '''

        results_file=results_file.replace('/classifications/','/feature_selection/classifications/')
        FileUtility.ensure_dir(results_file)
        file_name = results_file + '_LR'

        idxs=argsort(np.abs(clf_LR.coef_.tolist()[0]).tolist(),rev=True)[0:N]

        f = codecs.open(file_name, 'w')
        f.write('\t'.join(['feature', 'score']) + '\n')
        for idx in idxs:
            f.write('\t'.join([feature_names[idx], str(clf_LR.coef_.tolist()[0][idx])]) + '\n')
        f.close()


class RFClassifier:
    '''
        Random forest classifier
    '''

    def __init__(self, X, Y):
        '''
        :param X:
        :param Y:
        '''
        self.model = RandomForestClassifier(bootstrap=True, criterion='gini',
                                            min_samples_split=2, max_features='auto', min_samples_leaf=1,
                                            n_estimators=1000)
        self.X = X
        self.Y = Y
        RFClassifier.params_tuning = [{"n_estimators": [100, 200, 500, 1000],
                   "criterion": ["entropy"],  # "gini",
                   'max_features': ['auto'],  # 'auto',
                   'min_samples_split': [5],  # 2,5,10
                   'min_samples_leaf': [1]}]  # 'class_weight': ['balanced', None]}]

    def tune_and_eval(self, results_file, params=None, feature_names=None, njobs=50, kfold=10, optimized_for='f1_macro'):
        '''
        :param results_file:
        :param params:
        :param feature_names:
        :param njobs:
        :param kfold:
        :return:
        '''
        if params is None:
            params=RFClassifier.params_tuning
        self.CV = KFoldCrossVal(self.X, self.Y, folds=kfold)
        self.CV.tune_and_evaluate(self.model, parameters=params, score=optimized_for, file_name=results_file + '_RF',
                                  n_jobs=njobs)
        if feature_names is not None:
            [label_set, conf, label_set, best_score_, best_estimator_,
                              cv_results_, best_params_,  (cv_predictions_pred,cv_predictions_trues,isolates ), (Y_test_pred, Y_test) ] = FileUtility.load_obj(results_file + '_RF.pickle')
            self.generate_RF_important_features(best_estimator_, feature_names, results_file)

    def tune_and_eval_predefined(self, results_file, isolates, folds_file, test_file, params=None, feature_names=None, njobs=50, optimized_for='f1_macro'):
        '''
        :param results_file:
        :param isolates:
        :param folds:
        :param params:
        :param feature_names:
        :param njobs:
        :return:
        '''
        if params is None:
            params=RFClassifier.params_tuning
        self.CV = PredefinedFoldCrossVal(self.X, self.Y, isolates, folds_file, test_file)
        self.CV.tune_and_evaluate(self.model, parameters=params, score=optimized_for, file_name=results_file + '_RF',
                                  n_jobs=njobs)
        if feature_names is not None:
            [label_set, conf, label_set, best_score_, best_estimator_,
                              cv_results_, best_params_,  (cv_predictions_pred,cv_predictions_trues,isolates ), (Y_test_pred, Y_test) ]= FileUtility.load_obj(results_file + '_RF.pickle')
            self.generate_RF_important_features(best_estimator_, feature_names, results_file)

    def generate_RF_important_features(self, clf_random_forest, feature_names, results_file, N=1000):
        '''
        :param clf_random_forest:
        :param feature_names:
        :param results_file:
        :param N:
        :return:
        '''

        results_file=results_file.replace('/classifications/','/feature_selection/classifications/')
        FileUtility.ensure_dir(results_file)
        file_name = results_file + '_RF'
        clf_random_forest.fit(self.X, self.Y)
        std = np.std([tree.feature_importances_ for tree in clf_random_forest.estimators_], axis=0)

        scores = {feature_names[i]: (s, std[i]) for i, s in enumerate(list(clf_random_forest.feature_importances_)) if
                  not math.isnan(s)}
        scores = sorted(scores.items(), key=operator.itemgetter([1][0]), reverse=True)[0:N]
        f = codecs.open(file_name, 'w')
        f.write('\t'.join(['feature', 'score']) + '\n')
        for w, score in scores:
            #feature_array = self.X[:, feature_names.index(w)]
            #pos = [feature_array[idx] for idx, x in enumerate(self.Y) if x == 1]
            #neg = [feature_array[idx] for idx, x in enumerate(self.Y) if x == 0]
            f.write('\t'.join([str(w), str(score[0])]) + '\n')
        f.close()


class KNN:
    '''
        K-nearest neighbor classifier
    '''

    def __init__(self, X, Y):
        '''
        :param X:
        :param Y:
        '''
        self.model = KNeighborsClassifier(n_neighbors=3)
        self.X = X
        self.Y = Y
        KNN.parameter_tuning= [{"n_neighbors": [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 15, 20], 'weights': ['uniform', 'distance']}]

    def tune_and_eval(self, results_file, params=None, njobs=50, kfold=10, optimized_for='f1_macro'):
        '''
        :param results_file:
        :param params:
        :param njobs:
        :param kfold:
        :return:
        '''
        if params is None:
            params = KNN.parameter_tuning
        self.CV = KFoldCrossVal(self.X, self.Y, folds=kfold)
        self.CV.tune_and_evaluate(self.model, parameters=params, score=optimized_for, file_name=results_file + '_KNN',
                                  n_jobs=njobs)

    def tune_and_eval_predefined(self, results_file, isolates, folds_file, test_file, params=None, njobs=50, optimized_for='f1_macro'):
        '''
        :param results_file:
        :param isolates:
        :param folds:
        :param params:
        :param njobs:
        :return:
        '''
        if params is None:
            params = KNN.parameter_tuning
        self.CV = PredefinedFoldCrossVal(self.X, self.Y, isolates, folds_file, test_file)
        self.CV.tune_and_evaluate(self.model, parameters=params, score=optimized_for, file_name=results_file + '_KNN',
                                  n_jobs=njobs)
