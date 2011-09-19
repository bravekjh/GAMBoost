#include <math.h>

void pre_sum_index(double *x, int *n, int *p, 
                   int *n_total, int *n_index, int *p_index,
                   double *D, double *weights, double *presum) 
{
    int i, j;
    double actual_x;

    for (i = 0; i < *p; i++) {
        presum[i] = 0;
        for (j = 0; j < *n; j++) {
            actual_x = x[n_index[j] + (*n_total * p_index[i])];
            presum[i] += actual_x * actual_x * weights[j] * D[j];
        }
    }
}


void get_min_score_core(double *x, int *n, int *p,
                     int *n_total, int *n_index, int *p_index, 
                     double *y, double *actual_mu, double *weights, double *buffer0,
                     double *penalty, double *buffer,
                     int *min_index, double *score, double *best_beta_delta) 
{
    int i, j;
    double min_score, actual_beta_delta, actual_nominator;

    for (j = 0; j < *n; j++) buffer0[j] = weights[j] * (y[j] - actual_mu[j]);
    
    for (i = 0; i < *p; i++) {
        actual_nominator = 0;
        for (j = 0; j < *n; j++) actual_nominator += x[n_index[j] + (*n_total * p_index[i])] * buffer0[j];
        actual_beta_delta = actual_nominator / (buffer[i]+penalty[i]);
        score[i] = actual_nominator * actual_beta_delta;
        if (i == 0 || score[i] > min_score) {
            min_score = score[i];
            *best_beta_delta = actual_beta_delta;
            *min_index = i;
        } 
    }
}


void get_min_score_dev_binary_index(double *x, int *n, int *p,
                     int *n_total, int *n_index, int *p_index, 
                     double *y, double *actual_mu, double *actual_eta, double *D, double *weights, double *buffer0,
                     double *penalty, double *buffer,
                     int *min_index, double *score, double *min_dev) 
{
    int j;
    double mu;
    double best_beta_delta;

    pre_sum_index(x,n,p,n_total,n_index,p_index,D,weights,buffer);
    get_min_score_core(x,n,p,n_total,n_index,p_index,y,actual_mu,weights,buffer0,penalty,buffer,min_index,score,&best_beta_delta);

    *min_dev = 0;
    for (j = 0; j < *n; j++) {
        mu = 1/(1+exp(- (x[n_index[j] + (*n_total * p_index[*min_index])] * best_beta_delta + actual_eta[j])));
        *min_dev +=  -2*(y[j] == 1 ? log(mu) : log(1 - mu))*weights[j];
    }
}


void get_min_score_dev_general_index(double *x, int *n, int *p,
                     int *n_total, int *n_index, int *p_index, 
                     double *y, double *actual_mu, double *actual_eta, double *D, double *weights, double *buffer0,
                     double *penalty, double *buffer,
                     int *min_index, double *score, double *best_beta_delta) 
{
    pre_sum_index(x,n,p,n_total,n_index,p_index,D,weights,buffer);
    get_min_score_core(x,n,p,n_total,n_index,p_index,y,actual_mu,weights,buffer0,penalty,buffer,min_index,score,best_beta_delta);
}


void get_min_score_dev_gaussian_index(double *x, int *n, int *p, 
                     int *n_total, int *n_index, int *p_index,
                     double *y, double *actual_mu, double *actual_eta, double *D, double *weights, double *buffer0,
                     double *penalty, double *pre_sum,
                     int *min_index, double *score, double *min_dev) 
{
    int j;
    double mu;
    double best_beta_delta;

    get_min_score_core(x,n,p,n_total,n_index,p_index,y,actual_mu,weights,buffer0,penalty,pre_sum,min_index,score,&best_beta_delta);

    *min_dev = 0;
    for (j = 0; j < *n; j++) {
        mu = x[n_index[j] + (*n_total * p_index[*min_index])] * best_beta_delta + actual_eta[j];
        *min_dev +=  (y[j] - mu)*(y[j] - mu)*weights[j];
    }
}



void calc_dev_binary_index(double *x, int *n, int *p, 
                     int *n_total, int *n_index, int *p_index, 
                     double *y, double *actual_mu, double *actual_eta, double *D, double *weights, double *buffer0,
                     double *penalty, double *buffer, double *beta_delta_vec,
                     double *dev) 
{
    int i, j;
    double mu;

    pre_sum_index(x,n,p,n_total,n_index,p_index,D,weights,buffer);
    for (j = 0; j < *n; j++) buffer0[j] = weights[j] * (y[j] - actual_mu[j]);
    
    for (i = 0; i < *p; i++) {
        beta_delta_vec[i] = 0;
        for (j = 0; j < *n; j++) beta_delta_vec[i] += x[n_index[j] + (*n_total * p_index[i])] * buffer0[j];
        beta_delta_vec[i] *= 1/(buffer[i]+penalty[i]);
    }

    for (i = 0; i < *p; i++) {
        dev[i] = 0;
        for (j = 0; j < *n; j++) {
            mu = 1/(1+exp(- (x[n_index[j] + (*n_total * p_index[i])] * beta_delta_vec[i] + actual_eta[j])));
            dev[i] +=  -2*(y[j] == 1 ? log(mu) : log(1 - mu))*weights[j];
        }
    }
}
