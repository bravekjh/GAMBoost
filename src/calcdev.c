#include <math.h>

void pre_sum(double *x, int *n, int *p, double *D, double *weights, double *presum) 
{
    int i, j;

    for (i = 0; i < *p; i++) {
        presum[i] = 0;
        for (j = 0; j < *n; j++) {
            presum[i] += x[j + (*n * i)] * x[j + (*n * i)] * weights[j] * D[j];
        }
    }
}

void get_min_score_dev_binary(double *x, int *n, int *p, 
                     double *y, double *actual_mu, double *actual_eta, double *D, double *weights, double *buffer0,
                     double *penalty, double *buffer,
                     int *min_index, double *score, double *min_dev) 
{
    int i, j;
    double mu;
    double min_score, actual_beta_delta, actual_nominator;
    double best_beta_delta;

    pre_sum(x,n,p,D,weights,buffer);
    for (j = 0; j < *n; j++) buffer0[j] = weights[j] * (y[j] - actual_mu[j]);
    
    for (i = 0; i < *p; i++) {
        actual_nominator = 0;
        for (j = 0; j < *n; j++) actual_nominator += x[j + (*n * i)] * buffer0[j];
        actual_beta_delta = actual_nominator / (buffer[i]+penalty[i]);
        score[i] = actual_nominator * actual_beta_delta;
        if (i == 0 || score[i] > min_score) {
            min_score = score[i];
            best_beta_delta = actual_beta_delta;
            *min_index = i;
        } 
    }

    *min_dev = 0;
    for (j = 0; j < *n; j++) {
        mu = 1/(1+exp(- (x[j + (*n * *min_index)] * best_beta_delta + actual_eta[j])));
        *min_dev +=  -2*(y[j] == 1 ? log(mu) : log(1 - mu))*weights[j];
    }
}



void calc_dev_binary(double *x, int *n, int *p, 
                     double *y, double *actual_mu, double *actual_eta, double *D, double *weights, double *buffer0,
                     double *penalty, double *buffer, double *beta_delta_vec,
                     double *dev) 
{
    int i, j;
    double mu;

    pre_sum(x,n,p,D,weights,buffer);
    for (j = 0; j < *n; j++) buffer0[j] = weights[j] * (y[j] - actual_mu[j]);
    
    for (i = 0; i < *p; i++) {
        beta_delta_vec[i] = 0;
        for (j = 0; j < *n; j++) beta_delta_vec[i] += x[j + (*n * i)] * buffer0[j];
        beta_delta_vec[i] *= 1/(buffer[i]+penalty[i]);
    }

    for (i = 0; i < *p; i++) {
        dev[i] = 0;
        for (j = 0; j < *n; j++) {
            mu = 1/(1+exp(- (x[j + (*n * i)] * beta_delta_vec[i] + actual_eta[j])));
            dev[i] +=  -2*(y[j] == 1 ? log(mu) : log(1 - mu))*weights[j];
        }
    }
}
