function dist = JSDiv(p,q)

m = (p + q)/2;
zeros_idx = m==0;
p(zeros_idx) = [];
q(zeros_idx) = [];
m(zeros_idx) = [];

log_m = log(m);

log_p = log(p);
log_p(log_p==-Inf) = 0;

log_q = log(q);
log_q(log_q==-Inf) = 0;

KLDiv1 = p * log_p' - p * log_m';
KLDiv2 = q * log_q' - q * log_m';
dist = 0.5*KLDiv1 + 0.5*KLDiv2;