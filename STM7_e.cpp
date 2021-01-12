#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]

NumericMatrix bmv_modelQSS(int mult,
                           NumericVector inits, 
                           
                           double inc_factor, double beta, double risk_sort, 
                           double overest_factor, double theta, double kappa_neg,
                           double deltaID_low_pos,   double deltaLD_low_pos,  double deltaIL_low_pos,  double lambda_low_pos,
                           double deltaID_high_pos,  double deltaLD_high_pos, double deltaIL_high_pos, double lambda_high_pos,
                           double deltaID_low_neg,   double deltaLD_low_neg,  double deltaIL_low_neg,  double lambda_low_neg, 
                           double deltaID_high_neg,  double deltaLD_high_neg, double deltaIL_high_neg, double lambda_high_neg,
                           double HR, double sero_sort_pos, double sero_sort_neg,
                           NumericVector sigma_vector, NumericVector kappa_vector, NumericVector sigma_neg_vector,
                           NumericVector birth, NumericVector death, NumericVector migration){
  
  double Diag_first_low=0;
  double Diag_first_high=0;
  double Diag_second_low=0;
  double Diag_second_high=0;
  double Diag_first=0;
  double Diag_second=0;
  double Diag_pos=0;
  double Diag_neg=0;
  double Diag_all=0;
  
  double Diag_low_pos=0;
  double Diag_high_pos=0;
  double Diag_low_neg=0;
  double Diag_high_neg=0;  
  
  int num_years = 12;
  
  NumericMatrix variables(45,num_years*mult+1);
  variables(0,0) = inits(0);
  variables(1,0) = inits(1);
  variables(2,0) = inits(2);
  variables(3,0) = inits(3);
  variables(4,0) = inits(4);
  variables(5,0) = inits(5);
  variables(6,0) = inits(6);
  variables(7,0) = inits(7);
  
  variables(8,0)  = inits(8);
  variables(9,0)  = inits(9);
  variables(10,0) = inits(10);
  variables(11,0) = inits(11);
  variables(12,0) = inits(12);
  variables(13,0) = inits(13);
  variables(14,0) = inits(14);
  variables(15,0) = inits(15);
  
  variables(16,0) = inits(16);
  variables(17,0) = inits(17);
  variables(18,0) = inits(18);
  variables(19,0) = inits(19);
  variables(20,0) = inits(20);
  variables(21,0) = inits(21);
  variables(22,0) = inits(22);
  variables(23,0) = inits(23);
  
  variables(24,0) = inits(24);
  variables(25,0) = inits(25);
  variables(26,0) = inits(26);
  variables(27,0) = inits(27);
  variables(28,0) = inits(28);
  variables(29,0) = inits(29);
  variables(30,0) = inits(30);
  variables(31,0) = inits(31);
  
  variables(32,0) = Diag_first_low;
  variables(33,0) = Diag_second_low;
  variables(34,0) = Diag_first_high;
  variables(35,0) = Diag_second_high;
  
  variables(36,0) = Diag_first;
  variables(37,0) = Diag_second;
  variables(38,0) = Diag_pos;
  variables(39,0) = Diag_neg;
  variables(40,0) = Diag_all;  
  
  variables(41,0) = Diag_low_pos;
  variables(42,0) = Diag_high_pos;
  variables(43,0) = Diag_low_neg;
  variables(44,0) = Diag_high_neg;  
  
  double t_step = 1.0/mult;
  
  for (int i = 1; i <= num_years*mult; i += 1){
    
    int t = i/ mult;
    int curr_step = 0 + ((t) / 1);
    
    double sum_low_pos  = variables(0,i-1)  + variables(1,i-1)  + variables(2,i-1)  + variables(3,i-1)  + variables(4,i-1)  + variables(5,i-1)  + variables(6,i-1)  + variables(7,i-1);
    double sum_high_pos = variables(8,i-1)  + variables(9,i-1)  + variables(10,i-1) + variables(11,i-1) + variables(12,i-1) + variables(13,i-1) + variables(14,i-1) + variables(15,i-1);
    double sum_low_neg  = variables(16,i-1) + variables(17,i-1) + variables(18,i-1) + variables(19,i-1) + variables(20,i-1) + variables(21,i-1) + variables(22,i-1) + variables(23,i-1);
    double sum_high_neg = variables(24,i-1) + variables(25,i-1) + variables(26,i-1) + variables(27,i-1) + variables(28,i-1) + variables(29,i-1) + variables(30,i-1) + variables(31,i-1);
    //double sum_pos = sum_low_pos + sum_high_pos;
    //double sum_neg = sum_low_neg + sum_high_neg;
    //double sum_low = sum_low_neg + sum_low_pos;
    //double sum_high = sum_high_neg + sum_high_pos;
    //double sum_all = sum_low_pos + sum_high_pos + sum_low_neg + sum_high_neg;
    
    double inf_low_pos  = variables(1,i-1)  + variables(5,i-1);
    double inf_high_pos = variables(9,i-1)  + variables(13,i-1);
    double inf_low_neg  = variables(17,i-1) + variables(21,i-1);
    double inf_high_neg = variables(25,i-1) + variables(29,i-1);
    
    
    double   force_on_low_pos_by_low_pos   = beta*   ((  sero_sort_pos)*(  risk_sort)*inf_low_pos/sum_low_pos)  ;
    double   force_on_low_pos_by_high_pos  = beta*HR*((  sero_sort_pos)*(1-risk_sort)*inf_high_pos/sum_high_pos);
    double   force_on_low_pos_by_low_neg   = beta*   ((1-sero_sort_pos)*(  risk_sort)*inf_low_neg/sum_low_neg)  ;
    double   force_on_low_pos_by_high_neg  = beta*HR*((1-sero_sort_pos)*(1-risk_sort)*inf_high_neg/sum_high_neg);
    
    double   force_on_high_pos_by_low_pos  = beta*HR*((  sero_sort_pos)*(1-risk_sort)*inf_low_pos/sum_low_pos)  ;
    double   force_on_high_pos_by_high_pos = beta*HR*((  sero_sort_pos)*(  risk_sort)*inf_high_pos/sum_high_pos);
    double   force_on_high_pos_by_low_neg  = beta*HR*((1-sero_sort_pos)*(1-risk_sort)*inf_low_neg/sum_low_neg)  ;
    double   force_on_high_pos_by_high_neg = beta*HR*((1-sero_sort_pos)*(  risk_sort)*inf_high_neg/sum_high_neg);
    
    double   force_on_low_neg_by_low_pos   = beta*       ((1-sero_sort_neg)*(  risk_sort)*inf_low_pos/sum_low_pos);
    double   force_on_low_neg_by_high_pos  = beta*    HR*((1-sero_sort_neg)*(1-risk_sort)*inf_high_pos/sum_high_pos);
    double   force_on_low_neg_by_low_neg   = beta*       ((  sero_sort_neg)*(  risk_sort)*inf_low_neg/sum_low_neg);
    double   force_on_low_neg_by_high_neg  = beta*    HR*((  sero_sort_neg)*(1-risk_sort)*inf_high_neg/sum_high_neg);
    
    double   force_on_high_neg_by_low_pos  = beta*    HR*((1-sero_sort_neg)*(1-risk_sort)*inf_low_pos/sum_low_pos);
    double   force_on_high_neg_by_high_pos = beta*    HR*((1-sero_sort_neg)*(  risk_sort)*inf_high_pos/sum_high_pos);
    double   force_on_high_neg_by_low_neg  = beta*    HR*((  sero_sort_neg)*(1-risk_sort)*inf_low_neg/sum_low_neg);
    double   force_on_high_neg_by_high_neg = beta*    HR*((  sero_sort_neg)*(  risk_sort)*inf_high_neg/sum_high_neg);
    
    double   force_on_low_pos  = force_on_low_pos_by_low_pos +  force_on_low_pos_by_high_pos +  force_on_low_pos_by_low_neg +  force_on_low_pos_by_high_neg;
    double   force_on_high_pos = force_on_high_pos_by_low_pos + force_on_high_pos_by_high_pos + force_on_high_pos_by_low_neg + force_on_high_pos_by_high_neg;
    double   force_on_low_neg  = force_on_low_neg_by_low_pos +  force_on_low_neg_by_high_pos +  force_on_low_neg_by_low_neg +  force_on_low_neg_by_high_neg;
    double   force_on_high_neg = force_on_high_neg_by_low_pos + force_on_high_neg_by_high_pos + force_on_high_neg_by_low_neg + force_on_high_neg_by_high_neg;
    
    double total = variables(0,i-1) + variables(1,i-1) + variables(2,i-1) + variables(3,i-1) + variables(4,i-1) + variables(5,i-1) + variables(6,i-1) + variables(7,i-1) + variables(8,i-1) + variables(9,i-1) + variables(10,i-1) + variables(11,i-1) + variables(12,i-1) + variables(13,i-1) + variables(14,i-1) + variables(15,i-1) + variables(16,i-1) + variables(17,i-1) + variables(18,i-1) + variables(19,i-1) + variables(20,i-1) + variables(21,i-1) + variables(22,i-1) + variables(23,i-1) + variables(24,i-1) + variables(25,i-1) + variables(26,i-1) + variables(27,i-1) + variables(28,i-1) + variables(29,i-1) + variables(30,i-1) + variables(31,i-1);
    
    variables(0,i)  =    variables(0,i-1)  + t_step*(            -force_on_low_pos*variables(0,i-1)                                                                                                                                                                                                                                                                                                                                                                       - death[curr_step]*variables(0,i-1)/total  + migration[curr_step]*variables(0,i-1)/total          -sigma_vector[curr_step]*variables(0,i-1)              +kappa_vector[curr_step]*variables(8,i-1) );                   
    variables(1,i)  =    variables(1,i-1)  + t_step*(            +force_on_low_pos*variables(0,i-1)       - deltaIL_low_pos*variables(1,i-1)       - deltaID_low_pos*variables(1,i-1)                                                                                                                                                                                                                                                                                     - death[curr_step]*variables(1,i-1)/total  + migration[curr_step]*variables(1,i-1)/total          -sigma_vector[curr_step]*variables(1,i-1)              +kappa_vector[curr_step]*variables(9,i-1));         
    variables(2,i)  =    variables(2,i-1)  + t_step*(                                                       deltaIL_low_pos*variables(1,i-1)                                              - deltaLD_low_pos*variables(2,i-1)                                                                                                                                                                                                                                              - death[curr_step]*variables(2,i-1)/total  + migration[curr_step]*variables(2,i-1)/total          -sigma_vector[curr_step]*variables(2,i-1)              +kappa_vector[curr_step]*variables(10,i-1));         
    variables(3,i)  =    variables(3,i-1)  + t_step*(                                                                                                deltaID_low_pos*variables(1,i-1)     + deltaLD_low_pos*variables(2,i-1)                                                                                                                                                                                      -lambda_low_pos*variables(3,i-1)                        - death[curr_step]*variables(3,i-1)/total  + migration[curr_step]*variables(3,i-1)/total          -sigma_vector[curr_step]*variables(3,i-1)              +kappa_vector[curr_step]*variables(11,i-1));         
    variables(4,i)  =    variables(4,i-1)  + t_step*(                                                                                                                                                                                        -force_on_low_pos*variables(4,i-1)                                                                                                                                   +lambda_low_pos*(variables(3,i-1)+variables(7,i-1))     - death[curr_step]*variables(4,i-1)/total  + migration[curr_step]*variables(4,i-1)/total          -sigma_vector[curr_step]*variables(4,i-1)              +kappa_vector[curr_step]*variables(12,i-1));      
    variables(5,i)  =    variables(5,i-1)  + t_step*(                                                                                                                                                                                         force_on_low_pos*variables(4,i-1)       - deltaIL_low_pos*variables(5,i-1)           - deltaID_low_pos*variables(5,i-1)                                                                                                     - death[curr_step]*variables(5,i-1)/total  + migration[curr_step]*variables(5,i-1)/total          -sigma_vector[curr_step]*variables(5,i-1)              +kappa_vector[curr_step]*variables(13,i-1));      
    variables(6,i)  =    variables(6,i-1)  + t_step*(                                                                                                                                                                                                                                   deltaIL_low_pos*variables(5,i-1)                                                - deltaLD_low_pos*variables(6,i-1)                                                                - death[curr_step]*variables(6,i-1)/total  + migration[curr_step]*variables(6,i-1)/total          -sigma_vector[curr_step]*variables(6,i-1)              +kappa_vector[curr_step]*variables(14,i-1));      
    variables(7,i)  =    variables(7,i-1)  + t_step*(                                                                                                                                                                                                                                                                                deltaID_low_pos*variables(5,i-1)   + deltaLD_low_pos*variables(6,i-1)        -lambda_low_pos*variables(7,i-1)                        - death[curr_step]*variables(7,i-1)/total  + migration[curr_step]*variables(7,i-1)/total          -sigma_vector[curr_step]*variables(7,i-1)              +kappa_vector[curr_step]*variables(15,i-1));      
    
    variables(8,i)  =    variables(8,i-1)  + t_step*(            -force_on_high_pos*variables(8,i-1)                                                                                                                                                                                                                                                                                                                                                                      - death[curr_step]*variables(8,i-1)/total  + migration[curr_step]*variables(8,i-1)/total          +sigma_vector[curr_step]*variables(0,i-1)              -kappa_vector[curr_step]*variables(8,i-1)      +theta*variables(24,i-1));                    
    variables(9,i)  =    variables(9,i-1)  + t_step*(            +force_on_high_pos*variables(8,i-1)     - deltaIL_high_pos*variables(9,i-1)       - deltaID_high_pos*variables(9,i-1)                                                                                                                                                                                                                                                                                    - death[curr_step]*variables(9,i-1)/total  + migration[curr_step]*variables(9,i-1)/total          +sigma_vector[curr_step]*variables(1,i-1)              -kappa_vector[curr_step]*variables(9,i-1)      +theta*variables(25,i-1)); 
    variables(10,i) =    variables(10,i-1) + t_step*(                                                      deltaIL_high_pos*variables(9,i-1)                                              - deltaLD_high_pos*variables(10,i-1)                                                                                                                                                                                                                                            - death[curr_step]*variables(10,i-1)/total + migration[curr_step]*variables(10,i-1)/total         +sigma_vector[curr_step]*variables(2,i-1)              -kappa_vector[curr_step]*variables(10,i-1)     +theta*variables(26,i-1)); 
    variables(11,i) =    variables(11,i-1) + t_step*(                                                                                                deltaID_high_pos*variables(9,i-1)    + deltaLD_high_pos*variables(10,i-1)                                                                                                                                                                                    -lambda_high_pos*variables(11,i-1)                      - death[curr_step]*variables(11,i-1)/total + migration[curr_step]*variables(11,i-1)/total         +sigma_vector[curr_step]*variables(3,i-1)              -kappa_vector[curr_step]*variables(11,i-1)     +theta*variables(27,i-1)); 
    variables(12,i) =    variables(12,i-1) + t_step*(                                                                                                                                                                                        -force_on_high_pos*variables(12,i-1)                                                                                                                                 +lambda_high_pos*(variables(11,i-1)+variables(15,i-1))  - death[curr_step]*variables(12,i-1)/total + migration[curr_step]*variables(12,i-1)/total         +sigma_vector[curr_step]*variables(4,i-1)              -kappa_vector[curr_step]*variables(12,i-1)     +theta*variables(28,i-1)); 
    variables(13,i) =    variables(13,i-1) + t_step*(                                                                                                                                                                                         force_on_high_pos*variables(12,i-1)     - deltaIL_high_pos*variables(13,i-1)         - deltaID_high_pos*variables(13,i-1)                                                                                                   - death[curr_step]*variables(13,i-1)/total + migration[curr_step]*variables(13,i-1)/total         +sigma_vector[curr_step]*variables(5,i-1)              -kappa_vector[curr_step]*variables(13,i-1)     +theta*variables(29,i-1)); 
    variables(14,i) =    variables(14,i-1) + t_step*(                                                                                                                                                                                                                                   deltaIL_high_pos*variables(13,i-1)                                               - deltaLD_high_pos*variables(14,i-1)                                                             - death[curr_step]*variables(14,i-1)/total + migration[curr_step]*variables(14,i-1)/total         +sigma_vector[curr_step]*variables(6,i-1)              -kappa_vector[curr_step]*variables(14,i-1)     +theta*variables(30,i-1));
    variables(15,i) =    variables(15,i-1) + t_step*(                                                                                                                                                                                                                                                                                deltaID_high_pos*variables(13,i-1)  + deltaLD_high_pos*variables(14,i-1)     -lambda_high_pos*variables(15,i-1)                      - death[curr_step]*variables(15,i-1)/total + migration[curr_step]*variables(15,i-1)/total         +sigma_vector[curr_step]*variables(7,i-1)              -kappa_vector[curr_step]*variables(15,i-1)     +theta*variables(31,i-1)); 
    
    
    variables(16,i) =    variables(16,i-1) + t_step*(            -force_on_low_neg*variables(16,i-1)                                                                                                                                                                                                                                                                                                                                                                      - death[curr_step]*variables(16,i-1)/total + migration[curr_step]*variables(16,i-1)/total         -sigma_neg_vector[curr_step]*variables(16,i-1)         +kappa_neg*variables(24,i-1)                                                 +birth[curr_step]);                        
    variables(17,i) =    variables(17,i-1) + t_step*(            +force_on_low_neg*variables(16,i-1)     - deltaIL_low_neg*variables(17,i-1)         - deltaID_low_neg*variables(17,i-1)                                                                                                                                                                                                                                                                                  - death[curr_step]*variables(17,i-1)/total + migration[curr_step]*variables(17,i-1)/total         -sigma_neg_vector[curr_step]*variables(17,i-1)         +kappa_neg*variables(25,i-1));       
    variables(18,i) =    variables(18,i-1) + t_step*(                                                      deltaIL_low_neg*variables(17,i-1)                                              - deltaLD_low_neg*variables(18,i-1)                                                                                                                                                                                                                                             - death[curr_step]*variables(18,i-1)/total + migration[curr_step]*variables(18,i-1)/total         -sigma_neg_vector[curr_step]*variables(18,i-1)         +kappa_neg*variables(26,i-1));     
    variables(19,i) =    variables(19,i-1) + t_step*(                                                                                                  deltaID_low_neg*variables(17,i-1)  + deltaLD_low_neg*variables(18,i-1)                                                                                                                                                                                     -lambda_low_neg*variables(19,i-1)                       - death[curr_step]*variables(19,i-1)/total + migration[curr_step]*variables(19,i-1)/total         -sigma_neg_vector[curr_step]*variables(19,i-1)         +kappa_neg*variables(27,i-1));       
    variables(20,i) =    variables(20,i-1) + t_step*(                                                                                                                                                                                        -force_on_low_neg*variables(20,i-1)                                                                                                                                  +lambda_low_neg*(variables(19,i-1)+variables(23,i-1))   - death[curr_step]*variables(20,i-1)/total + migration[curr_step]*variables(20,i-1)/total         -sigma_neg_vector[curr_step]*variables(20,i-1)         +kappa_neg*variables(28,i-1));       
    variables(21,i) =    variables(21,i-1) + t_step*(                                                                                                                                                                                         force_on_low_neg*variables(20,i-1)      - deltaIL_low_neg*variables(21,i-1)         - deltaID_low_neg*variables(21,i-1)                                                                                                     - death[curr_step]*variables(21,i-1)/total + migration[curr_step]*variables(21,i-1)/total         -sigma_neg_vector[curr_step]*variables(21,i-1)         +kappa_neg*variables(29,i-1));       
    variables(22,i) =    variables(22,i-1) + t_step*(                                                                                                                                                                                                                                   deltaIL_low_neg*variables(21,i-1)                                               - deltaLD_low_neg*variables(22,i-1)                                                               - death[curr_step]*variables(22,i-1)/total + migration[curr_step]*variables(22,i-1)/total         -sigma_neg_vector[curr_step]*variables(22,i-1)         +kappa_neg*variables(30,i-1));       
    variables(23,i) =    variables(23,i-1) + t_step*(                                                                                                                                                                                                                                                                               deltaID_low_neg*variables(21,i-1)   + deltaLD_low_neg*variables(22,i-1)       -lambda_low_neg*variables(23,i-1)                       - death[curr_step]*variables(23,i-1)/total + migration[curr_step]*variables(23,i-1)/total         -sigma_neg_vector[curr_step]*variables(23,i-1)         +kappa_neg*variables(31,i-1));       
    
    variables(24,i) =    variables(24,i-1) + t_step*(           -force_on_high_neg*variables(24,i-1)                                                                                                                                                                                                                                                                                                                                                                      - death[curr_step]*variables(24,i-1)/total + migration[curr_step]*variables(24,i-1)/total         +sigma_neg_vector[curr_step]*variables(16,i-1)         -kappa_neg*variables(24,i-1)                   -theta*variables(24,i-1));                  
    variables(25,i) =    variables(25,i-1) + t_step*(           +force_on_high_neg*variables(24,i-1)     - deltaIL_high_neg*variables(25,i-1)       - deltaID_high_neg*variables(25,i-1)                                                                                                                                                                                                                                                                                  - death[curr_step]*variables(25,i-1)/total + migration[curr_step]*variables(25,i-1)/total         +sigma_neg_vector[curr_step]*variables(17,i-1)         -kappa_neg*variables(25,i-1)                   -theta*variables(25,i-1));
    variables(26,i) =    variables(26,i-1) + t_step*(                                                      deltaIL_high_neg*variables(25,i-1)                                             - deltaLD_high_neg*variables(26,i-1)                                                                                                                                                                                                                                            - death[curr_step]*variables(26,i-1)/total + migration[curr_step]*variables(26,i-1)/total         +sigma_neg_vector[curr_step]*variables(18,i-1)         -kappa_neg*variables(26,i-1)                   -theta*variables(26,i-1));
    variables(27,i) =    variables(27,i-1) + t_step*(                                                                                                 deltaID_high_neg*variables(25,i-1)  + deltaLD_high_neg*variables(26,i-1)                                                                                                                                                                                    -lambda_high_neg*variables(27,i-1)                      - death[curr_step]*variables(27,i-1)/total + migration[curr_step]*variables(27,i-1)/total         +sigma_neg_vector[curr_step]*variables(19,i-1)         -kappa_neg*variables(27,i-1)                   -theta*variables(27,i-1));
    variables(28,i) =    variables(28,i-1) + t_step*(                                                                                                                                                                                        -force_on_high_neg*variables(28,i-1)                                                                                                                                 +lambda_high_neg*(variables(27,i-1)+variables(31,i-1))  - death[curr_step]*variables(28,i-1)/total + migration[curr_step]*variables(28,i-1)/total         +sigma_neg_vector[curr_step]*variables(20,i-1)         -kappa_neg*variables(28,i-1)                   -theta*variables(28,i-1));
    variables(29,i) =    variables(29,i-1) + t_step*(                                                                                                                                                                                         force_on_high_neg*variables(28,i-1)     - deltaIL_high_neg*variables(29,i-1)        - deltaID_high_neg*variables(29,i-1)                                                                                                    - death[curr_step]*variables(29,i-1)/total + migration[curr_step]*variables(29,i-1)/total         +sigma_neg_vector[curr_step]*variables(21,i-1)         -kappa_neg*variables(29,i-1)                   -theta*variables(29,i-1));
    variables(30,i) =    variables(30,i-1) + t_step*(                                                                                                                                                                                                                                   deltaIL_high_neg*variables(29,i-1)                                              - deltaLD_high_neg*variables(30,i-1)                                                              - death[curr_step]*variables(30,i-1)/total + migration[curr_step]*variables(30,i-1)/total         +sigma_neg_vector[curr_step]*variables(22,i-1)         -kappa_neg*variables(30,i-1)                   -theta*variables(30,i-1));
    variables(31,i) =    variables(31,i-1) + t_step*(                                                                                                                                                                                                                                                                               deltaID_high_neg*variables(29,i-1)  + deltaLD_high_neg*variables(30,i-1)      -lambda_high_neg*variables(31,i-1)                      - death[curr_step]*variables(31,i-1)/total + migration[curr_step]*variables(31,i-1)/total         +sigma_neg_vector[curr_step]*variables(23,i-1)         -kappa_neg*variables(31,i-1)                   -theta*variables(31,i-1));
    
    //Diag_first_low, Diag_second_low, Diag_first_high, Diag_second_high
    variables(32,i) =    variables(32,i-1) + t_step*(      deltaID_low_pos*variables(1,i-1)    + deltaLD_low_pos*variables(2,i-1));
    variables(33,i) =    variables(33,i-1) + t_step*(      deltaID_low_pos*variables(5,i-1)    + deltaLD_low_pos*variables(6,i-1));
    variables(34,i) =    variables(34,i-1) + t_step*(      deltaID_high_pos*variables(9,i-1)   + deltaLD_high_pos*variables(10,i-1));
    variables(35,i) =    variables(35,i-1) + t_step*(      deltaID_high_pos*variables(13,i-1)  + deltaLD_high_pos*variables(14,i-1));
    
    //Diag_first, Diag_second in HIV
    variables(36,i) =    variables(36,i-1) + t_step*(      deltaID_low_pos*variables(1,i-1)   + deltaLD_low_pos*variables(2,i-1) + deltaID_high_pos*variables(9,i-1)    + deltaLD_high_pos*variables(10,i-1));
    variables(37,i) =    variables(37,i-1) + t_step*(      deltaID_low_pos*variables(5,i-1)   + deltaLD_low_pos*variables(6,i-1) + deltaID_high_pos*variables(13,i-1)   + deltaLD_high_pos*variables(14,i-1));
    
    //Diag_pos, Diag_neg
    variables(38,i) =    variables(38,i-1) + t_step*(      deltaID_low_pos*variables(1,i-1)   + deltaLD_low_pos*variables(2,i-1)  + deltaID_low_pos*variables(5,i-1)    + deltaLD_low_pos*variables(6,i-1)  + deltaID_high_pos*variables(9,i-1)    + deltaLD_high_pos*variables(10,i-1) + deltaID_high_pos*variables(13,i-1)   + deltaLD_high_pos*variables(14,i-1));
    variables(39,i) =    variables(39,i-1) + t_step*(      deltaID_low_neg*variables(17,i-1)  + deltaLD_low_neg*variables(18,i-1) + deltaID_low_neg*variables(21,i-1)   + deltaLD_low_neg*variables(22,i-1) + deltaID_high_neg*variables(25,i-1)   + deltaLD_high_neg*variables(26,i-1) + deltaID_high_neg*variables(29,i-1)   + deltaLD_high_neg*variables(30,i-1));
    
    variables(40,i) =    variables(40,i-1) + t_step*(      deltaID_low_pos*variables(1,i-1)   + deltaLD_low_pos*variables(2,i-1)  + deltaID_low_pos*variables(5,i-1)    + deltaLD_low_pos*variables(6,i-1)  + deltaID_high_pos*variables(9,i-1)    + deltaLD_high_pos*variables(10,i-1) + deltaID_high_pos*variables(13,i-1)   + deltaLD_high_pos*variables(14,i-1)+
      deltaID_low_neg*variables(17,i-1)  + deltaLD_low_neg*variables(18,i-1) + deltaID_low_neg*variables(21,i-1)   + deltaLD_low_neg*variables(22,i-1) + deltaID_high_neg*variables(25,i-1)   + deltaLD_high_neg*variables(26,i-1) + deltaID_high_neg*variables(29,i-1)   + deltaLD_high_neg*variables(30,i-1));
    
    //Diag_low_pos, Diag_high_pos, Diag_low_neg, Diag_high_neg
    variables(41,i) =    variables(41,i-1) + t_step*(      deltaID_low_pos*variables(1,i-1)    + deltaLD_low_pos*variables(2,i-1)   + deltaID_low_pos*variables(5,i-1)    + deltaLD_low_pos*variables(6,i-1));
    variables(42,i) =    variables(42,i-1) + t_step*(      deltaID_high_pos*variables(9,i-1)   + deltaLD_high_pos*variables(10,i-1) + deltaID_high_pos*variables(13,i-1)  + deltaLD_high_pos*variables(14,i-1));
    variables(43,i) =    variables(43,i-1) + t_step*(      deltaID_low_neg*variables(17,i-1)   + deltaLD_low_neg*variables(18,i-1)  + deltaID_low_neg*variables(21,i-1)   + deltaLD_low_neg*variables(22,i-1));
    variables(44,i) =    variables(44,i-1) + t_step*(      deltaID_high_neg*variables(25,i-1)  + deltaLD_high_neg*variables(26,i-1) + deltaID_high_neg*variables(29,i-1)  + deltaLD_high_neg*variables(30,i-1));
    
    
    
  }
  //NumericMatrix out(13,41);
  //for(int i = 0; i <= 12; i += 1) {
  //for (int j = 0; j <= 40; j+=1){
  //out(i,j) = variables(i*1*mult, j);
  //}
  //}
  //return(out);
  return(variables);
}