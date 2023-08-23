#ifndef LIONCPP_MATH_IPOPT_NLP_COMPLETE_OUTPUT_H    
#define LIONCPP_MATH_IPOPT_NLP_COMPLETE_OUTPUT_H    

namespace lioncpp {

    struct NLP_complete_output
    {
        Ipopt::Index number_of_variables;
        Ipopt::Index number_of_constraints;
        Ipopt::Index jacobian_nnz;
        Ipopt::Index hessian_nnz;
        std::vector<Ipopt::Index> col_jac;
        std::vector<Ipopt::Index> row_jac;
        std::vector<Ipopt::Index> col_hes;
        std::vector<Ipopt::Index> row_hes;
        Ipopt::Number fitness_function;
        std::vector<Ipopt::Number> constraints;
        std::vector<Ipopt::Number> grad_f;
        std::vector<Ipopt::Number> jacobian_constraints;
        std::vector<Ipopt::Number> hessian_lagrangian;
    };

}

#endif
