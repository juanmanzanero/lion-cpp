#ifndef LION_MATH_IPOPT_GENERIC_NLP_INTERFACE_H
#define LION_MATH_IPOPT_GENERIC_NLP_INTERFACE_H

//
// An auxiliary class that derives from Ipopt::TNLP and can be used as an interface with Ipopt
// This class is useful so that our own NLP class does not need to inherit from Ipopt::TNLP, which might
// be annoying, mainly because needs to instantiate its class as an Ipopt::SmartPtr<Ipopt::TNLP>.
//
// This class takes a pointer to our NLP, and overrides all the virtual functions that Ipopt calls
// to take a step with calls to our optimal control NLP transcription
//

namespace lioncpp {

    template<typename NLP_type>
    class Ipopt_generic_NLP_interface : public Ipopt::TNLP
    {
    public:

        Ipopt_generic_NLP_interface(NLP_type& nlp) : _nlp(&nlp) {}

        virtual bool get_nlp_info(Ipopt::Index& n, Ipopt::Index& m, Ipopt::Index& nnz_jac_g, Ipopt::Index& nnz_h_lag, Ipopt::TNLP::IndexStyleEnum& index_style) override
        {
            return _nlp->get_nlp_info(n, m, nnz_jac_g, nnz_h_lag, index_style);
        }

        virtual bool get_bounds_info(Ipopt::Index n, Ipopt::Number* x_l, Ipopt::Number* x_u, Ipopt::Index m, Ipopt::Number* g_l, Ipopt::Number* g_u) override
        {
            return _nlp->get_bounds_info(n, x_l, x_u, m, g_l, g_u);
        }

        virtual bool get_starting_point(Ipopt::Index n, bool init_x, Ipopt::Number* x, bool init_z, Ipopt::Number* z_L, Ipopt::Number* z_U, Ipopt::Index m, bool init_lambda, Ipopt::Number* lambda) override
        {
            return _nlp->get_starting_point(n, init_x, x, init_z, z_L, z_U, m, init_lambda, lambda);
        }

        virtual bool eval_f(Ipopt::Index n, const Ipopt::Number* x, bool new_x, Ipopt::Number& obj_value) override
        {
            return _nlp->eval_f(n, x, new_x, obj_value);
        }

        virtual bool eval_grad_f(Ipopt::Index n, const Ipopt::Number* x, bool new_x, Ipopt::Number* grad_f) override
        {
            return _nlp->eval_grad_f(n, x, new_x, grad_f);
        }

        virtual bool eval_g(Ipopt::Index n, const Ipopt::Number* x, bool new_x, Ipopt::Index m, Ipopt::Number* g) override
        {
            return _nlp->eval_g(n, x, new_x, m, g);
        }

        virtual bool eval_jac_g(Ipopt::Index n, const Ipopt::Number* x, bool new_x, Ipopt::Index m, Ipopt::Index nele_jac, Ipopt::Index* iRow, Ipopt::Index* jCol, Ipopt::Number* values) override
        {
            return _nlp->eval_jac_g(n, x, new_x, m, nele_jac, iRow, jCol, values);
        }

        virtual bool eval_h(Ipopt::Index n, const Ipopt::Number* x, bool new_x, Ipopt::Number obj_factor, Ipopt::Index m, const Ipopt::Number* lambda, bool new_lambda, Ipopt::Index nele_hess, Ipopt::Index* iRow, Ipopt::Index* jCol, Ipopt::Number* values) override
        {
            return _nlp->eval_h(n, x, new_x, obj_factor, m, lambda, new_lambda, nele_hess, iRow, jCol, values);
        }

        virtual void finalize_solution(Ipopt::SolverReturn status, Ipopt::Index n, const Ipopt::Number* x, const Ipopt::Number* z_L, const Ipopt::Number* z_U, Ipopt::Index m, const Ipopt::Number* g, const Ipopt::Number* lambda, Ipopt::Number obj_value, const Ipopt::IpoptData* ip_data, Ipopt::IpoptCalculatedQuantities* ip_cq) override
        {
            _nlp->finalize_solution(status, n, x, z_L, z_U, m, g, lambda, obj_value, ip_data, ip_cq);
        }

        virtual bool intermediate_callback(Ipopt::AlgorithmMode mode, Ipopt::Index iter, Ipopt::Number obj_value, Ipopt::Number inf_pr, Ipopt::Number inf_du, Ipopt::Number mu, Ipopt::Number d_norm, Ipopt::Number regularization_size, Ipopt::Number alpha_du, Ipopt::Number alpha_pr, Ipopt::Index ls_trials, const Ipopt::IpoptData* ip_data, Ipopt::IpoptCalculatedQuantities* ip_cq) override
        {
            return _nlp->intermediate_callback(mode, iter, obj_value, inf_pr, inf_du, mu, d_norm, regularization_size, alpha_du, alpha_pr, ls_trials, ip_data, ip_cq);
        }

    private:
        NLP_type* _nlp; // A pointer to the NLP
    };
}

#endif
