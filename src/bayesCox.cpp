#define R_PACKAGE

#include "TimeIndepCoxModel.h"
#include "TimeVaryingCoxModel.h"
#include "DynamicCoxModel.h"
#include "DynamicCoxModel_v2.h"
#include "GibbsSampler.h"

using namespace ir;

extern "C" {
  /* Bayesian Cox model */
  void bayesCox(double *p_LRX, int *p_N, int *p_nBeta,
                double *p_grid, int *p_K,
                char **p_out, int *p_id,
                double *p_bp1, double *p_bp2, double *p_cp1, double *p_cp2,
                int *p_iter, int *p_burn, int *p_thin, int *p_verbose, int *p_nReport,
                double *p_a0, double *p_eps0,
                double *p_lambda, double *p_beta, double *p_nu, int *p_jump,
                double *p_LPML, double *p_DHat, double *p_DBar, double *p_pD, double *p_DIC)
  {
    const Size N = p_N[0];
    const Size nBeta = p_nBeta[0];
    const Size K = p_K[0];

    ublas::matrix<double> LRX(N, nBeta + 2, 0.0);
    for (Size i = 0; i < LRX.size1(); ++i)
      for (Size j = 0; j < LRX.size2(); ++j)
        LRX(i, j) = p_LRX[i + N * j];

    ublas::vector<double> grid(K, 1.0);
    for (Size k = 0; k < grid.size(); ++k)
      grid(k) = p_grid[k];

    /* Pointer to data */
    boost::shared_ptr<IntRegData> pd(new IntRegData(LRX, grid));

    GetRNGstate();

    const int id = p_id[0];

    // TimeIndepCox + Gamma + Normal
    if (id == 11) {
      typedef CoxPrior<GammaPrior, NormalPrior> P;
      P prior(GammaPrior(p_bp1[0], p_bp2[0]), NormalPrior(p_cp1[0], p_cp2[0]));

      typedef TimeIndepCoxModel<P> M;
      boost::shared_ptr<M> pm(new M(pd));

      GibbsSampler<M> gs(pm, p_iter[0]);
      gs.runGibbs(prior, static_cast<bool>(p_verbose[0]), p_nReport[0]);

      TimeIndepCoxPar par(pm->initPar());
      gs.summaryFitR(p_burn[0], p_thin[0], par,
                     p_LPML[0], p_DHat[0], p_DBar[0], p_pD[0], p_DIC[0]);

      //par.print(std::cout);

      for (Size k = 0; k < K; ++k)
        p_lambda[k] = par.lambda(k);

      for (Size j = 0; j < nBeta; ++j)
        p_beta[j] = par.beta(j);

      std::string ss(p_out[0]);
      std::ofstream os(ss.c_str());
      gs.outputSample(os);
    }

    // TimeIndepCox + GammaProcess + Normal
    if (id == 12) {
      typedef CoxPrior<GammaProcessPrior, NormalPrior> P;
      P prior(GammaProcessPrior(p_bp1[0], p_bp2[0]), NormalPrior(p_cp1[0], p_cp2[0]));

      typedef TimeIndepCoxModel<P> M;
      boost::shared_ptr<M> pm(new M(pd));

      GibbsSampler<M> gs(pm, p_iter[0]);
      gs.runGibbs(prior, static_cast<bool>(p_verbose[0]), p_nReport[0]);

      TimeIndepCoxPar par(pm->initPar());
      gs.summaryFitR(p_burn[0], p_thin[0], par,
                     p_LPML[0], p_DHat[0], p_DBar[0], p_pD[0], p_DIC[0]);

      //par.print(std::cout);

      for (Size k = 0; k < K; ++k)
        p_lambda[k] = par.lambda(k);

      for (Size j = 0; j < nBeta; ++j)
        p_beta[j] = par.beta(j);

      std::string ss(p_out[0]);
      std::ofstream os(ss.c_str());
      gs.outputSample(os);
    }

    // TimeVaryingCox + Gamma + NormalProcess
    if (id == 21) {
      typedef CoxPrior<GammaPrior, NormalProcessPrior> P;
      P prior(GammaPrior(p_bp1[0], p_bp2[0]), NormalProcessPrior(p_cp1[0]));

      typedef TimeVaryingCoxModel<P> M;
      boost::shared_ptr<M> pm(new M(pd));
      pm->set_a0(p_a0[0]);

      GibbsSampler<M> gs(pm, p_iter[0]);
      gs.runGibbs(prior, static_cast<bool>(p_verbose[0]), p_nReport[0]);

      TimeVaryingCoxPar par(pm->initPar());
      gs.summaryFitR(p_burn[0], p_thin[0], par,
                     p_LPML[0], p_DHat[0], p_DBar[0], p_pD[0], p_DIC[0]);

      for (Size k = 0; k < K; ++k)
        p_lambda[k] = par.lambda(k);

      for (Size j = 0; j < nBeta; ++j)
        for (Size k = 0; k < K; ++k)
          p_beta[k + K * j] = par.beta(k, j);

      std::string ss(p_out[0]);
      std::ofstream os(ss.c_str());
      gs.outputSample(os);
    }

    // TimeVaryingCox + Gamma + NormalInvGammaProcess
    if (id == 22) {
      typedef CoxPrior<GammaPrior, NormalInvGammaProcessPrior> P;
      P prior(GammaPrior(p_bp1[0], p_bp2[0]),
              NormalInvGammaProcessPrior(p_cp1[0], p_cp2[0]));

      typedef TimeVaryingCoxModel<P> M;
      boost::shared_ptr<M> pm(new M(pd));
      pm->set_a0(p_a0[0]);

      GibbsSampler<M> gs(pm, p_iter[0]);
      gs.runGibbs(prior, static_cast<bool>(p_verbose[0]), p_nReport[0]);

      TimeVaryingCoxPar par(pm->initPar());
      gs.summaryFitR(p_burn[0], p_thin[0], par,
                     p_LPML[0], p_DHat[0], p_DBar[0], p_pD[0], p_DIC[0]);

      for (Size k = 0; k < K; ++k)
        p_lambda[k] = par.lambda(k);

      for (Size j = 0; j < nBeta; ++j)
        for (Size k = 0; k < K; ++k)
          p_beta[k + K * j] = par.beta(k, j);

      for (Size j = 0; j < nBeta; ++j)
        p_nu[j] = par.nu(j);

      std::string ss(p_out[0]);
      std::ofstream os(ss.c_str());
      gs.outputSample(os);
    }

    // TimeVaryingCox + GammaProcess + NormalProcess
    if (id == 23) {
      typedef CoxPrior<GammaProcessPrior, NormalProcessPrior> P;
      P prior(GammaProcessPrior(p_bp1[0], p_bp2[0]), NormalProcessPrior(p_cp1[0]));

      typedef TimeVaryingCoxModel<P> M;
      boost::shared_ptr<M> pm(new M(pd));
      pm->set_a0(p_a0[0]);

      GibbsSampler<M> gs(pm, p_iter[0]);
      gs.runGibbs(prior, static_cast<bool>(p_verbose[0]), p_nReport[0]);

      TimeVaryingCoxPar par(pm->initPar());
      gs.summaryFitR(p_burn[0], p_thin[0], par,
                     p_LPML[0], p_DHat[0], p_DBar[0], p_pD[0], p_DIC[0]);

      for (Size k = 0; k < K; ++k)
        p_lambda[k] = par.lambda(k);

      for (Size j = 0; j < nBeta; ++j)
        for (Size k = 0; k < K; ++k)
          p_beta[k + K * j] = par.beta(k, j);

      std::string ss(p_out[0]);
      std::ofstream os(ss.c_str());
      gs.outputSample(os);
    }

    // TimeVaryingCox + GammaProcess + NormalInvGammaProcess
    if (id == 24) {
      typedef CoxPrior<GammaProcessPrior, NormalInvGammaProcessPrior> P;
      P prior(GammaProcessPrior(p_bp1[0], p_bp2[0]),
              NormalInvGammaProcessPrior(p_cp1[0], p_cp2[0]));

      typedef TimeVaryingCoxModel<P> M;
      boost::shared_ptr<M> pm(new M(pd));
      pm->set_a0(p_a0[0]);

      GibbsSampler<M> gs(pm, p_iter[0]);
      gs.runGibbs(prior, static_cast<bool>(p_verbose[0]), p_nReport[0]);

      TimeVaryingCoxPar par(pm->initPar());
      gs.summaryFitR(p_burn[0], p_thin[0], par,
                     p_LPML[0], p_DHat[0], p_DBar[0], p_pD[0], p_DIC[0]);

      for (Size k = 0; k < K; ++k)
        p_lambda[k] = par.lambda(k);

      for (Size j = 0; j < nBeta; ++j)
        for (Size k = 0; k < K; ++k)
          p_beta[k + K * j] = par.beta(k, j);

      for (Size j = 0; j < nBeta; ++j)
        p_nu[j] = par.nu(j);

      std::string ss(p_out[0]);
      std::ofstream os(ss.c_str());
      gs.outputSample(os);
    }

    // Dynamic + Gamma + NormalProcess
    if (id == 31) {
      typedef CoxPrior<GammaPrior, NormalProcessPrior> P;
      P prior(GammaPrior(p_bp1[0], p_bp2[0]), NormalProcessPrior(p_cp1[0]));

      typedef DynamicCoxModel<P> M;
      boost::shared_ptr<M> pm(new M(pd));
      pm->set_a0(p_a0[0]);
      pm->set_eps0(p_eps0[0]);

      GibbsSampler<M> gs(pm, p_iter[0]);
      gs.runGibbs(prior, static_cast<bool>(p_verbose[0]), p_nReport[0]);

      DynamicCoxPar par(pm->initPar());
      gs.summaryFitR(p_burn[0], p_thin[0], par,
                     p_LPML[0], p_DHat[0], p_DBar[0], p_pD[0], p_DIC[0]);

      for (Size k = 0; k < K; ++k)
        p_lambda[k] = par.lambda(k);

      for (Size j = 0; j < nBeta; ++j) {
        for (Size k = 0; k < K; ++k) {
          p_beta[k + K * j] = par.beta(k, j);
          p_jump[k + K * j] = par.jump(k, j);
        }
      }

      std::string ss(p_out[0]);
      std::ofstream os(ss.c_str());
      gs.outputSample(os);
    }

    // DynamicCox + Gamma + NormalInvGammaProcess
    if (id == 32) {
      typedef CoxPrior<GammaPrior, NormalInvGammaProcessPrior> P;
      P prior(GammaPrior(p_bp1[0], p_bp2[0]),
              NormalInvGammaProcessPrior(p_cp1[0], p_cp2[0]));

      typedef DynamicCoxModel<P> M;
      boost::shared_ptr<M> pm(new M(pd));
      pm->set_a0(p_a0[0]);
      pm->set_eps0(p_eps0[0]);

      GibbsSampler<M> gs(pm, p_iter[0]);
      gs.runGibbs(prior, static_cast<bool>(p_verbose[0]), p_nReport[0]);

      DynamicCoxPar par(pm->initPar());
      gs.summaryFitR(p_burn[0], p_thin[0], par,
                     p_LPML[0], p_DHat[0], p_DBar[0], p_pD[0], p_DIC[0]);

      for (Size k = 0; k < K; ++k)
        p_lambda[k] = par.lambda(k);

      for (Size j = 0; j < nBeta; ++j) {
        for (Size k = 0; k < K; ++k) {
          p_beta[k + K * j] = par.beta(k, j);
          p_jump[k + K * j] = par.jump(k, j);
        }
      }

      for (Size j = 0; j < nBeta; ++j)
        p_nu[j] = par.nu(j);

      std::string ss(p_out[0]);
      std::ofstream os(ss.c_str());
      gs.outputSample(os);
    }

    // Dynamic + Const + NormalProcess
    if (id == 33) {
      typedef CoxPrior<ConstValuePrior, NormalProcessPrior> P;
      P prior(ConstValuePrior(1), NormalProcessPrior(p_cp1[0]));

      typedef DynamicCoxModel_v2<P> M;
      boost::shared_ptr<M> pm(new M(pd));
      pm->set_a0(p_a0[0]);
      pm->set_eps0(p_eps0[0]);

      GibbsSampler<M> gs(pm, p_iter[0]);
      gs.runGibbs(prior, static_cast<bool>(p_verbose[0]), p_nReport[0]);

      DynamicCoxPar par(pm->initPar());
      gs.summaryFitR(p_burn[0], p_thin[0], par,
                     p_LPML[0], p_DHat[0], p_DBar[0], p_pD[0], p_DIC[0]);

      for (Size k = 0; k < K; ++k)
        p_lambda[k] = par.lambda(k);

      for (Size j = 0; j < nBeta; ++j) {
        for (Size k = 0; k < K; ++k) {
          p_beta[k + K * j] = par.beta(k, j);
          p_jump[k + K * j] = par.jump(k, j);
        }
      }

      std::string ss(p_out[0]);
      std::ofstream os(ss.c_str());
      gs.outputSample(os);
    }

    // DynamicCox + Const + NormalInvGammaProcess
    if (id == 34) {
      typedef CoxPrior<ConstValuePrior, NormalInvGammaProcessPrior> P;
      P prior(ConstValuePrior(1), NormalInvGammaProcessPrior(p_cp1[0], p_cp2[0]));

      typedef DynamicCoxModel_v2<P> M;
      boost::shared_ptr<M> pm(new M(pd));
      pm->set_a0(p_a0[0]);
      pm->set_eps0(p_eps0[0]);

      GibbsSampler<M> gs(pm, p_iter[0]);
      gs.runGibbs(prior, static_cast<bool>(p_verbose[0]), p_nReport[0]);

      DynamicCoxPar par(pm->initPar());
      gs.summaryFitR(p_burn[0], p_thin[0], par,
                     p_LPML[0], p_DHat[0], p_DBar[0], p_pD[0], p_DIC[0]);

      for (Size k = 0; k < K; ++k)
        p_lambda[k] = par.lambda(k);

      for (Size j = 0; j < nBeta; ++j) {
        for (Size k = 0; k < K; ++k) {
          p_beta[k + K * j] = par.beta(k, j);
          p_jump[k + K * j] = par.jump(k, j);
        }
      }

      for (Size j = 0; j < nBeta; ++j)
        p_nu[j] = par.nu(j);

      std::string ss(p_out[0]);
      std::ofstream os(ss.c_str());
      gs.outputSample(os);
    }

    // Dynamic + GammaProcess + NormalProcess
    if (id == 35) {
      typedef CoxPrior<GammaProcessPrior, NormalProcessPrior> P;
      P prior(GammaProcessPrior(p_bp1[0], p_bp2[0]), NormalProcessPrior(p_cp1[0]));

      typedef DynamicCoxModel<P> M;
      boost::shared_ptr<M> pm(new M(pd));
      pm->set_a0(p_a0[0]);
      pm->set_eps0(p_eps0[0]);

      GibbsSampler<M> gs(pm, p_iter[0]);
      gs.runGibbs(prior, static_cast<bool>(p_verbose[0]), p_nReport[0]);

      DynamicCoxPar par(pm->initPar());
      gs.summaryFitR(p_burn[0], p_thin[0], par,
                     p_LPML[0], p_DHat[0], p_DBar[0], p_pD[0], p_DIC[0]);

      for (Size k = 0; k < K; ++k)
        p_lambda[k] = par.lambda(k);

      for (Size j = 0; j < nBeta; ++j) {
        for (Size k = 0; k < K; ++k) {
          p_beta[k + K * j] = par.beta(k, j);
          p_jump[k + K * j] = par.jump(k, j);
        }
      }

      std::string ss(p_out[0]);
      std::ofstream os(ss.c_str());
      gs.outputSample(os);
    }

    // DynamicCox + GammaProcess + NormalInvGammaProcess
    if (id == 36) {
      typedef CoxPrior<GammaProcessPrior, NormalInvGammaProcessPrior> P;
      P prior(GammaProcessPrior(p_bp1[0], p_bp2[0]),
              NormalInvGammaProcessPrior(p_cp1[0], p_cp2[0]));

      typedef DynamicCoxModel<P> M;
      boost::shared_ptr<M> pm(new M(pd));
      pm->set_a0(p_a0[0]);
      pm->set_eps0(p_eps0[0]);

      GibbsSampler<M> gs(pm, p_iter[0]);
      gs.runGibbs(prior, static_cast<bool>(p_verbose[0]), p_nReport[0]);

      DynamicCoxPar par(pm->initPar());
      gs.summaryFitR(p_burn[0], p_thin[0], par,
                     p_LPML[0], p_DHat[0], p_DBar[0], p_pD[0], p_DIC[0]);

      for (Size k = 0; k < K; ++k)
        p_lambda[k] = par.lambda(k);

      for (Size j = 0; j < nBeta; ++j) {
        for (Size k = 0; k < K; ++k) {
          p_beta[k + K * j] = par.beta(k, j);
          p_jump[k + K * j] = par.jump(k, j);
        }
      }

      for (Size j = 0; j < nBeta; ++j)
        p_nu[j] = par.nu(j);

      std::string ss(p_out[0]);
      std::ofstream os(ss.c_str());
      gs.outputSample(os);
    }

    PutRNGstate();
  }
}
