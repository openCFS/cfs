[1mdiff --git a/source/Materials/Models/EBHysteresis.cc b/source/Materials/Models/EBHysteresis.cc[m
[1mindex ffd435877..8701964af 100644[m
[1m--- a/source/Materials/Models/EBHysteresis.cc[m
[1m+++ b/source/Materials/Models/EBHysteresis.cc[m
[36m@@ -253,12 +253,13 @@[m [mDEFINE_LOG(eb, "EBHysteresis")[m
           if (numS_ > 1 ){ // hysteretic case[m
         if (jacobian_method_ == 1){ // use finite differences[m
           mu = EvaluateLocalMuFiniteDifferences(HVec, B_k, idx);[m
[31m-        } else if (jacobian_method_ == 2) {[m
[32m+[m[32m        } else if (jacobian_method_ == 2) { // use Broyden method[m
           mu = EvaluateLocalMuGBM(delta_H, delta_B, idx);[m
[32m+[m[32m        } else if (jacobian_method_ == 3) { // use simple finite differences[m
[32m+[m[32m          mu = EvaluateLocalMu(delta_H, delta_B, idx);[m
         } else {[m
           EXCEPTION("WRONG Jacobian_method!")[m
         }[m
[31m-        //mu = EvaluateLocalMuDFP(delta_H, delta_B, idx);[m
         M_dummy = Evaluate(HVec, true, idx);[m
       } else { // nonlinear case (only anhysteresis)[m
         mu = EvaluateLocalMuAnhystersisOnly(HVec, idx);[m
[36m@@ -266,9 +267,6 @@[m [mDEFINE_LOG(eb, "EBHysteresis")[m
       mu_[idx] = mu;[m
       M0_[idx] = M1_[idx];[m
       H0_[idx] = H1_[idx];[m
[31m-[m
[31m-[m
[31m-[m
       // mark this element as computed[m
       hasElemSolution_[idx] = true;[m
       return mu;[m
