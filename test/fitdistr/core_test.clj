(ns fitdistr.core-test
  (:require [clojure.test :refer [deftest is are]]
            [fitdistr.core :refer [fit bootstrap infer]]
            [fastmath.random :as r]
            [fastmath.core :as m]))

;; mle

;; N(0,1)

(def data [1.037568748877954 0.4771290201211284 0.6219697545985581 1.1796874270758235 -0.35424246293826306 0.3317089451323486 1.6567403486280303 0.09057886852474355 0.4114945430925313 -1.1966192741318902 0.7202179437969064 0.7733418244480172 -1.9929498516556228 -1.960234384631804 -1.4448575262428653 0.1162158786913258 -0.6439796931745334 -1.6104598614289811 -0.1882550424891307 -0.33565810632696885 -1.7493483701425638 -0.05088229211078229 -1.4982967525772628 -0.501140975142851 -0.17251619345061983 2.7322055506588683 0.36417712385979356 0.27965841926389345 -0.7911848532338075 -1.5564373350452367 -0.6960777645865691 0.078056255180603 0.5868709617607102 0.11312447661186921 -1.4323219812406904 0.6150418728809778 -1.1608835122327008 -2.4876492317214933 -0.13260008846789326 -1.6033868480261717 0.06911248425723576 -0.7058361774538193 0.05066357002647885 1.0973738282417742 -1.7522152096007264 0.1861645231455176 -1.105149890420801 -1.0501668294999578 -0.11267875137708652 1.4059311317279193 -1.6334991912177053 -2.4140577909577314 -0.45101965766172847 -1.4700585216817426 0.08550010652895919 0.2593915575837272 -1.5508506035148941 -0.5957653724332678 1.7540590976893686 -1.113940098647868 0.5218673641610219 1.0905296396484028 -0.5241530199073918 0.08299946117814815 -0.009270496365232757 -1.4766892071511633 -0.12044097182397016 2.8293248708536876 -0.5720481662740541 -0.21013603282610965 -1.2046901058501953 -0.7950521955836996 0.11795316635937708 -1.1551714062817526 -0.8398586392549506 0.8494787846547617 0.056149236987679915 -0.1730676445147461 0.5047079372171236 1.171590003300769 0.5523491484751555 1.3900012839451132 0.6329317071555962 1.5082756776920583 0.18962830379577802 0.5546328263299833 1.5302108160019807 0.39657041241776053 -1.897553548695323 1.836202282128761 0.40919407356835846 -2.8585552383408084 0.0752740478502172 1.7401770136743764 -1.1581844440384823 0.3247627534908861 0.7295941222461297 0.6570703756554964 1.0674541593726836 0.4919663039495628])

(deftest mle-norm
  (let [
        f (fit :mle :normal data)]
    (are [v k] (m/delta= v (get-in f k))
      -153.4408123294556 [:stats :mle]
      310.8816246589112 [:stats :aic]
      316.0919650308874 [:stats :bic]
      -0.12105211577857974 [:params :mu]
      1.1224003819345674 [:params :sd])
    (is (= :mle (:method f)))
    (is (= :normal (:distribution-name f)))))

(def usa-arrests-assault
  [236 263 294 190 276 204 110 238 335 211  46 120 249 113  56 115 109 249  83
   300 149 255  72 259 178 109 102 252  57 159 285 254 337  45 120 151 159 106
   174 279  86 188 201 120  48 156 145  81  53 161])

(deftest mle-data-fit-poisson
  (let [f (fit :mle :poisson usa-arrests-assault)]
    (is (m/delta= 170.76 (get-in f [:params :p])))
    (is (m/delta= -1211.7 (m/approx (get-in f [:stats :mle]))))))

(deftest mle-data-fit-poisson-gradient
  (let [f (fit :mle :poisson usa-arrests-assault {:optimizer :gradient})]
    (is (m/delta= 170.76 (m/approx (get-in f [:params :p]))))
    (is (m/delta= -1211.7 (m/approx (get-in f [:stats :mle]))))))

(deftest mle-data-fit-poisson-gradient-start
  (let [f (fit :mle :poisson usa-arrests-assault {:optimizer :gradient
                                                  :initial [2]})]
    (is (m/delta= 170 (int (get-in f [:params :p]))))
    (is (m/delta= -1211.7 (m/approx (get-in f [:stats :mle]))))))


(deftest mle-data-fit-nbinomial
  (let [f (fit :mle :negative-binomial usa-arrests-assault)]
    (is (m/delta= 3.82 (m/approx (get-in f [:params :r]) 2)))
    (is (m/delta= 0.02 (m/approx (get-in f [:params :p]) 2)))
    (is (m/delta= -290.3297 (m/approx (get-in f [:stats :mle]) 4)))))

;; mge

(def gb-serving
  [30.0  10.0  20.0  24.0  20.0  24.0  40.0  20.0  50.0  30.0  26.0  44.0
   25.0  20.0  40.0 100.0  30.0  80.0  50.0  40.0  20.0  50.0  50.0  50.0
   50.0  80.0  80.0  20.0  60.0  50.0  50.0  50.0  50.0  60.0  50.0  50.0
   50.0  50.0  17.0 120.0  80.0 100.0  80.0 100.0  50.0  70.0  50.0  50.0
   50.0  30.0  90.0 125.0  50.0  25.0  50.0  75.0 100.0  60.0 100.0  11.5
   100.0 100.0  50.0  40.0  50.0  50.0  50.0 150.0  50.0  80.0 100.0 100.0
   80.0  20.0  30.0 100.0  70.0 100.0  80.0  80.0  50.0  50.0  50.0  70.0
   25.0  40.0 100.0  25.0 130.0 100.0  60.0  60.0 100.0  50.0  50.0  50.0
   50.0  50.0  50.0 100.0  30.0  50.0  50.0  50.0 100.0  50.0  50.0  75.0
   25.0  20.0  75.0  50.0  50.0 100.0 100.0  80.0 130.0 200.0  40.0  64.0
   40.0  80.0  80.0  79.0  80.0 105.0  51.0 130.0 105.0 105.0  80.0  80.0
   40.0  40.0  50.0 105.0 105.0  80.0 105.0  40.0  50.0 130.0 105.0  80.0
   40.0  80.0  50.0 100.0 200.0  79.0 150.0  80.0  80.0  80.0  80.0  80.0
   65.0 130.0  40.0  80.0  80.0  80.0 130.0  51.0 130.0  80.0  80.0 130.0
   51.0  25.0  79.0  80.0  50.0  51.0  50.0 105.0  80.0  80.0 130.0 130.0
   80.0  51.0  50.0 130.0 115.0 115.0 200.0 160.0 105.0  80.0  80.0  80.0
   130.0 130.0  80.0  51.0  50.0 130.0 105.0 105.0  80.0 105.0  80.0  60.0
   51.0 105.0  80.0  51.0  40.0 130.0 130.0  40.0  32.0  80.0  25.0  80.0
   80.0  80.0  40.0 130.0 130.0  40.0  75.0  40.0  80.0  80.0 150.0  52.5
   130.0  80.0  80.0  75.0  40.0  25.5  51.0  50.0  80.0  51.0 130.0 130.0
   51.0 130.0  80.0 130.0  80.0 130.0 130.0  25.5 117.0 130.0  80.0  80.0
   80.0  80.0])

(deftest mge-cvm
  (let [f (fit :cvm :weibull gb-serving {:stats [:mle]})]
    (is (m/delta= 0.66 (m/approx (get-in f [:stats :cvm]))))
    (is (m/delta= -1255.63 (m/approx (get-in f [:stats :mle]))))
    (is (m/delta= 2.09 (m/approx (get-in f [:params :alpha]))))
    (is (m/delta= 82.66 (m/approx (get-in f [:params :beta]))))
    (is (= :cvm (:method f)))))

(deftest mge-ks
  (let [f (fit :ks :weibull gb-serving {:stats [:mle]})]
    (is (m/delta= 0.11 (m/approx (get-in f [:stats :ks]))))
    (is (m/delta= -1255.98 (m/approx (get-in f [:stats :mle]))))
    (is (m/delta= 2.07 (m/approx (get-in f [:params :alpha]))))
    (is (m/delta= 81.45 (m/approx (get-in f [:params :beta]))))
    (is (= :ks (:method f)))))

(deftest mge-ad
  (let [f (fit :ad :weibull gb-serving {:stats [:mle]})]
    (is (m/delta= 3.5 (m/approx (get-in f [:stats :ad]))))
    (is (m/delta= -1255.39 (m/approx (get-in f [:stats :mle]))))
    (is (m/delta= 2.13 (m/approx (get-in f [:params :alpha]))))
    (is (m/delta= 82.89 (m/approx (get-in f [:params :beta]))))
    (is (= :ad (:method f)))))

(deftest mge-adr
  (let [f (fit :adr :weibull gb-serving {:stats [:mle]})]
    (is (m/delta= 1.61 (m/approx (get-in f [:stats :adr]))))
    (is (m/delta= -1255.84 (m/approx (get-in f [:stats :mle]))))
    (is (m/delta= 2.07 (m/approx (get-in f [:params :alpha]))))
    (is (m/delta= 82.76 (m/approx (get-in f [:params :beta]))))
    (is (= :adr (:method f)))))

(deftest mge-adl
  (let [f (fit :adl :weibull gb-serving {:stats [:mle]})]
    (is (m/delta= 1.85 (m/approx (get-in f [:stats :adl]))))
    (is (m/delta= -1255.41 (m/approx (get-in f [:stats :mle]))))
    (is (m/delta= 2.2 (m/approx (get-in f [:params :alpha]))))
    (is (m/delta= 82.02 (m/approx (get-in f [:params :beta]))))
    (is (= :adl (:method f)))))

(deftest mge-ad2r
  (let [f (fit :ad2r :weibull gb-serving {:stats [:mle]})]
    (is (m/delta= 11.56 (m/approx (get-in f [:stats :ad2r]))))
    (is (m/delta= -1259.11 (m/approx (get-in f [:stats :mle]))))
    (is (m/delta= 1.9 (m/approx (get-in f [:params :alpha]))))
    (is (m/delta= 81.34 (m/approx (get-in f [:params :beta]))))
    (is (= :ad2r (:method f)))))

(deftest mge-ad2l
  (let [f (fit :ad2l :weibull gb-serving {:stats [:mle]})]
    (is (m/delta= 9.79 (m/approx (get-in f [:stats :ad2l]))))
    (is (m/delta= -1265.95 (m/approx (get-in f [:stats :mle]))))
    (is (m/delta= 2.48 (m/approx (get-in f [:params :alpha]))))
    (is (m/delta= 78.25 (m/approx (get-in f [:params :beta]))))
    (is (= :ad2l (:method f)))))

(deftest mge-ad2
  (let [f (fit :ad2 :weibull gb-serving {:stats [:mle]})]
    (is (m/delta= 26.95 (m/approx (get-in f [:stats :ad2]))))
    (is (m/delta= -1256.31 (m/approx (get-in f [:stats :mle]))))
    (is (m/delta= 2.08 (m/approx (get-in f [:params :alpha]))))
    (is (m/delta= 85.28 (m/approx (get-in f [:params :beta]))))
    (is (= :ad2 (:method f)))))

;; qme

(deftest qme-norm
  (let [f (fit :qme :normal data {:quantiles [m/THIRD m/TWO_THIRD]
                                  :stats [:mle]})]
    (are [v k] (m/delta= v (get-in f k))
      -153.49898608311454 [:stats :mle]
      7.139135716668963E-12 [:stats :qme]
      310.9979721662291 [:stats :aic]
      316.2083125382053 [:stats :bic]
      -0.08959004240194787 [:params :mu]
      1.138457318184922 [:params :sd])
    (is (= :qme (:method f)))))

(deftest qme-gb
  (let [f (fit :qme :weibull gb-serving {:stats [:mle]
                                         :quantiles 1000})]
    (is (m/delta= 39.77 (m/approx (get-in f [:stats :qme]))))
    (is (m/delta= -1255.25 (m/approx (get-in f [:stats :mle]))))
    (is (m/delta= 2.17 (m/approx (get-in f [:params :alpha]))))
    (is (m/delta= 82.87 (m/approx (get-in f [:params :beta]))))
    (is (= :qme (:method f)))))

;; mme

(deftest mme-norm
  (let [rng (r/rng :mersenne 1234)
        data (r/->seq (r/distribution :normal {:rng rng}) 100)
        f (fit :mme :normal data {:mse? false
                                  :stats [:mle]})]
    (is (m/delta= -153.4433291221637 (get-in f [:stats :mle])))
    (is (m/delta= 5.248972784199779E-11 (get-in f [:stats :mme])))
    (is (m/delta= 310.8866582443274 (get-in f [:stats :aic])))
    (is (m/delta= 316.0969986163036 (get-in f [:stats :bic])))
    (is (m/delta= -0.1210521157574994 (get-in f [:params :mu])))
    (is (m/delta= 1.1280548277428102 (get-in f [:params :sd])))
    (is (= :mme (:method f)))))

(deftest mme-gb
  (let [f (fit :mme :weibull gb-serving {:stats [:mle]
                                         :quantiles 1000})]
    (is (m/delta= 2.206E-11 (m/approx (get-in f [:stats :mme]) 14)))
    (is (m/delta= -1255.25 (m/approx (get-in f [:stats :mle]))))
    (is (m/delta= 2.16 (m/approx (get-in f [:params :alpha]))))
    (is (m/delta= 83.16 (m/approx (get-in f [:params :beta]))))
    (is (= :mme (:method f)))))

;; mps

(deftest mps-norm
  (let [rng (r/rng :mersenne 1234)
        data (r/->seq (r/distribution :normal {:rng rng}) 100)
        f (fit :mps :normal data {:mse? false
                                  :stats [:mle]})]
    (is (m/delta= -153.48470657973516 (get-in f [:stats :mle])))
    (is (m/delta= -5.187306023818127 (get-in f [:stats :mps])))
    (is (m/delta= 310.9694131594703 (get-in f [:stats :aic])))
    (is (m/delta= 316.1797535314465 (get-in f [:stats :bic])))
    (is (m/delta= -0.15431246746578386 (get-in f [:params :mu])))
    (is (m/delta= 1.1231921976781158 (get-in f [:params :sd])))
    (is (= :mps (:method f)))))

(deftest mps-gb
  (let [f (fit :mps :weibull gb-serving {:stats [:mle]})]
    (is (m/delta= -4.76 (m/approx (get-in f [:stats :mps]))))
    (is (m/delta= -1255.25 (m/approx (get-in f [:stats :mle]))))
    (is (m/delta= 2.19 (m/approx (get-in f [:params :alpha]))))
    (is (m/delta= 83.88 (m/approx (get-in f [:params :beta]))))
    (is (= :mps (:method f)))))


;; bootstrap

(deftest bootstrap-mle
  (let [f (bootstrap :mle :weibull gb-serving {:size 100
                                               :samples 2000
                                               :optimizer :gradient})]
    (is (= -1255 (int (get-in f [:stats :mle]))))
    (is (>= 2 (int (get-in f [:params :alpha]))))
    (is (< 80 (int (get-in f [:params :beta])) 85))
    (is (:bootstrap? f))
    (is (= :mle (:method f)))))

;; infer

(deftest infer-weibull
  (let [f (infer :weibull gb-serving)]
    (is (m/delta= -1257 (int (get-in f [:stats :mle]))))
    (is (m/delta= 2.385347199495471 (get-in f [:params :alpha])))
    (is (m/delta= 82.38093254929288 (get-in f [:params :beta])))
    (is (= :infer (:method f)))))

;; distribution

(deftest distribution-test
  (let [d (:distribution (infer :weibull gb-serving))]
    (is (m/delta= 73 (int (r/mean d))))
    (is (m/delta= 1061 (int (r/variance d))))))
