(ns fitdistr.core-test
  (:require [clojure.test :refer :all]
            [fitdistr.core :refer :all]
            [fastmath.random :as r]
            [fastmath.core :as m]))

;; mle

(deftest mle-norm
  (let [rng (r/rng :mersenne 1234)
        data (r/->seq (r/distribution :normal {:rng rng}) 100)
        f (fit :mle :normal data)]
    (is (= -132.75894977093353 (get-in f [:stats :mle])))
    (is (= 269.51789954186705 (get-in f [:stats :aic])))
    (is (= 274.72823991384325 (get-in f [:stats :bic])))
    (is (= -0.060565118544590624 (get-in f [:params :mu])))
    (is (= 0.9126990909183662 (get-in f [:params :sd])))
    (is (= :mle (:method f)))
    (is (= :normal (:distribution-name f)))))

(def usa-arrests-assault
  [236 263 294 190 276 204 110 238 335 211  46 120 249 113  56 115 109 249  83
   300 149 255  72 259 178 109 102 252  57 159 285 254 337  45 120 151 159 106
   174 279  86 188 201 120  48 156 145  81  53 161])

(deftest mle-data-fit-poisson
  (let [f (fit :mle :poisson usa-arrests-assault)]
    (is (= 170.76 (get-in f [:params :p])))
    (is (= -1211.7 (m/approx (get-in f [:stats :mle]))))))

(deftest mle-data-fit-poisson-gradient
  (let [f (fit :mle :poisson usa-arrests-assault {:optimizer :gradient})]
    (is (= 170.75999975700722 (get-in f [:params :p])))
    (is (= -1211.7 (m/approx (get-in f [:stats :mle]))))))

(deftest mle-data-fit-poisson-gradient-start
  (let [f (fit :mle :poisson usa-arrests-assault {:optimizer :gradient
                                                  :initial [2]})]
    (is (= 170 (int (get-in f [:params :p]))))
    (is (= -1211.7 (m/approx (get-in f [:stats :mle]))))))


(deftest mle-data-fit-nbinomial
  (let [f (fit :mle :negative-binomial usa-arrests-assault)]
    (is (= 3.82 (m/approx (get-in f [:params :r]) 2)))
    (is (= 0.02 (m/approx (get-in f [:params :p]) 2)))
    (is (= -290.3297 (m/approx (get-in f [:stats :mle]) 4)))))

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
    (is (= 0.66 (m/approx (get-in f [:stats :cvm]))))
    (is (= -1255.63 (m/approx (get-in f [:stats :mle]))))
    (is (= 2.09 (m/approx (get-in f [:params :alpha]))))
    (is (= 82.66 (m/approx (get-in f [:params :beta]))))
    (is (= :cvm (:method f)))))

(deftest mge-ks
  (let [f (fit :ks :weibull gb-serving {:stats [:mle]})]
    (is (= 0.11 (m/approx (get-in f [:stats :ks]))))
    (is (= -1255.98 (m/approx (get-in f [:stats :mle]))))
    (is (= 2.07 (m/approx (get-in f [:params :alpha]))))
    (is (= 81.45 (m/approx (get-in f [:params :beta]))))
    (is (= :ks (:method f)))))

(deftest mge-ad
  (let [f (fit :ad :weibull gb-serving {:stats [:mle]})]
    (is (= 3.5 (m/approx (get-in f [:stats :ad]))))
    (is (= -1255.39 (m/approx (get-in f [:stats :mle]))))
    (is (= 2.13 (m/approx (get-in f [:params :alpha]))))
    (is (= 82.89 (m/approx (get-in f [:params :beta]))))
    (is (= :ad (:method f)))))

(deftest mge-adr
  (let [f (fit :adr :weibull gb-serving {:stats [:mle]})]
    (is (= 1.61 (m/approx (get-in f [:stats :adr]))))
    (is (= -1255.84 (m/approx (get-in f [:stats :mle]))))
    (is (= 2.07 (m/approx (get-in f [:params :alpha]))))
    (is (= 82.76 (m/approx (get-in f [:params :beta]))))
    (is (= :adr (:method f)))))

(deftest mge-adl
  (let [f (fit :adl :weibull gb-serving {:stats [:mle]})]
    (is (= 1.85 (m/approx (get-in f [:stats :adl]))))
    (is (= -1255.41 (m/approx (get-in f [:stats :mle]))))
    (is (= 2.2 (m/approx (get-in f [:params :alpha]))))
    (is (= 82.02 (m/approx (get-in f [:params :beta]))))
    (is (= :adl (:method f)))))

(deftest mge-ad2r
  (let [f (fit :ad2r :weibull gb-serving {:stats [:mle]})]
    (is (= 11.56 (m/approx (get-in f [:stats :ad2r]))))
    (is (= -1259.11 (m/approx (get-in f [:stats :mle]))))
    (is (= 1.9 (m/approx (get-in f [:params :alpha]))))
    (is (= 81.34 (m/approx (get-in f [:params :beta]))))
    (is (= :ad2r (:method f)))))

(deftest mge-ad2l
  (let [f (fit :ad2l :weibull gb-serving {:stats [:mle]})]
    (is (= 9.79 (m/approx (get-in f [:stats :ad2l]))))
    (is (= -1265.95 (m/approx (get-in f [:stats :mle]))))
    (is (= 2.48 (m/approx (get-in f [:params :alpha]))))
    (is (= 78.25 (m/approx (get-in f [:params :beta]))))
    (is (= :ad2l (:method f)))))

(deftest mge-ad2
  (let [f (fit :ad2 :weibull gb-serving {:stats [:mle]})]
    (is (= 26.95 (m/approx (get-in f [:stats :ad2]))))
    (is (= -1256.31 (m/approx (get-in f [:stats :mle]))))
    (is (= 2.08 (m/approx (get-in f [:params :alpha]))))
    (is (= 85.28 (m/approx (get-in f [:params :beta]))))
    (is (= :ad2 (:method f)))))

;; qme

(deftest qme-norm
  (let [rng (r/rng :mersenne 1234)
        data (r/->seq (r/distribution :normal {:rng rng}) 100)
        f (fit :qme :normal data {:quantiles [m/THIRD m/TWO_THIRD]
                                  :stats [:mle]})]
    (is (= -133.00941301140114 (get-in f [:stats :mle])))
    (is (= 6.1693793635523286E-12 (get-in f [:stats :qme])))
    (is (= 270.0188260228023 (get-in f [:stats :aic])))
    (is (= 275.2291663947785 (get-in f [:stats :bic])))
    (is (= -0.001650866695727643 (get-in f [:params :mu])))
    (is (= 0.8960489369460864 (get-in f [:params :sd])))
    (is (= :qme (:method f)))))

(deftest qme-gb
  (let [f (fit :qme :weibull gb-serving {:stats [:mle]
                                         :quantiles 1000})]
    (is (= 39.77 (m/approx (get-in f [:stats :qme]))))
    (is (= -1255.25 (m/approx (get-in f [:stats :mle]))))
    (is (= 2.17 (m/approx (get-in f [:params :alpha]))))
    (is (= 82.87 (m/approx (get-in f [:params :beta]))))
    (is (= :qme (:method f)))))

;; mme

;; qme

(deftest mme-norm
  (let [rng (r/rng :mersenne 1234)
        data (r/->seq (r/distribution :normal {:rng rng}) 100)
        f (fit :mme :normal data {:mse? false
                                  :stats [:mle]})]
    (is (= -132.76146656358625 (get-in f [:stats :mle])))
    (is (= 4.921066679153263E-11 (get-in f [:stats :mme])))
    (is (= 269.5229331271725 (get-in f [:stats :aic])))
    (is (= 274.7332734991487 (get-in f [:stats :bic])))
    (is (= -0.060565118483739126 (get-in f [:params :mu])))
    (is (= 0.9172971003051852 (get-in f [:params :sd])))
    (is (= :mme (:method f)))))

(deftest mme-gb
  (let [f (fit :mme :weibull gb-serving {:stats [:mle]
                                         :quantiles 1000})]
    (is (= 2.206E-11 (m/approx (get-in f [:stats :mme]) 14)))
    (is (= -1255.25 (m/approx (get-in f [:stats :mle]))))
    (is (= 2.16 (m/approx (get-in f [:params :alpha]))))
    (is (= 83.16 (m/approx (get-in f [:params :beta]))))
    (is (= :mme (:method f)))))

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
    (is (= -1257 (int (get-in f [:stats :mle]))))
    (is (= 2.385347199495471 (get-in f [:params :alpha])))
    (is (= 82.38093254929288 (get-in f [:params :beta])))
    (is (= :infer (:method f)))))

;; distribution

(deftest distribution-test
  (let [d (:distribution (infer :weibull gb-serving))]
    (is (= 73 (int (r/mean d))))
    (is (= 1061 (int (r/variance d))))))
