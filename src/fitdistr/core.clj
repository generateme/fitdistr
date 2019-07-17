(ns fitdistr.core
  "Distribution fitting using MLE, MGE and QME methods."
  (:require [fastmath.core :as m]
            [fastmath.random :as r]
            [fastmath.stats :as stats]
            [fastmath.optimization :as o]
            [fitdistr.distributions :refer :all]))

(set! *warn-on-reflection* true)
(set! *unchecked-math* :warn-on-boxed)
(m/use-primitive-operators)

;; targets

(defn- log-likelihood [data] (fn [distr] (r/log-likelihood distr data)))

(defn- kolmogorov-smirnov
  [data]
  (let [n (count data)
        s (sort data)
        obsp (map #(/ (double %) n) (range 0 (inc n)))
        obspu (vec (rest obsp))
        obspl (vec (butlast obsp))]
    (fn [distr]
      (let [F (map (partial r/cdf distr) s)
            d (map (fn [^double u ^double f ^double l]
                     (max (m/abs (- u f))
                          (m/abs (- f l)))) obspu F obspl)]
        (stats/maximum d)))))

(defn- cramer-von-mises
  [data]
  (let [n (count data)
        n2 (double (+ n n))
        s (sort data)
        f12 (/ (* 12.0 n))
        obsp (mapv (fn [^long i] (/ (dec (* 2 i)) n2)) (range 1 (inc n)))]
    (fn [distr]
      (let [F (map (partial r/cdf distr) s)]
        (+ f12 ^double (reduce (fn [^double x ^double y] (+ x y)) 0.0 (map (fn [^double f ^double o]
                                                                            (m/sq (- f o))) F obsp)))))))

(defn- anderson-darling
  [data]
  (let [n (count data)
        s (sort data)
        obsp (mapv (fn [^long i] (dec (* 2 i))) (range 1 (inc n)))]
    (fn [distr]
      (let [F (map (partial r/cdf distr) s)
            res (- (- n) (stats/mean (map (fn [^double o ^double f ^double fr]
                                            (* o (+ (m/log f) (m/log (- 1.0 fr))))) obsp F (reverse F))))]
        (if (m/nan? res) ##Inf res)))))

(defn- uniform-quantilies
  [^long len]
  (cond
    (< len 2) [0.5]
    (== len 2) [0.2 0.8]
    (== len 3) [0.25 0.5 0.75]
    (== len 4) [0.2 0.4 0.6 0.8]
    (== len 5) [0.05 0.25 0.5 0.75 0.95]
    :else (mapv #(m/norm % 0 (dec len) 0.01 0.99) (range 0 len))))

(defn- qme
  [quantiles strategy measure data]
  (let [qs (cond
             (sequential? quantiles) quantiles
             (number? quantiles) (uniform-quantilies (long quantiles))
             :else (uniform-quantilies 50))
        measure (if (= :mse measure) m/sq m/abs)
        data-qs (stats/quantiles data qs strategy)]
    (fn [distr]
      (stats/mean (map (fn [^double q1 ^double q2]
                         (measure (- q1 q2))) data-qs (map (partial r/icdf distr) qs))))))

;;

(defn- method->fn
  [method data params]
  (case method 
    :ks (kolmogorov-smirnov data)
    :cvm (cramer-von-mises data)
    :ad (anderson-darling data)
    :qme (let [{:keys [quantiles strategy qmeasure]
                :or {quantiles 50 strategy :legacy qmeasure :mse}} params] (qme quantiles strategy qmeasure data))
    :mle (log-likelihood data)
    (throw (Exception. (str "Method " method " is not supported.")))))

(defn- aic-bic
  [^double ll ^long cnt ^long npar]
  (let [ll2 (* -2.0 ll)]
    {:aic (+ ll2 (* 2.0 npar))
     :bic (+ ll2 (* (m/log cnt) npar))}))

(defn- calc-stats
  [method distr data params result npar]
  (let [stats (disj (set (:stats params)) method)
        res (merge (when method {method result}) ;;step 1, fill in statistics
                   (->> stats
                        (map #(vector % ((method->fn % data params) distr)))
                        (into {})))
        res (if (contains? res :mle) ;; step 2, add aic and bic if log likelihood is available
              (merge res (aic-bic (:mle res) (count data) npar))
              res)]
    (when res {:stats res})))

(defn- method->opt-fn
  [method]
  (if (#{:ks :cvm :ad :qme} method) o/minimize o/maximize))

(defn- assert-values
  [data bounds validation {:keys [initial]}]
  (assert (validation data) "Data values do not fit required distribution")
  (when initial
    (assert (= (count initial)
               (count bounds)) "Initial values do not match required number of parameters")
    (assert (every? true? (map (fn [^double v [^double x1 ^double x2]]
                                 (<= x1 v x2)) initial bounds))
            "Initial values are out of bounds")))

(defn- fit-
  [target method distribution data {:keys [optimizer]
                                    :or {optimizer :nelder-mead}
                                    :as all}]
  (let [{:keys [param-names bounds inference] :as ddata} (distribution-data distribution)]
    (let [opt-fn (method->opt-fn method) ;; minimize or maximize?
          [pars result] (opt-fn optimizer target (merge {:initial (inference data)
                                                         :max-iters 1000} all {:bounds bounds 
                                                                               :bounded? true
                                                                               :stats? false})) ;; optimize!
          
          conf (zipmap param-names pars) ;; create final distribution parametrization...
          distr (r/distribution distribution conf)] ;; ...and distribution
      (merge (calc-stats method distr data all result (count param-names)) ;; return result and statistics
             {:params conf
              :distribution-name distribution
              :distribution distr
              :method method}))))

(defn- make-target
  [method distribution data all param-names]
  (let [raw-target (method->fn method data all)]
    (fn [& r]
      (let [d (r/distribution distribution (zipmap param-names r))]
        (raw-target d)))))

(defn fit
  "Fit distribution using given method

  * `:mle` - log likelihood
  * `:ad` - Anderson-Darling (default)
  * `:ks` - Kolmogorov-Smirnov
  * `:cvm` - Cramer-von-Mises
  * `:qme` - quantile matching estimation.

  For QME additional parameters can be provided:
  
  * `quantiles` - list of quantiles used to match or number of uniformly distributed (default: `50`).
  * `strategy` - quantile calculation strategy (default: `:legacy`).

  More about strategies:

  * https://generateme.github.io/fastmath/fastmath.stats.html#var-estimation-strategies-list
  * http://commons.apache.org/proper/commons-math/javadocs/api-3.6.1/org/apache/commons/math3/stat/descriptive/rank/Percentile.EstimationType.html

  Other parameters and distribution names see README"
  ([method distribution data] (fit method distribution data {}))
  ([method distribution data {:keys [assert?]
                              :or {assert? true}
                              :as all}]
   (let [{:keys [param-names bounds validation]} (distribution-data distribution)]
     (when assert? (assert-values data bounds validation all))
     (fit- (make-target method distribution data all param-names) method distribution data all))))

;; 

(defn infer
  "Infer parameters computationally"
  ([distribution data] (infer distribution data {}))
  ([distribution data {:keys [assert?]
                       :or {assert? true}
                       :as params}]
   (let [{:keys [param-names validation inference]} (distribution-data distribution)]
     (when (:assert? params) (assert (validation data) "Data values do not fit required distribution"))
     (let [conf (zipmap param-names (inference data))
           distr (r/distribution distribution conf)]
       (merge (calc-stats nil distr data (update params :stats conj :mle) nil (count param-names))
              {:params conf
               :distribution-name distribution
               :distribution distr})))))

;; bootstrap

(defn- ci->fn
  [ci-type]
  (case ci-type
    :stddev-mean stats/stddev-extent
    :mad-median stats/mad-extent
    :sem-mean stats/sem-extent
    :iqr-median stats/percentile-extent
    :adj-median stats/adjacent-values
    :ci stats/ci
    :min-max-mean stats/extent
    stats/mad-extent))

(defn bootstrap 
  "Bootstrapped version of fitting.

  Parameters are the same as [[fit]]
  
  Additional parameters:

  * `:size` - number of bootstrapped sequences (default: 100)
  * `:samples` - number of samples in each sequence (default: 10% of data, minimum 100, maximum 5000 samples)
  * `:ci-type` - confidence interval type (default: `:mad-median`)
  * `:all-params?` - fitted parameters for each sequence (default: false)"
  ([method distribution data] (bootstrap method distribution data {}))
  ([method distribution data {:keys [size samples ci-type all-params? assert?]
                              :or {size 100
                                   samples (m/constrain (* 0.1 (count data)) 100 5000)
                                   ci-type :mad-median
                                   all-params? false
                                   assert? true}
                              :as all}]
   (let [{:keys [param-names bounds validation inference]} (distribution-data distribution)]
     (when assert? (assert-values data bounds validation all))
     (let [bdata (stats/bootstrap data size samples) ;; create sequences of bootstrapped data
           res (pmap #(fit method distribution % (-> all
                                                     (dissoc :stats)
                                                     (assoc :assert? false))) bdata) ;; fit paralelly (do not calculate stats)
           all-params (map vals (map :params res)) ;; extract found parameters
           extents (map (ci->fn ci-type) (apply map vector all-params)) ;; calculate extents for each parameter
           params (map last extents) ;; extract mean or median
           param-names (:param-names (distribution-data distribution)) ;; extract parameter names
           conf (zipmap param-names params) ;; create distribution configuration
           distr (r/distribution distribution conf) ;; create distribution
           res (merge (calc-stats nil distr data (update all :stats conj method) nil (count param-names))
                      {ci-type (zipmap param-names extents)
                       :params conf
                       :distribution-name distribution
                       :distribution distr})]
       (if all-params? ;; maybe you want full list for each resampled data
         (assoc res :all-params all-params)
         res)))))

#_(do (def target (r/->seq (r/distribution :gumbel {:mu 10 :beta 3.21}) 100000))

      (take 10 target)

      (count (:bins (stats/histogram target :sqrt)))

      (time (bootstrap :mle :gumbel atv))

      (def atv [0.6 2.8 182.2 0.8 478.0 1.1 215.0 0.7 7.9 316.2 0.2 17780.0 7.8 100.0 0.9 180.0 0.3 300.9
                0.6 17.5 10.0 0.1 5.8 87.7 4.1 3.5 4.9 7060.0 0.2 360.0 100.8 2.3 12.3 40.0 2.3 0.1
                2.7 2.2 0.4 2.6 0.2 1.0 7.3 3.2 0.8 1.2 33.7 14.0 21.4 7.7 1.0 1.9 0.7 12.6
                3.2 7.3 4.9 4000.0 2.5 6.7 3.0 63.0 6.0 1.6 10.1 1.2 1.5 1.2 30.0 3.2 3.5 1.2
                0.2 1.9 0.7 17.0 2.8 4.8 1.3 3.7 0.2 1.8 2.6 5.9 2.6 6.3 1.4 0.8 670.0 810.0
                1890.0 1800.0 8500.0 21000.0 31.0 20.5 4370.0 1000.0 39891.8
                316.2 6400.0 1000.0 7400.0 31622.8]))

;;----------------------------------------

#_(do

    (def target (r/->seq (r/distribution :weibull {:alpha 0.5 :beta 2.2}) 10000))

    (fit :ad :weibull target {:stats [:mle]})
    ;; => {:stats
    ;;     {:ad 0.19749431207310408,
    ;;      :mle -19126.212671469282,
    ;;      :aic 38256.425342938564,
    ;;      :bic 38270.84602368252},
    ;;     :params {:alpha 0.5014214878565807, :beta 2.203213102262515},
    ;;     :distribution #object[org.apache.commons.math3.distribution.WeibullDistribution 0x430997b7 "org.apache.commons.math3.distribution.WeibullDistribution@430997b7"],
    ;;     :method :ad}

    (bootstrap :mle :weibull target {:stats #{:ad}
                                     :optimizer :nelder-mead})

    ;; => {:stats
    ;;     {:mle -19126.178345014738,
    ;;      :ad 0.35561024021990306,
    ;;      :aic 38256.356690029475,
    ;;      :bic 38270.77737077343},
    ;;     :mad-median
    ;;     {:alpha [0.4910043451070347 0.5056336146263343 0.4983189798666845],
    ;;      :beta [2.1185018316179285 2.326029409552982 2.222265620585455]},
    ;;     :params {:alpha 0.4983189798666845, :beta 2.222265620585455},
    ;;     :distribution #object[org.apache.commons.math3.distribution.WeibullDistribution 0x63a766b9 "org.apache.commons.math3.distribution.WeibullDistribution@63a766b9"]}

    (infer :weibull target {:stats #{:mle :ad}})
    ;; => {:stats
    ;;     {:mle -19126.13369575803,
    ;;      :ad 0.22838225327177497,
    ;;      :aic 38256.26739151606,
    ;;      :bic 38270.68807226002},
    ;;     :params {:alpha 0.5012938746206328, :beta 2.215448048490149},
    ;;     :distribution #object[org.apache.commons.math3.distribution.WeibullDistribution 0x3a0f2314 "org.apache.commons.math3.distribution.WeibullDistribution@3a0f2314"]}

    (def atv [0.6 2.8 182.2 0.8 478.0 1.1 215.0 0.7 7.9 316.2 0.2 17780.0 7.8 100.0 0.9 180.0 0.3 300.9
              0.6 17.5 10.0 0.1 5.8 87.7 4.1 3.5 4.9 7060.0 0.2 360.0 100.8 2.3 12.3 40.0 2.3 0.1
              2.7 2.2 0.4 2.6 0.2 1.0 7.3 3.2 0.8 1.2 33.7 14.0 21.4 7.7 1.0 1.9 0.7 12.6
              3.2 7.3 4.9 4000.0 2.5 6.7 3.0 63.0 6.0 1.6 10.1 1.2 1.5 1.2 30.0 3.2 3.5 1.2
              0.2 1.9 0.7 17.0 2.8 4.8 1.3 3.7 0.2 1.8 2.6 5.9 2.6 6.3 1.4 0.8 670.0 810.0
              1890.0 1800.0 8500.0 21000.0 31.0 20.5 4370.0 1000.0 39891.8
              316.2 6400.0 1000.0 7400.0 31622.8])

    (defn find-best
      [method ds]
      (let [selector (if (= method :mle) last first)]
        (dissoc (->> (map #(fit method % atv {:stats #{:mle :ad :ks :cvm}}) ds)
                     (sort-by (comp method :stats))
                     (selector))
                :distribution)))

    (find-best :mle [:weibull :log-normal :gamma :exponential :normal :pareto])
    ;; => {:stats
    ;;     {:mle -532.4052019871922,
    ;;      :cvm 0.6373592936482382,
    ;;      :ks 0.1672497620724005,
    ;;      :ad 3.4721179220009617,
    ;;      :aic 1068.8104039743844,
    ;;      :bic 1074.0991857726672},
    ;;     :params {:scale 2.553816262077493, :shape 3.147240361221695},
    ;;     :distribution-name :log-normal,
    ;;     :method :mle}

    (find-best :ad [:weibull :log-normal :gamma :exponential :normal :pareto])
    ;; => {:stats
    ;;     {:ad 3.0345123029861156,
    ;;      :cvm 0.4615381958965107,
    ;;      :ks 0.1332827771382316,
    ;;      :mle -532.9364810533066,
    ;;      :aic 1069.8729621066132,
    ;;      :bic 1075.161743904896},
    ;;     :params {:scale 2.2941800698596815, :shape 3.2934516278879205},
    ;;     :distribution-name :log-normal,
    ;;     :method :ad}

    (find-best :ks [:weibull :log-normal :gamma :exponential :normal :pareto])
    ;; => {:stats
    ;;     {:ks 0.07692307692307693,
    ;;      :cvm 0.11739378941793886,
    ;;      :mle ##-Inf,
    ;;      :ad ##Inf,
    ;;      :aic ##Inf,
    ;;      :bic ##Inf},
    ;;     :params {:scale 0.36510648416477365, :shape 0.2649915952623174},
    ;;     :distribution-name :pareto,
    ;;     :method :ks}

    (fit :mle :pareto atv {:stats [:mle :ad :ks :cvm :qme]
                           :optimizer :gradient
                           :max-iters 1000})



    (filter #(= :pareto (:distribution-name %)) (find-best :ad :weibull :log-normal :gamma :exponential :normal :pareto))

    (map (partial r/pdf (:distribution (find-best :ks :weibull :log-normal :gamma :exponential :normal :pareto))) atv)


    (r/pdf (r/distribution :pareto {:scale 1.06 :shape 0.34}) 0.5)

    (def target (r/->seq (r/distribution :pareto {:scale 0.1 :shape 4}) 2000))

    (fit :ks :pareto target {:initial [2 4]
                             :optimizer :multidirectional-simplex
                             })

    ((anderson-darling atv)
     (r/distribution :pareto {:scale 0.9 :shape 1}))

    (every? (complement m/inf?) (map #(m/log %) target))

    (defn- infer-pareto
      [data]
      (let [xs (double-array data)
            m (stats/mean xs)
            scale (m/prev-double (stats/minimum xs))
            shape (/ m (- m scale))]
        [scale shape]))

    (defn- infer-pareto
      [data]
      (let [xs (double-array data)
            m1 (stats/mean xs)
            x1 (stats/minimum xs)
            n (count data)
            nshape (/ (- (* n m1) x1)
                      (- m1 x1))
            shape (/ nshape n)
            scale (* x1 (/ (dec nshape) nshape))]
        [scale shape]))


    (infer-pareto atv)

    #_(do

        (def ^:private distribution-params
          {:beta {:params [:alpha :beta]
                  :bounds [[0 ##Inf] [0 ##Inf]]
                  :inference (constantly [1.0 1.0])
                  :validation (partial every? #(< 0 ^double % 1))}
           :cauchy {:params [:median :scale]
                    :bounds [[##-Inf ##Inf] [0 ##Inf]]
                    :inference (fn [xs] [(stats/median xs) (* 0.5 (stats/iqr xs))])
                    :validation (constantly true)}
           :chi-squared {:params [:degrees-of-freedom]
                         :bounds [[0 ##Inf]]
                         :inference (comp vec stats/mean)
                         :validation (partial every? #(pos? ^double %))}
           :exponential {:params [:mean]
                         :bounds [[0 ##Inf]]
                         :inference (comp vec stats/mean)
                         :validation (partial every? #(pos? ^double %))}
           :f {:params [:numerator-degrees-of-freedom :denominator-degrees-of-freedom]
               :bounds [[0 ##Inf] [0 ##Inf]]
               :inference #(let [m (stats/mean %)
                                 d2 (if (> m 2.0) (/ (+ m m) (dec m)) 1.0)
                                 p (/ d2 (+ d2 2.0))
                                 mode (stats/mode %)
                                 diff (- mode p)
                                 d1 (if (neg? diff) (/ (* -2.0 p) diff) 1.0)]
                             [d1 d2])
               :validation (partial every? (complement #(neg? ^double %)))}})

        (def goal-data (let [d (r/distribution :beta {:alpha 0.1 :beta 4})]
                         (take 10000 (r/->seq d))))


        goal-data

        (stats/iqr goal-data)

        ((partial every? #(< 0 ^double % 1)) [0.4])

        (stats/mean goal-data)



        (stats/mode goal-data)

        (let [t :beta
              xs goal-data
              {:keys [params bounds inference validation]} (distribution-params t)
              f (fn [& r]
                  (r/log-likelihood (r/distribution t (zipmap params r)) xs))]
          (when (validation xs)
            (maximize :nelder-mead f {:bounds bounds :initial (inference xs) :bounded? true})))


        (def atv [0.6     2.8   182.2     0.8   478.0     1.1   215.0     0.7     7.9
                  316.2     0.2 17780.0     7.8   100.0     0.9   180.0     0.3   300.9
                  0.6    17.5    10.0     0.1     5.8    87.7     4.1     3.5     4.9
                  7060.0     0.2   360.0   100.8     2.3    12.3    40.0     2.3     0.1
                  2.7     2.2     0.4     2.6     0.2     1.0     7.3     3.2     0.8
                  1.2    33.7    14.0    21.4     7.7     1.0     1.9     0.7    12.6
                  3.2     7.3     4.9  4000.0     2.5     6.7     3.0    63.0     6.0
                  1.6    10.1     1.2     1.5     1.2    30.0     3.2     3.5     1.2
                  0.2     1.9     0.7    17.0     2.8     4.8     1.3     3.7     0.2
                  1.8     2.6     5.9     2.6     6.3     1.4     0.8   670.0   810.0
                  1890.0  1800.0  8500.0 21000.0    31.0    20.5  4370.0  1000.0 39891.8
                  316.2  6400.0  1000.0  7400.0 31622.8])

        (def ll (r/distribution :log-normal {:scale 2.553816 :shape 3.147240}))

        (r/log-likelihood ll atv)


        (defn anderson-darling-l
          [data]
          (let [n (count data)
                s (sort data)
                f-15 (* -1.5 n)
                obsp (mapv #(dec (* 2 %)) (range 1 (inc n)))]
            (fn [distr]
              (let [F (mapv (partial r/cdf distr) s)]
                (- (+ f-15 (reduce + F)) (stats/mean (map #(* %1 (+ (m/log %2) (m/log (- 1.0 %3)))) obsp F)))))))


        ((anderson-darling atv) ll)

        (count atv)
        )
    )
