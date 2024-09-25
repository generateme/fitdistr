(ns fitdistr.core
  "Distribution fitting using MLE, MGE and QME methods."
  (:require [fastmath.core :as m]
            [fastmath.random :as r]
            [fastmath.stats :as stats]
            [fastmath.stats.bootstrap :as boot]
            [fastmath.optimization :as o]
            [fastmath.protocols :as prot]
            [fitdistr.distributions :refer [distribution-data]]))

(set! *warn-on-reflection* true)
(set! *unchecked-math* :warn-on-boxed)
(m/use-primitive-operators)

;; targets

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
        (+ f12 ^double (reduce m/+ 0.0 (map (fn [^double f ^double o]
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

(defn- anderson-darling-r
  [data]
  (let [n (count data)
        hn (/ n 2)
        s (sort data)
        obsp (mapv (fn [^long i] (dec (* 2 i))) (range 1 (inc n)))]
    (fn [distr]
      (let [F (map (partial r/cdf distr) s)
            res (- hn (* 2.0 ^double (reduce m/+ 0.0 F))
                   (stats/mean (map (fn [^double o ^double fr]
                                      (* o (m/log (- 1.0 fr)))) obsp (reverse F))))]
        (if (m/nan? res) ##Inf res)))))

(defn- anderson-darling-l
  [data]
  (let [n (count data)
        n32 (* 1.5 n)
        s (sort data)
        obsp (mapv (fn [^long i] (dec (* 2 i))) (range 1 (inc n)))]
    (fn [distr]
      (let [F (map (partial r/cdf distr) s)
            res (- (* 2.0 ^double (reduce m/+ 0.0 F)) n32
                   (stats/mean (map (fn [^double o ^double f]
                                      (* o (m/log f))) obsp F)))]
        (if (m/nan? res) ##Inf res)))))

(defn- anderson-darling-2r
  [data]
  (let [n (count data)
        s (sort data)
        obsp (mapv (fn [^long i] (dec (* 2 i))) (range 1 (inc n)))]
    (fn [distr]
      (let [F (map (partial r/cdf distr) s)
            res (+ (* 2.0 ^double (reduce m/+ 0.0 (map (fn [^double f]
                                                         (m/log (- 1.0 f))) F)))
                   (stats/mean (map (fn [^double o ^double fr]
                                      (/ o (- 1.0 fr))) obsp (reverse F))))]
        (if (m/nan? res) ##Inf res)))))

(defn- anderson-darling-2l
  [data]
  (let [n (count data)
        s (sort data)
        obsp (mapv (fn [^long i] (dec (* 2 i))) (range 1 (inc n)))]
    (fn [distr]
      (let [F (map (partial r/cdf distr) s)
            res (+ (* 2.0 ^double (reduce m/+ 0.0 (map (fn [^double f]
                                                         (m/log f)) F)))
                   (stats/mean (map (fn [^double o ^double f]
                                      (/ o f)) obsp F)))]
        (if (m/nan? res) ##Inf res)))))

(defn- anderson-darling-2
  [data]
  (let [n (count data)
        s (sort data)
        obsp (mapv (fn [^long i] (dec (* 2 i))) (range 1 (inc n)))]
    (fn [distr]
      (let [F (map (partial r/cdf distr) s)
            res (+ (* 2.0 ^double (reduce m/+ 0.0 (map (fn [^double f]
                                                         (+ (m/log f) (m/log (- 1.0 f)))) F)))
                   (stats/mean (map (fn [^double o ^double f ^double fr]
                                      (+ (/ o f) (/ o (- 1.0 fr)))) obsp F (reverse F))))]
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
  [quantiles strategy mse? data]
  (let [qs (cond
             (sequential? quantiles) quantiles
             (number? quantiles) (uniform-quantilies (long quantiles))
             :else (uniform-quantilies 50))
        measure (if mse? m/sq m/abs)
        data-qs (stats/quantiles data qs strategy)]
    (fn [distr]
      (stats/mean (map (fn [^double q1 ^double q2]
                         (measure (- q1 q2))) data-qs (map (partial r/icdf distr) qs))))))

(defn- mme
  [mse? data]
  (let [xs (m/seq->double-array data)
        mean (stats/mean xs)
        variance (stats/variance xs)
        measure (if mse? m/sq m/abs)]
    (fn [distr]
      (let [tmean (r/mean distr)
            tvariance (r/variance distr)]
        (if (or (m/nan? tmean)
                (m/nan? tvariance)
                (m/inf? tmean)
                (m/inf? tvariance))
          (throw (Exception. "Target distribution can't be used in mme method."))
          (* 0.5 (+ ^double (measure (- mean tmean))
                    ^double (measure (- variance tvariance)))))))))

(defn- mps
  [data]
  (let [c (count data)
        fac (/ (double c) (inc c))
        s (vec (sort (conj data ##-Inf ##Inf)))]
    (fn [distr]
      (let [starget (map (partial r/cdf distr) s)
            diffs (map (fn [^double d1 ^double d2 ^double d3]
                         (let [diff (- d2 d1)]
                           (if (< diff 1.0e-200)
                             (r/lpdf distr d3)
                             (m/log diff)))) starget (rest starget) data)]
        (* fac (stats/mean diffs))))))

;;

(defn- method->fn
  [method data params]
  (case method 
    :ks (kolmogorov-smirnov data)
    :cvm (cramer-von-mises data)
    :ad (anderson-darling data)
    :adr (anderson-darling-r data)
    :adl (anderson-darling-l data)
    :ad2r (anderson-darling-2r data)
    :ad2l (anderson-darling-2l data)
    :ad2 (anderson-darling-2 data)
    :qme (let [{:keys [quantiles strategy mse?]
                :or {quantiles 50 strategy :legacy mse? true}} params] (qme quantiles strategy mse? data))
    :mle (fn [distr] (r/log-likelihood distr data))
    :mme (mme (get params :mse? true) data)
    :mps (mps data)
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
  (if (#{:mle :mps} method) o/maximize o/minimize))

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
                                    :or {optimizer :lbfgsb}
                                    :as all}]
  (let [{:keys [param-names bounds inference]} (distribution-data distribution)
        opt-fn (method->opt-fn method) ;; minimize or maximize?
        [pars result] (opt-fn optimizer target (merge {:initial (inference data)
                                                       :gradient-h 1.0e-5
                                                       :max-iters 1000} all {:bounds bounds 
                                                                             :bounded? true
                                                                             :stats? false})) ;; optimize!
        
        conf (zipmap param-names pars) ;; create final distribution parametrization...
        distr (r/distribution distribution conf)] ;; ...and distribution
    (merge (calc-stats method distr data all result (count param-names)) ;; return result and statistics
           {:params conf
            :distribution-name distribution
            :distribution distr
            :method method})))

(defn- make-target
  [method distribution data all param-names]
  (let [raw-target (method->fn method data all)]
    (fn [& r]
      (let [d (r/distribution distribution (zipmap param-names r))]
        (try
          (raw-target d)
          (catch Exception e
            (.printStackTrace e)
            (throw e)))))))

(defn fit
  "Fit distribution using given method

  * `:mle` - log likelihood
  * `:ad` - Anderson-Darling
  * `:adr`, `:adl`, `:ad2r`, `:ad2l` and `:ad2` - Anderson-Darling variants
  * `:ks` - Kolmogorov-Smirnov
  * `:cvm` - Cramer-von-Mises
  * `:qme` - quantile matching estimation.
  * `:mme` - method of moments (modified)
  * `:mps` - maximum product of spacing estimation

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
     (when assert? (assert (validation data) "Data values do not fit required distribution"))
     (when assert? (assert (not (#{:frechet :johnson-sb} distribution))
                           (str distribution " returns wrong pararameters (use fit instead of infer)")))
     (let [conf (zipmap param-names (inference data))
           distr (r/distribution distribution conf)]
       (merge (calc-stats nil distr data (update params :stats conj :mle) nil (count param-names))
              {:params conf
               :method :infer
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
    :percentile-bc stats/percentile-bc-extent
    :percentile-bca stats/percentile-bca-extent
    :adj-median stats/adjacent-values
    :ci stats/ci
    :min-max-mean stats/extent
    stats/mad-extent))

(defn bootstrap 
  "Bootstrapped version of fitting.

  Parameters are the same as [[fit]]
  
  Additional parameters:

  * `:samples` - number of bootstrapped sequences (default: 100)
  * `:size` - number of samples in each sequence (default: 10% of data, minimum 100, maximum 5000 samples)
  * `:ci-type` - confidence interval type (default: `:mad-median`)
  * `:all-params?` - fitted parameters for each sequence (default: false)"
  ([method distribution data] (bootstrap method distribution data {}))
  ([method distribution data {:keys [size samples ci-type all-params? assert?]
                              :or {samples 100
                                   size (m/constrain (* 0.1 (count data)) 100 5000)
                                   ci-type :mad-median
                                   all-params? false
                                   assert? true}
                              :as all}]
   (let [{:keys [param-names bounds validation]} (distribution-data distribution)]
     (when assert? (assert-values data bounds validation all))
     (let [bdata (:samples (boot/bootstrap data nil {:samples samples :size size}))
           res (pmap #(fit method distribution % (-> all
                                                     (dissoc :stats)
                                                     (assoc :assert? false))) bdata) ;; fit paralelly (do not calculate stats)
           all-params (map vals (map :params res)) ;; extract found parameters
           extents (map (ci->fn ci-type) (apply map vector all-params)) ;; calculate extents for each parameter
           params (map last extents) ;; extract mean or median
           conf (zipmap param-names params) ;; create distribution configuration
           distr (r/distribution distribution conf) ;; create distribution
           res (merge (calc-stats nil distr data (update all :stats conj method) nil (count param-names))
                      {:ci (zipmap param-names extents)
                       :ci-type ci-type
                       :params conf
                       :distribution-name distribution
                       :distribution distr
                       :method method
                       :bootstrap? true})]
       (if all-params? ;; maybe you want full list for each resampled data
         (assoc res :all-params all-params)
         res)))))

;;

(def distribution r/distribution)

(defn ->distribution
  "Return distribution object"
  [distr]
  (cond
    (r/distribution? distr) distr
    (and (map? distr)
         (:distribution distr)) (:distribution distr)
    :else (throw (Exception. (str "Not a distribution: " distr)))))

(defn cdf
  "Cumulative probability."
  (^double [d v] (prot/cdf (->distribution d) v))
  (^double [d v1 v2] (prot/cdf (->distribution d) v1 v2)))

(defn pdf
  "Density"
  ^double [d v] (prot/pdf (->distribution d) v))

(defn lpdf
  "Log density"
  ^double [d v] (prot/lpdf (->distribution d) v))

(defn icdf
  "Inverse cumulative probability"
  [d ^double v] (prot/icdf (->distribution d) v))

(defn probability
  "Probability (PMF)"
  ^double [d v] (prot/probability (->distribution d) v))

(defn sample
  "Random sample"
  [d] (prot/sample (->distribution d)))

(defn observe1
  "Log of probability/density of the value. Alias for [[lpdf]]."
  ^double [d v]
  (prot/lpdf (->distribution d) v))

(defn log-likelihood
  "Log likelihood of samples"
  ^double [d vs]
  (let [d (->distribution d)]
    (reduce (fn [^double s ^double v] (if (m/invalid-double? s)
                                       (reduced s)
                                       (+ s v))) 0.0 (map #(prot/lpdf d %) vs))))

(defmacro observe
  "Log likelihood of samples. Alias for [[log-likelihood]]."
  [d vs]
  `(log-likelihood ~d ~vs))

(defn likelihood
  "Likelihood of samples"
  ^double [d vs]
  (m/exp (log-likelihood d vs)))

(defn mean
  "Distribution mean"
  ^double [d] (prot/mean (->distribution d)))

(defn variance
  "Distribution variance"
  ^double [d] (prot/variance (->distribution d)))

(defn lower-bound
  "Distribution lowest supported value"
  ^double [d] (prot/lower-bound (->distribution d)))

(defn upper-bound
  "Distribution highest supported value"
  ^double [d] (prot/upper-bound (->distribution d)))

(defn distribution-id
  "Distribution identifier as keyword."
  [d] (prot/distribution-id (->distribution d)))

(defn distribution-parameters
  "Distribution highest supported value.
  When `all?` is true, technical parameters are included, ie: `:rng` and `:inverser-cumm-accuracy`."
  ([d] (distribution-parameters d false))
  ([d all?]
   (if-not all?
     (-> (prot/distribution-parameters (->distribution d))
         (set)
         (disj :rng :inverse-cumm-accuracy)
         (vec))
     (prot/distribution-parameters (->distribution d)))))

(defn drandom
  "Return random double"
  ^double [d]
  (prot/drandom (->distribution d)))

(defn lrandom
  "Return random long"
  ^long [d]
  (prot/lrandom (->distribution d)))

(defn irandom
  "Return random integer"
  ^long [d]
  (prot/irandom (->distribution d)))

(defn ->seq
  "Return sequence of random values"
  ([d] (prot/->seq (->distribution d)))
  ([d n] (prot/->seq (->distribution d) n)))

(defn set-seed!
  "Set seed of the distribution, returns distribution object"
  [d seed]
  (prot/set-seed! (->distribution d) seed))

(m/unuse-primitive-operators)


;;;;;;; ignore

(comment
  (def target (r/->seq (r/distribution :cauchy {:median 2}) 1000))
  (take 10 target)
  (infer :rayleigh target)
  (time (fit :mps :normal target {:stats [:mle]
                                  :mse? false})))

(comment
  (def d (r/distribution :half-normal {:sigma 1.01}))
  (def samples (r/->seq d 100))
  (stats/maximum samples)
  (stats/minimum samples)
  (stats/mode samples :histogram {:bins (/ (count samples) 10)})
  (stats/mean samples)
  (r/variance d)
  (r/mean d)
  (infer :half-normal samples)
  (fit :mle :half-normal samples))

(comment
  (def in [21 35 55 100 134 13 15 17 16 30 34 45 60 34 55])
  (def in [16.5 48.5 2.5
         41.25 9.5 27.75 17 2.5 8.25 9.75 35.25 8.5
         3 14.5 11.5 12.25 7.5 2 25 17.5 33 5 9
         4 27 29 19.5 4.75 3 32.25 10.5 3.5 15
         6.75 15.25 16 2.25 1 1 1.5 4
         6.75 9.5 7 7.5 6.5 9 25.25 3
         1 1.5 1.25 4
         3
         3.25
         16.75
         6
         11.25
         4
         7.75
         16.25
         15.25
         2.75
         1.25
         15.5
         4.75
         16.75
         7.25
         7.25
         10.75
         6.5
         49.25
         10
         5.25
         8.5
         1])
  (def in [24
         16
         12
         12
         32
         12
         8
         12
         4
         8
         4
         12
         16
         8
         4
         16
         16
         12
         8
         8
         16
         4
         12
         8
         8
         2
         24
         8
         12
         4
         4
         32
         8
         4
         16
         1
         12
         16
         3
         2
         1
         1
         1
         2
         4
         8
         8
         24
         4
         8
         24
         1
         1
         1
         4
         2
         4
         8
         8
         12
         20
         4
         8
         40
         16
         16
         8
         5
         1
         1
         16
         8
         8
         8
         8
         12
         8
         24
         4
         12
         16
         8
         20
         12
         4
         4
         8
         1
         5
         16
         2
         8
         6
         2
         8
         16
         8
         16
         12
         3
         8
         12])
  (def in [13.0
         2.0
         5.0
         8.0
         13.0
         3.0
         3.0
         3.0
         3.0
         3.0
         13.0
         3.0
         3.0
         5.0
         3.0
         3.0
         3.0
         3.0
         3.0
         2.0
         5.0
         3.0
         3.0
         5.0
         5.0
         3.0
         3.0
         3.0
         2.0
         3.0
         8.0
         13.0
         5.0
         8.0
         8.0
         1.0
         13.0
         13.0
         3.0]
    (def in [3,
           4,
           3,
           13,
           6,
           30,
           29,
           2,
           12,
           40,
           16,
           382,
           180,
           36,
           6,
           49,
           36,
           13,
           55,
           12,
           20]))
  (stats/mean (map #(m/log %) in))
  (stats/stddev (map #(m/log %) in))
  (defn find-best
    ([method ds] (find-best method ds in))
    ([method ds data]
     (let [selector (if (= method :mle) last first)]
       (->> (map #(fit method % data {:stats #{:mle :ad :ks :cvm}}) ds)
            (remove (fn [v] (Double/isNaN (method (:stats v)))))
            (sort-by (comp method :stats))
            (selector)))))
  (def best (find-best :mle [:poisson :log-normal :gamma]))
  (def best-distr (:distribution best))
  (def means (repeatedly 2000 #(stats/mean (r/->seq best-distr 55))))
  (def sums-distr (r/distribution :real-discrete-distribution {:data means}))
  )
