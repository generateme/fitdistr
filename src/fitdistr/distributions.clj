(ns fitdistr.distributions
  "Distributions information necessary to infer parameters."
  (:require [fastmath.core :as m]
            [fastmath.random :as r]
            [fastmath.stats :as stats]))

(set! *warn-on-reflection* true)
(set! *unchecked-math* :warn-on-boxed)
(m/use-primitive-operators)

(defonce ^:private ^:const b0+ [1.0e-6 ##Inf])
(defonce ^:private ^:const b01 [0.0 1.0])
(defonce ^:private ^:const bfull [##-Inf ##Inf])

(def ^:private positives? (partial every? (fn [^double x] (pos? x))))
(def ^:private non-negatives? (partial every? (fn [^double x] (not (neg? x)))))
(def ^:private zero-one? (partial every? (fn [^double x] (< 0.0 x 1.0))))
(def ^:private all-accepted (constantly true))

(defn- log-data
  [data]
  (double-array (map #(m/log ^double %) data)))

(defn- continuous-mode
  ^double [data]
  (let [{:keys [^double step bins]} (stats/histogram data :sqrt)]
    (+ (* 0.5 step) ^double (ffirst (sort-by second (fn [^double a ^double b] (> a b)) bins)))))

(defmulti distribution-data (fn [k] k))

;; fitdistr.R
(defn- infer-weibull
  [data]
  (let [lx (log-data data)
        m (stats/mean lx)
        s (stats/stddev lx)
        shape (/ 1.28 s)]
    [shape (m/exp (+ m (/ m/GAMMA shape)))]))

(defmethod distribution-data :weibull [_]
  {:param-names [:alpha :beta]
   :bounds [b0+ b0+]
   :validation positives?
   :inference infer-weibull})

(defn- infer-normal
  [data]
  (let [xs (double-array data)
        m (stats/mean xs)
        s (stats/stddev xs)
        n (double (count data))
        f (m/sqrt (/ (dec n) n))]
    [m (* f s)]))

(defmethod distribution-data :normal [_]
  {:param-names [:mu :sd]
   :bounds [bfull b0+]
   :validation all-accepted
   :inference infer-normal})

(defn- infer-log-normal
  [data]
  (let [lxs (log-data data)
        m (stats/mean lxs)
        s (stats/stddev lxs)
        n (double (count data))
        f (m/sqrt (/ (dec n) n))]
    [m (* f s)]))

(defmethod distribution-data :log-normal [_]
  {:param-names [:scale :shape]
   :bounds [b0+ b0+]
   :validation positives?
   :inference infer-log-normal})

(defn- infer-negative-binomial
  [data]
  (let [xs (double-array data)
        m (stats/mean xs)
        n (double (count data))
        v (* (/ (dec n) n) (stats/variance xs))
        ^double s (if (> v m) (/ (* m m) (- v m)) 100)
        p (/ s (+ s m))]
    [s p]))

(defmethod distribution-data :negative-binomial [_]
  {:param-names [:r :p]
   :bounds [b0+ b01]
   :validation non-negatives?
   :inference infer-negative-binomial})

(defmethod distribution-data :poisson [_]
  {:param-names [:p]
   :bounds [b0+]
   :validation non-negatives?
   :inference (comp vector stats/mean)})

(defmethod distribution-data :exponential [_]
  {:param-names [:mean]
   :bounds [b0+]
   :validation positives?
   :inference (comp vector stats/mean)})

(defn- infer-gamma
  [data]
  (let [xs (double-array data)
        m (stats/mean xs)
        n (double (count data))
        v (* (/ (dec n) n) (stats/variance xs))
        vm (/ v m)]
    [vm (/ m vm)]))

(defmethod distribution-data :gamma [_]
  {:param-names [:scale :shape]
   :bounds [b0+ b0+]
   :validation positives?
   :inference infer-gamma})

(defn- infer-geometric
  [data]
  (let [m (stats/mean data)]
    [(if (pos? m) (/ (inc m)) 1.0)]))

(defmethod distribution-data :geometric [_]
  {:param-names [:p]
   :bounds [b01]
   :validation non-negatives?
   :inference infer-geometric})

(defn- infer-beta
  [data]
  (let [xs (double-array data)
        m (stats/mean xs)
        n (double (count data))
        v (* (/ (dec n) n) (stats/variance xs))
        aux (dec (/ (* m (- 1.0 m)) v))]
    [(* m aux) (* (- 1.0 m) aux)]))

(defmethod distribution-data :beta [_]
  {:param-names [:alpha :beta]
   :bounds [b0+ b0+]
   :validation zero-one?
   :inference infer-beta})

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

(defmethod distribution-data :pareto [_]
  {:param-names [:scale :shape]
   :bounds [b0+ b0+]
   :validation positives?
   :inference infer-pareto})

(defmethod distribution-data :cauchy [_]
  {:param-names [:median :scale]
   :bounds [bfull b0+]
   :validation all-accepted
   :inference #(let [xs (double-array %)]
                 [(stats/median xs) (* 0.5 (stats/iqr xs))])})

(defmethod distribution-data :chi-squared [_]
  {:param-names [:degrees-of-freedom]
   :bounds [b0+]
   :validation non-negatives?
   :inference (comp vector stats/mean)})

(defn- infer-f
  [data]
  (let [xs (double-array data)
        m (stats/mean xs)
        d2 (if (< 1.0 m 100.0) (/ (+ m m) (dec m)) 1.0)
        p (/ d2 (+ d2 2.0))
        mode (continuous-mode xs)
        diff (- mode p)
        d1 (if (neg? diff) (/ (* -2.0 p) diff) 1.0)]
    [d1 d2]))

(defmethod distribution-data :f [_]
  {:param-names [:numerator-degrees-of-freedom :denominator-degrees-of-freedom]
   :bounds [b0+ b0+]
   :validation non-negatives?
   :inference infer-f})

(def ^:private ^:const ^double gumbel-const (+ m/GAMMA (m/ln (m/ln 2))))

(defn- infer-gumbel
  [data]
  (let [xs (double-array data)
        m (stats/mean xs)
        median (stats/median xs)
        beta (/ (- m median) gumbel-const)
        mu (- m (* beta m/GAMMA))]
    [mu beta]))

(defmethod distribution-data :gumbel [_]
  {:param-names [:mu :beta]
   :bounds [bfull b0+]
   :validation all-accepted
   :inference infer-gumbel})