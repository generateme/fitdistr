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
(defonce ^:private ^:const b05+ [0.5 ##Inf])
(defonce ^:private ^:const b1+ [1.0 ##Inf])
(defonce ^:private ^:const bfull [##-Inf ##Inf])

(def ^:private positives (partial every? (fn [^double x] (pos? x))))
(def ^:private non-negatives (partial every? (fn [^double x] (not (neg? x)))))
(def ^:private zero-one (partial every? (fn [^double x] (< 0.0 x 1.0))))
(def ^:private zero-or-one (partial every? (fn [^double x] (bool-or (zero? x)
                                                                   (== x 1.0)))))
(def ^:private all-accepted (constantly true))

(defn- log-data
  [data]
  (double-array (map #(m/log ^double %) data)))

(defn- filter-outliers
  [data]
  (let [q1 (stats/percentile data 25)
        q3 (stats/percentile data 75)
        iqr15 (* 1.5 (- q3 q1))
        l (- q1 iqr15)
        u (+ q3 iqr15)]
    (filter (fn [^double v] (< l v u)) data)))

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
        s (stats/population-stddev lx)
        shape (/ 1.28 s)]
    [shape (m/exp (+ m (/ m/GAMMA shape)))]))

(defmethod distribution-data :weibull [_]
  {:param-names [:alpha :beta]
   :bounds [b0+ b0+]
   :validation positives
   :inference infer-weibull})

(defn- infer-normal
  [data]
  (let [xs (double-array data)
        m (stats/mean xs)
        s (stats/population-stddev xs)]
    [m s]))

(defmethod distribution-data :normal [_]
  {:param-names [:mu :sd]
   :bounds [bfull b0+]
   :validation all-accepted
   :inference infer-normal})

(defn- infer-log-normal
  [data]
  (let [lxs (log-data data)
        m (stats/mean lxs)
        s (stats/population-stddev lxs)]
    [m s]))

(defmethod distribution-data :log-normal [_]
  {:param-names [:scale :shape]
   :bounds [b0+ b0+]
   :validation positives
   :inference infer-log-normal})

(defn- infer-negative-binomial
  [data]
  (let [xs (double-array data)
        m (stats/mean xs)
        v (stats/population-variance xs)
        ^double s (if (> v m) (/ (* m m) (- v m)) 100)
        p (/ s (+ s m))]
    [s p]))

(defmethod distribution-data :negative-binomial [_]
  {:param-names [:r :p]
   :bounds [b0+ b01]
   :validation non-negatives
   :inference infer-negative-binomial})

(defmethod distribution-data :pascal [_]
  {:param-names [:r :p]
   :bounds [b0+ b01]
   :validation non-negatives
   :inference infer-negative-binomial})


(defmethod distribution-data :poisson [_]
  {:param-names [:p]
   :bounds [b0+]
   :validation non-negatives
   :inference (comp vector stats/mean)})

(defmethod distribution-data :exponential [_]
  {:param-names [:mean]
   :bounds [b0+]
   :validation positives
   :inference (comp vector stats/mean)})

(defn- infer-gamma
  [data]
  (let [xs (double-array data)
        m (stats/mean xs)
        v (stats/population-variance xs)
        vm (/ v m)]
    [vm (/ m vm)]))

(defmethod distribution-data :gamma [_]
  {:param-names [:scale :shape]
   :bounds [b0+ b0+]
   :validation positives
   :inference infer-gamma})

(defn- infer-geometric
  [data]
  (let [m (stats/mean data)]
    [(if (pos? m) (/ (inc m)) 1.0)]))

(defmethod distribution-data :geometric [_]
  {:param-names [:p]
   :bounds [b01]
   :validation non-negatives
   :inference infer-geometric})

(defn- infer-beta
  [data]
  (let [xs (double-array data)
        m (stats/mean xs)
        v (stats/population-variance xs)
        aux (dec (/ (* m (- 1.0 m)) v))]
    [(* m aux) (* (- 1.0 m) aux)]))

(defmethod distribution-data :beta [_]
  {:param-names [:alpha :beta]
   :bounds [b0+ b0+]
   :validation zero-one
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
   :validation positives
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
   :validation non-negatives
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
   :validation non-negatives
   :inference infer-f})

(defonce ^:private ^:const ^double gumbel-const (+ m/GAMMA (m/ln (m/ln 2))))

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

(defn- infer-laplace
  [data]
  (let [xs (double-array data)
        m (stats/mean data)
        s (stats/population-stddev data)]
    [m (/ s m/SQRT2)]))

(defmethod distribution-data :laplace [_]
  {:param-names [:mu :beta]
   :bounds [bfull b0+]
   :validation all-accepted
   :inference infer-laplace})

(defn- infer-levy
  [data]
  (let [vs (double-array data)
        m (stats/minimum vs)
        mode (continuous-mode (filter-outliers vs))
        c (* 3.0 (- mode m))] 
    [(- m (* 0.05 c)) c]))

(defmethod distribution-data :levy [_]
  {:param-names [:mu :c]
   :bounds [bfull b0+]
   :validation all-accepted
   :inference infer-levy})

(defn- infer-logistic
  [data]
  (let [vs (double-array data)
        m (stats/mean vs)
        v (stats/population-variance vs)]
    [m (/ (m/sqrt (* 3.0 v)) m/PI)]))

(defmethod distribution-data :logistic [_]
  {:param-names [:mu :s]
   :bounds [bfull b0+]
   :validation all-accepted
   :inference infer-logistic})

(defn- infer-nakagami
  [data]
  (let [o (stats/mean (map m/sq data))
        m (continuous-mode data)
        mu (* 0.5 (/ o (- o (* m m))))]
    [(max mu 1.0) o]))

(defmethod distribution-data :nakagami [_]
  {:param-names [:mu :omega]
   :bounds [b05+ b0+]
   :validation positives
   :inference infer-nakagami})

(defn- infer-t
  [data]
  (let [v (stats/population-variance data)]
    [(if (<= v 1.001) 2000 (/ (* (+ v v)) (dec v)))]))

(defmethod distribution-data :t [_]
  {:param-names [:degrees-of-freedom]
   :bounds [b0+]
   :validation all-accepted
   :inference infer-t})

(defn- infer-binomial
  [data]
  (let [vx (double-array data)
        m (stats/mean vx)
        v (stats/population-variance vx)
        p (if (zero? m)
            0.0
            (/ (m/abs (- m v)) m))]
    [p (if (zero? p) 20
           (/ m p))]))

(defmethod distribution-data :binomial [_]
  {:param-names [:p :trials]
   :bounds [b01 b1+]
   :validation non-negatives
   :inference infer-binomial})

(defn- infer-bernoulli
  [data]
  (let [vx (double-array data)
        m (stats/mean vx)
        v (stats/population-variance vx)
        p (if (zero? m)
            0.0
            (/ (m/abs (- m v)) m))]
    [p (if (zero? p) 20
           (/ m p))]))

(defmethod distribution-data :bernoulli [_]
  {:param-names [:p]
   :bounds [b01]
   :validation zero-or-one
   :inference (comp vector stats/mean)})


#_(def target (r/->seq (r/distribution :bernoulli {:p 0.2 :trials 100}) 10000))
#_(infer-binomial target)
