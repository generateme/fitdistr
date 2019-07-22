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
  (m/seq->double-array (map #(m/log ^double %) data)))

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
  (let [xs (m/seq->double-array data)
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
  (let [xs (m/seq->double-array data)
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
  (let [xs (m/seq->double-array data)
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
  (let [xs (m/seq->double-array data)
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
  (let [xs (m/seq->double-array data)
        m1 (stats/mean xs)
        x1 (- (stats/minimum xs) 0.01)
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
   :inference #(let [xs (m/seq->double-array %)]
                 [(stats/median xs) (* 0.5 (stats/iqr xs))])})

(defmethod distribution-data :chi-squared [_]
  {:param-names [:degrees-of-freedom]
   :bounds [b0+]
   :validation non-negatives
   :inference (comp vector stats/mean)})

(defn- infer-f
  [data]
  (let [xs (m/seq->double-array data)
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
  (let [xs (m/seq->double-array data)
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
  (let [xs (m/seq->double-array data)
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
  (let [vs (m/seq->double-array data)
        m (- (stats/minimum vs) 0.01)
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
  (let [vs (m/seq->double-array data)
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
  (let [vx (m/seq->double-array data)
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
  (let [vx (m/seq->double-array data)
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

(defn- infer-inverse-gamma
  [data]
  (let [m1 (stats/mean data)
        m2 (stats/mean (map m/sq data))
        m12 (m/sq m1)
        d (- m2 m12)
        shape (/ (- (* 2.0 m2) m12) d)
        scale (/ (* m1 m2) d)]
    [shape scale]))

(defmethod distribution-data :inverse-gamma [_]
  {:param-names [:alpha :beta]
   :bounds [b0+ b0+]
   :validation positives
   :inference infer-inverse-gamma})

(defonce ^:private ^:const ^double ll-const (* 2.0 (m/log 3.0)))

(defn- infer-log-logistic
  [data]
  (let [xs (m/seq->double-array data)
        beta (stats/median xs)
        [p25 p75] (stats/quantiles xs [0.25 0.75])
        shape (/ ll-const (- (m/log p75) (m/log p25)))]
    [shape beta]))

(defmethod distribution-data :log-logistic [_]
  {:param-names [:alpha :beta]
   :bounds [b0+ b0+]
   :validation positives
   :inference infer-log-logistic})

(defn- infer-chi
  [data]
  (let [xs (m/seq->double-array data)
        m2 (m/sq (stats/mean xs))
        sd (stats/population-stddev xs)]
    [(+ sd m2)]))

(defmethod distribution-data :chi [_]
  {:param-names [:nu]
   :bounds [b1+]
   :validation positives
   :inference infer-chi})

(defn- infer-chi-squared-noncentral
  [data]
  (let [xs (m/seq->double-array data)
        m (stats/mean xs)
        v (stats/population-variance xs)
        l (* 0.5 (- v m m))]
    [(- m l) l]))

(defmethod distribution-data :chi-squared-noncentral [_]
  {:param-names [:nu :lambda]
   :bounds [b0+ b0+]
   :validation positives
   :inference infer-chi-squared-noncentral})

(defn- infer-erlang
  [data]
  (let [xs (m/seq->double-array data)
        m (stats/mean xs)
        v (stats/population-variance xs)
        l (/ m v)]
    [(* m l) l]))

(defmethod distribution-data :erlang [_]
  {:param-names [:k :lambda]
   :bounds [b1+ b0+]
   :validation positives
   :inference infer-erlang})

(defn- infer-fatigue-life
  [data]
  (let [xs (m/seq->double-array data)
        mu (- (stats/minimum xs) 0.01)
        m (stats/mean xs)
        v (stats/population-variance xs)
        loc2 (m/sq (- m mu))
        a (* 0.25 (- v (* 5.0 loc2)))
        b (- v loc2)
        c v
        delta (- (* b b) (* 4.0 a c))
        gamma2 (/ (- (- b) (m/sqrt delta))
                  (+ a a))] 
    [mu (/ (- m mu)
           (inc (* 0.5 gamma2))) (m/sqrt gamma2)]))

(defmethod distribution-data :fatigue-life [_]
  {:param-names [:mu :beta :gamma]
   :bounds [bfull b0+ b0+]
   :validation all-accepted
   :inference infer-fatigue-life})

(defn- infer-frechet
  [data]
  (let [xs (m/seq->double-array data)
        mu (- (stats/minimum xs) 0.01)]
    (conj (vec (umontreal.ssj.probdist.FrechetDist/getMLE xs (alength xs) mu)) mu)))

(defmethod distribution-data :frechet [_]
  {:param-names [:alpha :beta :delta]
   :bounds [b0+ b0+ bfull]
   :validation all-accepted
   :inference infer-frechet})

(def ^:private infer-hyperbolic-secant infer-normal)

(defmethod distribution-data :hyperbolic-secant [_]
  {:param-names [:mu :sigma]
   :bounds [bfull b0+]
   :validation all-accepted
   :inference infer-hyperbolic-secant})

(defn- infer-inverse-gaussian
  [data]
  (let [m (stats/mean data)
        rm (/ m)
        s (/ (stats/mean (map (fn [^double x] (- (/ x) rm)) data)))]
    [m s]))

(defmethod distribution-data :inverse-gaussian [_]
  {:param-names [:mu :lambda]
   :bounds [b0+ b0+]
   :validation positives
   :inference infer-inverse-gaussian})

(defn- infer-johnson-su
  [data]
  (let [xs (m/seq->double-array data)
        m (stats/mean xs)
        s (stats/population-stddev xs)]
    [0.0 1.0 m (/ s 1.7)]))

(defmethod distribution-data :johnson-su [_]
  {:param-names [:gamma :delta :xi :lambda]
   :bounds [bfull b0+ bfull b0+]
   :validation all-accepted
   :inference infer-johnson-su})

(defn- infer-johnson-sb
  [data]
  (let [xs (m/seq->double-array data)
        [^double mn ^double mx] (stats/extent xs)
        xi (- mn 0.01)
        lambda (+ 0.01 (- mx xi))
        [gamma delta] (umontreal.ssj.probdist.JohnsonSBDist/getMLE xs (alength xs) xi lambda)]
    [gamma delta xi lambda]))

(defmethod distribution-data :johnson-sb [_]
  {:param-names [:gamma :delta :xi :lambda]
   :bounds [bfull b0+ bfull b0+]
   :validation all-accepted
   :inference infer-johnson-sb})

(defn- infer-johnson-sl
  [data]
  (vec (umontreal.ssj.probdist.JohnsonSLDist/getMLE (m/seq->double-array data) (count data))))

(defmethod distribution-data :johnson-sl [_]
  {:param-names [:gamma :delta :xi :lambda]
   :bounds [bfull b0+ bfull b0+]
   :validation all-accepted
   :inference infer-johnson-sl})

(defn- infer-pearson-6
  [data]
  (vec (umontreal.ssj.probdist.Pearson6Dist/getMLE (m/seq->double-array data) (count data))))

(defmethod distribution-data :pearson-6 [_]
  {:param-names [:alpha1 :alpha2 :beta]
   :bounds [b0+ b1+ b0+]
   :validation positives
   :inference infer-pearson-6})

(defn- infer-power
  [data]
  (let [xs (m/seq->double-array data)
        [^double mn ^double mx] (stats/extent xs)
        mn- (- mn 0.01)
        mx+ (+ mx 0.01)]
    (concat [mn- mx+] (umontreal.ssj.probdist.PowerDist/getMLE xs (alength xs) mn mx))))

(defmethod distribution-data :power [_]
  {:param-names [:a :b :c]
   :bounds [bfull bfull bfull]
   :validation all-accepted
   :inference infer-power})

(defn- infer-triangular
  [data]
  (let [xs (m/seq->double-array data)
        [^double mn ^double mx] (stats/extent xs)
        mn- (- mn 0.01)
        mx+ (+ mx 0.01)
        m (stats/mean xs)]
    [mn- (- (* 3.0 m) mn mx) mx+]))

(defmethod distribution-data :triangular [_]
  {:param-names [:a :c :b]
   :bounds [bfull bfull bfull]
   :validation all-accepted
   :inference infer-triangular})

(defn- infer-rayleigh
  [data]
  (let [xs (m/seq->double-array data)
        mn (stats/minimum xs)
        mn- (- mn 0.01)]
    (concat [mn-] (umontreal.ssj.probdist.RayleighDist/getMLE xs (alength xs) mn))))

(defmethod distribution-data :rayleigh [_]
  {:param-names [:a :beta]
   :bounds [bfull b0+]
   :validation all-accepted
   :inference infer-rayleigh})
