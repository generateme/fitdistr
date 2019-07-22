(ns gen-images
  (:require [cljplot.core :refer :all]
            [cljplot.build :as b]
            [fitdistr.core :as f]
            [fastmath.random :as r]
            [fastmath.stats :as stats]))

;; example 2

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
    (dissoc (->> (map #(f/fit method % atv {:stats #{:mle :ad :ks :cvm}}) ds)
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

(def pal [0x019FD7,0x9ECC3C,0xFFDB1A,0xCC3C3C,0x9D9D9D,:red])

(def f-log-normal (:distribution (f/fit :mle :log-normal atv)))
(def f-weibull (:distribution (f/fit :mle :weibull atv)))
(def f-gamma (:distribution (f/fit :mle :gamma atv)))
(def f-pareto (:distribution (f/fit :mle :pareto atv)))

(def cdf-conf {:domain (stats/extent atv)
               :margins {:x [0.0 0.01]
                         :y [0.0 0.01]}
               :samples 10000
               :stroke {:dash [10 5]}})

(def legend (map #(vector :line (name %1)
                          (if (= %2 5)
                            {:color (nth pal %2)}
                            {:color (nth pal %2)
                             :stroke (:stroke cdf-conf)}))
                 [:empirical :log-normal :weibull :gamma :pareto]
                 [5 0 1 3 4]))

(save (-> (xy-chart {:width 600 :height 500}
                    (b/series [:grid]
                              [:cdf atv {:margins {:x [0.0 0.01]
                                                   :y [0.0 0.01]}
                                         :samples 10000
                                         :color :red}]
                              [:function (partial r/cdf f-log-normal) (assoc cdf-conf :color (nth pal 0))]
                              [:function (partial r/cdf f-weibull) (assoc cdf-conf :color (nth pal 1))]
                              [:function (partial r/cdf f-gamma) (assoc cdf-conf :color (nth pal 3))]
                              [:function (partial r/cdf f-pareto) (assoc cdf-conf :color (nth pal 4))])
                    (b/update-scale :x :scale [:log])
                    (b/update-scale :y :scale [:log])
                    (b/update-scale :y :domain [0.15 1.0])
                    (b/update-scale :y :ticks [0.01 0.05 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0])
                    (b/add-legend "Series" legend)
                    (b/add-label :top "Maximum log likelihood, CDF comparison")
                    (b/add-label :bottom "log(data)")
                    (b/add-label :left "log(cdf)")
                    (b/add-axes :bottom)
                    (b/add-axes :left)))
      "utils/ex2-cdf.jpg")

;; example 3

(def gb [30.0 10.0 20.0 24.0 20.0 24.0 40.0 20.0 50.0 30.0 26.0 44.0 25.0 20.0 40.0 100.0 30.0 80.0 50.0 40.0 20.0 50.0 50.0 50.0
         50.0 80.0 80.0 20.0 60.0 50.0 50.0 50.0 50.0 60.0 50.0 50.0 50.0 50.0 17.0 120.0 80.0 100.0 80.0 100.0 50.0 70.0 50.0 50.0
         50.0 30.0 90.0 125.0 50.0 25.0 50.0 75.0 100.0 60.0 100.0 11.5 100.0 100.0 50.0 40.0 50.0 50.0 50.0 150.0 50.0 80.0 100.0 100.0
         80.0 20.0 30.0 100.0 70.0 100.0 80.0 80.0 50.0 50.0 50.0 70.0 25.0 40.0 100.0 25.0 130.0 100.0 60.0 60.0 100.0 50.0 50.0 50.0
         50.0 50.0 50.0 100.0 30.0 50.0 50.0 50.0 100.0 50.0 50.0 75.0 25.0 20.0 75.0 50.0 50.0 100.0 100.0 80.0 130.0 200.0 40.0 64.0
         40.0 80.0 80.0 79.0 80.0 105.0 51.0 130.0 105.0 105.0 80.0 80.0 40.0 40.0 50.0 105.0 105.0 80.0 105.0 40.0 50.0 130.0 105.0 80.0
         40.0 80.0 50.0 100.0 200.0 79.0 150.0 80.0 80.0 80.0 80.0 80.0 65.0 130.0 40.0 80.0 80.0 80.0 130.0 51.0 130.0 80.0 80.0 130.0
         51.0 25.0 79.0 80.0 50.0 51.0 50.0 105.0 80.0 80.0 130.0 130.0 80.0 51.0 50.0 130.0 115.0 115.0 200.0 160.0 105.0 80.0 80.0 80.0
         130.0 130.0 80.0 51.0 50.0 130.0 105.0 105.0 80.0 105.0 80.0 60.0 51.0 105.0 80.0 51.0 40.0 130.0 130.0 40.0 32.0 80.0 25.0 80.0
         80.0 80.0 40.0 130.0 130.0 40.0 75.0 40.0 80.0 80.0 150.0 52.5 130.0 80.0 80.0 75.0 40.0 25.5 51.0 50.0 80.0 51.0 130.0 130.0
         51.0 130.0 80.0 130.0 80.0 130.0 130.0 25.5 117.0 130.0 80.0 80.0 80.0 80.0])

(f/fit :qme :gamma gb {:stats [:mle]})
;; => {:stats
;;     {:qme 45.41350132338318,
;;      :mle -1253.652887342403,
;;      :aic 2511.305774684806,
;;      :bic 2518.3804432188426},
;;     :params {:scale 18.58013425207939, :shape 3.992585123164037},
;;     :distribution-name :gamma,
;;     :method :qme}

(f/fit :qme :weibull gb {:stats [:mle]})
;; => {:stats
;;     {:qme 52.75300844658587,
;;      :mle -1256.1208259372665,
;;      :aic 2516.241651874533,
;;      :bic 2523.3163204085704},
;;     :params {:alpha 2.052173279173565, :beta 83.1667239954531},
;;     :distribution-name :weibull,
;;     :method :qme}

(f/fit :qme :log-normal gb {:stats [:mle]})
;; => {:stats
;;     {:qme 58.022608223760685,
;;      :mle -1267.6351159421142,
;;      :aic 2539.2702318842285,
;;      :bic 2546.3449004182658},
;;     :params {:scale 4.204520160089395, :shape 0.46585471421274693},
;;     :distribution-name :log-normal,
;;     :method :qme}

(save (-> (xy-chart {:width 600 :height 400}
                    (b/series [:grid]
                              [:histogram gb {:bins 10
                                              :density? true}]
                              [:density gb {:color [0 0 0 50]
                                            :area? true}])
                    (b/add-axes :left)
                    (b/add-axes :bottom)
                    (b/add-label :top "Empirical density and histogram")
                    (b/add-label :left "Density")
                    (b/add-label :bottom "groundbeef servings")))
      "utils/ex3-hist.jpg")

(save (-> (xy-chart {:width 600 :height 400}
                    (b/series [:grid]
                              [:histogram gb {:bins 100
                                              :cumulative? true
                                              :density? true
                                              :type :lollipops}])
                    (b/add-axes :left)
                    (b/add-axes :bottom)
                    (b/add-label :top "Cumulative histogram")
                    (b/add-label :left "CDF")
                    (b/add-label :bottom "groundbeef servings")))
      "utils/ex3-cumulative.jpg")

(def ex2-gamma (f/fit :qme :gamma gb {:stats [:mle]}))
(def ex2-weibull (f/fit :qme :weibull gb {:stats [:mle]}))
(def ex2-log-normal (f/fit :qme :log-normal gb {:stats [:mle]}))
;; => {:stats
;;     {:qme 45.41350132338318,
;;      :mle -1253.652887342403,
;;      :aic 2511.305774684806,
;;      :bic 2518.3804432188426},
;;     :params {:scale 18.58013425207939, :shape 3.992585123164037},
;;     :distribution-name :gamma,
;;     :method :qme}

(def pdf-conf {:domain (stats/extent gb)
               :samples 500
               :points 150
               :size 8
               :stroke {:dash [10 5]}})

(def legend2 (map #(vector :line (name %1) {:color (nth pal %2)
                                            :stroke (:stroke pdf-conf)})
                  [:log-normal :weibull :gamma]
                  [3 1 0]))

(def legend2-c (map #(vector :circle (name %1) {:color (nth pal %2)
                                                :stroke (:stroke pdf-conf)})
                    [:log-normal :weibull :gamma]
                    [3 1 0]))

(save (-> (xy-chart {:width 600 :height 500}
                    (b/series [:grid]
                              [:histogram gb {:bins 10
                                              :density? true}]
                              [:function (partial r/pdf (:distribution ex2-gamma)) (assoc pdf-conf :color (nth pal 0))]
                              [:function (partial r/pdf (:distribution ex2-weibull)) (assoc pdf-conf :color (nth pal 1))]
                              [:function (partial r/pdf (:distribution ex2-log-normal)) (assoc pdf-conf :color (nth pal 3))])
                    (b/add-axes :left)
                    (b/add-axes :bottom)
                    (b/add-legend "PDF" legend2)
                    (b/add-label :top "Distribution densities and histogram")
                    (b/add-label :left "Density")
                    (b/add-label :bottom "groundbeef servings")))
      "utils/ex3-pdf.jpg")


(save (-> (xy-chart {:width 600 :height 400}
                    (b/series [:grid]
                              [:abline]
                              [:ppplot [(:distribution ex2-gamma) gb] (assoc pdf-conf :color (nth pal 0))]
                              [:ppplot [(:distribution ex2-weibull) gb] (assoc pdf-conf :color (nth pal 1))]
                              [:ppplot [(:distribution ex2-log-normal) gb] (assoc pdf-conf :color (nth pal 3))])
                    (b/add-axes :left)
                    (b/add-axes :bottom)
                    (b/add-legend "Distr" legend2-c)
                    (b/add-label :top "P-P plots")
                    (b/add-label :left "Empirical probabilities")
                    (b/add-label :bottom "Theoretical probabilities")))
      "utils/ex3-pp.jpg")

(save (-> (xy-chart {:width 600 :height 400}
                    (b/series [:grid]
                              [:abline]
                              [:qqplot [(:distribution ex2-gamma) gb] (assoc pdf-conf :color (nth pal 0))]
                              [:qqplot [(:distribution ex2-weibull) gb] (assoc pdf-conf :color (nth pal 1))]
                              [:qqplot [(:distribution ex2-log-normal) gb] (assoc pdf-conf :color (nth pal 3))])
                    (b/add-axes :left)
                    (b/add-axes :bottom)
                    (b/add-legend "Distr" legend2-c)
                    (b/add-label :top "Q-Q plots")
                    (b/add-label :left "Empirical quantiles")
                    (b/add-label :bottom "Theoretical quantiles")))
      "utils/ex3-qq.jpg")


(def b-gamma (f/bootstrap :qme :gamma gb {:all-params? true
                                          :ci-type :min-max-mean
                                          :samples 10000
                                          :size 100}))

(count (:all-params b-gamma))
;; => 10000

(:ci b-gamma)
;; => {:scale [9.436000626775108 28.586945660073017 17.967024829900268],
;;     :shape [2.5005015188673814 7.565181039020107 4.238506386828045]}

(defn show-scatter
  [b]
  (-> (xy-chart {:width 600 :height 600}
                (b/series [:grid]
                          [:scatter (:all-params b)]
                          [:density-2d (:all-params b) {:margins {:y [0 0] :x [0 0]}
                                                        :blur-kernel-size 0.1}]
                          [:scatter [((juxt :scale :shape) (:params b))] {:color :red
                                                                          :size 10}])
                (b/add-axes :left)
                (b/add-axes :bottom)
                (b/add-label :top "Bootstrapped parameters")
                (b/add-label :left "scale")
                (b/add-label :bottom "shape"))))

(save (show-scatter b-gamma)
      "utils/ex4-gamma.jpg")


(def b-gamma-lowq (f/bootstrap :qme :gamma gb {:all-params? true
                                               :quantiles 3
                                               :samples 10000
                                               :size 100}))

(save (show-scatter b-gamma) "utils/ex4-gamma-lowq.jpg")

;;

(show (xy-chart {:width 600 :height 600}
                (b/series [:grid]
                          [:function (partial r/pdf (r/distribution :watson-g {:n 40})) {:samples 600
                                                                                         :domain [-2 2]}]                          )
                (b/add-axes :left)
                (b/add-axes :bottom)
                (b/add-label :top "PDF")))

(show (xy-chart {:width 600 :height 600}
                (b/series [:grid]
                          [:histogram (r/->seq (r/distribution :levy {:mu 10 :c 10}) 1000) {:bins 10}])
                (b/add-axes :left)
                (b/add-axes :bottom)
                (b/add-label :top "PDF")))

;;
