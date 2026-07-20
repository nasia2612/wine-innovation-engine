
Computational Food Science portfolio project
# Wine Innovation Engine



-----Phase 1 — Fermentation digital twin-----

The model

The simulator solves a system of 5 coupled ODEs:

dX/dt  = (μ - k_d) · X                         # viable biomass
dXt/dt = μ · X                                  # total biomass  
dS/dt  = -(β·X / Y_es) - MNT·X                 # sugar consumption
dN/dt  = -μ·X / Y_xn                           # nitrogen depletion
dE/dt  = β·X + MNT·X                           # ethanol production


# ethanol turns out to be almost entirely non-growth-associated 
# sugar is not the limiting nutrient in grape juice. Sugar is 150–300 g/L — hugely abundant. What actually runs out first is nitrogen.

# Coleman 2007 uses an nth-order polynomial kd_prime(T) that rises sharply above 25°C, capturing thermal death at high temperatures.

We studied  also the stuck fermentation that is caused by low YAN <80 mg N/L (it might also cause sluggish fermentation,very slow fermantation).Stuck fermentation might also be caused by high temperature ->high kd

For phase 2 :
Boulton (1979) showed that wine fermentors can be modeled as well-mixed heat balances  an assumption valid for small white wine tanks. 
Miller & Block (2020) demonstrated that large-scale fermentations violate this assumption, exhibiting measurable temperature gradients. 
Coleman et al. (2007) established that these gradients directly alter local Monod kinetics, with kd rising sharply above 25°C. This module implements a 1D discretization of the Boulton heat equation coupled to the Coleman kinetic model, serving as a computationally lightweight alternative to full CFD (Miller 2019a).




#we have gradient because of the size of the tank not because of the cap (because we have white wine we dont have cap)

real tanks are not well-mixed:
Fermentation is exothermic — it generates heat
That heat is generated wherever yeast is eating sugar, which is mostly in the upper and middle liquid
The cooling jacket sits at the base and walls — pulling heat out from below
The top of the tank is open air — effectively insulated



#: the liquid develops a vertical temperature gradient. The top is warmer than the bottom.


lets say we use a Jacketed cylindrical white wine fermentor (we discretized it into N horizontal layers for CFD-lite)

Convection happens inside the liquid wine. During fermentation, yeast creates heat. The hot wine in the middle becomes lighter and floats to the top. When it hits the cold steel wall, it cools down, gets heavier, and sinks to the bottom. This creates a continuous 'mixing loop' called a convection current.

#Limitation: The 1D model captures conductive heat transport only. 

Molecular thermal diffusivity α = 1.33×10⁻⁷ m²/s (derived from  Bogard et al. 2020, Table 1) gives a conduction timescale 
τ = H²/α ≈ 10,000 h , approximately 40× longer than fermentation  duration. An effective diffusivity α_eff ≈ 1.33×10⁻⁴ m²/s is 
therefore used to approximate buoyancy-driven natural convection,  following the Nusselt-number scaling implicit in Bogard et al. 2020. 


At 14°C the simulated fermentation is slow (~500–700 h) and retains
modest residual sugar, consistent with Coleman et al. (2007), whose
measured 15°C fermentation (their Fig. 2a) similarly remained above
70 g/L sugar after 600 h. Cold white-wine fermentations are inherently
sluggish. μ_max was set to 0.05 /h following Coleman's reported value
near 11–15°C (their Fig. 3a).



#PHASE 3 RSM -OPTIMIZATION

What we want :
-final ethanol-> want higher 
-residual sugar->want low (we want it dry and to complete the fermentation)
-fermentation time -> want low (Economic/operational goal: less time in the tank = more production cycle, less risk of contamination)

We have conflicts Warmer finishes faster and drier but can overheat; more sugar means more ethanol but slower, stickier fermentation

we  know the screening factors(T,S0,N0) so we do not need a screening experiment.
The three factors were selected a priori based on the mechanistic structure of the Boulton ODE model (Arrhenius temperature dependence, dual-substrate Monod kinetics), not through empirical screening SO the DoE here serves to map the response surface around conditions that we already know have an effect, not to find out which ones do.

The DoE layer varies temperature, initial sugar (S₀), and initial assimilable nitrogen (N₀), but not initial yeast concentration (inoculum size, X₀).
. In this model the yeast population is not static: biomass grows logistically toward a carrying capacity set by the limiting nutrient (nitrogen, secondarily sugar), not by the starting cell count ,a small inoculum reaches the same population ceiling as a large one, just after a longer lag. Final ethanol is therefore governed by nutrient availability, not inoculum size.
 in a three-factor Box–Behnken optimisation of fruit-wine fermentation (Thuy et al., 2023, Food Sci. Technol), dry yeast concentration (0.15–0.25 g/L) was the least influential of the three factors on ethanol, with no further gain above ~0.2 g/L.
 there is no explicit enzyme state because grape must consists of directly fermentable hexoses (glucose, fructose) and needs no saccharification.
 So, enzymatic capacity is absorbed into βmax and the viable-biomass term X_V —> biomass is the model's proxy for enzyme concentration.

A three-factor Box–Behnken design (pyDOE2) is used over T, S0, N0. Factor levels
are coded (-1, 0, +1) with symmetric spacing so the quadratic terms are not
distorted by unequal step sizes.

Also  All levels are kept inside the validated domain of the Coleman et al. (2007) but The temperature range is deliberately kept away from thetested extremes (11 °C, 35 °C), since all four Coleman validation datasets resulted in stuck or incomplete fermentations.


Each model is evaluated by ANOVA; non-significant terms are dropped (model
reduction) so the retained surface reflects only statistically supported effects.
Model adequacy is judged on both R² and residual behaviour.
### Multi-response optimization (desirability)
Each fitted response is mapped to an individual desirability d(Yᵢ) ∈ [0,1] using
one-sided Derringer–Suich functions (`desirability` in R):

- i will add importance  weight to sugar because we care mostly for dryness/no stuck fermentation

- ethanol → `dMax`
- residual sugar → `dMin`
- fermentation time → `dMin`

The overall desirability is the geometric mean

    D = (d_ethanol · d_sugar · d_time)^(1/3)
Because the Phase 2 CFD-lite model shows a real spatial temperature gradient
(Dirichlet at tank base, Neumann at top), the scalar-T assumption is validated
**after** the optimization, not inside the DoE loop->so we follow these steps :
1. Take the optimum setpoint T* from the desirability optimization.
2. Run the Phase 2 spatial model **once** with T* as boundary condition.
3. Extract the resulting temperature gradient across tank height (ΔT, base vs top)
   over the fermentation window.



Removing the maintenance term from ethanol production had a larger-than-expected
effect on stuck/complete outcomes (5/15 → 3/15 stuck points), due to an indirect
feedback: lower ethanol accumulation reduces the ethanol-driven death rate (kd),
allowing viable biomass to persist longer and consume more sugar before dying
illustrating the coupled nature of the ODE system beyond the directly edited term.

##################################################################################
########## SUMMARY OF RSM for responce =peak rate #####################################
rsm(formula = peak_rate ~ SO(x1, x2, x3), data = design_coded)



# COEFFIENT TABLE(for each term,If I keep all the other 8 terms in the model, does this one term in particular help explain the peak_rate?)


               Estimate  Std. Error  t value  Pr(>|t|)    
(Intercept)  2.00459733  0.01883627 106.4222 1.389e-09 ***
x1           0.75680998  0.01153481  65.6109 1.557e-08 ***
x2           0.00133986  0.01153481   0.1162 0.9120484    
x3           0.80019664  0.01153481  69.3723 1.179e-08 ***
x1:x2        0.00041375  0.01631269   0.0254 0.9807460    
x1:x3        0.32579468  0.01631269  19.9719 5.816e-06 ***
x2:x3        0.00076199  0.01631269   0.0467 0.9645516    
x1^2         0.14144130  0.01697879   8.3305 0.0004075 ***
x2^2         0.00421467  0.01697879   0.2482 0.8138279    
x3^2        -0.03023040  0.01697879  -1.7805 0.1351113    
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Multiple R-squared:  0.9995,	Adjusted R-squared:  0.9985 
F-statistic:  1066 on 9 and 5 DF,  p-value: 1.153e-07

# Analysis of Variance Table (ANOVA ->ASKS FOR GROUPS OF TERMS)
# is it worth having ALL the category (all 3 terms together) in the model, or would #I do just as well without any of them?
#BUT WE NEED THE COEFFICIENT TABLE BECAUSE SOMETIMES TWO OF THE TERMS MATTER MORE
# FO = First Order x1, x2, x3.
# TWI = Two-Way Interactions  x1:x2, x1:x3, x2:x3.
# Q = Pure Quadratic  x1², x2², x3².




Response: peak_rate
                Df Sum Sq Mean Sq    F value    Pr(>F)
FO(x1, x2, x3)   3 9.7046  3.2349 3.0391e+03 1.433e-08
TWI(x1, x2, x3)  3 0.4246  0.1415 1.3296e+02 3.459e-05
PQ(x1, x2, x3)   3 0.0802  0.0267 2.5114e+01  0.001924
Residuals        5 0.0053  0.0011                     
Lack of fit      3 0.0053  0.0018 1.4393e+29 < 2.2e-16
Pure error       2 0.0000  0.0000                     

Stationary point of response surface:
         x1          x2          x3 
-2.48650992 -0.02208355 -0.16393906  # extrapolation 

Stationary point in original units:
       Ti      Sub0        N0 
 # 6.08094 264.66875 223.60609 

Eigenanalysis:
eigen() decomposition
$values
[1]  0.239734628  0.004214441 -0.128523497

$vectors
          [,1]          [,2]         [,3]
x1 0.856203253  0.0022131598  0.516634388
x2 0.001587816 -0.9999973743  0.001652349
x3 0.516636689 -0.0005944261 -0.856204519



# OBSERVATIONS
x1 and x3 (linear terms ) have significant results,second-order (quadratic) x1:x3 are significant 

#the 0.75 in x1 means that if x1 was increased from 0 to +1 (or 6 the step we have) the peak rate will be increased ~ ~0.757 g/L/h keeping the other the same 

#the standard errors are similar because of the model we choose (box-behnken) ,the design is balanced 

# t value 

estimate/std. error 

![alt text](image.png)

![alt text](image.png)

![alt text](image.png)

![alt text](image.png)

Call:
rsm(formula = time ~ SO(x1, x2, x3), data = design_coded)

            Estimate Std. Error  t value  Pr(>|t|)    
(Intercept)  232.616     20.112  11.5662 8.476e-05 ***
x1          -172.711     12.316 -14.0234 3.316e-05 ***
x2            47.586     12.316   3.8638  0.011834 *  
x3          -201.413     12.316 -16.3539 1.559e-05 ***
x1:x2        -52.776     17.417  -3.0301  0.029078 *  
x1:x3         94.047     17.417   5.3996  0.002943 ** 
x2:x3        -22.386     17.417  -1.2853  0.255005    
x1^2          27.826     18.129   1.5350  0.185386    
x2^2          53.464     18.129   2.9492  0.031916 *  
x3^2         104.740     18.129   5.7776  0.002185 ** 
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Multiple R-squared:  0.9911,	Adjusted R-squared:  0.9752 
F-statistic: 62.17 on 9 and 5 DF,  p-value: 0.0001345

Analysis of Variance Table

Response: time
                Df Sum Sq Mean Sq    F value    Pr(>F)
FO(x1, x2, x3)   3 581288  193763 1.5968e+02 2.201e-05
TWI(x1, x2, x3)  3  48525   16175 1.3330e+01  0.008048
PQ(x1, x2, x3)   3  49184   16395 1.3511e+01  0.007814
Residuals        5   6067    1213                     
Lack of fit      3   6067    2022 2.5036e+30 < 2.2e-16
Pure error       2      0       0                     

Stationary point of response surface:
        x1         x2         x3 
-17.761755  -7.508892   8.133282 

Stationary point in original units:
        Ti       Sub0         N0 
 -85.57053  152.36661 1053.32823 

Eigenanalysis:
eigen() decomposition
$values
[1] 133.010056  54.542174  -1.521715

$vectors
         [,1]       [,2]       [,3]
x1  0.4485044 -0.2019381  0.8706692
x2 -0.2687317  0.8985992  0.3468468
x3  0.8524242  0.3895387 -0.3487585







![alt text](image.png)