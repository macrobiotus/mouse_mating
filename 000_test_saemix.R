
# 2024-4-11 - Paul Czechowski

# Nonlinear Mixed-Effects Growth Models: A Tutorial Using 'saemix' in R
# Methodology, 2021, Vol. 17(4), 250–270, https://doi.org/10.5964/meth.7061

library("saemix")

# RQ1: Which of the Logistic, Gompertz, or Richards curves models best the
# growth in achievement?

# RQ2: Does Sex have an association with the total growth, rate of approach to
# the upper asymptote, or point of inflection (of the curve found optimal in
# answering the first research question?)

# RQ3: Does adding SES to the model in addition to Sex as a predictor of total
# growth, rate of approach to the upper asymptote, or point of inflection
# improve model fit? If so, how do SES and Sex relate to achievement growth?

# RQ4: Does Sex moderate the association between SES and total growth, rate of
# approach to the upper asymptote, or point of inflection?