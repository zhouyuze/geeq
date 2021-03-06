library(geepack)
library(geeM)

data(ohio)
fit.ohio.geepack <- geese(resp ~ age + smoke + age:smoke, id=id, data=ohio,
                          family=binomial, corstr="exch")
fit.ohio.geem <- geem(resp ~ age + smoke + age:smoke, id=id, data=ohio,
                      family=binomial, corstr="exchangeable")
fit.ohio.geeq <- geeq(resp ~ age + smoke + age:smoke, id=id, data=ohio,
                     method='gee', family=binomial(), corstr="exchangeable")

data(dietox)
dietox$Cu <- as.factor(dietox$Cu)
fit.dietox.geepack <- geeglm(Weight ~ Time + Cu + Cu * Time, id=Pig, data = dietox,
                              family=gaussian,corstr="exch")
fit.dietox.geem <- geem(Weight ~ Time + Cu + Cu * Time, id=Pig, data = dietox,
                        family=gaussian,corstr="exchangeable")
fit.dietox.geeq <- geeq(Weight ~ Time + Cu + Cu * Time, id=Pig, data = dietox,
                       method='gee', family=gaussian,corstr="exchangeable")


